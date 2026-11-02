from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.pipeline import *
import numpy as np
# Data import:

from ovito.data import DataCollection
import numpy
import os
import re


from ovito.data import CutoffNeighborFinder
def ShrinkWrapSimulationBoxFunction(frame: int, data: DataCollection):

    # There's nothing we can do if there are no input particles.
    if not data.particles or data.particles.count == 0:
        return

    # Compute min/max range of particle coordinates.
    coords_min = numpy.amin(data.particles.positions, axis=0)+0.1
    coords_max = numpy.amax(data.particles.positions, axis=0)-0.1

    # Build the new 3x4 cell matrix:
    #   (x_max-x_min  0            0            x_min)
    #   (0            y_max-y_min  0            y_min)
    #   (0            0            z_max-z_min  z_min)
    matrix = numpy.empty((3,4))
    matrix[:,:3] = numpy.diag(coords_max - coords_min)+ np.identity(3)*0.5
    matrix[:, 3] = coords_min

    # Assign the cell matrix - or create whole new SimulationCell object in
    # the DataCollection if there isn't one already.
    data.create_cell(matrix, (False, False, False))


def main(path_in, path_out, cluster):

    pipeline = import_file(path_in)


    pipeline.modifiers.append(IdentifyDiamondModifier())

    pipeline.modifiers.append(ExpressionSelectionModifier(expression = 'StructureType == 0'))

    pipeline.modifiers.append(ClusterAnalysisModifier(
        cutoff = 5, 
        only_selected = True, 
        sort_by_size = True, 
        compute_com = True, 
        compute_gyration = True, 
        cluster_coloring = True))


    data = pipeline.compute()
    pos = data.particles.positions[data.particles['Cluster'][...]==cluster]

    xx = max(pos[:,0]) - min(pos[:,0])
    yy = max(pos[:,1]) - min(pos[:,1])
    zz = max(pos[:,2]) - min(pos[:,2])

    x0 =  (max(pos[:,0]) + min(pos[:,0]))/2
    y0 =  (max(pos[:,1]) + min(pos[:,1]))/2
    z0 =  (max(pos[:,2]) + min(pos[:,2]))/2

    nx = round(round(xx/5.43+0.5)/2+0.5)+1
    ny = round(round(yy/5.43+0.5)/2+0.5)+1
    nz = round(round(zz/5.43+0.5)/2+0.5)+1


    xlo = x0 - nx*5.43
    ylo = y0 - ny*5.43
    zlo = z0 - nz*5.43
    xhi = x0 + nx*5.43
    yhi = y0 + ny*5.43
    zhi = z0 + nz*5.43

    pipeline.modifiers.append(ExpressionSelectionModifier(expression = f'Position.X > {xlo} && Position.X < {xhi} && Position.Y > {ylo} && Position.Y < {yhi} && Position.Z > {zlo} && Position.Z < {zhi}'))


    pipeline.modifiers.append(InvertSelectionModifier())

    pipeline.modifiers.append(DeleteSelectedModifier())

    pipeline.modifiers.append(ShrinkWrapSimulationBoxFunction)


    data = pipeline.compute()

    dx =  data.cell[0,0]
    dy =  data.cell[1,1]
    dz =  data.cell[2,2]
    o = data.cell[:,3]

    del pipeline

    pipeline = import_file(path_in)
    pipeline.modifiers.append(ExpressionSelectionModifier(expression = f'Position.X > {xlo} && Position.X < {xhi} && Position.Y > {ylo} && Position.Y < {yhi} && Position.Z > {zlo} && Position.Z < {zhi}'))

    pipeline.modifiers.append(InvertSelectionModifier())
    pipeline.modifiers.append(DeleteSelectedModifier())
    pipeline.modifiers.append(ShrinkWrapSimulationBoxFunction)
    pipeline.modifiers.append(ReplicateModifier(
        num_x = 1, 
        num_y = 1, 
        num_z = 1))
    pipeline.modifiers.append(IdentifyDiamondModifier())

    def trim_cell(trim):
        expression=  f'Position.X > {trim[0][0]} && Position.X < {trim[0][1]} && Position.Y > {trim[1,0]} && Position.Y < {trim[1][1]} && Position.Z > {trim[2][0]} && Position.Z < {trim[2][1]}'
        pipeline.modifiers[0].expression=expression
        attr =  pipeline.compute().attributes
        n_dia = attr['IdentifyDiamond.counts.CUBIC_DIAMOND_SECOND_NEIGHBOR'] #+ attr['IdentifyDiamond.counts.CUBIC_DIAMOND_FIRST_NEIGHBOR'] 

        return n_dia
       



    trim = np.array([[xlo,xhi],[ylo,yhi], [zlo,zhi]])

    r = [0,0,0]
    for dim in range(3):
        if dim == 0: pipeline.modifiers[-2].num_x =2
        elif dim == 1: pipeline.modifiers[-2].num_y =2
        elif dim == 2: pipeline.modifiers[-2].num_z =2
        
        r = [trim[dim, 0], trim[dim, 1], 1e8]
        for i in np.linspace(-1.25, 1.25,10):
            for j in np.linspace(-1.25, 1.25,10):
                print(j)
                 
                trim_test = trim.copy()
                trim_test[dim, 0]+=i
                trim_test[dim, 1]+=j
                n_dia = trim_cell(trim_test)
                if n_dia < r[-1]: 
                    r=[trim_test[dim,0], trim_test[dim,1], n_dia]

        trim[dim,0] = r[0]; trim[dim,1] = r[1] 
        pipeline.modifiers[-2].num_x =1;
        pipeline.modifiers[-2].num_y =1;
        pipeline.modifiers[-2].num_z =1;
        

        export_file(pipeline, f"outputdim{dim}.xyz", "xyz", columns = ["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z"])
       
    expression=  f'Position.X > {trim[0][0]} && Position.X < {trim[0][1]} && Position.Y > {trim[1,0]} && Position.Y < {trim[1][1]} && Position.Z > {trim[2][0]} && Position.Z < {trim[2][1]}'
    #print(expression)
    pipeline.modifiers[0].expression=expression

    del(pipeline.modifiers[-2])
    pipeline.modifiers.append(SelectTypeModifier(types = {2}))

    # Delete selected:
    pipeline.modifiers.append(DeleteSelectedModifier())



    export_file(pipeline, "output.xyz", "xyz", columns =
      ["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z"])

    del pipeline
    pipeline = import_file("output.xyz")

    overlap_distance = 1.0

    # The user-defined modifier function:
    def modify(frame, data):
        
        # Show this text in the status bar while the modifier function executes
        yield "Selecting overlapping particles"

        # Create 'Selection' output particle property
        selection = data.particles_.create_property('Selection')

        # Prepare neighbor finder
        finder = CutoffNeighborFinder(overlap_distance, data)

        # Iterate over all particles
        for index in range(data.particles.count):

            # Update progress display in the status bar
            yield (index / data.particles.count)
            # Iterate over all nearby particles around the current center particle
            for neigh in finder.find(index):
                # Once we find a neighbor which hasn't been marked yet,
                # mark the current center particle. This test is to ensure that we
                # always select only one of the particles in a close pair.
                if selection[neigh.index] == 0:
                    print(selection[neigh.index])
                    selection[index] = 1
                    break
    pipeline.modifiers.append(modify)

    pipeline.modifiers.append(WrapPeriodicImagesModifier())
    pipeline.modifiers.append(DeleteSelectedModifier())
    export_file(pipeline, path_out, "lammps/data")## columns =
    #["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z"])
    del pipeline

if __name__ == "__main__":
    """
    directory = "Sb_50keV"        # directory to scan
    cluster = 1            # example cluster identifier

    pattern = re.compile(r"data\.lastframe\.(\d+)$")

    for filename in os.listdir(directory):
        match = pattern.match(filename)
        if match:
            i = match.group(1)
            path_in = os.path.join(directory, filename)
            output_name = f"extracted_cluster.{i}_{cluster}"
            path_out = os.path.join(directory, output_name)
            print(f"Processing {filename} → {output_name}")
            main(path_in, path_out, cluster)
    """
    
    main("/home/fernbouv/out_9000_Bi.data", "aaaa", 1)
