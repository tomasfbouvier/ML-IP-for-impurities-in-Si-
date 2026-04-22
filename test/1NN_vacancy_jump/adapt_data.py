from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.pipeline import *


def adapt_data(file):

    # Data import:
    pipeline = import_file(file, atom_style = 'atomic')


    # Displacement vectors:
    mod = CalculateDisplacementsModifier()
    mod.reference = FileSource()
    pipeline.modifiers.append(mod)
    mod.reference.load('initial.data', atom_style = 'atomic')

    # Expression selection:
    #pipeline.modifiers.append(ExpressionSelectionModifier(expression = 'DisplacementMagnitude < 0.01'))
    #pipeline.modifiers.append(DeleteSelectedModifier())
    data=pipeline.compute()


    export_file(pipeline, "final.data", "lammps/dump",
        columns = ["Particle Identifier", "Position.X", "Position.Y", "Position.Z"])

    with open('final.data', 'r') as fr:
        lines = fr.readlines()
        ptr = 1
        with open('final.data', 'w') as fw:
            for line in lines:
                if ptr>9 or ptr==4:
                    fw.write(line)
                ptr += 1
    return
