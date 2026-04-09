#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from ase.calculators.vasp import Vasp
import ase.io
import shutil
import os


def write_quip_xyz(filename, atoms, energy, forces=None, stresses=None,
                   pbc=(1, 1, 1), writemode='w', config_type='gap_db'):
    """ Write xyz for input to GAP-training.
    Compared to ase.io.write: only includes needed data & "forces" --> "force" (gap_fit default)
    """
    with open(filename, mode=writemode) as xf:
        # 1st line
        print(len(atoms), file=xf)
        # 2nd line
        line2 = 'Lattice="'
        cell = atoms.get_cell()
        cell = np.ndarray.flatten(np.asarray(cell))
        cell = ' '.join(str(c) for c in cell)
        line2 += cell + '"'
        props = ' Properties=species:S:1:pos:R:3'
        if forces is not None:
            props += ':force:R:3'
        props += ':Z:I:1'
        line2 += props
        pbc_string = '"'
        for p in pbc:
            if p:
                pbc_string += 'T '
            else:
                pbc_string += 'F '
        pbc_string = pbc_string[:-1]  # delete last space
        pbc_string += '"'
        line2 += ' pbc=' + pbc_string
        line2 += ' config_type=' + config_type + ' energy=' + str(energy)
        if stresses is not None:
            line2 += ' virial="'
            stresses = np.ndarray.flatten(stresses)
            stresses = ' '.join(str(s) for s in stresses)
            line2 += stresses + '"'
        print(line2, file=xf)
        # atoms
        for i in range(len(atoms)):
            x, y, z = atoms[i].position
            if forces is None:
                print('{} {} {} {} {}'.
                      format(atoms[i].symbol, x, y, z,
                             atoms[i].number), file=xf)
            else:
                print('{} {} {} {} {} {} {} {}'.
                      format(atoms[i].symbol, x, y, z,
                             forces[i][0], forces[i][1], forces[i][2],
                             atoms[i].number), file=xf)


if __name__ == '__main__':

    ####################################################################
    ################## INPUT TO BE EDITED BEGINS  ######################
    ####################################################################

    ################## SET UP VASP CALCULATOR ##########################
    # run directory
    outdir = 'vasp'

    # INCAR and other parameters
    params = {'xc': 'pbe',            # PBE XC functional
              'pp': 'potpaw_PBE.54',  # PP folder
              'setups': {'W': '_sv',  # choose PAW potential, use _sv or_pv semi-core!
                         'Cr': '_pv',
                         'V': '_sv',
                         'Ta': '_pv'},
              'encut': 900,           # plane-wave energy cutoff
              ## note: careful with kpoints, for same "config" boxes with varying volumes,
              ## it's safer to specify kpts instead of kspacing, so that same k-grid is used for all volumes.
              ## (tip: to get kpts, launch test run with desired kspacing on a small/medium-volume box from config)
              'kspacing': 0.15,       # max spacing between k points (2pi/Å ??)
              #'kpts': (7, 7, 7),      # k-point grid
              'ismear': 0,            # Methfessel-Paxton smearing of order 1
                                      # recommended for metals
              'sigma': 0.03,           # sigma (eV) for smearing
              'ediff': 1.00e-06,          # convergence energy difference for SC loop
              'nelmin': 2,            # minimum number of SCF iterations
              'lasph': False,          # recommended for aspherical charge densities
              'ncore': 64,            # parallelisation parameter, should
                                      # be ~= number of cores per node (optimise!)
              'kpar': 8,              # k-point parallelisation, total n_cores should
                                      # be divisible by kpar, kpar ~= number of nodes
              'prec': 'Accurate',     # recommended for accurate forces/phonons
              'lwave': False,         # whether to write WAVECAR or not
              'nelm': 150,            # maximum number of SC steps
              #'lorbit': 10,           # lorbit >= 10 for magmoms in OUTCAR
              #'isym': 0,# symmetry treatment, turn off (isym=0) in case of problematic boxes (liquids)
              'algo': 'Fast',
              'gga': 'PE',
              'lreal': False
             }

    calc = Vasp(directory=outdir, **params)
    ##########################################################################


    ########################## (RE-)RUNNING INPUT ############################
    xyz_in = 'train.xyz'   # input (multi-frame) xyz file
    index = 300  # start from this frame (first frame = 0)
    end_index = 1000  # stop after this frame (put >= N_frames to do entire file)
    # what to save in output .xyz (energies are always written)
    get_forces = True
    get_stress = True
    # do spin-polarized run?
    do_spinpol = False
    initial_magmoms = 4.0
    ##########################################################################
    ################## INPUT TO BE EDITED ENDS  ##############################
    ##########################################################################

    # initialise output files/folders
    out_xyz = 'db.xyz'  # name of DFT-calculated xyz frames
    vasp_save = 'out_vasp'  # name of folder for saving vasp in/out
    runout_name = 'db.out' # name of file for printing run info
    runout = open(runout_name, 'a')  # always append
    if os.path.exists(out_xyz):
        print('# NOTE: output xyz {} already exists, appending new output'.
              format(out_xyz), flush=True, file=runout)
    if not os.path.exists(vasp_save):
        os.makedirs(vasp_save)
    else:
        print('# NOTE: output folder {} already exists, some output may be '
              'overwritten'.format(vasp_save), flush=True, file=runout)


    ##################### LOOP OVER FRAMES AND RUN VASP ######################
    while index <= end_index:
        print(index)
        try:
            atoms = ase.io.read(xyz_in, index=index)

            # save old energy, if it exists, before attaching calculator
            try:
                config_type = atoms.info['config_type']
            except:
                config_type = None

            try:
                old_energy = atoms.get_potential_energy()
            except:
                old_energy = 0.0

            atoms.set_calculator(calc)

            # magnetism?
            if do_spinpol:
                atoms.set_initial_magnetic_moments([initial_magmoms] * len(atoms))

            try:
                # runs VASP, get the free energy to be consistent with forces
                energy = atoms.get_potential_energy(force_consistent=True)
            except:
                print('# vasp run for frame index {} failed'.format(index),
                      file=runout, flush=True)
                shutil.copyfile(outdir + '/OUTCAR', vasp_save + '/OUTCAR-' + str(index))
                shutil.copyfile(outdir + '/INCAR', vasp_save + '/INCAR-' + str(index))
                shutil.copyfile(outdir + '/vasp.out', vasp_save + '/stdout-' + str(index))
                index += 1
                continue
            forces = None  # initialise in case we don't want them
            stress = None
            if get_forces:
                forces = atoms.get_forces()
            if get_stress:
                stress = atoms.get_stress(voigt=False)
                # QUIP wants the virial tensor in units of eV/Å^3*V_box = eV
                # (still calling it, slightly misleadingly, 'stress')
                # also note sign difference between ASE and QUIP (and VASP)
                # NOTE beware: this may or may not change in future ASE or QUIP!?
                stress = -stress * atoms.get_volume()

            # calculate old-new energy difference per atom
            ediff = 0.0
            if abs(old_energy) > 1e-4:
                ediff = (old_energy - energy) / len(atoms)

            # write output
            print('{:<6} {:<20} {:<20} {:<20} # force={} stress={} config={}'.
                  format(index, energy, old_energy, ediff, get_forces, get_stress,
                         config_type), flush=True, file=runout)
            write_quip_xyz(out_xyz, atoms, energy, forces, stress,
                           writemode='a', config_type=config_type)

            # save some of the vasp in/output
            shutil.copyfile(outdir + '/OUTCAR', vasp_save + '/OUTCAR-' + str(index))
            shutil.copyfile(outdir + '/INCAR', vasp_save + '/INCAR-' + str(index))
            shutil.copyfile(outdir + '/CONTCAR', vasp_save + '/CONTCAR-' + str(index))
            shutil.copyfile(outdir + '/vasp.out', vasp_save + '/stdout-' + str(index))
            #shutil.copyfile(outdir + '/CHG', vasp_save + '/CHG-' + str(index))
            shutil.copyfile(outdir + '/vasprun.xml', vasp_save + '/vasprun.xml-' + str(index))

            index += 1
        except StopIteration:
            print('# reached end of file {} for index {}'.format(xyz_in, index),
                  file=runout)
            exit()
