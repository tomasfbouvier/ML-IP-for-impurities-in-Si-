#!/usr/bin/python3
# -*- coding: utf-8 -*-
# -------------------------------------------------------------------------- #
# Filename:       run_gap_fit.py
# Author:         Jesper Byggmästar <jesper.byggmastar@helsinki.fi>
# Created:        27.03.2019 11:11
# Last modified:  22.01.2021  8:55
# -------------------------------------------------------------------------- #
import numpy as np
import argparse
import importlib
import subprocess
import ase.io
import os
import shutil
import time
import pathlib

# script that runs GAP-training with QUIP and calculates training and (optionally) test errors


def generate_run_cmd(external, general, descriptors):

    # executable path first
    try:
        cmd = external['gap_fit']
    except KeyError:
        print(" gap_fit executable not found!! Put the path in"
              " external['gap_fit'] in the input script")
        exit()

    for key, value in external.items():
        if key != 'gap_fit':
            cmd += ' ' + key + '=' + str(value)

    for key, value in general.items():
        cmd += ' ' + key + '=' + str(value)

    cmd += ' gap={'

    for des, desparam_group in descriptors.items():   # existing logic
        for desparam in desparam_group:          # new outer loop (list)
            cmd += des
            for key, value in desparam.items():
                cmd += ' ' + key + '=' + str(value)
            cmd += ' : '

    cmd = cmd[:-3]
    cmd += '}'

    return cmd


def run_quip(quip, xyzin, xyzout, potfile):

    # run quip
    run_cmd = quip + ' E=T F=T V=T local=T atoms_filename={}'.format(xyzin)
    run_cmd += ' param_filename={}'.format(potfile)
    with open('quip.out', 'w') as qout:
        subprocess.call(run_cmd, stdout=qout, shell=True)

    # get xyz file
    run_cmd = "grep AT quip.out | sed 's/AT//'"
    with open(xyzout, 'w') as qxyz:
        subprocess.call(run_cmd, stdout=qxyz, shell=True)


def gap_errors(gapxyz, dftxyz):
    """ calculate RMS errors from xyz files """

    gap_atoms = ase.io.read(gapxyz, ':')
    dft_atoms = ase.io.read(dftxyz, ':')
    E_gap = [atoms.get_potential_energy() / len(atoms) for atoms in gap_atoms]
    E_dft = [atoms.get_potential_energy() / len(atoms) for atoms in dft_atoms]
    E_gap = np.asarray(E_gap)
    E_dft = np.asarray(E_dft)

    error2 = (E_dft - E_gap)**2
    E_rms = np.sqrt(np.mean(error2)) * 1e3  # meV
    E_std = np.sqrt(np.var(error2)) * 1e3

    # NOTE all atoms/symbols (see validate_gap.py for specific symbol)
    F_gap = [atoms.arrays['force'][i] for atoms in gap_atoms for i in range(len(atoms))]
    F_dft = [atoms.arrays['force'][i] for atoms in dft_atoms for i in range(len(atoms))]
    F_gap = np.asarray(F_gap)
    F_dft = np.asarray(F_dft)

    error2 = (F_dft - F_gap)**2
    F_rms = np.sqrt(np.mean(error2))
    F_std = np.sqrt(np.var(error2))

    # check if we find virials
    s_rms = 999.999
    s_std = 999.999
    do_virials = True
    for atoms in gap_atoms:
        if 'virial' not in atoms.info:
            do_virials = False
            break
    for atoms in dft_atoms:
        if 'virial' not in atoms.info:
            do_virials = False
            break
    if do_virials:
        s_gap = [atoms.info['virial'] for atoms in gap_atoms]
        s_dft = [atoms.info['virial'] for atoms in dft_atoms]
        s_gap = np.asarray(s_gap)
        s_dft = np.asarray(s_dft)
        error2 = (s_dft - s_gap)**2
        s_rms = np.sqrt(np.mean(error2))
        s_std = np.sqrt(np.var(error2))

    return E_rms, E_std, F_rms, F_std, s_rms, s_std


def test_gap_db(quip, dbfile, xyzout, errorfile, potfile):
    run_quip(quip, dbfile, xyzout, potfile)
    errors = gap_errors(xyzout, dbfile)
    errors = [round(e, 3) for e in errors]
    with open(errorfile, 'a') as ef:
        print('{:<30} {} {} {} {} {} {}'.format(dbfile, *errors), file=ef,
              flush=True)


def teach_gap():
    default = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=default)
    parser.add_argument('--cmd', default=False, action='store_true',
                        help='only print gap_fit training command and exit')
    parser.add_argument('--nofit', default=False, action='store_true',
                        help='assume GAP exists and only get test/train errors')
    args = parser.parse_args()

    do_full_db = False  # calculate combined RMSEs on full training/test .xyz?

    # get input parameters
    from input_gap_fit import external, general, descriptors, validation

    # check that executables exist (possibly replace ~ with home dir path)
    if not pathlib.Path(external['gap_fit']).expanduser().is_file():
        print(f"gap_fit executable {external['gap_fit']} not found")
        exit()
    if args.cmd:
        run_cmd = generate_run_cmd(external, general, descriptors)
        print(run_cmd)
        exit()
    if not pathlib.Path(validation['quip']).expanduser().is_file():
        print(f"quip executable {validation['quip']} not found")
        exit()

    # prepare output folders
    out_train = 'out-train'
    out_test = 'out-test'
    if os.path.exists(out_train):
        shutil.rmtree(out_train)
    if os.path.exists(out_test):
        shutil.rmtree(out_test)
    os.makedirs(out_train)

    timefile = open('timing.out', 'w')

    # run gap_fit
    if not args.nofit:
        t0 = time.time()
        run_cmd = generate_run_cmd(external, general, descriptors)
        print(run_cmd)
        with open('gap_fit.out', 'w') as tout:
            subprocess.call(run_cmd, stdout=tout, shell=True)
        t1 = time.time()
        print(' gap_fit: {}'.format(t1 - t0), file=timefile, flush=True)

    # test trained gap
    # full training db
    if do_full_db:
        t0 = time.time()
        test_gap_db(validation['quip'], external['at_file'],
                    'gap_' + external['at_file'], 'accuracy.out', general['gp_file'])
        t1 = time.time()
        print(' QUIP GAP {}: {}'.format(external['at_file'], t1 - t0), file=timefile, flush=True)
        os.rename('gap_' + external['at_file'], out_train + '/gap_' + external['at_file'])

    # split training db
    all_frames = ase.io.read(external['at_file'], ':')
    configs = []
    for i, atoms in enumerate(all_frames):
        info = atoms.info
        config = info['config_type']
        if config not in configs:
            configs.append(config)
        xyzname = config + '.xyz'
        ase.io.write(xyzname, atoms, format='extxyz', append=True)
    for config in configs:
        xyz = config + '.xyz'
        t0 = time.time()
        test_gap_db(validation['quip'], xyz, 'gap_' + xyz, 'accuracy.out',
        general['gp_file'])
        t1 = time.time()
        print('  QUIP GAP {}: {}'.format(xyz, t1 - t0), file=timefile, flush=True)
        os.rename(xyz, out_train + '/db_' + xyz)
        os.rename('gap_' + xyz, out_train + '/gap_' + xyz)

    # move/copy training stuff to out folder
    os.rename('gap_fit.out', out_train + '/gap_fit.out')
    os.rename('accuracy.out', out_train + '/accuracy.out')
    shutil.copyfile(external['at_file'], out_train + '/' + external['at_file'])
    shutil.copyfile('input_gap_fit.py', out_train + '/' + \
                    'input_gap_fit.py')
    for filename in os.listdir('.'):
        if filename.startswith('sparse_indices') and filename.endswith('.out'):
            os.rename(filename, out_train + '/' + filename)
    # try:
        # os.rename('sparse_indices_2b.out', out_train + '/sparse_indices_2b.dat')
        # os.rename('sparse_indices_soap.out', out_train + '/sparse_indices_soap.dat')
    # except OSError:
        # pass


    # optional test db
    if 'testxyz' in validation:
        os.makedirs(out_test)
        if do_full_db:
            t0 = time.time()
            test_gap_db(validation['quip'], validation['testxyz'], 'gap_' + \
                        validation['testxyz'], 'accuracy.out', general['gp_file'])
            t1 = time.time()
            print(' QUIP GAP {}: {}'.format(validation['testxyz'], t1 - t0),
                  file=timefile, flush=True)
            os.rename('gap_' + validation['testxyz'],
                      out_train + '/gap_' + validation['testxyz'])

        # split training db
        all_frames = ase.io.read(validation['testxyz'], ':')
        configs = []
        for i, atoms in enumerate(all_frames):
            info = atoms.info
            config = info['config_type']
            if config not in configs:
                configs.append(config)
            xyzname = config + '.xyz'
            ase.io.write(xyzname, atoms, format='extxyz', append=True)
        for config in configs:
            xyz = config + '.xyz'
            t0 = time.time()
            test_gap_db(validation['quip'], xyz, 'gap_' + xyz, 'accuracy.out',
            general['gp_file'])
            t1 = time.time()
            print('  QUIP GAP {}: {}'.format(xyz, t1 - t0), file=timefile, flush=True)
            os.rename(xyz, out_test + '/db_' + xyz)
            os.rename('gap_' + xyz, out_test + '/gap_' + xyz)
        os.rename('accuracy.out', out_test + '/accuracy.out')

    # clean up
    os.remove('quip.out')
    dir = os.getcwd()
    files = os.listdir(dir)
    for file in files:
        if file.endswith('.idx'):
            os.remove(file)


if __name__ == '__main__':
    teach_gap()
