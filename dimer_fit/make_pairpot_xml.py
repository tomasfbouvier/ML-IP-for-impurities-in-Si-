#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -------------------------------------------------------------------------- #
# Filename:       make_pairpot_xml.py
# Author:         Jesper Byggmästar <jesper.byggmastar@helsinki.fi>
# Created:        21.02.2019 16:04
# Last modified:  06.02.2020 12:24
# -------------------------------------------------------------------------- #
import numpy as np
import argparse

# make QUIP-readable pair potential (as the pair part of a Glue model)

default = argparse.ArgumentDefaultsHelpFormatter
parser = argparse.ArgumentParser(formatter_class=default)
parser.add_argument('-f', nargs='+',
                    help='data file(s) with pair potential data')
parser.add_argument('-s', nargs='+', type=str,
                    help='element symbols of atom pair(s) (two per file)')
args = parser.parse_args()

# sanity check
if len(args.s) != len(args.f) * 2:
    print(' mismatch between number of files ({}) and number of element pairs'
          ' ({})'.format(len(args.f), len(args.s) / 2))
    exit()

# Z periodic table... (missing some rare high-Z elements)
Z = {'H': 1, 'He': 2,
     'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
     'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17,
     'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24,
     'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31,
     'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38,
     'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45,
     'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52,
     'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56,
     'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78,
     'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85,
     'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92,
    }
# periodic table other way around, i.e. {1: 'H'}
symbol = {v: k for k, v in Z.items()}

# initialise atom types
all_Z = [Z[s] for s in args.s]
sorted_Z = sorted(list(set(all_Z)))
at = {}
for i, zz in enumerate(sorted_Z):
    at[zz] = i + 1

# loop over Z and output xml and raw data
with open('pairpot.xml', 'w') as zf:
    # first lines
    label = 'pairpot'  # potential label for QUIP
    print('<{}>'.format(label), file=zf)
    print('<Potential label="{}" init_args="IP Glue"/>'.format(label), file=zf)
    print('<Glue_params n_types="{}">'.format(len(list(set(args.s)))),
                                                  file=zf)

    # per_type data (i.e. many-body part): garbage data that gives 0 for all r
    for i, Z1 in enumerate(sorted_Z):
        print(' <per_type_data atomic_num="{}" type="{}">'
              .format(Z1, at[Z1]), file=zf)
        print('  <!-- {}, Z: {}, many-body data ==> zero energy for'
              ' pair potential! -->'.format(symbol[Z1], Z1), file=zf)

        # density function
        # can be any trash data that doesn't make QUIP crash, like so:
        a = np.linspace(0.5, 5.0, 10)
        rho = [0.2] * len(a)
        edp = [0.0, 0.0]  # end-point derivatives
        print('  <density num_points="{}" density_y1="{}" density_yn="{}">'
              .format(len(a), edp[0], edp[1]), file=zf)
        for j in range(len(a)):
            print('   <point a="{}" rho="{}"/>'.format(a[j], rho[j]), file=zf)
        print('  </density>', file=zf)

        # embedding/glue function. Must give zero energies for pure ZBL
        rho = np.linspace(0, 1.0, 10)
        E = [0.0] * len(rho)  # ZERO!
        print('  <potential_density num_points="{}">'.format(len(rho)), file=zf)
        for j in range(len(rho)):
            print('   <point rho="{}" E="{}"/>'.format(rho[j], E[j]), file=zf)
        print('  </potential_density>', file=zf)
        print(' </per_type_data>', file=zf)

    # pair potential data
    for i in range(len(args.f)):
        s1 = args.s[2 * i]
        Z1 = Z[s1]
        s2 = args.s[2 * i + 1]
        Z2 = Z[s2]
        print(' <per_pair_data type1="{}" type2="{}">'
              .format(at[Z1], at[Z2]), file=zf)
        print('  <!-- {}-{}, Z: {}-{} pair potential -->'
              .format(s1, s2, Z1, Z2), file=zf)

        data = np.loadtxt(args.f[i])
        r = data[:, 0]
        E = data[:, 1]
        edp = [1e50, 0.0]  # end-point derivatives

        print('  <potential_pair num_points="{}" y1="{}" yn="{}">'
              .format(len(r), edp[0], edp[1]), file=zf)
        for j in range(len(r)):
            print('   <point r="{}" E="{}"/>'.format(r[j], E[j]), file=zf)
        print('  </potential_pair>', file=zf)
        print(' </per_pair_data>', file=zf)

    # last line
    print('</Glue_params>', file=zf)
    print('</{}>'.format(label), file=zf)
