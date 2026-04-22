import numpy as np
import argparse
from scipy.interpolate import splrep, splev
from ase.data import chemical_symbols, atomic_numbers

def cutoff_function(r, rmin, rmax):
    # Perriot polynomial cutoff
    if r < rmin:
        y = 1.0
    elif r > rmax:
        y = 0.0
    else:
        chi = (r - rmin) / (rmax - rmin)
        y = 1.0 - chi**3 * (6.0*chi**2 - 15*chi + 10.0)
    return y


def zbl(Z1, Z2, r):
    """ ZBL energy in eV and derivative in eV/Å, r input in Å """
    eps = 8.854187817e-12
    e = 1.60217657e-19

    a = 0.46848 / (Z1**0.23 + Z2**0.23)
    x= r/a
    phi = 0.18175 * np.exp(-3.1998 * x) + 0.50986 * np.exp(-0.94229 * x)
    phi += 0.28022 * np.exp(-0.4029 * x) + 0.02817 * np.exp(-0.20162 * x)
    r = r * 1e-10  # Å --> m
    E = Z1 * Z2 * e * phi / (4.0 * np.pi * eps * r)

    # derivative
    # (although we don't need it, PARCAS uses spline derivative nowadays too)
    Dphi = (-3.1998 * 0.18175 * np.exp(-3.1998 * x) - 0.94229*0.50986 * \
            np.exp(-0.94229 * x) - 0.4029 * 0.28022 * np.exp(-0.4029 * x) - \
            0.20162 * 0.02817 * np.exp(-0.20162 * x)) / (a * 1e-10)
    DE = -E / r + Z1 * Z2 * e * Dphi / (4.0 * np.pi * eps * r)
    DE = DE * 1e-10

    return E, DE


if __name__ == '__main__':

    default = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=default)
    parser.add_argument('z', nargs=2, help='Z1 and Z2 OR symbol1 symbol2')
    parser.add_argument('-rc', nargs=2, default=[None, None],
                        help='apply cutoff, give distance and width '
                        '(cutoff interval is [rc - width, rc])')
    parser.add_argument('--zbl', action='store_true', default=False,
                        help='also make ZBL data')
    parser.add_argument('--spline', default=False, action='store_true',
                        help='do and save cubic spline to dmol data')
    parser.add_argument('-s', default=0.0, type=float,
                        help='optional smoothing parameter for spline')
    parser.add_argument('-dd', default='/run/user/1290558/gvfs/smb-share:domain=ATKK,server=poiju1.ad.helsinki.fi,share=group,user=fernbouv/h528/acclab/fernbouv/my_old_machine/Desktop/dimer_fit/dmol/dmol-data',
                        help='path to dmol data folder')
    args = parser.parse_args()

    # support any stupid mix of inputs, e.g. Zr W, Zr 74, 40 74
    _1, _2 = args.z
    try:
        Z1 = int(_1)
        s1 = chemical_symbols[Z1]
    except ValueError:
        Z1 = atomic_numbers[_1]
        s1 = _1
    try:
        Z2 = int(_2)
        s2 = chemical_symbols[Z2]
    except ValueError:
        Z2 = atomic_numbers[_2]
        s2 = _2

    # read dmol data
    r_dmol = []
    E_dmol = []
    if Z1 > Z2:
        Z1, Z2 = Z2, Z1  # dmol data files are always Z1 < Z2
        s1, s2 = s2, s1
    dmolfile = args.dd + '/energies.' + str(Z1) + '.' + str(Z2)
    # inconsistent number of columns, so cannot use np.loadtxt
    with open(dmolfile, 'r') as file:
        for line in file:
            line = line.split()
            r_dmol.append(float(line[0]))
            E_dmol.append(float(line[1]))

    # subtract last energy offset
    E_dmol = [e - E_dmol[-1] for e in E_dmol]

    # remove last two points (r=100, 1000)
    while r_dmol[-1] > 10:
        del r_dmol[-1]
        del E_dmol[-1]

    # sometimes first point is E=0 for some reason, remove it
    if E_dmol[0] < 1.0:
        del r_dmol[0]
        del E_dmol[0]

    # do spline interpolation
    if args.spline:
        spl = splrep(r_dmol, E_dmol, k=3, s=args.s)
        x_spline = np.linspace(min(r_dmol), max(r_dmol), 1000)
        y_spline = splev(x_spline, spl)
        y_deriv = splev(x_spline, spl, der=1)

    # apply cutoff
    if args.rc[0] is not None:
        rmax = float(args.rc[0])
        rmin = rmax - float(args.rc[1])
        for i in range(len(r)):
            E_dmol[i] *= cutoff_function(r[i], rmin, rmax)
        if args.spline:
            y_spline = list(y_spline)
            for i in range(len(x_spline)):
                y_spline[i] *= cutoff_function(x_spline[i], rmin, rmax)
    else:
        rmax = 1e6

    # write data
    dmolout = f'dmol_{s1}-{s2}.dat'
    with open(dmolout, 'w') as file:
        for i in range(len(r_dmol)):
            print(r_dmol[i], E_dmol[i], file=file)

    if args.spline:
        splineout = f'spline_dmol_{s1}-{s2}.dat'
        np.savetxt(splineout, np.transpose([x_spline, y_spline, y_deriv]))

    # ZBL
    if args.zbl:
        r = np.linspace(r_dmol[0], r_dmol[-1], 1000)
        E, dEdr = zbl(Z1, Z2, r)
        zblout = f'zbl_{s1}-{s2}.dat'
        np.savetxt(zblout, np.transpose([r, E, dEdr]))
