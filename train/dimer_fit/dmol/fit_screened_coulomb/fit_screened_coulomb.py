import numpy as np
import argparse
import scipy.optimize as spo
import matplotlib.pyplot as plt


def zbl(r, Z1, Z2):
    eps = 8.854187817e-12
    e = 1.60217657e-19

    a = 0.46848 / (Z1**0.23 + Z2**0.23)
    x = r / a
    phi = 0.18175 * np.exp(-3.1998 * x) + 0.50986 * np.exp(-0.94229 * x)
    phi += 0.28022 * np.exp(-0.4029 * x) + 0.02817 * np.exp(-0.20162 * x)
    r = r * 1e-10  # Å --> m
    E = Z1 * Z2 * e * phi / (4.0 * np.pi * eps * r)
    return E


def screened_coulomb(r, c1, c2, c3, c4, c5, c6):
    eps = 8.854187817e-12
    e = 1.60217657e-19

    a = 0.46848 / (Z1**0.23 + Z2**0.23)
    x = r / a
    phi = c1 * np.exp(-c2 * x) + c3 * np.exp(-c4 * x) + c5 * np.exp(-c6 * x)
    r = r * 1e-10  # Å --> m
    E = Z1 * Z2 * e * phi / (4.0 * np.pi * eps * r)
    return E


def screened_coulomb_cutoff(r, c1, c2, c3, c4, c5, c6):
    eps = 8.854187817e-12
    e = 1.60217657e-19

    # apply smooth Perriot polynomial cutoff
    cutoff_array = []
    for rr in r:
        if rr > r2:
            cutoff = 0.0
        elif rr < r1:
            cutoff = 1.0
        else:
            chi = (rr - r1) / (r2 - r1)
            cutoff = 1.0 - chi**3 * (6.0*chi**2 - 15*chi + 10.0)
        cutoff_array.append(cutoff)
    cutoff_array = np.asarray(cutoff_array)

    a = 0.46848 / (Z1**0.23 + Z2**0.23)
    x = r / a
    phi = c1 * np.exp(-c2 * x) + c3 * np.exp(-c4 * x) + c5 * np.exp(-c6 * x)
    r = r * 1e-10  # Å --> m
    E = Z1 * Z2 * e * phi / (4.0 * np.pi * eps * r)
    return E * cutoff_array


def deriv_screened_coulomb(r, c1, c2, c3, c4, c5, c6):
    # NOTE NOTE does not include cutoff
    eps = 8.854187817e-12
    e = 1.60217657e-19

    a = 0.46848 / (Z1**0.23 + Z2**0.23)
    x = r / a
    phi = c1 * np.exp(-c2 * x) + c3 * np.exp(-c4 * x) + c5 * np.exp(-c6 * x)
    Dphi = (-c2 * c1 * np.exp(-c2 * x) - c4 * c3 * np.exp(-c4 * x) - \
            c6 * c5 * np.exp(-c6 * x)) / (a * 1e-10)
    r = r * 1e-10  # Å --> m
    E = Z1 * Z2 * e * phi / (4.0 * np.pi * eps * r)
    DE = -E / r + Z1 * Z2 * e * Dphi / (4.0 * np.pi * eps * r)
    DE = DE * 1e-10
    return DE


def fit_sc_reppot(Z1, Z2, reppotfile, plotdatafile=None, plot=False):

    print('\nfitting {} for {}-{} repulsion'.format(reppotfile, Z1, Z2))

    # prepare data, slice to fit desired interval
    data = np.loadtxt(reppotfile)
    x = list(data[:, 0])
    y = list(data[:, 1])
    for i in reversed(range(len(x))):
        if x[i] < args.rfit[0] or x[i] > args.rfit[1]:
            del x[i]
            del y[i]

    # fit screened Coulomb to data
    c1_guess = 0.2
    c2_guess = 1.0
    c3_guess = 0.2
    c4_guess = 1.0
    c5_guess = 0.2
    c6_guess = 1.0
    p0 = [c1_guess, c2_guess, c3_guess, c4_guess, c5_guess, c6_guess]
    # uncertainties (less weight on extremely large energies)
    sigma = range(len(x)+1, 1, -1)
    coeff, cov = spo.curve_fit(screened_coulomb, x, y, p0=p0, sigma=sigma,
                               maxfev=20000, bounds=(0, np.inf))
    coeff_err = np.sqrt(np.diag(cov))
    print(coeff)
    print(coeff_err)

    if plot:
        # plot log
        r = np.linspace(data[:, 0][0], data[:, 0][-1], 1000)
        plt.plot(r, zbl(r, Z1, Z2), label='ZBL')
        plt.plot(data[:, 0], data[:, 1], label='input data')
        plt.plot(x, y, '.-', label='sliced data for fitting')
        plt.plot(r, screened_coulomb(r, *coeff), label='fit')
        plt.plot(r, screened_coulomb_cutoff(r, *coeff), label='fit with cutoff')
        if plotdatafile is not None:
            for dfile in plotdatafile:
                data2 = np.loadtxt(dfile)
                plt.plot(data2[:, 0], data2[:, 1], 'o--', mfc='None', label=dfile)
        plt.yscale('log')
        plt.ylim([1, 1e8])
        plt.xlim([0, 3])
        plt.grid()
        plt.legend()

        # plot not log
        plt.figure()
        plt.plot(r, zbl(r, Z1, Z2), label='ZBL')
        plt.plot(data[:, 0], data[:, 1], label='input data')
        plt.plot(x, y, '.-', label='sliced data for fitting')
        plt.plot(r, screened_coulomb(r, *coeff), label='fit')
        plt.plot(r, screened_coulomb_cutoff(r, *coeff), label='fit with cutoff')
        if plotdatafile is not None:
            for dfile in plotdatafile:
                data2 = np.loadtxt(dfile)
                plt.plot(data2[:, 0], data2[:, 1], 'o--', mfc='None', label=dfile)
        plt.ylim([-10, 1000])
        plt.xlim([0, 5])
        plt.grid()
        plt.legend()
        plt.show()

    # save final data with cutoff
    outdata = 'screened_coulomb_' + str(Z1) + '-' + str(Z2) + '.dat'
    # NOTE upper limit hard-coded cutoff
    r = np.linspace(0.004, r2, 400)
    E = screened_coulomb_cutoff(r, *coeff)
    DE = deriv_screened_coulomb(r, *coeff)
    np.savetxt(outdata, np.transpose(np.array([r, E, DE])))
    print('fitted data in ', outdata)



default = argparse.ArgumentDefaultsHelpFormatter
parser = argparse.ArgumentParser(formatter_class=default)
parser.add_argument('data', help='file with data to fit reppot to')
parser.add_argument('-z', nargs=2, type=int, required=True,
                    help='atomic numbers Z1 and Z2')
parser.add_argument('-rcut', type=float, default=[1.0, 2.2], nargs=2,
                    help='cutoff interval')
parser.add_argument('-rfit', type=float, default=[0.04, 1.2], nargs=2,
                    help='only fit to data in this distance interval')
parser.add_argument('-d', nargs='*',
                    help='additional data files to plot (e.g. VASP data)')
args = parser.parse_args()

r1, r2 = args.rcut
Z1, Z2 = args.z

fit_sc_reppot(Z1, Z2, args.data, args.d, plot=True)
