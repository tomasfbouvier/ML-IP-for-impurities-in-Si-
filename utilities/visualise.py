import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import argparse
from scipy.interpolate import griddata
from tabgap.read import *
from cycler import cycler

# TODO support energy grids directly from tabgap class


def plot_2b(data, save_data=False):
    """ plot pair potential(s) """

    if len(data['energies_2b']) == 0:
        return

    fig, ax = plt.subplots(nrows=2, figsize=(5, 9))
    colormap = plt.cm.nipy_spectral
    colours = [colormap(i) for i in np.linspace(0, 0.96, len(data['energies_2b']))]
    epot_min = []
    for i, (pair, epot) in enumerate(data['energies_2b'].items()):
        ax[0].plot(data['distances_2b'], epot, label=pair, color=colours[i])
        ax[1].plot(data['distances_2b'], epot, label=pair, color=colours[i])
        epot_min.append(min(epot))
        if save_data:
            datfile = f'pair_potential_{"".join(pair)}.dat'
            np.savetxt(datfile, np.array([data['distances_2b'], epot]).T)
            print(f' saved {datfile}')

    ax[0].grid()
    ax[1].grid()
    ax[0].legend(fontsize=6)
    ax[1].legend(fontsize=6)
    ax[0].set_yscale('log')
    ax[0].set_ylim(ymin=0.1)
    ax[1].set_ylim(ymin=min(epot_min) - 1, ymax=10)
    ax[1].set_xlabel('Interatomic distance (Å)')
    ax[0].set_ylabel('Pair potential (eV)')
    ax[1].set_ylabel('Pair potential (eV)')

    #plt.savefig('pair_potentials.pdf', bbox_inches='tight')
    plt.savefig('pair_potentials.png', bbox_inches='tight', dpi=300)
    plt.show()
    print(' saved pair_potentials.pdf and .png')


def plot_eam(data, save_data=False):
    """ plot embedding energy and pair density functions """

    fig, ax = plt.subplots(nrows=2, figsize=(5, 9))

    colormap = plt.cm.nipy_spectral
    colours = [colormap(i) for i in np.linspace(0, 0.96, len(data['energies_eam']))]
    for i, (element, F) in enumerate(data['energies_eam'].items()):
        ax[0].plot(data['density'], F, label=element, color=colours[i])
        if save_data:
            datfile = f'embedding_energy_{element}.dat'
            np.savetxt(datfile, np.array([data['density'], F]).T)
            print(f' saved {datfile}')

    colours = [colormap(i) for i in np.linspace(0, 0.96, len(data['pair_densities']))]
    for i, (pair, dens) in enumerate(data['pair_densities'].items()):
        #if pair!="SiSb": continue
        ax[1].plot(data['distances_2b'], dens, label=pair, color=colours[i])
        if save_data:
            datfile = f'pair_density_{"".join(pair)}.dat'
            np.savetxt(datfile, np.array([data['distances_2b'], dens]).T)
            print(f' saved {datfile}')

    ax[0].grid()
    ax[1].grid()
    ax[0].legend(fontsize=6)
    ax[1].legend(fontsize=6)
    ax[0].set_xlabel('EAM density (arb. units)')
    ax[0].set_ylabel('Embedding energy (eV)')
    ax[1].set_xlabel('Interatomic distance (Å)')
    ax[1].set_ylabel('Pair density')

    #plt.savefig('eam_functions.pdf', bbox_inches='tight')
    plt.savefig('eam_functions.png', bbox_inches='tight', dpi=300)
    plt.show()
    print(' saved eam_functions.pdf and .png')


def plot_3b(tabgap_data, triplet=None, save_pngs=False,
            elev=20, azim=35):
    """ Make an animation of the pure 3b energy landscape """

    tg = tabgap_data
    if triplet is None:
        triplets = list(tg['energies_3b'].keys())
    else:
        raise NotImplementedError('...')

    # hard-coded parameters
    ngrid = 100  # keep RAM-use reasonable
    colourmap = 'inferno'
    shading = False

    for t, triplet in enumerate(triplets):
        print(f' rendering triplet {t+1}/{len(triplets)} {triplet}...')
        fig = plt.figure(figsize=plt.figaspect(1.))
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        mp4_file = 'threebody_potential_' + ''.join(triplet) + '.mp4'

        # read-in energy grid is n+2 due to spline filtering -- TODO check/fix rij min/max?
        rij_list = np.linspace(tg['r3b'][triplet][0], tg['r3b'][triplet][1], tg['n3b'][triplet][0]+2)
        dr = rij_list[-1] - rij_list[-2]

        def animate(i):
            rij = rij_list[i]
            data = []
            nstart = i * (tg['n3b'][triplet][1]+2) * (tg['n3b'][triplet][2]+2)
            # print(i, nstart)
            rik = []
            costheta = []
            energies = []
            c = 0
            for j in range(tg['n3b'][triplet][1]+2):
                for k in range(tg['n3b'][triplet][2]+2):
                    rik.append(tg['r3b'][triplet][0] + j * dr)
                    costheta.append(-1 + k * 2 / (tg['n3b'][triplet][2] + 2))  # range is 2 (-1 to 1)
                    energies.append(tg['energies_3b'][triplet][nstart + c])
                    c += 1
            theta = np.rad2deg(np.arccos(costheta))

            # actual grid to plot
            rik_grid = np.linspace(min(rik), max(rik), ngrid)
            theta_grid = np.linspace(min(theta), max(theta), ngrid)
            rik_grid, theta_grid = np.meshgrid(rik_grid, theta_grid)

            # interpolate energies on grid
            E_grid = griddata((rik, theta), energies, (rik_grid, theta_grid), method='linear')

            ax.clear()
            ax.plot_surface(rik_grid, theta_grid, E_grid,
                            rstride=5, cstride=5,
                            cmap=plt.get_cmap(colourmap),
                            edgecolors='k',
                            shade=shading)

            title = '-'.join(triplet)
            title += f': $r_{{ij}} = {round(rij, 2)}$ Å constant'
            ax.set_title(title)
            ax.set_xlabel(r'$r_{ik}$ $\mathrm{\AA}$')
            ax.set_ylabel(r'$\theta_{ijk}$ (degrees)')
            ax.set_zlabel(r'Energy of $i$ (eV)')

            elim = (-1.1,1.1) #(min(tg['energies_3b'][triplet])-1, max(tg['energies_3b'][triplet])+1)
            ax.set_zlim(elim)

            # adjust view angle
            # ax.view_init(elev=20, azim=315)
            ax.view_init(elev=elev, azim=azim)

            if save_pngs:
                ttt = ''.join(triplet)
                plt.savefig(f'threebody_{ttt}-{i}.png', dpi=300, bbox_inches='tight')

        ani = anim.FuncAnimation(fig, animate, frames=len(rij_list), interval=400, repeat=False)
        ani.save(mp4_file)
        print(f' saved {mp4_file}')
        # plt.show()


if __name__ == '__main__':

    default = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=default)
    parser.add_argument('file',
                        help='.tabgap or .eam.fs or .eam.alloy file to visualise/animate.')
    parser.add_argument('--data', action='store_true', default=False,
                        help='save plottable data in files (only pair and EAM data).')
    args = parser.parse_args()

    recognised_file_ext =['.tabgap', '.tabgap.h5', '.eam.fs', '.eam.alloy']

    if '.tabgap.h5' in args.file:
        print(f'Reading and plotting data from .tabgap.h5 file {args.file}')
        data = read_tabgap_potfile_h5(args.file)
        plot_2b(data, save_data=args.data)
        plot_3b(data)
    elif '.tabgap' in args.file:
        print(f'Reading and plotting data from .tabgap file {args.file}')
        data = read_tabgap_potfile(args.file)
        plot_2b(data, save_data=args.data)
        plot_3b(data)
    elif '.eam.fs' in args.file:
        print(f'Reading and plotting data from eam.fs file {args.file}')
        data = read_eam_potfile(args.file, style='eam/fs')
        plot_2b(data, save_data=args.data)
        plot_eam(data, save_data=args.data)
    elif '.eam.alloy' in args.file:
        print(f'Reading and plotting data from eam.alloy file {args.file}')
        data = read_eam_potfile(args.file, style='eam/alloy')
        plot_2b(data, save_data=args.data)
        plot_eam(data, save_data=args.data)
    else:
        print(f'Unrecognised file {args.file}, only following file extensions are recognised:')
        print(recognised_file_ext)
