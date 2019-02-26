import numpy as np
import matplotlib
import sys
import os
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'stixgeneral'
matplotlib.rcParams['figure.dpi'] = 150
import matplotlib.pyplot as plt
import yt
yt.mylog.setLevel("INFO")

from mpi4py import MPI
from mpl_toolkits.axes_grid1 import make_axes_locatable
from particles.particle_filters import *

comm = MPI.COMM_WORLD

yt.enable_parallelism(communicator=comm)

dir = './'

ls = {(0,10): ['solid', 2],
      (10,20): ['dotted', 2],
      (20,30): ['dashed', 1],
      (30,60): ['solid', 1],
      (60,100): ['dotted', 1]}

rmin, rmax = -10, 100

try:
    ind = int(sys.argv[1])
    ts = yt.DatasetSeries(os.path.join(dir,'data/*_hdf5_part_%04d' % ind), parallel=1)
except IndexError:
    ts = yt.DatasetSeries(os.path.join(dir,'data/*_hdf5_part_???0'), parallel=0)


maindir = os.path.join(dir, 'particles_dr_histogram/')

if yt.is_root():
    for subdir in [maindir]:
        if not os.path.exists(subdir):
            os.mkdir(subdir)

    f0 =  '/d/d8/ychen/2016_production_runs/1212_L45_M10_b1_h0_10Myr/data/MHD_Jet_10Myr_hdf5_part_0000'
    ds0 = yt.load(f0)
    ds0.add_particle_filter("metal")
    ad0 = ds0.all_data()
    arr0 = np.argsort(ad0[('metal', 'particle_tag')])
    rr0 = ad0[('metal','particle_position_spherical_radius')].in_units('kpc')[arr0]
else:
    rr0 = None

rr0 = comm.bcast(rr0,root=0)

for ds in ts.piter():
    print(ds.basename)
    ds.add_particle_filter("metal")
    ad = ds.all_data()
    arr = np.argsort(ad[('metal', 'particle_tag')])
    rr = ad[('metal','particle_position_spherical_radius')].in_units('kpc')[arr]
    drr = rr-rr0
    fig = plt.figure(figsize=(5,4))
    ax = fig.add_subplot(111)

    divider = make_axes_locatable(ax)
    ax2 = divider.append_axes("bottom", 1.1, pad=0.08, sharex=ax)

    for rlim in ls.keys():
        mask = np.logical_and(rr0 > rlim[0], rr0 < rlim[1])
        #ax.scatter(drr[mask], rr0[mask], s=1, lw=0)
        ax.plot(drr[mask], rr0[mask], 'o', markersize=0.5, alpha=0.7, rasterized=True)

        null = ax2.hist(drr[mask], bins=rmax-rmin, range=(rmin, rmax),
                        normed=True, histtype='step',
                        linestyle=ls[rlim][0], linewidth=ls[rlim][1],
                        cumulative=True,
                        label=r'%2i < $r_0$ < %2i' % rlim, alpha=0.9)

    ax.text(rmax-20, rmax-20, '%i Myr' % ds.current_time.in_units('Myr'))
    ax.set_ylabel(r'initial radius $r_0$ (kpc)')
    ax.set_ylim(0, rmax)
    ax.axvline(x=0, lw=1, ls=':', color='k')
    ax.tick_params(labelbottom=False, direction='in')

    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles[::-1], labels[::-1], fontsize=8, loc=4)
    ax2.axvline(x=0, lw=1, ls=':', color='k')
    ax2.grid(ls='--', alpha=0.5)

    ax2.set_xlabel(r'radial displacement $\Delta r$ (kpc)')
    ax2.set_xlim(rmin, rmax)
    ax2.set_ylabel('cumulative fraction')
    ax2.set_ylim(0.6, 1.0)
    plt.tight_layout()
    fig.savefig(os.path.join(maindir, 'particles_dr_%s.pdf' % ds.basename[-4:]), bbox_inches='tight')
    #fig.savefig(os.path.join(maindir, 'particles_dr_%s.png' % ds.basename[-4:]))
