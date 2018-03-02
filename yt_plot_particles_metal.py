#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import yt
import os
import sys
import util
import MPI_taskpull2
import logging
logging.getLogger('yt').setLevel(logging.ERROR)
import matplotlib.pyplot as plt
import numpy as np
from particles.particle_filters import *

dir = './'
regex = 'MHD_Jet*_hdf5_plt_cnt_[0-9][0-9][0-9]0'
files = None
figuredir = 'particles_metal/'
field = 'dr_hist'


def rescan(printlist=False):
    files = util.scan_files(dir, regex=regex, walk=True, printlist=printlist, reverse=True)
    return files

# Load the inital file to identify the inital locations of the metal particles
f0 =  util.scan_files(dir, regex='MHD_Jet*_hdf5_plt_cnt_0000', walk=True)[0].fullpath
ds = yt.load(f0)
ds.add_particle_filter("metal")
ad = ds.all_data()
arr = np.argsort(ad[('metal', 'particle_tag')])
rr0 = ad[('metal','particle_position_spherical_radius')].in_units('kpc')[arr]

# Filter the central slab
xx=ad[('metal', 'particle_posx')].in_units('kpc')[arr]
filtr = np.abs(xx)<30 if field == 'y_z_r0' else True

def worker_fn(file):
    ds = yt.load(file.fullpath)
    ds.add_particle_filter("metal")

    ad = ds.all_data()
    arr = np.argsort(ad[('metal', 'particle_tag')])
    yy=ad[('metal', 'particle_posy')].in_units('kpc')[arr]
    zz=ad[('metal', 'particle_posz')].in_units('kpc')[arr]

    plt.figure(figsize=(6,5))
    if field == 'y_z_r0':
        sc=plt.scatter(yy[filtr][::3], zz[filtr][::3], c=rr0[filtr][::3], s=1, lw=0.0, vmax=80, vmin=0, cmap='gnuplot2_r', alpha=0.9)
        cb=plt.colorbar(sc)
        cb.set_label('initial r (kpc)')
        plt.xlim(-50,50)
        plt.ylim(-50,50)
        plt.xlabel('y (kpc)')
        plt.ylabel('z (kpc)')
    elif field == 'r0_r':
        rr = ad[('metal','particle_position_spherical_radius')].in_units('kpc')[arr]
        rcyl = ad[('metal', 'particle_position_cylindrical_radius')].in_units('kpc')[arr]
        sc = plt.scatter(rr0[filtr][::3], rr[filtr][::3], c=rcyl[filtr][::3], s=1, lw=0.0, vmax=40, vmin=0)
        rline = np.linspace(0,100)
        plt.plot(rline, rline, ':', c='k', alpha=0.5)
        plt.xlim(0,100)
        plt.ylim(0,100)
        plt.xlabel('initial r (kpc)')
        plt.ylabel('r (kpc)')
    elif field == 'r0_dr':
        rr = ad[('metal','particle_position_spherical_radius')].in_units('kpc')[arr]
        rcyl = ad[('metal', 'particle_position_cylindrical_radius')].in_units('kpc')[arr]
        dr = rr-rr0
        sc = plt.scatter(rr0[filtr][::3], dr[filtr][::3], c=rcyl[filtr][::3], s=1, lw=0.0, vmax=40, vmin=0)
        plt.hlines(0, 0, 100, linestyles=':', alpha=0.5)
        plt.xlim(0,100)
        plt.ylim(-5,60)
        plt.xlabel('initial r (kpc)')
        plt.ylabel(r'$\Delta$ r (kpc)')
    elif field == 'dr_hist':
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        rr = ad[('metal','particle_position_spherical_radius')].in_units('kpc')[arr]
        drr = rr - rr0
        fig = plt.figure(figsize=(5,7))

        ax = fig.add_subplot(111)

        divider = make_axes_locatable(ax)
        ax2 = divider.append_axes("bottom", 1.3, pad=0.1, sharex=ax)

        for rlim in [(0,30), (30,60), (60,90)]:
            mask = np.logical_and(rr0 > rlim[0], rr0 < rlim[1])
            ax.scatter(drr[mask], rr0[mask], s=1, lw=0)

            null = ax2.hist(drr[mask], bins=65, range=(-5, 60),
                            normed=True, histtype='step',
                            label=r'%2i < $r_0$ < %2i' % rlim, alpha=0.9, log=True)
        ax.text(45, 85, '%.1f Myr' % ds.current_time.in_units('Myr'))
        ax.set_xlim(-5, 60)
        ax.set_ylim(0, 90)
        ax.set_ylabel(r'initial radius $r_0$ (kpc)')
        ax2.legend(loc=1)
        ax2.set_xlabel(r'$\Delta r$ (kpc)')
        ax2.set_ylim(1E-5, 1E0)

    plt.tight_layout()
    if file.pathname.split('/')[-1] == 'data':
        absfiguredir = os.path.join(os.path.dirname(file.pathname), figuredir)
    else:
        absfiguredir = os.path.join(file.pathname, figuredir)
    plt.savefig(os.path.join(absfiguredir,field,file.filename), dpi=150)

    return ds.basename[-4: ]


def tasks_gen(i0):
    #sys.stdout.write("Plotting %s... \n" % (file.filename))
    files = rescan(False)[i0:]
    for file in files:
        yield file


def init():
    if not os.path.exists(os.path.join(dir,figuredir)):
        os.mkdir(os.path.join(dir, figuredir))
    if not os.path.exists(os.path.join(dir,figuredir,field)):
        os.mkdir(os.path.join(dir,figuredir,field))

i0 = int(sys.argv[1]) if len(sys.argv) > 1 else 0
tasks = tasks_gen(i0)

results = MPI_taskpull2.taskpull(worker_fn, tasks, initialize=init, print_result=True)
