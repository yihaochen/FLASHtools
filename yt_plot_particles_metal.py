#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import yt
import os
import util
import MPI_taskpull2
import logging
logging.getLogger('yt').setLevel(logging.ERROR)
import matplotlib.pyplot as plt
import numpy as np

def metal(pfilter, data):
    filter = (data[("all", "particle_type")] == 2.0)
    return filter

yt.add_particle_filter("metal", function=metal, filtered_type='all', requires=["particle_type"])

dir = './'
regex = 'MHD_Jet*_hdf5_plt_cnt_[0-9][0-9][0-9]0'
files = None
figuredir = 'particles_metal/y_z_r0'


def rescan(printlist=False):
    files = util.scan_files(dir, regex=regex, walk=True, printlist=printlist, reverse=False)
    return files

# Load the inital file to identify the inital locations of the metal particles
f0 =  './MHD_Jet_hdf5_plt_cnt_0000'
ds = yt.load(f0)
ds.add_particle_filter("metal")
ad = ds.all_data()
arr = np.argsort(ad[('metal', 'particle_tag')])
rr0 = ad[('metal','particle_position_spherical_radius')].in_units('kpc')[arr]

# Filter the central slab
xx=ad[('metal', 'particle_posx')].in_units('kpc')[arr]
filtr = np.abs(xx)<10

def worker_fn(file):
    ds = yt.load(file.fullpath)
    ds.add_particle_filter("metal")

    ad = ds.all_data()
    arr = np.argsort(ad[('metal', 'particle_tag')])
    yy=ad[('metal', 'particle_posy')].in_units('kpc')[arr]
    zz=ad[('metal', 'particle_posz')].in_units('kpc')[arr]

    plt.figure(figsize=(12,10))
    sc=plt.scatter(yy[filtr][:], zz[filtr][:], c=rr0[filtr][:], s=8, lw=0.1, vmax=50, vmin=0, cmap='gnuplot2_r', alpha=0.9)
    cb=plt.colorbar(sc)
    cb.set_label('initial r (kpc)')
    plt.xlim(-40,40)
    plt.ylim(-40,40)
    plt.xlabel('y (kpc)')
    plt.ylabel('z (kpc)')

    plt.tight_layout()
    plt.savefig(os.path.join(dir,figuredir,file.filename), dpi=72)

    return ds.basename[-4: ]


def tasks_gen(i0):
    #sys.stdout.write("Plotting %s... \n" % (file.filename))
    files = rescan(False)[i0:]
    for file in files:
        yield file


def init():
    if not os.path.exists(os.path.join(dir,figuredir)):
        os.mkdir(os.path.join(dir, figuredir))

i0 = int(sys.argv[1]) if len(sys.argv) > 1 else 0
tasks = tasks_gen(i0)

results = MPI_taskpull2.taskpull(worker_fn, tasks, initialize=init, print_result=True)
