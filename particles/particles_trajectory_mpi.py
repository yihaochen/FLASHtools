#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['figure.dpi'] = 150
import matplotlib.pyplot as plt
import yt
import util
import MPI_taskpull2
import numpy as np
import h5py
import random
import pickle
import os
from mpi4py import MPI

from particles.particles_class import Particles
from tools import calcDen0_2015
from particle_filters import deltaden0

me = yt.utilities.physical_constants.mass_electron.v
c  = yt.utilities.physical_constants.speed_of_light.v
e  = yt.utilities.physical_constants.elementary_charge.v
Myr= yt.units.Myr.in_units('s').v
kpc= yt.units.kpc.in_units('cm').v

dir = './'
regex = '*_hdf5_part_????_updated_peak'
#regex = '*_hdf5_part_????'

# Index of the file in which random particles will be selected
rand_findex = 600
nparticles = 200

tags = None
tadds = None

#grid_fields = ["pressure", "grid_level", "entropy"]
grid_fields = []

def rescan(dir, printlist=False):
    files = util.scan_files(dir, regex=regex, walk=True, printlist=printlist)
    return files


if MPI_taskpull2.rank == 0:
    partfiles = rescan(dir)

    # Arbitrarily pick particles in this file
    h5file = h5py.File(partfiles[rand_findex].fullpath, 'r')
    tp = h5file['tracer particles']
    print('Choosing random particles in %s' % partfiles[rand_findex].filename)

    colname = [item[0].decode().strip() for item in h5file['particle names']]

    try:
        mask = tp[:,colname.index('type')] == 1.0
        tp = tp[mask,:]
    except:
        pass
    # Pick nparticles
    rints = sorted(random.sample(range(tp.shape[0]), nparticles))

    findices = {f: colname.index(f) for f in ['tag', 'tadd', 'shok']}

    tags =  tp[rints,colname.index('tag')]
    tadds = tp[rints,colname.index('tadd')]
    shoks = tp[rints,colname.index('shok')]

tags, tadds = MPI_taskpull2.comm.bcast((tags, tadds), root=0)

def init():
    for (tag, tadd, shok, rint) in zip(tags, tadds, shoks, rints):
        print(tag, tadd, shok, rint)


def worker_fn(filepath):
    parts = Particles(tags, tadds, grid_fields=grid_fields)
    parts.read_from_partfile(filepath)

    return parts

def tasks_gen():
    files = rescan(dir)
    for file in files[600:]:
        yield file.fullpath


picklename = "particles_trajectory.pickle"
if not os.path.exists(picklename):
    tasks = tasks_gen()
    results = MPI_taskpull2.taskpull(worker_fn, tasks, initialize=init)
    # Only on rank 0 since others return None
    if results:
        collected = {}
        #parts = Particles(tags, tadds, grid_fields=grid_fields)
        path, parts = results.popitem()
        for value in results.values():
            parts.combine(value)
        parts.convert_to_ndarray()

        pickle.dump(parts, open( picklename, "wb" ))
else:
    if MPI_taskpull2.rank == 0:
        print('%s found, unpickling...' % picklename)
        parts = pickle.load(open(picklename, 'rb'))
        print('Loaded %s' % picklename)

