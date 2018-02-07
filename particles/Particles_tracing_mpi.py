#!/usr/bin/env python
import yt
import util
import MPI_taskpull2
import logging
logging.getLogger('yt').setLevel(logging.WARNING)
import numpy as np
import h5py
import random
import pickle

from mpi4py import MPI
from particles.particles_class import Particles


dir = './'
regex = 'MHD_Jet*_hdf5_part_????_updated'

# Index of the file in which random particles will be selected
rand_findex = 600
nparticles = 100

tags = None
tadds = None

#grid_fields = ["pressure", "grid_level", "entropy"]
grid_fields = []

def rescan(dir, printlist=False):
    files = util.scan_files(dir, regex=regex, walk=True, reverse=True)
    return files


if MPI_taskpull2.rank == 0:
    partfiles = rescan(dir)

    # Arbitrarily pick particles in this file
    tp = h5py.File(partfiles[rand_findex].fullpath, 'r')['tracer particles']

    # Pick nparticles
    rints = sorted([random.randint(0, tp.shape[0]) for i in range(nparticles)])
    tags = tp[rints,15]
    tadds = tp[rints,8]
    shoks = tp[rints,7]

tags, tadds = MPI_taskpull2.comm.bcast((tags, tadds), root=0)

def init():
    for (tag, tadd, shok, ind) in zip(tags, tadds, shoks, rints):
        print(tag, tadd, shok, ind)


def worker_fn(filepath):
    parts = Particles(tags, tadds, grid_fields=grid_fields)
    parts.read_from_partfile(filepath)

    return parts

def tasks_gen():
    files = rescan(dir)
    for file in files[:]:
        yield file.fullpath

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

    picklename = "Particles_tracing.pickle"
    pickle.dump(parts, open( picklename, "wb" ))


