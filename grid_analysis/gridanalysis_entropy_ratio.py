#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import os.path as op
import util
import yt
import MPI_taskpull2
from yt_cluster_ratio_fields import *
yt.mylog.setLevel("ERROR")

# Scan for files
dirs = ['/home/ychen/data/0only_1022_h1_10Myr/']

regex = 'MHD_Jet*_hdf5_plt_cnt_????'
#regex = 'MHD_Jet*_hdf5_plt_cnt_[0-9][0-9][0-9][0-9]'

ratios = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
#ratios = [0.6]

def rescan(dir, printlist=False):
    files = util.scan_files(dir, regex=regex, walk=True, reverse=False)
    return files

def worker_fn(dirname, filepath):
    ds = yt.load(filepath)

    center = [0.0,0.0,0.0]
    radius = (200, 'kpc')
    sp = ds.sphere(center, radius)

    masses = []
    for ratio in ratios:
        low_entropy = sp.cut_region(["obj['entropy_ratio'] < %.1f" % ratio])
        masses.append(sum(low_entropy['cell_mass'].in_units('Msun')))

    return (int(ds.basename[-4:]), ds.current_time.in_units('Myr')) + tuple(masses)

def tasks_gen(dirs):
    for dir in dirs:
        files = rescan(dir)
        for file in reversed(files[:]):
            yield file.pathname, file.fullpath

tasks = tasks_gen(dirs)

results = MPI_taskpull2.taskpull(worker_fn, tasks, print_result=True)

if results:
    collected = {}
    for key, item in results.items():
        dirname, fname = key
        if 'restart' in dirname:
            dirname = op.dirname(dirname) + '/'
        if dirname in collected:
            collected[dirname].append(item)
        else:
            collected[dirname] = [item]
    for key, item in collected.items():
        collected[key] = sorted(item)

    #picklename = time.strftime("Bflux_table_%Y%m%d_%H%M%S.pickle")
    #pickle.dump(collected, open( picklename, "wb" ))

    fmt = '%04d %6.3f' + ' %e'*len(ratios)
    header = 'filenumber, t(Myr), ' + 'mass(Msun)'*len(ratios)
    for dirname in collected.keys():
        #print np.asarray(collected[dirname])
        if dirname.strip('/').split('/')[-1] == 'data':
            savedir = op.dirname(dirname)
        else:
            savedir = dirname
        np.savetxt(dirname+'/gridanalysis_entropy_ratio.txt', np.asarray(collected[dirname]), fmt=fmt, header=header)

#if MPI_taskpull2.rank == 0:
#    for key, item in results.items():
#        print item
