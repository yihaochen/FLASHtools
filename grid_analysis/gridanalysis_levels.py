#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import os.path as op
import h5py
import util
import yt
import MPI_taskpull2
yt.mylog.setLevel("ERROR")

# Scan for files
dirs = [
        '/d/d9/ychen/2018_production_runs/20180801_L430_rc10_beta07/data',\
        '/d/d9/ychen/2018_production_runs/20180802_L438_rc10_beta07/data',\
        '/d/d9/ychen/2018_production_runs/20180803_L446_rc10_beta07/data',\
        '/d/d9/ychen/2018_production_runs/20180824_L438_rc30_beta07/data',\
        '/d/d9/ychen/2018_production_runs/20180906_L430_rc30_beta07/data',\
        '/d/d9/ychen/2018_production_runs/20180907_L446_rc30_beta07/data'\
        ]

regex = '*_hdf5_plt_cnt_?[0,1]00'

maxlevel = 12

def rescan(dir, printlist=False):
    files = util.scan_files(dir, regex=regex, walk=True, reverse=False)
    return files

def worker_fn(dirname, filepath):
    h5f = h5py.File(filepath)
    current_time = h5f['real scalars'][0][1]/3.15576E13
    leaf_nodes = h5f['node type'][:] == 1
    zcoords = h5f['coordinates'][leaf_nodes,2]
    levels = h5f['refine level'][leaf_nodes]
    numgrids = [np.sum(np.logical_and(zcoords>0, levels==level)) for level in range(1,maxlevel+1)]+\
               [np.sum(np.logical_and(zcoords<0, levels==level)) for level in range(1,maxlevel+1)]

    return (int(filepath[-4:]), current_time, numgrids)

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

    fmt = '%04d %6.3f' + ' %6i'*len(maxlevel)*2
    header = 'filenumber, t(Myr), ' + 'upper domain ' + 'lev%2i,'*12 % tuple(range(1,maxlevel+1))+\
                                      'lower domain ' + 'lev%2i,'*12 % tuple(range(1,maxlevel+1))
    for dirname in collected.keys():
        #print np.asarray(collected[dirname])
        if dirname.strip('/').split('/')[-1] == 'data':
            savedir = op.dirname(dirname)
        else:
            savedir = dirname
        np.savetxt(savedir+'/gridanalysis_levels.txt', np.asarray(collected[dirname]), fmt=fmt, header=header)

#if MPI_taskpull2.rank == 0:
#    for key, item in results.items():
#        print item
