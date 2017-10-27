#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import os.path as op
import util
import yt
import MPI_taskpull2
import logging
logging.getLogger('yt').setLevel(logging.ERROR)

# Scan for files
dirs = ['/home/ychen/data/0only_0529_h1/',\
        '/home/ychen/data/0only_0605_h0/',\
        '/home/ychen/data/0only_0605_hinf/',\
        '/home/ychen/data/0only_0602_hydro/']

dirs = ['/home/ychen/data/0only_0602_hydro/']

#dirs = ['/home/ychen/data/0only_1022_h1_10Myr/',\
#        '/home/ychen/data/0only_0204_h0_10Myr/',\
#        '/home/ychen/data/0only_0204_hinf_10Myr/',\
#        '/home/ychen/data/0only_0413_hydro_10Myr/']

regex = 'MHD_Jet*_hdf5_plt_cnt_[0-9][0-9][0-9][0-9]'

def rescan(dir, printlist=False):
    files = util.scan_files(dir, regex=regex, walk=True, reverse=False)
    return files

def worker_fn(dirname, filepath):
    ds = yt.load(filepath)
    factor = 4 if '10Myr' in dirname else 8
    box = ds.box(ds.domain_left_edge/factor, ds.domain_right_edge/factor)
    enstrophy = sum(box['vorticity_squared']*box['cell_mass'])
    return int(ds.basename[-4:]), ds.current_time.in_units('Myr'), enstrophy

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
    print collected.keys()

    #picklename = time.strftime("Bflux_table_%Y%m%d_%H%M%S.pickle")
    #pickle.dump(collected, open( picklename, "wb" ))

    fmt = '%04d  %6.3f %8.3e'
    header = 'file t(Myr) enstrophy'
    for dirname in collected.keys():
        #print np.asarray(collected[dirname])
        np.savetxt(dirname+'/GridAnalysis_Enstrophy.txt', np.asarray(collected[dirname]), fmt=fmt, header=header)

#if MPI_taskpull2.rank == 0:
#    for key, item in results.items():
#        print item
