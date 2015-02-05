#!/usr/bin/env python
import time
t0 = time.time()
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import yt
try:
    from mpi4py import MPI
except:
    MPI = None
yt.enable_parallelism()
import logging
logging.getLogger().setLevel(logging.ERROR)
import pickle

t1 = time.time()


# Scan for files
fregex = 'MHD_Jet_hdf5_plt_cnt_00*'
tseries = yt.DatasetSeries(fregex, parallel=True)
storage = {}

for sto, ds in tseries.piter(storage=storage):
    if MPI:
        print 'Loading %s ... (%3i/%3i)' % (ds.basename, MPI.COMM_WORLD.rank+1, MPI.COMM_WORLD.size)
    else:
        print 'Loading %s ...' % (ds.basename)
    alldata = ds.all_data()
    Eintmin, Eintmax = alldata.quantities.extrema('eint')
    Pmin, Pmax = alldata.quantities.extrema('pressure')
    Rhomin, Rhomax = alldata.quantities.extrema('density')

    sto.result = (ds.basename, ds.current_time, Eintmin, Eintmax, Pmin, Pmax, Rhomin, Rhomax)

t2 = time.time()
array = None

if yt.is_root():

    names = ['basename', 'time', 'Eintmin', 'Eintmax', 'Pmin', 'Pmax', 'Rhomin', 'Rhomax']
    formats = ['S25',      'f8',   'f8'     , 'f8'     , 'f8'  , 'f8'  , 'f8'    , 'f8']
    dtype = dict(names=names, formats=formats)
    array = np.array(sorted(storage.values()), dtype=dtype)


    fig, ax0 = plt.subplots()
    lines = [None]*3

    lines[0] = ax0.plot(array['time'], array['Eintmin'], '*-', color='blue', label='min eint')[0]
    #ax0.set_xlabel('t')
    ax0.set_ylabel('Eint')
    ax0.semilogy()

    ax1 = ax0.twinx()
    lines[1] = ax1.plot(array['time'], array['Pmin'], 'o-', color='red', label='min pressure')[0]
    #ax1.set_ylabel('pressure')
    ax1.semilogy()

    ax2 = ax0.twinx()
    lines[2] = ax2.plot(array['time'], array['Rhomin'], '.-', color='green', label='min density')[0]
    #ax2.set_ylabel('density')
    ax2.semilogy()

    labels = [l.get_label() for l in lines]
    leg = ax0.legend(lines, labels, loc=3)
    for text, line in zip(leg.get_texts(), leg.get_lines()):
        text.set_color(line.get_color())

    plt.title(os.getcwd().split('/')[-1])
    #plt.xlim(0, 5.0E13)
    plt.savefig('min_t.png')

    pickle.dump(array, open('minarray.p', 'wb'))

    t3 = time.time()
    print 'Total time: %.2f s\n--\ninitialization: %.2f s\nfinding exterma: %.2f s\nploting: %.2f'\
           % (t3-t0, t1-t0, t2-t1, t3-t2)



