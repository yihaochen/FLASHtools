#!/usr/bin/env python
import time
t0 = time.time()
import yt
import pickle
import numpy as np
import logging
logging.getLogger().setLevel(logging.ERROR)
import multiprocessing
from plotSlices import plotSlices

Eint_threshold = 2.0E14


try:
    array = pickle.load(open('minarray.p', 'r'))
    if yt.is_root():
        print 'minarray.p loaded'
except:
    if yt.is_root():
        print 'minarray.p not found. Run yt_get_mins.py first'
        raise IOError

flag = array['Eintmin'] < Eint_threshold
toplot = array['basename'][flag]

def plotfile(fn):
    try:
        ds = yt.load(fn)
        print 'Finding minimum location in %s' % fn
        min, loc_min = ds.h.find_min('eint')
        print 'Plotting %s' % ds.basename
        plotSlices(ds, zoom_fac=8, center=loc_min, drawnozzle=False,\
                   markcenter=True, proj_axes=['x', 'z'],\
                   fields=['eint', 'pressure'])
    except KeyboardInterrupt:
        print 'KeyboardInterrupt catched...'
        raise Exception


nProcessor = multiprocessing.cpu_count()
pool = multiprocessing.Pool(nProcessor-1)

t1 = time.time()

pool.imap(plotfile, reversed(toplot))
pool.close()
pool.join()


t2 = time.time()
print 'Total time: %.2f s\n--\ninitialization: %.2f s\nploting: %.2f'\
       % (t2-t0, t1-t0, t2-t1)

