#!/usr/bin/env python
import time
t0 = time.time()
import yt
yt.enable_parallelism()
import pickle
import numpy as np
import logging
from plotSlices import plotSlices
from yt.utilities.parallel_tools.parallel_analysis_interface import\
        parallel_blocking_call

Eint_threshold = 2.0E14

t1 = time.time()

toplot = None
logging.getLogger().setLevel(logging.ERROR)

@parallel_blocking_call
def initialize():
    try:
        array = pickle.load(open('minarray.p', 'r'))
        if yt.is_root():
            print 'minarray.p loaded'
    except:
        if yt.is_root():
            print 'minarray.p not found. Run yt_get_mins.py first'
        return

    flag = array['Eintmin'] < Eint_threshold
    global toplot
    toplot = array['basename'][flag]


def yt_plot_mins():
    initialize()
    tseries_toplot = yt.DatasetSeries(toplot, parallel=True)
    logging.getLogger().setLevel(logging.INFO)
    for ds in tseries_toplot.piter():
        min, loc_min = ds.h.find_min('eint')
        print 'Plotting %s' % ds.basename
        plotSlices(ds, zoom_fac=4, center=loc_min, drawnozzle=False,\
                   markcenter=True, proj_axes=['x', 'y', 'z'],\
                   fields=['eint', 'pressure'])
    if yt.is_root():
        t2 = time.time()
        print 'Total time: %.2f s\n--\ninitialization: %.2f s\nploting: %.2f'\
               % (t2-t0, t1-t0, t2-t1)

if __name__ == '__main__':
    yt_plot_mins()
