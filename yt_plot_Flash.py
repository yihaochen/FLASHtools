#!/usr/bin/env python
import time
t0 = time.time()
import yt
import os
import sys
sys.path.append('/home/ychen/lib/util')
import util
#import matplotlib
import logging
logging.getLogger().setLevel(logging.ERROR)
import multiprocessing
from plotSlices import plotSlices

from tools import read_par, calcNozzleCoords

# Scan for files
#dir = '/d/d9/ychen/FLASH4/MHD_Jet_3D/1104_move_Heat'
dir = './'
pars = read_par(dir)

fields = ['density', 'temperature', 'pressure', 'velocity_z']
fields_grid = ['density', 'velocity_z']

def rescan(printlist=False):
    files = util.scan_files(dir, '*hdf5_plt_cnt_[0-9][0-9][0-9][0-9]', printlist=printlist)
    return files

def plotfile(file):
    try:
        ds = yt.load(file.fullpath)
        print "Plotting %s" % file.fullpath
        plotSlices(ds, zoom_fac=1, center="c", proj_axes=['x'], drawnozzle=True,\
                   plotgrid=fields_grid, savepath=os.path.join(dir, 'figures'))
    except KeyboardInterrupt:
        print 'KeyboardInterrupt catched...'
        raise Exception

files = rescan(True)
if not os.path.exists(os.path.join(files[-1].pathname,'figures/')):
    os.mkdir(os.path.join(files[-1].pathname,'figures/'))

t1 = time.time()

for file in files[:]:
    plotfile(file)

#nProcessor = multiprocessing.cpu_count()
#pool = multiprocessing.Pool(nProcessor)
#
#
#pool.imap(plotfile, reversed(files[:]))
#pool.close()
#pool.join()

t2 = time.time()
print 'Total time: %.2f s\n--\ninitialization: %.2f s\nploting: %.2f'\
       % (t2-t0, t1-t0, t2-t1)
