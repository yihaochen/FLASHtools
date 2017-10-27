#!/usr/bin/env python
import yt
import numpy as np
import matplotlib
matplotlib.rcParams['savefig.dpi'] = 150
import matplotlib.pyplot as plt
import os
import sys
sys.path.append('/home/ychen/lib/util')
import pickle
import time
import glob

import util

yt.enable_parallelism()

dirnames = [
    '/home/ychen/data/0only_0529_h1',
    '/home/ychen/data/0only_1022_h1_10Myr',
    '/home/ychen/data/0only_0605_hinf',
    '/home/ychen/data/0only_0204_hinf_10Myr',
    '/home/ychen/data/0only_0605_h0',
    '/home/ychen/data/0only_0204_h0_10Myr',
    ]


def get_Bs(fn):
    #print fn
    try:
        ds = yt.load(fn)
        ad = ds.all_data()
        Btor = sum(ad['magnetic_field_toroidal']**2*ad['cell_volume'])/(8.0*np.pi)
        Bpol = sum(ad['magnetic_field_poloidal']**2*ad['cell_volume'])/(8.0*np.pi)
        Emag = sum(ad['magnetic_pressure']*ad['cell_volume'])
        t = ds.current_time
        return t.in_units('s'), Btor.in_units('erg'), \
               Bpol.in_units('erg'), Emag.in_units('erg')
    except:
        print 'Error!'
        raise Exception

def get_Bflux(ds):
    ad = ds.all_data()
    Btor_flux = np.sum(ad['magnetic_field_toroidal']*(ad['dx']*ad['dy']))
    Bx_flux = np.sum(ad['magnetic_field_x']*ad['dy']*ad['dz'])
    By_flux = np.sum(ad['magnetic_field_y']*ad['dx']*ad['dz'])
    Bz_flux = np.sum(ad['magnetic_field_z']*ad['dx']*ad['dy'])
    t = ds.current_time
    return t.in_units('s'), Btor_flux.in_units('gauss*cm**2'),\
           Bx_flux.in_units('gauss*cm**2'), By_flux.in_units('gauss*cm**2'), Bz_flux.in_units('gauss*cm**2')
    #except:
    #    print 'Error!'
    #    raise Exception

def get_mags(dirname, parallel=1):
    #files = rescan(dirname, True)
    #Btors, Bpols, Emags = [], [], []
    #fnames = [f.fullpath for f in files[1:-1:5]]
    #print [fn.split('_')[-1] for fn in fnames]
    #if parallel:
    #    pool = multiprocessing.Pool(8)
    #print fnames
    #    result = pool.map(func, fnames)
    #else:
    #    result = map(func, fnames)
    storage = {}

    fnames = os.path.join(dirname,'*_hdf5_plt_cnt_???0')
    ts = yt.DatasetSeries(fnames, parallel=parallel)
    for sto, ds in ts.piter(storage=storage):
        sto.result = get_Bflux(ds)
        sto.result_id = str(ds)

    return storage



#files = glob.glob("Bflux_table_*.pickle")
#if files:
#    with open(files[-1], 'r') as f:
#        results = pickle.load(f)
#else:
#    results = {}
results = {}

t0 = time.time()
for dirname in dirnames:
    if dirname in results:
        next
    else:
        results[dirname] = get_mags(dirname, parallel=32)
        print 'Elapsed time:', time.time() - t0
        t0 = time.time()

picklename = time.strftime("Bflux_table_%Y%m%d_%H%M%S.pickle")
pickle.dump(results, open( picklename, "wb" ))

