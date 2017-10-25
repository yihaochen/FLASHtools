#!/usr/bin/env python
import yt
import os
import numpy as np
import glob
import time
import pickle
import MPI_taskpull2


def get_Bflux(dirname, fname):
    ds = yt.load(os.path.join(dirname,fname))
    t = ds.current_time
    #ad = ds.all_data()
    #ad = ds.sphere((0,0,0), (2, 'kpc'))
    ad = ds.disk([0,0,0], [0,0,1], (0.7, 'kpc'), (1, 'kpc'))
    Btor_flux = np.sum(ad['magnetic_field_toroidal']*(ad['dx']*ad['dy']))
    Bx_flux = np.sum(ad['magnetic_field_x']*ad['dy']*ad['dz'])
    By_flux = np.sum(ad['magnetic_field_y']*ad['dx']*ad['dz'])
    Bz_flux = np.sum(ad['magnetic_field_z']*ad['dx']*ad['dy'])
    #return t.in_units('s'), t.in_units('Gyr'), t.in_units('Myr'), t.in_units('kyr'), t.in_units('yr')
    return t.in_units('s'), Btor_flux.in_units('gauss*cm**2'),\
           Bx_flux.in_units('gauss*cm**2'), By_flux.in_units('gauss*cm**2'), Bz_flux.in_units('gauss*cm**2')

def worker(*args):
    import time
    time.sleep(0.5)
    return args[0]**2, args[1]**3

def tasks(n):
    for i in range(n):
        yield i, 2*i

def task_gen(dirnames):
    for dirname in dirnames:
        fnames = glob.glob(os.path.join(dirname, '*_hdf5_plt_cnt_???0'))
        print dirname, len(fnames)
        for fname in sorted(fnames, reverse=True):
            if '0000' in fname: next
            yield dirname, fname.split('/')[-1]

dirnames = [
    '/home/ychen/data/0only_1022_h1_10Myr',
    '/home/ychen/data/0only_0529_h1',
    '/home/ychen/data/0only_0204_hinf_10Myr',
    '/home/ychen/data/0only_0605_hinf',
    '/home/ychen/data/0only_0204_h0_10Myr',
    '/home/ychen/data/0only_0605_h0',
    ]

#try:
    #results = MPI_taskpull2.taskpull(worker, tasks(100))
#except KeyboardInterrupt:
#    print 'Interrupt Signal caught...'
#finally:

results = MPI_taskpull2.taskpull(get_Bflux, task_gen(dirnames))
if results:
    collected = {}
    for key, item in results.items():
        dirname, fname = key
        if dirname in collected:
            collected[dirname].append(item)
        else:
            collected[dirname] = [item]
    for key, item in collected.items():
        collected[key] = sorted(item)
    #print collected

    picklename = time.strftime("Bflux_table_%Y%m%d_%H%M%S.pickle")
    pickle.dump(collected, open( picklename, "wb" ))

    #results = sorted(results.items(), key=lambda (k,v):k)
    #print zip(*results)
    #for result in results:
    #    print result
    #for key, item in results.items():
    #    print key.split('/')[-2][6:], key.split('/')[-1], item
