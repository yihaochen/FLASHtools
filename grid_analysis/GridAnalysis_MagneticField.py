#!/usr/bin/env python
import yt
import numpy as np
import os.path as op
import sys
import util
import MPI_taskpull2
yt.mylog.setLevel('ERROR')


dirs = [
    '/home/ychen/data/0only_0529_h1',
    '/home/ychen/data/0only_1022_h1_10Myr',
    '/home/ychen/data/0only_0605_hinf',
    '/home/ychen/data/0only_0204_hinf_10Myr',
    '/home/ychen/data/0only_0605_h0',
    '/home/ychen/data/0only_0204_h0_10Myr',
    ]

regex = 'MHD_Jet*_hdf5_plt_cnt_[0-9][0-9][0-9][0-9]'

def rescan(dir, printlist=False):
    files = util.scan_files(dir, regex=regex, walk=True, reverse=False)
    return files

def worker_fn(dirname, filepath):
    ds = yt.load(filepath)
    kpc = yt.units.kpc.in_units('cm')
    leftedge =  [-100*kpc, -100*kpc, -100*kpc]
    rightedge = [ 100*kpc,  100*kpc,  100*kpc]
    box = ds.box(leftedge, rightedge)
    #EBtor = sum(box['magnetic_field_toroidal']**2*box['cell_volume'])/(8.0*np.pi)
    #EBpol = sum(box['magnetic_field_poloidal']**2*box['cell_volume'])/(8.0*np.pi)
    Emag  = sum(box['magnetic_pressure']*box['cell_volume'])

    return int(ds.basename[-4:]), ds.current_time.in_units('Myr'),\
            Emag.in_units('erg')
#            EBtor.in_units('erg'), EBpol.in_units('erg'), 

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
        if 'data' in dirname:
            dirname = op.dirname(dirname) + '/'
        if dirname in collected:
            collected[dirname].append(item)
        else:
            collected[dirname] = [item]
    for key, item in collected.items():
        collected[key] = sorted(item)
    #print collected

    #picklename = time.strftime("Bflux_table_%Y%m%d_%H%M%S.pickle")
    #pickle.dump(collected, open( picklename, "wb" ))

    fmt = '%04d %6.3f %e'# %e %e'
    #header = 'filenumber, t(Myr), EBtor(erg), EBpol(erg), Emag(erg)'
    header = 'filenumber, t(Myr), Emag(erg)'
    for dirname in collected.keys():
        #print np.asarray(collected[dirname])
        np.savetxt(dirname+'/GridAnalysis_MagneticFields.txt', np.asarray(collected[dirname]), fmt=fmt, header=header)
