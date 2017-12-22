#!/usr/bin/env python
import numpy as np
import os.path as op
import util
import yt
import MPI_taskpull2
yt.mylog.setLevel('ERROR')

# Scan for files
#dirs = ['/home/ychen/data/0only_0529_h1/',\
#        '/home/ychen/data/0only_0605_h0/',\
#        '/home/ychen/data/0only_0605_hinf/',\
#        '/home/ychen/data/0only_0602_hydro/']

dirs = ['/home/ychen/data/0only_1022_h1_10Myr/',\
        '/home/ychen/data/0only_0204_h0_10Myr/',\
        '/home/ychen/data/0only_0204_hinf_10Myr/',\
        '/home/ychen/data/0only_0413_hydro_10Myr/']

regex = 'MHD_Jet*_hdf5_plt_cnt_[0-9][0-9][0-9][0-9]'

metallicity = 0.5
emin = 0.1
emax = 100.0

def rescan(dir, printlist=False):
    files = util.scan_files(dir, regex=regex, walk=True, reverse=False)
    return files

def worker_fn(dirname, filepath):
    ds = yt.load(filepath)
    xray_fields = yt.add_xray_emissivity_field(ds, emin, emax, table_type='apec', metallicity=metallicity)
    emis, l, pemis = xray_fields
    kpc = yt.units.kpc.in_units('cm')
    leftedge =  [-100*kpc, -100*kpc, -100*kpc]
    rightedge = [ 100*kpc,  100*kpc,  100*kpc]
    box = ds.box(leftedge, rightedge)
    xray_lum = box.quantities.total_quantity(l)

    return int(ds.basename[-4:]), ds.current_time.in_units('Myr'), xray_lum

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

    fmt = '%04d %6.3f %e'
    header = 'filenumber, t(Myr), Lx_%.1f_%.1f (erg/s)' % (emin, emax)
    for dirname in collected.keys():
        #print np.asarray(collected[dirname])
        np.savetxt(dirname+'/GridAnalysis_Xray.txt', np.asarray(collected[dirname]), fmt=fmt, header=header)

