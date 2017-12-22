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
dirs = ['/home/ychen/data/0only_0330_h0_nojiggle',\
        '/home/ychen/data/0only_0518_hydro_nojiggle',\
        '/home/ychen/data/0only_0314_h1_nojiggle',\
        '/home/ychen/data/0only_0525_hinf_nojiggle']
#dirs = ['/home/ychen/data/0only_1106_M3_h1',\
#        '/home/ychen/data/0only_1204_M24_b01',\
#        '/home/ychen/data/0only_1110_h0_rerun']
#dirs = ['/home/ychen/data/0only_0529_h1/',\
#        '/home/ychen/data/0only_0605_h0/',\
#        '/home/ychen/data/0only_0605_hinf/',\
#        '/home/ychen/data/0only_0602_hydro/']

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
    kpc = yt.units.kpc.in_units('cm')
    #zlim = (ds.domain_left_edge[2]/factor, ds.domain_right_edge[2]/factor)
    #rlim = (0.0, ds.domain_right_edge[0]/factor)

    dx = ds.index.get_smallest_dx().in_units('cm')
    # Assuming the jet is propagating at velocity slower than 15 kpc/Myr
    zlim = max((15*ds.current_time.in_units('Myr').v*kpc)//dx, 100)*dx
    zlims = (-zlim, zlim)
    n_bins = int(2*zlim//dx)
    rlim = (0.0, ds.domain_right_edge[0]/8.0)

    ad = ds.all_data()
    prof = yt.create_profile(ad, 'z', ['jet '], logs={'z': False}, n_bins=n_bins, \
                             weight_field='cell_mass', extrema={'z': zlims})
    zmax = 0
    zmin = 0
    for (z, jet) in zip(reversed(prof.x), reversed(prof['jet '])):
        if jet > 1E-6:
            zmax = z
            break
    for (z, jet) in zip(prof.x, prof['jet ']):
        if jet > 1E-6:
            zmin = z
            break



    kpc = yt.units.kpc.in_units('cm')

    # Upper half of the simulation domain
    leftedge = [-30*kpc, -30*kpc, 0.0]
    rightedge = [ 30*kpc, 30*kpc, zmax.in_units('cm')*1.2]
    upbox = ds.box(leftedge, rightedge)
    upcenter = np.sum(upbox['jet ']*upbox['cell_mass']*upbox['z'])/np.sum(upbox['jet ']*upbox['cell_mass'])

    rmax_up = 0
    rprof = yt.create_profile(upbox, 'cylindrical_r', ['jet '], logs={'cylindrical_r': False}, n_bins=1024, \
                             weight_field='cell_mass', extrema={'cylindrical_r': rlim})
    for (r, jet) in zip(reversed(rprof.x), reversed(rprof['jet '])):
        if jet > 1E-6:
            rmax_up = r
            break

    # Lower half of the simulation domain
    leftedge = [-30*kpc, -30*kpc, zmin.in_units('cm')*1.2]
    rightedge = [ 30*kpc, 30*kpc, 0.0]
    bottombox = ds.box(leftedge, rightedge)
    bottomcenter = np.sum(bottombox['jet ']*bottombox['cell_mass']*bottombox['z'])/np.sum(bottombox['jet ']*bottombox['cell_mass'])

    rmax_low = 0
    rprof = yt.create_profile(bottombox, 'cylindrical_r', ['jet '], logs={'cylindrical_r': False}, n_bins=1024, \
                             weight_field='cell_mass', extrema={'cylindrical_r': rlim})
    for (r, jet) in zip(reversed(rprof.x), reversed(rprof['jet '])):
        if jet > 1E-6:
            rmax_low = r
            break

    return int(ds.basename[-4:]), ds.current_time.in_units('Myr'), zmax.in_units('kpc'), zmin.in_units('kpc'), \
           upcenter.in_units('kpc'), bottomcenter.in_units('kpc'), rmax_up.in_units('kpc'), rmax_low.in_units('kpc')

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
    #print collected

    #picklename = time.strftime("Bflux_table_%Y%m%d_%H%M%S.pickle")
    #pickle.dump(collected, open( picklename, "wb" ))

    fmt = '%04d  %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f'
    header = 'filenumber, t(Myr), zEdge+(kpc), zEdge-(kpc), zCenter+(kpc), zCenter-(kpc), rEdge+(kpc), rEdge-(kpc)'
    for dirname in collected.keys():
        #print np.asarray(collected[dirname])
        np.savetxt(dirname+'/GridAnalysis_LobeSize.txt', np.asarray(collected[dirname]), fmt=fmt, header=header)

#if MPI_taskpull2.rank == 0:
#    for key, item in results.items():
#        print item
