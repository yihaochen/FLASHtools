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
dirs = ['/d/d9/ychen/2018_production_runs/20180801_L430_rc10_beta07/',\
        '/d/d9/ychen/2018_production_runs/20180802_L438_rc10_beta07/',\
        '/d/d9/ychen/2018_production_runs/20180803_L446_rc10_beta07/',\
        '/d/d9/ychen/2018_production_runs/20180824_L438_rc30_beta07/']

regex = '*_hdf5_plt_cnt_[0-9][0-9][0-9][0-9]'

zoom_fac = 4

def rescan(dir, printlist=False):
    files = util.scan_files(dir, regex=regex, walk=True, reverse=False)
    return files

def worker_fn(dirname, filepath):
    ds = yt.load(filepath)
    kpc = yt.units.kpc.in_units('cm')
    #zlim = (ds.domain_left_edge[2]/factor, ds.domain_right_edge[2]/factor)
    #rlim = (0.0, ds.domain_right_edge[0]/factor)

    dx = ds.index.get_smallest_dx().v
    # Assuming the jet is propagating at velocity slower than 15 kpc/Myr
    #zlim = max((15*ds.current_time.in_units('Myr').v*kpc)//dx, 100)*dx
    zlim = ds.domain_width[2].v/zoom_fac/2.0
    zlims = (-zlim, zlim)
    n_bins = int(2*zlim//dx)
    rlim = (0.0, ds.domain_width[0]/zoom_fac/2.0)

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

    xlim = ds.domain_width[0]/zoom_fac/2.0
    ylim = ds.domain_width[1]/zoom_fac/2.0
    # Upper half of the simulation domain
    leftedge =  [-xlim, -ylim, 0.0]
    rightedge = [ xlim,  ylim, zmax.in_units('cm')*1.2]
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

collected = {}

def tasks_gen(dirs):
    analyzed_fnumbers = {}
    for dir in dirs:
        try:
            table = np.loadtxt(dir+'/GridAnalysis_LobeSize.txt')
            collected[dir] = list(table)
            analyzed_fnumbers[dir] = np.array(table[:,0], dtype=int)
        except:
            analyzed_fnumbers[dir] = []
    for dir in dirs:
        files = rescan(dir)
        for file in reversed(files[:]):
            if int(file.filename[-4:]) not in analyzed_fnumbers[dir]:
                yield file.pathname, file.fullpath

tasks = tasks_gen(dirs)

results = MPI_taskpull2.taskpull(worker_fn, tasks, print_result=True)

if results:
    for key, item in results.items():
        dirname, fname = key
        if dirname.split('/')[-1] == 'data':
            dirname = op.dirname(dirname) + '/'
        if dirname in collected:
            collected[dirname].append(np.array(item))
        else:
            collected[dirname] = [item]
    for key, item in collected.items():
        collected[key] = sorted(item, key=lambda arr: arr[0])
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
