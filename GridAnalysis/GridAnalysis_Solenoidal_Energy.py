#!/usr/bin/env python
import yt
import util
from solenoidal import solenoidal
import MPI_taskpull2
import logging
logging.getLogger('yt').setLevel(logging.WARNING)
import numpy as np


dirs = ['/d/d7/wattal/ACI/Bubble_beta_outsideMag_isolateBubble_ism_beta100']
regex = 'bubble_mhd_hdf5_plt_cnt_???0'

def rescan(dir, printlist=False):
    files = util.scan_files(dir, regex=regex, walk=True, reverse=False)
    return files

def worker_fn(dirname, filepath):
    ds = yt.load(filepath)
    #ad = ds.all_data()
    #xweight = (ad['x']*ad['jet ']*ad['cell_volume']).sum()/(ad['jet ']*ad['cell_volume']).sum()

    fields = ['velocity_x', 'velocity_y', 'velocity_z', 'density']
    max_level = 5
    low = (xweight-ds.quan(5.5, 'code_length'), -ds.quan(2., 'code_length'),-ds.quan(2., 'code_length'))
    dims = ds.arr([8.,4.,4.],'code_length')/(ds.domain_width/(ds.domain_dimensions*ds.refine_by**ds.max_level))
    cube = ds.covering_grid(ds.max_level, left_edge=low, dims=dims, fields=fields)
    vr = solenoidal(np.array([cube[field] for field in fields[:-1]]))
    kin_vr = (cube['cell_mass'].v*(vr**2).sum(axis=0)).sum()

    return int(ds.basename[-4:]), ds.current_time, kin_vr

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

    fmt = '%04d %8.4f %8.4fe'
    header = 'filenumber, t(s), E_vr(erg)'
    for dirname in collected.keys():
        #print np.asarray(collected[dirname])
        np.savetxt('./GridAnalysis_Solenoidal_Energy.txt', np.asarray(collected[dirname]), fmt=fmt, header=header)

