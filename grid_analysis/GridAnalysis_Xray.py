#!/usr/bin/env python
import yt
import numpy as np
import os.path as op
import sys
import util
import MPI_taskpull2
yt.mylog.setLevel('ERROR')
from synchrotron.yt_synchrotron_emissivity import _jet_volume_fraction

kpc = yt.units.kpc.in_units('cm')

dirs = [
    '/home/ychen/data/0only_0529_h1',
    '/home/ychen/data/0only_1022_h1_10Myr',
    '/home/ychen/data/0only_0605_hinf',
    '/home/ychen/data/0only_0602_hydro',
    '/home/ychen/data/0only_0413_hydro_10Myr',
    '/home/ychen/data/0only_0204_hinf_10Myr',
    '/home/ychen/data/0only_1110_h0_rerun',
    '/home/ychen/data/0only_1212_h0_10Myr_rerun',
    ]

regex = 'MHD_Jet*_hdf5_plt_cnt_[0-9][0-9][0-9][0-9]'

R = 100*kpc
metallicity = 0.5

def rescan(dir, printlist=False):
    files = util.scan_files(dir, regex=regex, walk=True, reverse=False)
    return files




def worker_fn(dirname, filepath):
    ds = yt.load(filepath)
    ds.add_field(('gas', 'jet_volume_fraction'), function=_jet_volume_fraction,
                 display_name="Jet Volume Fraction", sampling_type='cell')
    sp = ds.sphere([0,0,0], R)
    xray_fields = yt.add_xray_emissivity_field(ds, 0.1, 100, table_type='apec', metallicity=metallicity)
    emis, l, pemis = xray_fields
    def _xray_lum_icm(field, data):
        return data[l]*(1.0 - data['jet '])
    ds.add_field(('gas', l+'_icm'), function=_xray_lum_icm,
            units='erg/s',
            display_name='X-ray luminosity (icm)', sampling_type='cell')
    xray_lum = sp.quantities.total_quantity(l)

    xray_icm = sp.quantities.total_quantity(l+'_icm')

    return int(ds.basename[-4:]), ds.current_time.in_units('Myr'),\
            xray_lum.in_units('erg/s'), xray_icm.in_units('erg/s')
#            EBtor.in_units('erg'), EBpol.in_units('erg'), 

def tasks_gen(dirs):
    for dir in dirs:
        files = rescan(dir)
        for file in files[:]:
            yield file.pathname, file.fullpath

tasks = tasks_gen(dirs)

results = MPI_taskpull2.taskpull(worker_fn, tasks, print_result=True)

if results:
    collected = {}
    for key, item in results.items():
        dirname, fname = key
        if dirname.split('/')[-1] == 'data:
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

    fmt = '%04d %6.3f %e %e'# %e'
    #header = 'filenumber, t(Myr), EBtor(erg), EBpol(erg), Emag(erg)'
    header = 'filenumber, t(Myr), Lx(erg/s), Lx_icm(erg/s)'
    for dirname in collected.keys():
        #print np.asarray(collected[dirname])
        np.savetxt(dirname+'/GridAnalysis_Xray_%ikpc.txt' % R.in_units('kpc').v, np.asarray(collected[dirname]), fmt=fmt, header=header)
