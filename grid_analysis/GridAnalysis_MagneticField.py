#!/usr/bin/env python
import yt
import numpy as np
import os.path as op
import sys
import util
import MPI_taskpull2
yt.mylog.setLevel('ERROR')
from yt.units import dimensions
from yt.utilities.physical_constants import mu_0
from yt.utilities.math_utils import \
    get_cyl_r_component, \
    get_cyl_theta_component, \
    get_cyl_z_component
from synchrotron.yt_synchrotron_emissivity import _jet_volume_fraction

kpc = yt.units.kpc.in_units('cm')

dirs = [
    '/home/ychen/data/0only_0529_h1',
#    '/home/ychen/data/0only_1022_h1_10Myr',
    '/home/ychen/data/0only_0605_hinf',
    '/home/ychen/data/0only_0204_hinf_10Myr',
    '/home/ychen/data/0only_1110_h0_rerun',
    '/home/ychen/data/0only_1212_h0_10Myr_rerun',
#    '/home/ychen/data/0only_0330_h0_nojiggle',\
#    '/home/ychen/data/0only_0314_h1_nojiggle',\
#    '/home/ychen/data/0only_0525_hinf_nojiggle'
    ]

regex = 'MHD_Jet*_hdf5_plt_cnt_[0-9][0-9][0-9][0-9]'

R = 150*kpc

mag_factors = {dimensions.magnetic_field_cgs: 4.0*np.pi,
               dimensions.magnetic_field_mks: mu_0}
ftype = 'gas'

def _magnetic_toroidal_energy(field, data):
    normal = data.get_field_parameter("normal")
    d = data[ftype,'magnetic_field_x']
    Bfields = data.ds.arr(
                [data[ftype,'magnetic_field_x'],
                 data[ftype,'magnetic_field_y'],
                 data[ftype,'magnetic_field_z']],
                 d.units)

    theta = data["index", 'cylindrical_theta']
    B = get_cyl_theta_component(Bfields, theta, normal)
    return 0.5*B*B/mag_factors[B.units.dimensions]*data['cell_volume']
yt.add_field((ftype, "magnetic_toroidal_energy"), sampling_type="cell",
         function=_magnetic_toroidal_energy,
         units='erg')

def _magnetic_r_energy(field, data):
    normal = data.get_field_parameter("normal")
    d = data[ftype,'magnetic_field_x']
    Bfields = data.ds.arr(
                [data[ftype,'magnetic_field_x'],
                 data[ftype,'magnetic_field_y'],
                 data[ftype,'magnetic_field_z']],
                 d.units)

    theta = data["index", 'cylindrical_theta']
    B = get_cyl_r_component(Bfields, theta, normal)
    return 0.5*B*B/mag_factors[B.units.dimensions]*data['cell_volume']
yt.add_field((ftype, "magnetic_r_energy"), sampling_type="cell",
         function=_magnetic_r_energy,
         units='erg')

def _magnetic_z_energy(field, data):
    B = data[ftype,"magnetic_field_z"]
    return 0.5*B*B/mag_factors[B.units.dimensions]*data['cell_volume']
yt.add_field((ftype, "magnetic_z_energy"), sampling_type="cell",
         function=_magnetic_z_energy,
         units='erg')

def _magnetic_total_energy(field, data):
    B = data[ftype,"magnetic_field_strength"]
    return 0.5*B*B/mag_factors[B.units.dimensions]*data['cell_volume']
yt.add_field((ftype, "magnetic_total_energy"), sampling_type="cell",
         function=_magnetic_total_energy,
         units='erg')

yt.add_field(('gas', 'jet_volume_fraction'), function=_jet_volume_fraction,
         display_name="Jet Volume Fraction", sampling_type='cell')

def _jet_volume(field, data):
    return data[ftype, 'jet_volume_fraction']*data['index', 'cell_volume']
yt.add_field(('gas', 'jet_volume'), function=_jet_volume,
         display_name="Jet Volume", sampling_type='cell',
         units='cm**3')

def _jet_volume_thres(field, data):
    jet_thres = np.where(data['flash', 'jet '] > 1E-3, 1.0, 0.0)
    return jet_thres*data['index', 'cell_volume']
yt.add_field(('gas', 'jet_volume_thres'), function=_jet_volume_thres,
         display_name="Jet Volume (Threshold)", sampling_type='cell',
         units='cm**3')

def _jet_mass(field, data):
    return data['flash', 'jet ']*data['gas', 'cell_mass']
yt.add_field(('gas', 'jet_mass'), function=_jet_mass,
         display_name="Jet Mass", sampling_type='cell',
         units='g')

def rescan(dir, printlist=False):
    files = util.scan_files(dir, regex=regex, walk=True, reverse=False)
    return files

def worker_fn(dirname, filepath):
    ds = yt.load(filepath)
    sp = ds.sphere([0,0,0], R)
    EBtor = sp.quantities.total_quantity('magnetic_toroidal_energy')
    EBr = sp.quantities.total_quantity('magnetic_r_energy')
    EBz = sp.quantities.total_quantity('magnetic_z_energy')
    Emag = sp.quantities.total_quantity('magnetic_total_energy')
    jet_volume = sp.quantities.total_quantity('jet_volume')
    jet_volume_thres = sp.quantities.total_quantity('jet_volume_thres')
    jet_mass = sp.quantities.total_quantity('jet_mass')

    return int(ds.basename[-4:]), ds.current_time.in_units('Myr'),\
            EBtor.in_units('erg'), EBr.in_units('erg'),\
            EBz.in_units('erg'), Emag.in_units('erg'),\
            jet_volume.in_units('kpc**3'), jet_volume_thres.in_units('kpc**3'),\
            jet_mass.in_units('Msun')

def tasks_gen(dirs):
    for dir in dirs:
        files = rescan(dir)
        for file in files[:]:
        #for file in reversed(files[:]):
            yield file.pathname, file.fullpath

tasks = tasks_gen(dirs)

results = MPI_taskpull2.taskpull(worker_fn, tasks, print_result=True)

if results:
    collected = {}
    for key, item in results.items():
        dirname, fname = key
        if dirname.split('/')[-1] == 'data':
            dirname = op.dirname(dirname) + '/'
        if dirname in collected:
            collected[dirname].append(item)
        else:
            collected[dirname] = [item]
    for key, item in collected.items():
        collected[key] = sorted(item)

    fmt = '%04d %6.3f %e %e %e %e %e %e %e'
    header = 'filenumber, t(Myr), EBtor(erg), EBr(erg), EBz(erg), Emag(erg), V_jet(kpc^3), V_jet_thres(kpc^3), Mjet(Msun)'
    for dirname in collected.keys():
        #print np.asarray(collected[dirname])
        np.savetxt(dirname+'/GridAnalysis_MagneticFields_%ikpc.txt' % R.in_units('kpc').v, np.asarray(collected[dirname]), fmt=fmt, header=header)
