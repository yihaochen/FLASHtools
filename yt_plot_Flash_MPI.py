#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'stixgeneral'
import yt
import os
import sys
import util
import MPI_taskpull2
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
yt.mylog.setLevel('ERROR')

from plotSlices import plotSliceField, _entropy_ratio
from plotProjections import plotProjectionField
from particles.particle_filters import *
from synchrotron.yt_synchrotron_emissivity import setup_part_file

#dirs = ['/home/ychen/data/0only_1022_h1_10Myr']
#dirs = ['/home/ychen/data/0only_0204_hinf_10Myr',\
#        '/home/ychen/data/0only_0204_h0_10Myr']
dirs = ['./']
regex = '*_hdf5_plt_cnt_????'
#regex = 'MHD_Jet*_hdf5_plt_cnt_[0-9][0-9][0-9][0-9]'
files = None
zoom_facs = [2]
proj_axes= ['y']
figuredirtemplate = 'figures%s_zoom%i'
ptypes = ['lobe']


#annotate_particles = True if zoom_fac >= 2 else False
fields_velocity = None
fields_part = ['velocity_x', 'velocity_z']
fields = ['density', 'pressure', 'temperature', 'velocity_x', 'velocity_z', 'jet ',\
          'magnetic_field_y', 'magnetic_field_z', 'magnetic_pressure']
#fields += ['temperature_ratio', 'entropy_ratio']
#fields = ['density', 'shks']

def rescan(dir, printlist=False):
    files = util.scan_files(dir, regex=regex, walk=True, printlist=printlist, reverse=True)
    return files

e  = yt.utilities.physical_constants.elementary_charge #4.803E-10 esu
me = yt.utilities.physical_constants.mass_electron #9.109E-28
c  = yt.utilities.physical_constants.speed_of_light #2.998E10

def worker_fn(file, field, proj_axis, zoom_fac, ptype):
    ds = yt.load(file.fullpath)
    setup_part_file(ds)
    ds.add_particle_filter(ptype)
    #ds.periodicity = (True, True, True)

    #nozzleCoords = calcNozzleCoords(ds, proj_axis)
    #nozzleCoords = None if proj_axis=='z' or zoom_fac<4 else nozzleCoords
    nozzleCoords = None
    center = (0.0,0.0,7.714E22) if proj_axis=='z' else (0.0,0.0,0.0)
    fields_grid = ['density', 'pressure', 'velocity_y', 'velocity_z'] if zoom_fac >=16 \
             else ['velocity_z', 'density', 'shks']
    axis = '_%s' % proj_axis if proj_axis != 'x' else ''
    figuredir = figuredirtemplate % (axis, zoom_fac)
    if file.pathname.split('/')[-1] == 'data':
        absfiguredir = os.path.join(os.path.dirname(file.pathname), figuredir)
    else:
        absfiguredir = os.path.join(file.pathname, figuredir)

    #plotProjectionField(ds, zoom_fac=zoom_fac, center=center, proj_axis=proj_axis, field=field,\
    #               plotgrid=fields_grid, plotvelocity=fields_velocity, nozzleCoords=nozzleCoords, \
    #               annotate_particles=fields_part,annotate_part_info=False,\
    #               savepath=os.path.join(dir,figuredir,field))
    if file.pathname.split('/')[-1] == 'data':
        dirnamesplit = os.path.abspath(os.path.dirname(file.pathname)).split('_')
    else:
        dirnamesplit = os.path.abspath(file.pathname).split('_')

    if dirnamesplit[-1] in ['h1','hinf', 'h0'] and dirnamesplit[-2] in ['b1']:
        sim_name = dirnamesplit[-1]
    else:
        sim_name = dirnamesplit[-2] + '_' + dirnamesplit[-1]

    plotSliceField(ds, zoom_fac=zoom_fac, center=center, proj_axis=proj_axis, field=field,\
               plotgrid=fields_grid, plotvelocity=fields_velocity, \
               annotate_particles=fields_part,annotate_part_info=False, sim_name=sim_name,\
               savepath=os.path.join(absfiguredir,field.strip()))

    return sim_name, ds.basename[-4: ], field


def worker_test(file, field, proj_axis, zoom_fac):
    ds = yt.load(file.fullpath)
    return file.filename, field, proj_axis, zoom_fac


def tasks_gen(dirs, i0, fields, proj_axes, zoom_facs, ptypes):
    #sys.stdout.write("Plotting %s... \n" % (file.filename))
    for dir in dirs:
        files = rescan(dir, False)[i0:]
        for ptype in ptypes:
            for zoom_fac in zoom_facs:
                for proj_axis in proj_axes:
                    for file in files:
                        for field in fields:
                            yield file, field, proj_axis, zoom_fac, ptype


def init():
    for dir in dirs:
        for ptype in ptypes:
            for zoom_fac in zoom_facs:
                for proj_axis in proj_axes:
                    axis = '_%s' % proj_axis if proj_axis != 'x' else ''
                    figuredir = figuredirtemplate % (axis, zoom_fac)
                    if not os.path.exists(os.path.join(dir,figuredir)):
                        os.mkdir(os.path.join(dir, figuredir))
                    for field in fields:
                        if 'particle' in field:
                            fielddir = ptype+'_'+field.strip()
                        else:
                            fielddir = field.strip()
                        if not os.path.exists(os.path.join(dir,figuredir,fielddir)):
                            os.mkdir(os.path.join(dir,figuredir,fielddir))

i0 = int(sys.argv[1]) if len(sys.argv) > 1 else 0
tasks = tasks_gen(dirs, i0, fields, proj_axes, zoom_facs, ptypes)

results = MPI_taskpull2.taskpull(worker_fn, tasks, initialize=init, print_result=True)
