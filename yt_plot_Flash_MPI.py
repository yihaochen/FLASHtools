#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import yt
import os
import sys
#sys.path.append(os.getenv('HOME') + '/lib/util')
import util
import MPI_taskpull2
import logging
logging.getLogger('yt').setLevel(logging.ERROR)
import matplotlib.pyplot as plt
import numpy as np

from tools import calcNozzleCoords
from plotSlices import plotSliceField
from plotProjections import plotProjectionField

dir = './'
regex = 'MHD_Jet*_hdf5_plt_cnt_[0-9][0-9][0-9][0-9]'
files = None
zoom_facs = [8]
proj_axes= ['x']


#annotate_particles = True if zoom_fac >= 2 else False
fields_velocity = ['velocity_z']
fields_part = None
fields = ['density', 'pressure', 'temperature', 'velocity_y', 'velocity_z', 'jet ',\
          'magnetic_field_x', 'magnetic_field_y', 'magnetic_field_z', 'magnetic_pressure',\
          'particle_gamc', 'plasma_beta']
#fields = ['density', 'pressure','shok', 'velocity_y', 'velocity_z', 'jet ']
#fields = ['density', 'pressure', 'temperature', 'velocity_magnitude', 'velocity_z']
#fields = ['magnetic_pressure']
#fields = ['density', 'temperature']
#fields = ['plasma_beta']
#fields = ['magnetic_field_x', 'magnetic_field_y', 'magnetic_field_z', 'magnetic_pressure',\
#          'velocity_y', 'velocity_z']


def rescan(printlist=False):
    files = util.scan_files(dir, regex=regex, walk=True, printlist=printlist, reverse=False)
    return files


def worker_fn(file, field, proj_axis, zoom_fac):
    particle_path = file.fullpath.replace('plt_cnt', 'part')\
                    if 'plt_cnt' in file.fullpath\
                    else file.fullpath.replace('chk', 'part')
    if os.path.exists(particle_path):
        ds = yt.load(file.fullpath, particle_filename=particle_path)
    else:
        ds = yt.load(file.fullpath)
        global fields_part
        fields_part = False

    nozzleCoords = calcNozzleCoords(ds, proj_axis)
    #nozzleCoords = None if proj_axis=='z' or zoom_fac<4 else nozzleCoords
    nozzleCoords = None
    center = (0.0,0.0,7.714E22) if proj_axis=='z' else (0.0,0.0,0.0)
    #center = (6.125E+20,  1.414E+20, -8.953E+20) if zoom_fac==128 else (0.0,0.0,0.0)
    fields_grid = ['density', 'pressure', 'velocity_y', 'velocity_z'] if zoom_fac >=16 \
             else ['velocity_z', 'temperature', 'density']
    #fields_grid = True
    figuredir = 'figures_%s' % proj_axis if proj_axis != 'x' else 'figures'
    figuredir = figuredir + '_zoom%i' % zoom_fac

    #plotProjectionField(ds, zoom_fac=zoom_fac, center=center, proj_axis=proj_axis, field=field,\
    #               plotgrid=fields_grid, plotvelocity=fields_velocity, nozzleCoords=nozzleCoords, \
    #               annotate_particles=fields_part,annotate_part_info=False,\
    #               savepath=os.path.join(dir,figuredir,field))
    if field=='particle_gamc':
        ad = ds.all_data()
        filter = np.logical_and((ad['all', 'particle_tag'] % 20 == 0), (ad['all', 'particle_shok']==0))
        y = ad['all', 'particle_position_y'][filter]/3.08567758E21
        z = ad['all', 'particle_position_z'][filter]/3.08567758E21
        gamc = np.log10(np.abs(ad['all', 'particle_gamc'][filter]))
        #plt.hist(gamc, bins=50)
        fig=plt.figure(figsize=(5,8))
        ax=fig.add_subplot(111)
        ax.set_xlim(-30,30)
        ax.set_ylim(-60,60)
        ax.set_xlabel('kpc')
        ax.set_ylabel('kpc')
        ax.annotate('%6.3f Myr' % (float(ds.current_time)/3.15569E13),\
                    (1,0), xytext=(0.05, 0.96),  textcoords='axes fraction',\
                    horizontalalignment='left', verticalalignment='center')
        sc=ax.scatter(y,z,s=1,c=gamc,linewidth=0,cmap='algae',vmin=2,vmax=5,alpha=0.8)
        try:
            cb=plt.colorbar(sc)
            cb.set_label(u'log $\gamma_c$')
        except:
            pass
        plt.savefig(os.path.join(dir,figuredir,field,file.filename), dpi=150)


    else:
        dirnamesplit = os.path.abspath(file.pathname).split('_')
        if dirnamesplit[-1] in ['h1','hinf', 'h0']:
            sim_name = dirnamesplit[-1]
        else:
            sim_name = dirnamesplit[-2] + '_' + dirnamesplit[-1]
        plotSliceField(ds, zoom_fac=zoom_fac, center=center, proj_axis=proj_axis, field=field,\
                   plotgrid=fields_grid, plotvelocity=fields_velocity, nozzleCoords=nozzleCoords, \
                   annotate_particles=fields_part,annotate_part_info=False, sim_name=sim_name,\
                   savepath=os.path.join(dir,figuredir,field.strip()))
    return ds.basename[-4: ], field


def worker_test(file, field, proj_axis, zoom_fac):
    ds = yt.load(file.fullpath)
    return file.filename, field, proj_axis, zoom_fac


def tasks_gen(i0, fields, proj_axes, zoom_facs):
    #sys.stdout.write("Plotting %s... \n" % (file.filename))
    files = rescan(False)[i0:]
    for zoom_fac in zoom_facs:
        for proj_axis in proj_axes:
            for file in files:
                for field in fields:
                    yield file, field, proj_axis, zoom_fac


def init():
    for zoom_fac in zoom_facs:
        for proj_axis in proj_axes:
            figuredir = 'figures_%s' % proj_axis if proj_axis != 'x' else 'figures'
            figuredir = figuredir + '_zoom%i' % zoom_fac
            if not os.path.exists(os.path.join(dir,figuredir)):
                os.mkdir(os.path.join(dir, figuredir))
            for field in fields:
                if not os.path.exists(os.path.join(dir,figuredir,field.strip())):
                    os.mkdir(os.path.join(dir,figuredir,field.strip()))

i0 = int(sys.argv[1]) if len(sys.argv) > 1 else 0
tasks = tasks_gen(i0, fields, proj_axes, zoom_facs)

results = MPI_taskpull2.taskpull(worker_fn, tasks, initialize=init, print_result=True)
#print results
