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

from plotSlices import plotSliceField
from plotProjections import plotProjectionField
from particle_filters import *

#dirs = ['/home/ychen/data/0only_0314_h1_nojiggle',\
#        '/home/ychen/data/0only_0330_h0_nojiggle',\
#         '/home/ychen/data/0only_0525_hinf_nojiggle']
#dirs = ['/home/ychen/data/0only_0518_hydro_nojiggle']
#dirs = ['/home/ychen/data/0only_1022_h1_10Myr']
#dirs = ['/home/ychen/data/0only_0204_hinf_10Myr',\
#        '/home/ychen/data/0only_0204_h0_10Myr']
#dir = '/d/d9/ychen/FLASH4/stampede/0529_L45_M10_b1_h1/'
dirs = ['./']
#regex = 'MHD_Jet*_hdf5_plt_cnt_0[2-9][0,2,4,6,8]0'
regex = 'MHD_Jet*_hdf5_plt_cnt_[0-9][0-9][0-9][0-9]'
files = None
zoom_facs = [8]
proj_axes= ['x']
figuredirtemplate = 'figures%s_zoom%i'
ptype = 'jetp'


#annotate_particles = True if zoom_fac >= 2 else False
fields_velocity = None
fields_part = ['velocity_y']
fields = ['density', 'pressure', 'temperature', 'velocity_y', 'velocity_z', 'jet ',\
          'magnetic_field_x', 'magnetic_field_z', 'magnetic_pressure',\
          'plasma_beta', 'entropy', 'particle_gamc']
#fields = ['density', 'pressure', 'temperature', 'velocity_x', 'velocity_z', 'jet ',\
#          'entropy']
#fields = ['particle_nuc', 'particle_gamc', 'particle_age', 'particle_magp']
#fields = ['particle_gamc']
#fields = ['plasma_beta']


def rescan(dir, printlist=False):
    files = util.scan_files(dir, regex=regex, walk=False, printlist=printlist, reverse=True)
    return files

e  = yt.utilities.physical_constants.elementary_charge #4.803E-10 esu
me = yt.utilities.physical_constants.mass_electron #9.109E-28
c  = yt.utilities.physical_constants.speed_of_light #2.998E10

def worker_fn(file, field, proj_axis, zoom_fac):
    #particle_path = file.fullpath.replace('plt_cnt', 'part')\
    #                if 'plt_cnt' in file.fullpath\
    #                else file.fullpath.replace('chk', 'part')
    #if os.path.exists(particle_path):
    #    ds = yt.load(file.fullpath, particle_filename=particle_path)
    #else:
    #    ds = yt.load(file.fullpath)
    #    global fields_part
    #    fields_part = False
    ds = yt.load(file.fullpath)
    ds.add_particle_filter(ptype)
    ds.periodicity = (True, True, True)

    #nozzleCoords = calcNozzleCoords(ds, proj_axis)
    #nozzleCoords = None if proj_axis=='z' or zoom_fac<4 else nozzleCoords
    nozzleCoords = None
    center = (0.0,0.0,7.714E22) if proj_axis=='z' else (0.0,0.0,0.0)
    #center = (6.125E+20,  1.414E+20, -8.953E+20) if zoom_fac==128 else (0.0,0.0,0.0)
    fields_grid = ['density', 'pressure', 'velocity_y', 'velocity_z'] if zoom_fac >=16 \
             else ['velocity_z', 'density']
    #fields_grid = True
    axis = '_%s' % proj_axis if proj_axis != 'x' else ''
    figuredir = figuredirtemplate % (axis, zoom_fac)

    #plotProjectionField(ds, zoom_fac=zoom_fac, center=center, proj_axis=proj_axis, field=field,\
    #               plotgrid=fields_grid, plotvelocity=fields_velocity, nozzleCoords=nozzleCoords, \
    #               annotate_particles=fields_part,annotate_part_info=False,\
    #               savepath=os.path.join(dir,figuredir,field))
    dirnamesplit = os.path.abspath(file.pathname).split('_')
    if dirnamesplit[-1] in ['h1','hinf', 'h0'] and dirnamesplit[-2] in ['b1']:
        sim_name = dirnamesplit[-1]
    else:
        sim_name = dirnamesplit[-2] + '_' + dirnamesplit[-1]

    if field in ['particle_gamc', 'particle_age', 'particle_nuc', 'particle_magp']:
        ad = ds.all_data()
        #if proj_axis == 'x':
        #    xfield = (ptype, 'particle_position_y')
        #    yfield = (ptype, 'particle_position_z')
        #else:
        #    raise ValueError
        #cfield = (ptype, field)
        #plot = yt.ParticlePlot(ds, xfield, yfield, cfield)
        #plot.zoom(zoom_fac)
        #plot.set_zlim(cfield, 1E2, 1E5)
        #plot.set_cmap(cfield, 'algae')
        #plot.annotate_timestamp(corner='upper_left', time_format="{time:6.3f} {units}",
        #                        time_unit='Myr', text_args={'color':'k'})
        filter = ad[ptype, 'particle_tag'] % 5 == 0

        if field == 'particle_gamc':
            fdata = np.log10(np.abs(ad[ptype, 'particle_gamc'][filter]))
            vmin=2.5; vmax=4; cmap='algae'
            cblabel=u'log $\gamma_c$'
        elif field == 'particle_age':
            fdata = (ds.current_time.v - ad[ptype, 'particle_tadd'][filter])/ds.current_time.v
            vmin=0; vmax=1; cmap='algae_r'
            cblabel='normalized age'
        elif field == 'particle_nuc':
            B = np.sqrt(ad[(ptype, 'particle_magx')][filter]**2
                       +ad[(ptype, 'particle_magy')][filter]**2
                       +ad[(ptype, 'particle_magz')][filter]**2)*np.sqrt(4.0*np.pi)
            B = ad.apply_units(B, 'gauss')
            # Cutoff frequency
            fdata = np.log10(3.0*ad[(ptype, 'particle_gamc')][filter]**2*e*B/(4.0*np.pi*me*c))
            vmin=6; vmax=10; cmap='algae'
            cblabel=u'log $\\nu_c$'
        elif field == 'particle_magp':
            magp = (ad[(ptype, 'particle_magx')][filter]**2
                   +ad[(ptype, 'particle_magy')][filter]**2
                   +ad[(ptype, 'particle_magz')][filter]**2)/2
            magp = ad.apply_units(magp, 'erg/cm**3')
            fdata = np.log10(magp)
            vmin=-12; vmax=-10; cmap='algae'
            cblabel=u'log$P_B$'

        if proj_axis == 'x':
            xaxis = ad[ptype, 'particle_position_y'][filter]/3.08567758E21
            yaxis = ad[ptype, 'particle_position_z'][filter]/3.08567758E21
            fig=plt.figure(figsize=(8,12))
            xlim=ds.domain_width[1].in_units('kpc')/zoom_fac/2.0
            ylim=ds.domain_width[2].in_units('kpc')/zoom_fac/2.0
        elif proj_axis == 'y':
            xaxis = ad[ptype, 'particle_position_z'][filter]/3.08567758E21
            yaxis = ad[ptype, 'particle_position_x'][filter]/3.08567758E21
            fig=plt.figure(figsize=(12,6))
            xlim=ds.domain_width[2].in_units('kpc')/zoom_fac/2.0
            ylim=ds.domain_width[0].in_units('kpc')/zoom_fac/2.0
        elif proj_axis == 'z':
            xaxis = ad[ptype, 'particle_position_x'][filter]/3.08567758E21
            yaxis = ad[ptype, 'particle_position_y'][filter]/3.08567758E21
            fig=plt.figure(figsize=(8,7))
            xlim=ds.domain_width[0].in_units('kpc')/zoom_fac/2.0
            ylim=ds.domain_width[1].in_units('kpc')/zoom_fac/2.0

        ax=fig.add_subplot(111)
        ax.set_xlim(-xlim,xlim)
        ax.set_ylim(-ylim,ylim)
        ax.set_xlabel('kpc')
        ax.set_ylabel('kpc')
        ax.annotate('%6.3f Myr' % (float(ds.current_time)/3.15569E13),\
                    (1,0), xytext=(0.05, 0.96),  textcoords='axes fraction',\
                    horizontalalignment='left', verticalalignment='center')
        ax.annotate(sim_name, (1,0), xytext=(0.85, 0.96), textcoords='axes fraction',\
                    horizontalalignment='left', verticalalignment='center')

        sc=ax.scatter(xaxis,yaxis,s=1,c=fdata,linewidth=0,cmap=cmap,vmin=vmin,vmax=vmax,alpha=0.8)
        try:
            cb=plt.colorbar(sc)
            cb.set_label(cblabel)
        except:
            pass
        plt.tight_layout()
        plt.savefig(os.path.join(file.pathname,figuredir,ptype+'_'+field.strip(),file.filename), dpi=150)

        #plot.save(os.path.join(file.pathname,figuredir,ptype+'_'+field.strip()))

    else:
        plotSliceField(ds, zoom_fac=zoom_fac, center=center, proj_axis=proj_axis, field=field,\
                   plotgrid=fields_grid, plotvelocity=fields_velocity, \
                   annotate_particles=fields_part,annotate_part_info=False, sim_name=sim_name,\
                   savepath=os.path.join(file.pathname,figuredir,field.strip()))
    return sim_name, ds.basename[-4: ], field


def worker_test(file, field, proj_axis, zoom_fac):
    ds = yt.load(file.fullpath)
    return file.filename, field, proj_axis, zoom_fac


def tasks_gen(dirs, i0, fields, proj_axes, zoom_facs):
    #sys.stdout.write("Plotting %s... \n" % (file.filename))
    for dir in dirs:
        files = rescan(dir, False)[i0:]
        for zoom_fac in zoom_facs:
            for proj_axis in proj_axes:
                for file in files:
                    for field in fields:
                        yield file, field, proj_axis, zoom_fac


def init():
    for dir in dirs:
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
tasks = tasks_gen(dirs, i0, fields, proj_axes, zoom_facs)

results = MPI_taskpull2.taskpull(worker_fn, tasks, initialize=init, print_result=True)
#print results
