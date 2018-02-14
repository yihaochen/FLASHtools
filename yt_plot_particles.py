#!/usr/bin/env python
# Load particle information directly from hdf5 files
# Does not require plot files
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
from particles.particle_filters import *
from synchrotron.yt_synchrotron_emissivity import setup_part_file

dirs = ['./']
regex = '*hdf5_part_[0-9][0-9][0-9][0-9]_updated'
#regex = 'MHD_Jet*_hdf5_plt_cnt_[0-9][0-9][0-9][0-9]'
files = None
zoom_facs = [8]
proj_axes= ['y']
figuredirtemplate = 'figures%s_zoom%i'
ptypes = ['jetp']

fields = ['particle_gamc_dtau', 'particle_nuc_dtau']
fields += ['particle_den1']


def rescan(dir, printlist=False):
    files = util.scan_files(dir, regex=regex, walk=True, printlist=printlist, reverse=False)
    return files

e  = yt.utilities.physical_constants.elementary_charge #4.803E-10 esu
me = yt.utilities.physical_constants.mass_electron #9.109E-28
c  = yt.utilities.physical_constants.speed_of_light #2.998E10

def worker_fn(file, field, proj_axis, zoom_fac, ptype):
    ds = yt.load(file.fullpath)
    ds.add_particle_filter(ptype)

    axis = '_%s' % proj_axis if proj_axis != 'x' else ''
    figuredir = figuredirtemplate % (axis, zoom_fac)
    if file.pathname.split('/')[-1] == 'data':
        absfiguredir = os.path.join(os.path.dirname(file.pathname), figuredir)
    else:
        absfiguredir = os.path.join(file.pathname, figuredir)

    if file.pathname.split('/')[-1] == 'data':
        dirnamesplit = os.path.abspath(os.path.dirname(file.pathname)).split('_')
    else:
        dirnamesplit = os.path.abspath(file.pathname).split('_')
    if dirnamesplit[-1] in ['h1','hinf', 'h0'] and dirnamesplit[-2] in ['b1']:
        sim_name = dirnamesplit[-1]
    else:
        sim_name = dirnamesplit[-2] + '_' + dirnamesplit[-1]

    if 'particle' in field:
        ad = ds.all_data()
        filter = ad[ptype, 'particle_position_y'] < 0.24*yt.units.kpc.in_units('cm')

        if field in ['particle_gamc_dtau', 'particle_nuc_dtau', 'particle_gamc', 'particle_nuc']:

            if 'dtau' in field:
                gamc = (ad[ptype, 'particle_dens']/ad[ptype, 'particle_den0'])**(1./3.) \
                       / ad[ptype, 'particle_dtau']
            else:
                gamc = ad[ptype, 'particle_gamc']


            if 'particle_gamc' in field:
                fdata = np.log10(np.abs(gamc[filter]))
                vmin=2.5; vmax=4; cmap='arbre'
                cblabel=u'log $\gamma_c$'

            elif 'particle_nuc' in field:
                B = np.sqrt(ad[(ptype, 'particle_magx')][filter]**2
                           +ad[(ptype, 'particle_magy')][filter]**2
                           +ad[(ptype, 'particle_magz')][filter]**2)*np.sqrt(4.0*np.pi)
                B = ad.apply_units(B, 'gauss')
                # Cutoff frequency
                fdata = np.log10(3.0*gamc[filter]**2*e*B/(4.0*np.pi*me*c))
                vmin=7; vmax=10; cmap='arbre'
                cblabel=u'cutoff frequency $\\nu_c$'

            if 'dtau' in field:
                cblabel=cblabel+' (dtau)'

        elif field == 'particle_age':
            fdata = (ds.current_time.v - ad[ptype, 'particle_tadd'][filter])/ds.current_time.v
            vmin=0; vmax=1; cmap='arbre_r'
            cblabel='normalized age'
        elif field == 'particle_magp':
            magp = (ad[(ptype, 'particle_magx')][filter]**2
                   +ad[(ptype, 'particle_magy')][filter]**2
                   +ad[(ptype, 'particle_magz')][filter]**2)/2
            magp = ad.apply_units(magp, 'erg/cm**3')
            fdata = np.log10(magp)
            vmin=-12; vmax=-10; cmap='arbre'
            cblabel=u'log$P_B$'
        elif field == 'particle_den1':
            fdata = np.log10(ad[(ptype, 'particle_den1')][filter])
            vmin=-27; vmax=-25; cmap='arbre'
            cblabel=u'density when leaving jet $\\rho_1$ (g/cm$^3$)'

        if proj_axis == 'x':
            xaxis = ad[ptype, 'particle_position_y'][filter]/3.08567758E21
            yaxis = ad[ptype, 'particle_position_z'][filter]/3.08567758E21
            fig=plt.figure(figsize=(4,7))
            xlim=ds.domain_width[1].in_units('kpc')/zoom_fac/2.0
            ylim=ds.domain_width[2].in_units('kpc')/zoom_fac/2.0
            xlabel ='y'
            ylabel ='z'
        elif proj_axis == 'y':
            xaxis = ad[ptype, 'particle_position_z'][filter]/3.08567758E21
            yaxis = ad[ptype, 'particle_position_x'][filter]/3.08567758E21
            fig=plt.figure(figsize=(10,6))
            xlim=ds.domain_width[2].in_units('kpc')/zoom_fac/2.0
            ylim=15*yt.units.kpc
            xlabel ='z'
            ylabel ='x'
        elif proj_axis == 'z':
            xaxis = ad[ptype, 'particle_position_x'][filter]/3.08567758E21
            yaxis = ad[ptype, 'particle_position_y'][filter]/3.08567758E21
            fig=plt.figure(figsize=(8,7))
            xlim=ds.domain_width[0].in_units('kpc')/zoom_fac/2.0
            ylim=ds.domain_width[1].in_units('kpc')/zoom_fac/2.0
            xlabel ='x'
            ylabel ='y'

        ax=fig.add_subplot(111)
        sc=ax.scatter(xaxis,yaxis,s=1.5,c=fdata,linewidth=0,cmap=cmap,vmin=vmin,vmax=vmax,alpha=0.6)
        ax.set_xlim(-0,xlim)
        ax.set_ylim(-ylim,ylim)
        ax.set_xlabel(xlabel+' (kpc)')
        ax.set_ylabel(ylabel+' (kpc)')
        ax.set_aspect('equal', adjustable='box')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", 0.2, pad=0)
        cb=plt.colorbar(sc, cax=cax)
        cb.ax.tick_params(direction='in')
        cb.set_label(cblabel)
        if 'particle_nuc' in field:
            cb.set_ticks([7,8,9,10])
            cb.set_ticklabels(['10 MHz', '100 MHz', '1 GHz', '10 GHz'])
        ax.annotate('%6.3f Myr' % (float(ds.current_time)/3.15569E13),\
                    (1,0), xytext=(0.05, 0.9),  textcoords='axes fraction',\
                    horizontalalignment='left', verticalalignment='center')
        #ax.annotate(sim_name+'\n'+ptype, (1,0), xytext=(0.8, 0.96), textcoords='axes fraction',\
        #            horizontalalignment='left', verticalalignment='center')

        plt.tight_layout()
        plt.savefig(os.path.join(absfiguredir,ptype+'_'+field.strip(),file.filename+'.png'),\
                dpi=200, bbox_inches='tight')
        plt.close(fig)

    return sim_name, ds.basename[-12:-8], field


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
