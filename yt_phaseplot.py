#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import yt
import os
import sys
import util
import MPI_taskpull2
import logging
logging.getLogger('yt').setLevel(logging.ERROR)
import matplotlib.pyplot as plt
import numpy as np


dir = './'
regex = 'MHD_Jet*_hdf5_plt_cnt_[0-9][0-9][0-9][0-9]'
#regex = 'MHD_Jet*_hdf5_plt_cnt_00[0-3][0-9]'
files = None


def rescan(printlist=False):
    files = util.scan_files(dir, regex=regex, walk=True, printlist=printlist, reverse=False)
    return files

e  = yt.utilities.physical_constants.elementary_charge #4.803E-10 esu
me = yt.utilities.physical_constants.mass_electron #9.109E-28
c  = yt.utilities.physical_constants.speed_of_light #2.998E10

def worker_fn(file):
    particle_path = file.fullpath.replace('plt_cnt', 'part')\
                    if 'plt_cnt' in file.fullpath\
                    else file.fullpath.replace('chk', 'part')
    if os.path.exists(particle_path):
        ds = yt.load(file.fullpath, particle_filename=particle_path)
    else:
        ds = yt.load(file.fullpath)

    figuredir = 'phaseplot'

    ad = ds.all_data()
    v=3E9
    g=1.33333
    r=7.5E20
    bf=1.875E20
    M = 5.0+5.0*np.cos(0.5*np.pi*np.clip((ad['all', 'particle_tadd']-1.58E13)/1.58E13, -1.0, 0.0))
    den0 = 0.5*1E45/np.pi/v**3/( 0.5*r*r*(1.+1./M**2/(g-1.)) + r*bf*(0.3125+1./M**2/(g-1.))\
                             + bf*bf*(0.06056+0.29736/M**2/(g-1.)) )
    jetfil = np.logical_and((ad['all', 'particle_tag'] % 3 == 0),\
                            (ad['all', 'particle_shok']==0))
    corefil = np.logical_and((ad['all', 'particle_den0']>den0-3E-30),\
                            (ad['all', 'particle_den0']<den0+3E-30))
    #initfil = ad['all', 'particle_tadd'] < 1.58E13
    #filter = np.logical_or(initfil, np.logical_and(jetfil, corefil))
    filter = np.logical_and(jetfil, corefil)
    #filter = ad['all', 'particle_tag'] % 5 == 0

    x = ad['all', 'particle_tadd'][filter]/3.15569E13
    y = np.abs(ad['all', 'particle_posz'][filter])/3.08567758E21
    #y = np.log10(ad['all', 'particle_den0'][filter])
    #y = np.log10(ad['all', 'particle_dens'][filter]/ad['all', 'particle_den0'][filter])#/3.08567758E21
    #pos = ad['all', 'particle_position'][filter]
    #y = ds.find_field_values_at_points('jet ', pos)
    #z = np.log10(ad['all', 'particle_gamc'][filter])
    B = np.sqrt(ad[('all', 'particle_magx')][filter]**2
               +ad[('all', 'particle_magy')][filter]**2
               +ad[('all', 'particle_magz')][filter]**2)*np.sqrt(4.0*np.pi)
    B = ad.apply_units(B, 'gauss')
    #y = np.log10(B)
    # Cutoff frequency
    nuc = np.log10(3.0*ad[('all', 'particle_gamc')][filter]**2*e*B/(4.0*np.pi*me*c))

    #y = np.log10(ad['all', 'particle_magx'][filter]**2\
    #             +ad['all', 'particle_magy'][filter]**2\
    #             +ad['all', 'particle_magz'][filter]**2)
    fig=plt.figure(figsize=(12,8))
    sc=plt.scatter(x,y,c=nuc,s=1,linewidth=0,vmin=6,vmax=12)
    plt.xlim(-0.1,20)
    plt.ylim(-1,50)
    #plt.ylim(-25.78,-25.70)
    plt.xlabel('t (Myr)')
    plt.ylabel('z (kpc)')
    #plt.ylabel(u'log $\\rho_0$ (g/cm)')
    cb=plt.colorbar(sc)
    cb.set_label(u'log $\\nu_c$')

    plt.tight_layout()
    plt.savefig(os.path.join(dir,figuredir,file.filename), dpi=150)


    return ds.basename[-4: ]


def tasks_gen(i0):
    #sys.stdout.write("Plotting %s... \n" % (file.filename))
    files = rescan(False)[i0:]
    for file in files:
        yield file


def init():
    figuredir = 'phaseplot'
    if not os.path.exists(os.path.join(dir,figuredir)):
        os.mkdir(os.path.join(dir, figuredir))

i0 = int(sys.argv[1]) if len(sys.argv) > 1 else 0
tasks = tasks_gen(i0)

results = MPI_taskpull2.taskpull(worker_fn, tasks, initialize=init, print_result=True)
#print results
