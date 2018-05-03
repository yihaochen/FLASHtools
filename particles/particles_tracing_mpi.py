#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['figure.dpi'] = 150
import matplotlib.pyplot as plt
import yt
import util
import MPI_taskpull2
import numpy as np
import h5py
import random
import pickle
import os
from mpi4py import MPI

from particles.particles_class import Particles
from tools import calcDen0_2015
from particle_filters import deltaden0

me = yt.utilities.physical_constants.mass_electron.v
c  = yt.utilities.physical_constants.speed_of_light.v
e  = yt.utilities.physical_constants.elementary_charge.v
Myr= yt.units.Myr.in_units('s').v
kpc= yt.units.kpc.in_units('cm').v

dir = './'
regex = '*_hdf5_part_????'

# Index of the file in which random particles will be selected
rand_findex = 500
nparticles = 200

tags = None
tadds = None

#grid_fields = ["pressure", "grid_level", "entropy"]
grid_fields = []

def rescan(dir, printlist=False):
    files = util.scan_files(dir, regex=regex, walk=True)
    return files


if MPI_taskpull2.rank == 0:
    partfiles = rescan(dir)

    # Arbitrarily pick particles in this file
    h5file = h5py.File(partfiles[rand_findex].fullpath, 'r')
    tp = h5file['tracer particles']
    print('Choosing random particles in %s' % partfiles[rand_findex].filename)

    # Pick nparticles
    rints = sorted(random.sample(range(tp.shape[0]), nparticles))

    colname = [item[0].decode().strip() for item in h5file['particle names']]
    findices = {f: colname.index(f) for f in ['tag', 'tadd', 'shok']}

    tags =  tp[rints,colname.index('tag')]
    tadds = tp[rints,colname.index('tadd')]
    shoks = tp[rints,colname.index('shok')]

tags, tadds = MPI_taskpull2.comm.bcast((tags, tadds), root=0)

def init():
    for (tag, tadd, shok, rint) in zip(tags, tadds, shoks, rints):
        print(tag, tadd, shok, rint)


def worker_fn(filepath):
    parts = Particles(tags, tadds, grid_fields=grid_fields)
    parts.read_from_partfile(filepath)

    return parts

def tasks_gen():
    files = rescan(dir)
    for file in files[:]:
        yield file.fullpath


picklename = "Particles_tracing.pickle"
if not os.path.exists(picklename):
    tasks = tasks_gen()
    results = MPI_taskpull2.taskpull(worker_fn, tasks, initialize=init)
    # Only on rank 0 since others return None
    if results:
        collected = {}
        #parts = Particles(tags, tadds, grid_fields=grid_fields)
        path, parts = results.popitem()
        for value in results.values():
            parts.combine(value)
        parts.convert_to_ndarray()

        pickle.dump(parts, open( picklename, "wb" ))
else:
    if MPI_taskpull2.rank == 0:
        print('%s found, unpickling...' % picklename)
        parts = pickle.load(open(picklename, 'rb'))
        print('Loaded %s' % picklename)



################################################################
################################################################
v_thres = 0.5*3E9
figdir = os.path.join(dir, 'particle_tracing')
timeranges = ['alltime', 2]

def make_figdirs():
    subdirs = [figdir, os.path.join(figdir, 'alltime'),
               os.path.join(figdir, '%iMyr' % timeranges[1])]
    for subdir in subdirs:
        if not os.path.exists(subdir):
            os.mkdir(subdir)

def plot_particle_tracing(part, timerange='alltime'):
    if timerange == 'all':
        part.time -= part.time[0]
        part.time[0] += 1E-2
    fig = plt.figure(figsize=(10,6))
    #################### Density ################
    #rho_norm = 10**int(np.log10(max(part.dens)))
    #rho_norm = max(part.dens)/10
    rho_norm = 1E-28
    plt.fill_between(part.time, part.dens/rho_norm, lw=0, color='skyblue', alpha=0.7,
                     label=u'density/(%.0eg/cm$^3$)' % rho_norm)
    plt.fill_between(part.time, part.den1/rho_norm, lw=0, color='skyblue', alpha=0.5)

    #################### Velocity ################
    v_mag = np.sqrt(part.velx**2+part.vely**2+part.velz**2)
    v_norm = 6E8
    #plt.plot(part.time, v_mag/v_norm, '--', c='r', alpha=0.5, label='v/%.0e cm/s' % v_norm)
    plt.fill_between(part.time, v_mag/v_norm, color='r', lw=0, alpha=0.3, label=u'$|v| & v_z $/(%.0ecm/s)' % v_norm)
    plt.fill_between(part.time, np.abs(part.velz)/v_norm, color='r', lw=0, alpha=0.2)
    try:
        arg = np.where(v_mag<v_thres)[0][0]
        plt.axvline(part.time[arg], color='r', lw=1, alpha=0.3)
    except:
        pass
    #plt.hlines(v_thres/v_norm, part.time[0], part.time[-1], color='r', lw=1, alpha=0.3)

    #plt.twiny()

    #################### Magnetic Field ################
    B_mag = np.sqrt(part.magx**2+part.magy**2+part.magz**2)
    #B_norm = max(B_mag)/10
    B_norm = 1E-6
    plt.plot(part.time, B_mag/B_norm, '-', c='b', lw=0.8, alpha=0.5, label=u'B/(%.0f$\mu$G)' % (B_norm*1E6))

    #################### SHKS ################
    plt.plot(part.time, part.shks, '.-', ms=3, c='k', alpha=0.8, label=u'shks')
    plt.plot(part.time, part.whch, 'x', ms=3, c='k', alpha=0.8, label=u'whch')

    #################### TAU ################
    if part.tau1.size:
        tau_norm = max(part.tau1)/10
        plt.plot(part.time, part.tau1/tau_norm, '-',  c='g', label='tau1/%.0e' % (tau_norm), alpha=0.3)
        plt.plot(part.time, part.tau2/tau_norm, '--', c='g', label='tau2/%.0e' % (tau_norm), alpha=0.3)
        plt.plot(part.time, part.tau3/tau_norm, ':',  c='g', label='tau3/%.0e' % (tau_norm), alpha=0.3)
        #plt.plot(part.time, part.dtau/tau_norm, '--', c='g', label='dtau/%.0e' % (tau_norm))

        #################### IND ################
        plt.plot(part.time, part.ind1, '-',  c='purple', label='ind1', alpha=0.3)
        plt.plot(part.time, part.ind2, '--', c='purple', label='ind2', alpha=0.3)
        plt.plot(part.time, part.ind3, ':',  c='purple', label='ind3', alpha=0.3)

    #plt.hlines(0.00056966283093648636/tau_norm, part.time[0], part.time[-1], color='g', lw=0.5)
    #ind = np.argwhere(np.isclose(part.tau, 0.00056966283093648636))
    #plt.vlines(part.time[ind], 0, 12, color='g', lw=0.5)

    #################### Cutoff Gamma ################
    plt.plot(part.time, np.log10(part.gamc), '-.', c='c', label=u'log $\\gamma_{c}$', alpha=0.3)
    #gamc_dtau = (part.dens/part.den0)**(1./3.)/part.dtau
    #plt.plot(part.time, np.log10(gamc_dtau), '-', c='c', label=u'log $\\gamma_{c}$ (dtau)')


    #################### Cutoff frequency ################
    nuc = 3.0*part.gamc**2*e*B_mag/(4.0*np.pi*me*c)
    #nuc_dtau = 3.0*gamc_dtau**2*e*B_mag/(4.0*np.pi*me*c)
    plt.plot(part.time, np.log10(nuc), ':', c='gold', label=u'log $\\nu_{c}$/Hz', alpha=0.3)
    #plt.plot(part.time, np.log10(nuc_dtau), '-', c='gold', label=u'log $\\nu_{c}$/Hz (dtau)')

    #################### Cylindrical Radius ################
    #r_nozzle = 7.5E20
    #rr = np.sqrt(part.posx**2 + part.posy**2)/r_nozzle
    #plt.plot(part.time, rr, '-', c='indigo', label=u'$r_{cylindrical}$/(%ipc)' % (r_nozzle/kpc*1000))

    #################### Sampling Ticks ################
    plt.plot(part.time, [14]*len(part.time), '|', c='k')

    plt.ylim(0,14)
    plt.grid(axis='y', alpha=0.5)

    if timerange == 'alltime':
        plt.xlim(part.time[0], part.time[-1])
        plt.xlabel('time (Myr) - t0')
        plt.semilogx()
    else:
        plt.xlim(part.time[0], part.time[0]+timerange)
        plt.xlabel('time (Myr)')
    plt.legend(loc=1, ncol=2)
    plt.text(0.05, 0.92,
             't0=%.2f Myr, id=%i' % (part.tadd/Myr, part.tag),
             transform=plt.axes().transAxes)
    #plt.text(0.05, 0.87,
    #         r'$\rho_1$=%.2f*%.0e g/cm$^3$' % (part.den1/rho_norm, rho_norm),
    #         transform=plt.axes().transAxes)

    #coreden0 = calcDen0_2015(part.tadd)
    #plt.plot(part.time[0], coreden0/rho_norm, '+', c='skyblue', ms=5)
    if part.shok[0] == 1:
        plt.text(0.05, 0.82, 'shok', transform=plt.axes().transAxes)
    #elif part.den0 < coreden0 - deltaden0 or part.den0 > coreden0 + deltaden0:
    #    plt.text(0.05, 0.82, r'$\rho_0$=%.3e, $\rho_{0,core}=$%.3e' % (part.den0, coreden0), transform=plt.axes().transAxes)
    #    plt.text(0.05, 0.77, 'shelth', transform=plt.axes().transAxes)
    #else:

    #    plt.text(0.05, 0.82, r'$\rho_0$=%.3e, $\rho_{0,core}=$%.3e' % (part.den0, coreden0), transform=plt.axes().transAxes)

    if timerange == 'alltime':
        figfname = os.path.join(figdir, 'alltime/particle_tracing_%.2f_%i.png')
    else:
        figfname = os.path.join(figdir, '%iMyr/particle_tracing_%%.2f_%%i.png' % timerange)
    plt.tight_layout()
    plt.savefig(figfname % (part.tadd/Myr, part.tag))
    plt.close(fig)


def tasks_plot():
    for part in parts:
        for timerange in timeranges:
            yield part, timerange

tasks = tasks_plot()

MPI_taskpull2.comm.barrier()
results = MPI_taskpull2.taskpull(plot_particle_tracing, tasks, initialize=make_figdirs)
