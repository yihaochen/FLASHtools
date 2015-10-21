#!/usr/bin/env python
import time
t0 = time.time()
import matplotlib
matplotlib.use('Agg')
import yt
import os
import sys
sys.path.append(os.getenv('HOME') + '/lib/util')
import util
#from mpi_taskpull import taskpull
import logging
logging.getLogger('yt').setLevel(logging.ERROR)
from mpi4py import MPI
import matplotlib.pyplot as plt
import numpy as np

from tools import calcNozzleCoords, read_par
from plotSlices import plotSliceField
from plotProjections import plotProjectionField

def enum(*sequential, **named):
    """Handy way to fake an enumerated type in Python
    http://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python
    """
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

# Define MPI message tags
tags = enum('READY', 'DONE', 'EXIT', 'START')

# Initializations and preliminaries
comm = MPI.COMM_WORLD   # get MPI communicator object
size = comm.size        # total number of processes
rank = comm.rank        # rank of this process
status = MPI.Status()   # get MPI status object
name=MPI.Get_processor_name()
dir = './'
regex = 'MHD_Jet*_hdf5_plt_cnt_[0-9][0-9][0-9][0-9]'
#regex = 'MHD_Jet_10Myr_hdf5_plt_cnt_[0-9][0-9][0-9][0-9]'
#regex = 'MHD_Jet_hdf5_chk_001[0-9]'
files = None
zoom_facs = [8,16,128]
proj_axes= ['x']


#annotate_particles = True if zoom_fac >= 2 else False
fields_velocity = ['velocity_z']
fields_part = None
fields = ['density', 'pressure', 'temperature', 'velocity_y', 'velocity_z', 'jet ',\
          'magnetic_field_x', 'magnetic_field_y', 'magnetic_field_z', 'magnetic_pressure']
#fields = ['density', 'pressure','shok', 'velocity_y', 'velocity_z', 'jet ']
#fields = ['density', 'pressure', 'temperature', 'velocity_magnitude', 'velocity_z']
#fields = ['magnetic_pressure']
#fields = ['density', 'pressure', 'temperature', 'velocity_y', 'velocity_z']
#fields = ['pressure']
#fields = ['magnetic_field_z']
#fields = ['velocity_y']
#fields = ['shok']
#fields = ['plasma_beta']
#fields = ['particle_gamc']


def rescan(printlist=False):
    files = util.scan_files(dir, regex=regex, walk=True, printlist=printlist, reverse=False)
    return files

def worker_fn(file, field, proj_axis, zoom_fac):
    t_worker_0 = time.time()
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
             else ['velocity_z']
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
        plotSliceField(ds, zoom_fac=zoom_fac, center=center, proj_axis=proj_axis, field=field,\
                   plotgrid=fields_grid, plotvelocity=fields_velocity, nozzleCoords=nozzleCoords, \
                   annotate_particles=fields_part,annotate_part_info=False,\
                   savepath=os.path.join(dir,figuredir,field))
    t_worker_1 = time.time()
    elapsed_time = '%5.1f s' % (t_worker_1-t_worker_0)
    return ds.basename[-4: ], '%18s' % field, elapsed_time




def worker_test(file, field):
    ds = yt.load(file.fullpath)
    return (file.filename, field)

def tasks_gen(files, fields, proj_axes, zoom_facs):
    #sys.stdout.write("Plotting %s... \n" % (file.filename))
    for zoom_fac in zoom_facs:
        for proj_axis in proj_axes:
            for file in files:
                for field in fields:
                    yield (file, field, proj_axis, zoom_fac)


if rank == 0:
    i0 = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    files = rescan(True)[i0:]
    tasks = tasks_gen(files, fields, proj_axes, zoom_facs)
    for zoom_fac in zoom_facs:
        for proj_axis in proj_axes:
            figuredir = 'figures_%s' % proj_axis if proj_axis != 'x' else 'figures'
            figuredir = figuredir + '_zoom%i' % zoom_fac
            if not os.path.exists(os.path.join(dir,figuredir)):
                os.mkdir(os.path.join(dir, figuredir))
            for field in fields:
                if not os.path.exists(os.path.join(dir,figuredir,field)):
                    os.mkdir(os.path.join(dir,figuredir,field))
    t1 = time.time()
    sys.stdout.write("Timer started\n")

    # Master process executes code below
    task_index = 0
    num_workers = size - 1
    closed_workers = 0
    sys.stdout.write("Master starting with %d workers\n" % num_workers)
    while closed_workers < num_workers:
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        source = status.Get_source()
        tag = status.Get_tag()
        if tag == tags.READY:
            # Worker is ready, so send it a task
            try:
                comm.send(tasks.next(), dest=source, tag=tags.START)
                #print("Sending task to worker %03d" % (source))
            except StopIteration:
                comm.send(None, dest=source, tag=tags.EXIT)
            #if task_index < len(tasks):
            #    comm.send(tasks[task_index], dest=source, tag=tags.START)
            #    print("Sending task %d to worker %d" % (task_index, source))
            #    task_index += 1
            #else:
            #    comm.send(None, dest=source, tag=tags.EXIT)
        elif tag == tags.DONE:
            results = data
            sys.stdout.write("Worker %03d returned data: %s\n" % (source, str(data)))
        elif tag == tags.EXIT:
            sys.stdout.write("Worker %03d exited. (%3d/%3d)\n" % (source, num_workers-closed_workers-1, num_workers))
            closed_workers += 1

    print("Master finishing")
    t2 = time.time()
    print 'Total time: %.2f s\n--\ninitialization: %.2f s\nploting: %.2f'\
           %  (t2-t0, t1-t0, t2-t1)
else:
    # Worker processes execute code below
    name = MPI.Get_processor_name()
    #sys.stdout.write("I am a worker with rank %03d on %s.\n" % (rank, name))
    while True:
        comm.send(None, dest=0, tag=tags.READY)
        task = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
        tag = status.Get_tag()

        if tag == tags.START:
            #sys.stdout.write("Worker %03d on %s got a job.\n" % (rank, name))
            # Do the work here
            result = worker_fn(*task)
            comm.send(result, dest=0, tag=tags.DONE)
        elif tag == tags.EXIT:
            break

    comm.send(None, dest=0, tag=tags.EXIT)
