#!/bin/env python3
import re
import glob
import sys
import os

timelimit = int(sys.argv[1]) if len(sys.argv) > 1 else 48*60*60

plotfiles = glob.glob('./*_hdf5_plt_cnt_????') + glob.glob('./data/*_hdf5_plt_cnt_????')
partfiles = glob.glob('./*_hdf5_part_????') + glob.glob('./data/*_hdf5_part_????')
chkfiles = glob.glob('./*_hdf5_chk_????') + glob.glob('./data/*_hdf5_chk_????')

if len(chkfiles)*len(plotfiles) > 0:
    last_chkpoint_time = max([os.path.getmtime(f) for f in chkfiles])
    for f in sorted(plotfiles):
        # Don't consider plot files generated after the last checkpoint file
        if os.path.getmtime(f) > last_chkpoint_time:
            plotfiles.remove(f)
        # Don't consider forced plot files
        elif 'forced' in f:
            plotfiles.remove(f)

    for f in sorted(partfiles):
        if os.path.getmtime(f) > last_chkpoint_time:
            partfiles.remove(f)

    plotfileNumber = max([int(fn[-4:]) for fn in plotfiles]) + 1
    partfileNumber = max([int(fn[-4:]) for fn in partfiles]) + 1
    chkfileNumber = max([int(fn[-4:]) for fn in chkfiles])

    #assert plotFileNumber == partFileNumber

    restart = '.true.'
else:
    print('Setting flash.par for starting from scratch')
    restart = '.false.'
    plotfileNumber = 0
    partfileNumber = 0
    chkfileNumber = 0

timelimit -= max(plotfileNumber//10*4, 60)

with open('flash.par', 'r') as f:
    pars = f.read()

pars = re.sub(r'(wall_clock_time_limit[ \t]*?=)[ \t]*?(\d+)', r'\1 %i' % timelimit, pars)
pars = re.sub(r'(restart[ \t]*?=)[ \t]*?\..+?\.', r'\1 %s' % restart, pars)
pars = re.sub(r'(plotFileNumber[ \t]*?=)[ \t]*?(\d+)', r'\1 %i' % plotfileNumber, pars)
pars = re.sub(r'(particleFileNumber[ \t]*?=)[ \t]*?(\d+)', r'\1 %i' % partfileNumber, pars)
pars = re.sub(r'(checkpointFileNumber[ \t]*?=)[ \t]*?(\d+)', r'\1 %i' % chkfileNumber, pars)

with open('flash.par', 'w') as f:
    f.write(pars)
