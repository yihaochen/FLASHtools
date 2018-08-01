#!/bin/env python3
import re
import glob

try:
    plotfiles = glob.glob('./*_hdf5_plt_cnt_????')
    plotFileNumber = max([int(fn[-4:]) for fn in plotfiles]) + 1

    partfiles = glob.glob('./*_hdf5_part_????')
    partFileNumber = max([int(fn[-4:]) for fn in partfiles]) + 1

    chkfiles = glob.glob('./*_hdf5_chk_????')
    chkFileNumber = max([int(fn[-4:]) for fn in chkfiles])

    #assert plotFileNumber == partFileNumber

    with open('flash.par', 'r') as f:
        pars = f.read()

    pars = re.sub(r'(plotFileNumber[ \t]+=) (\d+)', r'\1 %i' % plotFileNumber, pars)
    pars = re.sub(r'(particleFileNumber[ \t]+=) (\d+)', r'\1 %i' % partFileNumber, pars)
    pars = re.sub(r'(checkpointFileNumber[ \t]+=) (\d+)', r'\1 %i' % chkFileNumber, pars)

    with open('flash.par', 'w') as f:
        f.write(pars)
except ValueError:
    pass

