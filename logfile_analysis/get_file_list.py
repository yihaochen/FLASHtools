#!/usr/bin/python
import sys
from logfile_parser import Logfile
import numpy as np

if __name__ == "__main__":
    fname = sys.argv[1]

    log = Logfile(fname)
    file_dict = log.parse_lines_hdf5_files()
    all_file = sorted(file_dict['checkpoint'] + file_dict['plotfile'] + file_dict['particles'], key=lambda x:x[2])
    np.savetxt(fname.replace('.log', '_file_list.txt').replace('/log/', '/'),
            np.array(all_file),
            fmt='%s %30s %7s %12s %6s',
            header='timestamp                filename                      nstep  time         time (Myr)')
