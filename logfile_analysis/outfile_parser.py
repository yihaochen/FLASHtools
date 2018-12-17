import numpy as np
import os
import re
from datetime import datetime
from collections import defaultdict
from glob import glob


class Outfiles():
    def __init__(self, fname_regex):
        self.__dict__ = defaultdict(list)
        #self.last_updated = None
        #self.last_parsed = None

        self.fnames = glob(fname_regex)
        self.fnames.sort(key=os.path.getmtime)

    def parse_lines(self, lines):
        # Regular expressions for retrieving information from lines
        sci_re = '(\s?-?\d\.\d+E\+\d+)'
        out_re = '\s+(\d+) %s %s  \(%s, %s, %s\) \|  %s %s\n' % tuple([sci_re]*7)

        step_dict = {}
        for line in lines:
            ma = re.match(out_re, line)
            if ma:
                step, t, dt, x, y, z, dt_hydro, dt_part = ma.groups()

                step_dict['step'] = int(step)
                step_dict['t'] = float(t)
                step_dict['dt'] = float(dt)
                step_dict['x'] = float(x)
                step_dict['y'] = float(y)
                step_dict['z'] = float(z)
                step_dict['dt_hydro'] = float(dt_hydro)
                step_dict['dt_part'] = float(dt_part)

                # Append the current step to the class dict
                for k, v in step_dict.items():
                    self.__dict__[k].append(v)

    def parse_outfiles(self, sizelimit=None):
        begin = datetime.now()
        # Last modified time of the file
        #self.last_updated = datetime.fromtimestamp(os.path.getmtime(fname))
        #if self.last_parsed and self.last_parsed > self.last_updated:
        #    return

        for k, v in self.__dict__.copy().items():
            if type(v) is np.ndarray:
                del self.__dict__[k]

        for fname in self.fnames:
            fsize = os.path.getsize(fname)/2**20
            print('Parsing ', fname, '%.2f MB' % fsize)
            if sizelimit and fsize > sizelimit:
                print('Skipped (> %.1f MB)' % sizelimit)
            with open(fname, 'r') as f:
                self.parse_lines(f.readlines())

        #self.last_parsed = datetime.now()

        # Convert lists to numpy arrays
        for k, v in self.__dict__.items():
            if type(v) is list:
                self.__dict__[k] = np.array(v)

        print('Parsed %6i steps (last n=%6i) in %.2f s' %
              (len(self.step), self.step[-1], (datetime.now()-begin).total_seconds()))

