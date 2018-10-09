import numpy as np
import os
import re
from datetime import datetime
from collections import defaultdict


class Logfile():
    def __init__(self, fname):
        self.__dict__ = defaultdict(list)
        self.sim_dict = defaultdict(list)
        self.fname = fname
        self.last_updated = None
        self.last_parsed = None
        self.dirname = fname.split('/')[-2]

    def append_sim_dict(self, startmark, endmark, startdtime, enddtime, nProc):
        # startmarks are the steps when starts (or restarts) of the simulation
        self.sim_dict['startmarks'].append(int(startmark))
        # endmarks are the steps just before the startmarks
        self.sim_dict['endmarks'].append(int(endmark))
        self.sim_dict['starttime'].append(startdtime)
        self.sim_dict['endtime'].append(enddtime)
        self.sim_dict['sim_nProc'].append(int(nProc))
        self.sim_dict['SUs'].append((enddtime - startdtime).total_seconds()/3600*nProc/48)

    def parse_lines(self, loglines):
        # Regular expressions for retrieving # of blocks
        leafblk_re = ' \[GRID amr_refine_derefine\] min leaf blks (\d+)[ ]*?max leaf blks (\d+)[ ]*?tot leaf blks[ ]*(\d+)'
        blk_re = ' \[GRID amr_refine_derefine\] min blks (\d+)[ ]*?max blks (\d+)[ ]*?tot blks (\d+)'
        step_re = ' \[ (.*) \] step: n=(\d+) t=(.+) dt=(.+)'

        strtime = None
        # Start of a simulation
        start = False
        # Continued from previous line
        continuation = False

        step_dict = {}
        for line in loglines:
            # leaf blk reporting is spread into two lines
            if continuation:
                # append the current line to the lines_buf of previous line
                lines_buf += line.lstrip()
                ma = re.match(leafblk_re, lines_buf)
                if ma:
                    step_dict['min_leaf_blks'], step_dict['max_leaf_blks'], step_dict['tot_leaf_blks'] = \
                        tuple(map(int, ma.groups()))
                continuation = False
            elif 'FLASH log file:' in line:
                # After at least one simulation
                if strtime:
                    #print(nProc, (dtime - startdtime).total_seconds()/3600*nProc/48)
                    dtime = datetime.strptime(strtime, '%m-%d-%Y  %H:%M:%S.%f')
                    self.append_sim_dict(startmark, step, startdtime, dtime, nProc)
                start = True
                startdtime = datetime.strptime(line, ' FLASH log file:  %m-%d-%Y %H:%M:%S.%f    Run number:  1\n')
            elif 'Number of MPI tasks:' in line:
                nProc = int(line.split()[-1])
                step_dict['nProc'] = nProc
            elif '[GRID amr_refine_derefine] min blks' in line:
                ma = re.match(blk_re, line)
                step_dict['min_blks'], step_dict['max_blks'], step_dict['tot_blks'] = tuple(map(int, ma.groups()))
            elif '[GRID amr_refine_derefine] min leaf blks' in line:
                lines_buf = line.strip('\n')
                continuation = True
            elif 'step: n=' in line:
                ma = re.match(step_re, line)
                strtime, step, t, dt = ma.groups()
                dtime = datetime.strptime(strtime, '%m-%d-%Y  %H:%M:%S.%f')
                step_dict['step'] = int(step)
                step_dict['t'] = float(t)
                step_dict['dt'] = float(dt)
                step_dict['step_dtime'] = dtime

                # Append the current step to the class dict
                for k, v in step_dict.items():
                    self.__dict__[k].append(v)
                if start:
                    startmark = step
                    start = False
            ma = re.match(' \[ (.*) \]', line)
            if ma:
                strtime, = ma.groups()
            # end of the for loop of lines

        # Append the information from the last simulation
        dtime = datetime.strptime(strtime, '%m-%d-%Y  %H:%M:%S.%f')
        self.append_sim_dict(startmark, step, startdtime, dtime, nProc)

    def parse_logfile(self):
        begin = datetime.now()
        # Last modified time of the logfile
        self.last_updated = datetime.fromtimestamp(os.path.getmtime(self.fname))
        if self.last_parsed and self.last_parsed > self.last_updated:
            return

        for k, v in self.__dict__.copy().items():
            if type(v) is np.ndarray:
                del self.__dict__[k]
        self.sim_dict = defaultdict(list)

        print('Parsing ', self.fname)
        with open(self.fname, 'r') as f:
            self.parse_lines(f.readlines())
        # Convert lists to numpy arrays
        for k, v in self.__dict__.items():
            if type(v) is list:
                self.__dict__[k] = np.array(v)
        # Add sim_dict to attributes
        for k, v in self.sim_dict.items():
            self.__dict__[k] = np.array(v)

        self.step_time = self.calculate_step_time()
        self.last_parsed = datetime.now()

        print('Parsed %6i steps (last n=%6i) in %.2f s' %
              (len(self.step), self.step[-1], (datetime.now()-begin).total_seconds()))

    def calculate_step_time(self):
        step_time = np.zeros(len(self.step))
        for i in np.arange(len(self.step)):
            if self.step[i] in self.startmarks:
                j = np.where(self.startmarks == self.step[i])[0][0]
                prv_time = self.starttime[j] if self.step_dtime[i] > self.starttime[j] else self.step_dtime[i-1]
            else:
                prv_time = self.step_dtime[i-1]
            step_time[i] = (self.step_dtime[i] - prv_time).total_seconds()

        return step_time
