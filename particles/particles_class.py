import numpy as np
import util
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['figure.dpi'] = 150
import matplotlib.pyplot as plt
import h5py
import yt

me = yt.utilities.physical_constants.mass_electron.v
c  = yt.utilities.physical_constants.speed_of_light.v
e  = yt.utilities.physical_constants.elementary_charge.v
Myr= yt.units.Myr.in_units('s').v

_fields_list = [
#    'tag', 'tadd', 'den0',
    'posx', 'posy', 'posz',
    'velx', 'vely', 'velz',
    'magx', 'magy', 'magz',
    'dens', 'gamc',
#    'den1', 'dtau',
    'shok',
    'tau0', 'cmb0', 'ict0',
    'ind1', 'tau1', 'den1', 'cmb1', 'ict1', 'tad1',
    'ind2', 'tau2', 'den2', 'cmb2', 'ict2', 'tad2',
    'ind3', 'tau3', 'den3', 'cmb3', 'ict3', 'tad3',
    'ind4', 'tau4', 'den4', 'cmb4', 'ict4', 'tad4',
    'whch', 'shks', 'jet'
]

def find_part_ind(h5file, tag, tadd):
    """
    Find the index of the particle matching tadd and tag.
    """
    tp = h5file['tracer particles']
    colname = [item[0].strip() for item in h5file['particle names']]
    itadd = colname.index('tadd')
    itag = colname.index('tag')
    ind_tadd = np.isclose(tp[:,itadd]-tadd, 0.0)
    ind_tag = np.in1d(tp[:,itag], tag)
    ind = np.where(np.logical_and(ind_tadd, ind_tag))[0]
    #print np.where(np.in1d(tp[:,findices['tadd']], tadd))
    if len(ind) == 0:
        return None
    elif len(ind) == 1:
        return ind
    else:
        raise IndexError


def find_part_indices(h5file, tags, tadds):
    """
    Find the indices of the particles matching tadds and tags. Return None if not found.
    """
    tp = h5file['tracer particles']
    colname = [item[0].decode().strip() for item in h5file['particle names']]
    itadd = colname.index('tadd')
    itag = colname.index('tag')

    masks_tadd = np.in1d(tp[:,itadd], tadds)
    masks_tag = np.in1d(tp[:,itag], tags)

    indices_temp = np.where(np.logical_and(masks_tadd, masks_tag))[0]

    # Make sure we find the same number of particles as we have
    assert len(indices_temp) <= len(tags),\
           'len(indices_temp) = %i, len(tags) = %i' % (len(indices_temp), len(tags))

    indices = [None]*len(tags)

    for i, (tag, tadd) in enumerate(zip(tags, tadds)):
        for ind in indices_temp:
            if tag == tp[ind,itag] and np.isclose(tadd, tp[ind,itadd]):
                indices[i] = ind
                break

    return indices

class Particle():
    def __init__(self, tag, tadd, grid_fields=[]):
        self.tadd = tadd
        self.tag = tag
        self.den0 = -1
        self.time = []
        #nfiles = len(partfiles)
        self.grid_fields = grid_fields

        if grid_fields:
            for grid_field in grid_fields:
                setattr(self, grid_field, [])

        for field in _fields_list:
            setattr(self, field, [])

    def read_from_h5file(self, ind, h5file, ds=None):
        colname = [item[0].decode().strip() for item in h5file['particle names']]
        try:
            findices = {f: colname.index(f) for f in _fields_list}
        except ValueError as err:
            print(err, colname)
            raise ValueError

        # If the supplied particle is found in this particle file
        if ind:
            self.time.append( h5file['real scalars'].value[0][1]/Myr )
            tp = h5file['tracer particles']
            # Assign den0 and den1 for the first time and make sure den0 are the same
            if self.den0 < 0:
                self.den0 = tp[ind, colname.index('den0')]
                #self.den1 = tp[ind, colname.index('den1')]
            else:
            # Compare den0 and den1 otherwise
                assert self.den0 == tp[ind, colname.index('den0')]
                #assert self.den1 == tp[ind, colname.index('den1')]

            for field in _fields_list:
                getattr(self, field).append( tp[ind,findices[field]] )
            if self.grid_fields and ds:
                pos = [self.posx[-1], self.posy[-1], self.posz[-1]]
                for grid_field in self.grid_fields:
                    #print ds, grid_field, pos
                    getattr(self, grid_field).append(ds.find_field_values_at_point(grid_field, pos)[0])

    def add(self, other_particle):
        assert self.tadd == other_particle.tadd
        assert self.tag == other_particle.tag
        self.time.extend(other_particle.time)
        if self.den0 < 0 and other_particle.den0 > 0:
            self.den0 = other_particle.den0
            self.den1 = other_particle.den1

        for field in _fields_list + self.grid_fields:
            getattr(self, field).extend( getattr(other_particle, field) )


class Particles():
    def __init__(self, tags, tadds, grid_fields=[]):
        assert len(tadds) == len(tags)
        self.tags = tags
        self.tadds = tadds
        self.nparticles = len(tadds)
        self.grid_fields = grid_fields

        # Initialize the particle
        self.data = [Particle(tag, tadd, grid_fields=grid_fields) for tag, tadd in zip(tags, tadds)]

    def __getitem__(self, key):
        return self.data[key]

    def convert_to_ndarray(self):
        # Convert the lists to numpy arrays
        for part in self.data:
            sortind = np.argsort(part.time)
            part.time = np.array(part.time)[sortind]
            for field in _fields_list:
                setattr(part, field, np.array(getattr(part, field))[sortind])
            if self.grid_fields:
                for grid_field in self.grid_fields:
                    setattr(part, grid_field, np.array(getattr(part, grid_field))[sortind])



    def read_from_partfile(self, filepath):
        h5file = h5py.File(filepath, 'r')
        indices = find_part_indices(h5file, self.tags, self.tadds)

        for i, ind in enumerate(indices):
            ds = yt.load(filepath.replace('part', 'plt_cnt')) if self.grid_fields else None
            self.data[i].read_from_h5file(ind, h5file, ds=ds)


    def read_from_partfiles(self, partfiles):
        for pf in partfiles:
            self.read_from_partfile(pf.fullpath)
        self.convert_to_ndarray()

    def combine(self, other_particles):
        for i, other_part in enumerate(other_particles.data):
            self.data[i].add(other_part)

