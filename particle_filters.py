import numpy as np
import yt
from tools import calcDen0

deltaden0 = 5E-31

@yt.particle_filter(name='jetp', requires=["particle_shok"], filtered_type="io")
def jetp(pfilter, data):
    den0 = calcDen0(data['io', 'particle_tadd'])
    jet = np.logical_and((data[("io", "particle_shok")] == 0),
                         (data[("io", "particle_gamc")] > 0.0))
    core = np.logical_and((data['io', 'particle_den0']>den0-deltaden0),\
                          (data['io', 'particle_den0']<den0+deltaden0))
    filter = np.logical_and(jet, core)
    return filter

@yt.particle_filter(name='shok', requires=["particle_shok"], filtered_type="io")
def shok(pfilter, data):
    shok = (data[("io", "particle_shok")] == 1)
    gamc = np.logical_and((data[("io", "particle_gamc")] > 0.0),\
                          (data[("io", "particle_gamc")] < 1E40))
    filter = np.logical_and(shok, gamc)
    return filter

@yt.particle_filter(name='jnsp', requires=["particle_shok"], filtered_type="io")
def jnsp(pfilter, data):
    den0 = calcDen0(data['io', 'particle_tadd'])
    jet = np.logical_and((data[("io", "particle_shok")] == 0),
                         (data[("io", "particle_gamc")] > 0.0))
    core = np.logical_and((data['io', 'particle_den0']>den0-deltaden0),\
                          (data['io', 'particle_den0']<den0+deltaden0))
    jetcore = np.logical_and(jet, core)
    shok = np.logical_and((data[("io", "particle_shok")] == 1),\
                          (data[("io", "particle_gamc")] > 0.0))
    filter = np.logical_or(jetcore, shok)
    return filter

@yt.particle_filter(name='lobe', requires=["particle_shok"], filtered_type="io")
def lobe(pfilter, data):
    c = yt.physical_constants.speed_of_light
    den0 = calcDen0(data['io', 'particle_tadd'])
    lobe = np.logical_and((data[("io", "particle_shok")] == 0),
                         (data["io", "particle_velocity_magnitude"] < 0.01*c))
    core = np.logical_and((data['io', 'particle_den0']>den0-deltaden0),\
                          (data['io', 'particle_den0']<den0+deltaden0))
    lobecore = np.logical_and(lobe, core)
    gamc = np.logical_and((data[("io", "particle_gamc")] > 0.0),\
                          (data[("io", "particle_gamc")] < 1E40))
    filter = np.logical_and(lobecore, gamc)
    return filter

@yt.particle_filter(name='lnsp', requires=["particle_shok"], filtered_type="io")
def lnsp(pfilter, data):
    c = yt.physical_constants.speed_of_light
    den0 = calcDen0(data['io', 'particle_tadd'])
    lobe = np.logical_and((data[("io", "particle_shok")] == 0),
                         (data["io", "particle_velocity_magnitude"] < 0.01*c))
    core = np.logical_and((data['io', 'particle_den0']>den0-deltaden0),\
                          (data['io', 'particle_den0']<den0+deltaden0))
    lobecore = np.logical_and(lobe, core)
    gamc = np.logical_and((data[("io", "particle_gamc")] > 0.0),\
                          (data[("io", "particle_gamc")] < 1E40))
    shok = (data[("io", "particle_shok")] == 1)
    filter = np.logical_and(np.logical_or(lobecore, shok), gamc)
    return filter
