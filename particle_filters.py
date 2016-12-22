import numpy as np
import yt
from tools import calcDen0

deltaden0 = 5E-31

@yt.particle_filter(name='metal', requires=["particle_type"], filtered_type="io")
def metal(pfilter, data):
    fil = np.isclose(data[("io", "particle_type")], 2.0)
    return fil


@yt.particle_filter(name='jet', requires=["particle_type", "particle_shok"], filtered_type="io")
def jet(pfilter, data):
    fil = np.logical_and(np.isclose(data[("io", "particle_shok")], 0.0),
                            np.isclose(data[("io", "particle_type")], 1.0))
    return fil

@yt.particle_filter(name='shok', requires=["particle_shok"], filtered_type="io")
def shok(pfilter, data):
    shok = np.isclose(data[("io", "particle_shok")], 1.0)
    gamc = np.logical_and((data[("io", "particle_gamc")] > 0.0),\
                          (data[("io", "particle_gamc")] < 1E40))
    fil = np.logical_and(shok, gamc)
    return fil

@yt.particle_filter(name='jetp', requires=["particle_den0", "particle_shok"], filtered_type="io")
def jetp(pfilter, data):
    # Leave only the particles injected at the core of the jet identified by den0
    den0 = calcDen0(data, ptype='io')
    try:
        jet = np.logical_and(np.isclose(data[("io", "particle_shok")], 0.0),
                             np.isclose(data[("io", "particle_type")], 1.0))
    except:
        jet = np.isclose(data[("io", "particle_shok")], 0.0)
    core = np.logical_and((data['io', 'particle_den0']>den0-deltaden0),\
                          (data['io', 'particle_den0']<den0+deltaden0))
    fil = np.logical_and(jet, core)
    return fil

@yt.particle_filter(name='jnsp', requires=["particle_shok"], filtered_type="io")
def jnsp(pfilter, data):
    den0 = calcDen0(data)
    jet = np.logical_and((data[("io", "particle_shok")] == 0),
                         (data[("io", "particle_gamc")] > 0.0))
    core = np.logical_and((data['io', 'particle_den0']>den0-deltaden0),\
                          (data['io', 'particle_den0']<den0+deltaden0))
    jetcore = np.logical_and(jet, core)
    shok = np.logical_and((data[("io", "particle_shok")] == 1),\
                          (data[("io", "particle_gamc")] > 0.0))
    fil = np.logical_or(jetcore, shok)
    return fil

@yt.particle_filter(name='lobe', requires=["particle_shok"], filtered_type="io")
def lobe(pfilter, data):
    c = yt.physical_constants.speed_of_light
    den0 = calcDen0(data, ptype='io')
    try:
        jet = np.logical_and(np.isclose(data[("io", "particle_shok")], 0.0),
                             np.isclose(data[("io", "particle_type")], 1.0))
    except:
        jet = np.isclose(data[("io", "particle_shok")], 0.0)
    core = np.logical_and((data['io', 'particle_den0']>den0-deltaden0),\
                          (data['io', 'particle_den0']<den0+deltaden0))
    jetcore = np.logical_and(jet, core)

    lobe = data["io", "particle_velocity_magnitude"] < 0.01*c
    gamc = np.logical_and((data[("io", "particle_gamc")] > 0.0),\
                          (data[("io", "particle_gamc")] < 1E40))
    lobegamc = np.logical_and(lobe, gamc)

    fil = np.logical_and(jetcore, lobegamc)
    return fil

@yt.particle_filter(name='lnsp', requires=["particle_shok"], filtered_type="io")
def lnsp(pfilter, data):
    c = yt.physical_constants.speed_of_light
    den0 = calcDen0(data)
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
