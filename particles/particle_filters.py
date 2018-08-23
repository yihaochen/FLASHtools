import numpy as np
import yt
from tools import calcDen0

# Density tolerence for selecting core and sheth particles
deltaden0 = 5E-31
c = yt.physical_constants.speed_of_light
lobe_v = 0.05*c

@yt.particle_filter(name='metal', requires=["particle_type"], filtered_type="io")
def metal(pfilter, data):
    metal = np.isclose(data[(pfilter.filtered_type, "particle_type")], 2.0)
    return metal


@yt.particle_filter(name='jet', requires=["particle_shok", "particle_dens"], filtered_type="io")
def jet(pfilter, data):
    try:
        jet = np.logical_and(np.isclose(data[(pfilter.filtered_type, "particle_shok")], 0.0),
                             np.isclose(data[(pfilter.filtered_type, "particle_type")], 1.0))
    except:
        jet = np.isclose(data[(pfilter.filtered_type, "particle_shok")], 0.0)
    finally:
        dens = data[(pfilter.filtered_type, "particle_dens")] > 0.0
        fil = np.logical_and(jet, dens)
    return fil


@yt.particle_filter(name='shok', requires=["particle_shok"], filtered_type="io")
def shok(pfilter, data):
    shok = np.isclose(data[(pfilter.filtered_type, "particle_shok")], 1.0)
    #gamc = np.logical_and((data[(pfilter.filtered_type, "particle_gamc")] > 0.0),\
    #                      (data[(pfilter.filtered_type, "particle_gamc")] < 1E40))
    #fil = np.logical_and(shok, gamc)
    return shok


# Lobe particles that include both particles injected at the core and the sheth
@yt.particle_filter(name='slow', requires=["particle_velocity_magnitude"], filtered_type="jet")
def slow(pfilter, data):
    c = yt.physical_constants.speed_of_light
    slow = data[pfilter.filtered_type, "particle_velocity_magnitude"] < lobe_v
    return slow

@yt.particle_filter(name='fast', requires=["particle_velocity_magnitude"], filtered_type="jet")
def fast(pfilter, data):
    c = yt.physical_constants.speed_of_light
    fast = data[pfilter.filtered_type, "particle_velocity_magnitude"] > lobe_v
    return fast


# Particles injected at the core of the jet
@yt.particle_filter(name='jetp', requires=["particle_den0"], filtered_type="jet")
def jetp(pfilter, data):
    # Keep only the particles injected at the core of the jet identified by den0
    den0 = calcDen0(data, ptype=pfilter.filtered_type)
    core = np.logical_and((data[pfilter.filtered_type, 'particle_den0']>den0-deltaden0),\
                          (data[pfilter.filtered_type, 'particle_den0']<den0+deltaden0))
    return core


# Lobe particles that include only particles injected at the core
@yt.particle_filter(name='lobe', requires=["particle_velocity_magnitude"], filtered_type="jetp")
def lobe(pfilter, data):
    lobe = data[pfilter.filtered_type, "particle_velocity_magnitude"] < lobe_v
    return lobe


# Particles injected at the outskirts of the jet
@yt.particle_filter(name='shth', requires=["particle_den0", "particle_shok"], filtered_type="jet")
def shth(pfilter, data):
    # Keep only the particles injected at the core of the jet identified by den0
    den0 = calcDen0(data, ptype=pfilter.filtered_type)
    sheth = np.logical_or((data[pfilter.filtered_type, 'particle_den0']<den0-deltaden0),\
                           (data[pfilter.filtered_type, 'particle_den0']>den0+deltaden0))
    return sheth


