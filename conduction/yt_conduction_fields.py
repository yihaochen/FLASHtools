import yt
import numpy as np
from yt.fields.derived_field import ValidateSpatial
from yt.funcs import just_one

jet_thresh = 1E-3

@yt.derived_field('xray_cooling_time', sampling_type='cell', units='Myr')
def _xray_cooling_time(field, data):
    xray_emis = ('gas', 'xray_emissivity_0.1_100_keV')
    if xray_emis not in data.ds.derived_field_list:
        xray_fields = yt.add_xray_emissivity_field(data.ds, 0.1, 100, table_type='apec', metallicity=0.5)
    gamma = 5./3.
    thermal_energy = 1./(gamma-1)*data['pressure']
    return thermal_energy/data[xray_emis]

@yt.derived_field('total_heating_rate_density', sampling_type='cell', units='erg/s/cm**3')
def _total_heating_rate_density(field, data):
    xray_emis = ('gas', 'xray_emissivity_0.1_100_keV')
    if xray_emis not in data.ds.derived_field_list:
        xray_fields = yt.add_xray_emissivity_field(data.ds, 0.1, 100, table_type='apec', metallicity=0.5)
    return -data['spitzer_heat_flux_divergence'] - data[xray_emis]

@yt.derived_field('spitzer_conduction_coefficient', sampling_type='cell', units='cm**2/s')
def _spitzer_conduction_coefficient(field, data):
    T1 = data['temperature'].to_equivalent('keV', 'thermal')/10/yt.units.keV
    n = data['H_nuclei_density']/1E-3/yt.units.cm**(-3)
    jetregion = data['jet '] > jet_thresh
    T1[jetregion] = T1[jetregion]/1E6
    return 0.1*4E32*T1**2.5/n*yt.units.cm**2*yt.units.s**(-1)

basename = 'spitzer_heat_flux'
xn = basename + '_x'
yn = basename + '_y'
zn = basename + '_z'

def _spitzer_heat_flux(ax):
    def func(field, data):
        return -data['spitzer_conduction_coefficient']*yt.physical_constants.kb*data['H_nuclei_density']*data['temperature_gradient_%s' % ax]
    return func

for ax in 'xyz':
    f = _spitzer_heat_flux(ax)
    yt.add_field('%s_%s' % (basename, ax), function=f, sampling_type='cell', units='erg/s/cm**2')


sl_left = slice(None, -2, None)
sl_right = slice(2, None, None)

div_fac = 2


@yt.derived_field(('gas', "%s_divergence" % basename), sampling_type="cell", units='erg/s/cm**3', validators=[ValidateSpatial(2)])
def _divergence(field, data):
    ds = div_fac * just_one(data["index", "dx"])
    f  = data[xn][sl_right,1:-1,1:-1]/ds
    f -= data[xn][sl_left ,1:-1,1:-1]/ds
    ds = div_fac * just_one(data["index", "dy"])
    f += data[yn][1:-1,sl_right,1:-1]/ds
    f -= data[yn][1:-1,sl_left ,1:-1]/ds
    ds = div_fac * just_one(data["index", "dz"])
    f += data[zn][1:-1,1:-1,sl_right]/ds
    f -= data[zn][1:-1,1:-1,sl_left ]/ds
    new_field = data.ds.arr(np.zeros(data[xn].shape, dtype=np.float64),
                            f.units)
    new_field[1:-1,1:-1,1:-1] = f
    return new_field



@yt.derived_field(('gas', 'spitzer_heating_rate'), sampling_type="cell", units='erg/s')
def _spitzer_heating_rate(field, data):
    return -data[('gas', 'spitzer_heat_flux_divergence')]*data[('gas', 'cell_volume')]


@yt.derived_field(('gas', 'total_heating_rate'), sampling_type="cell", units='erg/s')
def _total_heating_rate(field, data):
    return -data[('gas', 'total_heating_rate_density')]*data[('gas', 'cell_volume')]


@yt.derived_field(('gas', 'total_cooling_time'), sampling_type="cell", units='Myr')
def _total_cooling_time(field, data):
    gamma = 5./3.
    thermal_energy = 1./(gamma-1)*data['pressure']
    return thermal_energy/(-data['total_heating_rate_density'])

