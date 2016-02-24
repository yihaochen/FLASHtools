import numpy as np
import yt
from yt.data_objects.particle_filters import add_particle_filter
from yt.fields.derived_field import ValidateGridType
from yt.fields.field_detector import FieldDetector

def _jetp(pfilter, data):
    #filter = data[("all", "particle_shok")] == 0
    filter = np.logical_and((data[("io", "particle_shok")] == 0), 
                            (data[("io", "particle_gamc")] > 0.0))
    return filter

add_particle_filter("jetp", function=_jetp, filtered_type='io', requires=["particle_shok"])


me = yt.utilities.physical_constants.mass_electron #9.109E-28
c  = yt.utilities.physical_constants.speed_of_light #2.998E10
e  = yt.utilities.physical_constants.elementary_charge #4.803E-10 esu
gamma_min = yt.YTQuantity(10, 'dimensionless')

def _synchrotron_emissivity(field, data):
    '''
    Emissivity per energy. Needs to multiply by pressure to get emissivity per volumn.
    '''
    ptype = 'io'
    # To convert from FLASH "none" unit to cgs unit, times the B field from FLASH by sqrt(4*pi)
    B = np.sqrt(data[(ptype, 'particle_magx')]**2+data[(ptype, 'particle_magy')]**2+data[(ptype, 'particle_magz')]**2)\
        *np.sqrt(4.0*np.pi)
    B = data.apply_units(B, 'G')
    nuc = 3.0*data[(ptype, 'particle_gamc')]**2*e*B/(4.0*np.pi*me*c)
    nu = data.get_field_parameter("frequency", default=yt.YTQuantity(1.4, 'GHz'))
    fit_const = 5.8
    norm = 0.5*B**1.5*e**3.5/(c**2.5*me**1.5*(4.*np.pi)**0.5)
    N0 = 3.0/me/c/c/(np.log(np.abs(data[(ptype, 'particle_gamc')]/gamma_min)))

    return N0*norm*fit_const*nu**(-0.5)*np.exp(-nu/nuc)

def _deposit_average_filling(field, grid):
    #data._debug()
    ptype = 'jetp'
    deposit_field = 'particle_emissivity'
    uq = grid['gas', 'pressure'].uq
    if isinstance(grid, FieldDetector): return grid[ptype, deposit_field]*uq
    if len(grid[ptype, deposit_field]) > 0: return grid['gas', 'pressure']*grid[ptype, deposit_field].mean()
    elif grid.Level == 0: return grid['zeros']*uq
    else:
        pfield = np.zeros(0)
        for gr in grid.Parent.Children:
            pfield = np.concatenate([pfield, gr[ptype, deposit_field]])
        if len(pfield) == 0: 
            return grid['zeros']*uq
        else:
            return grid['gas', 'pressure']*pfield.mean()

def _intensity(field, data):
    '''
    Intensity per length. Integrate over line of sight to get intensity.
    '''
    return data['deposit', 'avgfill_emissivity']/yt.YTQuantity(4.*np.pi, 'sr')


def add_emissivity(ds, nu=yt.YTQuantity(1.4, 'GHz')):
    f1 = ds.add_field(('io', 'particle_emissivity'), function=_synchrotron_emissivity, particle_type=True,
                      display_name="%s Emissivity" % nu, force_override=True)

    filter = ds.add_particle_filter('jetp')

    f2 = ds.add_field(('deposit', 'avgfill_emissivity'), function=_deposit_average_filling, validators=[ValidateGridType()],
                      display_name="%s Emissivity" % nu, units='erg/s/cm**3/Hz', take_log=True, force_override=True)


    f3 = ds.add_field(('deposit', 'avgfill_intensity'), function=_intensity, display_name="%s Intensity" % nu,
                      units='Jy/cm/arcsec**2', take_log=True, force_override=True)

    return filter, f1, f2, f3, nu

