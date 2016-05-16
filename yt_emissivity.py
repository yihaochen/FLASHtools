import numpy as np
import yt
from yt.data_objects.particle_filters import add_particle_filter
from yt.data_objects.particle_unions import ParticleUnion
from yt.fields.derived_field import ValidateGridType
from yt.fields.field_detector import FieldDetector
from tools import calcDen0
from functools import partial

deltaden0 = 5E-31

@yt.particle_filter(requires=["particle_shok"], filtered_type="io")
def jetp(pfilter, data):
    den0 = calcDen0(data['io', 'particle_tadd'])
    jet = np.logical_and((data[("io", "particle_shok")] == 0),
                         (data[("io", "particle_gamc")] > 0.0))
    core = np.logical_and((data['io', 'particle_den0']>den0-deltaden0),\
                          (data['io', 'particle_den0']<den0+deltaden0))
    filter = np.logical_and(jet, core)
    return filter

@yt.particle_filter(requires=["particle_shok"], filtered_type="io")
def shok(pfilter, data):
    shok = (data[("io", "particle_shok")] == 1)
    gamc = np.logical_and((data[("io", "particle_gamc")] > 0.0),\
                          (data[("io", "particle_gamc")] < 1E40))
    filter = np.logical_and(shok, gamc)
    return filter

@yt.particle_filter(requires=["particle_shok"], filtered_type="io")
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

@yt.particle_filter(requires=["particle_shok"], filtered_type="io")
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

@yt.particle_filter(requires=["particle_shok"], filtered_type="io")
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



#add_particle_filter("jetp", function=_jetp, filtered_type='io', requires=["particle_shok"])
#add_particle_filter("shok", function=_shok, filtered_type='io', requires=["particle_shok"])


def _jet_volume_fraction(field, data):
    from yt.units import g, cm, Kelvin
    mH = yt.utilities.physical_constants.mass_hydrogen
    k  = yt.utilities.physical_constants.boltzmann_constant

    rhoCore = data.ds.parameters['sim_rhocore']*g/cm**3
    rCore   = data.ds.parameters['sim_rcore']*cm
    densitybeta = data.ds.parameters['sim_densitybeta']
    Tout    = data.ds.parameters['sim_tout']*Kelvin
    Tcore   = data.ds.parameters['sim_tcore']*Kelvin
    rCoreT  = data.ds.parameters['sim_rcoret']*cm
    gammaICM= data.ds.parameters['sim_gammaicm']
    mu      = data.ds.parameters['sim_mu']

    if not isinstance(data, FieldDetector):
        data.set_field_parameter('center', (0,0,0))
    r = data['index', 'spherical_radius']

    density0 = rhoCore*(1.0 + (r/rCore)**2)**(-1.5*densitybeta)
    T0 = Tout*(1.0+(r/rCoreT)**3)/(Tout/Tcore+(r/rCoreT)**3)
    P0 = density0/mu/mH*k*T0
    icm_mass_fraction = 1.0 - data['flash', 'jet ']
    P = data['gas', 'pressure']
    density = data['gas', 'density']

    icm_volume_fraction = (P0/P)**(1/gammaICM)*icm_mass_fraction*density/density0

    icm_volume_fraction = np.where(icm_volume_fraction < 1.0, icm_volume_fraction, 1.0)

    return 1.0 - icm_volume_fraction


def add_emissivity(ds, ptype='jnsp', nu=(1.4, 'GHz'), method='nearest'):
    me = yt.utilities.physical_constants.mass_electron #9.109E-28
    c  = yt.utilities.physical_constants.speed_of_light #2.998E10
    e  = yt.utilities.physical_constants.elementary_charge #4.803E-10 esu

    gamma_min = yt.YTQuantity(10, 'dimensionless')

    nu = yt.YTQuantity(*nu)
    nu_str = str(nu).replace(' ', '')
    fnames = []
    if ('gas', 'jet_volume_fraction') not in ds.derived_field_list:
        ds.add_field(('gas', 'jet_volume_fraction'), function=_jet_volume_fraction,
                     display_name="Jet Volume Fraction")

    def _synchrotron_spec(field, data):
        ptype = 'io'
        # To convert from FLASH "none" unit to cgs unit, times the B field from FLASH by sqrt(4*pi)
        B = np.sqrt(data[(ptype, 'particle_magx')]**2
                   +data[(ptype, 'particle_magy')]**2
                   +data[(ptype, 'particle_magz')]**2)*np.sqrt(4.0*np.pi)
        B = data.apply_units(B, 'gauss')
        # Cutoff frequency
        nuc = 3.0*data[(ptype, 'particle_gamc')]**2*e*B/(4.0*np.pi*me*c)
        #nu = data.get_field_parameter("frequency", default=yt.YTQuantity(1.4, 'GHz'))
        fit_const = 5.8
        norm = 0.5*e**3.5/(c**2.5*me**1.5*(4.*np.pi)**0.5)
        N0 = 3.0/me/c/c/(np.log(np.abs(data[(ptype, 'particle_gamc')]/gamma_min)))

        return N0*norm*fit_const*nu**(-0.5)*np.exp(-nu/nuc)

    fname1 =('io', 'particle_sync_spec_%s' % nu_str)
    ds.add_field(fname1, function=_synchrotron_spec, particle_type=True,
                 units='cm**(3/4)*s**(3/2)/g**(3/4)', force_override=True)
    ds.add_particle_filter('jetp')
    ds.add_particle_filter('shok')
    ds.add_particle_filter('jnsp')
    ds.add_particle_filter('lobe')
    ds.add_particle_filter('lnsp')
    #import pdb; pdb.set_trace()

    ###########################################################################
    ## Nearest Neighbor method
    ###########################################################################
    fname2 = ds.add_deposited_particle_field(
            (ptype, 'particle_sync_spec_%s' % nu_str), 'nearest')

    def _nn_emissivity(field, data):
        '''
        Emissivity using nearest neighbor. Integrate over line of sight to get intensity.
        '''
        PB =  data['gas', 'pressure']\
              *data['gas', 'magnetic_field_strength']**1.5\
              /yt.YTQuantity(4.*np.pi, 'sr')
        frac = data['gas', 'jet_volume_fraction']
        return PB*frac*data['deposit', '%s_nn_sync_spec_%s' % (ptype, nu_str)]


    #print ds.field_info[('jetp', 'particle_emissivity')].units
    #print ds.field_info[f4].units

    fname3 = ('deposit', 'nn_emissivity_%s_%s' % (ptype, nu_str))
    ds.add_field(fname3, function=_nn_emissivity,
                 display_name='%s NN Emissivity (%s)' % (nu_str, ptype),
                 units='Jy/cm/arcsec**2', take_log=True,
                 force_override=True)


    ###########################################################################
    ## Average Filling method
    ###########################################################################
    def _deposit_average_filling(field, grid):
        #data._debug()
        ptype = 'jetp'
        deposit_field = fname1[1]
        PB =  grid['gas', 'pressure']\
              *grid['gas', 'magnetic_field_strength']**1.5\
              /yt.YTQuantity(4.*np.pi, 'sr')
        frac = grid['gas', 'jet_volume_fraction']
        uq = PB.uq*grid[ptype, deposit_field].uq
        if isinstance(grid, FieldDetector): return grid[ptype, deposit_field]*PB
        if len(grid[ptype, deposit_field]) > 0: return PB*grid[ptype, deposit_field].mean()
        elif grid.Level == 0: return grid['zeros']*uq
        else:
            pfield = yt.YTArray([0], input_units=grid[ptype, deposit_field].units)
            for gr in grid.Parent.Children:
                pfield = np.concatenate([pfield, gr[ptype, deposit_field]])
            if len(pfield) == 0:
                return grid['zeros']*uq
            else:
                return PB*pfield.mean()*grid[ptype, deposit_field].uq

    fname4 =('deposit', 'avgfill_emissivity_%s' % nu_str)
    ds.add_field(fname4, function=_deposit_average_filling,
                 validators=[ValidateGridType()],
                 display_name="%s Emissivity" % nu,
                 units='Jy/cm/arcsec**2', take_log=True,
                 force_override=True)

    ###########################################################################
    ## Particle filter method
    ###########################################################################

    def _deposit_avg_pfilter(field, grid):
        ptype = 'jetp'
        p = grid.ds.all_data()[ptype, 'particle_position']
        deposit_field = fname1[1]
        pfield = grid.ds.all_data()[ptype, deposit_field]
        PB =  grid['gas', 'pressure']\
              *grid['gas', 'magnetic_field_strength']**1.5\
              /yt.YTQuantity(4.*np.pi, 'sr')
        if isinstance(grid, FieldDetector): return grid[ptype, deposit_field]*PB
        l = grid.LeftEdge
        r = grid.RightEdge
        c = (l+r)/2.
        filter = (p[:,0]>l[0])*(p[:,0]<r[0])
        for i in range(1,3):
            filter *= (p[:,i]>l[i])*(p[:,i]<r[i])
        #reg = grid.ds.region(c, c+(l-c)*fac, c+(r-c)*fac)
        expand = 0
        while len(pfield[filter]) == 0 and expand < 3:
            if grid.Level == 0: return grid['zeros']*PB.uq*grid[ptype, deposit_field].uq
            expand += 1
            l = c+(l-c)*2.
            r = c+(r-c)*2.
            filter = (p[:,0]>l[0])*(p[:,0]<r[0])
            for i in range(1,3):
                filter *= (p[:,i]>l[i])*(p[:,i]<r[i])
            #reg = grid.ds.region(c, c+(l-c)*fac, c+(r-c)*fac)
        return PB*pfield[filter].mean()

    fname5 = ('deposit', 'avgpf_emissivity_%s' %nu_str)
    ds.add_field(fname5, function=_deposit_avg_pfilter,
                 validators=[ValidateGridType()],
                 display_name='%s Avg (pfilter) Emissivity' % nu,
                 units='Jy/cm/arcsec**2', take_log=True,
                 force_override=True)



    return fname1, fname2, fname3, fname4, fname5, nu_str

