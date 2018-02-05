import os
import numpy as np
import yt
from yt.fields.field_detector import FieldDetector
#yt.enable_parallelism()
from tools import calcNozzleCoords

mu = 1.67E-24
k = 1.38E-16
T = 1.0E7
v = 0.1*3.0E10
mag = 4.0E-5

kpc = yt.units.kpc.in_units("cm")

Extrema = { 'density': (1.0E-4*mu, 1.0E-1*mu), 'pressure':(1.0E-12, 1.0E-8),\
            'temperature': (3.0E0*T, 6.0E0*T), 'entropy': (20, 60),\
            'temperature_ratio': (0.5, 2), 'entropy_ratio': (0.5, 2),\
            'magnetic_field_x': (-mag, mag), 'magnetic_field_y': (-mag, mag), 'magnetic_field_z':(-mag,mag),\
            'velocity_x':(-0.3333*v,0.3333*v), 'velocity_y': (-0.3333*v,0.3333*v), 'velocity_z':(-v,v),\
            #'velocity_x':(-0.1*v,0.1*v), 'velocity_y': (-1.0E7,1.0E7), 'velocity_z':(-1.0E5,1.0e5),\
            'velocity_magnitude':(1.0E3, 1.0E7),\
            'velocity_para':(-v,v), 'velocity_perp':(-0.05*v, 0.05*v), 'mach':(0.0, 30.0),\
            'plasma_beta':(1E-1, 1E3),\
            'shok': (0.0, 10.0),\
            'ism ': (0.0, 1.0), 'jet ': (0.0, 1.0), 'magnetic_pressure': (1.0E-4*mag*mag, 1.0E-1*mag*mag),\
            'dens': (1.0E-5*mu, 1.0E1*mu), 'pres':(1.0E-5*k*T, 1.0E0*k*T),'temp': (10.0*T, 1.0E4*T),\
            'magx': (-mag, mag), 'magy': (-mag, mag), 'magz':(-mag,mag),\
            'velx': (-v,v), 'vely': (-v,v), 'velz':(-v,v),\
            'magp': (1E-20, 1.0E-7), 'eint': (1.0E15, 1.0E18)\
}

logfield = { 'dens': True, 'pres': True, 'temp': True,\
             'velx': False, 'vely': False, 'velz': False,\
             'magx': False, 'magy': False, 'magz': False,\
             'ism ': False, 'jet ': False, 'magp': True,\
             'eint': True, 'shok': False,\
             'density': True, 'pressure': True, \
             'temperature': True, 'entropy': False,\
             'temperature_ratio': True, 'entropy_ratio': True,\
             'velocity_x': False, 'velocity_y': False, 'velocity_z': False,\
             'velocity_para':False, 'velocity_perp':False, 'mach': False,\
             'velocity_magnitude':True,\
             'magnetic_field_x': False, 'magnetic_field_y': False, 'magnetic_field_z': False,\
             'magnetic_pressure': True,\
             'plasma_beta': True}

fields_ = ['density', 'pressure', 'temperature', 'velocity_y', 'velocity_z', 'jet ',\
           'magnetic_field_x', 'magnetic_field_y', 'magnetic_field_z', 'magnetic_pressure']

axis = {'x': 0, 'y':1, 'z':2}

def _temperature_ratio(field, data):
    from yt.units import cm, Kelvin
    Tout    = data.ds.parameters['sim_tout']*Kelvin
    Tcore   = data.ds.parameters['sim_tcore']*Kelvin
    rCoreT  = data.ds.parameters['sim_rcoret']*cm

    if not isinstance(data, FieldDetector):
        data.set_field_parameter('center', (0,0,0))
    r = data['index', 'spherical_radius']
    T0 = Tout*(1.0+(r/rCoreT)**3)/(Tout/Tcore+(r/rCoreT)**3)

    return data['temperature']/T0


def _entropy_ratio(field, data):
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
    mw = data.get_field_parameter("mu")
    if mw is None: mw = 1.0
    mw *= mH
    density0 = rhoCore*(1.0 + (r/rCore)**2)**(-1.5*densitybeta)
    T0 = Tout*(1.0+(r/rCoreT)**3)/(Tout/Tcore+(r/rCoreT)**3)

    icm_entropy = k*T0*(density0/mw)**(-2./3.)

    return data['entropy'].in_units('keV*cm**2')/icm_entropy.in_units('keV*cm**2')

def plotSliceField(ds, proj_axis='x', field='density', center=(0.0,0.0,0.0),\
                zoom_fac=1, nozzleCoords=None, plotgrid=True, plotvelocity=False,\
                markcenter=False, savepath=None, sim_name=None, width=None,\
                show=False, annotate_particles=None, annotate_part_info=None,\
                buff_size=(400,800)):
    if field == 'entropy_ratio':
        ds.add_field('entropy_ratio', function=_entropy_ratio,
                display_name='Entropy Ratio', sampling_type='cell')
    if field == 'temperature_ratio':
        ds.add_field('temperature_ratio', function=_temperature_ratio,
                display_name='Temperature Ratio', sampling_type='cell')

    plot = yt.SlicePlot(ds, proj_axis, field, center=center, width=width)
    if plotvelocity:
        scale=10*v
        #scale = 5.0E5
        if type(plotvelocity) is list:
            if field in plotvelocity:
                plot.annotate_velocity(100, scale=scale)
        else:
            plot.annotate_velocity(100, scale=scale)
    plot.set_figure_size(10)
    #plot.set_font({'size': 36})
    plot.set_antialias(False)
    plot.set_buff_size((2**10,2**11))

    if field in ['ism ', 'jet ', 'shok']:
        plot.set_cmap(field, 'gist_heat')
    elif field in ['shok']:
        plot.set_cmap(field, 'gist_heat_r')
    elif field in ['plasma_beta']:
        plot.set_cmap(field, 'algae_r')
    elif field in ['magnetic_field_x', 'magnetic_field_y', 'magnetic_field_z',\
                'velocity_x', 'velocity_y', 'velocity_z', 'velocity_para', 'velocity_perp',\
                'temperature_ratio', 'entropy_ratio']:
        plot.set_cmap(field, 'seismic')
    else:
        plot.set_cmap(field, 'viridis')
    #if field in ['magnetic_pressure'] and proj_axis=='x':
    #    plot.annotate_line_integral_convolution('magnetic_field_y', 'magnetic_field_z', lim=(0.4,0.65), cmap='binary', alpha=0.8)

    plot.set_log(field, logfield[field])
    mi, ma = Extrema[field]
    plot.set_zlim(field, mi, ma)
    plot.zoom(zoom_fac)
    plot.annotate_timestamp(corner='upper_left', time_format="{time:6.3f} {units}",
                            time_unit='Myr', draw_inset_box=True)
    if sim_name:
        plot.annotate_text((0.85, 0.90), sim_name, coord_system='axis',
                           text_args={'color':'k'})
    if annotate_particles:
        try:
            #slab_width = ds.domain_width.value[axis[proj_axis]] - 2.0*center[axis[proj_axis]]
            slab_width = 0.5*kpc if field in ['shok'] else 2.0*kpc
            if type(annotate_particles) is list:
                if field in annotate_particles:
                    plot.annotate_particles(slab_width)
            else:
                plot.annotate_particles(slab_width)
        except:
            print('Cannot plot particles:', ds.basename)

    if annotate_part_info:
        ad = ds.all_data()
        part_tag = 7.0
        if part_tag in ad['particle_tag']:
            part_ind = abs(ad['particle_tag']-part_tag).argmin()
            textfmt = r'$\gamma_c = %.3e$'+'\n'+r'$\rho = %.3e$'+'\n'+'$|B| = %.3e$'
            B = np.sqrt(ad['particle_magx'][part_ind]**2\
                     + ad['particle_magy'][part_ind]**2\
                     + ad['particle_magz'][part_ind]**2)*np .sqrt(4.0*np.pi)
            text = textfmt % (ad['particle_gamc'][part_ind], ad['particle_dens'][part_ind], B)
            text_args = {'color': 'black', 'size':16, 'ha':'right'}

            marker_args = {'c':'white', 's':20}
            pos = (ad['particle_position_y'][part_ind], ad['particle_position_z'][part_ind])
            plot.annotate_marker(pos, marker='o', plot_args=marker_args)
            plot.annotate_text((0.95,0.89), text, text_args=text_args)
        part_tag = 46.0
        if part_tag in ad['particle_tag']:
            part_ind = abs(ad['particle_tag']-part_tag).argmin()
            textfmt = r'$\gamma_c = %.3e$'+'\n'+r'$\rho = %.3e$'+'\n'+'$|B| = %.3e$'
            B = np.sqrt(ad['particle_magx'][part_ind]**2\
                     + ad['particle_magy'][part_ind]**2\
                     + ad['particle_magz'][part_ind]**2)*np .sqrt(4.0*np.pi)
            text = textfmt % (ad['particle_gamc'][part_ind], ad['particle_dens'][part_ind], B)
            text_args = {'color': 'black', 'size':16, 'ha':'right'}

            marker_args = {'c':'cyan', 's':20}
            pos = (ad['particle_position_y'][part_ind], ad['particle_position_z'][part_ind])
            plot.annotate_marker(pos, marker='o', plot_args=marker_args)
            plot.annotate_text((0.45,0.02), text, text_args=text_args)
    if plotgrid:
        try:
            if type(plotgrid) is list:
                if field in plotgrid:
                    plot.annotate_grids(edgecolors="grey", linewidth=0.5, alpha=0.1)
            else:
                plot.annotate_grids(edgecolors="grey", linewidth=0.5, alpha=0.1)
        except:
            print('Cannot plot grids:', ds.basename)


    if markcenter:
        plot.annotate_marker(center, marker='s', plot_args={'color':'none',\
                             'edgecolor':'grey', 's':zoom_fac*8})
    if nozzleCoords:
        coordsets = np.array([coord for coord, size in nozzleCoords])
        sizes = (1.0+np.array([size for coord, size in nozzleCoords]))*zoom_fac*1
        plot.annotate_marker((coordsets[:,0], coordsets[:,1]), marker='o', plot_args={'color':'white',\
                             'edgecolor':'black', 's': sizes})


        #for pos, size in nozzleCoords:
        #    plot.annotate_marker(pos, marker='o', plot_args={'color':'white',\
        #                         'edgecolor':'black', 's':(1.0+size)*zoom_fac*4})
    if savepath:
        plot.save(savepath, mpl_kwargs={'dpi':200})
    #else:
    #    plot.save(os.getcwd())
    if show:
        plot.show()

    return plot


def plotSlices(ds, proj_axes=['x'], fields=fields_, drawnozzle=True, **kwargs):
    for proj_axis in proj_axes:
        if drawnozzle:
            nozzleCoords = calcNozzleCoords(ds, proj_axis=proj_axis)
        else:
            nozzleCoords = None
        for field in fields:
            plotSliceField(ds, proj_axis=proj_axis, field=field, nozzleCoords=nozzleCoords, **kwargs)
