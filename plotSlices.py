import os
import numpy as np
import yt
#yt.enable_parallelism()
from tools import read_par, calcNozzleCoords

mu = 1.67E-24
k = 1.38E-16
T = 1.0E7
v = 0.1*3.0E10
mag = 1.0E-5

Extrema = { 'density': (1.0E-4*mu, 1.0E-1*mu), 'pressure':(1.0E-12, 1.0E-9),'temperature': (1.0E0*T, 1.0E3*T),\
            'magnetic_field_x': (-mag, mag), 'magnetic_field_y': (-mag, mag), 'magnetic_field_z':(-mag,mag),\
            'velocity_x':(-0.3333*v,0.3333*v), 'velocity_y': (-0.3333*v,0.3333*v), 'velocity_z':(-v,v),\
            'ism ': (0.0, 1.0), 'jet ': (0.0, 1.0), 'magnetic_pressure': (1.0E-3*mag*mag, 1.0E-1*mag*mag),\
            'dens': (1.0E-5*mu, 1.0E1*mu), 'pres':(1.0E-5*k*T, 1.0E0*k*T),'temp': (10.0*T, 1.0E4*T),\
            'magx': (-mag, mag), 'magy': (-mag, mag), 'magz':(-mag,mag),\
            'velx': (-v,v), 'vely': (-v,v), 'velz':(-v,v),\
            'magp': (1E-20, 1.0E-7), 'eint': (1.0E15, 1.0E18)\
}

logfield = { 'dens': True, 'pres': True, 'temp': True,\
             'velx': False, 'vely': False, 'velz': False,\
             'magx': False, 'magy': False, 'magz': False,\
             'ism ': False, 'jet ': False, 'magp': True,\
             'eint': True,\
             'density': True, 'pressure': True, 'temperature': True,\
             'velocity_x': False, 'velocity_y': False, 'velocity_z': False,\
             'magnetic_field_x': False, 'magnetic_field_y': False, 'magnetic_field_z': False,\
             'magnetic_pressure': True}

fields_ = ['density', 'pressure', 'temperature', 'velocity_y', 'velocity_z', 'jet ',\
           'magnetic_field_x', 'magnetic_field_y', 'magnetic_field_z', 'magnetic_pressure']

axis = {'x': 0, 'y':1, 'z':2}

def plotSliceField(ds, proj_axis='x', field='density', center=(0.0,0.0,0.0),\
                zoom_fac=1, nozzleCoords=None, plotgrid=True, plotvelocity=False,\
                markcenter=False, savepath=None,\
                show=False, annotate_particles=None, annotate_part_info=None):
    plot = yt.SlicePlot(ds, proj_axis, field, center=center)#, width=(0.25,0.5))
    if plotvelocity:
        scale=9E10
        if type(plotvelocity) is list:
            if field in plotvelocity:
                plot.annotate_velocity(10, scale=scale)
        else:
            plot.annotate_velocity(10, scale=scale)
    #plot.set_figure_size(20)
    #plot.set_font({'size': 36})
    if proj_axis == 'z':
        plot.set_buff_size(256)
    else:
        plot.set_buff_size(256)

    if field in ['ism ', 'jet ']:
        plot.set_cmap(field, 'gist_heat')
    if field in ['magnetic_field_x', 'magnetic_field_y', 'magnetic_field_z',\
                'velocity_x', 'velocity_y', 'velocity_z']:
        plot.set_cmap(field, 'RdBu_r')
    plot.set_log(field, logfield[field])
    mi, ma = Extrema[field]
    plot.set_zlim(field, mi, ma)
    plot.zoom(zoom_fac)
    plot.annotate_timestamp(0.15, 0.95, normalized=True, format="{time:6.3f} {units}", fontsize=48)
    if annotate_particles:
        #slab_width = ds.domain_width.value[axis[proj_axis]] - 2.0*center[axis[proj_axis]]
        slab_width = 3E21
        if type(annotate_particles) is list:
            if field in annotate_particles:
                plot.annotate_particles(slab_width)
        else:
            plot.annotate_particles(slab_width)

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
        if type(plotgrid) is list:
            if field in plotgrid:
                plot.annotate_grids(edgecolors="grey")
        else:
            plot.annotate_grids(edgecolors="grey")


    if markcenter:
        plot.annotate_marker(center, marker='s', plot_args={'color':'none',\
                             'edgecolor':'grey', 's':zoom_fac*8})
    if nozzleCoords:
        coordsets = np.array([coord for coord, size in nozzleCoords])
        sizes = (1.0+np.array([size for coord, size in nozzleCoords]))*zoom_fac*4
        plot.annotate_marker((coordsets[:,0], coordsets[:,1]), marker='o', plot_args={'color':'white',\
                             'edgecolor':'black', 's': sizes})


        #for pos, size in nozzleCoords:
        #    plot.annotate_marker(pos, marker='o', plot_args={'color':'white',\
        #                         'edgecolor':'black', 's':(1.0+size)*zoom_fac*4})
    if savepath:
        plot.save(savepath) #, mpl_kwargs={'dpi':300})
    else:
        plot.save(os.getcwd())
    if show:
        plot.show()

    return plot


def plotSlices(ds, proj_axes=['x'], fields=fields_, drawnozzle=True, **kwargs):
    for proj_axis in proj_axes:
        if drawnozzle:
            nozzleCoords = calcNozzleCoords(read_par(os.getcwd()), ds.current_time.value,\
                                            proj_axis=proj_axis)
        else:
            nozzleCoords = None
        for field in fields:
            plotSliceField(ds, proj_axis=proj_axis, field=field, nozzleCoords=nozzleCoords, **kwargs)
