#!/usr/bin/env python
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'stixgeneral'
import matplotlib.pyplot as plt
import yt
yt.enable_parallelism()
yt.mylog.setLevel('INFO')
from synchrotron.yt_synchrotron_emissivity import\
        setup_part_file,\
        write_synchrotron_hdf5,\
        synchrotron_filename,\
        StokesFieldName
from scipy.ndimage import gaussian_filter

#dir = '/home/ychen/data/00only_0529_h1'
dir = '/d/d5/ychen/2015_production_runs/1022_h1_10Myr/'
#dir = '/home/ychen/data/00only_0605_hinf/'
#dir = '/home/ychen/data/00only_0605_h0/'
try:
    ind = int(sys.argv[1])
    ts = yt.DatasetSeries(os.path.join(dir,'data/*_hdf5_plt_cnt_%03d0' % ind), parallel=1, setup_function=setup_part_file)
except IndexError:
    ts = yt.DatasetSeries(os.path.join(dir,'data/*_hdf5_plt_cnt_???0'), parallel=2, setup_function=setup_part_file)

zoom_fac = 8
format = 'png'

#proj_axis = [1,0,2]
proj_axis = 'x'
#nus = [(150, 'MHz'), (233, 'MHz'), (325, 'MHz'), (610, 'MHz'), (1400, 'MHz')]
#nus = [(325, 'MHz'), (610, 'MHz'), (1400, 'MHz')]:
#nus = [(150, 'MHz'), (1400, 'MHz')]
nus = [(100, 'MHz')]
ptype = 'lobe'
gc = 32
maindir = os.path.join(dir, 'cos_synchrotron_QU_nn_%s/' % ptype)
if proj_axis != 'x':
    maindir = os.path.join(maindir, '%i_%i_%i' % tuple(proj_axis))
polline = os.path.join(maindir, 'polline')
spectral_index_dir = os.path.join(maindir, 'spectral_index')
if yt.is_root():
    for subdir in [maindir, polline, spectral_index_dir]:
        if not os.path.exists(subdir):
            os.mkdir(subdir)

for ds in ts.piter():
    if '0000' in ds.basename: continue

    projs, frb_I = {}, {}
    for nu in nus:
        norm = yt.YTQuantity(*nu).in_units('GHz').value**0.5

        width = ds.domain_width[1:]/zoom_fac
        #res = ds.domain_dimensions[1:]*ds.refine_by**ds.index.max_level//zoom_fac//2
        res = np.array([1200, 2400])
        stokes = StokesFieldName(ptype, nu, proj_axis, field_type='flash')

        write_synchrotron_hdf5(ds, ptype, nu, proj_axis, extend_cells=32)
        sync_fname = synchrotron_filename(ds, extend_cells=gc)
        if os.path.isfile(sync_fname):
            ds_sync = yt.load(sync_fname)
            # Need to build the field list for the field_into
            ds_sync.field_list
            for field in stokes.IQU:
                ds_sync.field_info[field].units = 'Jy/cm/arcsec**2'
                ds_sync.field_info[field].output_units = 'Jy/cm/arcsec**2'
        else:
        # We will calculate all things on the fly
            pars = add_synchrotron_dtau_emissivity(ds, ptype=ptype, nu=nu, proj_axis=proj_axis, extend_cells=gc)
            ds_sync = ds

        yt.mylog.info('Making projection plots from file: %s', sync_fname)
        if proj_axis == 'x':
            plot = yt.ProjectionPlot(ds_sync, proj_axis, stokes.IQU, center=[0,0,0], width=width)
        else:
            plot = yt.OffAxisProjectionPlot(ds_sync, proj_axis, stokes.IQU, center=[0,0,0], width=width,
                                            north_vector=[0,0,1])
        plot.set_buff_size(res)
        plot.set_axes_unit('kpc')

        # Setting up colormaps
        # Use "hot" for intensity plot and seismic for Q and U plots
        for field in stokes.IQU:
            if 'nn_emissivity_i' in field[1]:
                plot.set_zlim(field, 1E-3/norm, 1E1/norm)
                #cmap = plt.cm.get_cmap("algae")
                #cmap.set_bad((80./256., 0.0, 80./256.))
                cmap = plt.cm.hot
                cmap.set_bad('k')
                plot.set_cmap(field, cmap)

            elif 'nn_emissivity' in field[1]:
                cmap = plt.cm.seismic
                plot.set_cmap(field, cmap)
                plot.set_zlim(field, -1E0/norm, 1E0/norm)

        plot.annotate_timestamp(corner='upper_left', time_format="{time:6.3f} {units}",
                                time_unit='Myr', draw_inset_box=True)
        dirnamesplit = dir.strip('/').split('_')
        if dirnamesplit[-1] in ['h1','hinf', 'h0']:
            if dirnamesplit[-1] == 'h1':
                x = 0.80
                sim_name = 'helical'
            else:
                x = 0.85
                sim_name = dirnamesplit[-1]
        else:
            x = 0.80
            sim_name = dirnamesplit[-2] + '_' + dirnamesplit[-1]
        plot.annotate_text((x, 0.95), sim_name, coord_system='axis',
                            text_args={'color':'grey'})

        #plot.annotate_grids()

        # Saving annotated polarization lines images
        #if pol in ['q', 'u']:
        #    plot.annotate_polline(frb_I, frb_Q, frb_U, factor=factor)
        if yt.is_root():
            plot.save(maindir, suffix=format)
            projs[nu] = plot.data_source

            ###########################################################################
            ## Annotating polarization lines
            ###########################################################################

            # Binning pixels for annotating polarization lines
            factor = 1
            # Saving polarization lines image
            plot.set_buff_size(res//factor)
            plot._recreate_frb()
            frb_I[nu] = plot.frb.data[stokes.I].v
            frb_Q = plot.frb.data[stokes.Q].v
            frb_U = plot.frb.data[stokes.U].v
            #plot.annotate_clear(index=-1)
            plot.annotate_polline(frb_I[nu], frb_Q, frb_U, factor=25, scale=15)
            plot.save(polline, suffix=format)
        ds_sync.close()

    if yt.is_root():
        sigma = 1
        nu1, nu2 = nus[0], nus[-1]
        I1 = gaussian_filter(frb_I[nu1], sigma)
        I2 = gaussian_filter(frb_I[nu2], sigma)
        alpha = np.log10(I2/I1)/np.log10(1400/150)
        alpha = np.ma.masked_where(I2<1E-3, np.array(alpha))
        ext = ds.arr([-0.5*width[0], 0.5*width[0], -0.5*width[1], 0.5*width[1]])
        plt.figure(figsize=(8,12), dpi=150)
        cmap = plt.cm.jet
        cmap.set_bad('navy')
        plt.imshow(alpha, cmap=cmap, vmin=-1.5, vmax=-0.5, extent=ext.in_units('kpc'), origin='lower', aspect='equal')
        plt.xlabel('y (kpc)')
        plt.ylabel('z (kpc)')
        plt.axes().tick_params(direction='in')
        cb = plt.colorbar(fraction=0.10, pad=0, aspect=50)
        cb.set_label('Spectral Index (%s) (1.4GHz/150MHz)' % ptype)
        cb.ax.tick_params(direction='in')
        #pickle.dump(projs, open(dir+'projs/%s_projs.pickle' % ds.basename, 'wb'))
        #projs[(1.4, 'GHz')].save_object('proj_1.4GHz', dir+'projs/'+ds.basename+'_projs.cpkl')
        #projs[(150, 'MHz')].save_object('proj_150MHz', dir+'projs/'+ds.basename+'_projs.cpkl')

        dirnamesplit = dir.strip('/').split('_')
        if dirnamesplit[-1] in ['h1','hinf', 'h0']:
            if dirnamesplit[-1] == 'h1':
                x = 0.80
                sim_name = 'helical'
            else:
                x = 0.85
                sim_name = dirnamesplit[-1]
        else:
            x = 0.80
            sim_name = dirnamesplit[-2] + '_' + dirnamesplit[-1]
        plt.annotate(sim_name, (1,1), xytext=(x, 0.96),  textcoords='axes fraction',\
                    horizontalalignment='left', verticalalignment='center')

        plt.annotate('%6.3f Myr' % (float(ds.current_time)/3.15569E13),\
                    (0,1), xytext=(0.05, 0.96),  textcoords='axes fraction',\
                    horizontalalignment='left', verticalalignment='center')
        plt.tight_layout()
        sifname = ds.basename
        if proj_axis == 'x':
            sifname += '_x_'
        else:
            sifname += '%i_%i_%i_' % tuple(proj_axis)
        sifname += 'proj_spectral_index'
        plt.savefig(os.path.join(spectral_index_dir, sifname), format=format)

    # Done with this dataset
    yt.mylog.info('Closing file: %s', ds)
    ds.close()
    del projs
