#!/usr/bin/env python
import pdb
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yt
from yt_synchrotron_emissivity import *
yt.enable_parallelism()
import logging
from yt.utilities.fits_image import FITSProjection
logging.getLogger('yt').setLevel(logging.INFO)

#nu = yt.YTQuantity(500, 'MHz')


dir = '/home/ychen/data/0only_1110_h0_rerun/'
#dir = '/home/ychen/data/0only_0605_hinf/'
#dir = '/home/ychen/data/0only_1022_h1_10Myr/'
#dir = '/d/d8/ychen/MHD_Jet/0314_L45_M10_b1_h1_nojiggle'
ts = yt.DatasetSeries(os.path.join(dir,'*_hdf5_plt_cnt_1[0-3]??'), parallel=20)
#ts = yt.DatasetSeries(os.path.join(dir,'*_hdf5_plt_cnt_1410'), parallel=1)
#ts = yt.DatasetSeries(os.path.join(dir,'*_hdf5_plt_cnt_0600'), parallel=8)

zoom_fac = 8

ptype = 'lobe'
maindir = os.path.join(dir, 'pol_synchrotron_QU_nn_%s/' % ptype)
lowres = os.path.join(maindir, 'lowres')
spectral_index_dir = os.path.join(maindir, 'spectral_index')
#emisdir = os.path.join(maindir, 'emissivity')
if yt.is_root():
#    for subdir in [maindir, spectral_index_dir, emisdir]:
#    for subdir in [maindir, lowres, spectral_index_dir]:
    for subdir in [maindir, lowres]:
        if not os.path.exists(subdir):
            os.mkdir(subdir)

for ds in ts.piter():
    #proj = yt.ProjectionPlot(ds, proj_axis, ('gas', 'density'), center=(0,0,0))
    ##proj.set_zlim(('gas', 'density'), 1E-5, 1E-2)
    ##proj.set_cmap(('flash', 'jet '), 'gist_heat')
    #proj.annotate_timestamp(corner='upper_left', time_format="{time:6.3f} {units}", time_unit='Myr', draw_inset_box=True)
    #proj.annotate_text((0.85, 0.95), dir.split('b1_')[-1].strip('/'), coord_system='axis', text_args={'color':'k'})
    #proj.zoom(10)
    #savefn = 'Projection_%s_density_%s.png' % (proj_axis, str(ds).split('_')[-1])
    #proj.save(os.path.join(figuredir,savefn))

    projs = {}
    proj_axis = 'x'
    for nu in [(150, 'MHz')]:
        norm = yt.YTQuantity(*nu).in_units('GHz').value**0.5
#        pars = add_synchrotron_emissivity(ds, ptype=ptype, nu=nu)
#        field = ('deposit', ('nn_emissivity_%s_%%.1f%%s' % ptype) % nu)
#        #field = ('deposit', 'jetp_nn_sync_spec_%.1f%s' % nu)
#        #field = [('deposit', 'avgfill_emissivity_%.1f%s' % nu),\
#        #         ('deposit', 'avgpf_emissivity_%.1f%s' % nu)]
#        #projs[nu] = ds.proj(field, proj_axis, center=[0,0,0])
#
#
#        ###########################################################################
#        ## Slices
#        ###########################################################################
#        #center = yt.YTArray([0.25,0,20], input_units='kpc')
#        #plot = yt.SlicePlot(ds, proj_axis, field, center=center, width=(10,'kpc'))
#        #plot.set_zlim(field, 1E-25/norm, 1E-21/norm)
#        #plot.set_zlim(field, 1E-20/norm, 1E-16/norm)
#        #plot.annotate_particles((0.5,'kpc'), p_size=4, marker='.', ptype='jetp')
#
#        ###########################################################################
#        ## Projection
#        ###########################################################################
#        plot = yt.ProjectionPlot(ds, proj_axis, field, center=[0,0,0])
#        plot.set_zlim(field, 1E-3/norm, 1E1/norm)
#
#        ###########################################################################
#        ## Plot Settings
#        ###########################################################################
#        #plot.set_cmap(field, 'kamae_r')
#        cmap = plt.cm.hot
#        cmap.set_bad('k')
#        #cmap = plt.cm.get_cmap("algae")
#        #cmap.set_bad((80./256., 0.0, 80./256.))
#        plot.set_cmap(field, cmap)
#        plot.annotate_timestamp(corner='upper_left', time_format="{time:6.3f} {units}",
#                                time_unit='Myr', draw_inset_box=True)
#        dirnamesplit = dir.split('_')
#        if dirnamesplit[-1] in ['h1','hinf', 'h0']:
#            x = 0.85
#            sim_name = dirnamesplit[-1]
#        else:
#            x = 0.80
#            sim_name = dirnamesplit[-2] + '_' + dirnamesplit[-1]
#        plot.annotate_text((x, 0.95), sim_name, coord_system='axis',
#                            text_args={'color':'grey'})
#        ##plot.annotate_grids()
#        plot.zoom(10)
#        plot.save(emisdir)
#        #projs[nu] = plot.data_source

        ###########################################################################
        ## Polarizations
        ###########################################################################

        pars = add_synchrotron_pol_emissivity(ds, ptype=ptype, nu=nu, proj_axis=proj_axis)
        fields = []
        for pol in ['i', 'q', 'u']:
            #figuredir = os.path.join(maindir, 'emissivity_%s' % pol)
            #if yt.is_root():
            #    if not os.path.exists(figuredir):
            #        os.mkdir(figuredir)
            field = ('deposit', ('nn_emissivity_%s_%s_%%.1f%%s' % (pol, ptype)) % nu)
            fields.append(field)


        width = ds.domain_width[1:]/zoom_fac
        res = ds.domain_dimensions[1:]*ds.refine_by**ds.index.max_level/zoom_fac/2
        plot = yt.ProjectionPlot(ds, proj_axis, fields, center=[0,0,0], width=width)
        plot.set_buff_size(res)
        plot.set_axes_unit('kpc')
        frb_I = plot.frb.data[fields[0]].v
        frb_Q = plot.frb.data[fields[1]].v
        frb_U = plot.frb.data[fields[2]].v
        #plot.save(maindir)

        for field in fields:
            #print plot.get_log(field)
            if 'nn_emissivity_i' in field[1]:
                plot.set_zlim(field, 1E-3/norm, 1E1/norm)
                #cmap = plt.cm.get_cmap("algae")
                #cmap.set_bad((80./256., 0.0, 80./256.))
                cmap = plt.cm.hot
                cmap.set_bad('k')
                plot.set_cmap(field, cmap)

            else:
                cmap = plt.cm.seismic
                #cmap.set_bad('k')
                plot.set_cmap(field, cmap)
                plot.set_zlim(field, -1E0/norm, 1E0/norm)

        plot.annotate_timestamp(corner='upper_left', time_format="{time:6.3f} {units}",
                                time_unit='Myr', draw_inset_box=True)
        dirnamesplit = dir.split('_')
        if dirnamesplit[-1] in ['h1','hinf', 'h0']:
            x = 0.85
            sim_name = dirnamesplit[-1]
        else:
            x = 0.80
            sim_name = dirnamesplit[-2] + '_' + dirnamesplit[-1]
        plot.annotate_text((x, 0.95), sim_name, coord_system='axis',
                            text_args={'color':'grey'})

        ##plot.annotate_grids()

        # Binning pixels for annotating polarization lines
        factor = 16
        # Saving annotated polarization lines images
        plot.annotate_polline(frb_I, frb_Q, frb_U, factor=factor)
        plot.save(maindir)
        projs[nu] = plot.data_source

        # Saving low resolution image
        plot.set_buff_size(res/factor)
        plot._recreate_frb()
        frb_I = plot.frb.data[fields[0]].v
        frb_Q = plot.frb.data[fields[1]].v
        frb_U = plot.frb.data[fields[2]].v
        plot.annotate_clear(index=-1)
        plot.annotate_polline(frb_I, frb_Q, frb_U, factor=1)
        plot.save(lowres)

    #if yt.is_root():
        #pickle.dump(projs, open(dir+'projs/%s_projs.pickle' % ds.basename, 'wb'))
        #projs[(1.4, 'GHz')].save_object('proj_1.4GHz', dir+'projs/'+ds.basename+'_projs.cpkl')
        #projs[(150, 'MHz')].save_object('proj_150MHz', dir+'projs/'+ds.basename+'_projs.cpkl')

        #ext = ds.arr([-7.72E22, 7.72E22, -1.544E23, 1.544E23], input_units='code_length')
        #frb1 = projs[(1.4, 'GHz')].to_frb(ext[1]-ext[0], (512,1024), height=(ext[3]-ext[2]))
        #frb2 = projs[(150, 'MHz')].to_frb(ext[1]-ext[0], (512,1024), height=(ext[3]-ext[2]))
        #S1 = frb1[('deposit', 'nn_emissivity_%s_1.4GHz' % ptype)]
        #S2 = frb2[('deposit', 'nn_emissivity_%s_150.0MHz' % ptype)]
        #alpha = np.log(S1/S2)/np.log(1400/150)
        #alpha = np.ma.masked_where(S1<0.001, np.array(alpha))

        #fig = plt.figure(figsize=(8,12), dpi=150)
        #ims = plt.imshow(alpha, vmin=-2, vmax=-0.5, extent=ext.in_units('kpc'), origin='lower', aspect='equal')
        #plt.xlabel('z (kpc)')
        #plt.ylabel('x (kpc)')
        #cb = plt.colorbar()
        #cb.set_label('Spectral Index (%s) (1.4GHz/150MHz)' % ptype)

        #dirnamesplit = dir.split('_')
        #if dirnamesplit[-1] in ['h1','hinf', 'h0']:
        #    sim_name = dirnamesplit[-1]
        #else:
        #    sim_name = dirnamesplit[-2] + '_' + dirnamesplit[-1]
        #x = 0.80
        #plt.annotate(sim_name, (1,1), xytext=(x, 0.96),  textcoords='axes fraction',\
        #            horizontalalignment='left', verticalalignment='center')

        #plt.annotate('%6.3f Myr' % (float(ds.current_time)/3.15569E13),\
        #            (0,1), xytext=(0.05, 0.96),  textcoords='axes fraction',\
        #            horizontalalignment='left', verticalalignment='center')
        #plt.tight_layout()
        #plt.savefig(spectral_index_dir+'/%s_proj_spectral_index.png' % ds.basename)

    # Done with this dataset
    #del projs
