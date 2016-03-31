#!/usr/bin/env python
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yt
from yt_emissivity import *
yt.enable_parallelism()
import logging
logging.getLogger('yt').setLevel(logging.INFO)

#nu = yt.YTQuantity(500, 'MHz')

#dir = '/home/ychen/d9/FLASH4/stampede/0529_L45_M10_b1_h1/'
dir = '/d/d8/ychen/MHD_Jet/0314_L45_M10_b1_h1_nojiggle'
ts = yt.DatasetSeries(os.path.join(dir,'*_hdf5_plt_cnt_0??0'), parallel=20)
#ts = yt.DatasetSeries(os.path.join(dir,'*_hdf5_plt_cnt_1410'), parallel=1)
#ts = yt.DatasetSeries(os.path.join(dir,'*_hdf5_plt_cnt_0600'), parallel=8)

figuredir = os.path.join(dir, 'synchrotron_nn/')
if yt.is_root():
    if not os.path.exists(os.path.join(dir,figuredir)):
        os.mkdir(os.path.join(dir, figuredir))

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
    for nu in [(150, 'MHz'), (1.4, 'GHz')]:
        norm = yt.YTQuantity(*nu).in_units('GHz').value**0.5
        pars = add_emissivity(ds, nu=nu)
        field = ('deposit', 'nn_emissivity_%.1f%s' % nu)
        #field = ('deposit', 'jetp_nn_sync_spec_%.1f%s' % nu)
        #field = [('deposit', 'avgfill_emissivity_%.1f%s' % nu),\
        #         ('deposit', 'avgpf_emissivity_%.1f%s' % nu)]


        ###########################################################################
        ## Slices
        ###########################################################################
        #center = yt.YTArray([0.25,0,20], input_units='kpc')
        #plot = yt.SlicePlot(ds, proj_axis, field, center=center, width=(10,'kpc'))
        #plot.set_zlim(field, 1E-25/norm, 1E-21/norm)
        #plot.set_zlim(field, 1E-20/norm, 1E-16/norm)
        #plot.annotate_particles((0.5,'kpc'), p_size=4, marker='.', ptype='jetp')

        ###########################################################################
        ## Projection
        ###########################################################################

        plot = yt.ProjectionPlot(ds, proj_axis, field, center=[0,0,0])
        plot.set_zlim(field, 1E-3/norm, 1E1/norm)

        ###########################################################################
        ## Plot Settings
        ###########################################################################
        #plot.set_cmap(field, 'kamae_r')
        plot.annotate_timestamp(corner='upper_left', time_format="{time:6.3f} {units}",
                                time_unit='Myr', draw_inset_box=True)
        plot.annotate_text((0.85, 0.95), dir.split('b1_')[-1].strip('/'), coord_system='axis',
                            text_args={'color':'k'})
        #plot.annotate_grids()
        plot.zoom(10)
        plot.save(figuredir)

    #if yt.is_root():
        #pickle.dump(projs, open(dir+'projs/%s_projs.pickle' % ds.basename, 'wb'))
        #projs[(1.4, 'GHz')].save_object('proj_1.4GHz', dir+'projs/'+ds.basename+'_projs.cpkl')
        #projs[(150, 'MHz')].save_object('proj_150MHz', dir+'projs/'+ds.basename+'_projs.cpkl')
        #frb1 = projs[(1.4, 'GHz')].to_frb((1.544E23, 'code_length'), (1600,800), height=(3.088E23, 'code_length'))
        #frb2 = projs[(150, 'MHz')].to_frb((1.544E23, 'code_length'), (1600,800), height=(3.088E23, 'code_length'))
        #S1 = frb1[('deposit', 'nn_emissivity_1.4GHz')]
        #S2 = frb2[('deposit', 'nn_emissivity_150.0MHz')]
        #alpha = np.log(S1/S2)/np.log(1400/150)

        #fig = plt.figure(figsize=(8,12), dpi=150)
        #ext = np.array([-7.72E22, 7.72E22, -1.544E23, 1.544E23])/3.086E21
        #ims = plt.imshow(np.array(alpha), vmin=-2, vmax=-0.5, extent=ext, origin='lower', aspect='equal')
        #plt.xlabel('y (kpc)')
        #plt.ylabel('z (kpc)')
        #cb = plt.colorbar()
        #cb.set_label('Spectral Index (1.4GHz/150MHz)')
        #plt.savefig(dir+'projs/%s_proj_spectral_index.png' % ds.basename)
