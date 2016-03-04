#!/usr/bin/env python
import os
import yt
from yt_emissivity import *
yt.enable_parallelism()
import logging
logging.getLogger('yt').setLevel(logging.INFO)

#nu = yt.YTQuantity(500, 'MHz')

dir = '/home/ychen/d9/FLASH4/stampede/0529_L45_M10_b1_h1/'
ts = yt.DatasetSeries(os.path.join(dir,'*_hdf5_plt_cnt_*00'), parallel=7)

figuredir = os.path.join(dir, 'synchrotron')
if yt.is_root():
    if not os.path.exists(os.path.join(dir,figuredir)):
        os.mkdir(os.path.join(dir, figuredir))

for ds in ts.piter():
    #proj_axis = 'x'
    #proj = yt.ProjectionPlot(ds, proj_axis, ('gas', 'density'), center=(0,0,0))
    ##proj.set_zlim(('gas', 'density'), 1E-5, 1E-2)
    ##proj.set_cmap(('flash', 'jet '), 'gist_heat')
    #proj.annotate_timestamp(corner='upper_left', time_format="{time:6.3f} {units}", time_unit='Myr', draw_inset_box=True)
    #proj.annotate_text((0.85, 0.95), dir.split('b1_')[-1].strip('/'), coord_system='axis', text_args={'color':'k'})
    #proj.zoom(10)
    #savefn = 'Projection_%s_density_%s.png' % (proj_axis, str(ds).split('_')[-1])
    #proj.save(os.path.join(figuredir,savefn))

    for nu in [(150, 'MHz'), (1.4, 'GHz')]:
        norm = yt.YTQuantity(*nu).in_units('GHz').value**0.5
        pars = add_emissivity(ds, nu=nu)
        #print pars

        #field = pars[2]
        #proj_axis = 'x'
        #slice = yt.SlicePlot(ds, proj_axis, ('deposit', field), center=(0,0,0),
        #                        field_parameters={'frequency': nu})
        #slice.set_zlim(field, 1E-36/norm, 1E-32/norm)
        #slice.annotate_timestamp(corner='upper_left', time_format="{time:6.3f} {units}", time_unit='Myr', draw_inset_box=True)
        #slice.annotate_text((0.85, 0.95), dir.split('b1_')[-1].strip('/'), coord_system='axis', text_args={'color':'k'})
        #slice.annotate_grids()
        #slice.zoom(10)
        #savefn = 'Slice_%s_emissivity_%s_%s.png' % (proj_axis, str(nu).replace(' ',''), str(ds).split('_')[-1])
        #slice.save(os.path.join(figuredir,savefn))

        field = pars[4]
        proj_axis = 'x'
        proj = yt.ProjectionPlot(ds, proj_axis, field, center=(0,0,0),
                                        field_parameters={'frequency': nu})
        proj.set_zlim(field, 1E-3/norm, 1E1/norm)
        proj.annotate_timestamp(corner='upper_left', time_format="{time:6.3f} {units}", time_unit='Myr', draw_inset_box=True)
        proj.annotate_text((0.85, 0.95), dir.split('b1_')[-1].strip('/'), coord_system='axis', text_args={'color':'k'})
        proj.zoom(10)
        proj.save(figuredir)
        #savefn = 'Projection_%s_%s_%s_%s.png' % (proj_axis, field, str(nu).replace(' ',''), str(ds).split('_')[-1])
        #proj.save(os.path.join(figuredir,savefn))
