#!/usr/bin/env python
import os
import yt
from yt_emissivity import *
yt.enable_parallelism()
import logging
logging.getLogger('yt').setLevel(logging.INFO)


#nu = yt.YTQuantity(500, 'MHz')

dir = '/home/ychen/d9/FLASH4/stampede/1022_L45_M10_b1_h1_10Myr/'
#fname = '/home/ychen/d9/FLASH4/stampede/0529_L45_M10_b1_h1/MHD_Jet_hdf5_plt_cnt_0320'
#ds = yt.load(fname, particle_filename=fname.replace('plt_cnt', 'part'))
ts = yt.DatasetSeries(os.path.join(dir,'*_hdf5_plt_cnt_*50'), parallel=4)

figuredir = os.path.join(dir, 'synchrotron')
if yt.is_root():
    if not os.path.exists(os.path.join(dir,figuredir)):
        os.mkdir(os.path.join(dir, figuredir))

for ds in ts.piter():
    for nu in [yt.YTQuantity(150, 'MHz'), yt.YTQuantity(1.4, 'GHz')]:
        add_emissivity(ds, nu=nu)

        #field = 'avgfill_emissivity'
        #slice = yt.SlicePlot(ds, 'x', ('deposit', field), center=(0,0,0),
        #                        field_parameters={'frequency': nu})
        #slice.set_zlim(('deposit', field), 1E-36/nu.in_units('GHz').value**0.5,
        #                   1E-32/nu.in_units('GHz').value**0.5)
        #slice.annotate_grids()
        #slice.zoom(10)
        #slice.save('Slice_emissivty_%s.png' % nu)

        field = 'avgfill_intensity'
        proj_axis = 'x'
        proj = yt.ProjectionPlot(ds, proj_axis, ('deposit', field), center=(0,0,0),
                                        field_parameters={'frequency': nu})
        proj.set_zlim(('deposit', field), 1E-3/nu.in_units('GHz').value**0.5, 
                              1E1/nu.in_units('GHz').value**0.5)
        proj.annotate_timestamp(corner='upper_left', time_format="{time:6.3f} {units}", time_unit='Myr', draw_inset_box=True)
        proj.annotate_text((0.85, 0.95), dir.split('b1_')[-1].strip('/'), coord_system='axis', text_args={'color':'k'})
        proj.zoom(10)
        #proj.save(dir+'Synchrotron/')
        savefn = 'Projection_%s_intensity_%s_%s.png' % (proj_axis, str(nu).replace(' ',''), str(ds).split('_')[-1])
        proj.save(os.path.join(figuredir,savefn))
