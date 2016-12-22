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
from yt.utilities.fits_image import FITSProjection, FITSSlice
logging.getLogger('yt').setLevel(logging.INFO)

#nu = yt.YTQuantity(500, 'MHz')


#dir = '/home/ychen/data/0only_0529_h1/'
#dir = '/home/ychen/data/0only_0605_hinf/'
dir = '/home/ychen/data/0only_1022_h1_10Myr/'
#dir = '/d/d8/ychen/MHD_Jet/0314_L45_M10_b1_h1_nojiggle'
ts = yt.DatasetSeries(os.path.join(dir,'*_hdf5_plt_cnt_1060'), parallel=1)
#ts = yt.DatasetSeries(os.path.join(dir,'*_hdf5_plt_cnt_1410'), parallel=1)
#ts = yt.DatasetSeries(os.path.join(dir,'*_hdf5_plt_cnt_0600'), parallel=8)

zoom_fac = 4

ptype = 'lobe'
maindir = os.path.join(dir, 'pol_synchrotron_QU_nn_%s/' % ptype)
fitsdir = os.path.join(maindir, 'fits/')
if yt.is_root():
    for subdir in [maindir, fitsdir]:
        if not os.path.exists(subdir):
            os.mkdir(subdir)

for ds in ts.piter():
    proj_axis = 'x'
    #fields = []
    width = ds.domain_width[1:]/zoom_fac
    res = ds.domain_dimensions[1:]*ds.refine_by**ds.index.max_level/zoom_fac
    fits_slice = FITSSlice(ds, proj_axis, 'magnetic_field_strength', center=[0,0,0], width=width, image_res=res)
    slicefitsfname = os.path.join(fitsdir, 'magnetic_field_strength_%s.fits' % ds.basename[-4:])
    fits_slice.writeto(slicefitsfname)
    #for nu in [(150, 'MHz'), (233, 'MHz'), (325, 'MHz'), (610, 'MHz'), (1400, 'MHz')]:
    ##for nu in [(233, 'MHz'), (325, 'MHz'), (610, 'MHz'), (1.4, 'GHz')]:
    #    ###########################################################################
    #    ## Polarizations
    #    ###########################################################################

    #    pars = add_synchrotron_pol_emissivity(ds, ptype=ptype, nu=nu, proj_axis=proj_axis)
    #    #for pol in ['i', 'q', 'u']:
    #    for pol in ['i']:
    #        field = ('nn_emissivity_%s_%s_%%.1f%%s' % (pol, ptype)) % nu
    #        #fields.append(field)


    #    fits_proj = FITSProjection(ds, proj_axis, field, center=[0,0,0], width=width, image_res=res)

    #    fitsfname = os.path.join(fitsdir, ('synchrotron_i_%s_%s_%%i%%s.fits' % (ptype, ds.basename[-4:])) % nu)
    #    fits_proj.writeto(fitsfname)