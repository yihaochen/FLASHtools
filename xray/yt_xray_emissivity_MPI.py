#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import yt
#yt.enable_parallelism()
from yt.analysis_modules.spectral_integrator.api import \
     add_xray_emissivity_field
#from yt.visualization.plot_window import get_oblique_window_parameters, OffAxisProjectionDummyDataSource, PWViewerMPL
#from yt.visualization.fixed_resolution import OffAxisProjectionFixedResolutionBuffer
import logging
logging.getLogger('yt').setLevel(logging.WARNING)
import numpy as np
import matplotlib.pyplot as plt
import MPI_taskpull2


#center = [0.0, 0.0, 0.0]
#normal = [0.0, 1.0, 0.5]
#kpc = yt.YTQuantity(1, 'kpc')
#width = [40*kpc, 40*kpc, 80*kpc]
#Npix = 512
#field = 'density'
#
#proj = yt.ProjectionPlot(ds, 'x', field)
#image = yt.off_axis_projection(ds, center, normal, width, Npix, field, num_threads=4)
#yt.write_image(np.log10(image), "%s_offaxis_projection.png" % ds)


def xray_3color(proj_angle=30.):
    dir ='/home/ychen/d9/FLASH4/stampede/0529_L45_M10_b1_h1/'
    fname = dir+'MHD_Jet_hdf5_plt_cnt_0620'
    #fname = '~/d9/yt_testdata/GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0150'

    energy_ranges = [(0.3, 1.5), (1.5, 3.5), (3.5, 7.0)]
    fields = ['xray_emissivity_%.1f_%.1f_keV' % er for er in energy_ranges]

    ds = yt.load(fname)

    for energy_range in energy_ranges:
        print add_xray_emissivity_field(ds, *energy_range, \
                with_metals=False, constant_metallicity=0.02,\
                filename='/d/d9/ychen/yt_testdata/apec_emissivity.h5')

    normal = (0.0, np.tan(proj_angle/180.*np.pi), 1.0)
    north = (0.0, 0.1, 1.0)
    max_lv = 6

    plot = yt.OffAxisProjectionPlot(ds, normal, fields, width=((40,'kpc'), (80,'kpc')),\
                                    depth=(80, 'kpc'), max_level=max_lv, north_vector=north)
    #plot.zoom(10)
    #plot.frb.export_fits(fname+'.fits', fields=fname)
    #plot.save()

    ext = ds.arr([-20, 20, -40, 40], input_units='kpc')

    r = plot['xray_emissivity_0.3_1.5_keV'].image._A
    g = plot['xray_emissivity_1.5_3.5_keV'].image._A
    b = plot['xray_emissivity_3.5_7.0_keV'].image._A
    image = np.array([r,g,b]).transpose([1,2,0])
    image = image/image.max()
    plt.figure(figsize = (10,20))
    plt.imshow(image, origin='lower', extent=ext)
    plt.xlabel('(kpc)')
    #plt.ylabel('(kpc)')
    plt.savefig('xray_apec002_3color%i_lv%i.png' % (int(proj_angle), int(max_lv)))

tasks = iter(range(0, 90, 10))

results = MPI_taskpull2.taskpull(xray_3color, tasks, print_result=True)
