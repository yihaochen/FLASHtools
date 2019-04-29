#!/usr/bin/env python
import os
import numpy as np
import yt
from yt_cluster_ratio_fields import *
yt.enable_parallelism()

dir = './'
ts = yt.DatasetSeries(os.path.join(dir,'data/*_hdf5_plt_cnt_1364'), parallel=1)

figuredir = os.path.join(dir, 'volume_rendering_entropy_ratio_300Myr_perspective')
tfdir = os.path.join(figuredir, 'transfer_function')

if yt.is_root():
    for subdir in [figuredir, tfdir]:
        if not os.path.exists(subdir):
            os.mkdir(subdir)

#fname = '/home/ychen/d9/FLASH4/stampede/0529_L45_M10_b1_h1/MHD_Jet_hdf5_plt_cnt_0620'

bounds = (0.4, 2)

# Since this rendering is done in log space, the transfer function needs
# to be specified in log space.
tf = yt.ColorTransferFunction(np.log10(bounds))

tf.sample_colormap(-0.35, 0.001, alpha=0.9,  colormap="arbre")
tf.sample_colormap(-0.2 , 0.001, alpha=0.3,  colormap="arbre")
tf.sample_colormap(-0.1 , 0.001, alpha=0.2,  colormap="arbre")
tf.sample_colormap( 0.2 , 0.001, alpha=0.4,  colormap="arbre")
tf.sample_colormap( 0.3 , 0.0005, alpha=0.7,  colormap="arbre")

for ds in ts[::-1].piter():
    figfname = figuredir+'/'+ds.basename + '.png'
    #if os.path.exists(figfname): continue
    sp = ds.sphere([0,0,0], (150, 'kpc'))
    sc = yt.create_scene(sp, field='entropy_ratio', lens_type='perspective')

    render_source = sc.get_source(0)
    render_source.transfer_function = tf
    render_source.tfh.tf = tf
    render_source.tfh.bounds = bounds
    #render_source.set_use_ghost_zones(True)

    render_source.tfh.plot(tfdir+'/%s_render_transfer_function.png' % ds.basename,
                           profile_field='density')

    cam = sc.camera
    cam.resolution = (1920, 1080)
    cam.width = ds.arr([512, 288, 512], 'kpc')
    cam.position = ds.arr([256, 0, 0], 'kpc')
    cam.switch_orientation(normal_vector=[-1, 0, 0],
    north_vector = [0, 1, 0])
    #cam.yaw(np.pi/180*70, [0,0,0])
    for i in range(0, 360):
        figfname = figuredir+'/'+ds.basename+'_%03ideg.png'
        sc.render()
        sc.save(figfname % (i*1), sigma_clip=4.0)
        cam.yaw(np.pi/180*1, [0,0,0])
    #sc.save_annotated(figuredir+'/'+ds.basename+'_annotated', sigma_clip=4.0,
    #        text_annotate=[(0.05,0.9), '%5.1f Myr' % ds.current_time.in_units('Myr'),\
    #            {'color': 'grey'}])
