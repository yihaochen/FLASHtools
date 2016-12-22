#!/usr/bin/env python
import os
import numpy as np
import yt
yt.enable_parallelism()

dir = '/home/ychen/data/0605_L45_M10_b1_h0/'
ts = yt.DatasetSeries(os.path.join(dir,'*_hdf5_plt_cnt_??[0,5]0'), parallel=22)

figuredir = os.path.join(dir, 'volume_rendering_plane_parallel')
tfdir = os.path.join(figuredir, 'transfer_function')

if yt.is_root():
    for subdir in [figuredir, tfdir]:
        if not os.path.exists(subdir):
            os.mkdir(subdir)

#fname = '/home/ychen/d9/FLASH4/stampede/0529_L45_M10_b1_h1/MHD_Jet_hdf5_plt_cnt_0620'
#ds = yt.load(fname)

bounds = (1e-28, 1e-25)

# Since this rendering is done in log space, the transfer function needs
# to be specified in log space.
tf = yt.ColorTransferFunction(np.log10(bounds))

tf.sample_colormap(np.log10(1E-25), 0.005, alpha=0.03, colormap="arbre")
#tf.sample_colormap(np.log10(2E-26), 0.005, alpha=0.05, colormap="arbre")
tf.sample_colormap(np.log10(1E-26), 0.005, alpha=0.2,  colormap="arbre")
tf.sample_colormap(np.log10(3E-27), 0.005, alpha=0.5,  colormap="arbre")
tf.sample_colormap(np.log10(1E-27), 0.005, alpha=0.9,  colormap="arbre")

for ds in ts.piter():
    sc = yt.create_scene(ds, field='density', lens_type='plane-parallel')

    render_source = sc.get_source(0)
    render_source.transfer_function = tf
    render_source.tfh.tf = tf
    render_source.tfh.bounds = bounds

    render_source.tfh.plot(tfdir+'/%s_render_transfer_function.png' % ds.basename,
                           profile_field='density')

    cam = sc.add_camera(ds, lens_type='plane-parallel')
    cam.resolution = (1600, 1600)
    cam.width = ds.quan(80, 'kpc')
    cam.position = ds.arr([50, 0, 0], 'kpc')
    cam.switch_orientation(normal_vector=[-1, 0, 0],
    north_vector=[0, 0, 1])

    sc.render()

    sc.save(figuredir+'/'+ds.basename, sigma_clip=4.0)
