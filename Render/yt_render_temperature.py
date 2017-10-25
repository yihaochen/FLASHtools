#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import os
import sys
import numpy as np
import yt
yt.enable_parallelism()

#dir = '/home/ychen/data/0only_0529_h1/'
dir = '/d/d11/ychen/MHD_jet/0517_L45_M10_b1_h1_20Myr'
try:
    ind = int(sys.argv[1])
    ts = yt.DatasetSeries(os.path.join(dir,'*_hdf5_plt_cnt_%02d?0' % ind), parallel=10)
except IndexError:
    ts = yt.DatasetSeries(os.path.join(dir,'*_hdf5_plt_cnt_???0'), parallel=10)


figuredir = os.path.join(dir, 'volume_rendering_temperature_clip4')
annotateddir = os.path.join(figuredir, 'annotated')
tfdir = os.path.join(figuredir, 'transfer_function')

if yt.is_root():
    for subdir in [figuredir, tfdir, annotateddir]:
        if not os.path.exists(subdir):
            os.mkdir(subdir)

#fname = '/home/ychen/d9/FLASH4/stampede/0529_L45_M10_b1_h1/MHD_Jet_hdf5_plt_cnt_0620'
#ds = yt.load(fname)

bounds = (2e7, 3e9)

# Since this rendering is done in log space, the transfer function needs
# to be specified in log space.
tf = yt.ColorTransferFunction(np.log10(bounds))

tf.sample_colormap(np.log10(2E9), 0.005, alpha=0.1, colormap="arbre")
#tf.sample_colormap(np.log10(1E9), 0.005, alpha=0.5, colormap="arbre")
#tf.sample_colormap(np.log10(6E8), 0.005, alpha=0.1,  colormap="arbre")
tf.sample_colormap(np.log10(7.5E8), 0.005, alpha=0.1,  colormap="arbre")
tf.sample_colormap(np.log10(2.5E8), 0.005, alpha=0.1,  colormap="arbre")
tf.sample_colormap(np.log10(1E8), 0.005, alpha=0.1,  colormap="arbre")
tf.sample_colormap(np.log10(3.5E7), 0.005, alpha=0.01,  colormap="arbre")

for ds in ts.piter():
    sc = yt.create_scene(ds, field='temperature', lens_type='plane-parallel')

    render_source = sc.get_source(0)
    render_source.transfer_function = tf
    render_source.tfh.tf = tf
    render_source.tfh.bounds = bounds

    render_source.tfh.plot(tfdir+'/%s_render_transfer_function.png' % ds.basename,
                           profile_field='density')

    cam = sc.add_camera(ds, lens_type='plane-parallel')
    cam.resolution = (1920, 1080)
    cam.width = ds.arr([100,56.25,56.25], 'kpc')
    cam.position = ds.arr([50, 0, 0], 'kpc')
    cam.switch_orientation(normal_vector=[-1, 0, 0], north_vector=[0, 1, 0])

    sc.render()


    # save an annotated version of the volume rendering including a representation
    # of the transfer function and a nice label showing the simulation time.
    text_string = "T = {:6.2f} Myr".format(float(ds.current_time.to('Myr')))
    #text_kwargs = {'color': 'grey'}
    sc.save_annotated(annotateddir+'/'+ds.basename, sigma_clip=4,
                              text_annotate=[[(0.05, 0.9), text_string]])


    #sc.save(figuredir+'/'+ds.basename, sigma_clip=4.0)
    sc.save(figuredir+'/'+ds.basename, sigma_clip=4)

