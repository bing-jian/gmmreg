#!/usr/bin/env python
#coding=utf-8

import os, subprocess, time
import numpy as np

TMP_PATH = './tmp'

if not os.path.exists(TMP_PATH):
    os.makedirs(TMP_PATH)


BINARY_DIR = '../C++/build'
GMMREG_BINARY = {
        'nt' : r'gmmreg_demo.exe',
        'posix': r'gmmreg_demo'
}
BINARY_FULLPATH = os.path.join(BINARY_DIR, GMMREG_BINARY[os.name])


# TODO(bing-jian): use a better downsample method.
def downsampled_subset(x, num):
    downsample_rate = max(int(np.floor(len(x) / num * 1.0)), 1)
    return x[::downsample_rate]

# Run pairwise rigid registration using specified binary and configuration.
def run_rigid_pairwise(gmmreg_exe, model, scene, f_config):
    if type(model) == type('abc'):
        model = np.loadtxt(model)
        scene = np.loadtxt(scene)
    if model.shape[0] > 5000:
        model = downsampled_subset(model, 5000)
    if scene.shape[0] > 5000:
        scene = downsampled_subset(scene, 5000)
    print(model.shape)
    print(scene.shape)
    np.savetxt(os.path.join(TMP_PATH, 'model.txt'), model)
    np.savetxt(os.path.join(TMP_PATH, 'scene.txt'), scene)
    cmd = '%s %s %s'%(gmmreg_exe, f_config, 'rigid')

    t1 = time.time()
    subprocess.call(cmd, shell=True)
    t2 = time.time()
    print("Run time : %s seconds" % (t2 - t1))
    param = np.loadtxt(os.path.join(TMP_PATH, 'final_rigid.txt'))
    matrix = np.loadtxt(os.path.join(TMP_PATH, 'final_rigid_matrix.txt'))
    elasped_time_in_ms = np.loadtxt(os.path.join(TMP_PATH, 'elapsed_time_in_ms.txt'))

    return param, matrix, elasped_time_in_ms.item(), t2-t1


# http://qianyi.info/scenedata.html
def convert_depth_to_pcloud(depth):
    fx = 525.0; fy = 525.0; # default focal length
    cx = 319.5; cy = 239.5; # default optical center

    ydim, xdim = depth.shape
    depth = depth.flatten()
    valid_ind = np.where(depth > 0)
    xx = range(xdim)
    yy = range(ydim)
    u, v = np.meshgrid(xx, yy)
    # translation from depth pixel (u,v,d) to a point (x,y,z)
    u = u.flatten()
    v = v.flatten()
    u = u[valid_ind]
    v = v[valid_ind]
    depth = depth[valid_ind]
    z = depth / 1000.0;  # convert the unit from mm to meter.
    x = (u - cx) * z / fx;
    y = (v - cy) * z / fy;
    return np.c_[x, y, z]
