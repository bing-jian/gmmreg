#!/usr/bin/env python
#coding=utf-8

"""
This Python script can be used to test the gmmreg algorithm on the Stanford
"dragon stand" dataset as described in Section 6.1 of the gmmreg PAMI paper.
"""

import os, subprocess, time
import numpy as np

from common_utils import *

# https://github.com/bistromath/gr-air-modes/blob/master/python/Quaternion.py
from Quaternion import Quat, normalize

DATA_PATH = '../data/dragon_stand'
CONFIG_FILE = './dragon_stand.ini'


def load_dragon_conf(data_path):
    conf = os.path.join(data_path, 'dragonStandRight.conf')
    o = open(conf, 'r').read().splitlines()
    pos = {}
    for line in o[2::]:
        words = line.split()
        fname = words[1]
        param = [float(x) for x in words[2::]]
        pos[fname] = param
    return pos

GT_POS = load_dragon_conf(DATA_PATH)

def get_all_plyfiles(data_path):
    ind = range(0, 360, 24)
    return [os.path.join(data_path,
                         'dragonStandRight_%d.ply'%j) for j in ind]


PLY_FILES = get_all_plyfiles(DATA_PATH)



def ply2txt(plyfile):
    return '%s.txt' % os.path.splitext(plyfile)[0]


def run_rigid_batch(step, res):
    for i in range(15):
        j = (i + step) % 15
        model_ply = PLY_FILES[i]
        scene_ply = PLY_FILES[j]
        model_txt = ply2txt(model_ply)
        scene_txt = ply2txt(scene_ply)
        print(model_txt, scene_txt)
        res[(i, j)] = run_rigid_pairwise(BINARY_FULLPATH, model_txt, scene_txt, CONFIG_FILE)
    return res


def lookup_ground_truth(i, j):
    model_ply = PLY_FILES[i]
    scene_ply = PLY_FILES[j]
    pos_i = GT_POS[os.path.basename(model_ply)]
    pos_j = GT_POS[os.path.basename(scene_ply)]
    q_i = Quat(normalize(np.array(pos_i[3::])))
    q_j = Quat(normalize(np.array(pos_j[3::])))
    q_ji = q_i / q_j
    q_ij = q_j / q_i
    return q_ij, q_ji


# Distance between rotations
# http://www.boris-belousov.net/2016/12/01/quat-dist/
def compute_accuracy(reg_result):
    error = {} # error in degrees.
    total_run_time = {}
    core_run_time = {}
    for i, j in reg_result:
        q_ij, q_ji = lookup_ground_truth(i, j)
        param = reg_result[(i,j)][0]
        q = Quat(normalize(param[0:4]))
        similarity = np.abs(np.dot(q.q, q_ji.q))
        error[(i, j)] = 2 * np.arccos(similarity)*180 / np.pi
        total_run_time[(i, j)] = reg_result[(i,j)][-1]
        core_run_time[(i, j)] = reg_result[(i,j)][-2]
    return error, total_run_time, core_run_time


# Run pair-wise registrations, record errors and run time.
def main():
    errors = {}
    total_run_times = {}
    core_run_times = {}
    for step in [1]:
        reg_result = {}
        reg_result = run_rigid_batch(step, reg_result)
        reg_result = run_rigid_batch(-step, reg_result)
        error, total_run_time, core_run_time = compute_accuracy(reg_result)
        errors[step] = np.array([error[ij] for ij in error])
        total_run_times[step] = np.array([total_run_time[ij] for ij in total_run_time])
        core_run_times[step] = np.array([core_run_time[ij] for ij in core_run_time])

    for step in errors:
        print("-----------")
        print("Input pairs with pose difference ~= %d degrees" % (step * 24))
        error = errors[step]
        print("<avg_err, min_err, max_err, median_err>: %f, %f, %f, %f (in degrees)" % (
            error.mean(), error.min(), error.max(), np.median(error)))
        print("<# of small errors (<3 degree)>: %d out of %d)" % (
            len(np.where(error < 3)[0]), len(error)))
        run_time = total_run_times[step]
        print("Total run time: <avg_time, min_time, max_time, median_time>: %f, %f, %f, %f (in seconds)" % (
            run_time.mean(), run_time.min(), run_time.max(), np.median(run_time)))
        print("Please note that the run time reported above includes time spend on reading/writing files.")
        core_run_time = core_run_times[step]
        print("Registration time: <avg_time, min_time, max_time, median_time>: %f, %f, %f, %f (in milliseconds)" % (
            core_run_time.mean(), core_run_time.min(), core_run_time.max(), np.median(core_run_time)))


# Please comment the code using open3d (www.open3d.org) if it is not installed.
def visualize_registration(i, j):
    try:
        from open3d.geometry import PointCloud
        from open3d.utility import Vector3dVector
        from open3d.visualization import draw_geometries
    except:
        from open3d import PointCloud, Vector3dVector, draw_geometries
    model_ply = PLY_FILES[i]
    scene_ply = PLY_FILES[j]
    model_txt = ply2txt(model_ply)
    scene_txt = ply2txt(scene_ply)
    model = np.loadtxt(model_txt)
    scene = np.loadtxt(scene_txt)
    print(model_txt, scene_txt)
    pcloud_model = PointCloud()
    pcloud_model.points = Vector3dVector(model)
    pcloud_model.paint_uniform_color([1, 0, 0]) # red
    pcloud_scene = PointCloud()
    pcloud_scene.points = Vector3dVector(scene)
    pcloud_scene.paint_uniform_color([0, 1, 0]) # green
    draw_geometries([pcloud_model, pcloud_scene])
    res = run_rigid_pairwise(BINARY_FULLPATH, model, scene, CONFIG_FILE)
    print(res)
    transformed = np.loadtxt(os.path.join(TMP_PATH, 'transformed_model.txt'))
    pcloud_transformed = PointCloud()
    pcloud_transformed.points = Vector3dVector(transformed)
    pcloud_transformed.paint_uniform_color([0, 0, 1]) # blue
    draw_geometries([pcloud_transformed, pcloud_scene])


import sys
if __name__ == "__main__":
    # Register just one pair, visualize point clouds before/after alignment.
    try:
        import open3d
        visualize_registration(0, 1)
    except:
        pass
    # Run registration on many more specified pairs and collect metrics.
    main()
