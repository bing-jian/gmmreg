#!/usr/bin/env python
#coding=utf-8
"""
This Python script can be used to test the gmmreg algorithm on the Stanford
"lounge" dataset.  For how to get this dataset, please refer to
http://qianyi.info/scenedata.html

Dataset can be downloaded with:
  wget https://github.com/isl-org/open3d_downloads/releases/download/20220301-data/LoungeRGBDImages.zip
"""

import argparse
import copy
import os
import subprocess
import time

import cv2
import numpy as np
import open3d as o3d

from common_utils import *

_HERE = os.path.dirname(os.path.abspath(__file__))
BINARY_DIR = os.path.join(_HERE, '..', 'C++', 'build')
GMMREG_BINARY = {'nt': r'gmmreg_demo.exe', 'posix': r'gmmreg_demo'}

BINARY_FULLPATH = os.path.join(BINARY_DIR, GMMREG_BINARY[os.name])
TMP_PATH = os.path.join(_HERE, 'tmp')
CONFIG_FILE = os.path.join(_HERE, 'lounge.ini')

if not os.path.exists(TMP_PATH):
    os.makedirs(TMP_PATH)

# Populated by _init_data() before any registration calls.
GT_POS = None
DEPTH_FILES = None


def _init_data(data_path):
    global GT_POS, DEPTH_FILES
    trajectory_path = os.path.join(data_path, 'lounge_trajectory.log')
    GT_POS = parse_ground_truth(trajectory_path)
    DEPTH_FILES = get_all_depth_files(data_path)


def parse_ground_truth(trajectory_conf):
    lines = open(trajectory_conf).readlines()
    matrices = {}
    for i in range(3000):
        ll = lines[5 * i + 1:5 * i + 5]
        ss = ''.join(ll)
        mat = np.fromstring(ss, dtype=float, sep='\t').reshape(4, 4)
        matrices[i] = mat
    return matrices


def get_all_depth_files(data_path):
    ind = range(1, 3001)
    return [os.path.join(data_path, 'depth/%06d.png' % j) for j in ind]


def draw_registration_result(source, target, transformation):
    source_temp = copy.deepcopy(source)
    target_temp = copy.deepcopy(target)
    source_temp.paint_uniform_color([1, 0.706, 0])
    target_temp.paint_uniform_color([0, 0.651, 0.929])
    source_temp.transform(transformation)
    o3d.visualization.draw_geometries([source_temp, target_temp])


# Please comment the code using open3d (www.open3d.org) if it is not installed.
def run_pairwise_registration(i, j, visualize=False, icp_refine=False):
    model_depth = cv2.imread(DEPTH_FILES[i],
                             cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)
    scene_depth = cv2.imread(DEPTH_FILES[j],
                             cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)
    model = convert_depth_to_pcloud(model_depth)
    scene = convert_depth_to_pcloud(scene_depth)
    pcloud_model = o3d.geometry.PointCloud()
    pcloud_model.points = o3d.utility.Vector3dVector(model)
    pcloud_scene = o3d.geometry.PointCloud()
    pcloud_scene.points = o3d.utility.Vector3dVector(scene)

    if visualize:
        pcloud_model.paint_uniform_color([1, 0, 0])  # red
        pcloud_scene.paint_uniform_color([0, 1, 0])  # green
        o3d.visualization.draw_geometries([pcloud_model, pcloud_scene])

    down_model = pcloud_model.voxel_down_sample(voxel_size=0.065)
    down_scene = pcloud_scene.voxel_down_sample(voxel_size=0.065)
    model_pts = np.array(down_model.points)
    scene_pts = np.array(down_scene.points)

    res = run_rigid_pairwise(BINARY_FULLPATH, model_pts, scene_pts, CONFIG_FILE)
    # http://www.boris-belousov.net/2016/12/01/quat-dist/
    print("Transformation estimated by gmmreg:")
    print(res[1])
    gt = np.dot(np.linalg.inv(GT_POS[j]), GT_POS[i])
    print("Transformation from ground truth:")
    print(gt)
    theta_before = np.arccos((np.trace(gt[:3, :3]) - 1) * 0.5)
    print("pose difference (in degrees) before alignment:",
          theta_before * 180 / np.pi)
    R = np.dot(gt[:3, :3].T, res[1][:3, :3])
    theta_after = np.arccos((np.trace(R) - 1) * 0.5)
    print("pose difference (in degrees) after alignment:",
          theta_after * 180 / np.pi)
    transformed = np.loadtxt(os.path.join(TMP_PATH, 'transformed_model.txt'))
    pcloud_transformed = o3d.geometry.PointCloud()
    pcloud_transformed.points = o3d.utility.Vector3dVector(transformed)
    pcloud_transformed.paint_uniform_color([0, 0, 1])  # blue
    matrix = np.loadtxt(os.path.join(TMP_PATH, 'final_rigid_matrix.txt'))
    if visualize:
        o3d.visualization.draw_geometries(
            [pcloud_transformed, down_model, down_scene])
        transformed = np.dot(model, matrix[:3, :3].T) + matrix[:3, 3].T
        pcloud_transformed.points = o3d.utility.Vector3dVector(transformed)
        pcloud_transformed.paint_uniform_color([0, 0, 1])  # blue
        o3d.visualization.draw_geometries([pcloud_transformed, pcloud_scene])

    if icp_refine:
        print("Apply point-to-point ICP")
        threshold = 0.02
        trans_init = matrix
        t1 = time.time()
        reg_p2p = o3d.pipelines.registration.registration_icp(
            down_model, down_scene, threshold, trans_init,
            o3d.pipelines.registration.TransformationEstimationPointToPoint())
        t2 = time.time()
        print("ICP Run time : %s seconds" % (t2 - t1))
        print(reg_p2p)
        print("Transformation by ICP is:")
        print(reg_p2p.transformation)
        print("")
        R = np.dot(gt[:3, :3].T, reg_p2p.transformation[:3, :3])
        theta = np.arccos((np.trace(R) - 1) * 0.5)
        print("pose difference (in degrees) after icp-refinement:",
              theta * 180 / np.pi)
        if visualize:
            draw_registration_result(pcloud_model, pcloud_scene,
                                     reg_p2p.transformation)
    return res, theta_before * 180 / np.pi, theta_after * 180 / np.pi


# Run pair-wise registrations, record errors and run time.
def main():
    o = []
    for i in range(0, 2995):
        j = i + 5
        res, theta_before, theta_after = run_pairwise_registration(
            i, j, visualize=False, icp_refine=False)
        o.append((theta_before, theta_after, res[-2]))
    o = np.array(o)
    print("-----------")
    print("Input pairs with pose difference ~= %f degrees" % (np.mean(o[:, 0])))
    error = o[:, 1]
    print(
        "<avg_err, min_err, max_err, median_err>: %f, %f, %f, %f (in degrees)" %
        (np.nanmean(error), np.nanmin(error), np.nanmax(error),
         np.nanmedian(error)))
    print("<# of small errors (<3 degree)>: %d out of %d)" %
          (len(np.where(error < 3)[0]), len(error)))
    core_run_time = o[:, 2]
    print(
        "Registration time: <avg_time, min_time, max_time, median_time>: %f, %f, %f, %f (in milliseconds)"
        % (core_run_time.mean(), core_run_time.min(), core_run_time.max(),
           np.median(core_run_time)))
    np.savetxt('./tmp/loung_expts_results.txt', o)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Testing script for running gmmreg on the Lounge RGB-D dataset.')
    parser.add_argument(
        '--data_path',
        default='../data/lounge',
        help='Directory containing the Lounge dataset (depth/ subfolder and lounge_trajectory.log).')
    args = parser.parse_args()

    _init_data(args.data_path)

    # Register just one pair, visualize point clouds before/after alignment.
    #run_pairwise_registration(1, 6, visualize=True)
    #run_pairwise_registration(20, 26, visualize=True)
    run_pairwise_registration(1, 11, visualize=True)
    #run_pairwise_registration(900, 910)
    main()
