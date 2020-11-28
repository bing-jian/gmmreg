#!/usr/bin/env python
#coding=utf-8
"""
This Python script can be used to test the gmmreg algorithm on the Stanford
"dragon stand" dataset as described in Section 6.1 of the gmmreg PAMI paper.
"""

import collections, copy, os, subprocess, time

import numpy as np
from scipy.spatial.transform import Rotation
import open3d as o3d

try:
    from configparser import ConfigParser
except ImportError:
    from ConfigParser import ConfigParser  # ver. < 3.0

DATA_PATH = '../data/dragon_stand'
CONFIG_FILE = './dragon_stand.ini'
BINARY_DIR = '../C++/build'
GMMREG_BINARY = {'nt': r'gmmreg_demo.exe', 'posix': r'gmmreg_demo'}
BINARY_FULLPATH = os.path.join(BINARY_DIR, GMMREG_BINARY[os.name])


def load_dragon_conf(data_path):
    conf = os.path.join(data_path, 'dragonStandRight.conf')
    lines = open(conf, 'r').read().splitlines()
    poses = {}
    for line in lines[2::]:  # skip the first 2 lines
        """
        bmesh dragonStandRight_0.ply -0.00042573013 -3.6194921e-05 -0.00096614892 0.00068114395 -0.00029109037 -0.00073466112 0.99999946
        """
        words = line.split()
        basename = words[1][:-4]
        idx = int(os.path.basename(basename).split('_')[1]) / 24
        param = [float(x) for x in words[2::]]  # tx, ty, tz, q0, q1, q2, q3
        poses[idx] = param
    return poses


def lookup_ground_truth(gt_poses, i, j):
    pose_i = gt_poses[i]
    pose_j = gt_poses[j]
    q_i = Rotation.from_quat(pose_i[3::])
    q_j = Rotation.from_quat(pose_j[3::])
    q_ji = (q_i * q_j.inv()).as_quat()
    q_ij = (q_i * q_i.inv()).as_quat()
    return q_ij, q_ji


def batch_downsample(data_dir, output_dir, voxel_size=0.005):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    txt_files = []
    for k in range(0, 360, 24):
        ply = os.path.join(data_dir, 'dragonStandRight_%d.ply' % k)
        pcd = o3d.io.read_point_cloud(ply)
        pcd_new = pcd.voxel_down_sample(voxel_size=voxel_size)
        txt = os.path.join(output_dir, '%d.txt' % (k))
        np.savetxt(txt, np.array(pcd_new.points))
        txt_files.append(txt)
    return txt_files


def _get_cfg_and_result_files(working_dir, tuuid):
    config_file = os.path.join(working_dir, '%s.ini' % tuuid)
    output_dir = os.path.join(working_dir, tuuid)
    final_rigid = os.path.join(output_dir, 'final_rigid.txt')
    final_rigid_matrix = os.path.join(output_dir, 'final_rigid_matrix.txt')
    transformed_model = os.path.join(output_dir, 'transformed_model.txt')
    elapsed_time_in_ms = os.path.join(output_dir, 'elapsed_time_in_ms.txt')
    return config_file, final_rigid, final_rigid_matrix, transformed_model, elapsed_time_in_ms


def update_config(config_tmpl, txt_dir, i, j, working_dir):
    config = ConfigParser()
    config.read(config_tmpl)
    txt_i = os.path.join(txt_dir, '%d.txt' % (i * 24))
    txt_j = os.path.join(txt_dir, '%d.txt' % (j * 24))
    pts_i = np.loadtxt(txt_i)
    pts_j = np.loadtxt(txt_j)
    model_scene_switched = False
    if pts_i.shape[0] > pts_j.shape[0]:
        config.set('FILES', 'model', txt_j)
        config.set('FILES', 'scene', txt_i)
        model_scene_switched = True
    else:
        config.set('FILES', 'model', txt_i)
        config.set('FILES', 'scene', txt_j)

    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    tuuid = '{}_vs_{}'.format(i * 24, j * 24)
    output_dir = os.path.join(working_dir, tuuid)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    config_file, final_rigid, final_rigid_matrix, transformed_model, elapsed_time_in_ms = _get_cfg_and_result_files(
        working_dir, tuuid)
    config.set('FILES', 'final_rigid', final_rigid)
    config.set('FILES', 'final_rigid_matrix', final_rigid_matrix)
    config.set('FILES', 'transformed_model', transformed_model)
    config.set('FILES', 'elapsed_time_in_ms', elapsed_time_in_ms)
    with open(config_file, 'w') as fp:
        config.write(fp)
    return config_file, model_scene_switched


class BatchProcessor(object):

    def __init__(self, data_dir, config_tmpl, working_dir):
        self.data_dir_orig = data_dir
        self.gt_poses = load_dragon_conf(self.data_dir_orig)
        self.config_tmpl = config_tmpl
        self.working_dir = working_dir
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)
        self.data_dir = None
        self.reg_results = {}
        self.run_time = {}
        self.eval_results = {}
        print('data_dir_orig: {}'.format(data_dir))
        print('working_dir: {}'.format(working_dir))
        print('config_tmpl: {}'.format(config_tmpl))

    def set_gmmreg_exe(self, gmmreg_exe_path):
        self.gmmreg_exe = gmmreg_exe_path

    def downsample(self, voxel_size):
        self.voxel_size = voxel_size
        output_dir = os.path.join(self.working_dir,
                                  'data_voxel_size_%.3f' % voxel_size)
        batch_downsample(self.data_dir_orig, output_dir)
        self.data_dir = output_dir
        print('data_dir: {}'.format(self.data_dir))

    def process(self, i, j):
        f_config, model_scene_switched = update_config(self.config_tmpl,
                                                       self.data_dir, i, j,
                                                       self.working_dir)
        cmd = '%s %s %s' % (self.gmmreg_exe, f_config, 'rigid')
        print(cmd)
        t1 = time.time()
        subprocess.call(cmd, shell=True)
        t2 = time.time()
        print("Run time : %s seconds" % (t2 - t1))
        tuuid = '{}_vs_{}'.format(i * 24, j * 24)
        _, final_rigid, final_rigid_matrix, transformed_model, elapsed_time_in_ms = _get_cfg_and_result_files(
            self.working_dir, tuuid)
        rigid_param = np.loadtxt(final_rigid)
        matrix = np.loadtxt(final_rigid_matrix)
        elapsed_time_in_ms = np.loadtxt(elapsed_time_in_ms).item()
        if model_scene_switched:
            matrix = np.linalg.inv(matrix)
            rigid_param[0:3] = -rigid_param[0:3]
            rigid_param[4:7] = matrix[:3, 3].T
        self.reg_results[(i, j)] = rigid_param, matrix
        self.run_time[(i, j)] = elapsed_time_in_ms, t2 - t1

    def evaluate(self, i, j):
        assert ((i, j) in self.reg_results)
        q_ij, q_ji = lookup_ground_truth(self.gt_poses, i, j)
        param = self.reg_results[(i, j)][0]
        q = Rotation.from_quat(param[0:4]).as_quat()
        # Distance between rotations
        # http://www.boris-belousov.net/2016/12/01/quat-dist/
        self.eval_results[(i, j)] = np.abs(np.dot(q, q_ji))

    def run_fixed_step(self, step):
        for i in range(15):
            j = (i + step) % 15
            self.process(i, j)
            self.evaluate(i, j)

    def summarize(self):
        print("===========================================")
        scores = collections.defaultdict(list)
        for i, j in self.eval_results:
            step = abs(i - j)
            if step > 7:
                step = 15 - step
            scores[step].append(self.eval_results[(i, j)])
        for step in scores:
            score = np.array(scores[step])
            print("\nSummary stats for registering pairs that are %d degrees away from each other:" % (step * 24))
            print("<avg_score, min_score, max_score, median_score>: %f, %f, %f, %f" %
                  (score.mean(), score.min(), score.max(), np.median(score)))
            print("<# of high accuracy scores (>.99)>: %d out of %d)" %
                  (len(np.where(score > .99)[0]), len(score)))
            error = 2 * np.arccos(score) * 180 / np.pi
            print(
                "<avg_err, min_err, max_err, median_err>: %f, %f, %f, %f (in degrees)"
                % (error.mean(), error.min(), error.max(), np.median(error)))
            print("<# of small errors (<3 degree)>: %d out of %d)" %
                  (len(np.where(error < 3)[0]), len(error)))

    def visualize(self):
        for i, j in self.eval_results:
            if self.eval_results[(i, j)] < 0.99:
                continue
            rigid_param, matrix = self.reg_results[(i, j)]
            model_ply = os.path.join(self.data_dir_orig,
                                     'dragonStandRight_%d.ply' % (i * 24))
            scene_ply = os.path.join(self.data_dir_orig,
                                     'dragonStandRight_%d.ply' % (j * 24))
            pcd_model = o3d.io.read_point_cloud(model_ply)
            pcd_scene = o3d.io.read_point_cloud(scene_ply)
            matrix[:3, 3] = 0  # clear translation
            pcd_model_copy = copy.deepcopy(pcd_model)
            pcd_model_copy.transform(matrix)
            pcd_model.paint_uniform_color([1, 0, 0])  # red
            pcd_scene.paint_uniform_color([0, 1, 0])  # green
            pcd_model_copy.paint_uniform_color([0, 0, 1])  # blue
            o3d.visualization.draw_geometries(
                [pcd_model, pcd_scene, pcd_model_copy])
            break


# https://stackoverflow.com/questions/4248399/page-range-for-printing-algorithm
def parse_range(astr):
    result = set()
    for part in astr.split(','):
        x = part.split('-')
        result.update(range(int(x[0]), int(x[-1]) + 1))
    return sorted(result)


# Run pair-wise registrations, record errors and run time.
def main(args):
    if not args.working_dir:
        tuuid = base64.urlsafe_b64encode(
            uuid.uuid4().bytes).rstrip(b'=').decode('ascii')
        working_dir = '/tmp/{}'.format(tuuid)
    else:
        working_dir = args.working_dir
    batch_processor = BatchProcessor(args.data_dir, args.config_tmpl,
                                     working_dir)
    batch_processor.downsample(args.voxel_size)
    if not args.gmmreg_exe:
        batch_processor.set_gmmreg_exe(BINARY_FULLPATH)
    else:
        batch_processor.set_gmmreg_exe(args.gmmreg_exe)
    steps = set(parse_range(args.intervals))
    steps.intersection_update({1, 2, 3, 4})
    for step in steps:
        batch_processor.run_fixed_step(step)
        batch_processor.run_fixed_step(-step)
    batch_processor.summarize()
    if args.visualize:
        batch_processor.visualize()


import argparse
import uuid, base64
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Testing scrpt for running gmmreg on Dragonstand dataset.')

    parser.add_argument(
        '--data_dir',
        default='../data/dragon_stand',
        required=False,
        help=
        'Directory containing original dragonstand dataset (ply files) and ground truth poses.'
    )

    parser.add_argument('--config_tmpl',
                        default='./dragon_stand.ini',
                        required=False,
                        help='Template config file for running gmmreg.')
    parser.add_argument('--working_dir',
                        required=False,
                        help='Output directory',
                        default='')
    parser.add_argument('--voxel_size',
                        default=0.005,
                        type=float,
                        required=False,
                        help='Run voxel_down_sample(voxel_size=...).')
    parser.add_argument('--gmmreg_exe',
                        required=False,
                        default='',
                        help='Path to the gmmreg executable.')
    parser.add_argument(
        '--intervals',
        required=False,
        default='1,2,3',
        help=
        'Specify the subsets of pairs to be tested. If number n (taking value from {1,2,3,4} is included, '
        'then all pairs whose pose difference is nx24 degrees will be tested.')
    parser.add_argument(
        '--visualize',
        required=False,
        action='store_true',
        help='If true, an example registration result will be visualized.')
    args = parser.parse_args()
    main(args)
