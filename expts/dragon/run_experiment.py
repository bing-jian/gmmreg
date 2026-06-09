#!/usr/bin/env python3
"""Step 2: run pairwise rigid registrations on the prepared dragon stand data.

Reads the downsampled point clouds and GT poses produced by prepare_data.py,
runs gmmreg_demo for each pair, evaluates rotation accuracy, and prints a
summary.

Usage
-----
    python run_experiment.py \\
        --prepared_dir /path/to/prepared \\
        --output_dir   /path/to/results  \\
        [--config_tmpl dragon_stand.ini] \\
        [--gmmreg_exe  /path/to/gmmreg_demo] \\
        [--intervals   1,2,3]
"""

import argparse
import collections
import os
import subprocess
import time

import numpy as np
from configparser import ConfigParser
from scipy.spatial.transform import Rotation

_HERE = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_CONFIG = os.path.join(_HERE, 'dragon_stand.ini')
_DEFAULT_EXE = os.path.join(_HERE, '..', '..', 'C++', 'build', 'gmmreg_demo')

NUM_VIEWS = 15
DEGREE_STEP = 24


def load_gt_poses(gt_dir):
    """Load GT poses saved by prepare_data.py → {0-based idx: [tx,ty,tz,q0,q1,q2,q3]}."""
    poses = {}
    for i in range(NUM_VIEWS):
        poses[i] = np.loadtxt(os.path.join(gt_dir, f'gt_pose_{i:02d}.txt')).tolist()
    return poses


def lookup_ground_truth(gt_poses, i, j):
    q_i = Rotation.from_quat(gt_poses[i][3:])
    q_j = Rotation.from_quat(gt_poses[j][3:])
    return (q_j * q_i.inv()).as_quat(), (q_i * q_j.inv()).as_quat()


def _result_files(working_dir, tag):
    out = os.path.join(working_dir, tag)
    return (os.path.join(working_dir, f'{tag}.ini'),
            os.path.join(out, 'final_rigid.txt'),
            os.path.join(out, 'final_rigid_matrix.txt'),
            os.path.join(out, 'transformed_model.txt'),
            os.path.join(out, 'elapsed_time_in_ms.txt'))


def update_config(config_tmpl, pcloud_dir, i, j, working_dir):
    config = ConfigParser()
    config.read(config_tmpl)
    txt_i = os.path.join(pcloud_dir, f'{i * DEGREE_STEP}.txt')
    txt_j = os.path.join(pcloud_dir, f'{j * DEGREE_STEP}.txt')
    pts_i = np.loadtxt(txt_i)
    pts_j = np.loadtxt(txt_j)
    switched = pts_i.shape[0] > pts_j.shape[0]
    config.set('FILES', 'model', txt_j if switched else txt_i)
    config.set('FILES', 'scene', txt_i if switched else txt_j)

    tag = f'{i * DEGREE_STEP}_vs_{j * DEGREE_STEP}'
    os.makedirs(os.path.join(working_dir, tag), exist_ok=True)
    cfg_file, final_rigid, final_rigid_matrix, transformed_model, elapsed_ms = \
        _result_files(working_dir, tag)
    config.set('FILES', 'final_rigid', final_rigid)
    config.set('FILES', 'final_rigid_matrix', final_rigid_matrix)
    config.set('FILES', 'transformed_model', transformed_model)
    config.set('FILES', 'elapsed_time_in_ms', elapsed_ms)
    with open(cfg_file, 'w') as fp:
        config.write(fp)
    return cfg_file, switched


class BatchProcessor:

    def __init__(self, prepared_dir, config_tmpl, gmmreg_exe, output_dir):
        self.pcloud_dir = os.path.join(prepared_dir, 'pclouds')
        self.gt_poses = load_gt_poses(os.path.join(prepared_dir, 'gt_poses'))
        self.config_tmpl = config_tmpl
        self.gmmreg_exe = gmmreg_exe
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        self.reg_results = {}
        self.run_time = {}
        self.eval_results = {}

    def process(self, i, j):
        cfg, switched = update_config(
            self.config_tmpl, self.pcloud_dir, i, j, self.output_dir)
        cmd = f'{self.gmmreg_exe} {cfg} rigid'
        print(cmd)
        t1 = time.time()
        subprocess.call(cmd, shell=True)
        elapsed_wall = time.time() - t1
        print(f'  wall time: {elapsed_wall:.2f}s')

        tag = f'{i * DEGREE_STEP}_vs_{j * DEGREE_STEP}'
        _, final_rigid, final_rigid_matrix, _, elapsed_ms_file = \
            _result_files(self.output_dir, tag)
        rigid_param = np.loadtxt(final_rigid)
        matrix = np.loadtxt(final_rigid_matrix)
        elapsed_ms = np.loadtxt(elapsed_ms_file).item()
        if switched:
            matrix = np.linalg.inv(matrix)
            rigid_param[0:3] = -rigid_param[0:3]
            rigid_param[4:7] = matrix[:3, 3].T
        self.reg_results[(i, j)] = rigid_param, matrix
        self.run_time[(i, j)] = elapsed_ms, elapsed_wall

    def evaluate(self, i, j):
        assert (i, j) in self.reg_results
        _, q_ji = lookup_ground_truth(self.gt_poses, i, j)
        q = Rotation.from_quat(self.reg_results[(i, j)][0][:4]).as_quat()
        self.eval_results[(i, j)] = np.abs(np.dot(q, q_ji))

    def run_fixed_step(self, step):
        for i in range(NUM_VIEWS):
            self.process(i, (i + step) % NUM_VIEWS)
            self.evaluate(i, (i + step) % NUM_VIEWS)

    def summarize(self):
        print('===========================================')
        scores = collections.defaultdict(list)
        for i, j in self.eval_results:
            step = abs(i - j)
            if step > NUM_VIEWS // 2:
                step = NUM_VIEWS - step
            scores[step].append(self.eval_results[(i, j)])
        for step in sorted(scores):
            score = np.array(scores[step])
            error = 2 * np.arccos(score) * 180 / np.pi
            print(f'\nPairs {step * DEGREE_STEP}° apart:')
            print(f'  score  <avg,min,max,median>: '
                  f'{score.mean():.4f}, {score.min():.4f}, {score.max():.4f}, {np.median(score):.4f}')
            print(f'  high accuracy (>.99): {(score > .99).sum()} / {len(score)}')
            print(f'  error  <avg,min,max,median> (deg): '
                  f'{error.mean():.2f}, {error.min():.2f}, {error.max():.2f}, {np.median(error):.2f}')
            print(f'  small errors (<3°): {(error < 3).sum()} / {len(error)}')


def parse_range(astr):
    result = set()
    for part in astr.split(','):
        x = part.split('-')
        result.update(range(int(x[0]), int(x[-1]) + 1))
    return sorted(result)


def main():
    parser = argparse.ArgumentParser(
        description='Run pairwise rigid registrations on prepared dragon stand data.')
    parser.add_argument('--prepared_dir', required=True,
                        help='Directory from prepare_data.py (contains pclouds/ and gt_poses/)')
    parser.add_argument('--output_dir', required=True,
                        help='Directory for per-pair configs and registration results')
    parser.add_argument('--config_tmpl', default=_DEFAULT_CONFIG,
                        help='Template INI file with [GMMREG_OPT] section')
    parser.add_argument('--gmmreg_exe', default=_DEFAULT_EXE,
                        help='Path to gmmreg_demo executable')
    parser.add_argument('--intervals', default='1,2,3',
                        help='Step sizes to test (subset of {1,2,3,4}); default: 1,2,3')
    args = parser.parse_args()

    processor = BatchProcessor(
        args.prepared_dir, args.config_tmpl, args.gmmreg_exe, args.output_dir)
    for step in set(parse_range(args.intervals)) & {1, 2, 3, 4}:
        processor.run_fixed_step(step)
        processor.run_fixed_step(-step)
    processor.summarize()


if __name__ == '__main__':
    main()
