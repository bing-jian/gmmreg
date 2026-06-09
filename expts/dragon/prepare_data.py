#!/usr/bin/env python3
"""Step 1: extract per-frame downsampled point clouds and GT poses.

Reads the Stanford dragon stand dataset, voxel-downsamples each PLY frame,
and saves GT poses from dragonStandRight.conf as individual text files.

Outputs
-------
<output_dir>/pclouds/<k>.txt          -- downsampled point cloud (N×3), k=0,24,...,336
<output_dir>/gt_poses/gt_pose_NN.txt  -- 7-element pose [tx,ty,tz,q0,q1,q2,q3], NN=00..14

Usage
-----
    python prepare_data.py \\
        --data_dir  /path/to/dragon_stand \\
        --output_dir /path/to/output \\
        [--voxel_size 0.005] \\
        [--overwrite]
"""

import argparse
import os
import sys

import numpy as np
import open3d as o3d

NUM_VIEWS = 15
DEGREE_STEP = 24


def parse_gt_poses(data_dir):
    """Parse dragonStandRight.conf → {0-based index: [tx,ty,tz,q0,q1,q2,q3]}."""
    conf = os.path.join(data_dir, 'dragonStandRight.conf')
    poses = {}
    for line in open(conf).read().splitlines()[2:]:
        words = line.split()
        basename = words[1][:-4]
        k = int(os.path.basename(basename).split('_')[1])
        poses[k // DEGREE_STEP] = [float(x) for x in words[2:]]
    return poses


def save_gt_poses(poses, gt_dir, overwrite=False):
    """Write per-frame GT pose files. Returns count of files written."""
    os.makedirs(gt_dir, exist_ok=True)
    n_written = 0
    for idx, param in poses.items():
        path = os.path.join(gt_dir, f'gt_pose_{idx:02d}.txt')
        if not overwrite and os.path.exists(path):
            continue
        np.savetxt(path, param, fmt='%.10f')
        n_written += 1
    return n_written


def main():
    parser = argparse.ArgumentParser(
        description='Prepare downsampled point clouds and GT poses from the dragon stand dataset.')
    parser.add_argument('--data_dir', required=True,
                        help='Dragon stand dataset root (PLY files + dragonStandRight.conf)')
    parser.add_argument('--output_dir', required=True,
                        help='Destination for pclouds/ and gt_poses/ subdirectories')
    parser.add_argument('--voxel_size', type=float, default=0.005,
                        help='Voxel down-sample size in metres (default: 0.005)')
    parser.add_argument('--overwrite', action='store_true',
                        help='Re-process frames even if output files already exist')
    args = parser.parse_args()

    pcloud_dir = os.path.join(args.output_dir, 'pclouds')
    gt_dir = os.path.join(args.output_dir, 'gt_poses')
    os.makedirs(pcloud_dir, exist_ok=True)

    # GT poses
    print(f'Parsing GT poses: {os.path.join(args.data_dir, "dragonStandRight.conf")}')
    poses = parse_gt_poses(args.data_dir)
    n_gt_written = save_gt_poses(poses, gt_dir, overwrite=args.overwrite)
    n_gt_skipped = len(poses) - n_gt_written
    print(f'GT poses: wrote={n_gt_written}  skipped={n_gt_skipped}  → {gt_dir}')

    # Point clouds
    n_saved = n_skipped = n_failed = 0
    for k in range(0, NUM_VIEWS * DEGREE_STEP, DEGREE_STEP):
        out_path = os.path.join(pcloud_dir, f'{k}.txt')
        if not args.overwrite and os.path.exists(out_path):
            n_skipped += 1
            continue
        ply = os.path.join(args.data_dir, f'dragonStandRight_{k}.ply')
        pcd = o3d.io.read_point_cloud(ply)
        if not pcd.has_points():
            print(f'  WARNING: could not read {ply}', file=sys.stderr)
            n_failed += 1
            continue
        pts = np.asarray(pcd.voxel_down_sample(voxel_size=args.voxel_size).points)
        np.savetxt(out_path, pts, fmt='%.6f')
        n_saved += 1
        print(f'  k={k:3d}  ({pts.shape[0]:5d} points)')

    print(f'\nDone.  saved={n_saved}  skipped={n_skipped}  failed={n_failed}')
    print(f'Output directory: {args.output_dir}')


if __name__ == '__main__':
    main()
