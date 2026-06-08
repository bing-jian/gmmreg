#!/usr/bin/env python3
"""Step 1: extract per-frame point clouds and ground-truth poses.

Reads the lounge RGB-D dataset, converts each depth image to a 3-D point
cloud, voxel-downsamples it, and writes everything to a structured output
directory that later steps can consume without any dataset dependency.

Outputs
-------
<output_dir>/pclouds/pcloud_XXXXX.txt   -- downsampled point cloud (N×3)
<output_dir>/gt_poses/gt_pose_XXXXX.txt -- 4×4 absolute pose matrix

Frame indices are 1-based (XXXXX = 00001 … 03000) to match depth filenames.

Usage
-----
    python prepare_data.py \\
        --data_dir /path/to/lounge \\
        --output_dir /path/to/output \\
        [--voxel_size 0.065] \\
        [--overwrite]
"""

import argparse
import os
import sys

import cv2
import numpy as np
import open3d as o3d


# Lounge dataset camera intrinsics (http://qianyi.info/scenedata.html).
_FX, _FY = 525.0, 525.0
_CX, _CY = 319.5, 239.5

NUM_FRAMES = 3000


# ---------------------------------------------------------------------------
# Core helpers
# ---------------------------------------------------------------------------

def depth_to_pointcloud(depth_img: np.ndarray) -> np.ndarray:
    """Convert a 16-bit depth image to an (N, 3) float64 point cloud."""
    ydim, xdim = depth_img.shape
    depth_flat = depth_img.flatten().astype(np.float64)
    valid = np.where(depth_flat > 0)[0]
    u, v = np.meshgrid(np.arange(xdim), np.arange(ydim))
    u = u.flatten()[valid]
    v = v.flatten()[valid]
    z = depth_flat[valid] / 1000.0          # mm → metres
    x = (u - _CX) * z / _FX
    y = (v - _CY) * z / _FY
    return np.column_stack([x, y, z])


def parse_trajectory(trajectory_path: str) -> dict:
    """Return {0-based frame index: 4×4 pose matrix} for all NUM_FRAMES frames."""
    lines = open(trajectory_path).readlines()
    poses = {}
    for i in range(NUM_FRAMES):
        block = ''.join(lines[5 * i + 1: 5 * i + 5])
        poses[i] = np.fromstring(block, dtype=float, sep='\t').reshape(4, 4)
    return poses


def downsample(pts: np.ndarray, voxel_size: float) -> np.ndarray:
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(pts)
    return np.asarray(pcd.voxel_down_sample(voxel_size=voxel_size).points)


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def save_gt_poses(poses: dict, gt_dir: str, overwrite: bool = False) -> int:
    """Write per-frame GT pose files. Returns count of files written."""
    os.makedirs(gt_dir, exist_ok=True)
    n_written = 0
    for idx, mat in poses.items():
        path = os.path.join(gt_dir, f'gt_pose_{idx + 1:05d}.txt')
        if not overwrite and os.path.exists(path):
            continue
        np.savetxt(path, mat, fmt='%.10f')
        n_written += 1
    return n_written


def process_frame(depth_path: str, voxel_size: float):
    """Load depth image, convert and downsample. Returns (N,3) array or None."""
    depth = cv2.imread(depth_path, cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)
    if depth is None:
        return None
    pts = depth_to_pointcloud(depth)
    return downsample(pts, voxel_size)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description='Prepare per-frame point clouds and GT poses from the lounge RGB-D dataset.')
    parser.add_argument('--data_dir', required=True,
                        help='Lounge dataset root (contains depth/ and lounge_trajectory.log)')
    parser.add_argument('--output_dir', required=True,
                        help='Destination for pclouds/ and gt_poses/ subdirectories')
    parser.add_argument('--voxel_size', type=float, default=0.065,
                        help='Voxel down-sample size in metres (default: 0.065)')
    parser.add_argument('--overwrite', action='store_true',
                        help='Re-process frames even if output files already exist')
    args = parser.parse_args()

    pcloud_dir = os.path.join(args.output_dir, 'pclouds')
    gt_dir     = os.path.join(args.output_dir, 'gt_poses')
    os.makedirs(pcloud_dir, exist_ok=True)

    # --- Ground-truth poses (parse once, write all) -------------------------
    trajectory_path = os.path.join(args.data_dir, 'lounge_trajectory.log')
    print(f'Parsing ground-truth trajectory: {trajectory_path}')
    poses = parse_trajectory(trajectory_path)
    n_gt_written = save_gt_poses(poses, gt_dir, overwrite=args.overwrite)
    n_gt_skipped = len(poses) - n_gt_written
    print(f'GT poses: wrote={n_gt_written}  skipped={n_gt_skipped}  → {gt_dir}')

    # --- Point clouds -------------------------------------------------------
    depth_dir = os.path.join(args.data_dir, 'depth')
    n_saved = n_skipped = n_failed = 0

    for i in range(NUM_FRAMES):
        frame_no = i + 1                    # 1-based, matches depth filenames
        out_path = os.path.join(pcloud_dir, f'pcloud_{frame_no:05d}.txt')

        if not args.overwrite and os.path.exists(out_path):
            n_skipped += 1
            continue

        depth_path = os.path.join(depth_dir, f'{frame_no:06d}.png')
        pts = process_frame(depth_path, args.voxel_size)
        if pts is None:
            print(f'  WARNING: could not read {depth_path}', file=sys.stderr)
            n_failed += 1
            continue

        np.savetxt(out_path, pts, fmt='%.6f')
        n_saved += 1

        if frame_no % 100 == 0:
            print(f'  {frame_no:4d}/{NUM_FRAMES}  '
                  f'({pts.shape[0]:5d} points after voxel_size={args.voxel_size})')

    print(f'\nDone.  saved={n_saved}  skipped={n_skipped}  failed={n_failed}')
    print(f'Output directory: {args.output_dir}')


if __name__ == '__main__':
    main()
