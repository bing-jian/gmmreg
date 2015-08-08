#!/usr/bin/env python
#coding=utf-8

import os, numpy, subprocess, time
from Quaternion import Quat, normalize

gmmreg_binary = {'nt' : r'gmmreg_demo.exe',
                 'posix': r'gmmreg_demo'}

data_path = '../data/dragon_stand'
tmp_path = './tmp'
config_file = './dragon_stand.ini'

if not os.path.exists(tmp_path):
    os.makedirs(tmp_path)

def run_rigid(gmmreg_exe, model, scene, f_config, method = 'rigid'):
    if type(model) == type('abc'):
        model = numpy.loadtxt(model)
        scene = numpy.loadtxt(scene)
    if model.shape[0] > 10000:
        model = model[0::10]
        scene = scene[0::10]
    numpy.savetxt(os.path.join(tmp_path, 'model.txt'), model)
    numpy.savetxt(os.path.join(tmp_path, 'scene.txt'), scene)
    cmd = '%s %s %s'%(gmmreg_exe, f_config, 'rigid')
    t1 = time.time()
    subprocess.call(cmd, shell=True)
    t2 = time.time()
    print "Run time : %s seconds" % (t2-t1)
    param = numpy.loadtxt(os.path.join(tmp_path, 'final_rigid.txt'))
    return param, t2-t1

def ply2txt(plyfile):
    return '%s.txt' % os.path.splitext(plyfile)[0]

def pairwise_rigid3d(model_ply, scene_ply, f_config):
    model_txt = ply2txt(model_ply)
    scene_txt = ply2txt(scene_ply)
    print model_txt, scene_txt
    binary = os.path.join(os.getcwd(), gmmreg_binary[os.name])
    return run_rigid(binary, model_txt, scene_txt, f_config)

def ground_truth(model_ply, scene_ply, gt_pos):
    pos_i = gt_pos[os.path.basename(model_ply)]
    pos_j = gt_pos[os.path.basename(scene_ply)]
    q_i = Quat(normalize(numpy.array(pos_i[3::])))
    q_j = Quat(normalize(numpy.array(pos_j[3::])))
    q_ji = q_i / q_j
    q_ij = q_j / q_i
    return q_ij, q_ji

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

def get_all_plyfiles(data_path):
    ind = range(0, 360, 24)
    return [os.path.join(data_path,
                         'dragonStandRight_%d.ply'%j) for j in ind]

def run_rigid_batch(ply_files, step, config_file, res):
    for i in range(15):
        j = (i + step) % 15
        model_ply = ply_files[i]
        scene_ply = ply_files[j]
        print model_ply, scene_ply
        res[(i, j)] = pairwise_rigid3d(model_ply, scene_ply, config_file)
    return res

def compare_rotation(ply_files, reg_result, gt_pos):
    res = {}
    for i, j in reg_result:
        model_ply = ply_files[i]
        scene_ply = ply_files[j]
        q_ij, q_ji = ground_truth(model_ply, scene_ply, gt_pos)
        param = reg_result[(i,j)][0]
        q = Quat(normalize(param[0:4]))
        res[(i, j)] = abs(numpy.dot(q.q, q_ji.q))
    return res

def run_single_step(step, reg_result):
    ply_files = get_all_plyfiles(data_path)
    reg_result = run_rigid_batch(ply_files, step, config_file, reg_result)
    return reg_result

def run_all_steps():
    reg_result = {}
    for step in [-1, 1, -2, 2, -3, 3]:
        reg_result = run_single_step(step, reg_result)
    return reg_result

def compute_accuracy(reg_result):
    gt_pos = load_dragon_conf(data_path)
    ply_files = get_all_plyfiles(data_path)
    accuracy = compare_rotation(ply_files, reg_result, gt_pos)
    return accuracy
