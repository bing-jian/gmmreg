#!/usr/bin/env python
#coding=utf-8

##==============================================================================
## $Author: bingjian $
## $Date: 2013-01-04 01:39:25 -0500 (Fri, 04 Jan 2013) $
## $Revision: 145 $
## $URL: http://gmmreg.googlecode.com/svn/trunk/Python/_core.py $
##==============================================================================

import ConfigParser
import time
from math import cos,sin,log,exp,sqrt
from numpy import loadtxt,arange,array,dot,delete,reshape,kron,eye,ones,zeros,trace,s_,r_,c_,squeeze
from numpy.linalg import svd,qr,norm
from scipy.optimize import fmin_bfgs, fmin_l_bfgs_b
#from _extension import *
#import _plotting


def normalize(x):
    """
    Translate and scale a point set to make it have zero mean and unit variance. Return the point set after normalization, its original centroid and scale.
    """
    centroid = x.mean(0)
    x = x - centroid
    scale = norm(x,'fro')/sqrt(x.shape[0])
    x = x/scale
    return x,centroid,scale


def denormalize(x,centroid,scale):
    """Denormalize a point set from saved centroid and scale."""
    x = x*scale + centroid
    return x


def L2_distance(model, scene, scale):
    """
    Compute the L2 distance between the two Gaussian mixture densities constructed from a moving 'model' point set and a fixed 'scene' point set at a given 'scale'. The term that only involves the fixed 'scene' is excluded from the returned distance.  The gradient with respect to the 'model' is calculated and returned as well.
    """
    f1, g1 = gauss_transform(model, model, scale)
    f2, g2 = gauss_transform(model, scene, scale)
    # print f1,f2
    f =  f1 - 2*f2
    g = 2*array(g1) - 2*array(g2)
    return f,g


def correlation(model, scene, scale):
    """
    Compute the normalized correlation between the two Gaussian mixture densities constructed from a moving 'model' point set and a fixed 'scene' point set at a given 'scale'. The term that only involves the fixed 'scene' is excluded from the returned distance.  The gradient with respect to the 'model' is calculated and returned as well.
    """
    f1, g1 = gauss_transform(model, model, scale)
    f2, g2 = gauss_transform(model, scene, scale)
    #f = -f2/sqrt(f1)
    #g = array(g1)*f2/(f1*sqrt(f1)) - array(g2)/sqrt(f1)
    #below is the version using -log(correlation), supposely slightly faster because no square root
    #print f1,f2
    #epsilon = 1e-10
    #f = log(f1+epsilon) - 2*log(f2+epsilon)
    #g = 2*array(g1)/(f1+epsilon) - 2*array(g2)/(f2+epsilon)
    #squared correlation
    f = -f2*f2/f1
    g = 2*array(g1)*f2*f2/(f1*f1) - 2*array(g2)*f2/f1
    return f,g


def init_param(n, d, opt_affine=True):
    """
    Construct the initial parameters of thin-plate splines.
    """
    init_tps = [0.0]*(d*n-d*(d+1))
    # init_tps is the vectorized version of the (n-d-1)-by-d coefficient matrix
    init_affine = ([0.0]*d+[1.0])*d
    # init_affine is the vectorized version of the (d+1)-by-d affine matrix
    if opt_affine:
        # init_affine is part of parameter and will be updated
        init_param = init_affine + init_tps
    else:
        # init_affine is fixed during the optimization
        init_param = init_tps
    return array(init_param)
    #return init_param


def transform_points(param, basis):
    """
    Perform the thin-plate spline transform.
    """
    (nL,n) = basis.shape
    d = param.shape[0]/n
    affine_param = param[0:d*(d+1)].reshape(d+1,d)
    tps_param = param[d*(d+1):d*n].reshape(n-d-1,d)
    after_tps = dot(basis,r_[affine_param,tps_param])
    return after_tps

def compute_GRBF(ctrl_pts, landmarks, sigma):
    """
    Compute the kernel matrix for GRBF splines.
    """
    def kernel_func(r, sigma):
        return exp(-r*r/(sigma*sigma))

    K = array([kernel_func(norm(x-y), sigma) for x in ctrl_pts for y in ctrl_pts])
    K = K.reshape(len(ctrl_pts),len(ctrl_pts))
    U = array([kernel_func(norm(x-y), sigma) for x in landmarks for y in ctrl_pts])
    U = U.reshape(len(landmarks),len(ctrl_pts))
    return K,U


def compute_TPS_K(ctrl_pts, landmarks = None, _lambda = 0):
    """
    Compute the kernel matrix for thin-plate splines.
    Reference:
      Landmark-based Image Analysis, Karl Rohr, p195
    """

    #kernel_func = [lambda r,_lambda=0: 0 if r==0 else r*r*log(r), lambda r,_lambda=0: -r]
    #the above definition is not used because the if else syntax is not supported in Python2.4
    def kernel_func_2d(r, _lambda=0):
        #_lambda reserved for regularization
        if r == 0:
            return 0
        else:
            return r*r*log(r)
    def kernel_func_3d(r, _lambda=0):
        #_lambda reserved for regularization
        return -r
    kernel_func = (kernel_func_2d, kernel_func_3d)

    [n,d] = ctrl_pts.shape
    K = [kernel_func[d-2](norm(ctrl_pts[i]-ctrl_pts[j]), _lambda) for i in arange(n) for j in arange(n)]
    K = array(K).reshape(n,n)
    if landmarks is not None:
        [m,d] = landmarks.shape  # assert (d,d) equal
        U = [kernel_func[d-2](norm(landmarks[i]-ctrl_pts[j]), _lambda) for i in arange(m) for j in arange(n)]
        U = array(U).reshape(m,n)
    else:
        U = None
    return K,U

def prepare_TPS_basis(landmarks,ctrl_pts):
    """
    Return the basis for performing TPS transforms and the kernel matrix for computing the bending energy.
    """
    [m,d] = landmarks.shape
    [n,d] = ctrl_pts.shape
    [K,U] = compute_TPS_K(ctrl_pts,landmarks)
    Pm = c_[ones((m,1)),landmarks]
    Pn = c_[ones((n,1)),ctrl_pts]
    u,s,vh = svd(Pn)
    PP = u[:,d+1:]
    TPS_basis = c_[Pm,dot(U,PP)]
    TPS_kernel = dot(PP.T,dot(K,PP))
    return TPS_basis,TPS_kernel


def obj_TPS(dist_func, param, basis, kernel, scene, scale, _lambda): #, init_affine=None):
    """
    The cost function based on TPS model
    Todo: enable the option to fix the affine parameter
    """
    nL,n = basis.shape     # (control-pts, landmarks)
    d = scene.shape[1]
    affine_param = param[0:d*(d+1)].reshape(d+1,d)
    tps_param = param[d*(d+1):d*n].reshape(n-d-1,d)
    after_tps = dot(basis,r_[affine_param,tps_param])
    bending = trace(dot(tps_param.T,dot(kernel,tps_param)))
    distance, grad = dist_func(after_tps, scene, scale)
    energy = distance + _lambda * bending
    grad = dot(basis.T, grad)
    grad[d+1:n] += 2*_lambda*dot(kernel,tps_param)
    grad = grad.reshape(d*n)
    return energy, grad

def obj_L2_TPS(param, basis, kernel, scene, scale, _lambda):
    #return obj_TPS(L2_distance, param, basis, kernel, scene, scale, alpha, beta)
    nL,n = basis.shape     # (control-pts, landmarks)
    d = scene.shape[1]
    affine_param = param[0:d*(d+1)].reshape(d+1,d)
    tps_param = param[d*(d+1):d*n].reshape(n-d-1,d)
    after_tps = dot(basis,r_[affine_param,tps_param])
    bending = trace(dot(tps_param.T,dot(kernel,tps_param)))
    distance, grad = L2_distance(after_tps, scene, scale)
    energy = distance + _lambda * bending
    grad = dot(basis.T, grad)
    grad[d+1:n] += 2*_lambda*dot(kernel,tps_param)
    grad = grad.reshape(d*n)
    return energy, grad


def obj_KC_TPS(param, basis, kernel, scene, scale, alpha, beta):
    #return obj_TPS(correlation, param, basis, kernel, scene, scale, alpha, beta)
    nL,n = basis.shape     # (control-pts, landmarks)
    d = scene.shape[1]
    affine_param = param[0:d*(d+1)].reshape(d+1,d)
    tps_param = param[d*(d+1):d*n].reshape(n-d-1,d)
    after_tps = dot(basis,r_[affine_param,tps_param])
    bending = trace(dot(tps_param.T,dot(kernel,tps_param)))
    distance, grad = correlation(after_tps, scene, scale)
    energy = alpha*distance + beta * bending
    grad = alpha*dot(basis.T, grad)
    grad[d+1:n] += 2*beta*dot(kernel,tps_param)
    grad = grad.reshape(d*n)
    return energy, grad


def run_multi_level(model,scene,ctrl_pts,level,scales,lambdas,iters):
    [n,d] = ctrl_pts.shape
    x0 = init_param(n,d)
    [basis, kernel] = prepare_TPS_basis(model, ctrl_pts)
    for i in range(level):
        x = fmin_l_bfgs_b(obj_L2_TPS, x0, None, args=(basis,kernel,scene,scales[i],lambdas[i]),maxfun=iters[i])
        x0 = x[0]
    after_tps = transform_points(x0,basis)
    return after_tps


def run_ini(f_config):

    section_common = 'FILES'
    section_option = 'GMMREG_OPT'

    c = ConfigParser.ConfigParser()
    c.read(f_config)
    model_file = c.get(section_common,'model')
    scene_file = c.get(section_common,'scene')
    model = loadtxt(model_file)
    scene = loadtxt(scene_file)
    try:
        ctrl_pts_file = c.get(section_common,'ctrl_pts')
        ctrl_pts = loadtxt(ctrl_pts_file)
    except:
        ctrl_pts = model
    level = int(c.get(section_option,'level'))
    option_str = c.get(section_option,'sigma')
    scales = [float(s) for s in option_str.split(' ')]
    option_str = c.get(section_option,'lambda')
    lambdas = [float(s) for s in option_str.split(' ')]

    option_str = c.get(section_option,'max_function_evals')
    iters = [int(s) for s in option_str.split(' ')]

    normalize_flag = int(c.get(section_option,'normalize'))
    #print normalize_flag
    if normalize_flag==1:
        [model, c_m, s_m] = normalize(model)
        [scene, c_s, s_s] = normalize(scene)
        [ctrl_pts, c_c, s_c] = normalize(ctrl_pts)
    t1 = time.time()
    after_tps = run_multi_level(model,scene,ctrl_pts,level,scales,lambdas,iters)
    if normalize_flag==1:
        model = denormalize(model,c_m,s_m)
        scene = denormalize(scene,c_s,s_s)
        after_tps = denormalize(after_tps,c_s,s_s)
    t2 = time.time()
    print "Elasped time is %s seconds"%(t2-t1)
    #_plotting.displayABC(model,scene,after_tps)
    #_plotting.display2Dpointsets(after_tps,scene)
    return model,scene,after_tps

