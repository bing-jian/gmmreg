#!/usr/bin/env python
#coding=utf-8

##====================================================
## $Author: bing.jian $
## $Date: 2009-02-10 02:13:49 -0500 (Tue, 10 Feb 2009) $
## $Revision: 121 $
## $URL: http://gmmreg.googlecode.com/svn/trunk/Python/_demo.py $
##====================================================

import time
import subprocess
import ConfigParser

from numpy import loadtxt

import _core
import _plotting

def test(f_config, display = True):
    model,scene,after_tps = _core.run_ini(f_config)
    if display:
        _plotting.displayABC(model,scene,after_tps)


def run_executable(gmmreg_exe, f_config, method, display = True):
    cmd = '%s %s %s'%(gmmreg_exe, f_config, method)
    t1 = time.time()
    subprocess.call(cmd,shell=True)
    t2 = time.time()
    print "Elasped time is %s seconds"%(t2-t1)
    if display:
        display_pts(f_config)

def display_pts(f_config):
    c = ConfigParser.ConfigParser()
    c.read(f_config)
    section_common = 'FILES'
    mf = c.get(section_common,'model')
    sf = c.get(section_common,'scene')
    tf = c.get(section_common,'transformed_model')

    m = loadtxt(mf)
    s = loadtxt(sf)
    t = loadtxt(tf)
    _plotting.displayABC(m,s,t)




