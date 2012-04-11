#!/usr/bin/env python
# encoding: utf-8
"""
A few dictionaries which to the link between the type of computation and the basic file names 

 - directory
 - mesh file
 - args file

:author; psicsic
:date: 2011/04/11
"""

JOB_ID = {
    'bar':'bar{0.dim:d}d-EPS{1.epsilon:.3f}-h{0.meshsize:.4f}-TMAX{0.tmax:.3f}-NS{0.numsteps:d}-AT{1.atnum}',
    'PM':'PM-xc{0.xc:.1f}-theta{0.thetac}-TMAX{0.tmax:.1f}-NS{0.numsteps:d}',
    }

MESH_FILE = {
    'bar':"Bar{0.dim:d}D-Lx{0.length:.1f}-Ly{0.width:.2f}-mesh{0.meshsize:.4f}.gen",
    'PM':"PM-mesh{0.meshsize:.4f}.gen",
    }

ARGS_FILE = {
    'bar':"Bar{0.dim:d}D.args",
    'PM':"PM-Layers-Scaling.args",
    }


