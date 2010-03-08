#!/usr/bin/env python
import cubit
import numpy as np
import pymef90

cubit.init([""])
cubit.cmd("reset")
###
###  DCB Basic Geometry
###
lx=8.
ly1=4
ly=7.
lz=1.
epscrack=1/16.
lcrack=lx/2.
###
### Grips
###
ri=1/4.
xi=.5
yi=1/2.
###
### Layers
###
loffset=0.1
llayers=3.
bb=[lx/2.-llayers/2., lx/2.+llayers/2.+loffset*2., -ly/2., ly/2., -lz/2., lz/2.]
theta1=.3
theta2=.5
alpha=15
depth=lz/3.
###
### Create the DCB shape
###
DCB_3D = pymef90.DCBCreate(lx, ly1, ly, lz, epscrack, lcrack)
### 
### create the top grip
###
cubit.cmd('create cylinder radius %f height %f' % (ri, lz*2.))
tmp_ID=cubit.get_last_id("volume")
cubit.cmd('move volume %i X %f Y %f' % (tmp_ID, xi, yi))
(DCB_3D, GRIP1_3D) = WebcutTool(DCB_3D, tmp_ID, delete=True)
### 
### create the bottom grip
###
cubit.cmd('create cylinder radius %f height %f' % (ri, lz*2.))
tmp_ID=cubit.get_last_id("volume")
cubit.cmd('move volume %i X %f Y %f' % (tmp_ID, xi, -yi))
(DCB_3D, Grip2_3D) = WebcutTool(DCB_3D, tmp_ID, delete=True)
###
### Reserve the center section for layers
###
cubit.cmd('create brick X %f Y %f Z %f' % (llayers, ly, lz*2.))
tmp_ID=cubit.get_last_id("volume")
cubit.cmd('move volume %i X %f' % (tmp_ID, lx/2.))
(DCB_3D, LAYER1_3D) = pymef90.WebcutTool(DCB_3D, tmp_ID, delete=True)
### 
### Layer
###
(LAYER1_3D, LAYER2_3D) = pymef90.MilledLayer(LAYER1_3D, bb, alpha, theta1, theta2, depth)