#!/usr/bin/env python
import cubit
import numpy as np
import pymef90
import os
from optparse import OptionParser
import ConfigParser

### 
### Parse input
###
parser = OptionParser()
parser.add_option("-g", "--geometry", dest="geometry", help="geometry cfg file")
parser.add_option("-m", "--meshsize", dest="meshsize", help="mesh size")
parser.add_option("-o", "--output", dest="genfile", help="exodus mesh file")

(options, args) = parser.parse_args()
if options.geometry == None:
  parser.error("must specify geometry cfg file with -g or --geometry")
if not options.geometry[0] == "/":
    options.geometry = os.getcwd() + '/' + options.geometry
if options.meshsize == None:
  parser.error("must specify mesh size with -m or --meshsize")
if options.genfile == None:
  parser.error("must specify output mesh file with -o or --output")

DCBconfig = ConfigParser.ConfigParser()
DCBconfig.read(options.geometry)

lx       = eval(DCBconfig.get('geometry', 'lx'))
ly1      = eval(DCBconfig.get('geometry', 'ly1'))
ly       = eval(DCBconfig.get('geometry', 'ly'))
lz       = eval(DCBconfig.get('geometry', 'lz'))
epscrack = eval(DCBconfig.get('geometry', 'epscrack'))
lcrack   = eval(DCBconfig.get('geometry', 'lcrack'))
ri       = eval(DCBconfig.get('geometry', 'ri'))
xi       = eval(DCBconfig.get('geometry', 'xi'))
yi       = eval(DCBconfig.get('geometry', 'yi'))

loffset  = eval(DCBconfig.get('geometry', 'loffset'))
llayers  = eval(DCBconfig.get('geometry', 'llayers'))
theta1   = eval(DCBconfig.get('geometry', 'theta1'))
theta2   = eval(DCBconfig.get('geometry', 'theta2'))
alpha    = eval(DCBconfig.get('geometry', 'alpha'))
bb       = eval(DCBconfig.get('geometry', 'bb'))
depth    = eval(DCBconfig.get('geometry', 'depth'))
scheme3d   = DCBconfig.get('mesh', 'scheme3d')

cubit.init([""])
cubit.cmd("reset")
cubit.cmd("set journal off")
cubit.cmd("set echo off")

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
(DCB_3D, GRIP1_3D) = pymef90.WebcutTool(DCB_3D, tmp_ID, delete=True)
### 
### create the bottom grip
###
cubit.cmd('create cylinder radius %f height %f' % (ri, lz*2.))
tmp_ID=cubit.get_last_id("volume")
cubit.cmd('move volume %i X %f Y %f' % (tmp_ID, xi, -yi))
(DCB_3D, GRIP2_3D) = pymef90.WebcutTool(DCB_3D, tmp_ID, delete=True)
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
(LAYER_3D) = pymef90.MilledLayer(LAYER1_3D, bb, alpha, theta1, theta2, depth)

###
### imprint and merge
###
cubit.cmd("imprint all")
cubit.cmd("merge all")
###
### create groups
###
pymef90.GroupAddVolList("GRIP1_3D", GRIP1_3D)
pymef90.GroupAddVolList("GRIP2_3D", GRIP2_3D)
pymef90.GroupAddVolList("DCB_3D", DCB_3D)
pymef90.GroupAddVolList("LAYER1_3D", LAYER_3D)
###
### mesh
###
cubit.cmd("volume all size %f" % float(options.meshsize))
cubit.cmd("volume all scheme %s" % scheme3d)
cubit.cmd("mesh volume all")
###
### BC
###
cubit.cmd('nodeset 1 volume in GRIP1_3D')
cubit.cmd("nodeset 2 volume in GRIP2_3D")
###
cubit.cmd("block 1 volume in GRIP1_3D")
cubit.cmd("block 2 volume in GRIP2_3D")
cubit.cmd("block 3 volume in DCB_3D")
cubit.cmd("block 4 volume in LAYER1_3D")
###
### export mesh
###
cubit.cmd('export mesh "%s" dimension 3 overwrite' % options.genfile)
