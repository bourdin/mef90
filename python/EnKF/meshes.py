#!/usr/bin/env python
# encoding: utf-8


"""
:author: Paul Sicsic
:date: 2011/04/11
"""


import os, sys
try : 
    import cubit
except ImportError as err: 
    print "\n Unable to import cubit"
    print "Check PYTHONPATH"
    print "ImportError : {0} \n".format(err)

def cubit_init():
    """init commands for generating a file thourgh python commands"""
    cubit.init([""])
    cubit.cmd('reset')
    cubit.cmd('set journal on')
    cubit.cmd('set echo off')

def geom_bar_2d(lx, ly):
    """ define the size of the 2d bar """
    cubit.cmd("create surface rectangle width %f height %f zplane"%(lx, ly))    
    return cubit.get_last_id("surface")


def geom_bar_3d(lx, ly, lz):
    """ define size of the 3d bar"""
    cubit.cmd("brick x %f y %f z %f"%(lx, ly,  lz))  
    return cubit.get_last_id("volume")


def bar2d(lx, ly, meshsize, meshfile="bar_2d_python.gen"):
    """ 
    Generate a file for a 2d bar in exodusII format
    :param lx: length of the bar : uniaxial direcion in which the force
    is applied
    :type lx: float

    """

    cubfile = meshfile.split('.gen')[0]+'.cub'
    cubit_init()

    s1id = geom_bar_2d(lx, ly)      
    cubit.cmd('surface %i scheme tridelaunay'%s1id)
    cubit.cmd('mesh surface %i'%s1id)

    cubit.cmd('block 1 surface %i'%s1id)
#Rque plus robuste utilise la poistion des noeuds .... 
#    cubit.cmd('nodeset 1 node with x_coord = %f'%(-lx/2))                                                                                                                               
#    cubit.cmd('nodeset 2 node with x_coord = %f'%(lx/2))   

    cubit.cmd('nodeset 1 curve 2')
    cubit.cmd('nodeset 2 curve 4') 

    cubit.cmd('export mesh "%s" dimension 2 overwrite'%meshfile)
    cubit.cmd('save as "%s" overwrite'%cubfile)


def bar3d(lx, ly, lz, meshsize, meshfile="bar_3d_python.gen"):
    """ 
    Generate a file for 3d 
    
    :param lx: length of the bar : uniaxial direcion in which the force
    is applied
    :type lx: float
    """
    cubfile = meshfile.split('.gen')[0]+'.cub'
    cubit_init()

    vol_id  = geom_bar_3d(lx, ly, lz)
    cubit.cmd('volume %i size %f scheme Tetmesh'%(vol_id, meshsize))
    cubit.cmd('imprint all')
#    cubit.cmd('merge all')

    cubit.cmd('mesh volume %i'%vol_id)

    cubit.cmd('block 1 volume %i'%vol_id)
#Rque plus robuste utilise la poistion des noeuds .... 
#    cubit.cmd('nodeset 1 node with x_coord = %f'%(-lx/2))                                                                                                                               
#    cubit.cmd('nodeset 2  node with x_coord = %f'%(lx/2))   
    cubit.cmd('nodeset 1 surface 4')
    cubit.cmd('nodeset 2 surface 6')

    cubit.cmd('export mesh "%s" dimension 3 overwrite'%meshfile)
    cubit.cmd('save as "%s" overwrite'%cubfile)
    
def mesh_selection(options):
    """
    A function  generating from computation options
    """
    #from pymef.mesh.bar import bar2d, bar3d
    #from pymef.mesh.pied import pied2d, pied3d_real
    #from pymef.mesh.fiber_matrix import fiber2d
    #from pymef.mesh.bending import bending_3pts_2d
    #from pymef.mesh.bending import bending_4pts_2d, bending_sen_2d
    #from pymef.mesh.masonry import french_pantheon 
    #from pymef.mesh.therm import therm2d
    #from pymef.mesh.pm_inclusions import PM2d
    #from pymef.mesh.rings import peach, ring
    if options.type == 'bar' :
        if options.dim == 2:
            bar2d(options.length, options.width, options.meshsize, options.meshfile)
        elif options.dim == 3:
            bar3d(options.length, options.width, options.height,
                    options.meshsize, options.meshfile)
        else:
            raise ValueError("Dimension can only be 2 or 3")
    
    elif options.type == 'therm':
        therm2d(options.length, options.width, options.meshsize, options.meshfile)
    elif options.type == 'PM':
        PM2d(options.r, options.lx, options.ly, options.lz,
                 options.theta1, options.theta2, options.alpha,
                 options.thetac, options.xc, options.meshsize, options.meshfile)
    else :
        print options.type 
        raise ValueError("The job type you have asked for is not implemented \
                         pymef only handles : bar, pied, fiber, 3pbt, PM")






