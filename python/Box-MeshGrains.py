import os
import sys
import pymef90

if __name__ == '__main__':
    import cubit
    from optparse import OptionParser
    debug=True
    ### 
    ### Parse input
    ###
    parser = OptionParser()
    parser.add_option("-l", dest="l", help="box size")
    parser.add_option("-m", "--meshsize", dest="meshsize", help="mesh size")
    parser.add_option("-o", "--output",      dest="genfile",    help="exodus mesh file")
    parser.add_option("-c", "--cub",             dest="cubfile",    help="cubit geometry file")
    parser.add_option("-n", "--ngrains",     dest="ngrains",    help= "number of grains")
    ###
    (options, args) = parser.parse_args()
    ###
    ngrains = eval(options.ngrains)
    ###
    l = eval(options.l)
    lz = 1.
    if not options.genfile == None:
        if options.meshsize == None:
            parser.error("must specify mesh size with -m or --meshsize")
        else:
            meshsize = eval(options.meshsize)
    genfile  = options.genfile
    cubfile  = options.cubfile
    ###
    cubit.init([""])
else:
    debug=True
    l    = 1.
    lz = 1.
    ngrains = 10.
    genfile = 'TEST.gen'
    cubfile = 'TEST.cub'
    meshsize = .1
    ###
cubit.cmd("reset")
cubit.cmd("set journal off")
cubit.cmd("set echo off")
###
### generate domain and grains
###
GRAINS_3D = pymef90.BoxCrystalCreate(l, ngrains)
###
### imprint and merge
###
cubit.cmd("imprint all")
cubit.cmd("merge all")
###
### create element blocks
###
index=1
for i, grain in enumerate(GRAINS_3D):
    if len(grain) > 0:
        for vol in grain:
            cubit.cmd("block %i surface with z_min>0 in volume %i" % (index, vol))
        index += 1
###
### save geometry
###
if not cubfile == None:
    cubit.cmd('save as "%s" overwrite' % cubfile)
###
### mesh
###
if not genfile == None:
    cubit.cmd("surface with z_min>0 size %f" %meshsize)
    cubit.cmd("Tridelaunay point placement gq")
    cubit.cmd("surface with z_min>0 scheme tridelaunay")
    cubit.cmd("mesh surface with z_min>0")
    cubit.cmd("surface with z_min>0 smooth scheme laplacian free")
    ###
    ### create nodesets
    ###
    eps=.005
    cubit.cmd("nodeset 1 node with z_coord > 0 and y_max <= %f and x_max > %f" % (-l/2. * (1.-eps), -l/2. ))
    if debug:
        print("### curves in nodeset 1:")
        cubit.cmd("list curve in nodeset 1 ids")
    cubit.cmd("nodeset 2 node with z_coord > 0. and y_min >= %f" % ( l/2. * (1.-eps)))
    if debug:
        print("### curves in nodeset 2:")
        cubit.cmd("list curve in nodeset 2 ids")
    cubit.cmd("nodeset 3 vertex with x_coord = %f and y_coord = %f and z_coord > 0." % (-l/2., -l/2.))
    if debug:
        print("### vertices in nodeset 3:")
        cubit.cmd("list vertex in nodeset 3 ids")
    cubit.cmd('export mesh "%s" dimension 2 overwrite' % genfile)
