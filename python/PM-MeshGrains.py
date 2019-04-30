import os
import sys
import pymef90

if __name__ == '__main__':
    import cubit
    from optparse import OptionParser
    import ConfigParser
    debug=True
    ### 
    ### Parse input
    ###
    parser = OptionParser()
    parser.add_option("-g", "--geometry",    dest="geometry",    help="geometry cfg file")
    parser.add_option("-m", "--meshsize1", dest="meshsize1", help="mesh size of the center coin")
    parser.add_option("-M", "--meshsize2", dest="meshsize2", help="mesh size of the outer shell")
    parser.add_option("-o", "--output",      dest="genfile",    help="exodus mesh file")
    parser.add_option("-c", "--cub",             dest="cubfile",    help="cubit geometry file")
    parser.add_option("-n", "--ngrains",     dest="ngrains",    help= "number of grains")
    ###
    (options, args) = parser.parse_args()
    ###
    if options.geometry == None:
        parser.error("must specify geometry cfg file with -g or --geometry")
    if not options.geometry[0] == "/":
            options.geometry = os.getcwd() + '/' + options.geometry
    if options.meshsize1 == None:
        options.meshsize1 = .1
    if options.meshsize2 == None:
        options.meshsize2 = options.meshsize1
    ###
    CBconfig = ConfigParser.ConfigParser()
    CBconfig.read(options.geometry)
    ###
    R                = eval(CBconfig.get('geometry', 'ro'))
    r                = eval(CBconfig.get('geometry', 'ri'))
    lz           = eval(CBconfig.get('geometry', 'lz'))
    thetac   = eval(CBconfig.get('geometry', 'thetac'))
    xc           = eval(CBconfig.get('geometry', 'xc'))
    if options.ngrains == None:
        ngrains = eval(CBconfig.get('geometry', 'ngrains'))
    else:
        ngrains = eval(options.ngrains)
    ###
    scheme2d = CBconfig.get('mesh', 'scheme2d')
    ###
    meshsize1 = eval(options.meshsize1)
    meshsize2 = eval(options.meshsize2)
    genfile  = options.genfile
    cubfile  = options.cubfile
    ###
    cubit.init([""])
else:
    debug=False
    R    = 5.
    r = 2.
    lz = 1
    thetac = 15.
    xc       = -.1
    ngrains = 10.
    genfile = 'TEST.gen'
    cubfile = 'TEST.cub'
    meshsize1 = .1
    meshsize2 = .5
    scheme2d = 'tridelaunay'
    ###
cubit.cmd("reset")
cubit.cmd("set journal off")
cubit.cmd("set echo off")
###
### generate domain and grains
###
cubit.cmd("Graphics Pause")
(OUTSIDE_3D, GRAINS_3D) = pymef90.PacmanCoinCrystalCreate(R, r, lz, thetac, xc, ngrains)
cubit.cmd("Display")
###
### imprint and merge
###
#for i in range(len(GRAINS_3D)):
#    cubit.cmd("imprint volume in OUTSIDE_3D with volume in Grain%i" % (i-1))
cubit.cmd("imprint all")
cubit.cmd("merge all")
###
### find edges of the domain and create nodesets
###
allcurves = []
allcurves.extend(cubit.get_relatives("volume", OUTSIDE_3D, "curve"))
for grain in GRAINS_3D:
    for vol in grain:
        allcurves.extend(cubit.get_relatives("volume", vol, "curve"))
alledges = [c for c in allcurves if len(cubit.get_relatives("curve", c, "volume"))==1]
for e in alledges:
    cubit.cmd('group "ALLEDGES" add curve %i'%e)
cubit.cmd("nodeset 1 curve with z_min>0 in ALLEDGES")
cubit.cmd("delete ALLEDGES")
cubit.cmd("nodeset 1 remove curve with z_min>0 and x_max=%f"%R)
cubit.cmd("nodeset 2 curve with z_min>0 and x_max=%f"%R)
###
### create element blocks
###
index=1
for i, grain in enumerate(GRAINS_3D):
    if len(grain) > 0:
        for vol in grain:
            cubit.cmd("block %i surface with z_min>0 in volume %i" % (index, vol))
        index += 1
cubit.cmd("block %i surface with z_min>0 in volume %i" % (index, OUTSIDE_3D))
###
### save geometry
###
if not cubfile == None:
    cubit.cmd('save as "%s" overwrite' % cubfile)
###
### mesh
###
if not genfile == None:
    index=1
    cubit.cmd("Tridelaunay point placement gq")
    for i, grain in enumerate(GRAINS_3D):
        if len(grain) > 0:
            for vol in grain:
                cubit.cmd("surface with z_min>0 in volume %i size %f" % (vol, meshsize1))
                cubit.cmd("surface with z_min>0 in volume %i scheme %s" % (vol, scheme2d))
                cubit.cmd("mesh surface with z_min>0 in volume %i" % vol)
            index += 1    
    cubit.cmd('surface with z_min>0 in volume in OUTSIDE_3D name "OUTSIDE_2D"')
    cubit.cmd('group "CoinEdge" add curve with x_min > %f in OUTSIDE_2D' % (-(r+R)/2.))
    cubit.cmd("curve with z_min>0 and x_max=%f size %f" % (R, meshsize2))
    cubit.cmd("OUTSIDE_2D sizing function type bias start curve in CoinEdge finish curve with z_min>0 and x_max=%f" % R)
    cubit.cmd("OUTSIDE_2D scheme %s" % scheme2d)
    cubit.cmd("mesh OUTSIDE_2D")
    #cubit.cmd("mesh surface %i"%(ngrains+1))
    cubit.cmd('delete CoinEdge')
    #
    cubit.cmd("surface with z_min>0 smooth scheme laplacian free")
    cubit.cmd("smooth surface with z_min>0")
    cubit.cmd('export mesh "%s" dimension 2 overwrite' % genfile)
