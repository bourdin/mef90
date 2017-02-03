#!/usr/bin/env python

import exodus as exo
import sys
import os.path
import numpy as np

f = sys.argv[-1]
if not os.path.isfile(f) or not (len(sys.argv) == 2):
    print "usage: {0} <filename>".format(sys.argv[0])
    exit
else:
    fout = f.split(".gen")[0]+"_fixed.gen"
    if os.path.exists(fout):
        print "{0} exists. Delete and run again".format(fout)
        exit

e = exo.exodus(f,"r")
dim       = e.num_dimensions()
nVertices = e.num_nodes()
nCells    = e.num_elems()
nCSets    = e.num_blks()
nVSets    = e.num_node_sets()

title     = "{0} fixed by {1}".format(e.title(),sys.argv[0])

eout = exo.exodus(fout,'w',numDims=dim,title=title,numNodes=nVertices-4,numElems=nCells,numBlocks=nCSets,numNodeSets=nVSets,numSideSets=0)

print "number of dimensions: {0}".format(dim)
print "number of vertices: {0}".format(nVertices)
print "number of cells: {0}".format(nCells)
print "number of cell sets: {0}".format(nCSets)
print "number of vertex sets: {0}".format(nVSets)

print "removing first 4 vertices from mesh"
X,Y,Z=e.get_coords()
eout.put_coord_names(["x","y","z"])
eout.put_coords(X[4:],Y[4:],Z[4:])


for set in range(nCSets):
    setID = e.get_elem_blk_ids()[set]
    setName = e.get_elem_blk_name(setID)
    cellType,numCells,numVertexPerCell,numAttr = e.elem_blk_info(setID)
    print "Cell set {0} {1}".format(setID,setName)
    print "\tNumber of cells: {0}".format(numCells)
    print "\tCell type: {0}".format(cellType)
    print "\tNumber of vertex per cell: {0}".format(numVertexPerCell)

    eout.put_elem_blk_info(setID,cellType,numCells,numVertexPerCell,1)
    eout.put_elem_blk_name(setID,setID)

    connect = e.get_elem_connectivity(setID)
    newConnect = np.array(connect[0],dtype=exo.c_int)-4
    eout.put_elem_connectivity(setID,newConnect)

eout.close()

