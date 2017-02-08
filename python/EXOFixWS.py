#!/usr/bin/env python

import exodus as exo
import sys
import os.path
import numpy as np

cell1D = ("BAR","BAR2","BEAM2","BAR3","BEAM3")
cell2D = ("TRI","TRI3","TRIANGLE","TRISHELL","TRISHELL3","TRI6","TRISHELL6","QUAD","QUAD4","SHELL","SHELL4","QUAD9","SHELL9")
cell3D = ("TETRA","TETRA4","TETRA10","HEX","HEX8","HEX27")

f = sys.argv[-1]
if not os.path.isfile(f) or not (len(sys.argv) == 2):
    print "usage: {0} <filename>".format(sys.argv[0])
    exit
else:
    fout = f.split(".gen")[0]+"_fixed.gen"
    if os.path.exists(fout):
        print "{0} exists. Delete and run again".format(fout)
        exit

print "Opening {0}".format(f)
e = exo.exodus(f,"r")
dim       = e.num_dimensions()
nVertices = e.num_nodes()
nCells    = e.num_elems()
nCSets    = e.num_blks()
nVSets    = e.num_node_sets()

print "\tnumber of dimensions: {0}".format(dim)
print "\tnumber of vertices: {0}".format(nVertices)
print "\tnumber of cells: {0}".format(nCells)
print "\tnumber of cell sets: {0}".format(nCSets)
print "\tnumber of vertex sets: {0}".format(nVSets)



### Reorder cell sets 
### List cells of coDIm0 first then cells of coDim 1
###
print "Reordering cell sets by increasing co-dimension"
cellsType = []
for set in range(nCSets):
    setID = e.get_elem_blk_ids()[set]
    cellType,numCells,numVertexPerCell,numAttr = e.elem_blk_info(setID)
    cellsType.append(cellType)

blocksOrder = []
if dim == 2:
    cellCoDim0 = cell2D
    cellCoDim1 = cell1D
else:
    cellCoDim0 = cell3D
    cellCoDim1 = cell2D

for i in range(len(cellsType)):
    if cellsType[i].upper() in cellCoDim0:
        blocksOrder.append(i)
for i in range(len(cellsType)):
    if cellsType[i].upper() in cellCoDim1:
        blocksOrder.append(i)


### Find hanging nodes
###
### 1. find largest vertex index in all connectivity tables
print "Removing hanging nodes"
minV = 100000000
maxV = 0
for set in range(nCSets):
    setID = e.get_elem_blk_ids()[set]
    connect = e.get_elem_connectivity(setID)
    for v in connect[0]:
        maxV = max(maxV,v)
        minV = min(minV,v)
print "\tVertex range: {0}/{1}".format(minV,maxV)

listedVertices = [False,]*maxV
for set in range(nCSets):
    setID = e.get_elem_blk_ids()[set]
    connect = e.get_elem_connectivity(setID)
    for v in connect[0]:
        listedVertices[v-1] = True

numMissing = 0 
vertexReordering = np.array(range(maxV),dtype=int)
for v in range(len(listedVertices)):
    if not listedVertices[v]:
        print "\tvertex {0} is missing".format(v)
        vertexReordering[v:] = vertexReordering[v:]-1
        vertexReordering[v]  = -1
        numMissing += 1

### Create fixed file
title     = "{0} fixed by {1}".format(e.title(),sys.argv[0])[:exo.MAX_LINE_LENGTH]
nVertices -= numMissing
eout = exo.exodus(fout,'w',numDims=dim,title=title,numNodes=nVertices,numElems=nCells,numBlocks=nCSets,numNodeSets=nVSets,numSideSets=0)
### Write coordinates
###
X,Y,Z=e.get_coords()
eout.put_coord_names(["x","y","z"])
fixedX = np.empty(nVertices,dtype=exo.c_double)
fixedY = np.empty(nVertices,dtype=exo.c_double)
fixedZ = np.empty(nVertices,dtype=exo.c_double)

for i in range(nVertices):
    if listedVertices[i]:
        fixedX[i] = X[vertexReordering[i]]
        fixedY[i] = Y[vertexReordering[i]]
        fixedZ[i] = Z[vertexReordering[i]]
eout.put_coords(fixedX,fixedY,fixedZ)

### Write cell sets
###
for set in blocksOrder:
    setID = e.get_elem_blk_ids()[set]
    setName = e.get_elem_blk_name(setID)
    cellType,numCells,numVertexPerCell,numAttr = e.elem_blk_info(setID)
    print "Cell set {0} {1} {2}".format(set,setID,setName)
    print "\tNumber of cells: {0}".format(numCells)
    print "\tCell type: {0}".format(cellType)
    print "\tNumber of vertex per cell: {0}".format(numVertexPerCell)

    eout.put_elem_blk_info(setID,cellType,numCells,numVertexPerCell,1)
    eout.put_elem_blk_name(setID,setID)

    connect = np.array(e.get_elem_connectivity(setID)[0],dtype=exo.c_int)
    for i in range(len(connect)):
        connect[i] = 1+vertexReordering[connect[i]-1]
    eout.put_elem_connectivity(setID,connect)

print "Ignoring nodesets for now"
eout.close()
