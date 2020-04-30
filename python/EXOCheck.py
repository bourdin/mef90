#!/usr/bin/env python
import sys
if sys.version_info.major == 3:
    import exodus3 as exo
else:
    import exodus2 as exo
import os.path
import numpy as np

f = sys.argv[-1]
if not os.path.isfile(f) or not (len(sys.argv) == 2):
    print "usage: {0} <filename>".format(sys.argv[0])
    exit

e = exo.exodus(f,"r")
dim       = e.num_dimensions()
nVertices = e.num_nodes()
nCells    = e.num_elems()
nCSets    = e.num_blks()
nVSets    = e.num_node_sets()

print "number of dimensions: {0}".format(dim)
print "number of vertices: {0}".format(nVertices)
print "number of cells: {0}".format(nCells)
print "number of cell sets: {0}".format(nCSets)
print "number of vertex sets: {0}".format(nVSets)

#listedVertices = np.empty([nVertices,],dtype=bool)
listedVertices = [False,]*nVertices
maxV = 0
minV = 2*nVertices
for set in range(nCSets):
	setID = e.get_elem_blk_ids()[set]
	print "Cell set {0}".format(setID)
	connect = e.get_elem_connectivity(setID)
	print "\tNumber of cells: {0}".format(connect[1])
	for v in connect[0]:
		if v > nVertices:
			print "PROBLEM: vertex {0} is out of range {1}".format(v,nVertices)
			exit
		else:
			listedVertices[v-1] = True
			minV = min(minV,v)
			maxV = max(maxV,v)

print "vertex range: ({0}/{1})".format(minV,maxV)
for v in range(nVertices):
	if not listedVertices[v]:
		print "PROBLEM: vertex {0} is not part of any cell".format(v+1)

for set in range(nVSets):
	setID = e.get_elem_blk_ids()[nVSets]
	print "Vertex set {0}".format(setID)