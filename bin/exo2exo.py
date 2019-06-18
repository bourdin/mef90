#!/usr/bin/env python
import exodus as exo
import sys
import os.path
import numpy as np
import pymef90


def exo2exo(fin,fout):
    import warnings
    warnings.filterwarnings('ignore','.*buffer.*',)

    cell1D = ("BAR","BAR2","BEAM2","BAR3","BEAM3")
    cell2D = ("TRI","TRI3","TRIANGLE","TRISHELL","TRISHELL3","TRI6","TRISHELL6","QUAD","QUAD4","SHELL","SHELL4","QUAD9","SHELL9")
    cell3D = ("TETRA","TETRA4","TETRA10","HEX","HEX8","HEX27")

    print "Opening {0}".format(fin)
    e = exo.exodus(fin,"r")
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
            print "\tvertex {0} is missing.".format(v)
            vertexReordering[v:] = vertexReordering[v:]-1
            vertexReordering[v]  = -1
            numMissing += 1
    if numMissing > 0:
        print ("{0} vertices are not referenced in the input file. Was the mesh renumbered?".format(numMissing))

    ### Create fixed file
    title     = "{0} fixed by {1}".format(e.title(),sys.argv[0])[:exo.MAX_LINE_LENGTH]
    nVertices -= numMissing
    eout = exo.exodus(fout,'w',numDims=dim,title=title,numNodes=nVertices,numElems=nCells,numBlocks=nCSets,numNodeSets=nVSets,numSideSets=0)
    ### Write coordinates
    ###
    X,Y,Z=e.get_coords()
    eout.put_coord_names(["x","y","z"][0:dim])
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
    maxID = 0
    usedID = [0,]
    for set in blocksOrder:
        setID = e.get_elem_blk_ids()[set]
        setName = e.get_elem_blk_name(setID)
        try:
            setFixedID = int(setName)
            if setFixedID in usedID:
                setFixedID =  max(usedID)+1    
            usedID.append(setFixedID)
        except ValueError:
            setFixedID =  max(usedID)+1
            usedID.append(setFixedID)
        cellType,numCells,numVertexPerCell,numAttr = e.elem_blk_info(setID)
        print('Assigning ID {0:4d} to cell set "{1}". \tmef90/vDef name will be cs{0:04}'.format(setFixedID,setName))
        print "\tNumber of cells: {0}".format(numCells)
        print "\tCell type: {0}".format(cellType)
        print "\tNumber of vertex per cell: {0}".format(numVertexPerCell)

        eout.put_elem_blk_info(setFixedID,cellType,numCells,numVertexPerCell,1)
        eout.put_elem_blk_name(setFixedID,setName)

        connect = np.array(e.get_elem_connectivity(setID)[0],dtype=exo.c_int)
        for i in range(len(connect)):
            connect[i] = 1+vertexReordering[connect[i]-1]
        eout.put_elem_connectivity(setFixedID,connect)

    ### Write node sets
    ###
    setIDs = e.get_node_set_ids()
    usedID = [0,]
    for set in setIDs:
        setName = e.get_node_set_name(set)
        try:
            setFixedID = int(setName)
            if setFixedID in usedID:
                setFixedID =  max(usedID)+1    
            usedID.append(setFixedID)
        except ValueError:
            setFixedID =  max(usedID)+1
            usedID.append(setFixedID)

        print('Assigning ID {0:4d} to vertex set "{1}". \tmef90/vDef name will be vs{0:04}'.format(setFixedID,setName))
        num_ns_nodes, num_ns_dist_facts = e.get_node_set_params(set)

        if num_ns_nodes == 0:
            print("\tset is empty and will NOT be written to output file")
        else:
            print("\tnumber of nodes: {1}".format(set,num_ns_nodes,setName))
            eout.put_node_set_params(setFixedID,num_ns_nodes, num_ns_dist_facts)
            eout.put_node_set(setFixedID,e.get_node_set_nodes(set))
            eout.put_node_set_name(setFixedID,setName)
    eout.close()
    e.close()

if __name__ == '__main__':
    if not (len(sys.argv) == 3):
        print "usage: {0} <input filename> <output filename>".format(sys.argv[0])
        sys.exit(-1)
    fin = sys.argv[-2]
    fout = sys.argv[-1]
    if not os.path.isfile(fin):
        print "usage: {0} <input filename> <output filename>".format(sys.argv[0])
        sys.exit(-1)
    if os.path.exists(fout):
        if pymef90.confirm("ExodusII file {0} already exists. Overwrite?".format(fout)):
            os.remove(fout)
        else:
            print "bye!"
            sys.exit(-1)

    exo2exo(fin,fout)
