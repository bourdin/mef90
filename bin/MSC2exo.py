#!/usr/bin/env python
import numpy as np
import exodus as exo
import argparse
import io
import pymef90

def MSCImporter(filename):
    TetFaces = [[1, 2, 4],
                [2, 3, 4],
                [1, 3, 2],
                [1, 4, 3]]
    f = io.open(filename, 'r') 
    nCell = 0
    nVert = 0
    cellSet = {}
    vertexSet = {}
    cellSetOffset = 0
    faceSetOffset = 10
    vertexSetOffset = 100
    cellOffset = 0
    line=f.readline().strip()
    while line!="":
        if line.startswith('sizing'):
            nCell = int(line.split()[2])
            nVert = int(line.split()[3])
            print ('Number of cells: {0}, number of vertices: {1}'.format(nCell,nVert))
        elif line.startswith('connectivity'):
            f.readline() # skip one line
            connect = []
            for cell in range(nCell):
                connect.append([int(x) for x in f.readline().strip().split()[2:]])
        elif line.startswith('coordinates'):
            f.readline() # skip one line
            coord = np.empty([nVert,3])
            for vert in range(nVert):
                line = f.readline().strip().replace('+', 'e+').replace('-','e-')
                coord[vert,:] = [float(x) for x in line.split()[1:]]
        elif line.startswith('define'):
            if line.strip().split()[1].startswith('element'):
                cellSetOffset += 1
                print 'Cell set {0} is now cs{1:04d}'.format(line.split()[-1],cellSetOffset)
                cellSet[cellSetOffset] = {}
                cellSet[cellSetOffset]['numVPE'] = 4
                cellSet[cellSetOffset]['elemType'] = 'TETRA4'
                line=f.readline()
                csStart = int(line.strip().split()[0])
                csEnd   = int(line.strip().split()[2])
                cellSet[cellSetOffset]['connect'] = []
                for cell in range(csStart,csEnd+1):
                    for v in connect[cell-1]:
                        cellSet[cellSetOffset]['connect'].append(v)
                cellOffset += csEnd - csStart + 1
                #line=f.readline()
            elif line.strip().split()[1].startswith('facemt'):
                faceSetOffset += 1
                print 'Cell set {0} is now cs{1:04d}'.format(line.split()[-1],faceSetOffset)
                cellSet[faceSetOffset] = {}
                cellSet[faceSetOffset]['numVPE'] = 3
                cellSet[faceSetOffset]['elemType'] = 'TRI3'
                cellSet[faceSetOffset]['connect'] = []
                line=f.readline()
                while line.strip()[0] in '0123456789':
                    for ff in line.strip().split():
                        cell = int(ff.split(':')[0])
                        face = int(ff.split(':')[1])
                        for v in TetFaces[face-1]:
                            cellSet[faceSetOffset]['connect'].append(connect[cell-1][v-1])
                    line=f.readline()
            if line.strip().split()[1].startswith('ndsq'):
                print 'Vertex set {0} is now vs{1:04d}'.format(line.split()[-1],vertexSetOffset)
                vertexSet[vertexSetOffset] = []

                line=f.readline()
                while line.strip()[0] in '0123456789':
                    for v in line.strip().split():
                        vertexSet[vertexSetOffset].append(int(v))
                    line=f.readline()
                vertexSetOffset += 1
        line=f.readline().strip()
    return coord,vertexSet,cellSet

def FixOrientation(connect,coord):
    e0 = (coord[1][0] - coord[0][0],coord[1][1] - coord[0][1])
    e1 = (coord[2][0] - coord[0][0],coord[2][1] - coord[0][1])
    det = e0[0] * e1[1] - e0[1]*e1[0]
    if det <= 0.:
        return (connect[0],connect[2],connect[1])
    else:
        return connect

#------Function for writing to exo format
def exoWriter(coords,vertexSets,cellSets,exoFile):
    cell1D = ("BAR","BAR2","BEAM2","BAR3","BEAM3")
    cell2D = ("TRI","TRI3","TRIANGLE","TRISHELL","TRISHELL3","TRI6","TRISHELL6","QUAD","QUAD4","SHELL","SHELL4","QUAD9","SHELL9")
    cell3D = ("TETRA","TETRA4","TETRA10","HEX","HEX8","HEX27")

    X = coords[:,0]         #set of all X coords
    Y = coords[:,1]         #set of all Y coords
    Z = coords[:,2]         #set of all Z coords
    if max(Z) == min(Z):        #two dimensional case if Z is static
        numDim = 2
    else:                       #otherwise we assume we have a three dimensional case
        numDim = 3

    numElem = 0
    for k in cellSets.keys():       #finding number of elements
        numElem += len(cellSets[k]['connect'])/cellSets[k]['numVPE']

    filename = exoFile           #this is the file specified when function is called

    #setting up exo file for writing
    e=exo.exodus(filename, mode='w',title='title',numDims=numDim,numNodes=len(X), 
                  numElems=numElem,numBlocks=len(cellSets),numNodeSets=len(vertexSets),numSideSets=0)
    #coordinates
    if numDim == 3:
        e.put_coord_names(["x","y","z"])    #name of each coordinate
    else:
        e.put_coord_names(["x","y"])    #name of each coordinate
    e.put_coords(X,Y,Z)                 #actual coordinates

    #block info and connectivity
    ###
    ### Sort cell sets in order to write boundary cell sets last
    ###
    blocksOrder = []
    if numDim == 2:
        cellCoDim0 = cell2D
        cellCoDim1 = cell1D
    else:
        cellCoDim0 = cell3D
        cellCoDim1 = cell2D

    for setID in cellSets.keys():
        if cellSets[setID]['elemType'].upper() in cellCoDim0:
            blocksOrder.append(setID)
    for setID in cellSets.keys():
        if cellSets[setID]['elemType'].upper() in cellCoDim1:
            blocksOrder.append(setID)

    #for setID in cellSets.keys():
    for setID in blocksOrder:
            ###---setID, elemType, num elems, num nodes per elem, num attributes per elem
            e.put_elem_blk_info(setID,cellSets[setID]['elemType'],len(cellSets[setID]['connect'])/cellSets[setID]['numVPE'],cellSets[setID]['numVPE'],0)
            ###---setID, connectivity table
            e.put_elem_connectivity(setID,cellSets[setID]['connect'])
    #node set info
    for setID in vertexSets.keys():
        ###---setID, num nodes, num distribution factors in a node set
        e.put_node_set_params(setID,len(vertexSets[setID]),0)
        ###---setID, nodes
        e.put_node_set(setID,vertexSets[setID])
    e.close()



#------Function for reordering connectivity table
def reorder(celltype,connect):
    ### other celltypes look okay, but we can add more here
    orderFix = {
               # 10-node tetrahedron
               11:[1,4,2,3,8,10,5,7,9,6],   
               #27 node hex
               12:[1,5,8,4,2,6,7,3,11,18,16,10,9,17,20,14,13,19,15,12,27,23,24,21,26,22,25]
               }
    #computes reordering
    if celltype in orderFix:
        reordered = [connect[i-1] for i in orderFix[celltype]]
    else:
        reordered = connect
    return reordered                    #returns fixed connect table


#------Main Function
def main():
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument("MSCFile", help = "The name of the mesh file to be parsed.", type = str)
    parser.add_argument("exoFile", help = "The name of the exodus file to be written.", type = str)
    parser.add_argument("--force",action="store_true",default=False,help="Overwrite existing files without prompting")
    args = parser.parse_args()

    if  os.path.exists(args.exoFile):
        if args.force:
            os.remove(args.exoFile)
        else:
            if pymef90.confirm("ExodusII file {0} already exists. Overwrite?".format(args.exoFile)):
                os.remove(args.exoFile)
            else:
                print '\n\t{0} was NOT generated from {1}\n'.format(args.exoFile,args.MSCFile)
                return -1
    try:
        (coord,vertexSet,cellSet) = MSCImporter(args.MSCFile)
    except TypeError:
        print "Cannot read {0}".format(args.MSCFile)
        return -1
    exoWriter(coord,vertexSet,cellSet,args.exoFile)

if __name__ == '__main__':
    main()

