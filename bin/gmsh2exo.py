#!/usr/bin/env python
import numpy as np
import exodus as exo
import argparse
import io
import pymef90

#--------Function for parsing mesh file
def GMSHImporter(filename):
    elemDim = {1:2,   # 2 node line
               2:3,   # 3 node triangle
               3:4,   # 4 node quad
               4:4,   # 4-node tetrahedron
               5:8,   # 8 node hexahedron
               8:3,   # 3 node line
               9:6,   # 6 node triangle
               10:9,  # 9 node quad
               11:10, # 10 node tetrahedron
               12:27, # 27 node hexahedron
               15:1   # vertex
              }
    elemType = {1:"BAR2",       #2 node line
                2:"TRI3",       #3 node triangle
                3:"QUAD4",      #4 node quad
                4:"TETRA4",     #4 node tet
                5:"HEX8",       #8 node hexahedron
                8:"BAR3",       #3 node line
                9:"TRI6",       #6 node triangle
                10:"QUAD9",     #9 node quad
                11:"TETRA10",   #10 node tet
                12:"HEX27",     #27 node hexahedron
                15:""           #vertex
                }

    cl = 0
    # Opening and reading mesh file
    f = io.open(filename, 'r') 
    if not f.readline() == '$MeshFormat\n':     #checking if mesh file
        raise(RuntimeError,"Unable to parse file {0}".format(filename))
    buffer = f.readline()
    if not buffer.lstrip().startswith('2.'):
        raise(RuntimeError,"gmsh2exo.py can only read msh2 format, but got format {0}.\nRegenerate your mesh with the option -format msh2".format(buffer.rstrip()))
    if not f.readline() == '$EndMeshFormat\n':
        raise(RuntimeError,"Unable to parse file {0}: $EndMeshFormat not found when expected".format(filename))

    #checking to see if we are at the correct location in mesh file again
    if not f.readline() == '$Nodes\n':
        raise(RuntimeError,"Unable to parse file {0}:  $Nodes not found when expected".format(filename))
    nVert = int(f.readline())                   #number of vertices

    #creating and filling list of coordinates from mesh file
    coord = []
    for v in range(nVert):
        coord.append([float(x) for x in f.readline().strip().split()[1:]])
    coord = np.array(coord)         #making coord list into numpy array
    
    #checking to see if we are at the correct location in mesh file
    if not f.readline() == '$EndNodes\n':
        raise(RuntimeError,"Unable to parse file {0}:  $EndNodes not found when expected".format(filename))

    #checking to see if we are at the correct location in mesh file again
    if not f.readline() == '$Elements\n':
        raise(RuntimeError,"Unable to parse file {0}:  $Elements not found when expected".format(filename))

    #finding number of cells in mesh file
    nCell = int(f.readline())
    
    #creating and filling vertexSet and cellSet dictionaries
    vertexSet = {}
    cellSet   = {}
    for c in range(nCell):
        tmp = [int(x) for x in f.readline().strip().split()[1:]]
        cellType = tmp[0]       #first number is cell type
        if cellType == 15:      #this is when element is a node
            ### vertex sets:
            tag = tmp[2]        #third number in tmp is a tag
            if not tag in vertexSet.keys():         #adding new tags
                vertexSet[tag] = []
            vertexSet[tag].append(tmp[-1])          #setting element to tag
        else:
            ### cell sets:
            tag = tmp[2]                    #third number in tmp is a tag
            if not tag in cellSet.keys():   #looking for tag in cellSet  
                cellSet[tag] = {}           #creating new tag if doesn't exist
                cellSet[tag]['connect'] = []
            cellConnect = tmp[-elemDim[cellType]:]
            cellConnect = reorder(cellType,cellConnect)
            if elemType[cellType] == 'TRI3':
                locCoord = (coord[cellConnect[0]-1],coord[cellConnect[1]-1],coord[cellConnect[2]-1])
                cellConnect = FixOrientation(cellConnect,locCoord)
            cellSet[tag]['connect'] += cellConnect  #reordering and adding to array        #adding connectivity table to tag
            cellSet[tag]['numVPE']  = elemDim[cellType]     #adding number of vertex per element to tag
            cellSet[tag]['elemType'] = elemType[cellType]   #adding element type to tag
                        
    #checking to make sure we are at end of file
    if not f.readline() == '$EndElements\n':    
        raise(RuntimeError,"Unable to parse file {0}:  $EndElements not found when expected".format(filename))

    #block info and connectivity
    ###
    ### Sort cell sets in order to write boundary cell sets last
    ###
    X = coord[:,0]
    Y = coord[:,1]
    Z = coord[:,2]
    if max(Z) == min(Z):
        numDim = 2
    else:
        numDim = 3

    return coord,vertexSet,cellSet,numDim

def FixOrientation(connect,coord):
    e0 = (coord[1][0] - coord[0][0],coord[1][1] - coord[0][1])
    e1 = (coord[2][0] - coord[0][0],coord[2][1] - coord[0][1])
    det = e0[0] * e1[1] - e0[1]*e1[0]
    if det <= 0.:
        return (connect[0],connect[2],connect[1])
    else:
        return connect

#------Function for writing to exo format
def exoWriter(coords,vertexSets,cellSets,numDim,exoFile):
    cell1D = ("BAR","BAR2","BEAM2","BAR3","BEAM3")
    cell2D = ("TRI","TRI3","TRIANGLE","TRISHELL","TRISHELL3","TRI6","TRISHELL6","QUAD","QUAD4","SHELL","SHELL4","QUAD9","SHELL9")
    cell3D = ("TETRA","TETRA4","TETRA10","HEX","HEX8","HEX27")

    # Reorder blocks with cells of coDimension 0 first
    blocksOrder = []
    if numDim == 2:
        cellCoDim0 = cell2D
        cellCoDim1 = cell1D
        cellCoDim2 = ()
    else:
        cellCoDim0 = cell3D
        cellCoDim1 = cell2D
        cellCoDim2 = cell1D

    #in 3D, cellCoDim2 cells need to be converted into a vertex set
    for setID in cellSets.keys():
        if cellSets[setID]['elemType'].upper() in cellCoDim2:
            vs = []
            for v in cellSets[setID]['connect']:
                if v not in vs:
                    vs.append(v)
            cellSets.pop(setID,'None')
            if setID in vertexSets.keys():
                print("Codimension 2 cell set {0} renamed vertex set {1} so as not to clash with existing vertex set".format(setID,max(vertexSets.keys())+1))
                setID = max(vertexSets.keys())+1
            vertexSets[setID]=vs


    for setID in cellSets.keys():
        if cellSets[setID]['elemType'].upper() in cellCoDim0:
            blocksOrder.append(setID)
    for setID in cellSets.keys():
        if cellSets[setID]['elemType'].upper() in cellCoDim1:
            blocksOrder.append(setID)

    #Writting exodusII file
    numElem = 0
    for k in cellSets.keys():       #finding number of elements
        numElem += len(cellSets[k]['connect'])/cellSets[k]['numVPE']
    numNodes = len(coords[:,0])
    e=exo.exodus(exoFile, mode='w',title='title',numDims=numDim,numNodes=numNodes, 
                  numElems=numElem,numBlocks=len(cellSets),numNodeSets=len(vertexSets),numSideSets=0)
    #coordinates
    if numDim == 3:
        e.put_coord_names(["x","y","z"])
    else:
        e.put_coord_names(["x","y"])
    e.put_coords(coords[:,0],coords[:,1],coords[:,2])

    for setID in blocksOrder:
            ###---setID, elemType, num elems, num nodes per elem, num attributes per elem
            e.put_elem_blk_info(setID,cellSets[setID]['elemType'],len(cellSets[setID]['connect'])/cellSets[setID]['numVPE'],cellSets[setID]['numVPE'],0)
            ###---setID, connectivity table
            e.put_elem_connectivity(setID,cellSets[setID]['connect'])
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
    parser.add_argument("gmeshFile", help = "The name of the mesh file to be parsed.", type = str)
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
                print '\n\t{0} was NOT generated from {1}\n'.format(args.exoFile,args.gmeshFile)
                return -1
    try:
        (coord,vertexSet,cellSet,numDim) = GMSHImporter(args.gmeshFile)
    except TypeError:
        print "Cannot read {0} is it in gmsh 2 format?".format(args.gmeshFile)
        return -1
    exoWriter(coord,vertexSet,cellSet,numDim,args.exoFile)

if __name__ == '__main__':
    main()

