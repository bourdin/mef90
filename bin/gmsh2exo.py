#!/usr/bin/env python
import numpy as np
import exodus as exo
import argparse
import io
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

    # Opening and reading mesh file
    f = io.open(filename, 'r') 
    if not f.readline() == '$MeshFormat\n':     #checking if mesh file
        print "Unknown file format"
        return -1
    buffer = f.readline()
    if not buffer.lstrip().startswith('2.2'):
        print "gmsh2exo.py can only read msh2 format, but got format {0}.\nRegenerate your mesh with the option -format msh2".format(buffer.rstrip())
        return -1
    if not f.readline() == '$EndMeshFormat\n':
        print "Something weird is happening here at EndMeshFormat"
        return -1

    #checking to see if we are at the correct location in mesh file again
    if not f.readline() == '$Nodes\n':
        print "Something weird is happening here at Nodes"
        return -1
    nVert = int(f.readline())                   #number of vertices

    #creating and filling list of coordinates from mesh file
    coord = []
    for v in range(nVert):
        coord.append([float(x) for x in f.readline().strip().split()[1:]])
    coord = np.array(coord)         #making coord list into numpy array
    
    #checking to see if we are at the correct location in mesh file
    if not f.readline() == '$EndNodes\n':
        print "Something weird is happening here at EndNodes"
        return -1

    #checking to see if we are at the correct location in mesh file again
    if not f.readline() == '$Elements\n':
        print "Something weird is happening here at Elements"
        return -1

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
        print "Something weird is happening here at End Elements"
        return -1
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

def confirm(prompt=None, resp=False):
    """prompts for yes or no response from the user. Returns True for yes and
    False for no.

    'resp' should be set to the default value assumed by the caller when
    user simply types ENTER.

    >>> confirm(prompt='Create Directory?', resp=True)
    Create Directory? [y]|n: 
    True
    >>> confirm(prompt='Create Directory?', resp=False)
    Create Directory? [n]|y: 
    False
    >>> confirm(prompt='Create Directory?', resp=False)
    Create Directory? [n]|y: y
    True

    """
    
    if prompt is None:
        prompt = 'Confirm'

    if resp:
        prompt = '%s [%s]|%s: ' % (prompt, 'y', 'n')
    else:
        prompt = '%s [%s]|%s: ' % (prompt, 'n', 'y')
        
    while True:
        ans = raw_input(prompt)
        if not ans:
            return resp
        if ans not in ['y', 'Y', 'n', 'N']:
            print 'please enter y or n.'
            continue
        if ans == 'y' or ans == 'Y':
            return True
        if ans == 'n' or ans == 'N':
            return False


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
            if confirm("ExodusII file {0} already exists. Overwrite?".format(args.exoFile)):
                os.remove(args.exoFile)
            else:
                print '\n\t{0} was NOT generated from {1}\n'.format(args.exoFile,args.gmeshFile)
                return -1
    try:
        (coord,vertexSet,cellSet) = GMSHImporter(args.gmeshFile)
    except TypeError:
        print "Cannot read {0} is it in gmsh 2 format?".format(args.gmeshFile)
        return -1
    exoWriter(coord,vertexSet,cellSet,args.exoFile)

if __name__ == '__main__':
    main()

