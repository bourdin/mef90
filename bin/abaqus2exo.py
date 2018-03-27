#!/usr/bin/env python
import numpy as np
import exodus as exo
import argparse
import time

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



#--------Function for parsing mesh file
def ABAQUSImporter(filename):
    import sys
    elemDim = {
               "S3R":3,      # 3 node triangle
               "S4R":4,      # 4 node quad
               "C3D8R":8,    # 8 node hex
               "C3D10":10,   # 10 node tet
               "C3D4":4,     # 4 node tet
               "STRI3":3,    # 3 node triangle
               "B21":2,      # 2 node line
               "T2D2":2,     # 2 node line
               "CPE3":3,     # 3 node triangle
               "CPE6":6,     # 6 node triangle
               "CPE6S":6,    # 6 node triangle
               "CPE4":4,     # 4 node quad
               "CPE4S":4,    # 4 node quad
               "CPE4R":4,    # 4 node quad
               "CPE8":8,     # 8 node quad
               "CPE8S":8,    # 8 node quad
               "CPS3":3,     # 3 node tri
               "CPS6":6,     # 6 node tri
               "CPS4":4,     # 4 node quad
               "CPS4S":4,    # 4 node quad
               "CPS4R":4,    # 4 node quad
               "C3D4":4,     # 4 node tet
               "C3D10":10,   # 10 node tet
               "C3D10S":10,  # 10 node tet
               "C3D8":8,     # 8 node hex
               "C3D20":20,   # 20 node hex
               "C3D20S":20,  # 20 node hex
              }
    elemType = {
                "S3R":"TRI3",       # 3 node triangle
                "S4R":"QUAD4",      # 4 node quad
                "C3D8R":"HEX8",     # 8 node hex
                "C3D10":"TETRA10",  # 10 node tet
                "C3D4":"TETRA4",    # 4 node tet
                "STRI3":"TRI3",     # 3 node triangle
                "B21":"BAR2",       # 2 node line
                "T2D2":"BAR2",      # 2 node line
                "CPE3":"TRI3",      # 3 node triangle
                "CPE6":"TRI6",      # 6 node triangle
                "CPE6S":"TRI6",     # 6 node triangle
                "CPE4":"QUAD4",     # 4 node quad
                "CPE4S":"QUAD4",    # 4 node quad
                "CPE4R":"QUAD4",    # 4 node quad
                "CPS3":"TRI3",      # 3 node triangle
                "CPS6":"TRI6",      # 6 node triangle
                "CPS4":"QUAD4",     # 4 node quad
                "CPS4S":"QUAD4",    # 4 node quad
                "CPS4R":"QUAD4",    # 4 node quad
                "CPE8":"QUAD8",     # 8 node quad
                "CPE8S":"QUAD8",    # 8 node quad
                "C3D4":"TETRA4",    # 4 node tet
                "C3D10":"TETRA10",  # 10 node tet
                "C3D10S":"TETRA10", # 10 node tet
                "C3D8":"HEX8",      # 8 node hex
                "C3D20":"HEX20",    # 20 node hex
                "C3D20S":"HEX20",   # 20 node hex
                }

    # Opening and reading abaqus file
    f = open(filename, 'r') 
    line = ' '
    cellSet = {}
    vertexSet = {}
    nodeNames = {}
    while line != '':
        if line.startswith('*NODE'):
            sys.stdout.write( "Parsed {0} \n".format(line.strip()))
            (coord,nodeID,line,f,order) = readCoords(f,line)
        elif line.startswith('*ELEM'):
            sys.stdout.write( "Parsed {0} ".format(line.strip()))
            #sys.stdout.write( "elemType: {0}".format(elemType))
            (cellSet,elemNames,line,f) = readElems(f,line,nodeID,elemDim,elemType,cellSet,order)
            #print cellSet[cellName]['elemType'],'\n'
        elif line.startswith('*NSET'):
            sys.stdout.write( "Parsed {0} \n".format(line.strip()))
            (vertexSet,nodeNames,line,f) = readNodes(f,line,nodeID,vertexSet,nodeNames,order)
        else:
            line = f.readline()
    return coord,vertexSet,cellSet,elemNames,nodeNames
    
def readCoords(f,line):
    #creating and filling list of coordinates from abaqus file
    coord = []
    nodeID = []
    line = f.readline()
    #skip lines describing nodeset, and coords list ends on '**'
    while not line.startswith('*'):
        coord.append([float(x) for x in line.strip().split(", ")[1:]])
        nodeID.append(int(line.strip().split(", ")[0]))
        line = f.readline()
    order = [-1]*max(nodeID)
    coord = np.array(coord)         #making coord list into numpy array
    for i in range(len(nodeID)):
        order[nodeID[i]-1] = i+1
    return coord,nodeID,line,f,order 
    
def readElems(f,line,nodeID,elemDim,elemType,cellSet,order):
    #creating and filling cellSet dictionary
    #get rid of word '*ELEMENT'
    elemTmp = [x for x in line.strip().split(",")[1:]]
    #sometimes type comes first in line, sometimes name comes first
    if elemTmp[0].startswith('TYPE') or elemTmp[0].startswith(' TYPE'):
        cellType = elemTmp[0].split('=')[1]
        cellName = elemTmp[1].split('=')[1]
    else:
        cellName = elemTmp[0].split('=')[1]
        cellType = elemTmp[1].split('=')[1] 
    line = f.readline()     #continue to actual data        
    while not line.startswith('*') and not line == '':
        #get rid of index, this is connectivity table
        #this condition is for reading lines which have a comma at the end
        #splitting at "," will return an additional string with nothing in it
        if line.strip().split(",")[len(line.strip().split(",")[0:])-1] == '':
            cellConnect = [order[int(x)-1] for x in line.strip().split(",")[1:len(line.strip().split(","))-1]]
            line = f.readline()
            if line.strip().split(",")[len(line.strip().split(",")[0:])-1] == '':
                tmp = [order[int(x)-1] for x in line.strip().split(",")[0:len(line.strip().split(","))-1]]
                for i in range(len(tmp)):
                    cellConnect.append(tmp[i])
                line = f.readline()
                tmp = [order[int(x)-1] for x in line.strip().split(",")[0:]]
                for i in range(len(tmp)):
                    cellConnect.append(tmp[i])
            else:
                tmp = [order[int(x)-1] for x in line.strip().split(",")[0:]]
                for i in range(len(tmp)):
                    cellConnect.append(tmp[i])
        else:
            cellConnect = [order[int(x)-1] for x in line.strip().split(",")[1:]]    
        ### cell sets:
        if not cellName in cellSet.keys():   #looking for tag in cellSet  
            cellSet[cellName] = {}           #creating new tag if doesn't exist
            cellSet[cellName]['connect'] = []
            cellSet[cellName]['numVPE']  = elemDim[cellType]     #adding number of vertex per element to array
            cellSet[cellName]['elemType'] = elemType[cellType]   #adding element type to array
            print ("set type {0} converted to {1}\n".format(cellType,cellSet[cellName]['elemType']))
        cellSet[cellName]['connect'] += cellConnect  #reordering and adding to array
        line = f.readline() 

    #create dictionary linking cellSet number to its name
    cellNames = [x for x in cellSet]
    elemNames = {}
    for i in range(len(cellNames)):
        elemNames[i] = cellNames[i] 
    return cellSet,elemNames,line,f

def readNodes(f,line,nodeID,vertexSet,nodeNames,order):
    #now filling the node set information 
    #tmp is the node set name
    tmp = line.strip().split(",")[1]
    tmp = tmp.split('=')[1]
    #continue to next line which has actual node numbers
    line = f.readline()
    #add nodes to list
    #sometimes the node set will end with '**', and the next one will begin right after
    while not line.startswith('*'):
        #this condition is for reading lines which have a comma at the end
        #splitting at "," will return an additional string with nothing in it
        if line.strip().split(",")[len(line.strip().split(",")[0:])-1] == '':
            nodes = [order[int(x)-1] for x in line.strip().split(",")[0:len(line.strip().split(","))-1]]
        else:
            nodes = [order[int(x)-1] for x in line.strip().split(",")[0:]]
        #if the node set name isnt already in vertexSet
        if not tmp in vertexSet.keys():
            vertexSet[tmp] = []
        #add the actual node to vertexSet
        for i in range(len(nodes)):
            list_node = nodes[i]
            vertexSet[tmp].append(list_node)
        #move onto next line
        line = f.readline()
        #this condition is for the odd declaration of element sets that
        #sometimes show up between node sets
        if line.startswith('*ELSET') or line == '':
            while not line.startswith('*NSET') and not line.startswith('**'):
                line = f.readline()
                #check to see if we've reached the end of the file
                if line == '':
                    break
        #check to see if we've reached the end of the file
        if line == '':
            break
    #create dictionary linking the node set number to its name
    vertNames = [x for x in vertexSet]
    for i in range(len(vertNames)):
        nodeNames[i] = vertNames[i]
    return vertexSet,nodeNames,line,f 
    
    
#------Function for writing to exo format
def exoWriter(coords,vertexSets,cellSets,filename,elemNames,nodeNames):
    cell1D = ("BAR","BAR2","BEAM2","BAR3","BEAM3")
    cell2D = ("TRI","TRI3","TRIANGLE","TRISHELL","TRISHELL3","TRI6","TRISHELL6","QUAD","QUAD4","SHELL","SHELL4","QUAD9","SHELL9")
    cell3D = ("TETRA","TETRA4","TETRA10","HEX","HEX8","HEX27")
    X = coords[:,0]         #set of all X coords
    Y = coords[:,1]         #set of all Y coords
    if len(coords[1])>2 and max(coords[1]) == min(coords[1]):
       Z = coords[:,2]         #set of all Z coords
       numDim = 3
    else:                       #otherwise we assume we have a three dimensional case
       Z = [0]*len(X) 
       numDim = 2

    numElem = 0
    for k in cellSets.keys():       #finding number of elements
        numElem += len(cellSets[k]['connect'])/cellSets[k]['numVPE']

    #setting up exo file for writing
    e=exo.exodus(filename, mode='w',title='title',numDims=numDim,numNodes=len(X), 
                  numElems=numElem,numBlocks=len(cellSets),numNodeSets=len(vertexSets),numSideSets=0)
    e.put_coord_names(["x","y","z"][0:numDim])    #name of each coordinate
    e.put_coords(X,Y,Z)                 #actual coordinates
    
    print("\n")
    elemIDs = {}
    usedID = [0,]

    #block info and connectivity
    ###
    ### Sort cell sets in order to write boundary cell sets last
    ###
    blocksOrder = []
    c0blocks = []
    c1blocks = []
    if numDim == 2:
        cellCoDim0 = cell2D
        cellCoDim1 = cell1D
    else:
        cellCoDim0 = cell3D
        cellCoDim1 = cell2D

    for setID in cellSets.keys():
        if cellSets[setID]['elemType'].upper() in cellCoDim0:
            c0blocks.append(setID)
    c0blocks.sort()
    for setID in cellSets.keys():
        if cellSets[setID]['elemType'].upper() in cellCoDim1:
            c1blocks.append(setID)
    c1blocks.sort()
    blocksOrder = c0blocks + c1blocks
        
    for bn in blocksOrder:
        try:
            setFixedID = int(bn)
            if setFixedID in usedID:
                setFixedID =  max(usedID)+1    
            usedID.append(setFixedID)
        except ValueError:
            setFixedID =  max(usedID)+1
            usedID.append(setFixedID)
        elemIDs[bn] = setFixedID
        print("Assigning ID {0:4d} to ELSET {1}. \tmef90/vDef name will be cs{0:04}".format(elemIDs[bn],bn))

    #block info and connectivity
    for setID in blocksOrder:#cellSets.keys():
        ###---setID, elemType, num elems, num nodes per elem, num attributes per elem
        e.put_elem_blk_info(elemIDs[setID],cellSets[setID]['elemType'],len(cellSets[setID]['connect'])/cellSets[setID]['numVPE'],cellSets[setID]['numVPE'],0)
        ###---setID, connectivity table
        e.put_elem_connectivity(elemIDs[setID],cellSets[setID]['connect'])

    ###---assign string names of sets to index
    for i in range(len(elemNames)):
        e.put_elem_blk_name(elemIDs[elemNames[i]],elemNames[i])
    
    nodeIDs = {}
    usedID = [0,]
    for i in range(len(nodeNames)):
        try:
            setFixedID = int(nodeNames[i])
            if setFixedID in usedID:
                setFixedID =  max(usedID)+1    
            usedID.append(setFixedID)
        except ValueError:
            setFixedID =  max(usedID)+1
            usedID.append(setFixedID)
        nodeIDs[nodeNames[i]] = setFixedID
        print("Assigning ID {0:4d} to NSET {1}. \tmef90/vDef name will be vs{0:04}".format(nodeIDs[nodeNames[i]],nodeNames[i]))

    #node set info
    for setID in vertexSets.keys():
        ###---setID, num nodes, num distribution factors in a node set
        e.put_node_set_params(nodeIDs[setID],len(vertexSets[setID]),0)
        ###---setID, nodes
        e.put_node_set(nodeIDs[setID],vertexSets[setID])

    for i in range(len(nodeNames)):
        e.put_node_set_name(nodeIDs[nodeNames[i]],nodeNames[i])
    e.close()


#------Main Function
def Main():
    import sys
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument("abaqusFile", help = "The name of the ABAQUS file to be parsed.", type = str)
    parser.add_argument("exoFile", help = "The name of the exodus file to be written.", type = str)
    parser.add_argument("--force",action="store_true",default=False,help="Overwrite existing files without prompting")
    parser.add_argument("--keepordering",action="store_true",default=False,help="Keep the original ELSET ordering")

    #parser.add_argument("--time_min",type=float,help="first time step",default=0.)
    args = parser.parse_args()
    if  os.path.exists(args.exoFile):
        if args.force:
            os.remove(args.exoFile)
        else:
            if confirm("ExodusII file {0} already exists. Overwrite?".format(args.exoFile)):
                os.remove(args.exoFile)
            else:
                print '\n\t{0} was NOT generated from {1}\n'.format(args.exoFile,args.abaqusFile)
                return -1
    (coord,vertexSet,cellSet,elemNames,nodeNames) = ABAQUSImporter(args.abaqusFile)
    exoWriter(coord,vertexSet,cellSet,args.exoFile,elemNames,nodeNames)

if __name__ == '__main__':
    Main()
