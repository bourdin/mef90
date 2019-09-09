import numpy as np
import io
def GMSHread(filename):
    """
    Reads a mesh in gmsh v2. format
    Parameters
    ----------
    filename : name of the mesh

    returns:
       coord:     a numpy array of vertex coordinates
       vertexSet: a dictionary of vertex sets properties
       cellSet:   a dictionary of cell sets properties, inclusing connectivity table
       numDim:    dimensionality of the mesh
    """
    def FixOrientation(connect,coord):
        e0 = (coord[1][0] - coord[0][0],coord[1][1] - coord[0][1])
        e1 = (coord[2][0] - coord[0][0],coord[2][1] - coord[0][1])
        det = e0[0] * e1[1] - e0[1]*e1[0]
        if det <= 0.:
            return (connect[0],connect[2],connect[1])
        else:
            return connect

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
    line = f.readline().strip()
    cellSet = {}
    vertexSet = {}
    physicalNames = {}

    while line != '':
        #print line
        if line.startswith('$MeshFormat'):
            line = f.readline().strip()
            if not line.startswith('2.'):
                raise(RuntimeError,"gmsh2exo.py can only read msh2 format, but got format {0}.\nUse the option -format msh2".format(line.strip()))
            while not line.startswith('$EndMeshFormat'):
                line = f.readline()
        if line.startswith('$PhysicalNames'):
            nNames = int(f.readline())
            for v in range(nNames):
                line = f.readline()
                physicalNames[int(line.split()[1])] = line.split()[2].strip('"')
            line = f.readline().strip()
            if not line.startswith('$EndPhysicalNames'):
                raise(RuntimeError,'something went wrong. Was expecting $EndPhysicalNames, and got {0} instead.'.format(line))
        if line.startswith('$Nodes'):
            nVert = int(f.readline())
            coord = []
            for v in range(nVert):
                coord.append([float(x) for x in f.readline().strip().split()[1:]])
            coord = np.array(coord)
            line = f.readline().strip()
            if not line.startswith('$EndNodes'):
                raise(RuntimeError,'something went wrong. Was expecting $EndNodes, and got {0} instead.'.format(line))
        if line.startswith('$Elements'):
            nCell = int(f.readline())
            vertexSet = {}
            cellSet   = {}
            for c in range(nCell):
                tmp = [int(x) for x in f.readline().strip().split()[1:]]
                cellType = tmp[0]       #first number is cell type
                if cellType == 15:      #this is when element is a node
                    ### vertex sets:
                    tag = tmp[2]        #third number in tmp is a tag
                    if not tag in vertexSet.keys():         #adding new tags
                        vertexSet[tag] = {}
                        vertexSet[tag]['vertex'] = []
                    vertexSet[tag]['vertex'].append(tmp[-1])          #setting element to tag
                else:
                    ### cell sets:
                    tag = tmp[2]                    #third number in tmp is a tag
                    if not tag in cellSet.keys():   #looking for tag in cellSet  
                        cellSet[tag] = {}           #creating new tag if doesn't exist
                        cellSet[tag]['connect'] = []
                        cellSet[tag]['elemType'] = elemType[cellType]   #adding element type to tag
                    elif not cellSet[tag]['elemType'] == elemType[cellType]:
                        raise(RuntimeError,'ID {0} is shared by more than 1 cell set.\n This is probably due to Physical volume, surface, lines sharing the same ID '.format(tag))
                    cellConnect = tmp[-elemDim[cellType]:]
                    cellConnect = reorder(cellType,cellConnect)
                    if elemType[cellType] == 'TRI3':
                        locCoord = (coord[cellConnect[0]-1],coord[cellConnect[1]-1],coord[cellConnect[2]-1])
                        cellConnect = FixOrientation(cellConnect,locCoord)
                    cellSet[tag]['connect'] += cellConnect  #reordering and adding to array        #adding connectivity table to tag

            line = f.readline().strip()
            if not line.startswith('$EndElements'):
                raise(RuntimeError,'something went wrong. Was expecting $EndElements, and got {0} instead.'.format(line))
        line = f.readline().strip()

    for set in vertexSet:
        if set in physicalNames:
            vertexSet[set]['name'] = physicalNames[set]
        else:
            vertexSet[set]['name'] = ''
    for set in cellSet:
        if set in physicalNames:
            cellSet[set]['name'] = physicalNames[set]
        else:
            cellSet[set]['name'] = ''

    if max(coord[:,2]) == min(coord[:,2]):
        numDim = 2
    else:
        numDim = 3
    return coord,vertexSet,cellSet,numDim

