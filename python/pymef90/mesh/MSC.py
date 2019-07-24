import EXODUS
import numpy as np
import io

def MSCread(filename):
    def MSCIntervalParse(str):
        I = []
        for s in str.upper().split('AND'):
            if 'TO' in s:
                ssplit = s.split('TO')
                if len(ssplit) == 2:
                    I += list(range(int(ssplit[0]),int(ssplit[1])+1))
                else:
                    print('error parsing string {0}'.format(s))
            else:
                for ss in s.split(','):
                    I += [int(ss),]
        return I

    def MSCcellType(cellType):
        if cellType in (1, 5, 9, 39):
            # class 1
            return "BAR2"
        elif cellType in (9, 39):
            # class 2
            return "BAR2"
        elif cellType in (2, 6, 37, 38, 50, 196, 201, 228, 229):
            # class 3
            return "TRI3"
        elif cellType in (3, 10, 11, 18, 19, 20, 39, 40, 80, 81, 82, 83, 111, 112, 160, 161, 162, 198, 230, 231):
            # class 4
            return "QUAD4"
        elif cellType in (7, 43, 84, 113, 163):
            # class 5
            return "HEX8"
        elif cellType in (64, 65):
            # class 6
            return "BAR3"
        elif cellType in (26, 27, 28, 29, 30, 32, 33, 34, 41, 42, 62, 63, 66, 67, 199, 234, 235):
            # class 7
            return "QUAD8"
        elif cellType in (53, 54, 55, 56, 58, 59, 60, 69, 70, 73, 74):
            # class 8
            return "QUAD8"
        elif cellType in (21, 35, 44, 236):
            # class 9
            return "HEX20"
        elif cellType in (57, 61, 71):
            # class 10
            return "HEX20"
        elif cellType in (13, 14, 25, 52):
            # class 13
            return "BAR2"
        elif cellType in (124, 125, 126, 128, 129, 131, 132, 197, 200, 231, 232):
            # class 14
            return "TRI6"
        elif cellType in (127, 130, 133, 233):
            # class 15
            return "TETRA10"
        elif cellType in (114, 115, 116, 118, 119, 121, 122):
            # class 16
            return "QUAD4"
        elif cellType in (117, 120, 123):
            # class 17
            return "HEX8"
        elif cellType in (134, 135, 164):
            # class 18
            return "TETRA4"
        else:
            return "UNKNOWN"

    def MSCcellReorder(cellType):
        if cellType == "BAR3":
            return [0,2,1]
        elif cellType == "HEX20":
            return [0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15]
        else:
            # other cell type are numbered consitently with EXODUS
            return  list(range(EXODUS.EXOCellSize(cellType)))

    def MSCcellGetFace(cellType,faceNum):
        if cellType == "TRI3":
            faceList = [[0,1],
                        [1,2],
                        [2,0]]
        elif cellType == "TRI6":
            faceList = [[0,1,3],
                        [1,2,4],
                        [0,2,5]]
        elif cellType == "TETRA4":
            faceList = [[0,1,2],
                        [0,1,3],
                        [1,2,3],
                        [0,2,3]]
            faceList = [[0, 1, 3],
                        [1, 2, 3],
                        [0, 1, 2],
                        [0, 2, 3]]
        elif cellType == "TETRA10":
            faceList = [[0,1,2,4,5,6],
                        [0,1,3,4,8,7],
                        [1,2,3,8,9,5],
                        [0,2,3,6,0,7]]
        else:
            print("Unknown face list for cell type {0}".format(cellType))
            return [],''
        return faceList[faceNum]

    def MSCcellGetFaceType(cellType):
        if cellType == "TRI3":
            return "BAR2"
        elif cellType == "TRI6":
            return "BAR3"
        elif cellType == "TETRA4":
            return "TRI3"
        elif cellType == "TETRA10":
            return "TRI6"
        else:
            print("Unknown face list for cell type {0}".format(cellType))
            return ""

    def MSCparseFloat(s):
        return float(s.replace('+', 'e+').replace('-','e-'))

    def MSCreadline(f):
        line = f.readline().strip()
        fullLine = ''
        while line.endswith('c'):
            fullLine += line.rstrip('c')
            line = f.readline().strip()
        fullLine += line
        return fullLine

    f = io.open(filename, 'r') 
    nCell = 0
    nVert = 0
    cellSet = {}
    vertexSet = {}
    usedID = [0,]

    line=f.readline().strip()
    while line!="":
        if line.startswith('sizing'):
            nCell = int(line.split()[2])
            nVert = int(line.split()[3])
            print ('Number of cells: {0}, number of vertices: {1}'.format(nCell,nVert))
        elif line.startswith('connectivity'):
            f.readline() # skip one line
            connect = []
            cellType = np.empty(nCell)
            for cell in range(nCell):
                line = MSCreadline(f).split()
                cellType[cell]= int(line[1])
                connect.append([int(c) for c in line[2:]])
        elif line.startswith('coordinates'):
            f.readline() # skip one line
            coord = np.empty([nVert,3])
            for vert in range(nVert):
                line = f.readline()[10:].strip().replace('+', 'e+').replace('-','e-')
                coord[vert,:] = [float(x) for x in line.split()[:]]
        elif line.startswith('define'):
            if line.strip().split()[1].startswith('element'):
                setName = line.split()[-1].strip()
                try:
                    setID = int(setName)
                    if setID in usedID:
                        setID =  max(usedID)+1    
                except ValueError:
                    setID =  max(usedID)+1
                usedID.append(setID)
                cellSet[setID] = {}
                cellSet[setID]['name'] = setName
                cellList = MSCIntervalParse(f.readline().strip())
                cellSet[setID]['elemType'] = MSCcellType(cellType[cellList[0]-1])
                cellOrdering = MSCcellReorder(cellSet[setID]['elemType'])
                cellSet[setID]['connect']  = []
                for cell in cellList:
                    for c in cellOrdering:
                        cellSet[setID]['connect'].append(connect[cell-1][c])
            elif line.strip().split()[1].startswith('facemt'):
                setName = line.split()[-1].strip()
                try:
                    setID = int(setName)
                    if setID in usedID:
                        setID =  max(usedID)+1    
                except ValueError:
                    setID =  max(usedID)+1
                usedID.append(setID)
                cellSet[setID] = {}
                cellSet[setID]['name'] = setName
                print 'Face set {0} is now cell set cs{1:04d}'.format(cellSet[setID]['name'],setID)
                cellSet[setID]['connect']  = []
                line = MSCreadline(f).split()
                setCellType = MSCcellType(cellType[int(line[0].split(':')[0])-1])
                cellSet[setID]['elemType'] = MSCcellGetFaceType(setCellType)
                for face in line:
                    faceConnect = MSCcellGetFace(setCellType,int(face.split(':')[1])-1)
                    for c in faceConnect:
                        cellSet[setID]['connect'].append(connect[int(face.split(':')[0])-1][c])
            elif line.strip().split()[1].startswith('ndsq'):
                setName = line.split()[-1].strip()
                try:
                    setID = int(setName)
                    if setID in usedID:
                        setID =  max(usedID)+1    
                except ValueError:
                    setID =  max(usedID)+1
                usedID.append(setID)
                vertexSet[setID] = {}
                vertexSet[setID]['name'] = setName
                print 'Vertex set {0} is now vertex set vs{1:04d}'.format(vertexSet[setID]['name'],setID)
                line = MSCreadline(f)
                vertexSet[setID]['vertex'] = [int(v) for v in line.split()]
        line=f.readline().strip()
    if max(coord[:,2]) == min(coord[:,2]):
        numDim = 2
    else:
        numDim = 3
    return coord,vertexSet,cellSet,numDim
