import numpy as np
import io

def ABAQUSread(filename):
    '''
    Import an abaqus file
    '''
    def ABAQUScellType(cellType):
        if cellType.upper() in ("B21", "B21H", "T2D2"):
            return "BAR2"
        elif cellType.upper() in ("S3R","STRI3","CPE3", "CPE3R", "CPE3S", "CPS3", "CPS3R", "CPS3S"):
            return "TRI3"
        elif cellType.upper() in ("CPE6", "CPE6R", "CPE6S", "CPS6", "CPS6R", "CPS6S"):
            return "TRI6"
        elif cellType.upper() in ("S4R", "CPE4", "CPE4S", "CPE4R", "CPS4", "CPS4S", "CPS4R"):
            return "QUAD4"
        elif cellType.upper() in ("CPE8", "CPE8S", "CPE8R", "CPS8", "CPS8S", "CPS8R"):
            return "QUAD8"
        elif cellType.upper() in ("C3D4"):
            return "TETRA4"
        elif cellType.upper() in ("C3D10", "C3D10R", "C3D10S"):
            return "TETRA10"
        elif cellType.upper() in ("C3D8", "C3D8R", "C3D8S"):
            return "HEX8"
        elif cellType.upper() in ("C3D20", "C3D20R", "C3D20S"):
            return "HEX20"

    def ABAQUSreadCoords(f,line):
        #creating and filling list of coordinates from abaqus file
        coord = []
        line = f.readline()
        #skip lines describing nodeset, and coords list ends on '**'
        while not line.startswith('*'):
            coord.append([float(x) for x in line.strip().split(", ")[1:]])
            line = f.readline()
        coord = np.array(coord)         #making coord list into numpy array
        if max(coord[:,2]) == min(coord[:,2]):
            numDim = 2
        else:
            numDim = 3
        return coord,numDim,line
        
    def ABAQUSreadELSET(f,cellSet,line):
        #creating and filling cellSet dictionary
        #get rid of word '*ELEMENT'
        for s in line.split(",")[1:]:
            if s.upper().startswith('TYPE'):
                cellType = s.split('=')[1].strip()
            if s.upper().startswith('ELSET'):
                setName = s.split('=')[1].strip()
            else:
                setName = ''

        usedID = list(cellSet.keys())+[0,]
        try:
            setID = int(setName)
            if setID in usedID:
                setID =  max(usedID)+1    
        except ValueError:
            setID =  max(usedID)+1

        line = f.readline()     #continue to actual data        
        while not line.startswith('*') and not line == '':
            if line.strip().split(",")[len(line.strip().split(",")[0:])-1] == '':
                cellConnect = [int(x)-1 for x in line.strip().split(",")[1:len(line.strip().split(","))-1]]
                line = f.readline()
                if line.strip().split(",")[len(line.strip().split(",")[0:])-1] == '':
                    tmp = [int(x)-1 for x in line.strip().split(",")[0:len(line.strip().split(","))-1]]
                    for i in range(len(tmp)):
                        cellConnect.append(tmp[i])
                    line = f.readline()
                    tmp = [int(x)-1 for x in line.strip().split(",")[0:]]
                    for i in range(len(tmp)):
                        cellConnect.append(tmp[i])
                else:
                    tmp = [int(x)-1 for x in line.strip().split(",")[0:]]
                    for i in range(len(tmp)):
                        cellConnect.append(tmp[i])
            else:
                cellConnect = [int(x)-1 for x in line.strip().split(",")[1:]]    
            ### create cell sets:
            line = f.readline() 
        cellSet[setID]             = {}
        cellSet[setID]['connect']  = []
        cellSet[setID]['elemType'] = ABAQUScellType(cellType)
        cellSet[setID]['name']     = setName.strip()
        print("Assigning ID {0:4d}, type {2:s} to ELSET {1:s}. mef90/vDef name will be cs{0:04}".format(setID,setName,cellSet[setID]['elemType']))
        cellSet[setID]['connect'] += cellConnect  #reordering and adding to array
        return cellSet,line

    def ABAQUSreadNSET(f,vertexSet,line):
        setName = line.strip().split(",")[1]
        setName = setName.split('=')[1]
        usedID = list(vertexSet.keys())+[0,]
        try:
            setID = int(setName)
            if setID in usedID:
                setID =  max(usedID)+1    
        except ValueError:
            setID =  max(usedID)+1
        vertexSet[setID] = {}
        vertexSet[setID]['name'] = setName
        vertexSet[setID]['vertex'] = []

        #continue to next line which has actual node numbers
        line = f.readline()
        #add nodes to list
        #sometimes the node set will end with '**', and the next one will begin right after
        while not line.startswith('*') or line == '':
            #this condition is for reading lines which have a comma at the end
            #splitting at "," will return an additional string with nothing in it
            if line.strip().split(",")[len(line.strip().split(",")[0:])-1] == '':
                vertexSet[setID]['vertex'] += [int(x)-1 for x in line.strip().split(",")[0:len(line.strip().split(","))-1]]
            else:
                vertexSet[setID]['vertex'] += [int(x)-1 for x in line.strip().split(",")[0:]]
            #move onto next line
            line = f.readline()
        return vertexSet,line

    # Opening and reading abaqus file
    f = io.open(filename, 'r') 
    line = f.readline()
    cellSet = {}
    vertexSet = {}
    while line != '':
        if line.upper().startswith('*NODE'):
            (coord,numDim,line) = ABAQUSreadCoords(f,line)
        elif line.upper().startswith('*ELEM'):
            (cellSet,line) = ABAQUSreadELSET(f,cellSet,line)
        elif line.upper().startswith('*NSET'):
            (vertexSet,line) = ABAQUSreadNSET(f,vertexSet,line)
        else:
            line = f.readline()
    f.close()
    return coord,vertexSet,cellSet,numDim
    
