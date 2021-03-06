import sys
if sys.version_info.major == 3:
    import exodus3 as exo
else:
    import exodus2 as exo

def EXOCellSize(cellType):
    if cellType.upper() in ("BAR","BAR2","BEAM"):
        return 2
    elif cellType.upper() in ("BAR3","BEAM3"):
        return 3
    elif cellType.upper() in ("TRI","TRI3","TRIANGLE","TRISHELL","TRISHELL3"):
        return 3
    elif cellType.upper() in ("TRI6","TRISHELL6"):
        return 6
    elif cellType.upper() in ("QUAD","QUAD4","SHELL","SHELL4"):
        return 4
    elif cellType.upper() in ("QUAD8","SHELL8"):
        return 8
    elif cellType.upper() in ("QUAD9","SHELL9"):
        return 9
    elif cellType.upper() in ("TETRA","TETRA4"):
        return 4
    elif cellType.upper() in ("TETRA10",):
        return 10
    elif cellType.upper() in ("HEX","HEX8"):
        return 8
    elif cellType.upper() in ("HEX20",):
        return 20
    elif cellType.upper() in ("HEX27",):
        return 27
    else:
        print("Unknown cell type {0}".format(cellType))
        return -1

def EXODUSwrite(coords,vertexSets,cellSets,numDim,exoFile):
    #''' 
    #   Writes an exodus file
    #'''
    cell1D = ("BAR","BAR2","BEAM2","BAR3","BEAM3")
    cell2D = ("TRI","TRI3","TRIANGLE","TRISHELL","TRISHELL3","TRI6","TRISHELL6","QUAD","QUAD4","SHELL","SHELL4","QUAD9","SHELL9")
    cell3D = ("TETRA","TETRA4","TETRA10","HEX","HEX8","HEX27")



    # Reorder blocks with cells of coDimension 0 first
    if numDim == 2:
        cellCoDim0 = cell2D
        cellCoDim1 = cell1D
        cellCoDim2 = ()
    else:
        cellCoDim0 = cell3D
        cellCoDim1 = cell2D
        cellCoDim2 = cell1D

    #in 3D, cellCoDim2 cells need to be converted into a vertex set
    convertedSets = []
    for setID in cellSets.keys():
        if cellSets[setID]['elemType'].upper() in cellCoDim2:
            print("Converting set {0} of codimension 2 into a vertex set".format(setID))
            vs = []
            for v in cellSets[setID]['connect']:
                if v not in vs:
                    vs.append(v)
            if setID in vertexSets.keys():
                print("Codimension 2 cell set {0} renamed vertex set {1} so as not to clash with existing vertex set".format(setID,max(vertexSets.keys())+1))
                setID = max(vertexSets.keys())+1
            vertexSets[setID] = {}
            vertexSets[setID]['vertex'] = vs
            vertexSets[setID]['name']   = cellSets[setID]['name']
            convertedSets.append(setID)
    for setID in convertedSets:
        cellSets.pop(setID,'None')

    # reorder cells so that sets of codimension 0 are written first in the mesh
    blocksOrder = []
    for setID in cellSets.keys():
        if cellSets[setID]['elemType'].upper() in cellCoDim0:
            blocksOrder.append(setID)
    for setID in cellSets.keys():
        if cellSets[setID]['elemType'].upper() in cellCoDim1:
            blocksOrder.append(setID)

    #Writting exodusII file
    numElems = 0
    for k in cellSets.keys():       #finding number of elements
        numElems += len(cellSets[k]['connect'])//EXOCellSize(cellSets[k]['elemType'])
    numNodes = len(coords[:,0])
    ex_pars = exo.ex_init_params(num_dim=numDim, num_nodes=numNodes,
        num_elem=numElems, num_elem_blk=len(cellSets), num_node_sets=len(vertexSets),num_side_sets=0)
    e=exo.exodus(exoFile,array_type='numpy',mode='w',init_params=ex_pars)
    #coordinates
    if numDim == 3:
        e.put_coord_names(["x","y","z"])
    else:
        e.put_coord_names(["x","y"])
    e.put_coords(coords[:,0],coords[:,1],coords[:,2])

    for setID in blocksOrder:
        elemType = cellSets[setID]['elemType']
        ###---setID, elemType, num elems, num nodes per elem, num attributes per elem
        e.put_elem_blk_info(setID,elemType,len(cellSets[setID]['connect'])//EXOCellSize(elemType),EXOCellSize(elemType),0)
        ###---setID, connectivity table
        e.put_elem_connectivity(setID,cellSets[setID]['connect'])
        if not cellSets[setID]['name'] == '':
            e.put_elem_blk_name(setID,cellSets[setID]['name'])
    for setID in vertexSets.keys():
        ###---setID, num nodes, num distribution factors in a node set
        e.put_node_set_params(setID,len(vertexSets[setID]['vertex']),0)
        ###---setID, nodes
        e.put_node_set(setID,vertexSets[setID]['vertex'])
        if not vertexSets[setID]['name'] == '':
            e.put_node_set_name(setID,vertexSets[setID]['name'])
    ### Adding a QA record, needed until visit fixes its exodus reader
    import datetime
    import os.path
    import sys
    QA_rec_len = 32
    QA = [os.path.basename(sys.argv[0]),os.path.basename(__file__),datetime.date.today().strftime('%Y%m%d'),datetime.datetime.now().strftime("%H:%M:%S")]
    e.put_qa_records([[ q[0:31] for q in QA],])
    e.close()
