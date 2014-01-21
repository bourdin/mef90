EXODUS_PY_COPYRIGHT_AND_LICENSE = """ 

exodus.py v 1.0 (beta) is a python wrapper of some of the exodus II library

Copyright (c) 2013 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the Sandia Corporation nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

EXODUS_PY_VERSION = "1.0 (beta)"

EXODUS_PY_COPYRIGHT = """
You are using exodus.py v 1.0 (beta), a python wrapper of some of the exodus II library.
Copyright (c) 2013 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.
"""

EXODUS_PY_CONTACTS = """
Authors:  Timothy Shelton  (trshelt@sandia.gov)
          Michael Veilleux (mgveill@sandia.gov)
          David Littlewood (djlittl@sandia.gov)
"""

import sys
sys.dont_write_bytecode = True

oneline = "Gather from or export to Exodus II files using the Exodus II library"

# Imports and external programs
from ctypes import *
import os

def basename(file):
  """
    Extract base name from file.
    basename("test.e") -> "test"
  """
  fileParts = file.split(".")
  basename = ".".join(fileParts[:-1])
  return basename

def getExodusVersion():
  """
  Parse the exodusII.h header file and return the version number or 0 if not
  found.
  """
  for line in open(getAccessPath() + '/include/exodusII.h'):
    fields = line.split()
    if (len(fields) == 3 and
        fields[0] == '#define' and
        fields[1] == 'EX_API_VERS_NODOT'):
      return int(fields[2])
  return 0

def getAccessPath():
  """
  Return the base path of the exodus installation.
  """
  accessPth = os.path.join(os.environ.get('PETSC_DIR'),os.environ.get('PETSC_ARCH'))
  return accessPth

def findLibExoPath():
  # FIXME: this is a hack to try to find the ExodusII and NetCDF shared
  #        libraries in an external ExodusII distribution; this could be made
  #        portable if the shared libraries were built in the Sierra distribution
  accessPth = getAccessPath()
  trialShLibDirs = [ '/shlib/', '/lib/shared/', '/lib/' ]
  for shLibDir in trialShLibDirs:
    libExoPth = accessPth + shLibDir
    if os.path.isdir(libExoPth): return libExoPth
  raise Exception, "Shared library sub-directory not found in $ACCESS path: %s" % accessPth

EXODUS_PATH = findLibExoPath()
if os.uname()[0] == 'Darwin':
  NETCDF_SO = EXODUS_PATH + "libnetcdf.dylib"
  EXODUS_SO = EXODUS_PATH + "libexodus.dylib"
else:
  NETCDF_SO = EXODUS_PATH + "libnetcdf.so"
  EXODUS_SO = EXODUS_PATH + "libexodus.so"
NETCDF_LIB = cdll.LoadLibrary(NETCDF_SO)
EXODUS_LIB = cdll.LoadLibrary(EXODUS_SO)
  
# set exodus error output option
exErrPrintMode = c_int(1) # VERBOSE
EXODUS_LIB.ex_opts(exErrPrintMode)

EX_API_VERSION_NODOT = getExodusVersion()
EX_NOCLOBBER         = 0x0004  # does not overwrite existing exodus file
EX_CLOBBER           = 0x0008  # overwrites existing exodus file
MAX_STR_LENGTH       = 32      # match exodus default
MAX_LINE_LENGTH      = 80      # match exodus default
EX_MAPS_INT64_API    = 0x2000  # all maps (id, order, ...) store int64_t values
EX_IDS_INT64_API     = 0x4000  # all entity ids (sets, blocks, maps) are int64_t
EX_BULK_INT64_API    = 0x8000  # all integer bulk data (not ids) are int64_t
EX_INQ_INT64_API     = 0x10000 # integers passed to/from ex_inquire are int64_t

def ex_inquiry(inquiry):
  # create dictionary for return types
  inquiry_dictionary = {
    'EX_INQ_FILE_TYPE': 1,                    # inquire EXODUS II file type
    'EX_INQ_API_VERS': 2,                     # inquire API version number
    'EX_INQ_DB_VERS': 3,                      # inquire database version number
    'EX_INQ_TITLE': 4,                        # inquire database title
    'EX_INQ_DIM': 5,                          # inquire number of dimensions
    'EX_INQ_NODES': 6,                        # inquire number of nodes
    'EX_INQ_ELEM': 7,                         # inquire number of elements
    'EX_INQ_ELEM_BLK': 8,                     # inquire number of element blocks
    'EX_INQ_NODE_SETS': 9,                    # inquire number of node sets
    'EX_INQ_NS_NODE_LEN': 10,                 # inquire length of node set node list
    'EX_INQ_SIDE_SETS': 11,                   # inquire number of side sets
    'EX_INQ_SS_NODE_LEN': 12,                 # inquire length of side set node list
    'EX_INQ_SS_ELEM_LEN': 13,                 # inquire length of side set element list
    'EX_INQ_QA': 14,                          # inquire number of QA records
    'EX_INQ_INFO': 15,                        # inquire number of info records
    'EX_INQ_TIME': 16,                        # inquire number of time steps in the database
    'EX_INQ_EB_PROP': 17,                     # inquire number of element block properties
    'EX_INQ_NS_PROP': 18,                     # inquire number of node set properties
    'EX_INQ_SS_PROP': 19,                     # inquire number of side set properties
    'EX_INQ_NS_DF_LEN': 20,                   # inquire length of node set distribution factor list
    'EX_INQ_SS_DF_LEN': 21,                   # inquire length of side set distribution factor list
    'EX_INQ_LIB_VERS': 22,                    # inquire API Lib vers number
    'EX_INQ_EM_PROP': 23,                     # inquire number of element map properties
    'EX_INQ_NM_PROP': 24,                     # inquire number of node map properties
    'EX_INQ_ELEM_MAP': 25,                    # inquire number of element maps
    'EX_INQ_NODE_MAP': 26,                    # inquire number of node maps
    'EX_INQ_EDGE': 27,                        # inquire number of edges
    'EX_INQ_EDGE_BLK': 28,                    # inquire number of edge blocks
    'EX_INQ_EDGE_SETS': 29,                   # inquire number of edge sets
    'EX_INQ_ES_LEN': 30,                      # inquire length of concat edge set edge list
    'EX_INQ_ES_DF_LEN': 31,                   # inquire length of concat edge set dist factor list
    'EX_INQ_EDGE_PROP': 32,                   # inquire number of properties stored per edge block
    'EX_INQ_ES_PROP': 33,                     # inquire number of properties stored per edge set
    'EX_INQ_FACE': 34,                        # inquire number of faces
    'EX_INQ_FACE_BLK': 35,                    # inquire number of face blocks
    'EX_INQ_FACE_SETS': 36,                   # inquire number of face sets
    'EX_INQ_FS_LEN': 37,                      # inquire length of concat face set face list
    'EX_INQ_FS_DF_LEN': 38,                   # inquire length of concat face set dist factor list
    'EX_INQ_FACE_PROP': 39,                   # inquire number of properties stored per face block
    'EX_INQ_FS_PROP': 40,                     # inquire number of properties stored per face set
    'EX_INQ_ELEM_SETS': 41,                   # inquire number of element sets
    'EX_INQ_ELS_LEN': 42,                     # inquire length of concat element set element list
    'EX_INQ_ELS_DF_LEN': 43,                  # inquire length of concat element set dist factor list
    'EX_INQ_ELS_PROP': 44,                    # inquire number of properties stored per elem set
    'EX_INQ_EDGE_MAP': 45,                    # inquire number of edge maps
    'EX_INQ_FACE_MAP': 46,                    # inquire number of face maps
    'EX_INQ_COORD_FRAMES': 47,                # inquire number of coordinate frames
    'EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH': 48,  # inquire size of MAX_NAME_LENGTH dimension on database
    'EX_INQ_DB_MAX_USED_NAME_LENGTH': 49,     # inquire size of MAX_NAME_LENGTH dimension on database
    'EX_INQ_READ_NAME_LENGTH': 50,            # inquire client-specified max size of returned names
    'EX_INQ_DB_FLOAT_SIZE': 51                # inquire size of floating-point values stored on database
    }
  # search dictionary for the requested code
  if inquiry in inquiry_dictionary:
    return inquiry_dictionary[inquiry]
  # none found, must be invalid
  return -1 # EX_INQ_INVALID

def ex_entity_type(varType):
  entity_dictionary = {
    'EX_ELEM_BLOCK': 1,   #   element block property code
    'EX_NODE_SET': 2,     #   node set property code
    'EX_SIDE_SET': 3,     #   side set property code
    'EX_ELEM_MAP': 4,     #   element map property code
    'EX_NODE_MAP': 5,     #   node map property code
    'EX_EDGE_BLOCK': 6,   #   edge block property code
    'EX_EDGE_SET': 7,     #   edge set property code
    'EX_FACE_BLOCK': 8,   #   face block property code
    'EX_FACE_SET': 9,     #   face set property code
    'EX_ELEM_SET': 10,    #   face set property code
    'EX_EDGE_MAP': 11,    #   edge map property code
    'EX_FACE_MAP': 12,    #   face map property code
    'EX_GLOBAL': 13,      #   global 'block' for variables
    'EX_NODAL': 14,       #   nodal 'block' for variables
    'EX_NODE_BLOCK': 14,  #   alias for EX_NODAL
    'EX_COORDINATE': 15,  #   kluge so some internal wrapper functions work
    }
  # search dictionary for the requested code
  if varType in entity_dictionary:
    return entity_dictionary[varType]
  # none found, must be invalid`
  return -1 ##EX_INVALID

# Class definition

class exodus:
  """
    e = exodus('file.e',mode='r') -> open existing exodus file for data extraction
    e = exodus('file.e',mode='a') -> open existing exodus file for data insertion/extraction
    e = exodus('file.e',mode='w',title=title,
               numDims=numDims,numNodes=numNodes
               numElems=numElems,numBlocks=numBlocks,
               numNodeSets=numNodeSets,
               numSideSets=numSideSets) -> open new exodus file for data insertion
  """
  
  # --------------------------------------------------------------------

  def __init__(self,file,mode=None,title=None,numDims=None,
               numNodes=None,numElems=None,numBlocks=None,
               numNodeSets=None,numSideSets=None,io_size=0):
    print EXODUS_PY_COPYRIGHT
    if mode == None:
      mode = 'r'
    self.EXODUS_LIB = EXODUS_LIB
    self.fileName = file
    self.basename = basename(file)
    self.modeChar = mode
    self.__open(io_size=io_size)
    if mode.lower() == 'w':
      info = [title,numDims,numNodes,numElems,numBlocks,
              numNodeSets,numSideSets]
      assert None not in info
      self.__ex_put_info(info)
      self.numTimes = c_int(0)
    else:
      self.__ex_get_info()
      self.numTimes = c_int(self.__ex_inquire_int(ex_inquiry("EX_INQ_TIME")))

  #
  # copy to a new database
  #
  # --------------------------------------------------------------------

  def copy(self,fileName):
    """
      e.copy(fileName) -> copies Exodus database to fileName and returns
                          this copy as a new exodus object
    """
    new = exodus( fileName,
                  "w",
                  self.title(),
                  self.num_dimensions(),
                  self.num_nodes(),
                  self.num_elems(),
                  self.num_blks(),
                  self.num_node_sets(),
                  self.num_side_sets() )
    self.__copy_file(new)
    return new

  #
  # general info
  #
  # --------------------------------------------------------------------

  def title(self):
    """ 
      e.title() -> title in exodus file
    """
    return self.Title.value

  # --------------------------------------------------------------------

  def version_num(self):
    """
      e.version_num() -> string representation of exodus version number
    """
    return "%1.2f" % self.version.value

  # --------------------------------------------------------------------
  
  def put_info(self,Title,numDim,numNodes,numElem,numElemBlk,numNodeSets,numSideSets):
    """
      e.put_info(self,Title,numDim,numNodes,numElem,numElemBlk,numNodeSets,numSideSets)
      -> put initialization information into exodus file
    """
    self.__ex_put_info([Title,numDim,numNodes,numElem,numElemBlk,numNodeSets,numSideSets])
    return True

  # --------------------------------------------------------------------

  def get_qa_records(self):
    """
      e.get_qa_records() -> get a list of QA records from the exodus database;
                            each QA record is a length-4 tuple of strings:
                               1) the analysis code name
                               2) the analysis code descriptor, e.g. version
                               3) the analysis data
                               4) the analysis time
    """
    return self.__ex_get_qa()

  # --------------------------------------------------------------------

  def put_qa_records(self,records):
    """
      e.put_qa_records(records) -> put a list of QA records into the exodus database;
                                   each QA record must be a length-4 tuple of strings:
                                     1) the analysis code name
                                     2) the analysis code descriptor, e.g. version
                                     3) the analysis data
                                     4) the analysis time
    """
    for rec in records:
      assert len(rec) == 4
      for recEntry in rec:
        assert len(str(recEntry)) < MAX_STR_LENGTH
    if self.__ex_put_qa(records):
      return True
    else:
      return False

  # --------------------------------------------------------------------
  def num_info_records(self):
    return int(self.__ex_inquire_int(ex_inquiry("EX_INQ_INFO")))

  # --------------------------------------------------------------------
  def get_info_records(self):
    """
      e.get_info_records() -> get a list info records from the exodus database;
                              each entry in the list is one line of info, e.g.
                              a line of an input deck
    """
    info_recs = self.__ex_get_info_recs()
    return info_recs

  # --------------------------------------------------------------------

  def put_info_records(self,info):
    """
      e.put_info_records(info) -> put a list of info records into an exodus database;
                                  each entry in the list is one line of info, e.g.
                                  a line of an input deck
    """
    for rec in info:
      if len(str(rec)) > MAX_LINE_LENGTH:
        print "WARNING: max line length reached for one or more info records;"
        print "         info stored to exodus file is incomplete for these records"
        break
    if self.__ex_put_info_recs(info):
      return True
    else:
      return False

  # --------------------------------------------------------------------

  def get_sierra_input(self,inpFileName=None):
    """
      e.get_sierra_input(inpFileName) -> get sierra input deck from exodus database;
                                         if inpFileName is passed, input is written
                                         to this file, e.g. myInputDeck.i; otherwise
                                         a list of file lines is returned
    """
    info_recs = self.__ex_get_info_recs()
    sierra_inp = []
    begin = False
    for rec in info_recs:
      vals = rec.split()
      if not begin: # have not reached Sierra block
        if len(vals) >= 2 and vals[0].lower() == 'begin' and vals[1].lower() == "sierra":
          begin = True
      if begin: # inside Sierra block
        sierra_inp.append(rec)
        if len(rec) > MAX_LINE_LENGTH:
          print "WARNING: max line length reached for one or more input lines;"
          print "         input data might be incomplete for these lines"
          break
        if len(vals) >= 2 and vals[0].lower() == "end" and vals[1].lower() == "sierra":
          break # end of Sierra block
    if inpFileName:
      fd = open(inpFileName,"w")
      for fileLine in sierra_inp:
        print >> fd, fileLine
      fd.close()
      return ""
    else:
      return sierra_inp

  #
  # time steps
  #
  # --------------------------------------------------------------------

  def num_times(self):
    """ 
      e.num_times() -> number of time steps in exodus file
    """
    return self.numTimes.value

  # --------------------------------------------------------------------

  def get_times(self):
    """ 
      e.get_times() -> get list of times in exodus file
    """
    if self.numTimes.value == 0:
      self.times = []
    else:
      self.__ex_get_all_times()
    return self.times

  # --------------------------------------------------------------------

  def put_time(self,step,value):
    """ 
      e.put_time(step,value) -> put time step and value into exodus file
    """
    self.__ex_put_time(step,value)
    return True

  #
  # coordinate system
  #
  # --------------------------------------------------------------------

  def num_dimensions(self):
    """ 
      e.num_dimensions() -> number of dimensions in exodus file
    """
    return self.numDim.value

  # --------------------------------------------------------------------

  def get_coord_names(self):
    """ 
      e.get_coord_names() -> get list of coordinate names in exodus file
    """
    names = self.__ex_get_coord_names()
    return names

  # --------------------------------------------------------------------

  def put_coord_names(self,names):
    """ 
      e.put_coord_names() -> set list of coordinate names in exodus file
    """
    self.__ex_put_coord_names(names)

  #
  # nodes
  #
  # --------------------------------------------------------------------

  def num_nodes(self):
    """ 
      e.num_nodes() -> number of nodes in exodus file
    """
    return self.numNodes.value

  # --------------------------------------------------------------------

  def get_coords(self):
    """ 
      e.get_coords() -> get tuple of lists of coordinates (X,Y,Z) in exodus file
    """
    self.__ex_get_coord()
    return (self.coordsX,self.coordsY,self.coordsZ)

  # --------------------------------------------------------------------

  def get_coord(self,i):
    """
      e.get_coord(i) -> get (x,y,z) of i^th node in exodus file
    """
    listX,listY,listZ = self.__ex_get_n_coord(i,1)
    return (listX[0],listY[0],listZ[0])

  # --------------------------------------------------------------------

  def put_coords(self,xCoords,yCoords,zCoords):
    """ 
      e.put_coords() -> put coordinates (xCoords,yCoords,zCoords) in exodus file
    """
    self.__ex_put_coord(xCoords,yCoords,zCoords)
    return True

  # --------------------------------------------------------------------

  def get_node_num_map(self):
    """ 
      e.get_node_num_map() -> **DEPRECATED** 
                              use: e.ex_get_node_id_map()

                              get list mapping local node id to global id 
                              in exodus file
    """
    nodeNumMap = self.__ex_get_node_num_map()
    return nodeNumMap

  # --------------------------------------------------------------------

  def get_node_id_map(self):
    """ 
      e.get_node_id_map() -> get list mapping local node id to global id 
                             in exodus file
    """
    objType = ex_entity_type("EX_NODE_MAP")
    inqType = ex_inquiry("EX_INQ_NODES")
    nodeIdMap = self.__ex_get_id_map(objType,inqType)
    return nodeIdMap

  # --------------------------------------------------------------------

  def put_node_id_map(self,map):
    """ 
      e.put_node_id_map(map) -> put list mapping local node id to global id 
                                into exodus file
    """
    objType = ex_entity_type("EX_NODE_MAP")
    inqType = ex_inquiry("EX_INQ_NODES")
    return self.__ex_put_id_map(objType,inqType,map)

  # --------------------------------------------------------------------

  def get_node_variable_names(self):
    """ 
      e.get_node_variable_names() -> get list of node variable names in exodus file
    """
    if self.__ex_get_var_param('n').value == 0:
      return []
    return self.__ex_get_var_names("n")

  # --------------------------------------------------------------------

  def get_node_variable_number(self):
    """ 
      e.get_node_variable_number() -> get number of node variables in exodus file
    """
    ndType = ex_entity_type("EX_NODAL")
    num = self.__ex_get_variable_param(ndType)
    return num.value

  # --------------------------------------------------------------------

  def set_node_variable_number(self,number):
    """ 
      e.set_node_variable_number(number) -> set number of node variables in exodus file
    """
    ndType = ex_entity_type("EX_NODAL")
    self.__ex_put_variable_param(ndType,number)
    return True

  # --------------------------------------------------------------------

  def put_node_variable_name(self,name,index):
    """ 
      e.put_node_variable_name("name",index) -> node variable with name at index into exodus file
    """
    ndType = ex_entity_type("EX_NODAL")
    NDvarNames = self.get_node_variable_names()
    if name in NDvarNames:
      print "WARNING:node variable \"", name, "\" already exists."
    if index > len(NDvarNames):
      raise Exception, ("ERROR: variable index out of range.")
    self.__ex_put_variable_name(ndType,index,name)
    return True

  # --------------------------------------------------------------------

  def get_node_variable_values(self,name,step):
    """ 
      e.get_node_variable_values("name",step) -> get list of node variable values
      for a step within exodus file
    """
    names = self.get_node_variable_names()
    var_id = names.index(name) + 1
    ndType = ex_entity_type("EX_NODAL")
    numVals = self.num_nodes()
    values =  self.__ex_get_var(step,ndType,var_id,0,numVals)
    return values

  # --------------------------------------------------------------------

  def put_node_variable_values(self,name,step,values):
    """ 
      e.put_node_variable_values("name",step,values) -> put node values into
      variable name at step into exodus file
    """
    names = self.get_node_variable_names()
    var_id = names.index(name) + 1
    ndType = ex_entity_type("EX_NODAL")
    numVals = self.num_nodes()
    self.__ex_put_var(step,ndType,var_id,0,numVals,values)
    return True

  #
  # elements
  #
  # --------------------------------------------------------------------

  def num_elems(self):
    """ 
      e.num_elems() -> number of elements in exodus file
    """
    return self.numElem.value

  # --------------------------------------------------------------------

  def get_elem_num_map(self):
    """ 
      e.get_elem_num_map() -> **DEPRECATED** 
                              use: e.ex_get_elem_id_map()

                              get list mapping local element id to 
                              global id in exodus file
    """
    elemNumMap = self.__ex_get_elem_num_map()
    return elemNumMap

  # --------------------------------------------------------------------

  def get_elem_id_map(self):
    """ 
      e.get_elem_id_map() -> get list mapping local elem id to global id 
                             in exodus file
    """
    objType = ex_entity_type("EX_ELEM_MAP")
    inqType = ex_inquiry("EX_INQ_ELEM")
    elemIdMap = self.__ex_get_id_map(objType,inqType)
    return elemIdMap

  # --------------------------------------------------------------------

  def put_elem_id_map(self,map):
    """ 
      e.put_elem_id_map(map) -> put list mapping local elem id to global id 
                                into exodus file
    """
    objType = ex_entity_type("EX_ELEM_MAP")
    inqType = ex_inquiry("EX_INQ_ELEM")
    return self.__ex_put_id_map(objType,inqType,map)

  # --------------------------------------------------------------------

  def get_elem_order_map(self):
    """ 
      e.get_elem_order_map() -> get list of optimized element ordering
    """

    elemOrderMap = self.__ex_get_elem_order_map()
    return elemOrderMap

  #
  # element blocks
  #
  # --------------------------------------------------------------------

  def num_blks(self):
    """ 
      e.num_blks() -> number of element blocks in exodus file
    """
    return self.numElemBlk.value

  # --------------------------------------------------------------------

  def get_elem_blk_ids(self):
    """ 
      e.get_elem_blk_ids() -> get list of element block ids in exodus file
    """
    self.__ex_get_elem_blk_ids()
    return self.elemBlkIds

  # --------------------------------------------------------------------

  def get_elem_blk_name(self,id):
    """ 
      e.get_elem_blk_name(id) -> get element block name for block "id" in 
                                 exodus file
    """
    objType = ex_entity_type("EX_ELEM_BLOCK")
    elemBlkName = self.__ex_get_name(objType,id)
    return elemBlkName

  # --------------------------------------------------------------------

  def put_elem_blk_name(self,id,name):
    """ 
      e.put_elem_blk_name(id,name) -> write element block name for block 
                                      "id" in exodus file
    """
    objType = ex_entity_type("EX_ELEM_BLOCK")
    self.__ex_put_name(objType,id,name)

  # --------------------------------------------------------------------

  def get_elem_blk_names(self):
    """ 
      e.get_elem_blk_names() -> get list of element block names in 
                                exodus file
    """
    objType = ex_entity_type("EX_ELEM_BLOCK")
    inqType = ex_inquiry("EX_INQ_ELEM_BLK")
    elemBlkNames = self.__ex_get_names(objType,inqType)
    return elemBlkNames

  # --------------------------------------------------------------------

  def put_elem_blk_names(self,names):
    """ 
      e.put_elem_blk_names(names) -> write list of element block names to 
                                     exodus file
    """
    objType = ex_entity_type("EX_ELEM_BLOCK")
    inqType = ex_inquiry("EX_INQ_ELEM_BLK")
    self.__ex_put_names(objType,inqType,names)

  # --------------------------------------------------------------------

  def elem_blk_info(self,id):
    """
    e.elem_blk_info(id) -> element block info for block "id" in exodus file
    """
    (elemType,numElem,nodesPerElem,numAttr) = self.__ex_get_elem_block(id)
    return elemType.value, numElem.value, nodesPerElem.value, numAttr.value

  # --------------------------------------------------------------------

  def put_elem_blk_info(self,id,elemType,numElems,
                        numNodesPerElem,numAttrsPerElem):
    """
    e.put_elem_blk_info(id) -> sets ID, element type string (all caps), 
                               number of elements, number of nodes per element,
                               and number of attributes per element in an
                               element block
    """
    self.__ex_put_elem_block(id,elemType,numElems,
                               numNodesPerElem,numAttrsPerElem)

  # --------------------------------------------------------------------

  def put_concat_elem_blk(self,elemBlkIDs, elemType, numElemThisBlk,\
                          numNodesPerElem,numAttr,defineMaps):
    """
      e.put_concat_elem_blk(self,elemBlkIDs, elemType, numElemThisBlk,... \
                       numNodesPerElem,numAttr,defineMaps)
    """
    self.__ex_put_concat_elem_blk(elemBlkIDs,elemType,numElemThisBlk, \
               numNodesPerElem,numAttr,defineMaps)
    return True

  # --------------------------------------------------------------------

  def get_elem_connectivity(self,id):
    """ 
      e.get_elem_connectivity(id) -> get tuple of one list and two reals 
      (elem_block_connectivity,num_elem_this_blk,num_nodes_per_elem)
    """
    (elem_block_connectivity,num_elem_this_blk,num_nodes_per_elem) = self.__ex_get_elem_conn(id);
    return (elem_block_connectivity,num_elem_this_blk.value,num_nodes_per_elem.value)

  # --------------------------------------------------------------------

  def put_elem_connectivity(self,id,connectivity):
    """ 
      e.put_elem_connectivity() -> set the element connectivity array
                                   for all elements in a block.  The
                                   information for this block must have
                                   already been set by calling
                                   e.put_elem_blk_info() on this block.
    """
    d1,numBlkElems,numNodesPerElem,d2 = self.elem_blk_info(id)
    assert len(connectivity) == (numBlkElems * numNodesPerElem)
    self.__ex_put_elem_conn(id,connectivity)

  # --------------------------------------------------------------------

  def get_elem_attr(self,elemBlkID):
    """
      e.get_elm_attr(elemBlkID) -> get element block attributes
    """
    attribute = self.__ex_get_elem_attr(elemBlkID)
    return attribute

  # --------------------------------------------------------------------

  def put_elem_attr(self,elemBlkID,Attr):
    """
      e.put_elem_attr(elemBlkID,Attr) -> write element attribute lists for each block
    """
    self.__ex_put_elem_attr(elemBlkID,Attr)

  # --------------------------------------------------------------------

  def elem_type(self,id):
    """ 
      e.elem_type(id) -> type of element, i.e "QUAD", "TETRA", etc.
    """
    (elemType,numElem,nodesPerElem,numAttr) = self.__ex_get_elem_block(id)
    return elemType.value

  # --------------------------------------------------------------------

  def num_attr(self,id):
    """ 
      e.num_attr(id) -> number of element attributes
    """
    (elemType,numElem,nodesPerElem,numAttr) = self.__ex_get_elem_block(id)
    return numAttr.value

  # --------------------------------------------------------------------

  def num_elems_in_blk(self,id):
    """ 
      e.num_elems_in_blk(id) -> number of elements in block "id" in exodus file
    """
    (elemType,numElem,nodesPerElem,numAttr) = self.__ex_get_elem_block(id)
    return numElem.value

  # --------------------------------------------------------------------

  def num_nodes_per_elem(self,id):
    """ 
      e.num_nodes_per_elem(id) -> number of nodes per element in block "id" in exodus file
    """
    (elemType,numElem,nodesPerElem,numAttr) = self.__ex_get_elem_block(id)
    return nodesPerElem.value

  # --------------------------------------------------------------------

  def get_element_variable_truth_table(self,blockId=None):
    """
      e.get_element_variable_truth_table(blockId) -> gets a truth table indicating which variables are
                                                     defined for a block.  If blockId is not passed, then
                                                     a concatenated truth table for all blocks is returned
                                                     with variable index cycling faster than blockId
    """
    truthTable = self.__ex_get_elem_var_tab()
    if blockId != None:
      self.get_elem_blk_ids()
      assert blockId in list(self.elemBlkIds)
      indx = list(self.elemBlkIds).index(blockId)
      numVars = self.__ex_get_var_param("e").value
      start,stop = (indx * numVars, (indx + 1) * numVars)
      return truthTable[start:stop]
    return truthTable

  # --------------------------------------------------------------------

  def set_element_variable_truth_table(self,table):
    """
      e.set_element_variable_truth_table(table) -> sets a truth table indicating which variables are
                                                   defined for each block in the model.  
                                                   Element variable index cycles faster than the 
                                                   element block index
    """
    self.get_elem_blk_ids()
    numBlks = len(self.elemBlkIds)
    numVars = int(self.__ex_get_var_param("e").value)
    assert len(table) == (numBlks * numVars)
    return self.__ex_put_elem_var_tab(table)

  # --------------------------------------------------------------------

  def get_element_variable_values(self,blockId,name,step):
    """ 
      e.get_element_variable_values(blockId,"name",step) -> get list of element variable values
      for an element block id at step within exodus file
    """
    names = self.get_element_variable_names()
    var_id = names.index(name) + 1
    ebType = ex_entity_type("EX_ELEM_BLOCK")
    numVals = self.num_elems_in_blk(blockId)
    values =  self.__ex_get_var(step,ebType,var_id,blockId,numVals)
    return values

  # --------------------------------------------------------------------

  def put_element_variable_values(self,blockId,name,step,values):
    """ 
      e.put_element_variable_values(blockId,"name",step,values) -> put values into element block id
      and variable name at step into exodus file
    """
    names = self.get_element_variable_names()
    var_id = names.index(name) + 1
    ebType = ex_entity_type("EX_ELEM_BLOCK")
    numVals = self.num_elems_in_blk(blockId)
    self.__ex_put_var(step,ebType,var_id,blockId,numVals,values)
    return True

  # --------------------------------------------------------------------

  def get_element_variable_number(self):
    """ 
      e.get_element_variable_number() -> get number of element variables in exodus file
    """
    ebType = ex_entity_type("EX_ELEM_BLOCK")
    num = self.__ex_get_variable_param(ebType)
    return num.value

  # --------------------------------------------------------------------

  def set_element_variable_number(self,number):
    """ 
      e.set_element_variable_number(number) -> set number of element variables in exodus file
    """
    ebType = ex_entity_type("EX_ELEM_BLOCK")
    self.__ex_put_variable_param(ebType,number)
    return True

  # --------------------------------------------------------------------

  def get_element_variable_names(self):
    """ 
      e.get_element_variable_names() -> get list of element variable names in exodus file
    """
    if self.__ex_get_var_param("e").value == 0:
      return []
    return self.__ex_get_var_names("e")

  # --------------------------------------------------------------------

  def put_element_variable_name(self,name,index):
    """ 
      e.put_element_variable_name("name",index) -> element variable with name at index into exodus file
    """
    ebType = ex_entity_type("EX_ELEM_BLOCK")
    EBvarNames = self.get_element_variable_names()
    if name in EBvarNames:
      print "WARNING:element variable \"", name, "\" already exists."
    if index > len(EBvarNames):
      print "index", index, "len", len(EBvarNames)
      raise Exception, ("ERROR: variable index out of range.")
    self.__ex_put_variable_name(ebType,index,name)
    return True

  # --------------------------------------------------------------------

  def get_element_attribute_names(self,blkId):
    """ 
      e.get_element_attribute_names(blId) -> get list of element attribute names for block in exodus file
    """
    names = self.__ex_get_elem_attr_names(blkId)
    return list(names)

  # --------------------------------------------------------------------

  def put_element_attribute_names(self,blkId,names):
    """ 
      e.put_element_attribute_names(blId,names) -> set element attribute names for block in exodus file
    """
    return self.__ex_put_elem_attr_names(blkId,names)

  # --------------------------------------------------------------------

  def get_element_property_names(self):
    """ 
      e.get_element_property_names() -> get list of element property names in exodus file
    """
    names = []
    ebType = ex_entity_type("EX_ELEM_BLOCK")
    inqType = "EX_INQ_EB_PROP"
    names = self.__ex_get_prop_names(ebType,inqType)
    return list(names)

  # --------------------------------------------------------------------

  def get_element_property_value(self,id,name):
    """ 
      e.get_element_property_value(id,name) -> get element property value in exodus file
                                               for an element block ID and property name
    """
    ebType = ex_entity_type("EX_ELEM_BLOCK")
    propVal = self.__ex_get_prop(ebType,id,name)
    return int(propVal)

  # --------------------------------------------------------------------

  def put_element_property_value(self,id,name,value):
    """
      e.put_element_property_value(id,name,value) -> put an element property name and
                                                     its integer value for a block ID
                                                     into an exodus file
    """
    ebType = ex_entity_type("EX_ELEM_BLOCK")
    if self.__ex_put_prop(ebType,id,name,value):
      return True
    else:
      return False

  #
  # nodesets
  #
  # --------------------------------------------------------------------

  def num_node_sets(self):
    """ 
      e.num_node_sets() -> number of node sets in exodus file
    """
    return self.numNodeSets.value

  # --------------------------------------------------------------------

  def get_node_set_ids(self):
    """ 
      e.get_node_set_ids() -> get list of node set ids in exodus file
    """
    self.__ex_get_node_set_ids()
    return list(self.nodeSetIds)

  # --------------------------------------------------------------------

  def get_node_set_name(self,id):
    """ 
      e.get_node_set_name(id) -> get node set name for node set "id" in 
                                 exodus file
    """
    objType = ex_entity_type("EX_NODE_SET")
    nodeSetName = self.__ex_get_name(objType,id)
    return nodeSetName

  # --------------------------------------------------------------------

  def put_node_set_name(self,id,name):
    """ 
      e.put_node_set_name(id,name) -> write node set name for node set 
                                      "id" in exodus file
    """
    objType = ex_entity_type("EX_NODE_SET")
    self.__ex_put_name(objType,id,name)

  # --------------------------------------------------------------------

  def get_node_set_names(self):
    """ 
      e.get_node_set_names() -> get list of node set names in 
                                exodus file
    """
    objType = ex_entity_type("EX_NODE_SET")
    inqType = ex_inquiry("EX_INQ_NODE_SETS")
    nodeSetNames = self.__ex_get_names(objType,inqType)
    return nodeSetNames

  # --------------------------------------------------------------------

  def put_node_set_names(self,names):
    """ 
      e.put_node_set_names(names) -> write list of node set names to 
                                     exodus file
    """
    objType = ex_entity_type("EX_NODE_SET")
    inqType = ex_inquiry("EX_INQ_NODE_SETS")
    self.__ex_put_names(objType,inqType,names)

  # --------------------------------------------------------------------

  def num_nodes_in_node_set(self,id):
    """ 
      e.num_nodes_in_node_set(id) -> number of nodes in node set id
    """
    node_set_nodes = self.get_node_set_nodes(id)
    return len(node_set_nodes)

  # --------------------------------------------------------------------

  def get_node_set_nodes(self,id):
    """
      e.get_node_set_nodes(id) -> list of node id's in node set
    """
    node_set_ids = self.get_node_set_ids()
    assert id in node_set_ids
    node_set_nodes = self.__ex_get_node_set(id)
    return list(node_set_nodes)

  # --------------------------------------------------------------------

  def put_node_set(self,id,nodeSetNodes):
    """
      e.put_node_set() -> set the id and node list for a node set (creates the node set).
    """
    self.__ex_put_node_set(id,nodeSetNodes)

  # --------------------------------------------------------------------

  def get_node_set_dist_facts(self,id):
    """
      e.get_node_set_nodes(id) -> list of distribution factors for nodes in node set
    """
    node_set_dfs = self.__ex_get_node_set_dist_fact(id)
    return list(node_set_dfs)

  # --------------------------------------------------------------------

  def put_node_set_dist_fact(self,id,nodeSetDistFact):
    """
      e.put_node_set_dist_fact() -> sets the distribution factors for a node set.
    """
    self.__ex_put_node_set_dist_fact(id,nodeSetDistFact)

  # --------------------------------------------------------------------

  def get_node_set_variable_number(self):
    """ 
      e.get_node_set_variable_number() -> get number of node set variables in exodus file
    """
    nsType = ex_entity_type("EX_NODE_SET")
    num = self.__ex_get_variable_param(nsType)
    return num.value

  # --------------------------------------------------------------------

  def set_node_set_variable_number(self,number):
    """ 
      e.set_node_set_variable_number(number) -> set number of node set variables in exodus file
    """
    nsType = ex_entity_type("EX_NODE_SET")
    self.__ex_put_variable_param(nsType,number)
    return True

  # --------------------------------------------------------------------

  def get_node_set_variable_truth_table(self,nodeSetId=None):
    """
      e.get_node_set_variable_truth_table(nodeSetId) -> gets a truth table indicating which variables are
                                                     defined for a node set.  If nodeSetId is not passed, then
                                                     a concatenated truth table for all node sets is returned
                                                     with variable index cycling faster than nodeSetId
    """
    truthTable = self.__ex_get_nset_var_tab()
    if nodeSetId != None:
      self.get_node_set_ids()
      assert nodeSetId in list(self.nodeSetIds)
      indx = list(self.nodeSetIds).index(nodeSetId)
      numVars = self.__ex_get_var_param("m").value
      start,stop = (indx * numVars, (indx + 1) * numVars)
      return truthTable[start:stop]
    return truthTable

  # --------------------------------------------------------------------

  def set_node_set_variable_truth_table(self,table):
    """
      e.set_node_set_variable_truth_table(table) -> sets a truth table indicating which variables are
                                                   defined for each node set in the model.  
                                                   Node set variable index cycles faster than the 
                                                   node set index
    """
    self.get_node_set_ids()
    numBlks = len(self.nodeSetIds)
    numVars = int(self.__ex_get_var_param("m").value)
    assert len(table) == (numBlks * numVars)
    return self.__ex_put_nset_var_tab(table)

  # --------------------------------------------------------------------

  def get_node_set_variable_names(self):
    """ 
      e.get_node_set_variable_names() -> get list of node set variable names in exodus file
    """
    names = []
    nsType = ex_entity_type("EX_NODE_SET")
    num_vars = self.__ex_get_variable_param(nsType)
    for varid in range(num_vars.value):
      varid += 1
      name = self.__ex_get_variable_name(nsType,varid)
      names.append(name.value)
    return names

  # --------------------------------------------------------------------

  def put_node_set_variable_name(self,name,index):
    """ 
      e.put_node_set_variable_name("name",index) -> put node set variable with name at index into exodus file
    """
    nsType = ex_entity_type("EX_NODE_SET")
    NSvarNames = self.get_node_set_variable_names()
    if name in NSvarNames:
      print "WARNING: Node set variable \"", name, "\" already exists."
    if index > len(NSvarNames):
      raise Exception, ("ERROR: variable index out of range.")
    self.__ex_put_variable_name(nsType,index,name)
    return True

  # --------------------------------------------------------------------

  def get_node_set_variable_values(self,id,name,step):
    """ 
      e.get_node_set_variable_names(id,"name",step) -> get list of node set variable values 
      for node set id at step within exodus file
    """
    names = self.get_node_set_variable_names()
    var_id = names.index(name) + 1
    values =  self.__ex_get_nset_var(step,var_id,id)
    return values

  # --------------------------------------------------------------------

  def put_node_set_variable_values(self,id,name,step,values):
    """ 
      e.put_node_set_variable_values(id,"name",step,values) -> put values into node set id 
      and variable name at step into exodus file
    """
    names = self.get_node_set_variable_names()
    var_id = names.index(name) + 1
    self.__ex_put_nset_var(step,var_id,id,values)
    return True

  # --------------------------------------------------------------------

  def get_all_node_set_params(self):
    """
      e.get_all_node_set_params() -> get total number of nodes and distribution factors
                                     combined in all node sets
    """
    self.__ex_get_node_set_ids()
    totNumSetNodes, totNumSetDistFacts = 0, 0
    for nodeSetId in self.nodeSetIds: 
      (numSetNodes,numSetDistFacts) = self.__ex_get_node_set_param(int(nodeSetId))
      totNumSetNodes += numSetNodes
      totNumSetDistFacts += numSetDistFacts
    return (totNumSetNodes, totNumSetDistFacts)

  # --------------------------------------------------------------------

  def get_node_set_params(self,id):
    """
      e.get_node_set_params() -> get number of nodes and distribution factors
                                     in a node set
    """
    (numSetNodes,numSetDistFacts) = self.__ex_get_node_set_param(int(id))
    return (numSetNodes, numSetDistFacts)

  # --------------------------------------------------------------------

  def put_node_set_params(self,id,numSetNodes,numSetDistFacts=None):
    """
      e.put_node_set_params() -> sets the ID, number of nodes, and number of
                                 distribution factors in a node set.  If the
                                 last argument is not given, then the number
                                 of distribution factors is set equal to the
                                 number of nodes
    """
    if numSetDistFacts == None: numSetDistFacts = numSetNodes
    assert numSetDistFacts == 0 or numSetDistFacts == numSetNodes
    self.__ex_put_node_set_param(id,numSetNodes,numSetDistFacts)

  # --------------------------------------------------------------------

  def get_node_set_property_names(self):
    """ 
      e.get_node_set_property_names() -> get list of nodeset property names in exodus file
    """
    names = []
    nsType = ex_entity_type("EX_NODE_SET")
    inqType = "EX_INQ_NS_PROP"
    names = self.__ex_get_prop_names(nsType,inqType)
    return list(names)

  # --------------------------------------------------------------------

  def get_node_set_property_value(self,id,name):
    """ 
      e.get_node_set_property_value(id,name) -> get nodeset property value in exodus file
                                                for an nodeset ID and property name
    """
    nsType = ex_entity_type("EX_NODE_SET")
    propVal = self.__ex_get_prop(nsType,id,name)
    return int(propVal)

  # --------------------------------------------------------------------

  def put_node_set_property_value(self,id,name,value):
    """
      e.put_node_set_property_value(id,name,value) -> put a nodeset property name and
                                                      its integer value for a nodeset
                                                      ID into an exodus file
    """
    nsType = ex_entity_type("EX_NODE_SET")
    if self.__ex_put_prop(nsType,id,name,value):
      return True
    else:
      return False

  #
  # sidesets
  #
  # --------------------------------------------------------------------

  def num_side_sets(self):
    """ 
      e.num_side_sets() -> number of side sets in exodus file
    """
    return self.numSideSets.value

  # --------------------------------------------------------------------

  def get_side_set_ids(self):
    """ 
      e.get_side_set_ids() -> get list of side set ids in exodus file
    """
    self.__ex_get_side_set_ids()
    return self.sideSetIds

  # --------------------------------------------------------------------

  def get_side_set_name(self,id):
    """ 
      e.get_side_set_name(id) -> get side set name for side set "id" in 
                                 exodus file
    """
    objType = ex_entity_type("EX_SIDE_SET")
    sideSetName = self.__ex_get_name(objType,id)
    return sideSetName

  # --------------------------------------------------------------------

  def put_side_set_name(self,id,name):
    """ 
      e.put_side_set_name(id,name) -> write side set name for side set 
                                      "id" in exodus file
    """
    objType = ex_entity_type("EX_SIDE_SET")
    self.__ex_put_name(objType,id,name)

  # --------------------------------------------------------------------

  def get_side_set_names(self):
    """ 
      e.get_side_set_names() -> get list of side set names in 
                                exodus file
    """
    objType = ex_entity_type("EX_SIDE_SET")
    inqType = ex_inquiry("EX_INQ_SIDE_SETS")
    sideSetNames = self.__ex_get_names(objType,inqType)
    return sideSetNames

  # --------------------------------------------------------------------

  def put_side_set_names(self,names):
    """ 
      e.put_side_set_names(names) -> write list of side set names to 
                                     exodus file
    """
    objType = ex_entity_type("EX_SIDE_SET")
    inqType = ex_inquiry("EX_INQ_SIDE_SETS")
    self.__ex_put_names(objType,inqType,names)

  # --------------------------------------------------------------------

  def num_faces_in_side_set(self,id):
    """ 
      e.num_faces_in_side_set(id) -> number of faces in side set id
    """
    ssids = self.get_side_set_ids()
    if ( id not in ssids ): 
      print "WARNING: queried side set ID does not exist in database"
      return 0
    (num_side_in_set,num_dist_fact_in_set) = self.__ex_get_side_set_param(id)
    return num_side_in_set

  # --------------------------------------------------------------------

  def get_all_side_set_params(self):
    """
      e.get_all_side_set_params() -> get total number of elements, nodes, and distribution
                                     factors combined in all side sets
    """
    self.__ex_get_side_set_ids()
    totNumSetSides, totNumSetDistFacts = 0, 0 # totNumSetDistFacts = totNumSetNodes
    for sideSetId in self.sideSetIds: 
      (numSetSides,numSetDistFacts) = self.__ex_get_side_set_param(int(sideSetId))
      totNumSetSides += numSetSides
      totNumSetDistFacts += numSetDistFacts
    totNumSetNodes = totNumSetDistFacts
    return (totNumSetSides, totNumSetNodes, totNumSetDistFacts)

  # --------------------------------------------------------------------

  def get_side_set_params(self,id):
    """
      e.get_side_set_params() -> get number of sides, and distribution factors
                                 in a side set
    """
    (numSetSides,numSetDistFacts) = self.__ex_get_side_set_param(int(id))
    return (numSetSides, numSetDistFacts)

  # --------------------------------------------------------------------

  def put_side_set_params(self,id,numSetSides,numSetDistFacts):
    """
      e.put_side_set_params() -> set ID, num elements, and num nodes of a sideset
    """
    self.__ex_put_side_set_param(id,numSetSides,numSetDistFacts)

  # --------------------------------------------------------------------

  def get_side_set(self,id):
    """ 
      e.get_side_set(id) -> get tuple of two lists (side_set_elem_list,
      side_set_side_list) the side_set_elem_list contains the elements in the
      side set, and the side_set_side_list contains the side ID for each of the
      elements in the side_set_elem_list (identifies the side onto which the load
      is applied).
    """
    (side_set_elem_list,side_set_side_list) = self.__ex_get_side_set(id)
    return (side_set_elem_list,side_set_side_list)

  # --------------------------------------------------------------------

  def put_side_set(self,id,sideSetElements,sideSetSides):
    """
      e.put_side_set() -> set the id, element ids, and side ids for
                          a side set (creates the side set).
    """
    self.__ex_put_side_set(id,sideSetElements,sideSetSides)

  # --------------------------------------------------------------------

  def get_side_set_dist_fact(self,id):
    """ 
      e.get_side_set_dist_fact(id) ->  list of distribution factors for sides in side set
    """
    side_set_dfs = self.__ex_get_side_set_dist_fact(id)
    return list(side_set_dfs)

  # --------------------------------------------------------------------

  def put_side_set_dist_fact(self,id,sideSetDistFact):
    """
      e.put_side_set_dist_fact() -> sets the distribution factors for a side set.
                                    A distribution factor is given for each node
                                    on each face in the side set.
    """
    self.__ex_put_side_set_dist_fact(id,sideSetDistFact)

  # --------------------------------------------------------------------  

  def get_side_set_node_list(self,id):
    """ 
      e.get_side_set_nodes(id) -> get tuple of two lists (side_set_node_cnt_list,
      side_set_node_list) the side_set_node_cnt_list is the number of nodes on each
      face and the side_set_node_list is the list of local node ids for the side set 
      in exodus file
    """
    (side_set_node_cnt_list,side_set_node_list) = self.__ex_get_side_set_node_list(id)
    return (side_set_node_cnt_list,side_set_node_list)

  # --------------------------------------------------------------------

  def get_side_set_variable_truth_table(self,sideSetId=None):
    """
      e.get_side_set_variable_truth_table(sideSetId) -> gets a truth table indicating which variables are
                                                     defined for a side set.  If sideSetId is not passed, then
                                                     a concatenated truth table for all side sets is returned
                                                     with variable index cycling faster than sideSetId
    """
    truthTable = self.__ex_get_sset_var_tab()
    if sideSetId != None:
      self.get_side_set_ids()
      assert sideSetId in list(self.sideSetIds)
      indx = list(self.sideSetIds).index(sideSetId)
      numVars = self.__ex_get_var_param("s").value
      start,stop = (indx * numVars, (indx + 1) * numVars)
      return truthTable[start:stop]
    return truthTable

  # --------------------------------------------------------------------

  def set_side_set_variable_truth_table(self,table):
    """
      e.set_side_set_variable_truth_table(table) -> sets a truth table indicating which variables are
                                                   defined for each side set in the model.  
                                                   Side set variable index cycles faster than the 
                                                   side set index
    """
    self.get_side_set_ids()
    numBlks = len(self.sideSetIds)
    numVars = int(self.__ex_get_var_param("s").value)
    assert len(table) == (numBlks * numVars)
    return self.__ex_put_sset_var_tab(table)

  # --------------------------------------------------------------------

  def get_side_set_variable_number(self):
    """ 
      e.get_side_set_variable_number() -> get number of side set variables in exodus file
    """
    ssType = ex_entity_type("EX_SIDE_SET")
    num = self.__ex_get_variable_param(ssType)
    return num.value

  # --------------------------------------------------------------------

  def set_side_set_variable_number(self,number):
    """ 
      e.set_side_set_variable_number(number) -> set number of side set variables in exodus file
    """
    ssType = ex_entity_type("EX_SIDE_SET")
    self.__ex_put_variable_param(ssType,number)
    return True

  # --------------------------------------------------------------------

  def get_side_set_variable_names(self):
    """ 
      e.get_side_set_variable_names() -> get list of side set variable names in exodus file
    """
    names = []
    ssType = ex_entity_type("EX_SIDE_SET")
    num_vars = self.__ex_get_variable_param(ssType)
    for varid in range(num_vars.value):
      varid += 1
      name = self.__ex_get_variable_name(ssType,varid)
      names.append(name.value)
    return names

  # --------------------------------------------------------------------

  def put_side_set_variable_name(self,name,index):
    """ 
      e.put_side_set_variable_name("name",index) -> put side set variable with name at index into exodus file
    """
    ssType = ex_entity_type("EX_SIDE_SET")
    SSvarNames = self.get_side_set_variable_names()
    if name in SSvarNames:
      print "WARNING:Side set variable \"", name, "\" already exists."
    if index > len(SSvarNames):
      raise Exception, ("ERROR: variable index out of range.")
    self.__ex_put_variable_name(ssType,index,name)
    return True

  # --------------------------------------------------------------------

  def get_side_set_variable_values(self,id,name,step):
    """ 
      e.get_side_set_variable_names(id,"name",step) -> get list of side set variable values 
      for side set id at step within exodus file
    """
    names = self.get_side_set_variable_names()
    var_id = names.index(name) + 1
    values =  self.__ex_get_sset_var(step,var_id,id)
    return values

  # --------------------------------------------------------------------

  def put_side_set_variable_values(self,id,name,step,values):
    """ 
      e.put_side_set_variable_values(id,"name",step,values) -> put values into side set id 
      and variable name at step into exodus file
    """
    names = self.get_side_set_variable_names()
    var_id = names.index(name) + 1
    self.__ex_put_sset_var(step,var_id,id,values)
    return True

  # --------------------------------------------------------------------

  def get_side_set_property_names(self):
    """ 
      e.get_side_set_property_names() -> get list of sideset property names in exodus file
    """
    names = []
    ssType = ex_entity_type("EX_SIDE_SET")
    inqType = "EX_INQ_SS_PROP"
    names = self.__ex_get_prop_names(ssType,inqType)
    return list(names)

  # --------------------------------------------------------------------

  def get_side_set_property_value(self,id,name):
    """ 
      e.get_side_set_property_value(id,name) -> get sideset property value in exodus file
                                                for an sideset ID and property name
    """
    ssType = ex_entity_type("EX_SIDE_SET")
    propVal = self.__ex_get_prop(ssType,id,name)
    return int(propVal)

  # --------------------------------------------------------------------

  def put_side_set_property_value(self,id,name,value):
    """
      e.put_side_set_property_value(id,name,value) -> put a sideset property name and
                                                      its integer value for a sideset
                                                      ID into an exodus file
    """
    ssType = ex_entity_type("EX_SIDE_SET")
    if self.__ex_put_prop(ssType,id,name,value):
      return True
    else:
      return False

  #
  # global variables
  #
  # --------------------------------------------------------------------

  def get_global_variable_number(self):
    """ 
      e.get_global_variable_number() -> get number of global variables in exodus file
    """
    gbType = ex_entity_type("EX_GLOBAL")
    num = self.__ex_get_variable_param(gbType)
    return num.value

  # --------------------------------------------------------------------

  def set_global_variable_number(self,number):
    """ 
      e.set_global_variable_number(number) -> set number of global variables in exodus file
    """
    gbType = ex_entity_type("EX_GLOBAL")
    self.__ex_put_variable_param(gbType,number)
    return True

  # --------------------------------------------------------------------

  def get_global_variable_names(self):
    """ 
      e.get_global_variable_names() -> get list of global variable names in exodus file
    """
    if self.get_global_variable_number() == 0:
      return []
    return self.__ex_get_var_names("g")

  # --------------------------------------------------------------------

  def put_global_variable_name(self,name,index):
    """ 
      e.put_global_variable_name("name",index) -> put global variable with name at index into exodus file
    """
    gbType = ex_entity_type("EX_GLOBAL")
    GlobVarNames = self.get_global_variable_names()
    if name in GlobVarNames:
      print "WARNING: global variable \"", name, "\" already exists."
    if index > len(GlobVarNames):
      print "index", index, "len", len(GlobVarNames)
      raise Exception, ("ERROR: variable index out of range.")
    self.__ex_put_variable_name(gbType,index,name)
    return True

  # --------------------------------------------------------------------

  def get_global_variable_value(self,name,step):
    """ 
      e.get_global_variable_value("name",step) -> get global variable value
      at time step within exodus file
    """
    names = self.get_global_variable_names()
    var_id = names.index(name)
    gbType = ex_entity_type("EX_GLOBAL")
    num = self.__ex_get_variable_param(gbType)
    gvalues =  self.__ex_get_var(step,gbType,0,1,num.value)
    return gvalues[var_id]

  # --------------------------------------------------------------------

  def get_all_global_variable_values(self,step):
    """ 
      e.get_all_global_variable_values(step) -> get all global variable values
      at time step within exodus file
    """
    gbType = ex_entity_type("EX_GLOBAL")
    num = self.__ex_get_variable_param(gbType)
    gvalues =  self.__ex_get_var(step,gbType,0,1,num.value)
    values = []
    for i in xrange(num.value):
      values.append(gvalues[i])
    return values

  # --------------------------------------------------------------------

  def put_global_variable_value(self,name,step,value):
    """ 
      e.put_global_variable_value("name",step,value) -> put global variable value
      and variable name at time step into exodus file
    """
    # we must write all values at once, not individually
    names   = self.get_global_variable_names()
    # get all values
    numVals = self.get_global_variable_number()
    values = (c_double * numVals)()
    for i in xrange(numVals):
      values[i] = c_double(self.get_global_variable_value(names[i], step))
    # adjust one of them
    values[names.index(name)] = c_double(value)
    # write them all
    EXODUS_LIB.ex_put_glob_vars(self.fileId,
                                c_int(step),
                                c_int(numVals),
                                values)
    return True

  # --------------------------------------------------------------------

  def put_all_global_variable_values(self,step,values):
    """ 
      e.put_all_global_variable_values(step,values) -> put all global variable values
      at time step into exodus file
    """
    numVals = self.get_global_variable_number()
    gvalues = (c_double * numVals)()
    for i in xrange(numVals):
      gvalues[i] = c_double(values[i])
    EXODUS_LIB.ex_put_glob_vars(self.fileId,
                                c_int(step),
                                c_int(numVals),
                                gvalues)
    return True

  # --------------------------------------------------------------------

  def get_global_variable_values(self,name):
    """ 
      e.get_global_variable_values("name",step) -> get global variable values
      within exodus file
    """
    names = self.get_global_variable_names()
    var_id = names.index(name)
    gbType = ex_entity_type("EX_GLOBAL")
    num = self.__ex_get_variable_param(gbType)
    values = []
    for i in range(self.numTimes.value):
      gvalues =  self.__ex_get_var(i+1,gbType,0,1,num.value)
      values.append( gvalues[var_id] )
    return values

  # --------------------------------------------------------------------

  def close(self):
    """ 
      e.close() -> close the exodus file.
    """
    print "Closing exodus file: " + self.fileName
    errorInt = EXODUS_LIB.ex_close(self.fileId)
    if errorInt != 0:
      raise Exception, ("ERROR: Closing file " + self.fileName + " had problems.")

  # --------------------------------------------------------------------
  #
  # Exodus API calls
  #
  # --------------------------------------------------------------------

  def __open(self, io_size=0):
    print "Opening exodus file: " + self.fileName
    self.mode = 0
    if self.modeChar.lower() == "a": self.mode = 1
    if self.modeChar.lower() in ["a","r"] and not os.path.isfile(self.fileName):
      raise Exception, ("ERROR: Cannot open " + self.fileName + " for read. Does not exist.")
    elif self.modeChar.lower() == "w" and os.path.isfile(self.fileName):
      raise Exception, ("ERROR: Cowardly not opening " + self.fileName + \
                        " for write. File already exists.")
    elif self.modeChar.lower() not in ["a","r","w"]:
      raise Exception, ("ERROR: File open mode " + self.modeChar + " unrecognized.")
    self.comp_ws = c_int(8)
    self.io_ws = c_int(io_size)
    self.version = c_float(0.0)
    if self.modeChar.lower() in ["a","r"]: # open existing file
      self.fileId = EXODUS_LIB.ex_open_int(self.fileName,self.mode,
                                           byref(self.comp_ws),
                                           byref(self.io_ws),
                                           byref(self.version),
                                           EX_API_VERSION_NODOT)
    else: # create file
      if io_size == 0:
        io_size = 8
        self.io_ws = c_int(io_size)
      self.__create()

  # --------------------------------------------------------------------

  def __create(self):
    cMode = c_int(EX_NOCLOBBER)
    self.fileId = EXODUS_LIB.ex_create_int(self.fileName,cMode,
                                           byref(self.comp_ws),
                                           byref(self.io_ws),
                                           EX_API_VERSION_NODOT)

  # --------------------------------------------------------------------

  def __copy_file(self,other):
    EXODUS_LIB.ex_copy(self.fileId,other.fileId)

  # --------------------------------------------------------------------

  def __ex_get_info(self):
    self.Title       = create_string_buffer(MAX_LINE_LENGTH+1)
    if (EXODUS_LIB.ex_int64_status(self.fileId) & EX_BULK_INT64_API):
      self.numDim      = c_longlong(0)
      self.numNodes    = c_longlong(0)
      self.numElem     = c_longlong(0)
      self.numElemBlk  = c_longlong(0)
      self.numNodeSets = c_longlong(0)
      self.numSideSets = c_longlong(0)
    else:
      self.numDim      = c_int(0)
      self.numNodes    = c_int(0)
      self.numElem     = c_int(0)
      self.numElemBlk  = c_int(0)
      self.numNodeSets = c_int(0)
      self.numSideSets = c_int(0)
    EXODUS_LIB.ex_get_init(self.fileId,self.Title,byref(self.numDim),byref(self.numNodes),
                           byref(self.numElem),byref(self.numElemBlk),byref(self.numNodeSets),
                           byref(self.numSideSets))

  # --------------------------------------------------------------------

  def __ex_put_info(self,info):
    self.Title       = create_string_buffer(info[0],MAX_LINE_LENGTH+1)
    self.numDim      = c_longlong(info[1])
    self.numNodes    = c_longlong(info[2])
    self.numElem     = c_longlong(info[3])
    self.numElemBlk  = c_longlong(info[4])
    self.numNodeSets = c_longlong(info[5])
    self.numSideSets = c_longlong(info[6])
    EXODUS_LIB.ex_put_init(self.fileId,self.Title,self.numDim,self.numNodes,self.numElem,
                           self.numElemBlk,self.numNodeSets,self.numSideSets)
    self.version = self.__ex_inquire_float(ex_inquiry("EX_INQ_DB_VERS"))

  # --------------------------------------------------------------------

  def __ex_put_concat_elem_blk(self,elemBlkIDs, elemType, numElemThisBlk,\
                                 numNodesPerElem,numAttr,defineMaps):
    if (EXODUS_LIB.ex_int64_status(self.fileId) & EX_IDS_INT64_API):
      elem_blk_ids = (c_longlong * len(elemBlkIDs))()
      elem_blk_ids[:] = elemBlkIDs
      num_elem_this_blk = (c_longlong * len(elemBlkIDs))()
      num_elem_this_blk[:] = numElemThisBlk
      num_nodes_per_elem = (c_longlong * len(elemBlkIDs))()
      num_nodes_per_elem[:] = numNodesPerElem
      num_attr = (c_longlong * len(elemBlkIDs))()
      num_attr[:] = numAttr
    else:
      elem_blk_ids = (c_int * len(elemBlkIDs))()
      elem_blk_ids[:] = elemBlkIDs
      num_elem_this_blk = (c_int * len(elemBlkIDs))()
      num_elem_this_blk[:] = numElemThisBlk
      num_nodes_per_elem = (c_int * len(elemBlkIDs))()
      num_nodes_per_elem[:] = numNodesPerElem
      num_attr = (c_int * len(elemBlkIDs))()
      num_attr[:] = numAttr
    elem_type = (c_char_p * len(elemBlkIDs))()
    elem_type[:] = elemType
    define_maps = c_int(defineMaps)
    EXODUS_LIB.ex_put_concat_elem_block(self.fileId,elem_blk_ids,elem_type, \
               num_elem_this_blk,num_nodes_per_elem,num_attr,define_maps)
    


  # --------------------------------------------------------------------

  def __ex_get_qa(self):
    num_qa_recs = c_int(self.__ex_inquire_int(ex_inquiry("EX_INQ_QA")))
    qa_rec_ptrs = ((POINTER(c_char * (MAX_STR_LENGTH+1)) * 4) * num_qa_recs.value)()
    for i in range(num_qa_recs.value):
      for j in range(4):
        qa_rec_ptrs[i][j] = pointer(create_string_buffer(MAX_STR_LENGTH+1))
    if num_qa_recs.value:
      EXODUS_LIB.ex_get_qa(self.fileId,byref(qa_rec_ptrs))
    qa_recs = []
    for qara in qa_rec_ptrs:
      qa_rec_list = []
      for ptr in qara:
        qa_rec_list.append(ptr.contents.value)
      qa_rec_tuple = tuple(qa_rec_list)
      assert len(qa_rec_tuple) == 4
      qa_recs.append(qa_rec_tuple)
    return qa_recs

  # --------------------------------------------------------------------

  def __ex_put_qa(self,qaRecs):
    num_qa_recs = c_int(len(qaRecs))
    qa_rec_ptrs = ((POINTER(c_char * (MAX_STR_LENGTH+1)) * 4) * num_qa_recs.value)()
    for i in range(num_qa_recs.value):
      for j in range(4):
        qa_rec_ptrs[i][j] = pointer(create_string_buffer(str(qaRecs[i][j]),MAX_STR_LENGTH+1))
    EXODUS_LIB.ex_put_qa(self.fileId,num_qa_recs,byref(qa_rec_ptrs))
    return True

  # --------------------------------------------------------------------

  def _ex_get_info_recs_quietly(self):
    num_infos = c_int(self.__ex_inquire_int(ex_inquiry("EX_INQ_INFO")))
    info_ptrs = (POINTER(c_char * (MAX_LINE_LENGTH+1)) * num_infos.value)()
    for i in range(num_infos.value):
      info_ptrs[i] = pointer(create_string_buffer(MAX_LINE_LENGTH+1))
    if num_infos.value:
      EXODUS_LIB.ex_get_info(self.fileId,byref(info_ptrs))
    info_recs = []
    for irp in info_ptrs:
      info_recs.append(irp.contents.value)
    return info_recs

  # --------------------------------------------------------------------

  def __ex_get_info_recs(self):
    num_infos = c_int(self.__ex_inquire_int(ex_inquiry("EX_INQ_INFO")))
    info_ptrs = (POINTER(c_char * (MAX_LINE_LENGTH+1)) * num_infos.value)()
    for i in range(num_infos.value):
      info_ptrs[i] = pointer(create_string_buffer(MAX_LINE_LENGTH+1))
    EXODUS_LIB.ex_get_info(self.fileId,byref(info_ptrs))
    info_recs = []
    for irp in info_ptrs:
      info_recs.append(irp.contents.value)
    for rec in info_recs:
      if len(rec) > MAX_LINE_LENGTH:
        print "WARNING: max line length reached for one or more info records;"
        print "         info might be incomplete for these records"
        break
    return info_recs

  # --------------------------------------------------------------------

  def __ex_put_info_recs(self,infoRecs):
    num_infos = c_int(len(infoRecs))
    info_ptrs = (POINTER(c_char * (MAX_LINE_LENGTH+1)) * num_infos.value)()
    for i in range(num_infos.value):
      info_ptrs[i] = pointer(create_string_buffer(str(infoRecs[i]),MAX_LINE_LENGTH+1))
    EXODUS_LIB.ex_put_info(self.fileId,num_infos,None) ## Define Number of Info
    EXODUS_LIB.ex_put_info(self.fileId,num_infos,byref(info_ptrs)) ## Set Info
    return True

  # --------------------------------------------------------------------

  def __ex_inquire_float(self,id):
    val        = c_int(0)
    dummy_char = create_string_buffer(MAX_LINE_LENGTH+1)
    ret_float  = c_float(0.0)
    if (EXODUS_LIB.ex_int64_status(self.fileId) & EX_INQ_INT64_API):
      dummy_int = c_longlong(0)
    else:
      dummy_int = c_int(0)
    val = EXODUS_LIB.ex_inquire(self.fileId,id,byref(dummy_int),byref(ret_float),dummy_char)
    if val < 0:
      raise Exception, ("ERROR: ex_inquire(" + str(id) + ") failed on " + self.fileName)
    return ret_float

  # --------------------------------------------------------------------

  def __ex_inquire_int(self,id):
    val = c_longlong(0)
    val = EXODUS_LIB.ex_inquire_int(self.fileId,id)
    if val < 0:
      raise Exception, ("ERROR: ex_inquire_int(" + str(id) + ") failed on " + self.fileName)
    return val

  # --------------------------------------------------------------------

  def __ex_get_coord_names(self):
    coord_name_ptrs = (POINTER(c_char * (MAX_STR_LENGTH+1)) * self.numDim.value)()
    for i in range(self.numDim.value):
      coord_name_ptrs[i] = pointer(create_string_buffer(MAX_STR_LENGTH+1))
    EXODUS_LIB.ex_get_coord_names(self.fileId,byref(coord_name_ptrs))
    coord_names = []
    for cnp in coord_name_ptrs: 
      coord_names.append(cnp.contents.value)
    return coord_names

  # --------------------------------------------------------------------

  def __ex_put_coord_names(self,names):
    coord_name_ptrs = (POINTER(c_char * (MAX_STR_LENGTH+1)) * self.numDim.value)()
    assert len(names) == self.numDim.value
    for i in range(self.numDim.value):
      coord_name_ptrs[i] = pointer(create_string_buffer(names[i],MAX_STR_LENGTH+1))
    EXODUS_LIB.ex_put_coord_names(self.fileId,byref(coord_name_ptrs))

  # --------------------------------------------------------------------

  def __ex_get_all_times(self):
    self.times = (c_double * self.numTimes.value)()
    EXODUS_LIB.ex_get_all_times(self.fileId,byref(self.times))

  # --------------------------------------------------------------------

  def __ex_get_time(self,timeStep):
    time_step = c_int(timeStep)
    time_val = c_double(0.0)
    EXODUS_LIB.ex_get_time(self.fileId,time_step,byref(time_val))
    return time_val.value()

  # --------------------------------------------------------------------

  def __ex_put_time(self,timeStep,timeVal):
    time_step = c_int(timeStep)
    time_val = c_double(timeVal)
    EXODUS_LIB.ex_put_time(self.fileId,time_step,byref(time_val))
    return True

  # --------------------------------------------------------------------

  def __ex_get_name(self,objType,objId):
    obj_type = c_int(objType)
    obj_id = c_int(objId)
    obj_name = create_string_buffer(MAX_STR_LENGTH+1)
    EXODUS_LIB.ex_get_name(self.fileId,obj_type,obj_id,byref(obj_name))
    return obj_name.value

  # --------------------------------------------------------------------

  def __ex_put_name(self,objType,objId,objName):
    obj_type = c_int(objType)
    obj_id = c_int(objId)
    obj_name = create_string_buffer(objName,MAX_STR_LENGTH+1)
    EXODUS_LIB.ex_put_name(self.fileId,obj_type,obj_id,obj_name)

  # --------------------------------------------------------------------

  def __ex_get_names(self,objType,inqType):
    obj_type = c_int(objType)
    num_objs = c_int(self.__ex_inquire_int(inqType))
    numObjs = num_objs.value
    obj_name_ptrs = (POINTER(c_char * (MAX_STR_LENGTH+1)) * numObjs)()
    for i in range(numObjs):
      obj_name_ptrs[i] = pointer(create_string_buffer(MAX_STR_LENGTH+1))
    EXODUS_LIB.ex_get_names(self.fileId,obj_type,byref(obj_name_ptrs))
    obj_names = []
    for onp in obj_name_ptrs: 
      obj_names.append(onp.contents.value)
    return obj_names

  # --------------------------------------------------------------------

  def __ex_put_names(self,objType,inqType,objNames):
    num_objs = c_int(self.__ex_inquire_int(inqType))
    numObjs = num_objs.value
    assert numObjs == len(objNames)
    obj_name_ptrs = (POINTER(c_char * (MAX_STR_LENGTH+1)) * numObjs)()
    obj_type = c_int(objType)
    for i in range(numObjs):
      obj_name_ptrs[i] = pointer(create_string_buffer(objNames[i],MAX_STR_LENGTH+1))
    EXODUS_LIB.ex_put_names(self.fileId,obj_type,byref(obj_name_ptrs))

  # --------------------------------------------------------------------

  def __ex_get_elem_blk_ids(self):
    if (EXODUS_LIB.ex_int64_status(self.fileId) & EX_IDS_INT64_API):
      self.elemBlkIds = (c_longlong * self.numElemBlk.value)()
    else:
      self.elemBlkIds = (c_int * self.numElemBlk.value)()
    if self.numElemBlk.value > 0:
      EXODUS_LIB.ex_get_elem_blk_ids(self.fileId,byref(self.elemBlkIds))

  # --------------------------------------------------------------------

  def __ex_get_side_set_ids(self):
    if (EXODUS_LIB.ex_int64_status(self.fileId) & EX_IDS_INT64_API):
      self.sideSetIds = (c_longlong * self.numSideSets.value)()
    else:
      self.sideSetIds = (c_int * self.numSideSets.value)()
    if self.num_side_sets() > 0:
      EXODUS_LIB.ex_get_side_set_ids(self.fileId,byref(self.sideSetIds))

  # --------------------------------------------------------------------

  def __ex_get_node_set_ids(self):
    if (EXODUS_LIB.ex_int64_status(self.fileId) & EX_IDS_INT64_API):
      self.nodeSetIds = (c_longlong * self.numNodeSets.value)()
    else:
      self.nodeSetIds = (c_int * self.numNodeSets.value)()
    if self.num_node_sets() > 0:
      EXODUS_LIB.ex_get_node_set_ids(self.fileId,byref(self.nodeSetIds))

  # --------------------------------------------------------------------

  def __ex_get_node_set_param(self, nodeSetId):
    node_set_id = c_longlong(nodeSetId)
    if (EXODUS_LIB.ex_int64_status(self.fileId) & EX_BULK_INT64_API):
      num_set_nodes      = c_longlong(0)
      num_set_dist_facts = c_longlong(0)
    else:
      num_set_nodes      = c_int(0)
      num_set_dist_facts = c_int(0)
    EXODUS_LIB.ex_get_node_set_param(self.fileId,node_set_id,byref(num_set_nodes),byref(num_set_dist_facts))
    return (int(num_set_nodes.value),int(num_set_dist_facts.value))

  # --------------------------------------------------------------------

  def __ex_put_node_set_param(self,nodeSetId,numNodes,numDistFacts):
    node_set_id        = c_longlong(nodeSetId)
    num_set_nodes      = c_longlong(numNodes)
    num_set_dist_facts = c_longlong(numDistFacts)
    EXODUS_LIB.ex_put_node_set_param(self.fileId,node_set_id,num_set_nodes,num_set_dist_facts)

  # --------------------------------------------------------------------

  def __ex_get_node_set(self, nodeSetId):
    node_set_id = c_longlong(nodeSetId)
    num_node_set_nodes = self.__ex_get_node_set_param(nodeSetId)[0]
    if num_node_set_nodes == 0:
      return []
    if (EXODUS_LIB.ex_int64_status(self.fileId) & EX_BULK_INT64_API):
      set_nodes = (c_longlong * num_node_set_nodes)()
    else:
      set_nodes = (c_int * num_node_set_nodes)()
    EXODUS_LIB.ex_get_node_set(self.fileId,node_set_id,byref(set_nodes))
    return set_nodes

  # --------------------------------------------------------------------

  def __ex_put_node_set(self,nodeSetId,nodeSetNodes):
    node_set_id    = c_longlong(nodeSetId)
    if (EXODUS_LIB.ex_int64_status(self.fileId) & EX_BULK_INT64_API):
      node_set_nodes = (c_longlong * len(nodeSetNodes))()
      for i in range(len(nodeSetNodes)):
        node_set_nodes[i] = c_longlong(nodeSetNodes[i])
    else:
      node_set_nodes = (c_int * len(nodeSetNodes))()
      for i in range(len(nodeSetNodes)):
        node_set_nodes[i] = c_int(nodeSetNodes[i])
    EXODUS_LIB.ex_put_node_set(self.fileId,node_set_id,node_set_nodes)

  # --------------------------------------------------------------------

  def __ex_get_node_set_dist_fact(self, nodeSetId):
    node_set_id = c_longlong(nodeSetId)
    num_node_set_nodes = self.__ex_get_node_set_param(nodeSetId)[0]
    set_dfs = (c_double * num_node_set_nodes)()
    EXODUS_LIB.ex_get_node_set_dist_fact(self.fileId,node_set_id,byref(set_dfs))
    return set_dfs

  # --------------------------------------------------------------------

  def __ex_put_node_set_dist_fact(self,nodeSetId,nodeSetDistFact):
    node_set_id        = c_longlong(nodeSetId)
    node_set_dist_fact = (c_double * len(nodeSetDistFact))()
    for i in range(len(nodeSetDistFact)):
      node_set_dist_fact[i] = c_double(nodeSetDistFact[i])
    EXODUS_LIB.ex_put_node_set_dist_fact(self.fileId,node_set_id,node_set_dist_fact)

  # --------------------------------------------------------------------

  def __ex_get_nset_var(self,timeStep,varId,id):
    step        = c_int(timeStep)
    var_id      = c_int(varId)
    node_set_id = c_longlong(id)
    (numNodeInSet,numDistFactInSet) = self.__ex_get_node_set_param(id)
    num_node_in_set = c_longlong(numNodeInSet)
    ns_var_vals = (c_double * numNodeInSet)()
    EXODUS_LIB.ex_get_nset_var(self.fileId,step,var_id,node_set_id,num_node_in_set,ns_var_vals)
    return list(ns_var_vals)

  # --------------------------------------------------------------------

  def __ex_get_nset_var_tab(self):
    self.__ex_get_node_set_ids()
    node_set_count = c_int(len(self.nodeSetIds))
    variable_count = self.__ex_get_var_param("m")
    truth_table = (c_int * (node_set_count.value * variable_count.value))()
    EXODUS_LIB.ex_get_nset_var_tab(self.fileId,
                                   node_set_count,
                                   variable_count,
                                   byref(truth_table))
    truthTab = []
    for val in truth_table:
      if val:
        truthTab.append(True)
      else:
        truthTab.append(False)
    return truthTab

  # --------------------------------------------------------------------

  def __ex_put_nset_var_tab(self,truthTab):
    self.__ex_get_node_set_ids()
    num_blks = c_int(len(self.nodeSetIds))
    num_vars = self.__ex_get_var_param("m")
    truth_tab = (c_int * (num_blks.value*num_vars.value))()
    for i in xrange(len(truthTab)):
      boolVal = truthTab[i]
      if boolVal: 
        truth_tab[i] = c_int(1)
      else: 
        truth_tab[i] = c_int(0)
    EXODUS_LIB.ex_put_nset_var_tab(self.fileId,num_blks,num_vars,truth_tab)
    return True

  # --------------------------------------------------------------------

  def __ex_put_nset_var(self,timeStep,varId,id,values):
    step        = c_int(timeStep)
    var_id      = c_int(varId)
    node_set_id = c_longlong(id)
    (numNodeInSet,numDistFactInSet) = self.__ex_get_node_set_param(id)
    num_node_in_set = c_longlong(numNodeInSet)
    ns_var_vals = (c_double * numNodeInSet)()
    for i in range(numNodeInSet):
      ns_var_vals[i] = float(values[i])
    EXODUS_LIB.ex_put_nset_var(self.fileId,step,var_id,node_set_id,num_node_in_set,ns_var_vals)
    return True

  # --------------------------------------------------------------------

  def __ex_get_coord(self):
    self.coordsX = (c_double * self.numNodes.value)()
    self.coordsY = (c_double * self.numNodes.value)()
    self.coordsZ = (c_double * self.numNodes.value)()
    EXODUS_LIB.ex_get_coord(self.fileId,byref(self.coordsX),byref(self.coordsY),byref(self.coordsZ))

  # --------------------------------------------------------------------

  def __ex_put_coord(self,xCoords,yCoords,zCoords):
    self.coordsX = (c_double * self.numNodes.value)()
    self.coordsY = (c_double * self.numNodes.value)()
    self.coordsZ = (c_double * self.numNodes.value)()
    for i in range(self.numNodes.value):
      self.coordsX[i] = float(xCoords[i])
      self.coordsY[i] = float(yCoords[i])
      self.coordsZ[i] = float(zCoords[i])
    EXODUS_LIB.ex_put_coord(self.fileId,byref(self.coordsX),byref(self.coordsY),byref(self.coordsZ))

  # --------------------------------------------------------------------

  def __ex_get_n_coord(self,startNodeId,numNodes):
    start_node_num = c_longlong(startNodeId)
    num_nodes      = c_longlong(numNodes)
    coordsX        = (c_double * numNodes)()
    coordsY        = (c_double * numNodes)()
    coordsZ        = (c_double * numNodes)()
    EXODUS_LIB.ex_get_n_coord(self.fileId,start_node_num,num_nodes,byref(coordsX),byref(coordsY),byref(coordsZ))
    return list(coordsX),list(coordsY),list(coordsZ)

  # --------------------------------------------------------------------

  def __ex_get_id_map(self,objType,inqType):
    obj_type = c_int(objType)
    num_objs = c_int(self.__ex_inquire_int(inqType))
    numObjs = num_objs.value
    id_map = []
    if (EXODUS_LIB.ex_int64_status(self.fileId) & EX_IDS_INT64_API):
      id_map = (c_longlong * numObjs)()
    else:
      id_map = (c_int * numObjs)()
    EXODUS_LIB.ex_get_id_map(self.fileId,obj_type,byref(id_map))
    idMap = []
    for i in xrange(numObjs):
      idMap.append(id_map[i])
    return idMap

  # --------------------------------------------------------------------

  def __ex_put_id_map(self,objType,inqType,map):
    obj_type = c_int(objType)
    num_objs = c_int(self.__ex_inquire_int(inqType))
    numObjs = num_objs.value
    assert numObjs == len(map)
    id_map = []
    if (EXODUS_LIB.ex_int64_status(self.fileId) & EX_IDS_INT64_API):
      id_map = (c_longlong * numObjs)()
      for i in xrange(numObjs):
        id_map[i] = c_longlong( map[i] )
    else:
      id_map = (c_int * numObjs)()
      for i in xrange(numObjs):
        id_map[i] = c_int( map[i] )
    EXODUS_LIB.ex_put_id_map(self.fileId,obj_type,byref(id_map))
    return True

  # --------------------------------------------------------------------

  def __ex_get_elem_num_map(self):
    if (EXODUS_LIB.ex_int64_status(self.fileId) & EX_MAPS_INT64_API):
      elemNumMap = (c_longlong * self.numElem.value)()
    else:
      elemNumMap = (c_int * self.numElem.value)()
    EXODUS_LIB.ex_get_elem_num_map(self.fileId,byref(elemNumMap))
    return elemNumMap

  # --------------------------------------------------------------------

  def __ex_get_node_num_map(self):
    if (EXODUS_LIB.ex_int64_status(self.fileId) & EX_MAPS_INT64_API):
      nodeNumMap = (c_longlong * self.numNodes.value)()
    else:
      nodeNumMap = (c_int * self.numNodes.value)()
    EXODUS_LIB.ex_get_node_num_map(self.fileId,byref(nodeNumMap))
    return nodeNumMap

  # --------------------------------------------------------------------

  def __ex_get_elem_order_map(self):
    if (EXODUS_LIB.ex_int64_status(self.fileId) & EX_MAPS_INT64_API):
      elemOrderMap = (c_longlong * self.numElem.value)()
    else:
      elemOrderMap = (c_int * self.numElem.value)()
    EXODUS_LIB.ex_get_map(self.fileId,byref(elemOrderMap))
    return elemOrderMap

  # --------------------------------------------------------------------

  def __ex_get_elem_block(self,id):
    elem_block_id = c_longlong(id)
    elem_type = create_string_buffer(MAX_STR_LENGTH+1) 
    if (EXODUS_LIB.ex_int64_status(self.fileId) & EX_BULK_INT64_API):
      num_elem_this_blk  = c_longlong(0)
      num_nodes_per_elem = c_longlong(0)
      num_attr           = c_longlong(0)
    else:
      num_elem_this_blk  = c_int(0)
      num_nodes_per_elem = c_int(0)
      num_attr           = c_int(0)
    EXODUS_LIB.ex_get_elem_block(self.fileId,elem_block_id,elem_type,
                                 byref(num_elem_this_blk),byref(num_nodes_per_elem),\
                                 byref(num_attr))
    return(elem_type,num_elem_this_blk,num_nodes_per_elem,num_attr)

  # --------------------------------------------------------------------

  def __ex_put_elem_block(self,id,eType,numElems,numNodesPerElem,numAttrsPerElem):
    elem_block_id      = c_longlong(id)
    elem_type          = create_string_buffer(eType.upper(),MAX_STR_LENGTH+1)
    num_elem_this_blk  = c_longlong(numElems)
    num_nodes_per_elem = c_longlong(numNodesPerElem)
    num_attr           = c_longlong(numAttrsPerElem)
    EXODUS_LIB.ex_put_elem_block(self.fileId,elem_block_id,elem_type,
                                 num_elem_this_blk,num_nodes_per_elem,
                                 num_attr)

  # --------------------------------------------------------------------

  def __ex_get_elem_conn(self,id):
    (elem_type,num_elem_this_blk,num_nodes_per_elem,num_attr) = self.__ex_get_elem_block(id)
    elem_block_id = c_longlong(id)
    if (EXODUS_LIB.ex_int64_status(self.fileId) & EX_BULK_INT64_API):
      elem_block_connectivity = (c_longlong * (num_elem_this_blk.value * num_nodes_per_elem.value))()
    else:
      elem_block_connectivity = (c_int * (num_elem_this_blk.value * num_nodes_per_elem.value))()
    EXODUS_LIB.ex_get_elem_conn(self.fileId,elem_block_id,byref(elem_block_connectivity))
    return (elem_block_connectivity,num_elem_this_blk,num_nodes_per_elem)

  # --------------------------------------------------------------------

  def __ex_put_elem_conn(self,id,connectivity):
    (elem_type,num_elem_this_blk,num_nodes_per_elem,num_attr) = self.__ex_get_elem_block(id)
    elem_block_id = c_longlong(id)
    if (EXODUS_LIB.ex_int64_status(self.fileId) & EX_BULK_INT64_API):
      elem_block_connectivity = (c_longlong * (num_elem_this_blk.value * num_nodes_per_elem.value))()
      for i in range(num_elem_this_blk.value * num_nodes_per_elem.value):
        elem_block_connectivity[i] = c_longlong(connectivity[i])
    else:
      elem_block_connectivity = (c_int * (num_elem_this_blk.value * num_nodes_per_elem.value))()
      for i in range(num_elem_this_blk.value * num_nodes_per_elem.value):
        elem_block_connectivity[i] = c_int(connectivity[i])
    EXODUS_LIB.ex_put_elem_conn(self.fileId,elem_block_id,elem_block_connectivity)

  # --------------------------------------------------------------------

  def __ex_put_elem_attr(self,elemBlkID,Attr):
    elem_blk_id = c_longlong(elemBlkID)
    attrib = (c_double * len(Attr))()
    for i in range(len(Attr)):
      attrib[i] = c_double(Attr[i])
    EXODUS_LIB.ex_put_elem_attr(self.fileId,elem_blk_id,attrib)
  
  # --------------------------------------------------------------------

  def __ex_get_elem_attr(self,elemBlkID):
    elem_blk_id = c_longlong(elemBlkID)
    numAttrThisBlk = self.num_attr(elemBlkID)
    numElemsThisBlk = self.num_elems_in_blk(elemBlkID)
    totalAttr = numAttrThisBlk*numElemsThisBlk
    attrib = (c_double * totalAttr)()
    EXODUS_LIB.ex_get_elem_attr(self.fileId,elem_blk_id,byref(attrib))
    return attrib

  # --------------------------------------------------------------------

  def __ex_get_var_param(self,varChar):
    assert varChar.lower() in 'ngems'
    var_char = c_char(varChar)
    num_vars = c_int()
    EXODUS_LIB.ex_get_var_param(self.fileId,byref(var_char),byref(num_vars))
    return num_vars

  # --------------------------------------------------------------------

  def __ex_get_var_names(self,varChar):
    assert varChar.lower() in 'ngems'
    var_char = c_char(varChar)
    num_vars = self.__ex_get_var_param(varChar)
    var_name_ptrs = (POINTER(c_char * (MAX_STR_LENGTH+1)) * num_vars.value)()
    for i in range(num_vars.value):
      var_name_ptrs[i] = pointer(create_string_buffer(MAX_STR_LENGTH+1))
    EXODUS_LIB.ex_get_var_names(self.fileId,byref(var_char),num_vars,byref(var_name_ptrs))
    var_names = []
    for vnp in var_name_ptrs: var_names.append(vnp.contents.value)
    return var_names

  # --------------------------------------------------------------------

  def __ex_get_elem_var_tab(self):
    self.__ex_get_elem_blk_ids()
    num_blks = c_int(len(self.elemBlkIds))
    num_vars = self.__ex_get_var_param("e")
    truth_tab = (c_int * (num_blks.value * num_vars.value))()
    EXODUS_LIB.ex_get_elem_var_tab(self.fileId, num_blks, num_vars, byref(truth_tab))
    truthTab = []
    for val in truth_tab:
      if val:
        truthTab.append(True)
      else:
        truthTab.append(False)
    return truthTab

  # --------------------------------------------------------------------

  def __ex_put_elem_var_tab(self,truthTab):
    self.__ex_get_elem_blk_ids()
    num_blks = c_int(len(self.elemBlkIds))
    num_vars = self.__ex_get_var_param("e")
    truth_tab = (c_int * (num_blks.value*num_vars.value))()
    for i in xrange(len(truthTab)):
      boolVal = truthTab[i]
      if boolVal: 
        truth_tab[i] = c_int(1)
      else: 
        truth_tab[i] = c_int(0)
    EXODUS_LIB.ex_put_elem_var_tab(self.fileId,num_blks,num_vars,truth_tab)
    return True

  # --------------------------------------------------------------------

  def __ex_get_var(self,timeStep,varType,varId,blkId,numValues):
    step = c_int(timeStep)
    var_type = c_int(varType)
    var_id   = c_int(varId)
    block_id = c_longlong(blkId)
    num_values = c_longlong(numValues)
    var_vals = (c_double * num_values.value)()
    EXODUS_LIB.ex_get_var(self.fileId,step,var_type,var_id,block_id,num_values,var_vals)
    return var_vals

  # --------------------------------------------------------------------

  def __ex_put_var(self,timeStep,varType,varId,blkId,numValues,values):
    step = c_int(timeStep)
    var_type = c_int(varType)
    var_id   = c_int(varId)
    block_id = c_longlong(blkId)
    num_values = c_longlong(numValues)
    var_vals = (c_double * num_values.value)()
    for i in range(num_values.value):
      var_vals[i] = float(values[i])
    EXODUS_LIB.ex_put_var(self.fileId,step,var_type,var_id,block_id,num_values,var_vals)
    return True

  # --------------------------------------------------------------------

  def __ex_get_side_set_node_list_len(self,id):
    side_set_id = c_longlong(id)
    if (EXODUS_LIB.ex_int64_status(self.fileId) & EX_BULK_INT64_API):
      side_set_node_list_len = c_longlong(0)
    else:
      side_set_node_list_len = c_int(0)
    EXODUS_LIB.ex_get_side_set_node_list_len(self.fileId,side_set_id,byref(side_set_node_list_len))
    return side_set_node_list_len

  # --------------------------------------------------------------------

  def __ex_get_side_set_param(self,id):
    side_set_id = c_longlong(id)
    if (EXODUS_LIB.ex_int64_status(self.fileId) & EX_BULK_INT64_API):
      num_side_in_set      = c_longlong(0)
      num_dist_fact_in_set = c_longlong(0)
    else:
      num_side_in_set      = c_int(0)
      num_dist_fact_in_set = c_int(0)
    EXODUS_LIB.ex_get_side_set_param(self.fileId,side_set_id,byref(num_side_in_set),\
                                     byref(num_dist_fact_in_set))
    return (int(num_side_in_set.value),int(num_dist_fact_in_set.value))

  # --------------------------------------------------------------------

  def __ex_put_side_set_param(self,id,numSides,numDistFacts):
    side_set_id          = c_longlong(id)
    num_side_in_set      = c_longlong(numSides)
    num_dist_fact_in_set = c_longlong(numDistFacts)
    EXODUS_LIB.ex_put_side_set_param(self.fileId,side_set_id,num_side_in_set,num_dist_fact_in_set)
    return True

  # --------------------------------------------------------------------

  def __ex_get_side_set(self,sideSetId):
    side_set_id = c_longlong(sideSetId)
    (num_side_in_set,num_dist_fact_in_set) = self.__ex_get_side_set_param(sideSetId)
    if num_side_in_set == 0:
      return ([], [])
    if (EXODUS_LIB.ex_int64_status(self.fileId) & EX_BULK_INT64_API):
      side_set_elem_list = (c_longlong * num_side_in_set)()
      side_set_side_list = (c_longlong * num_side_in_set)()
    else:
      side_set_elem_list = (c_int * num_side_in_set)()
      side_set_side_list = (c_int * num_side_in_set)()
    EXODUS_LIB.ex_get_side_set(self.fileId,side_set_id,\
                               byref(side_set_elem_list),\
                               byref(side_set_side_list) )
    return (side_set_elem_list,side_set_side_list)

  # --------------------------------------------------------------------

  def __ex_put_side_set(self,id,sideSetElements,sideSetSides):
    side_set_id = c_longlong(id)
    if (EXODUS_LIB.ex_int64_status(self.fileId) & EX_BULK_INT64_API):
      side_set_elem_list = (c_longlong * len(sideSetElements))()
      side_set_side_list = (c_longlong * len(sideSetSides))()
      for i in range(len(sideSetElements)):
        side_set_elem_list[i] = c_longlong(sideSetElements[i])
        side_set_side_list[i] = c_longlong(sideSetSides[i])
    else:
      side_set_elem_list = (c_int * len(sideSetElements))()
      side_set_side_list = (c_int * len(sideSetSides))()
      for i in range(len(sideSetElements)):
        side_set_elem_list[i] = c_int(sideSetElements[i])
        side_set_side_list[i] = c_int(sideSetSides[i])
    EXODUS_LIB.ex_put_side_set(self.fileId,side_set_id,side_set_elem_list,side_set_side_list)
    return True

  # --------------------------------------------------------------------

  def __ex_get_sset_var_tab(self):
    self.__ex_get_side_set_ids()
    side_set_count = c_int(len(self.sideSetIds))
    variable_count = self.__ex_get_var_param("s")
    truth_table = (c_int * (side_set_count.value * variable_count.value))()
    EXODUS_LIB.ex_get_sset_var_tab(self.fileId,
                                   side_set_count,
                                   variable_count,
                                   byref(truth_table))
    truthTab = []
    for val in truth_table:
      if val:
        truthTab.append(True)
      else:
        truthTab.append(False)
    return truthTab

  # --------------------------------------------------------------------

  def __ex_put_sset_var_tab(self,truthTab):
    self.__ex_get_side_set_ids()
    num_blks = c_int(len(self.sideSetIds))
    num_vars = self.__ex_get_var_param("s")
    truth_tab = (c_int * (num_blks.value*num_vars.value))()
    for i in xrange(len(truthTab)):
      boolVal = truthTab[i]
      if boolVal: 
        truth_tab[i] = c_int(1)
      else: 
        truth_tab[i] = c_int(0)
    EXODUS_LIB.ex_put_sset_var_tab(self.fileId,num_blks,num_vars,truth_tab)
    return True

  # --------------------------------------------------------------------

  def __ex_get_side_set_dist_fact(self, sideSetId):
    side_set_id = c_longlong(sideSetId)
    side_set_node_list_len = self.__ex_get_side_set_node_list_len(sideSetId)
    set_dfs = (c_double * side_set_node_list_len.value)()
    EXODUS_LIB.ex_get_side_set_dist_fact(self.fileId,side_set_id,byref(set_dfs))
    return set_dfs

  # --------------------------------------------------------------------

  def __ex_put_side_set_dist_fact(self,sideSetId,sideSetDistFact):
    side_set_id = c_longlong(sideSetId)
    side_set_dist_fact = (c_double * len(sideSetDistFact))()
    for i in range(len(sideSetDistFact)):
      side_set_dist_fact[i] = c_double(sideSetDistFact[i])
    EXODUS_LIB.ex_put_side_set_dist_fact(self.fileId,side_set_id,side_set_dist_fact)

  # --------------------------------------------------------------------

  def __ex_get_side_set_node_list(self,id):
    side_set_id = c_longlong(id)
    side_set_node_list_len = self.__ex_get_side_set_node_list_len(id)
    (num_side_in_set,num_dist_fact_in_set) = self.__ex_get_side_set_param(id)
    if (EXODUS_LIB.ex_int64_status(self.fileId) & EX_BULK_INT64_API):
      side_set_node_cnt_list = (c_longlong * num_side_in_set)()
      side_set_node_list     = (c_longlong * side_set_node_list_len.value)()
    else:
      side_set_node_cnt_list = (c_int * num_side_in_set)()
      side_set_node_list     = (c_int * side_set_node_list_len.value)()
    EXODUS_LIB.ex_get_side_set_node_list(self.fileId,side_set_id,\
                                         byref(side_set_node_cnt_list),\
                                         byref(side_set_node_list) )
    return (side_set_node_cnt_list,side_set_node_list)

  # --------------------------------------------------------------------

  def __ex_get_sset_var(self,timeStep,varId,id):
    step        = c_int(timeStep)
    var_id      = c_int(varId)
    side_set_id = c_longlong(id)
    (numSideInSet,numDistFactInSet) = self.__ex_get_side_set_param(id)
    ss_var_vals = (c_double * numSideInSet)()
    num_side_in_set = c_longlong(numSideInSet)
    EXODUS_LIB.ex_get_sset_var(self.fileId,step,var_id,side_set_id,num_side_in_set,ss_var_vals)
    return ss_var_vals

  # --------------------------------------------------------------------

  def __ex_put_sset_var(self,timeStep,varId,id,values):
    step        = c_int(timeStep)
    var_id      = c_int(varId)
    side_set_id = c_longlong(id)
    (numSideInSet,numDistFactInSet) = self.__ex_get_side_set_param(id)
    num_side_in_set = c_longlong(numSideInSet)
    ss_var_vals = (c_double * numSideInSet)()
    for i in range(numSideInSet):
      ss_var_vals[i] = float(values[i])
    EXODUS_LIB.ex_put_sset_var(self.fileId,step,var_id,side_set_id,num_side_in_set,ss_var_vals)
    return True

  # --------------------------------------------------------------------

  def __ex_get_variable_param(self,varType):
    var_type = c_int(varType)
    num_vars = c_int(0)
    EXODUS_LIB.ex_get_variable_param(self.fileId,var_type,byref(num_vars))
    return num_vars

  # --------------------------------------------------------------------

  def __ex_put_variable_param(self,varType,numVars):
    var_type = c_int(varType)
    num_vars = c_int(numVars)
    current_num = self.__ex_get_variable_param(varType)
    if current_num.value == num_vars.value:
      ##print "value already set"
      return True
    errorInt = EXODUS_LIB.ex_put_variable_param(self.fileId,var_type,num_vars)
    if errorInt != 0:
      print "ERROR code =", errorInt
      raise Exception, ("ERROR: ex_put_variable_param had problems. This can only be called once per varType.")
    return True

  # --------------------------------------------------------------------

  def __ex_get_variable_name(self,varType,varId):
    var_type = c_int(varType)
    var_id   = c_int(varId)
    name = create_string_buffer(MAX_STR_LENGTH+1)
    EXODUS_LIB.ex_get_variable_name(self.fileId,var_type,var_id,name)
    return name

  # --------------------------------------------------------------------

  def __ex_put_variable_name(self,varType,varId,varName):
    var_type = c_int(varType)
    var_id   = c_int(varId)
    name = create_string_buffer(varName,MAX_STR_LENGTH+1)
    EXODUS_LIB.ex_put_variable_name(self.fileId,var_type,var_id,name)
    return True

  # --------------------------------------------------------------------

  def __ex_get_elem_attr_names(self,blkId):
    object_id = c_int(blkId)
    num_attr = c_int(self.num_attr(blkId))
    len_name = self.__ex_inquire_int(ex_inquiry("EX_INQ_READ_NAME_LENGTH"))
    attr_name_ptrs = (POINTER(c_char * (len_name+1)) * num_attr.value)() 
    for i in range(num_attr.value): 
      attr_name_ptrs[i] = pointer(create_string_buffer(len_name+1)) 
    EXODUS_LIB.ex_get_elem_attr_names(self.fileId,object_id,byref(attr_name_ptrs)) 
    attr_names = []
    for cnp in attr_name_ptrs: attr_names.append(cnp.contents.value)
    return attr_names

  # --------------------------------------------------------------------

  def __ex_put_elem_attr_names(self,blkId,varNames):
    object_id = c_int(blkId)
    num_attr = c_int(self.num_attr(blkId))
    len_name = self.__ex_inquire_int(ex_inquiry("EX_INQ_READ_NAME_LENGTH"))
    attr_name_ptrs = (POINTER(c_char * (len_name+1)) * num_attr.value)()
    assert len(varNames) == num_attr.value
    for i in range(num_attr.value):
      attr_name_ptrs[i] = pointer(create_string_buffer(varNames[i],len_name+1))
    EXODUS_LIB.ex_put_elem_attr_names(self.fileId,object_id,byref(attr_name_ptrs))
    return True

  # --------------------------------------------------------------------

  def __ex_get_prop_names(self,varType,inqType):
    var_type = c_int(varType)
    num_props = c_int(self.__ex_inquire_int(ex_inquiry(inqType)))
    prop_name_ptrs = (POINTER(c_char * (MAX_STR_LENGTH+1)) * num_props.value)() 
    for i in range(num_props.value): 
      prop_name_ptrs[i] = pointer(create_string_buffer(MAX_STR_LENGTH+1)) 
    EXODUS_LIB.ex_get_prop_names(self.fileId,var_type,byref(prop_name_ptrs)) 
    prop_names = []
    for cnp in prop_name_ptrs: prop_names.append(cnp.contents.value)
    return prop_names

  # --------------------------------------------------------------------

  def __ex_get_prop(self,objType,objId,propName):
    obj_type = c_int(objType)
    obj_id = c_longlong(objId)
    prop_name = create_string_buffer(propName,MAX_STR_LENGTH+1)
    if (EXODUS_LIB.ex_int64_status(self.fileId) & EX_IDS_INT64_API):
      prop_val = c_longlong(0)
    else:
      prop_val = c_int(0)
    EXODUS_LIB.ex_get_prop(self.fileId,obj_type,obj_id,byref(prop_name),byref(prop_val))
    return prop_val.value

  # --------------------------------------------------------------------

  def __ex_put_prop(self,objType,objId,propName,propVal):
    obj_type  = c_int(objType)
    obj_id    = c_longlong(objId)
    prop_name = create_string_buffer(propName,MAX_STR_LENGTH+1)
    prop_val  = c_longlong(propVal)
    EXODUS_LIB.ex_put_prop(self.fileId,obj_type,obj_id,byref(prop_name),prop_val)
    return True

  # --------------------------------------------------------------------

  def __ex_update(self):
    EXODUS_LIB.ex_update(self.fileId)
    return True

# --------------------------------------------------------------------
# Utility Functions
# --------------------------------------------------------------------

def collectElemConnectivity(exodusHandle,connectivity):
  """
    This function generates a list of lists that represent the element connectivity.

    Usage:

    exodusHandle = exodus("file.g","r")
    connectivity = []
    collectElemConnectivity(exodusHandle,connectivity)
    
    exodusHandle.close()
  """

  if type(connectivity) is not list:
    raise Exception, ("ERROR: connectivity is not a list in call to collectElemConnectivity().")
  if connectivity:
    raise Exception, ("ERROR: connectivity is not empty in call to collectElemConnectivity().")

  blockIds = exodusHandle.get_elem_blk_ids()
  for blId in blockIds:
    (elem_block_conn,num_elem,num_nodes) = exodusHandle.get_elem_connectivity(blId)
    for k in range(num_elem):
      i = k * num_nodes
      j = i + num_nodes
      local_elem_conn = elem_block_conn[i:j]
      connectivity.append( local_elem_conn )

# --------------------------------------------------------------------

def collectLocalNodeToLocalElems(exodusHandle,connectivity,localNodeToLocalElems):
  """
    This function generates a list of lists to go from local node id
    to local elem id.

    Usage:

    exodusHandle = exodus("file.g","r")
    connectivity = [] ## If this is not empty it will assume it is already filled.
    localNodeToLocalElems = []
    collectLocalNodeToLocalElems(exodusHandle,connectivity,localNodeToLocalElems)
    
    exodusHandle.close()
  """

  if type(connectivity) is not list:
    raise Exception, ("ERROR: connectivity is not a list in call to collectLocalNodeToLocalElems().")
  if type(localNodeToLocalElems) is not list:
    raise Exception, ("ERROR: localNodeToLocalElems is not a list in call to collectLocalNodeToLocalElems().")
  if localNodeToLocalElems:
    raise Exception, ("ERROR: localNodeToLocalElems is not empty in call to collectLocalNodeToLocalElems().")

  if not connectivity:
    collectElemConnectivity(exodusHandle,connectivity)

  numNodes = exodusHandle.num_nodes()
  for i in range(numNodes+1):
    localNodeToLocalElems.append([])
  localElemId = 0
  for local_elem_conn in connectivity:
    for n in local_elem_conn:
      localNodeToLocalElems[n].append(localElemId)
    localElemId = localElemId + 1

# --------------------------------------------------------------------

def collectLocalElemToLocalElems(exodusHandle,connectivity,localNodeToLocalElems,localElemToLocalElems):
  """
    This function generates a list of lists to go from local elem id
    to connected local elem ids.

    Usage:

    exodusHandle = exodus("file.g","r")
    connectivity = [] ## If this is not empty it will assume it is already filled.
    localNodeToLocalElems = [] ## If this is not empty it will assume it is already filled.
    localElemToLocalElems = []
    collectLocalElemToLocalElems(exodusHandle,connectivity,localNodeToLocalElems,localElemToLocalElems)
    
    exodusHandle.close()
  """

  if type(connectivity) is not list:
    raise Exception, ("ERROR: connectivity is not a list in call to collectLocalElemToLocalElems().")
  if type(localNodeToLocalElems) is not list:
    raise Exception, ("ERROR: localNodeToLocalElems is not a list in call to collectLocalElemToLocalElems().")
  if type(localElemToLocalElems) is not list:
    raise Exception, ("ERROR: localElemToLocalElems is not a list in call to collectLocalElemToLocalElems().")
  if localElemToLocalElems:
    raise Exception, ("ERROR: localElemToLocalElems is not empty in call to collectLocalElemToLocalElems().")

  if not connectivity:
    collectElemConnectivity(exodusHandle,connectivity)
  if not localNodeToLocalElems:
    collectLocalNodeToLocalElems(exodusHandle,connectivity,localNodeToLocalElems)

  numElems = exodusHandle.num_elems()
  for i in range(numElems):
    localElemToLocalElems.append([])
  for localElemId in range(numElems):
    nodeList = list(connectivity[localElemId])
    newConnectedElems = []
    for n in nodeList:
      for elem in localNodeToLocalElems[n]:
        newConnectedElems.append( elem )
    localElemToLocalElems[localElemId] = list( set(newConnectedElems) )

# --------------------------------------------------------------------

def copyTransfer(fromFileName,toFileName,additionalGlobalVariables=[],additionalNodalVariables=[],additionalElementVariables=[]):
  """
    This function creates an exodus file toFileName and copies everything from exodus file fromFileName returning
    a file handle to toFileName.

    Additional space is allocated for additionalGlobalVariables, additionalNodalVariables and additionalElementVariables if specified.

    additionalGlobalVariables: list of global variable names to add.

    additionalNodalVaraibles: list of nodal variable names to add.

    additionalElementVariables: should be a list of element variable names to add to all blocks or
                                tuples ( name, blkIds ) where name is the element variable to add
                                and blkIds is a list of blkIds to add it to.
    Usage:

    fromFileName = "input.e"
    toFileName = "output.e"
    addGlobalVariables = [] ## Do not add any new global variables
    addNodeVariables = ["node_dummy1","node_dummy2"] ## Add node_dummy1 and node_dummy2 as new node variables
    addElementVariables = [ ("elem_dummy1",[1,2,3]), "elem_dummy2" ] ## Add elem_dummy1 on blkIds 1,2,3 and elem_dummy2 on all blocks

    toFileHandle = copyTranfer(fromFileName,toFileName,addGlobalVariables,addNodeVariables,addElementVariables)

    ## Fill in new variables
    
    toFileHandle.close()

  """

  debugPrint = False

  if type(additionalGlobalVariables) is not list:
    raise Exception, ("ERROR: additionalGlobalVariables is not a list in call to copyTransfer().")
  if type(additionalNodalVariables) is not list:
    raise Exception, ("ERROR: additionalNodalVariables is not a list in call to copyTransfer().")
  if type(additionalElementVariables) is not list:
    raise Exception, ("ERROR: additionalElementVariables is not a list in call to copyTransfer().")
  if ( os.path.isfile(toFileName) ):
    raise Exception, ("ERROR: ", toFileName, " file already exists cowardly exiting instead of overwriting in call to copyTransfer().")
  
  exoFrom = exodus(fromFileName,"r")
  
  if debugPrint: print "Initialize New File"
  title = exoFrom.title()
  numDim = exoFrom.num_dimensions()
  numNodes = exoFrom.num_nodes()
  numElems = exoFrom.num_elems()
  numBlks = exoFrom.num_blks()
  numNodeSets = exoFrom.num_node_sets()
  numSideSets = exoFrom.num_side_sets()
  exoTo = exodus( toFileName, "w", title, numDim, numNodes, numElems, numBlks, numNodeSets, numSideSets )
  
  if debugPrint: print "Transfer QA records"
  qaRecords = exoFrom.get_qa_records()
  exoTo.put_qa_records( qaRecords )
  
  if debugPrint: print "Transfer Info records"
  numInfoRecs = exoFrom.num_info_records()
  if numInfoRecs > 0:
    infoRecs = exoFrom.get_info_records()
    exoTo.put_info_records( infoRecs )
  
  if debugPrint: print "Transfer Nodal Coordinates and Names"
  exoTo.put_coord_names( exoFrom.get_coord_names() )
  (xCoords,yCoords,zCoords) = exoFrom.get_coords()
  exoTo.put_coords(xCoords,yCoords,zCoords)
  
  if debugPrint: print "Transfer Node Id Map"
  nodeIdMap = exoFrom.get_node_id_map()
  exoTo.put_node_id_map(nodeIdMap)
  
  if debugPrint: print "Transfer Element Data"
  blkIds = exoFrom.get_elem_blk_ids()
  for blkId in blkIds:
    (elemType,numElem,nodesPerElem,numAttr) = exoFrom.elem_blk_info(blkId)
    exoTo.put_elem_blk_info(blkId,elemType,numElem,nodesPerElem,numAttr)
    (connectivity,numElem,nodesPerElem) = exoFrom.get_elem_connectivity(blkId)
    exoTo.put_elem_connectivity(blkId,connectivity)
    if numAttr > 0:
      attrNames = exoFrom.get_element_attribute_names(blkId)
      exoTo.put_element_attribute_names(blkId,attrNames)
      exoTo.put_elem_attr(blkId, exoFrom.get_elem_attr(blkId))
    elemProps = exoFrom.get_element_property_names()
    for elemProp in elemProps:
      propVal = exoFrom.get_element_property_value(blkId,elemProp)
      if elemProp == "ID" and propVal == blkId:
        continue
      else:
        exoTo.put_element_property_value(blkId,elemProp,propVal)
    blockName = exoFrom.get_elem_blk_name(blkId)
    exoTo.put_elem_blk_name(blkId,blockName)
  
  if debugPrint: print "Transfer Element Id Map"
  elemIdMap = exoFrom.get_elem_id_map()
  exoTo.put_elem_id_map(elemIdMap)
  
  if debugPrint: print "Transfer Node Sets"
  if numNodeSets > 0:
    nodeSetIds = exoFrom.get_node_set_ids()
    for nsId in nodeSetIds:
      (numSetNodes,numSetDistFacts) = exoFrom.get_node_set_params(nsId)
      exoTo.put_node_set_params(nsId,numSetNodes,numSetDistFacts)
      nsNodes = exoFrom.get_node_set_nodes(nsId)
      exoTo.put_node_set(nsId,nsNodes)
      if numSetDistFacts > 0:
        nsDF = exoFrom.get_node_set_dist_facts(nsId)
        exoTo.put_node_set_dist_fact(nsId,nsDF)
      nodeSetName = exoFrom.get_node_set_name(nsId)
      exoTo.put_node_set_name(nsId,nodeSetName)
    nodeSetProps = exoFrom.get_node_set_property_names()
    for nodeSetProp in nodeSetProps:
      propVal = exoFrom.get_node_set_property_value(nsId,nodeSetProp)
      if nodeSetProp == "ID" and propVal == nsId:
        continue
      else:
        exoTo.put_node_set_property_value(nsId,nodeSetProp,propVal)
  
  if debugPrint: print "Transfer Side Sets"
  if numSideSets > 0:
    sideSetIds = exoFrom.get_side_set_ids()
    for ssId in sideSetIds:
      (numSetSides,numSetDistFacts) = exoFrom.get_side_set_params(ssId)
      exoTo.put_side_set_params(ssId,numSetSides,numSetDistFacts)
      (elemList,sideList) = exoFrom.get_side_set(ssId)
      exoTo.put_side_set(ssId,elemList,sideList)
      if numSetDistFacts > 0:
        ssDF = exoFrom.get_side_set_dist_fact(ssId)
        exoTo.put_side_set_dist_fact(ssId,ssDF)
      sideSetName = exoFrom.get_side_set_name(ssId)
      exoTo.put_side_set_name(ssId,sideSetName)
    sideSetProps = exoFrom.get_side_set_property_names()
    for sideSetProp in sideSetProps:
      propVal = exoFrom.get_side_set_property_value(ssId,sideSetProp)
      if sideSetProp == "ID" and propVal == ssId:
        continue
      else:
        exoTo.put_side_set_property_value(ssId,sideSetProp,propVal)
  
  if debugPrint: print "Transfer time values"
  nSteps = exoFrom.num_times()
  timeVals = exoFrom.get_times()
  for step in xrange(nSteps):
    exoTo.put_time( step+1, timeVals[step] )
  
  if debugPrint: print "Transfer Global Variables"
  nNewGlobalVars = len(additionalGlobalVariables)
  nGlobalVars = exoFrom.get_global_variable_number() + nNewGlobalVars
  defaultNewVarVals = []
  for i in xrange(nNewGlobalVars):
    defaultNewVarVals.append(0.0)
  if nGlobalVars > 0:
    exoTo.set_global_variable_number(nGlobalVars)
    gVarNames = exoFrom.get_global_variable_names()
    gVarNames.extend( additionalGlobalVariables )
    for nameIndex in xrange(nGlobalVars):
      globalVarName = gVarNames[nameIndex]
      exoTo.put_global_variable_name(globalVarName,nameIndex+1)
    for step in xrange(nSteps):
      gValues = exoFrom.get_all_global_variable_values(step+1)
      gValues.extend( defaultNewVarVals )
      exoTo.put_all_global_variable_values(step+1,gValues)
  
  if debugPrint: print "Transfer Nodal Variables"
  nNewNodalVars = len(additionalNodalVariables)
  nOrigNodalVars = exoFrom.get_node_variable_number()
  nNodalVars = nOrigNodalVars + nNewNodalVars
  if nNodalVars > 0:
    exoTo.set_node_variable_number(nNodalVars)
    nVarNames = exoFrom.get_node_variable_names()
    nVarNames.extend( additionalNodalVariables )
    for nameIndex in xrange(nNodalVars):
      nodalVarName = nVarNames[nameIndex]
      exoTo.put_node_variable_name(nodalVarName,nameIndex+1)
      if nameIndex < nOrigNodalVars:
        for step in xrange(nSteps):
          nValues = exoFrom.get_node_variable_values(nodalVarName,step+1)
          exoTo.put_node_variable_values(nodalVarName,step+1,nValues)
  
  if debugPrint: print "Construct Truth Table for additionalElementVariables"
  newElemVariableNames = []
  newElemVariableBlocks = []
  for item in additionalElementVariables:
    if type(item) is tuple:
      newElemVariableNames.append( item[0] )
      inBlks = []
      for blkId in item[1]:
        if blkId in blkIds:
          inBlks.append(blkId)
      newElemVariableBlocks.append( inBlks )
    elif type(item) is str:
      newElemVariableNames.append( item )
      newElemVariableBlocks.append( blkIds )
    else:
      print "Warning additionalElementVariable item ", item, " is not right type to add."
      print "should be a string or tuple, skipping"
  
  if debugPrint: print "Transfer Element Variables"
  nNewElemVars = len(newElemVariableNames)
  nOrigElemVars = exoFrom.get_element_variable_number()
  nElemVars = nOrigElemVars + nNewElemVars
  if nElemVars > 0:
    exoTo.set_element_variable_number(nElemVars)
    origElemVarNames = exoFrom.get_element_variable_names()
    eVarNames =  exoFrom.get_element_variable_names()
    eVarNames.extend( newElemVariableNames )
    truthTable = []
    if nOrigElemVars > 0:
      truthTable = exoFrom.get_element_variable_truth_table()
    if nNewElemVars > 0:
      newTruth = []
      for j in xrange(numBlks):
        for k in xrange(nOrigElemVars):
          index = j*nOrigElemVars + k
          newTruth.append( truthTable[index] )
        for m in xrange(nNewElemVars):
          if blkIds[j] in newElemVariableBlocks[m]:
            newTruth.append(True)
          else:
            newTruth.append(False)
      truthTable = newTruth
    exoTo.set_element_variable_truth_table(truthTable)
    for nameIndex in xrange(nElemVars):
      elemVarName = eVarNames[nameIndex]
      exoTo.put_element_variable_name(elemVarName,nameIndex+1)
    truthIndex = 0
    for blkId in blkIds:
      for eVarName in origElemVarNames:
        if truthTable[truthIndex]:
          for step in xrange(nSteps):
            eValues = exoFrom.get_element_variable_values(blkId,eVarName,step+1)
            exoTo.put_element_variable_values(blkId,eVarName,step+1,eValues)
        truthIndex = truthIndex + 1
      truthIndex = truthIndex + nNewElemVars
  
  ## TODO: Transfer Nodeset Variables
  
  ## TODO: Transfer Sideset Variables
  
  exoFrom.close()
  return exoTo

