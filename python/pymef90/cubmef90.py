###
### geometry tools
###
def WebcutTool(body_list_in, tool_ID, delete=False):
  import cubit
  ### delete group 'webcut_group' if it exists
  cubit.delete_group(cubit.get_id_from_name('webcut_group'))
  ### webcut 
  cmd='webcut volume '
  for i in body_list_in:
    cmd += '%i ' % i
  cmd += 'tool volume %i group_results' % tool_ID
  cubit.cmd(cmd)
  ### identify the grain and the outside
  ### I need a list of the volume IDs in webcut_group, but cubit.get_group_volumes returns a t-uple
  tmp_IDs = cubit.get_id_from_name('webcut_group')
  webcut_group_IDs = [vol_id for vol_id in cubit.get_group_volumes(tmp_IDs) ]
  webcut_group_IDs.append(tool_ID)
  ovlp = cubit.get_overlapping_volumes(webcut_group_IDs) 
  ### ovlp now contains the overlap and the tool
  ### remove tool_id to get tool_list_out
  tmp_list = [vol_id for vol_id in ovlp if not vol_id == tool_ID]
  tool_list_out = tmp_list # just so that I can make a call with tool_ID = tool_list_out
  ### remove tool_list_out and tool_ID from ovlp
  body_list_out = [vol_id for vol_id in webcut_group_IDs if vol_id not in ovlp]
  ### delete the tool of requested
  if delete:
    cubit.cmd('delete volume %i' % tool_ID)
  ### Be nice and clean up
  cubit.delete_group(cubit.get_id_from_name('webcut_group'))
  return body_list_out, tool_list_out

def GrainCreate(vertices, LZ=10):
  import cubit
  ### Create a grain from a vertex list
  ###
  ### Create the vertices
  vert_ID=[]
  for X in vertices:
    cubit.cmd("create vertex X %f Y %f Z %f" % (X[0], X[1], -LZ/2.0))
    vert_ID.append(cubit.get_last_id("vertex"))
  ###
  ### Create the edges
  curv_ID=[]
  for i in range(len(vert_ID)-1):
    cubit.cmd ("create curve vertex %i vertex %i" % (vert_ID[i], vert_ID[i+1]))
    curv_ID.append(cubit.get_last_id("curve"))
  cubit.cmd ("create curve vertex %i vertex %i" % (vert_ID[-1], vert_ID[0]))
  curv_ID.append(cubit.get_last_id("curve"))
  ###
  ### Create the surface
  cmd="create surface curve"
  for edge in curv_ID:
    cmd+=" %i" % edge
  cubit.cmd(cmd)
  surf_ID=cubit.get_last_id("surface")
  cubit.cmd("sweep surface %i direction Z distance %f" % (surf_ID, LZ))
  grain_ID=cubit.get_last_id("volume")
  return grain_ID

def IsInBox(X, box):
### Check if x is in the quad / hex box
  result=True
  for j in range(len(X)):
    if not box[2*j] <= X[j] <= box[2*j+1]:
      result=False
  return result


def SitesGenRect2D(N, rect, BB):
### box is given in the form xmin, xmax, ymin, ymax
### BB is given in the form xmin. xmax, ymin, ymax
  import numpy as np

  sites=np.empty([N+4,2])
  sites[:N,0]  = np.random.uniform(rect[0], rect[1], N)
  sites[:N,1]  = np.random.uniform(rect[2], rect[3], N)
  sites[N,:]   = np.array([BB[0], BB[2]])
  sites[N+1,:] = np.array([BB[1], BB[2]])
  sites[N+2,:] = np.array([BB[0], BB[3]])
  sites[N+3,:] = np.array([BB[1], BB[3]])
  return sites

def SitesGenCircle2D(N, circle, BB):
### circle is given in the form xc, yc, R
### BB is given in the form xmin, xmax, ymin, ymax
  import numpy as np

  sites=np.empty([N+4,2])
  R = circle[0]*(np.random.power(2, size=N))
  Theta = np.random.uniform(0, 2*np.pi, N)
  sites[:N,0]  = circle[1]+R*np.cos(Theta)
  sites[:N,1]  = circle[2]+R*np.sin(Theta)
  sites[N,:]   = np.array([BB[0], BB[2]])
  sites[N+1,:] = np.array([BB[1], BB[2]])
  sites[N+2,:] = np.array([BB[0], BB[3]])
  sites[N+3,:] = np.array([BB[1], BB[3]])
  return sites

def SitesRead(filename):
  import numpy as np

  inputfile=open(filename, 'r')
  ### read file header
  buffer=inputfile.readline()
  numsites=int(buffer.split()[3])
  numdim=int(buffer.split()[0])

  ### Throw away line 2
  buffer=inputfile.readline()

  coords=np.empty([numsites, numdim])
  ### Read coordinates

  for i in range(numsites):
    buffer=inputfile.readline()
    coords[i,:] = buffer.split()
  inputfile.close()
  return coords

def SitesWrite(filename, sitecoords):
  import numpy as np

  sitefile=open(filename, 'w')
  sitefile.write('%i  rbox D%i %i\n' % (sitecoords.shape[0], sitecoords.shape[0], sitecoords.shape[1]))
  sitefile.write('%i\n' % sitecoords.shape[0])
  for site in sitecoords:
    for coord in site:
      sitefile.write('%f   ' % coord)
    sitefile.write('\n')
  sitefile.close()

def SitesToStr(sitecoords):
  import numpy as np

  sitestr='%i  rbox D%i %i\n%i\n' % (sitecoords.shape[0], sitecoords.shape[0], sitecoords.shape[1], sitecoords.shape[1])
  for site in sitecoords:
    for coord in site:
      sitestr += '%f   ' % coord
    sitestr += '\n '
  return(sitestr)

def SitesPlot(sites, pattern=None):
  import matplotlib.pyplot as plt
  import numpy as np

  if not pattern:
    pattern= 'ro'
  for site in sites:
    plt.plot(site[0], site[1], pattern)

def VoronoiRead():
  verbose = False

  inputfile=open(filename, 'r')
  ### read file header
  f = open(inputfile, mode = 'r')
  buffer=f.readline()
  #if dim!= 2:
  #  parser.error(option.inputfile + 'is not a 2d voronoi file')
  buffer=f.readline()
  numvertices = int(buffer.split()[0])
  numcells    = int(buffer.split()[1])

  if verbose:
    print( "number of vertices: %i \nnumber of cells:    %i\n" % (numvertices, numcells))

  #### Read vertices
  vertices=np.empty([numvertices,2])
  for i in range(numvertices):
    buffer = f.readline()
    vertices[i,:]= buffer.split()

  if verbose:
    print("vertices:\n")
    print(vertices)

  ### Read cells, and filter through
  cells=[]
  for i in range(numcells):
    buffer = f.readline()
    tmpcell = buffer.split()
    is_infinity=False
    for j in tmpcell:
      if j=='0':
        is_infinity=True
    if verbose:
      print(i, tmpcell, is_infinity)
    if not is_infinity:
      cells.append(tmpcell)

  numcells=len(cells)
  if verbose:
    print('number of interior cells: %i\n' % numcells)
    print(cells)
  inputfile.close()

  return vertices, cells

def VoronoiPlot(vertices, cells, pattern=None, legend=False):
  import matplotlib
  import numpy as np
  from matplotlib.patches import Polygon
  from matplotlib.collections import PatchCollection
  import pylab

  ax=matplotlib.pyplot.subplot(111, aspect='equal')

###$  ### draw vertices
###$  if not pattern:
###$    pattern='b*'
###$  for vertex in vertices[1:]:
###$    matplotlib.pyplot.plot(vertex[0], vertex[1], pattern)
###$  
  ### draw cells
  patches = []
  for cell in cells:
    patches.append(Polygon(vertices[cell], True))
  colors = np.random.rand(len(patches))
  p = PatchCollection(patches, cmap=matplotlib.cm.hsv, alpha=0.4)
  p.set_array(pylab.array(colors))
  ax.add_collection(p)

  if legend:
    for i in range(len(cells)):
      cell = cells[i]
      X = sum(vertices[cell,0]) / len(cell)
      Y = sum(vertices[cell,1]) / len(cell)
      matplotlib.pyplot.text(X, Y, str(i), size=12)

def VoronoiGen(sites):
  import subprocess
  import numpy as np

  verbose = False
  qvoronoi = subprocess.Popen('qvoronoi o',
                          shell=True,
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          )
  qvoronoi.stdin.write(SitesToStr(sites))
  voronoistr=qvoronoi.communicate()[0].split('\n')

  if verbose:
    print(voronoistr)

  numdim=int(voronoistr[0])
  (numvertices, numcells, junk) = [int(s) for s in voronoistr[1].split()]
  if verbose:
    print('Number of vertices: \t%i\nNumber of cells: \t%i\n' % (numvertices, numcells))

  vertices=np.empty([numvertices, numdim])
  for i in range(numvertices):
    vertices[i] = [float(s) for s in voronoistr[i+2].split()]

  if verbose:
    print("vertices:\n")
    print(vertices)

  ### Read cells, and filter through
  cells=[]
  for i in range(numcells):
    if not int(voronoistr[numvertices+2+i].split()[0]) == 0:
      cells.append([int(s) for s in voronoistr[numvertices+2+i].split()[1:]])

  return vertices, cells

def VoronoiRemoveInfiniteCells(cells):
### Read cells, and filter through
  realcells=[]
  for cell in cells:
    if 0 not in cell:
      realcells.append(cell)
  return realcells

def SitesSmoother(sites, vertices, cells):
### replace sites with center of their voronoi cells
### leave sites associated with infinite cells unmoved
  import numpy as np

  newsites=np.empty(sites.shape)
  newsites[:]=sites

  for i in range(sites.shape[0]):
    cell = cells[i]
    if 0 not in cell:
      for j in range(sites.shape[1]):
        newsites[i,j] = sum(vertices[cell,j]) / len(cell)
  return newsites

def SitesSmootherBox(sites, vertices, cells, box):
### replace sites with center of their voronoi cells
### leave sites associated with infinite cells unmoved
### leave sites outside of the box unchanged
### does not move sites outside of the box
  import numpy as np

  newsites=np.empty(sites.shape)
  newsites[:]=sites

  print("number of sites: %i vertices: %i cells: %i" % (len(sites), len(vertices), len(cells)) )
  for i in range(sites.shape[0]):
    cell = cells[i]
    if (0 not in cell) and IsInBox(sites[i,:], box):
      for j in range(sites.shape[1]):
        newsites[i,j] = min(max(sum(vertices[cell,j]) / len(cell), box[2*j]), box[2*j+1])
  return newsites

def SitesSmootherSphere(sites, vertices, cells, sphere):
  ### replace sites with center of their voronoi cells
  ### leave sites associated with infinite cells unmoved
  ### leave sites outside of the box unchanged
  ### does not move sites outside of the box
    import numpy as np

    newsites=np.empty(sites.shape)
    newsites[:]=sites

    for i in range(sites.shape[0]):
      cell = cells[i]
      if 0 not in cell:
        R=np.sqrt(sum(pow(sites[i,:]-sphere[1:],2)))
        if R<=sphere[0]:
          for j in range(sites.shape[1]):
            newsites[i,j] = sum(vertices[cell,j]) / len(cell)
          R=np.sqrt(sum(pow(sites[i,:]-sphere[1:],2)))
          if R>=sphere[0]:
            newsites[i,:] = sphere[1:] + (newsites[i,:]-sphere[1:])/R * sphere[0]
    return newsites

def BoxDraw(box):
  import matplotlib
  from matplotlib.patches import Polygon
  ax=matplotlib.pyplot.subplot(111, aspect='equal')
  rect=Polygon([ [box[0], box[2]], [box[1], box[2]], [box[1], box[3]], [box[0], box[3]] ])
  rect.set_fill(False)
  rect.set_linewidth(2)
  ax.add_patch(rect)

def CircleDraw(circle):
  import matplotlib
  from matplotlib.patches import Circle

  ax=matplotlib.pyplot.subplot(111, aspect='equal')
  circ = Circle((circle[1], circle[2]), circle[0])
  circ.set_fill(False)
  circ.set_linewidth(2)
  ax.add_patch(circ)

###
### Microstructure constructors
###
def Layer(Body_IDs, BB, Alpha, Theta1, Theta2, Xoffset=.5):
  import cubit
  import numpy as np
  ###
  ###    Body_ID: a single volume ID
  ###    BB:      bounding box in the form [xmin, xmax, ymin, ymax, zmin, zmax]
  ###    Alpha:   lamination angle
  ###    Theta1:  thickness of layer 1
  ###    Theta2:  thickness of layer 2
  ###    Xoffset: distance along the x-axis between the center of the bounding box
  ###             and the first interface l1-l2 if specified, in fraction of Theta1
  ###             i.e. 0 -> interface l2-l1 at center
  ###                  1 -> interface l1-l2 at center
  ###             otherwise, the center of BB corresponds to the center of a Layer1
  ###             i.e. Xoffset = .5
  ###                                            
  YL = np.sqrt((BB[1]-BB[0])**2+(BB[3]-BB[2])**2)
  XC = (BB[0] + BB[1]) / 2.
  YC = (BB[2] + BB[3]) / 2.
  ZC = (BB[4] + BB[5]) / 2.
  Theta = Theta1 + Theta2
  if Alpha>0:
    l1 = (BB[0]-XC) * np.cos(np.radians(Alpha)) + (BB[2]-YC) * np.sin(np.radians(Alpha))
    l2 = (BB[1]-XC) * np.cos(np.radians(Alpha)) + (BB[3]-YC) * np.sin(np.radians(Alpha))
  else:
    l1 = (BB[0]-XC) * np.cos(np.radians(Alpha)) + (BB[3]-YC) * np.sin(np.radians(Alpha))
    l2 = (BB[1]-XC) * np.cos(np.radians(Alpha)) + (BB[2]-YC) * np.sin(np.radians(Alpha))    
  n1 = int(np.ceil(l1/Theta))
  n2 = int(np.ceil(l2/Theta))
  LAYER1_IDs = Body_IDs[:]
  LAYER2_IDs = []
  ###
  ### offset layers
  ###
  for i in range(n1,n2):
    X = XC + i * np.cos( np.radians(Alpha) ) * Theta + np.cos( np.radians(Alpha) ) * Theta1 * (Xoffset - .5)
    Y = YC + i * np.sin( np.radians(Alpha) ) * Theta
    cubit.cmd('create brick X %f Y %f Z %f' % (Theta1, YL, BB[5]-BB[4]))
    tmp_ID=cubit.get_last_id("volume")
    cubit.cmd('move volume %i X %f Y %f Z %f' % (tmp_ID, X, Y, ZC))
    if not Alpha == 0.:
      ### cubit doesn't like rotations will 0 angle...
      cubit.cmd('create vertex X %f Y %f Z %f' % (X, Y, BB[4]))
      v1_ID=cubit.get_last_id("vertex")
      cubit.cmd('create vertex X %f Y %f Z %f' % (X, Y, BB[5]))
      v2_ID=cubit.get_last_id("vertex")
      cubit.cmd('rotate volume %i about vertex %i vertex %i angle %f' % (tmp_ID, v1_ID, v2_ID, Alpha))
      cubit.cmd('delete vertex %i' % v1_ID)
      cubit.cmd('delete vertex %i' % v2_ID)
    (LAYER1_IDs, tmp_layer) = WebcutTool(LAYER1_IDs, tmp_ID, delete=True)
    for l in tmp_layer:
      LAYER2_IDs.append(l)
  return (LAYER1_IDs, LAYER2_IDs)

def MilledLayer(Body_IDs, BB, Alpha, Theta1, Theta2, secmin, Xoffset=.5):
  import cubit
  import numpy as np
  ###
  ###    Body_ID: a single volume ID
  ###    BB:      bounding box in the form [xmin, xmax, ymin, ymax, zmin, zmax]
  ###    Alpha:   lamination angle
  ###    Theta1:  thickness of layer 1
  ###    Theta2:  thickness of layer 2
  ###    secmin:  the cross section after machining
  ###    Xoffset: distance along the x-axis between the center of the bounding box
  ###             and the first interface l1-l2 if specified, in fraction of Theta1
  ###             i.e. 0 -> interface l2-l1 at center
  ###                  1 -> interface l1-l2 at center
  ###             otherwise, the center of BB corresponds to the center of a Layer1
  ###             i.e. Xoffset = .5
  ###                                            
  YL = np.sqrt((BB[1]-BB[0])**2+(BB[3]-BB[2])**2)
  ZL = BB[5] - BB[4]
  XC = (BB[0] + BB[1]) / 2.
  YC = (BB[2] + BB[3]) / 2.
  ZC = (BB[4] + BB[5]) / 2.
  Theta = Theta1 + Theta2
  if Alpha>0:
    l1 = (BB[0]-XC) * np.cos(np.radians(Alpha)) + (BB[2]-YC) * np.sin(np.radians(Alpha))
    l2 = (BB[1]-XC) * np.cos(np.radians(Alpha)) + (BB[3]-YC) * np.sin(np.radians(Alpha))
  else:
    l1 = (BB[0]-XC) * np.cos(np.radians(Alpha)) + (BB[3]-YC) * np.sin(np.radians(Alpha))
    l2 = (BB[1]-XC) * np.cos(np.radians(Alpha)) + (BB[2]-YC) * np.sin(np.radians(Alpha))    
  n1 = int(np.ceil(l1/Theta))
  n2 = int(np.ceil(l2/Theta))
  for i in range(n1,n2):
    X = XC + i * np.cos( np.radians(Alpha) ) * Theta + np.cos( np.radians(Alpha) ) * Theta1 * (Xoffset - .5)
    Y = YC + i * np.sin( np.radians(Alpha) ) * Theta
    cubit.cmd('create brick X %f Y %f Z %f' % (Theta1, YL, ZL))
    tmp_ID=cubit.get_last_id("volume")
    cubit.cmd('move volume %i X %f Y %f Z %f' % (tmp_ID, X, Y, ZL/2.+secmin[1]))
    if not Alpha == 0.:
      ### cubit doesn't like rotations will 0 angle...
      cubit.cmd('create vertex X %f Y %f Z %f' % (X, Y, BB[4]))
      v1_ID=cubit.get_last_id("vertex")
      cubit.cmd('create vertex X %f Y %f Z %f' % (X, Y, BB[5]))
      v2_ID=cubit.get_last_id("vertex")
      cubit.cmd('rotate volume %i about vertex %i vertex %i angle %f' % (tmp_ID, v1_ID, v2_ID, Alpha))
      cubit.cmd('delete vertex %i' % v1_ID)
      cubit.cmd('delete vertex %i' % v2_ID)
    (Body_IDs, tmp_layer) = WebcutTool(Body_IDs, tmp_ID, delete=False)
    cmd="delete volume "
    for v in tmp_layer:
      cmd+=" %i" % v
    cubit.cmd(cmd)
    cubit.cmd("move volume %i Z %f" % (tmp_ID, -ZL-secmin[1]+secmin[0]))
    (Body_IDs, tmp_layer) = WebcutTool(Body_IDs, tmp_ID, delete=True)
    cmd="delete volume "
    for v in tmp_layer:
      cmd+=" %i" % v
    cubit.cmd(cmd)
  return (Body_IDs)

### 
### Various geometry constructors
###

def DCBCreate(LX, LY1, LY, LZ, EPSCRACK, LCRACK):
  import cubit
  import numpy as np
  ###
  vertices=np.empty([7,3])
  vertices[0,:] = [LCRACK, 0, -LZ/2.]
  vertices[1,:] = [0, EPSCRACK/2., -LZ/2.]
  vertices[2,:] = [0, LY1/2., -LZ/2.]
  vertices[3,:] = [LX, LY/2., -LZ/2.]
  vertices[4,:] = [LX, -LY/2., -LZ/2.]
  vertices[5,:] = [0, -LY1/2., -LZ/2.]
  vertices[6,:] = [0, -EPSCRACK/2., -LZ/2.]
  ###
  ### Create vertices
  ###
  v_ID=[]
  for v in vertices: 
    cubit.cmd("create vertex X %f Y %f Z %f" % (v[0], v[1], v[2]))
    v_ID.append(cubit.get_last_id("vertex"))
  ###
  ### Create curves
  ### 
  c_ID=[]
  for v in range(len(v_ID)):
    cubit.cmd("create curve vertex %i vertex %i" % (v_ID[v], v_ID[(v+1)%len(v_ID)]))
    c_ID.append(cubit.get_last_id("curve"))
  ### 
  ### create surface
  ###
  cmd="create surface curve"
  for c in c_ID:
    cmd+=" %i" % c
  cubit.cmd(cmd)
  s_ID=cubit.get_last_id("surface")
  ###
  ### sweep volume
  ###
  cubit.cmd("sweep surface %i direction Z distance %f" % (s_ID, LZ))
  DCB_3D=[cubit.get_last_id("volume")]
  return DCB_3D
  
def CylCrackCreate(lx, ly, lz, thetacrack):
   import cubit
   import numpy as np
   import pymef90
   ###
   ### Create Main body
   ###
   cubit.cmd ("Create cylinder height %f major radius %f minor radius %f" % (lz, lx, ly)) 
   CYL = cubit.get_last_id("volume")
   ###
   ### Create wedge
   ###
   w_ID=[]
   cubit.cmd("create vertex X 0 Y 0 Z %f" % (-lz))
   w_ID.append(cubit.get_last_id("vertex"))
   X = 1.1 * np.max(lx, ly) * np.cos(np.radians(thetacrack/2.))
   Y = 1.1 * np.max(lx, ly) * np.sin(np.radians(thetacrack/2.))
   cubit.cmd("create vertex X %f Y %f Z %f" % (-X, Y, -lz))
   w_ID.append(cubit.get_last_id("vertex"))
   cubit.cmd("create vertex X %f Y %f Z %f" % (-X,  -Y, -lz))
   w_ID.append(cubit.get_last_id("vertex"))
   cubit.cmd("create surface vertex %i %i %i" % (w_ID[0], w_ID[1], w_ID[2]))
   cubit.cmd("sweep surface %i perpendicular distance %f" % (cubit.get_last_id("surface"), 2*lz))
   WEDGE=cubit.get_last_id("volume")
   ###
   ### remove crack
   ###
   cubit.cmd("subtract volume %i from volume %i" % (WEDGE, CYL))
   ###
   ### Return the ID of the cylinder
   ###
   return CYL

def CylCrackCreateLayered(lx, ly, lz, alpha, theta1, theta2, thetacrack, lcrack=.5, ):
  import cubit
  import numpy as np
  import pymef90
  ###
  ### Create Main body
  ###
  cubit.cmd ("Create cylinder height %f major radius %f minor radius %f" % (lz, lx, ly)) 
  CYL = cubit.get_last_id("volume")
  ###
  ### Create wedge
  ###
  w_ID=[]
  cubit.cmd("create vertex X %f Y 0 Z %f" % ( (2. * lcrack - 1.) *lx, -lz))
  #  cubit.cmd("create vertex X 0 Y 0 Z %f" % ( -lz))
  w_ID.append(cubit.get_last_id("vertex"))
  X = 1.1 * np.max(lx, ly) * np.cos(np.radians(thetacrack/2.))
  Y = 1.1 * np.max(lx, ly) * np.sin(np.radians(thetacrack/2.))
  cubit.cmd("create vertex X %f Y %f Z %f" % (-X, Y, -lz))
  w_ID.append(cubit.get_last_id("vertex"))
  cubit.cmd("create vertex X %f Y %f Z %f" % (-X,  -Y, -lz))
  w_ID.append(cubit.get_last_id("vertex"))
  cubit.cmd("create surface vertex %i %i %i" % (w_ID[0], w_ID[1], w_ID[2]))
  cubit.cmd("sweep surface %i perpendicular distance %f" % (cubit.get_last_id("surface"), 2*lz))
  WEDGE=cubit.get_last_id("volume")
  ###
  ### remove crack
  ###
  cubit.cmd("subtract volume %i from volume %i" % (WEDGE, CYL))
  CYL_ID=[]
  CYL_ID.append(CYL)
  ###
  ### create layers
  ###
  bb=[-lx, lx, -ly, ly, -2.*lz, 2.*lz]
  print('CYL_ID, bb, alpha, theta1, theta2', CYL_ID, bb, alpha, theta1, theta2)
  (LAYER1_3D, LAYER2_3D) = pymef90.Layer(CYL_ID, bb, alpha, theta1, theta2, 0.)
  ###
  ### imprint and merge
  ###
  cubit.cmd("imprint all")
  cubit.cmd("merge all")
  ###
  ### groups
  ###
  pymef90.GroupAddVolList("LAYER1_3D", LAYER1_3D)
  pymef90.GroupAddVolList("LAYER2_3D", LAYER2_3D)
  ###
  ### return LAYER1_3D and LAYER2_3D
  ###
  return LAYER1_3D, LAYER2_3D
  
def CylCrackCreateLayeredCoin(R, lx, ly, lz, alpha, theta1, theta2, thetacrack, lcrack=.5, ):
  import cubit
  import numpy as np
  import pymef90
  ###
  ### Create Main body
  ###
  OUTSIDE_3D=[]
  cubit.cmd ("Create cylinder height %f radius %f" % (lz, R)) 
  OUTSIDE_3D.append(cubit.get_last_id("volume"))
  ###
  ### Create wedge
  ###
  w_ID=[]
  cubit.cmd("create vertex X %f Y 0 Z %f" % ( (2. * lcrack - 1.) * R, -lz))
  w_ID.append(cubit.get_last_id("vertex"))
  X = 1.1 * R * np.cos(np.radians(thetacrack/2.))
  Y = 1.1 * R * np.sin(np.radians(thetacrack/2.))
  cubit.cmd("create vertex X %f Y %f Z %f" % (-X, Y, -lz))
  w_ID.append(cubit.get_last_id("vertex"))
  cubit.cmd("create vertex X %f Y %f Z %f" % (-X,  -Y, -lz))
  w_ID.append(cubit.get_last_id("vertex"))
  cubit.cmd("create surface vertex %i %i %i" % (w_ID[0], w_ID[1], w_ID[2]))
  cubit.cmd("sweep surface %i perpendicular distance %f" % (cubit.get_last_id("surface"), 2*lz))
  WEDGE=cubit.get_last_id("volume")
  ###
  ### remove crack
  ###
  cubit.cmd("subtract volume %i from volume %i" % (WEDGE, OUTSIDE_3D[0]))
  ###
  ### Create center coin
  ###
  cubit.cmd ("Create cylinder height %f major radius %f minor radius %f" % (lz, lx, ly)) 
  CYL = cubit.get_last_id("volume")
  (OUTSIDE_3D, CYL_ID) = pymef90.WebcutTool(OUTSIDE_3D, CYL)
  ###
  ### create layers
  ###
  bb=[-lx, lx, -ly, ly, -2.*lz, 2.*lz]
  print('CYL_ID, bb, alpha, theta1, theta2', CYL_ID, bb, alpha, theta1, theta2)
  (LAYER1_3D, LAYER2_3D) = pymef90.Layer(CYL_ID, bb, alpha, theta1, theta2, 0.)
  ###
  ### imprint and merge
  ###
  cubit.cmd("imprint all")
  cubit.cmd("merge all")
  ###
  ### groups
  ###
  pymef90.GroupAddVolList("LAYER1_3D",  LAYER1_3D)
  pymef90.GroupAddVolList("LAYER2_3D",  LAYER2_3D)
  pymef90.GroupAddVolList("OUTSIDE_3D", OUTSIDE_3D)
  ###
  ### return LAYER1_3D and LAYER2_3D
  ###
  return OUTSIDE_3D, LAYER1_3D, LAYER2_3D

def GroupAddVolList(groupname, vol_IDs):
  import cubit
  cmd='group "%s" add volume ' % groupname
  for vol in vol_IDs:
    cmd += ' %i' % vol
  cubit.cmd(cmd)