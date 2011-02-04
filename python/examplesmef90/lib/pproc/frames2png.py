#        Visit Script to generate and save PNG images of a mef90-sieve output
#
#                  How to use
#
# 1) set the environmental variable FILE_PREFIX to set the file name:
#    e.g., for the file set ShrinkingSlab2D-X1-Y2-nc1-lc0.05-h0.025-XXXX.gen
#    "export FILE_PREFIX=ShrinkingSlab2D-X1-Y2-nc1-lc0.05-h0.025"
#    
# 1) run VisIt IN THE DIRECTORY WHERE THE FILE IS LOCATED with Command Line Interface and No Windows:
#    "visit -cli -nowin -s Visit_FracturePNG.py"
#
# Note: the script identify the database name from the .e2c filename 
# To do: consider the case qith more than 1 .e2c file
#
import os
import sys
import glob
# --------------------------------------------------
# Step 0; Files and directories
# --------------------------------------------------
## Get the working directory
#Directory=os.environ.get('WORKDIR')+'/'
Directory=os.environ.get('PWD') + '/'
## Get the prefix for the file name
#############
Prefix=glob.glob('*.e2c')[0].split('.e2c')[0]
#Prefix=os.environ.get('GENPREFIX')
## Create directory to store images
#ImageDir=Directory+'Images_'+Prefix+'/'
#if not os.path.isdir(ImageDir):
#  os.mkdir(ImageDir)
##  
## Step 1: Open the database
##
MyDatabase= Prefix+'-*.gen database'
OpenDatabase(MyDatabase,0)
# --------------------------------------------------
# Step 2: Add plots and set their properties
# --------------------------------------------------
AddPlot("Mesh", "Mesh", 1, 1)
Mesh = MeshAttributes()
Mesh.legendFlag = 1
Mesh.lineStyle = MeshAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
Mesh.lineWidth = 0
Mesh.meshColor = (255, 255, 255, 255)
Mesh.outlineOnlyFlag = 0
Mesh.errorTolerance = 0.01
Mesh.meshColorSource = MeshAtts.MeshCustom  # Foreground, MeshCustom
Mesh.opaqueColorSource = MeshAtts.Background  # Background, OpaqueCustom
Mesh.opaqueMode = MeshAtts.Auto  # Auto, On, Off
Mesh.pointSize = 0.05
Mesh.opaqueColor = (255, 255, 255, 255)
Mesh.smoothingLevel = MeshAtts.None  # None, Fast, High
Mesh.pointSizeVarEnabled = 0
Mesh.pointSizeVar = "default"
Mesh.pointType = MeshAtts.Point  # Box, Axis, Icosahedron, Point, Sphere
Mesh.showInternal = 0
Mesh.pointSizePixels = 2
Mesh.opacity = 0.305882
SetPlotOptions(Mesh)

AddPlot("Pseudocolor", "Fracture", 1, 1)
Frac = PseudocolorAttributes()
Frac.legendFlag = 1
Frac.lightingFlag = 1
Frac.minFlag = 1
Frac.maxFlag = 1
Frac.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
Frac.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
Frac.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
Frac.min = 0
Frac.max = 1
Frac.pointSize = 0.05
Frac.pointType = PseudocolorAtts.Point  # Box, Axis, Icosahedron, Point, Sphere
Frac.skewFactor = 1
Frac.opacity = 1
Frac.colorTableName = "grayalpha"
Frac.smoothingLevel = 0
Frac.pointSizeVarEnabled = 0
Frac.pointSizeVar = "default"
Frac.pointSizePixels = 2
Frac.lineStyle = PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
Frac.lineWidth = 0
Frac.opacityType = PseudocolorAtts.ColorTable  # Explicit, ColorTable
SetPlotOptions(Frac)
AddPlot("Pseudocolor", "Delamination", 1, 1)
Delam = PseudocolorAttributes()
Delam.legendFlag = 1
Delam.lightingFlag = 1
Delam.minFlag = 1
Delam.maxFlag = 1
Delam.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
Delam.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
Delam.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
Delam.min = 0
Delam.max = 1
Delam.pointSize = 0.05
Delam.pointType = PseudocolorAtts.Point  # Box, Axis, Icosahedron, Point, Sphere
Delam.skewFactor = 1
Delam.opacity = 1
Delam.colorTableName = "hot"
Delam.smoothingLevel = 0
Delam.pointSizeVarEnabled = 0
Delam.pointSizeVar = "default"
Delam.pointSizePixels = 2
Delam.lineStyle = PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
Delam.lineWidth = 0
Delam.opacityType = PseudocolorAtts.Explicit  # Explicit, ColorTable
SetPlotOptions(Delam)
DrawPlots()

# Set the view
#v = GetView2D()
#windowCoords = (-0.776642, 1.70354, -0.425681, 0.626574)
#viewportCoords = (0, 2, 0, 1)
#SetView2D(v)
AnnotationAtts = AnnotationAttributes()
AnnotationAtts.axes2D.visible = 0
AnnotationAtts.userInfoFlag = 0
AnnotationAtts.databaseInfoFlag = 0
AnnotationAtts.legendInfoFlag = 0
AnnotationAtts.axesArray.visible = 0
SetAnnotationAttributes(AnnotationAtts)
# --------------------------------------------------
# Step 3: Set Save file parameters 
# --------------------------------------------------
s = SaveWindowAttributes()
s.format = s.PNG
s.saveTiled = 0
s.family=0
s.width = 600
s.height = 200
s.screenCapture = 0
s.progressive = 1
forceMerge = 1
s.quality = 80
#
# Save images of all time steps and add each image filename to a list.
#
names = []
for state in range(1,TimeSliderGetNStates()):
   SetTimeSliderState(state)
   # Save the image
   s.fileName = Directory+ Prefix+'-%04i'%state
   SetSaveWindowAttributes(s)
   n = SaveWindow()
   names = names + [n]
   if opts.debug:
	print "Processing image for frame %s" %state
   #print "Processing image for frame %s" %n
#print names
exit()
