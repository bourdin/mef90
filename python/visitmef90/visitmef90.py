###############################################################################
# Visit Script to generate and save PNG images of a mef90-sieve output
#
# How to use:
#  
#  run VisIt IN THE DIRECTORY WHERE THE FILE IS LOCATED with Command Line Interface and No Windows:
#    "visit -cli -nowin -s $MEF_HOME/python/visitmef90/Visit_FracturePNG2.py"
#
# Note: the script identify the database name from the .e2c filename 
#
# Author: Corrado Maurini: cmaurini@gmail.com
###############################################################################

import os; import sys; import glob; 
import pymef90;  from visit import * 

def FigureFracture(Prefix,state=1):
    print Prefix

    ##  
    ## Step 1: Open the database
    ##
    if os.path.exists(Prefix+'-0001.gen'):
      MyDatabase= Prefix+'-*.gen database'
    else:
      MyDatabase= Prefix+'-0000.gen'
    OpenDatabase(MyDatabase,0)
    # --------------------------------------------------
    # Step 2: Add plots and set their properties
    # --------------------------------------------------
    AddPlot('Pseudocolor', 'Fracture')
    #AddPlot('Mesh','mesh')
    p = PseudocolorAttributes()
    p.lightingFlag = 1
    p.centering = p.Natural  # Natural, Nodal, Zonal
    p.scaling = p.Linear  # Linear, Log, Skew
    p.limitsMode = p.OriginalData  # OriginalData, CurrentPlot
    p.pointSize = 0.05
    p.pointType = p.Point  # Box, Axis, Icosahedron, Point, Sphere
    p.skewFactor = 1
    p.opacity = 1
    p.colorTableName = "hot"
    p.invertColorTable = 1
    p.smoothingLevel = 0
    p.pointSizeVarEnabled = 0
    p.pointSizeVar = "default"
    p.pointSizePixels = 2
    p.lineStyle = p.SOLID  # SOLID, DASH, DOT, DOTDASH
    p.lineWidth = 0
    p.opacityType = p.Explicit  # Explicit, ColorTable
    # Set the min/max values
    p.minFlag = 1
    p.maxFlag = 1
    p.min=0.0
    p.max=1.0
    p.legendFlag=0
    SetPlotOptions(p)
    # Step 3: Draw the plots
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
    SetTimeSliderState(state)

      
def ExportSingleFigure(prefix,state=1):
    SetTimeSliderState(state)
    s = SaveWindowAttributes()
    s.fileName = prefix
    s.format = s.PNG
    s.saveTiled = 0
    s.family=0
    s.width = 600
    s.height = 200
    s.screenCapture = 0
    s.progressive = 1
    s.quality = 80
    SetSaveWindowAttributes(s)    
    print "Saving Image %s.png" %prefix
    n = SaveWindow()

def ExportTimeFigure(prefix):
    names = []
    for state in range(0,TimeSliderGetNStates()):
       SetTimeSliderState(state)
       # Save the image
       FigureName = prefix+'-%04i'%state      
       ExportSingleFigure(FigureName,state)
       n = SaveWindow()
       names = names + [n]
       print "Processing image for frame %s" %state
