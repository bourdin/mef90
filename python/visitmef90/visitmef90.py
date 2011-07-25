###############################################################################
# Visit Script to generate and save PNG images of a mef90-sieve output
#
# How to use:
#  
#  run VisIt IN THE DIRECTORY WHERE THE FILE IS LOCATED with Command Line Interface and No Windows:
#
# Author: Corrado Maurini: cmaurini@gmail.com
###############################################################################

import os; import sys; import glob;    import math;
import pymef90;  from visit import * 


def FigureFracture(Prefix,ImageOptions=None,state=1):
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
    if not ImageOptions == None:
        # Set the view
        v = GetView2D()
        v.windowCoords = (ImageOptions['x0'], ImageOptions['lx'],ImageOptions['y0'],ImageOptions['ly'])
        v.viewportCoords = (ImageOptions['x0v'],ImageOptions['lxv'],ImageOptions['y0v'],ImageOptions['lyv'])
        SetView2D(v)
    AnnotationAtts = AnnotationAttributes()
    AnnotationAtts.axes2D.visible = 0
    AnnotationAtts.userInfoFlag = 0
    AnnotationAtts.databaseInfoFlag = 0
    AnnotationAtts.legendInfoFlag = 0
    AnnotationAtts.axesArray.visible = 0
    SetAnnotationAttributes(AnnotationAtts)
    # --------------------------------------------------
    SetTimeSliderState(state)

      
def ExportSingleFigure(prefix,ImageOptions,state=1):
    aspectratio=ImageOptions['ly']/ImageOptions['lx']
    SetTimeSliderState(state)
    s = SaveWindowAttributes()
    s.fileName = prefix
    s.format = s.PNG
    s.saveTiled = 1
    s.family=0
    s.width =  ImageOptions['res']
    s.height =  math.ceil(ImageOptions['res']*aspectratio)
    s.screenCapture = 0
    s.progressive = 0
    s.quality = 80
    s.resConstraint = 1
    SetSaveWindowAttributes(s)    
    print "Saving Image %s.png" %prefix
    n = SaveWindow()
    
    
def MakeFiguresLastStep(ImageOptions=None):
    energyfiles=glob.glob('*.ener')
    
    if len(energyfiles)==0:
      print "It seems there are not files to process"
    else:
      for energyfilename in energyfiles:
            # gen prefix
            prefix=energyfilename.split('.ener')[0]
            
            # get last step
            nsteps = pymef90.getlaststep(prefix+'.ener')
            
            # Represt Fracture 
            FigureFracture(prefix,ImageOptions,nsteps-1) 
            
            
            # Save image of a single time step in the JOB directory  
            ExportSingleFigure(prefix,ImageOptions,nsteps- 1) 
            


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
       


def MakeFiguresAllSteps(savedir='PNG'):
    energyfiles=glob.glob('*.ener')
    
    if len(energyfiles)==0:
      print "It seems there are not files to process"
    else:
      for energyfilename in energyfiles:
            # gen prefix
            prefix=energyfilename.split('.ener')[0]
            
            # get last step
            nsteps = pymef90.getlaststep(prefix+'.ener')
            
            # Represt Fracture 
            FigureFracture(prefix,nsteps-1) 
            
            savedir='PNG'
            saveprefix=savedir+'/'+prefix
                        
            # Delete images if not up to date
            if os.path.isdir(savedir):
                if os.path.getmtime(savedir)<os.path.getmtime(energyfilename):
                    ExportTimeFigure(saveprefix)
                else:
                    print ' --------------------------------------------------------------'
                    print "| The figures for all time steps are up to date. Nothing to do |'" 
                    print ' --------------------------------------------------------------'              
            else:
                    os.mkdir(savedir)  
                    ExportTimeFigure(saveprefix)
