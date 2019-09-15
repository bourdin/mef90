from visit import *
def savePNG(filename,geometry=[4096,4096]):
    SaveWindowAtts = SaveWindowAttributes()
    SaveWindowAtts.outputToCurrentDirectory = 1
    SaveWindowAtts.fileName = filename
    SaveWindowAtts.family = 0
    SaveWindowAtts.format = SaveWindowAtts.PNG
    SaveWindowAtts.width = geometry[0]
    SaveWindowAtts.height = geometry[1]
    SaveWindowAtts.resConstraint = SaveWindowAtts.NoConstraint
    SaveWindowAtts.forceMerge = 1
    SetSaveWindowAttributes(SaveWindowAtts)
    return SaveWindow()

def SetAnnotations():
    import visit
    AnnotationAtts = AnnotationAttributes()
    AnnotationAtts.axes3D.visible = 0
    AnnotationAtts.axes3D.autoSetTicks = 1
    AnnotationAtts.axes3D.autoSetScaling = 1
    AnnotationAtts.axes3D.lineWidth = 0
    AnnotationAtts.axes3D.tickLocation = AnnotationAtts.axes3D.Inside  # Inside, Outside, Both
    AnnotationAtts.axes3D.axesType = AnnotationAtts.axes3D.ClosestTriad  # ClosestTriad, FurthestTriad, OutsideEdges, StaticTriad, StaticEdges
    AnnotationAtts.axes3D.triadFlag = 1
    AnnotationAtts.axes3D.bboxFlag = 0
    AnnotationAtts.axes3D.xAxis.title.visible = 0
    AnnotationAtts.axes3D.xAxis.title.font.font = AnnotationAtts.axes3D.xAxis.title.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axes3D.xAxis.title.font.scale = 1
    AnnotationAtts.axes3D.xAxis.title.font.useForegroundColor = 1
    AnnotationAtts.axes3D.xAxis.title.font.color = (0, 0, 0, 255)
    AnnotationAtts.axes3D.xAxis.title.font.bold = 0
    AnnotationAtts.axes3D.xAxis.title.font.italic = 0
    AnnotationAtts.axes3D.xAxis.title.userTitle = 0
    AnnotationAtts.axes3D.xAxis.title.userUnits = 0
    AnnotationAtts.axes3D.xAxis.title.title = "X-Axis"
    AnnotationAtts.axes3D.xAxis.title.units = ""
    AnnotationAtts.axes3D.xAxis.label.visible = 0
    AnnotationAtts.axes3D.xAxis.label.font.font = AnnotationAtts.axes3D.xAxis.label.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axes3D.xAxis.label.font.scale = 1
    AnnotationAtts.axes3D.xAxis.label.font.useForegroundColor = 1
    AnnotationAtts.axes3D.xAxis.label.font.color = (0, 0, 0, 255)
    AnnotationAtts.axes3D.xAxis.label.font.bold = 0
    AnnotationAtts.axes3D.xAxis.label.font.italic = 0
    AnnotationAtts.axes3D.xAxis.label.scaling = 0
    AnnotationAtts.axes3D.xAxis.tickMarks.visible = 1
    AnnotationAtts.axes3D.xAxis.tickMarks.majorMinimum = 0
    AnnotationAtts.axes3D.xAxis.tickMarks.majorMaximum = 1
    AnnotationAtts.axes3D.xAxis.tickMarks.minorSpacing = 0.02
    AnnotationAtts.axes3D.xAxis.tickMarks.majorSpacing = 0.2
    AnnotationAtts.axes3D.xAxis.grid = 0
    AnnotationAtts.axes3D.yAxis.title.visible = 0
    AnnotationAtts.axes3D.yAxis.title.font.font = AnnotationAtts.axes3D.yAxis.title.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axes3D.yAxis.title.font.scale = 1
    AnnotationAtts.axes3D.yAxis.title.font.useForegroundColor = 1
    AnnotationAtts.axes3D.yAxis.title.font.color = (0, 0, 0, 255)
    AnnotationAtts.axes3D.yAxis.title.font.bold = 0
    AnnotationAtts.axes3D.yAxis.title.font.italic = 0
    AnnotationAtts.axes3D.yAxis.title.userTitle = 0
    AnnotationAtts.axes3D.yAxis.title.userUnits = 0
    AnnotationAtts.axes3D.yAxis.title.title = "Y-Axis"
    AnnotationAtts.axes3D.yAxis.title.units = ""
    AnnotationAtts.axes3D.yAxis.label.visible = 0
    AnnotationAtts.axes3D.yAxis.label.font.font = AnnotationAtts.axes3D.yAxis.label.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axes3D.yAxis.label.font.scale = 1
    AnnotationAtts.axes3D.yAxis.label.font.useForegroundColor = 1
    AnnotationAtts.axes3D.yAxis.label.font.color = (0, 0, 0, 255)
    AnnotationAtts.axes3D.yAxis.label.font.bold = 0
    AnnotationAtts.axes3D.yAxis.label.font.italic = 0
    AnnotationAtts.axes3D.yAxis.label.scaling = 0
    AnnotationAtts.axes3D.yAxis.tickMarks.visible = 1
    AnnotationAtts.axes3D.yAxis.tickMarks.majorMinimum = 0
    AnnotationAtts.axes3D.yAxis.tickMarks.majorMaximum = 1
    AnnotationAtts.axes3D.yAxis.tickMarks.minorSpacing = 0.02
    AnnotationAtts.axes3D.yAxis.tickMarks.majorSpacing = 0.2
    AnnotationAtts.axes3D.yAxis.grid = 0
    AnnotationAtts.axes3D.zAxis.title.visible = 0
    AnnotationAtts.axes3D.zAxis.title.font.font = AnnotationAtts.axes3D.zAxis.title.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axes3D.zAxis.title.font.scale = 1
    AnnotationAtts.axes3D.zAxis.title.font.useForegroundColor = 1
    AnnotationAtts.axes3D.zAxis.title.font.color = (0, 0, 0, 255)
    AnnotationAtts.axes3D.zAxis.title.font.bold = 0
    AnnotationAtts.axes3D.zAxis.title.font.italic = 0
    AnnotationAtts.axes3D.zAxis.title.userTitle = 0
    AnnotationAtts.axes3D.zAxis.title.userUnits = 0
    AnnotationAtts.axes3D.zAxis.title.title = "Z-Axis"
    AnnotationAtts.axes3D.zAxis.title.units = ""
    AnnotationAtts.axes3D.zAxis.label.visible = 0
    AnnotationAtts.axes3D.zAxis.label.font.font = AnnotationAtts.axes3D.zAxis.label.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axes3D.zAxis.label.font.scale = 1
    AnnotationAtts.axes3D.zAxis.label.font.useForegroundColor = 1
    AnnotationAtts.axes3D.zAxis.label.font.color = (0, 0, 0, 255)
    AnnotationAtts.axes3D.zAxis.label.font.bold = 0
    AnnotationAtts.axes3D.zAxis.label.font.italic = 0
    AnnotationAtts.axes3D.zAxis.label.scaling = 0
    AnnotationAtts.axes3D.zAxis.tickMarks.visible = 1
    AnnotationAtts.axes3D.zAxis.tickMarks.majorMinimum = 0
    AnnotationAtts.axes3D.zAxis.tickMarks.majorMaximum = 1
    AnnotationAtts.axes3D.zAxis.tickMarks.minorSpacing = 0.02
    AnnotationAtts.axes3D.zAxis.tickMarks.majorSpacing = 0.2
    AnnotationAtts.axes3D.zAxis.grid = 0
    AnnotationAtts.axes3D.setBBoxLocation = 1
    #AnnotationAtts.axes3D.bboxLocation = (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5)
    AnnotationAtts.userInfoFlag = 0
    AnnotationAtts.userInfoFont.font = AnnotationAtts.userInfoFont.Arial  # Arial, Courier, Times
    AnnotationAtts.userInfoFont.scale = 1
    AnnotationAtts.userInfoFont.useForegroundColor = 1
    AnnotationAtts.userInfoFont.color = (0, 0, 0, 255)
    AnnotationAtts.userInfoFont.bold = 0
    AnnotationAtts.userInfoFont.italic = 0
    AnnotationAtts.databaseInfoFlag = 1
    AnnotationAtts.timeInfoFlag = 1
    AnnotationAtts.databaseInfoFont.font = AnnotationAtts.databaseInfoFont.Arial  # Arial, Courier, Times
    AnnotationAtts.databaseInfoFont.scale = 1
    AnnotationAtts.databaseInfoFont.useForegroundColor = 1
    AnnotationAtts.databaseInfoFont.color = (0, 0, 0, 255)
    AnnotationAtts.databaseInfoFont.bold = 0
    AnnotationAtts.databaseInfoFont.italic = 0
    AnnotationAtts.databaseInfoExpansionMode = AnnotationAtts.File  # File, Directory, Full, Smart, SmartDirectory
    AnnotationAtts.databaseInfoTimeScale = 1
    AnnotationAtts.databaseInfoTimeOffset = 0
    AnnotationAtts.legendInfoFlag = 1
    AnnotationAtts.backgroundColor = (255, 255, 255, 255)
    AnnotationAtts.foregroundColor = (0, 0, 0, 255)
    AnnotationAtts.gradientBackgroundStyle = AnnotationAtts.Radial  # TopToBottom, BottomToTop, LeftToRight, RightToLeft, Radial
    AnnotationAtts.gradientColor1 = (0, 0, 255, 255)
    AnnotationAtts.gradientColor2 = (0, 0, 0, 255)
    AnnotationAtts.backgroundMode = AnnotationAtts.Solid  # Solid, Gradient, Image, ImageSphere
    AnnotationAtts.backgroundImage = ""
    AnnotationAtts.imageRepeatX = 1
    AnnotationAtts.imageRepeatY = 1
    AnnotationAtts.axesArray.visible = 1
    AnnotationAtts.axesArray.ticksVisible = 1
    AnnotationAtts.axesArray.autoSetTicks = 1
    AnnotationAtts.axesArray.autoSetScaling = 1
    AnnotationAtts.axesArray.lineWidth = 0
    AnnotationAtts.axesArray.axes.title.visible = 1
    AnnotationAtts.axesArray.axes.title.font.font = AnnotationAtts.axesArray.axes.title.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axesArray.axes.title.font.scale = 1
    AnnotationAtts.axesArray.axes.title.font.useForegroundColor = 1
    AnnotationAtts.axesArray.axes.title.font.color = (0, 0, 0, 255)
    AnnotationAtts.axesArray.axes.title.font.bold = 0
    AnnotationAtts.axesArray.axes.title.font.italic = 0
    AnnotationAtts.axesArray.axes.title.userTitle = 0
    AnnotationAtts.axesArray.axes.title.userUnits = 0
    AnnotationAtts.axesArray.axes.title.title = ""
    AnnotationAtts.axesArray.axes.title.units = ""
    AnnotationAtts.axesArray.axes.label.visible = 1
    AnnotationAtts.axesArray.axes.label.font.font = AnnotationAtts.axesArray.axes.label.font.Arial  # Arial, Courier, Times
    AnnotationAtts.axesArray.axes.label.font.scale = 1
    AnnotationAtts.axesArray.axes.label.font.useForegroundColor = 1
    AnnotationAtts.axesArray.axes.label.font.color = (0, 0, 0, 255)
    AnnotationAtts.axesArray.axes.label.font.bold = 0
    AnnotationAtts.axesArray.axes.label.font.italic = 0
    AnnotationAtts.axesArray.axes.label.scaling = 0
    AnnotationAtts.axesArray.axes.tickMarks.visible = 1
    AnnotationAtts.axesArray.axes.tickMarks.majorMinimum = 0
    AnnotationAtts.axesArray.axes.tickMarks.majorMaximum = 1
    AnnotationAtts.axesArray.axes.tickMarks.minorSpacing = 0.02
    AnnotationAtts.axesArray.axes.tickMarks.majorSpacing = 0.2
    AnnotationAtts.axesArray.axes.grid = 0
    visit.SetAnnotationAttributes(AnnotationAtts)
    return 0

def getlaststep(fname):
  ### open file
  f=open(fname)
  ### Read last line in a string
  lastline = f.readlines()[-1]
  laststep = lastline.rsplit()[0] 
  return(int(laststep))

def drawCrack(displacementScaling=.1,damageThreshold=.99):
    ##
    ## Add pseudocolor plot of fracture field
    ##

    AddPlot('Pseudocolor', 'Damage')
    PseudocolorAtts = PseudocolorAttributes()
    PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
    PseudocolorAtts.skewFactor = 1
    PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
    PseudocolorAtts.minFlag = 1
    PseudocolorAtts.min = 0
    PseudocolorAtts.maxFlag = 1
    PseudocolorAtts.max = 1
    PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
    PseudocolorAtts.colorTableName = "hot"
    PseudocolorAtts.invertColorTable = 0
    PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
    PseudocolorAtts.renderSurfaces = 1
    PseudocolorAtts.renderWireframe = 0
    PseudocolorAtts.renderPoints = 0
    PseudocolorAtts.smoothingLevel = 0
    PseudocolorAtts.legendFlag = 1
    PseudocolorAtts.lightingFlag = 0
    SetPlotOptions(PseudocolorAtts)

    AddOperator("Displace", 1)
    SetActivePlots(0)
    SetActivePlots(0)
    DisplaceAtts = DisplaceAttributes()
    DisplaceAtts.factor = displacementScaling
    DisplaceAtts.variable = "Displacement"
    SetOperatorOptions(DisplaceAtts, 1)
    AddOperator("Isovolume", 1)
    IsovolumeAtts = IsovolumeAttributes()
    IsovolumeAtts.lbound = -1e+37
    IsovolumeAtts.ubound = damageThreshold
    IsovolumeAtts.variable = "Damage"
    SetOperatorOptions(IsovolumeAtts, 1)
    DrawPlots()
    Query("SpatialExtents", use_actual_data=1)
    BB = GetQueryOutputValue() 
    #SetView(BB)

    return BB
    

def setBGBlack():
    AnnotationAtts = GetAnnotationAttributes()
    AnnotationAtts.backgroundColor = (0, 0, 0, 255)
    AnnotationAtts.foregroundColor = (255, 255, 255, 255)
    SetAnnotationAttributes(AnnotationAtts)
    return 0

def setBGWhite():
    AnnotationAtts = GetAnnotationAttributes()
    AnnotationAtts.backgroundColor = (255, 255, 255, 255)
    AnnotationAtts.foregroundColor = (0, 0, 0, 255)
    SetAnnotationAttributes(AnnotationAtts)
    return 0

def setView3DZup(thetadeg = 30,phideg = 45,perspective=True):
    import math
    
    theta = math.radians(thetadeg)
    phi = math.radians(phideg)
    View3DAtts = GetView3D()
    View3DAtts.SetViewNormal(math.cos(theta)*math.cos(phi),math.sin(theta)*math.cos(phi),math.sin(phi))
    View3DAtts.SetViewUp(0,0,1)
    View3DAtts.SetPerspective(perspective)
    SetView3D(View3DAtts)
    return 0

def setView3DXup(thetadeg = 30,phideg = 45,zoom=1,perspective=True):
    import math
    
    theta = math.radians(thetadeg)
    phi   = math.radians(phideg)
    View3DAtts = GetView3D()
    n = (math.cos(phi),
         math.cos(theta)*math.sin(phi),
         math.sin(theta)*math.sin(phi))
    u = (math.sin(phi),
        -math.cos(theta)*math.cos(phi),
        -math.sin(theta)*math.cos(phi))
    View3DAtts.viewNormal = n
    View3DAtts.viewUp = u
    View3DAtts.focus = (0,0,0)
    View3DAtts.perspective = perspective
    View3DAtts.imageZoom = zoom
    SetView3D(View3DAtts)
    return 0

def setView3DYup(thetadeg = 30,phideg = 45,zoom=1,perspective=True):
    import math
    
    theta = math.radians(thetadeg)
    phi   = math.radians(phideg)
    View3DAtts = GetView3D()
    n = (math.sin(theta)*math.sin(phi),
         math.cos(phi),
         math.cos(theta)*math.sin(phi))
    u = (-math.sin(theta)*math.cos(phi),
          math.sin(phi),
         -math.cos(theta)*math.cos(phi))
    View3DAtts.viewNormal = n
    View3DAtts.viewUp = u
    View3DAtts.focus = (0,0,0)
    View3DAtts.perspective = perspective
    View3DAtts.imageZoom = zoom
    SetView3D(View3DAtts)
    return 0

def setView3DZup(thetadeg = 30,phideg = 45,zoom=1,perspective=True):
    import math
    
    theta = math.radians(thetadeg)
    phi   = math.radians(phideg)
    View3DAtts = GetView3D()
    n = (math.cos(theta)*math.sin(phi),
         math.sin(theta)*math.sin(phi),
         math.cos(phi))
    u = (-math.cos(theta)*math.cos(phi),
         -math.sin(theta)*math.cos(phi),
          math.sin(phi))
    View3DAtts.viewNormal = n
    View3DAtts.viewUp = u
    View3DAtts.focus = (0,0,0)
    View3DAtts.perspective = perspective
    View3DAtts.imageZoom = zoom
    SetView3D(View3DAtts)
    return 0

def plot(options):
    import json
    import os
    import os.path
    import shutil
    import math
    
    prefix,ext = os.path.splitext(options.inputfile)

    laststep = 1000000
    enerfile = prefix+'.ener'
    if options.step == 0:
        if os.path.exists(enerfile):
            laststep = getlaststep(enerfile)
        else:
            enerfile = prefix.split('_out')[0]+'.ener'
            if os.path.exists(enerfile):
                laststep = getlaststep(enerfile)
            else:
                print "unable to find step to plot."
                return -1
        step = laststep
    else:
        step = options.step

    ##  
    ## Open the database
    ##
    MyDatabase = options.inputfile
      
    print ('Trying to load {0}'.format(options.inputfile))
    status = OpenDatabase(options.inputfile, step-1)       
    if not status:
        print ("unable to open database {0}".format(options.inputfile))
        return -1

    BB = drawCrack(options.displacementScaling,options.damageThreshold)
    SetAnnotations()
    DrawPlots()

    if options.geometry == None:
        geometry = [2048,1536]
    else: 
        geometry = options.geometry

    if options.output != None:
        filename=os.path.splitext(options.output)[0]
    else:
        filename = '{basename}-{step:04d}'.format(basename = prefix,step=step)
    if options.bg == 'white':
        setBGWhite()
    else:
        setBGBlack()
    status = savePNG(filename,geometry)

    DeleteAllPlots()
    CloseDatabase(MyDatabase)
    return 0

def parse(args=None):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--step',type=int,default=0)
    parser.add_argument('--bg',choices=['white','black'],default='white')
    parser.add_argument('--displacementScaling',type=float,default=0)
    parser.add_argument('--damageThreshold',type=float,default=.99)
    parser.add_argument('--output',default=None)
    parser.add_argument('--geometry',default=None)
    parser.add_argument('inputfile',help='input file')
    return parser.parse_args()

if __name__ == "__main__":
    import sys  
    import os.path

    options = parse()
    if os.path.exists(options.inputfile):   
        print('processing {0}'.format(options.inputfile)) 
        plot(options)   
        sys.exit(0)
    else:
        print('Unable to find input file {0}'.format(options.inputfile))
        sys.exit(-1)

