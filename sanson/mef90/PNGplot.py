#!/usr/bin/env python
import infotxt
import sys
import os

def getlaststep(fname):
  ### open file
  f=open(fname)
  ### Read last line in a string
  lastline = f.readlines()[-1]
  laststep = lastline.rsplit()[0] 
  return(int(laststep))

def SavePNG(prefix,geometry=[1920,1080]):
    SaveWindowAtts = SaveWindowAttributes()
    SaveWindowAtts.outputToCurrentDirectory = 1
    SaveWindowAtts.fileName = prefix
    SaveWindowAtts.family = 1
    SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
    SaveWindowAtts.width = geometry[0]
    SaveWindowAtts.height = geometry[1]
    SaveWindowAtts.resConstraint = SaveWindowAtts.NoConstraint
    SaveWindowAtts.forceMerge = 1
    SetSaveWindowAttributes(SaveWindowAtts)
    name = SaveWindow()
    return name

def SetAnnotations():
    # Logging for SetAnnotationObjectOptions is not implemented yet.
    AnnotationAtts = AnnotationAttributes()
    AnnotationAtts.axes3D.yAxis.label.visible = 0
    AnnotationAtts.userInfoFlag = 0
    AnnotationAtts.databaseInfoFlag = 0
    AnnotationAtts.legendInfoFlag = 0
    AnnotationAtts.axesArray.visible = 0
    SetAnnotationAttributes(AnnotationAtts)

def SetAnnotations3D():
    # Logging for SetAnnotationObjectOptions is not implemented yet.
    AnnotationAtts = AnnotationAttributes()
    AnnotationAtts.axes3D.xAxis.label.visible = 0
    AnnotationAtts.axes3D.yAxis.label.visible = 0
    AnnotationAtts.axes3D.zAxis.label.visible = 0
    AnnotationAtts.userInfoFlag = 0
    AnnotationAtts.databaseInfoFlag = 0
    AnnotationAtts.legendInfoFlag = 0
    AnnotationAtts.axesArray.visible = 0
    AnnotationAtts.databaseInfoExpansionMode = AnnotationAtts.SmartDirectory
    AnnotationAtts.axes3D.visible = 0
    SetAnnotationAttributes(AnnotationAtts)
    
def main():
    import os.path
    import shutil
    if os.path.exists('00_INFO.txt'):
        Param = infotxt.Dictreadtxt('00_INFO.txt')
    
        ##  
        ## Open the database
        ##
        MyDatabase = str(Param['JOBID'])+'.xmf'
        status = OpenDatabase(MyDatabase,0)
        if not status:
            MyDatabase = str(Param['JOBID'])+'.*.xmf database'
            status = OpenDatabase(MyDatabase,0)
            if not status:
                print "Cannot open database, exiting"
                return 0
        ##
        ## Add pseudocolor plot of fracture field
        ##
        AddPlot('Pseudocolor', 'Fracture')
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
        p.smoothingLevel = 2
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
        
        
        
        laststep = TimeSliderGetNStates()-1
        SetTimeSliderState(laststep-1)
        SetPlotOptions(p)
    

        if Param['NY'] == 2:
            dim = 2
        else:
            dim = 3

        if dim == 2:
            SetAnnotations()
            View3DAtts = View3DAttributes()
            View3DAtts.viewNormal = (0,-1.,0.)
            View3DAtts.focus = (4, 0.005, 4)
            View3DAtts.viewUp = (0,0,1)
            View3DAtts.viewAngle = 30
            View3DAtts.parallelScale = 5.65686
            View3DAtts.nearPlane = -11.3137
            View3DAtts.farPlane = 11.3137
            View3DAtts.imagePan = (0, 0)
            View3DAtts.imageZoom = 1
            View3DAtts.perspective = 1
            View3DAtts.eyeAngle = 2
            View3DAtts.centerOfRotationSet = 0
            View3DAtts.centerOfRotation = (4, 0.005, 4)
            View3DAtts.axis3DScaleFlag = 0
            View3DAtts.axis3DScales = (1, 1, 1)
            View3DAtts.shear = (0, 0, 1)
            SetView3D(View3DAtts)
        
        
            geometry = [0,0]
            if Param['NX'] > Param['NZ']:
                geometry[0] = 2000;
                geometry[1] = int(2000 * Param['NZ'] / Param['NX'])
            else:
                geometry[0] = int(2000 * Param['NX'] / Param['NZ'])
                geometry[1] = 2000;
        else:
            SetAnnotations3D()
            View3DAtts = View3DAttributes()
            View3DAtts.viewNormal = (0.657417, 0.711041, 0.249448)
            View3DAtts.focus = (4, 4, 4)
            View3DAtts.viewUp = (-0.166863, -0.18545, 0.968383)
            View3DAtts.viewAngle = 15
            View3DAtts.parallelScale = 6.9282
            View3DAtts.nearPlane = -13.8564
            View3DAtts.farPlane = 13.8564
            View3DAtts.imagePan = (0, 0)
            View3DAtts.imageZoom = 1
            View3DAtts.perspective = 1
            View3DAtts.eyeAngle = 2
            View3DAtts.centerOfRotationSet = 0
            View3DAtts.centerOfRotation = (4, 4, 4)
            View3DAtts.axis3DScaleFlag = 0
            View3DAtts.axis3DScales = (1, 1, 1)
            View3DAtts.shear = (0, 0, 1)
            SetView3D(View3DAtts)
            ViewAxisArrayAtts = ViewAxisArrayAttributes()
            ViewAxisArrayAtts.domainCoords = (0, 1)
            ViewAxisArrayAtts.rangeCoords = (0, 1)
            ViewAxisArrayAtts.viewportCoords = (0.15, 0.9, 0.1, 0.85)
            SetViewAxisArray(ViewAxisArrayAtts)

            AddOperator("Isosurface", 1)
            SetActivePlots(0)
            IsosurfaceAtts = IsosurfaceAttributes()
            IsosurfaceAtts.contourNLevels = 10
            IsosurfaceAtts.contourValue = (0.1)
            IsosurfaceAtts.contourPercent = ()
            IsosurfaceAtts.contourMethod = IsosurfaceAtts.Value  # Level, Value, Percent
            IsosurfaceAtts.minFlag = 0
            IsosurfaceAtts.min = 0
            IsosurfaceAtts.maxFlag = 0
            IsosurfaceAtts.max = 1
            IsosurfaceAtts.scaling = IsosurfaceAtts.Linear  # Linear, Log
            IsosurfaceAtts.variable = "Fracture"
            SetOperatorOptions(IsosurfaceAtts, 1)

        InvertBackgroundColor()        
        DrawPlots()
        pngname = SavePNG(str(Param['JOBID']))
        shutil.move(pngname,Param['JOBID']+'.png')


import sys  
if __name__ == "__main__":
    #main()
    sys.exit(main())
