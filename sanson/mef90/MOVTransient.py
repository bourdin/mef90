#!/usr/bin/env python
import infotxt
import sys
import os

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
    AnnotationAtts.databaseInfoExpansionMode = AnnotationAtts.SmartDirectory
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
    import math
    if os.path.exists('00_INFO.txt'):
        Param = pymef90.Dictreadtxt('00_INFO.txt')
        prefix = str(Param['prefix'])
        enerfile = prefix+'.ener'
        laststep = pymef90.energies.getlaststep(enerfile)
    
        ##  
        ## Open the database
        ##
        if os.path.exists(prefix+'-0001.gen'):
          MyDatabase = prefix+'-*.gen database'
        else:
          MyDatabase = prefix+'-0000.gen'

        status = OpenDatabase(MyDatabase,0)
        print MyDatabase, status
        
        if not status:
            print "unable to open database %s"%MyDatabase
            return -1        

        dim = GetMetaData(MyDatabase).GetMeshes(0).spatialDimension
        laststep = TimeSliderGetNStates()-1
        SetTimeSliderState(laststep-1)

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
        
        SetPlotOptions(p)
    
        bb = GetBB(0)
        
        if dim == 2:
            pass
        else:
            SetAnnotations3D()
            SetActivePlots(0)
            AddOperator("Isovolume", 0)
            IsovolumeAtts = IsovolumeAttributes()
            IsovolumeAtts.lbound = -1e+37
            IsovolumeAtts.ubound = 0.1
            IsovolumeAtts.variable = "Fracture"
            SetOperatorOptions(IsovolumeAtts, 0)

            phi = math.radians(20)
            theta = math.radians(30)
            View3DAtts = View3DAttributes()
            View3DAtts.viewNormal = (math.cos(theta)*math.cos(phi),math.sin(theta)*math.cos(phi),math.sin(phi))
            # use View3DAtts.viewNormal = (math.cos(theta)*math.cos(phi),math.sin(phi),math.sin(theta)*math.cos(phi))
            # for Brick Y-up computations
            View3DAtts.focus = ((bb[0]+bb[1])/2.,(bb[2]+bb[3])/2.,(bb[4]+bb[5])/2.)
            View3DAtts.viewUp = (0,0,1)

            View3DAtts.viewAngle = 39 #horizontal angle of view for a 50mm 24x36 camera
            View3DAtts.parallelScale = 1
            View3DAtts.nearPlane = -1.5*math.sqrt(bb[0]*bb[0] + bb[2]*bb[2] + bb[4]*bb[4])
            View3DAtts.farPlane = 1.5*math.sqrt(bb[1]*bb[1] + bb[3]*bb[3] + bb[5]*bb[5])
            View3DAtts.imagePan = (0, 0)
            View3DAtts.imageZoom = 1
            View3DAtts.perspective = 1
            View3DAtts.eyeAngle = 2
            View3DAtts.centerOfRotationSet = 0
            View3DAtts.centerOfRotation = ((bb[0]+bb[1])/2.,(bb[2]+bb[3])/2.,(bb[4]+bb[5])/2.)
            View3DAtts.axis3DScaleFlag = 0
            View3DAtts.axis3DScales = (1, 1, 1)
            View3DAtts.shear = (0, 0, 1)
            SetView3D(View3DAtts)
            ViewAxisArrayAtts = ViewAxisArrayAttributes()
            ViewAxisArrayAtts.domainCoords = (0, 1)
            ViewAxisArrayAtts.rangeCoords = (0, 1)
            SetViewAxisArray(ViewAxisArrayAtts)



        InvertBackgroundColor()        

        tmpdir = tempfile.mkdtemp()
        ### generate individual frames
        for step in range(laststep):
            SetTimeSliderState(step)
            DrawPlots()
            pngname = SavePNG(os.path.join(tmpdir,prefix+"-",geometry)
    
        ### use ffmpeg to generate animation
        pattern = os.path.join(tmpdir,PREFIX)+"-%04d.png"
        cmd = "ffmpeg -y -i %s -vcodec mjpeg -qscale 0 %s-Transient.avi"%(pattern,prefix)
        print "Now running %s"%cmd
        os.system(cmd)
        shutil.rmtree(tmpdir)


import sys  
if __name__ == "__main__":
    #main()
    sys.exit(main())
