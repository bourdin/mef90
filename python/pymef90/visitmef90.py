from visit import *

def SetAnnotations():
    AnnotationAtts = visit.GetAnnotationAttributes()
    AnnotationAtts.axes3D.xAxis.label.visible = 0
    AnnotationAtts.axes3D.yAxis.label.visible = 0
    AnnotationAtts.axes3D.zAxis.label.visible = 0
    AnnotationAtts.userInfoFlag = 0
    AnnotationAtts.databaseInfoFlag = 0
    AnnotationAtts.legendInfoFlag = 0
    AnnotationAtts.axesArray.visible = 0
    AnnotationAtts.databaseInfoExpansionMode = AnnotationAtts.SmartDirectory
    AnnotationAtts.axes3D.visible = 0
    AnnotationAtts.backgroundColor = (0, 0, 0, 255)
    AnnotationAtts.foregroundColor = (255, 255, 255, 255)
    visit.SetAnnotationAttributes(AnnotationAtts)
    return 0

def setView3DXup(thetadeg = 30,phideg = 45):
    import math
    
    theta = math.radians(thetadeg)
    phi = math.radians(phideg)
    View3DAtts = GetView3D()
    View3DAtts.SetViewNormal(math.sin(phi),math.cos(theta)*math.cos(phi),math.sin(theta)*math.cos(phi))
    View3DAtts.SetViewUp(1,0,0)
    SetView3D(View3DAtts)
    return 0

def setView3DYup(thetadeg = 30,phideg = 45):
    import math
    
    theta = math.radians(thetadeg)
    phi = math.radians(phideg)
    View3DAtts = GetView3D()
    View3DAtts.SetViewNormal(math.sin(theta)*math.cos(phi),math.sin(phi),math.cos(theta)*math.cos(phi))
    View3DAtts.SetViewUp(0,1,0)
    SetView3D(View3DAtts)
    return 0

def setView3DZup(thetadeg = 30,phideg = 45):
    import math
    
    theta = math.radians(thetadeg)
    phi = math.radians(phideg)
    View3DAtts = GetView3D()
    View3DAtts.SetViewNormal(math.cos(theta)*math.cos(phi),math.sin(theta)*math.cos(phi),math.sin(phi))
    View3DAtts.SetViewUp(0,0,1)
    SetView3D(View3DAtts)
    return 0

def SavePNG(prefix,geometry=[1024,768]):
    import shutil
    import os.path
    
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
    pngname = SaveWindow()
    if os.path.exists(pngname):
        print pngname, prefix+'.png'
        shutil.move(pngname,prefix+'.png')
        return 0
    else:
        return -1

