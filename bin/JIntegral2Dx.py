from visit import *

def parse(args=None):
    import argparse
    ### Get options from the command line
    parser = argparse.ArgumentParser(description='Compute the value of G_theta')
    parser.add_argument('-i','--inputfile',help='rootdir and prefix of the file')
    #parser.add_argument("--time_min",type=int,help="first time step",default=0)
    #parser.add_argument("--time_max",type=int,help="last time step",default=1)
    parser.add_argument("--E",type=float,help="Youngs modulus",default=1.)
    parser.add_argument("--nu",type=float,help="Poisson Ratio",default=0.)
    parser.add_argument('--bb',type=float,nargs=4,default=[-1,-1,2,2],help='bounding box of the layered area (xmin ymin lx ly)')
    parser.add_argument('--nint',type=int,nargs=2,default=[10,10],help='Number of grid points along axes of the opts.bb')

    return parser.parse_args()


def mytrapz(xy):
    trapz = 0
    for i in range(len(xy)/2-1):
        trapz += (xy[2*i+2]-xy[2*i])*(xy[2*i+1]+xy[2*i+3])
    trapz = trapz*.5
    return trapz

def plot(opts):
    import json
    import os
    import os.path
    import shutil
    import math
    
    ##  
    ## Open the database
    ##
    MyDatabase = opts.inputfile
    print "Trying to open database %s"%MyDatabase
    
    status = OpenDatabase(MyDatabase,0)       
    
    DefineScalarExpression("X", 'coord(Mesh)[0]')
    DefineScalarExpression("Y", 'coord(Mesh)[1]')
    DefineScalarExpression("E", '{0}* (1.-Damage) * (1.-Damage)'.format(opts.E))
    DefineScalarExpression("v", '{0}'.format(opts.nu))
    DefineVectorExpression("Displacement_2D", "{Displacement_X,Displacement_Y}")
    DefineScalarExpression("Strain_2D", "{{gradient(Displacement_X)[0],0.5*(gradient(Displacement_X)[1]+gradient(Displacement_Y)[0])},{0.5*(gradient(Displacement_X)[1]+gradient(Displacement_Y)[0]),gradient(Displacement_Y)[1]}} ")
    DefineTensorExpression("Grad_U", "{{gradient(Displacement_2D[0])[0],gradient(Displacement_2D[0])[1]},{gradient(Displacement_2D[1])[0],gradient(Displacement_2D[1])[1]}}")
    DefineTensorExpression("Stress_2D", "E/((1+v)*(1-2*v))*{{(1-v)*Strain_2D[0][0] , (1-2*v)*Strain_2D[0][1]},{(1-2*v)*Strain_2D[1][0] , (1-v)*Strain_2D[1][1]}}")
    DefineVectorExpression("dUx","gradient(<Displacement_X>)")
    DefineVectorExpression("dUy","gradient(<Displacement_Y>)")
    DefineScalarExpression("EnergyDensity","Stress_2D[0][0] * Grad_U[0][0] + Stress_2D[1][1] * Grad_U[1][1] + Stress_2D[0][1] * Grad_U[0][1] ")
    DefineScalarExpression("C11","EnergyDensity - Stress_2D[0][0] * Grad_U[0][0] - Stress_2D[0][1] * Grad_U[0][1] * 2.")
    DefineScalarExpression("C21","              - Stress_2D[0][1] * Grad_U[0][0] - Stress_2D[1][1] * Grad_U[0][0] * 2.")

    opts.bb = opts.bb
    AddPlot("Curve", "operators/Lineout/C11")
    DrawPlots()
    SuppressQueryOutputOn()

    filename = opts.inputfile+'_Jint.txt'
    f=open(filename,'w')
    f.write("#t   J, Jleft, Jtop, Jright, Jbottom\n")
    for s in range(0,TimeSliderGetNStates(),1):
        SetTimeSliderState(s)
        Query("Time")
        t = GetQueryOutputValue()

        ChangeActivePlotsVar("C11")
        lineAtts = LineoutAttributes()
        lineAtts.point1= (opts.bb[0],opts.bb[1]+opts.bb[3], 0)
        lineAtts.point2= (opts.bb[0],opts.bb[1], 0)
        lineAtts.interactive = 0
        lineAtts.ignoreGlobal = 1
        lineAtts.samplingOn = 1
        lineAtts.numberOfSamplePoints = opts.nint[1]
        SetOperatorOptions(lineAtts)
        Cleft = -mytrapz(GetPlotInformation()["Curve"])

        lineAtts.point1= (opts.bb[0]+opts.bb[2],opts.bb[1], 0)
        lineAtts.point2= (opts.bb[0]+opts.bb[2],opts.bb[1]+opts.bb[3], 0)
        lineAtts.interactive = 0
        lineAtts.ignoreGlobal = 1
        lineAtts.samplingOn = 1
        lineAtts.numberOfSamplePoints = opts.nint[1]
        SetOperatorOptions(lineAtts)
        Cright = mytrapz(GetPlotInformation()["Curve"])
        
        ChangeActivePlotsVar("C21")
        lineAtts.point1= (opts.bb[0]+opts.bb[2],opts.bb[1]+opts.bb[3], 0)
        lineAtts.point2= (opts.bb[0],       opts.bb[1]+opts.bb[3], 0)
        lineAtts.interactive = 0
        lineAtts.ignoreGlobal = 1
        lineAtts.samplingOn = 1
        lineAtts.numberOfSamplePoints = opts.nint[0]
        SetOperatorOptions(lineAtts)
        Ctop = mytrapz(GetPlotInformation()["Curve"])

        lineAtts.point1= (opts.bb[0],      opts.bb[1], 0)
        lineAtts.point2= (opts.bb[0]+opts.bb[2],opts.bb[1], 0)
        lineAtts.interactive = 0
        lineAtts.ignoreGlobal = 1
        lineAtts.samplingOn = 1
        lineAtts.numberOfSamplePoints = opts.nint[0]
        SetOperatorOptions(lineAtts)
        Cbot = -mytrapz(GetPlotInformation()["Curve"])
        
        print "****** step %i load = %e, Jint = %e  ( = %e %+e %+e %+e)"%(s,t,Cleft+Ctop+Cright+Cbot,Cleft,Ctop,Cright,Cbot)
        f.write("%e \t%e \t%e \t%e \t%e \t%e\n"%(t,Cleft+Ctop+Cright+Cbot,Cleft,Ctop,Cright,Cbot))
        f.flush()
        os.fsync(f)
    f.close()
    DeleteAllPlots()
    CloseDatabase(MyDatabase)
    return 0

if __name__ == "__main__":
    import sys  
    import os.path

    opts = parse()
    plot(opts)
    exit()

    
