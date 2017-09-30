from visit import *

def parse(args=None):
    import argparse
    ### Get options from the command line
    parser = argparse.ArgumentParser(description='Compute the value of J')
    parser.add_argument('-i','--inputfile',help='rootdir and prefix of the file')
    parser.add_argument("--step_min",type=int,help="first time step",default=None)
    parser.add_argument("--step_max",type=int,help="last time step",default=None)
    parser.add_argument("--E",type=float,help="Youngs modulus",default=1.)
    parser.add_argument("--nu",type=float,help="Poisson Ratio",default=0.)
    parser.add_argument('--bb',type=float,nargs=4,default=[-1,-1,2,2],help='bounding box of the layered area (xmin ymin lx ly)')
    parser.add_argument('--nint',type=int,nargs=2,default=[10,10],help='Number of grid points along axes of the opts.bb')

    return parser.parse_args()

def Jintegraldx2(BB,nintx=100,ninty=100):
    def mytrapz(xy):
        trapz = 0
        for i in range(len(xy)/2-1):
            trapz += (xy[2*i+2]-xy[2*i])*(xy[2*i+1]+xy[2*i+3])
        trapz = trapz*.5
        return trapz

    ### -C e1 . n left edge
    SetActiveWindow(1)
    ChangeActivePlotsVar("C11")
    DrawPlots()
    Query("Lineout", start_point=(BB[0],BB[1]+BB[3], 0), 
                       end_point=(BB[0],BB[1], 0), 
                    use_sampling=1,num_samples=ninty)
    SetActiveWindow(2)
    Cleft = -mytrapz(GetPlotInformation()["Curve"])
    DeleteActivePlots()


    ### C e1 . n right edge
    SetActiveWindow(1)
    Query("Lineout", start_point=(BB[0]+BB[2],BB[1], 0), 
                       end_point=(BB[0]+BB[2],BB[1]+BB[3], 0), 
                    use_sampling=1,num_samples=ninty)
    SetActiveWindow(2)
    Cright = mytrapz(GetPlotInformation()["Curve"])
    DeleteActivePlots()

    ### C e1 . n top edge
    SetActiveWindow(1)
    ChangeActivePlotsVar("C21")
    DrawPlots()
    Query("Lineout", start_point=(BB[0]+BB[2],BB[1]+BB[3], 0), 
                       end_point=(BB[0],      BB[1]+BB[3], 0), 
                    use_sampling=1,num_samples=nintx)
    SetActiveWindow(2)
    Ctop = mytrapz(GetPlotInformation()["Curve"])
    DeleteActivePlots()

    ### C e1 . n bottom edge
    SetActiveWindow(1)
    Query("Lineout", start_point=(BB[0],      BB[1], 0), 
                       end_point=(BB[0]+BB[2],BB[1], 0), 
                    use_sampling=1,num_samples=nintx)
    SetActiveWindow(2)
    Cbot = -mytrapz(GetPlotInformation()["Curve"])
    DeleteActivePlots()
    return Cleft,Ctop,Cright,Cbot

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
    
    DefineVectorExpression("dUx","gradient(<Displacement_X>)")
    DefineVectorExpression("dUy","gradient(<Displacement_Y>)")
    DefineScalarExpression("EnergyDensity","(<Stress_XX> * dUx[0] + <Stress_YY> * dUy[1] + <Stress_XY> * (dUx[1]+dUy[0])) * .5 ")
    DefineScalarExpression("C11","EnergyDensity - <Stress_XX> * dUx[0] - <Stress_XY> * dUy[0]")
    DefineScalarExpression("C21","              - <Stress_XY> * dUx[0] - <Stress_YY> * dUy[0]")

    #DefineScalarExpression("C11","EnergyDensity - <Stress_XY> * dUx[0] - <Stress_YY> * dUy[0]")
    #DefineScalarExpression("C21","              - <Stress_XX> * dUx[0] - <Stress_XY> * dUy[0]")

    DefineScalarExpression("C11","<Stress_XX>")
    DefineScalarExpression("C21","<Stress_YY>")

    opts.bb = opts.bb
    AddPlot("Pseudocolor","Damage")
    DrawPlots()
    SuppressQueryOutputOn()

    filename = opts.inputfile+'_Jint.txt'
    f=open(filename,'w')
    f.write("#t   J, Jleft, Jtop, Jright, Jbottom\n")


    BB = opts.bb
    nintx = opts.nint[0]
    ninty = opts.nint[1]
    if opts.step_max:
        laststep = opts.step_max
    else:
        laststep = TimeSliderGetNStates()
    if opts.step_min:
        firststep = opts.step_min
    else:
        firststep = 0
    for s in range(firststep,laststep,1):
        SetActiveWindow(1)
        SetTimeSliderState(s)
        Query("Time")
        t = GetQueryOutputValue()
        DrawPlots()
        Cleft,Ctop,Cright,Cbot = Jintegraldx2(BB,nintx,ninty)    

        print "****** step %i load = %e, Jint = %e  ( = %e %+e %+e %+e)"%(s,t,Cleft+Ctop+Cright+Cbot,Cleft,Ctop,Cright,Cbot)
        f.write("%e \t%e \t%e \t%e \t%e \t%e\n"%(t,Cleft+Ctop+Cright+Cbot,Cleft,Ctop,Cright,Cbot))
        f.flush()
        os.fsync(f)
    f.close()
    SetActiveWindow(1)
    DeleteAllPlots()
    SetActiveWindow(2)
    DeleteAllPlots()
    DeleteWindow()
    CloseDatabase(MyDatabase)
    return 0

if __name__ == "__main__":
    import sys  
    import os.path

    opts = parse()
    plot(opts)
    exit()

    
