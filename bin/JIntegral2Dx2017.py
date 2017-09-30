from visit import *
def parse(args=None):
    import argparse
    ### Get options from the command line
    parser = argparse.ArgumentParser(description='Compute the value of J')
    parser.add_argument('-i','--inputfile',help='rootdir and prefix of the file')
    parser.add_argument("--stepmin",type=int,help="first time step to process",default=0)
    parser.add_argument("--stepmax",type=int,help="last time step to process",default=None)
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
    
    #DefineVectorExpression("dUx","gradient(<Displacement_X>)")
    #DefineVectorExpression("dUy","gradient(<Displacement_Y>)")
    #
    #DefineScalarExpression("e11","dUx[0]")
    #DefineScalarExpression("e12","(dUx[1]+dUy[0])*0.5")
    #DefineScalarExpression("e22","dUy[1]")
    #DefineScalarExpression("sigma11","e11")
    #DefineScalarExpression("sigma12","e12")
    #DefineScalarExpression("sigma22","e22")
    #
    #
    #DefineScalarExpression("EnergyDensity","(sigma11 * e11 + 2. * sigma12 * e12 + sigma22 * e22) * 0.5")
    #DefineScalarExpression("C11","EnergyDensity - sigma11 * dUx[0] - sigma12 * dUy[0]")
    #DefineScalarExpression("C21","              - sigma12 * dUx[0] - sigma22 * dUy[0]")


    DefineVectorExpression("dUx","gradient(<Displacement_X>)")
    DefineVectorExpression("dUy","gradient(<Displacement_Y>)")
    DefineScalarExpression("EnergyDensity","(<Stress_XX> * dUx[0] + <Stress_YY> * dUy[1] + <Stress_XY> * (dUx[1]+dUy[0])) * .5 ")
    DefineScalarExpression("C11","EnergyDensity - <Stress_XX> * dUx[0] - <Stress_XY> * dUy[0]")
    DefineScalarExpression("C21","              - <Stress_XY> * dUx[0] - <Stress_YY> * dUy[0]")


    opts.bb = opts.bb
    AddPlot("Pseudocolor","Damage")
    DrawPlots()
    SuppressQueryOutputOn()

    filename = os.path.splitext(opts.inputfile)[0]+'_Jint.txt'
    f=open(filename,'w')
    f.write("#t   J, Jleft, Jtop, Jright, Jbottom\n")


    BB = opts.bb
    nintx = opts.nint[0]
    ninty = opts.nint[1]



    Cleft  = []
    Ctop   = []
    Cright = []
    Cbot   = []
    Jint   = []
    load   = []

    print "Computing J integral on 4 edges of a box:"
    ### -C e1 . n left edge
    ChangeActivePlotsVar("C11")
    DrawPlots()
    Lineout(start_point = (BB[0],BB[1]+BB[3],0), 
            end_point   = (BB[0],BB[1],      0)) 
    SetActiveWindow(2)
    firststep = opts.stepmin
    if opts.stepmax == None:
        laststep = TimeSliderGetNStates()
    else:
        laststep = opts.stepmin+opts.stepmax+1
    
    for s in range(firststep,laststep):
        SetTimeSliderState(s)
        Query("Time")
        t = GetQueryOutputValue()
        info = GetPlotInformation()
        Cleft.append(-mytrapz(info["Curve"]))
        load.append(t)
        sys.stdout.write("-")
        sys.stdout.flush()
    DeleteActivePlots()
    print "Done with left edge"

    ### C e1 . n right edge
    SetActiveWindow(1)
    Lineout(start_point = (BB[0]+BB[2],BB[1],      0), 
            end_point   = (BB[0]+BB[2],BB[1]+BB[3],0)) 
    SetActiveWindow(2)
    for s in range(firststep,laststep):
        SetTimeSliderState(s)
        info = GetPlotInformation()
        #print info["Curve"]
        #print '**** ', mytrapz(info["Curve"]), '***\n\n'
        Cright.append(mytrapz(info["Curve"]))
        sys.stdout.write("-")
        sys.stdout.flush()
    DeleteActivePlots()
    print "Done with right edge"

    ### C e1 . n top edge
    SetActiveWindow(1)
    ### -C e1 . n left edge
    ChangeActivePlotsVar("C21")
    Lineout(start_point = (BB[0]+BB[2],BB[1]+BB[3],0), 
            end_point   = (BB[0]      ,BB[1]+BB[3],0))
    SetActiveWindow(2)
    for s in range(firststep,laststep):
        SetTimeSliderState(s)
        info = GetPlotInformation()
        Ctop.append(mytrapz(info["Curve"]))
        sys.stdout.write("-")
        sys.stdout.flush()
    DeleteActivePlots()
    print "Done with top edge"

    ### C e1 . n bottom edge
    SetActiveWindow(1)
    Lineout(start_point = (BB[0]      ,BB[1],0), 
            end_point   = (BB[0]+BB[2],BB[1],0)) 
    SetActiveWindow(2)
    for s in range(firststep,laststep):
        SetTimeSliderState(s)
        info = GetPlotInformation()
        Cbot.append(-mytrapz(info["Curve"]))
        sys.stdout.write("-")
        sys.stdout.flush()
    DeleteActivePlots()
    print "Done with bottom edge"

    for s in range(firststep,laststep):
        Jint.append(Cleft[s]+Cright[s]+Ctop[s]+Cbot[s])
        print "****** step {0:d} load = {1:e}, Jint = {2:e}  ( = {3:e} + {4:e} + {5:e} + {6:e})".format(s,load[s],Jint[s],Cleft[s],Ctop[s],Cright[s],Cbot[s])
        f.write("{0:e}\t {1:e} \t{2:e} \t{3:e} \t{4:e} \t{5:e}\n".format(load[s],Jint[s],Cleft[s],Ctop[s],Cright[s],Cbot[s]))
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
    
