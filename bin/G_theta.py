#visit -cli -nowin -s G_theta.py 
#visit -cli -nowin -np 12 -l srun -s G_theta.py

#! /usr/bin/env python
from visit import *


def parse(args=None):
    import argparse
    ### Get options from the command line
    parser = argparse.ArgumentParser(description='Compute the value of G_theta')
    parser.add_argument('-i','--inputfile',help='rootdir and prefix of the file')
    parser.add_argument("--R_inn",type=float,help="interior radius of the control area",default=.5)
    parser.add_argument("--R_out",type=float,help="outer radius of the control area",default=1.)
    parser.add_argument('--Angle_crack_propagation',type=float,help="crack propagation angle made with the x direction in degree",default=0.)
    parser.add_argument('--Crack_position',type=float,nargs=2,default=[0,0],help='postion of the crack in X,Y')
    #parser.add_argument('--cs',type=int,nargs='*',default=1,help='set of blocks where G theta is applied, check on visit to get the number of the set')
    return parser.parse_args()




def Get_G(opts):
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
    status = OpenDatabase(MyDatabase)
    
    print "Assumptions: \n"
    print "(i) 2D domain with a straigth crack (can be extended to a 3D, but it is not supported in this version)\n"
    print "(ii) homogenous material and cracks stay in the same material property. Do not apply if the crack cross different material properties \n"
    print "(iii) free stress at crack lips + no lips interpentration (it can be expressed with pressure applied on crack lips, but it is not supported in this version)\n"
    print "(iv) no force density (it can be expressed with force density, but it is not supported in this version )\n"

    #### Define displacement, stress, and strain
    DefineVectorExpression("Displacement_2D", "{Displacement_X,Displacement_Y}")
    DefineTensorExpression("Stress_Tensor_2D", "{{Stress_XX,Stress_XY},{Stress_XY,Stress_YY}}")
    DefineTensorExpression("GradU", "{{gradient(Displacement_2D[0])[0],gradient(Displacement_2D[0])[1]},{gradient(Displacement_2D[1])[0],gradient(Displacement_2D[1])[1]}}")

    #### define the controlled area for the 
    DefineVectorExpression("Crack_Pos", '{'+str(opts.Crack_position[0])+','+str(opts.Crack_position[1])+'}')
    DefineVectorExpression("Coords", "coord(Mesh)-Crack_Pos")
    DefineScalarExpression("R", "sqrt( (Coords[0])^2.+(Coords[1])^2.)")
    DefineVectorExpression("InnerOuterRadius", '{'+str(opts.R_inn)+','+str(opts.R_out)+'}')

    #### adding terms for force density = -Crack_Pressure*gradient(Damage)
    #DefineVectorExpression("ThetaVector", "-Crack_Pressure*gradient(Damage)")
    #DefineVectorExpression("GradFdotThetaVector","{gradient(CrackPressureForce[0])[0]*ThetaVector[0]+gradient(CrackPressureForce[0])[1]*ThetaVector[0],gradient(CrackPressureForce[1])[0]*ThetaVector[0]+gradient(CrackPressureForce[1])[1]*ThetaVector[1]}")
    #DefineScalarExpression("G_theta_Density","trace(Stress_Tensor_2D*(GradU*GradThetaVector))-.5*trace(Stress_Tensor_2D*GradU)*divergence(ThetaVector)+dot(CrackPressureForce,Displacement_2D)*divergence(ThetaVector)+dot(Displacement_2D,GradFdotThetaVector)")

    #### adding terms when crack is pressurized
    #### G_theta_Density + int_{\Gamma(Rout)-\Gamma(Rint])} Gamma_pressure_crack ds
    # DefineScalarExpression("Gamma_pressure_crack","dot(normale_crack*Displacement_2D)*divergence(ThetaVector)")

    #### Define the vector theta and its gradients
    DefineScalarExpression("Theta_function", "if( le(R, InnerOuterRadius[0]), 1 , if(ge(R, InnerOuterRadius[1]), 0, (R - InnerOuterRadius[1] )/( InnerOuterRadius[0] - InnerOuterRadius[1]) ))")
    DefineVectorExpression("ThetaVector", '{Theta_function*cos('+str(opts.Angle_crack_propagation*180/math.pi)+'),Theta_function*sin('+str(opts.Angle_crack_propagation*180/math.pi)+')}')
    DefineTensorExpression("GradThetaVector", "{{gradient(ThetaVector[0])[0],gradient(ThetaVector[0])[1]},{gradient(ThetaVector[1])[0],gradient(ThetaVector[1])[1]}}")
    
    DefineScalarExpression("G_theta_Density", "trace(Stress_Tensor_2D*(GradU*GradThetaVector))-.5*trace(Stress_Tensor_2D*GradU)*divergence(ThetaVector)")



    AddPlot("Pseudocolor","G_theta_Density")
    #silr = SILRestriction()
    #silr.SuspendCorrectnessChecking()
    #silr.TurnOffAll()
    #for var_cs in opts.cs:
    #    silr.TurnOnSet(var_cs)
    #silr.EnableCorrectnessChecking()
    #SetPlotSILRestriction(silr ,1)

    DrawPlots()
    SuppressQueryOutputOn()

    filename = opts.inputfile+'_Gtheta.txt'
    f=open(filename,'w')
    f.write("#t             load            G            for  R_inner = %e   R_outer = %e  \n"%(opts.R_inn,opts.R_out))
    print 'options ',opts
    for s in range(0,TimeSliderGetNStates(),1):
        SetActiveWindow(1)
        SetTimeSliderState(s)


        Query("Time")
        t=GetQueryOutputValue()

        Query("Weighted Variable Sum")
        G= GetQueryOutputValue()

        print "****** step %i load = %e, G_theta = %e "%(s,t,G)
        f.write("%e \t%e \t%e \n"%(s,t,G))
        f.flush()
        os.fsync(f)
    f.close()
    SetActiveWindow(1)
    DeleteAllPlots()
    DeleteWindow()
    CloseDatabase(MyDatabase)
    return 0

if __name__ == "__main__":
    import sys  
    import os.path

    opts = parse()
    Get_G(opts)
    exit()
