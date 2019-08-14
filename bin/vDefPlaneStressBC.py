#!/usr/bin/env python
import sys

def parse(args=None):
    import argparse
    ### Get options from the command line
    parser = argparse.ArgumentParser(description='Compute boundary displacement for a surfing computation')
    parser.add_argument('-i','--inputfile',help='input file',default=None)
    parser.add_argument('-o','--outputfile',help='output file',default=None)
    parser.add_argument('--E',type=float,help='Youngs Modulus',default=1)
    parser.add_argument('--nu',type=float,help='Poisson Ratio',default=0)    
    parser.add_argument('--sigma',type=float,help='Applied stress (sigma_11, sigma_22, sigma_12)',nargs=3,default=[1,0,0])
    parser.add_argument("--cs",type=int,nargs='*',help="list of cell sets where the beam is applied",default=[])
    parser.add_argument("--vs",type=int,nargs='*',help="list of vertex sets where the beam is applied",default=[])
    parser.add_argument("--force",action="store_true",default=False,help="Overwrite existing files without prompting")
    parser.add_argument("--time_min",type=float,default=1.,help='Start time')
    parser.add_argument("--time_max",type=float,default=1.,help='End time')
    parser.add_argument("--time_numstep",type=int,default=1,help='Number of time steps')
    return parser.parse_args()
    
def exoformat(e,plasticity=False):
    if plasticity:
        global_variable_name = ["Elastic Energy","Work","Surface Energy","Total Energy","Dissipation Plastic"]
        if e.num_dimensions() == 2: 
            node_variable_name  = ["Temperature","Damage","Displacement_X","Displacement_Y"]
            element_variable_name   = ["External_Temperature","Heat_Flux","Pressure_Force",
                                       "Force_X","Force_Y",
                                       "Stress_XX","Stress_YY","Stress_XY",
                                       "Cumulated_Plastic_Energy","plasticStrain_XX","plasticStrain_YY","plasticStrain_XY"]
        else:
            node_variable_name  = ["Temperature","Damage","Displacement_X","Displacement_Y","Displacement_Z"]
            element_variable_name   = ["External_Temperature","Heat_Flux","Pressure_Force",
                                       "Force_X","Force_Y","Force_Z",
                                       "Stress_XX","Stress_YY","Stress_ZZ","Stress_YZ","Stress_XZ","Stress_XY",
                                       "Cumulated_Plastic_Energy","plasticStrain_XX","plasticStrain_YY","plasticStrain_ZZ","plasticStrain_XY","plasticStrain_YZ","plasticStrain_XZ","plasticStrain_XY"]
    else:
        global_variable_name = ["Elastic Energy","Work","Surface Energy","Total Energy"]
        if e.num_dimensions() == 2: 
            node_variable_name  = ["Temperature","Damage","Displacement_X","Displacement_Y"]
            element_variable_name   = ["External_Temperature","Heat_Flux","Pressure_Force",
                                       "Force_X","Force_Y",
                                       "Stress_XX","Stress_YY","Stress_XY"]
        else:
            node_variable_name  = ["Temperature","Damage","Displacement_X","Displacement_Y","Displacement_Z"]
            element_variable_name   = ["External_Temperature","Heat_Flux","Pressure_Force",
                                       "Force_X","Force_Y","Force_Z",
                                       "Stress_XX","Stress_YY","Stress_ZZ","Stress_YZ","Stress_XZ","Stress_XY"]
    e.set_global_variable_number(0)
    e.set_node_variable_number(len(node_variable_name))
    for i in range(len(node_variable_name)):
        e.put_node_variable_name(node_variable_name[i],i+1)
    e.set_element_variable_number(len(element_variable_name))
    for i in range(len(element_variable_name)):
        e.put_element_variable_name(element_variable_name[i],i+1)
    e.set_element_variable_truth_table([True] * e.numElemBlk.value * len(element_variable_name))
    return(0)

def displacementBC(e,t,options):
    import exodus as exo
    import numpy as np
    
    dim = e.num_dimensions()
    if not dim == 2:
        print('This function only makes sense in 2D')
        exit(-1)

    e11 = (options.sigma[0] - options.nu * options.sigma[1]) / options.E
    e22 = (options.sigma[1] - options.nu * options.sigma[0]) / options.E
    e12 = options.sigma[2] * (1. + options.nu) / options.E

    X,Y,Z=e.get_coords()
    U = np.zeros([2,len(X)],dtype=exo.c_double)
    
    csoffset = [e.elem_blk_info(set)[1] for set in options.cs]        
    for set in options.cs:
        connect = e.get_elem_connectivity(set)
        for cid in range(connect[1]):
            vertices = [connect[0][cid*connect[2]+c] for c in range(connect[2])]
            for v in vertices:
                U[0,v-1] = t * (e11 * X[v-1] + e12 * Y[v-1])
                U[1,v-1] = t * (e12 * X[v-1] + e22 * Y[v-1])
        
    for set in options.vs:
        for v in e.get_node_set_nodes(set):
            U[0,v-1] = t * (e11 * X[v-1] + e12 * Y[v-1])
            U[1,v-1] = t * (e12 * X[v-1] + e22 * Y[v-1])
    return U


def main():
    import exodus as exo
    import numpy as np
    import os
    import pymef90
    options = parse()
    
    if  os.path.exists(options.outputfile):
        if options.force:
            os.remove(options.outputfile)
        else:
            if pymef90.confirm("ExodusII file {0} already exists. Overwrite?".format(options.outputfile)):
                os.remove(options.outputfile)
            else:
                print ('\n\t{0} was NOT generated.\n'.format(options.outputfile))
                return -1
    exoin  = exo.exodus(options.inputfile,mode='r')
    exoout = exoin.copy(options.outputfile)
    exoin.close()
    exoformat(exoout)
    
    if not  exoout.num_dimensions() == 2:
        print("This program only makes sense in 2D")
        return (-1)

    T = np.linspace(options.time_min,options.time_max,options.time_numstep)
    for step in range(options.time_numstep):
        t = T[step]
        print "writing step",step+1,t
        exoout.put_time(step+1,t)
        U = displacementBC(exoout,t,options)
        exoout.put_node_variable_values("Displacement_X",step+1,U[0,:])
        exoout.put_node_variable_values("Displacement_Y",step+1,U[1,:])
    exoout.close()
    return (0)
    
if __name__ == "__main__":
        sys.exit(main())
