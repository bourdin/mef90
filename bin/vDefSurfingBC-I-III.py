#!/usr/bin/env python
import sys
if sys.version_info.major == 3:
    import exodus3 as exo
else:
    import exodus2 as exo

def parse(args=None):
    import argparse
    ### Get options from the command line
    parser = argparse.ArgumentParser(description='Compute boundary displacement for a surfing computation')
    parser.add_argument('-i','--inputfile',help='input file',default=None)
    parser.add_argument('-o','--outputfile',help='output file',default=None)
    parser.add_argument("--time_min",type=float,help="first time step",default=0.)
    parser.add_argument("--time_max",type=float,help="last time step",default=1.)
    parser.add_argument("--time_numstep",type=int,help="number of time step",default=11)
    parser.add_argument("--E",type=float,help="Youngs modulus",default=1.)
    parser.add_argument("--nu",type=float,help="Poisson Ratio",default=.3)
    parser.add_argument("--ampl",type=float,help="Amplification",default=1.)  
    parser.add_argument("--mixity",type=float,help="theta = K_I/K_III: mode mixity",default=0.)        
    parser.add_argument("--initialpos",type=float,nargs=3,help="initial crack tip postion",default=[0.,0,0])
    parser.add_argument("--cs",type=int,nargs='*',help="list of cell sets where surfing boundary displacements are applied",default=[])
    parser.add_argument("--vs",type=int,nargs='*',help="list of vertex sets where surfing boundary displacements are applied",default=[])
    parser.add_argument("--force",action="store_true",default=False,help="Overwrite existing files without prompting")
    return parser.parse_args()
    
def exoformat(e):
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

def cart2polar(x, y):
    import numpy as np
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return r, theta
        
def surfingBC(e,t,Xc,cslist,vslist,E,nu,ampl,mixity):
    import numpy as np
    
    kappa = (3.0-nu)/(1.0+nu)
    mu = E / (1. + nu) * .5

    dim = e.num_dimensions()
    X,Y,Z=e.get_coords()
    U = np.zeros([3,len(X)])
    
    csoffset = [e.elem_blk_info(set)[1] for set in cslist]        
    for set in cslist:
        connect = e.get_elem_connectivity(set)
        for cid in range(connect[1]):
            vertices = [connect[0][cid*connect[2]+c] for c in range(connect[2])]
            for v in vertices:
                r,theta = cart2polar(X[v-1]-Xc[0]-t,Y[v-1]-Xc[1])
                z = Z[v-1]-Xc[2]
                U[0,v-1] = mixity * ampl * np.sqrt(r / np.pi * .5) / mu * .5 * np.cos(theta * .5) * (kappa - np.cos(theta))
                U[1,v-1] = mixity * ampl * np.sqrt(r / np.pi * .5) / mu * .5 * np.sin(theta * .5) * (kappa - np.cos(theta))
                U[2,v-1] = np.sqrt(2.*r/np.pi) * np.sin(theta * .5) / mu 
        
    for set in vslist:
        for v in e.get_node_set_nodes(set):
            r,theta = cart2polar(X[v-1]-xc-t,Y[v-1]-yc)
            z = Z[v-1]
            U[0,v-1] = mixity * ampl * np.sqrt(r / np.pi * .5) / mu * .5 * np.cos(theta * .5) * (kappa - np.cos(theta))
            U[1,v-1] = mixity * ampl * np.sqrt(r / np.pi * .5) / mu * .5 * np.sin(theta * .5) * (kappa - np.cos(theta))
            U[2,v-1] = np.sqrt(2.*r/np.pi) * np.sin(theta * .5) / mu
    return U

def main():
    import numpy as np
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
    exoout.close()
    exoout  = exo.exodus(options.outputfile,mode='a',array_type='numpy')
    ### Adding a QA record, needed until visit fixes its exodus reader
    import datetime
    import os.path
    import sys
    QA_rec_len = 32
    QA = [os.path.basename(sys.argv[0]),os.path.basename(__file__),datetime.date.today().strftime('%Y%m%d'),datetime.datetime.now().strftime("%H:%M:%S")]
    exoout.put_qa_records([[ q[0:31] for q in QA],])

    exoformat(exoout)
    
    dim = exoout.num_dimensions()
    step = 0
    for t in np.linspace(options.time_min,options.time_max,options.time_numstep):
        print "writing step",step+1,t
        exoout.put_time(step+1,t)
        U = surfingBC(exoout,t,options.initialpos,options.cs,options.vs,options.E,options.nu,options.ampl,options.mixity)
        X,Y,Z=exoout.get_coords()
        exoout.put_node_variable_values("Displacement_X",step+1,U[0,:])
        exoout.put_node_variable_values("Displacement_Y",step+1,U[1,:])
        if dim == 3:
            exoout.put_node_variable_values("Displacement_Z",step+1,U[2,:])
        step += 1
    exoout.close()
    ### 
    ### compute boundary displacement at vertex sets
    ###
    
if __name__ == "__main__":
        sys.exit(main())

