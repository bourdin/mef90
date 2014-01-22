#!/usr/bin/env python
import sys

def parse(args=None):
    import argparse
    ### Get options from the command line
    parser = argparse.ArgumentParser(description='Compute boundary displacement for a surfing computation')
    parser.add_argument('-i','--inputfile',help='input file',default=None)
    parser.add_argument('-o','--outputfile',help='output file',default=None)
    parser.add_argument("--time_min",type=float,help="first time step",default=0.)
    parser.add_argument("--time_max",type=float,help="last time step",default=1.)
    parser.add_argument("--time_numstep",type=int,help="number of time step",default=11)
    parser.add_argument("--velocity",type=float,help="surfing velocity",default=1.)
    parser.add_argument("--surfingcs",type=int,nargs='?',help="list of cell sets where surfing boundary displacements are applied",default=[])
    parser.add_argument("--surfingvs",type=int,nargs='?',help="list of vertex sets where surfing boundary displacements are applied",default=[])
    return parser.parse_args()
    
def exoformat(e):
    global_variable_name = ["Elastic Energy","Work","Surface Energy","Total Energy"]
    if e.num_dimensions() == 2: 
        node_variable_name  = ["Temperature","Damage","Displacement_X","Displacement_Y"]
        element_variable_name   = ["empty","empty","Pressure_Force","Force_X","Force_Y"]
    else:
        node_variable_name  = ["Temperature","Damage","Displacement_X","Displacement_Y","Displacement_Z"]
        element_variable_name   = ["empty","empty","Pressure_Force","Force_X","Force_Y","Force_Y"]
    e.set_global_variable_number(0)
    e.set_node_variable_number(len(node_variable_name))
    for i in range(len(node_variable_name)):
        e.put_node_variable_name(node_variable_name[i],i+1)
    e.set_element_variable_number(len(element_variable_name))
    return(0)

def cart2polar(x, y):
    import numpy as np
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return r, theta
    
def surfingBC(e,t,xc,yc,cslist,vslist,E,nu):
    import exodus as exo
    import numpy as np
    
    dim = e.num_dimensions()
    X,Y,Z=e.get_coords()
    #Ux = X
    #Uy = Y
    #Uz = Z
    #for i in range(len(Ux)):
    #    Ux[i] = 0.
    #    Uy[i] = 0.
    #    Uz[i] = 0.
    U = np.zeros([3,len(X)],dtype=exo.c_double)
        
    for set in vslist:
        for v in e.get_node_set_nodes(set):
            r,theta = cart2polar(X[v-1]-xc-t,Y[v-1]-yc)
            z = Z[v-1]
            #Ux[v-1] = r
            #Uy[v-1] = theta
            #Uz[v-1] = z
            U[0,v-1] = r
            U[1,v-1] = theta 
            U[2,v-1] = z
    return U

def main():
    import exodus as exo
    import numpy as np
    options = parse()
    
    print options.inputfile
    
    exoin  = exo.exodus(options.inputfile,mode='r')
    exoout = exoin.copy(options.outputfile)
    exoin.close()
    exoformat(exoout)
    
    E  = 1
    nu = .3
    step = 0
    for t in np.linspace(options.time_min,options.time_max,options.time_numstep):
        print "writing step",step+1,t
        exoout.put_time(step+1,t)
        U = surfingBC(exoout,t,0,0.,exoout.get_elem_blk_ids,exoout.get_node_set_ids(),E,nu)
        X,Y,Z=exoout.get_coords()
        #print U[0,:],len(U[0,:])
        exoout.put_node_variable_values("Displacement_X",step+1,U[0,:])
        exoout.put_node_variable_values("Displacement_Y",step+1,U[1,:])
        #exoout.put_node_variable_values("Displacement_X",step,X)
        #exoout.put_node_variable_values("Displacement_Y",step,Y)
        step += 1
    exoout.close()
    ### 
    ### compute boundary displacement at vertex sets
    ###
    
if __name__ == "__main__":
        sys.exit(main())

