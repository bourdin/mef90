#!/usr/bin/env python
import sys
if sys.version_info.major == 3:
    import exodus3 as exo
else:
    import exodus2 as exo

def parse(args=None):
    import argparse
    ### Get options from the command line
    parser = argparse.ArgumentParser(description='Creates an inputfile for TestHeatXferCtx.F90')
    parser.add_argument('-i','--inputfile',help='input file',default=None)
    parser.add_argument('-o','--outputfile',help='output file',default=None)
    parser.add_argument("--time_min",type=float,help="first time step",default=0.)
    parser.add_argument("--time_max",type=float,help="last time step",default=1.)
    parser.add_argument("--time_numstep",type=int,help="number of time step",default=11)
    parser.add_argument("--force",action="store_true",default=False,help="Overwrite existing files without prompting")
    return parser.parse_args()
    
def exoformat(e):
    node_variable_name  = ["Temperature",]
    element_variable_name   = ["Flux",]
    e.set_global_variable_number(0)
    e.set_node_variable_number(len(node_variable_name))
    for i in range(len(node_variable_name)):
        e.put_node_variable_name(node_variable_name[i],i+1)
    e.set_element_variable_number(len(element_variable_name))
    for i in range(len(element_variable_name)):
        e.put_element_variable_name(element_variable_name[i],i+1)
    e.set_element_variable_truth_table([True] * e.numElemBlk.value * len(element_variable_name))
    return(0)

def main():
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
    
    dim          = exoout.num_dimensions()
    # numCells     = exoout.num_elems()
    # numVertices  = exoout.num_nodes()
    X,Y,Z        = exoout.get_coords()
    # Temperature  = np.zeros(numVertices)
    # Flux         = np.zeros(num_cells)

    exoout.put_node_variable_values("Temperature",1,X)
    exoout.put_node_variable_values("Temperature",2,Y)
    if (dim == 3):
        exoout.put_node_variable_values("Temperature",3,Z)

    exoout.close()
    
if __name__ == "__main__":
        sys.exit(main())

