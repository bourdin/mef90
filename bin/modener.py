#!/usr/bin/env python
def parse(args=None):
    import argparse
    import sys
    ### Get options from the command line
    parser = argparse.ArgumentParser(description='Plot energy evolution for VarFracQS.')
    parser.add_argument("-f","--forces",default=False,action="store_true",help="displays forces work")
    parser.add_argument('inputfile',nargs='?',type=argparse.FileType('r'),help='Input file',default=sys.stdin)
    parser.add_argument('outputfile',nargs='?',type=argparse.FileType('r'),help='(optional) output file',default=sys.stdout)
    parser.add_argument("-d","--debug",action="store_true",default=False,help="Display useless debugging information")
    parser.add_argument("-M","--stepmax",type=int,help="last time step")
    parser.add_argument("-o","--old",action="store_true",default=False,help="old style energy file (no forces)")
    parser.add_argument("-r","--relative",action="store_true",default=False,help="offset surface energy")
    return parser.parse_args()


def main():
    import pymef90
    import matplotlib
    import numpy as np
    options = parse()

    if options.debug:
      print("Option table: {0}".format(options))
    if options.old:
      energies_old=np.loadtxt(options.inputfile)
      energies=np.zeros( (energies_old.shape[0], 6) )
      energies[:,:3]=energies_old[:,:3]
      energies[:,-2:]=energies_old[:,-2:]
    else:
      energies=np.loadtxt(options.inputfile)

    if options.debug:
        print('size of energies: {0}'.format(energies.shape))
        print('Energies: {0}'.format(energies))

    energiesBT=pymef90.energies.fixBT(energies, options.stepmax)
    
    pymef90.energies.save(options.outputfile, energiesBT)

if __name__ == "__main__":
    import sys
    sys.exit(main())
