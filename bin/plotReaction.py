#!/usr/bin/env python3
import sys

def parse(args=None):
    import argparse
    ### Get options from the command line
    parser = argparse.ArgumentParser(description='Plot reaction force for vDef.')
    parser.add_argument('inputfile',type=argparse.FileType('r'),nargs='?',help='Input file',default=sys.stdin)
    parser.add_argument('-o','--outputfile',help='output file',default=None)
    parser.add_argument("-d","--debug",action="store_true",default=False,help="Display useless debugging information")
    parser.add_argument("-m","--stepmin",type=int,help="first time step")
    parser.add_argument("-M","--stepmax",type=int,help="last time step")
    parser.add_argument("-e","--Emin",type=float,help="Force min")
    parser.add_argument("-E","--Emax",type=float,help="Force max")
    parser.add_argument("--old",action="store_true",default=False,help="old style energy file (no forces)")
    parser.add_argument("-r","--relative",action="store_true",default=False,help="offset surface energy")
    parser.add_argument("--size",type=float,nargs=2,default=None,help="Figure size")
    parser.add_argument("--title",default=None,help="Figure title")
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
      energies=np.zeros( (energies_old.shape[0],6) )
      energies[:,:3]=energies_old[:,:3]
      energies[:,-2:]=energies_old[:,-2:]
    else:
      energies=np.loadtxt(options.inputfile)
      
    if options.stepmin == None:
      tmin = 0
    else:
      tmin = int(options.stepmin)
    if options.stepmax == None:
      tmax = energies.shape[0]
    else:
      tmax = int(options.stepmax)
    
    if options.debug:
      print('size of energies: {0}'.format(energies.shape))
      print('requested slice: {0}:{1}'.format(tmin,tmax))
      print('Energies: {0}'.format(energies[tmin:tmax,:]))
    
    if options.outputfile != None:
      matplotlib.use('Agg')
      useTex=True
    else:
     useTex = False
    import matplotlib.pyplot as plt

    
    #pymef90.matplotlibdefaults('medium',useTex)
    fig = plt.figure(figsize=options.size)
    
    if options.relative:
      energies[:,-2] -= energies[tmin,-2]
      energies[:,-1] -= energies[tmin,-1]
    ### plot
    plt.plot(energies[tmin:tmax,1],energies[tmin:tmax,2]/energies[tmin:tmax,1]/2,'-')
    plt.grid()
    plt.xlabel('t')
    plt.ylabel('Force')
    if options.title:
        plt.title(options.title)
    else:
        plt.title('Reaction force vs normalized time')
    pymef90.setspines()

    BB = plt.axis()
    bb = [b for b in BB]
    if options.Emin:
        bb[2] = options.Emin      
    if options.Emax:
        bb[3] = options.Emax      
    plt.axis(bb)

    ### export plot if needed
    if options.outputfile != None:
      fig.tight_layout(pad=0.1)
      plt.savefig(options.outputfile)
    else:
      plt.show()

if __name__ == "__main__":
        sys.exit(main())
