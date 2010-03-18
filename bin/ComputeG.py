#!/usr/bin/env python 

def computeG(t, Eel, Es):
    import numpy as np

### Compute G/Gc one more time
    G = -(Eel[1:]/pow(t[1:],2)-Eel[:-1]/pow(t[:-1],2)) / (Es[1:]-Es[:-1])
    #G = -(Eel[1:]-Eel[:-1]) / (Es[1:]-Es[:-1])

    return G
    
def plotG(G, l):
    import matplotlib.pyplot as plt

    plt.plot(Es[1:], G, 'r*')
    plt.grid()
    plt.xlabel('Surface energy')
    plt.ylabel('$G/G_c(1,l)$')
    plt.title('Energy release rate')

    return 0

def plotGriffith(G,t):
    import matplotlib.pyplot as plt

    plt.plot(t[1:], pow(t[1:],2)*G, 'b*')
    plt.title('Griffith criterion')
    plt.xlabel('$t$')
    plt.ylabel('$t^2G/G_c$')
    plt.grid()


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np
    from optparse import OptionParser
  
### Get options from the command line
    parser = OptionParser()
    parser.add_option("-d", "--debug", dest="debug", action="store_true", default=False, help="display useless debugging information")
    parser.add_option("-i", "--input", dest="inputfile", help="energy file")
    parser.add_option("-o", "--output", dest="outputfile", help="output file name")
    parser.add_option("-m", "--stepmin", dest="stepmin", help="first time step")
    parser.add_option("-M", "--stepmax", dest="stepmax", help="last time step")
    parser.add_option("--old", dest="old", action="store_true", default=False, help="old style energy file (no forces)")
    (options, args) = parser.parse_args()
    if options.debug:
        print(options)
    if options.inputfile == None:
        parser.error("must specify input file with -i or --input")

### Read input file and crop the energies if requested
    if options.old:
        energies_old=np.loadtxt(options.inputfile)
        energies=np.zeros( (energies_old.shape[0], 6) )
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

    t   = energies[tmin:tmax,1]
    Eel = energies[tmin:tmax,2]
    Es  = energies[tmin:tmax,4]

    G = computeG(t, Eel, Es)
    #plt.subplot(211)
    plotG(G, Eel)
    
    #plt.subplot(212)
    #plotGriffith(G,t)

    plt.show()
