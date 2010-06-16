#!/usr/bin/env python
import pymef90
import matplotlib
import numpy as np
from optparse import OptionParser

### Get options from the command line
parser = OptionParser()
parser.add_option("-d", "--debug", dest="debug", action="store_true", default=False, help="display useless debugging information")
parser.add_option("-i", "--input", dest="inputfile", help="energy file")
parser.add_option("-f", "--forces", dest="forces", action="store_true", help="displays forces work", default=False)
parser.add_option("-o", "--output", dest="outputfile", help="output file name")
parser.add_option("-m", "--stepmin", dest="stepmin", help="first time step")
parser.add_option("-M", "--stepmax", dest="stepmax", help="last time step")
parser.add_option("--old", dest="old", action="store_true", default=False, help="old style energy file (no forces)")
parser.add_option("-r", "--relative", dest="relative", action="store_true", default=False, help="offset surface energy")

(options, args) = parser.parse_args()
if options.debug:
  print("Option table: {0}".format(options))
if options.inputfile == None:
  parser.error("must specify input file with -i or --input")

### Read input file
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
  print('Energies: {0}'.format(energies[tmin:tmax,:]))


if options.outputfile != None:
  matplotlib.use('Agg')
import matplotlib.pyplot as plt

if options.relative:
  energies[:,4] -= energies[tmin,4]
### plot
pymef90.energies.plot(energies[tmin:tmax,:], options.forces)
  
### export plot if needed
if options.outputfile != None:
  plt.savefig(options.outputfile)
else:
  plt.show()
