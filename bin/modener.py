#!/usr/bin/env python
import pymef90
import numpy as np
from optparse import OptionParser

### Get options from the command line
parser = OptionParser()
parser.add_option("-d", "--debug", dest="debug", action="store_true", default=False, help="display useless debugging information")
parser.add_option("-i", "--input", dest="inputfile", help="energy file")
parser.add_option("-o", "--output", dest="outputfile", help="output file name")
parser.add_option("--old", dest="old", action="store_true", default=False, help="old style energy file (no forces)")
parser.add_option("-M", "--stepmax", dest="stepmax", default=None, help="step to parse until" )


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

if options.debug:
  print('size of energies: {0}'.format(energies.shape))
  print('Energies: {0}'.format(energies))

energiesBT=pymef90.energies.fixBT(energies, int(options.stepmax))

pymef90.energies.save(options.outputfile, energiesBT)
