#!/usr/bin/env python
import argparse
import pymef90.mesh
import sys
import os

parser = argparse.ArgumentParser()
parser.add_argument("gmeshFile", help = "The name of the mesh file to be parsed.", type = str)
parser.add_argument("exoFile", help = "The name of the exodus file to be written.", type = str)
parser.add_argument("--force",action="store_true",default=False,help="Overwrite existing files without prompting")
args = parser.parse_args()

if  os.path.exists(args.exoFile):
    if args.force:
        os.remove(args.exoFile)
    else:
        if pymef90.confirm("ExodusII file {0} already exists. Overwrite?".format(args.exoFile)):
            os.remove(args.exoFile)
        else:
            print ('\n\t{0} was NOT generated from {1}\n'.format(args.exoFile,args.gmeshFile))
            sys.exit()
(coord,vertexSet,cellSet,numDim) = pymef90.mesh.GMSHread(args.gmeshFile)
pymef90.mesh.EXODUSwrite(coord,vertexSet,cellSet,numDim,args.exoFile)

