#
# run as:  visit -cli -s example_visitmef90.py
# Author corrado.maurini@upmc.fr
# 
from visit import *
import os; import sys; import glob; 
import pymef90;
import visitmef90 as vm

# gen prefix
prefix='Slab2DNG-QTri-X100-Y60-h3'

# get last step
nsteps = pymef90.getlaststep(prefix+'.ener')

# Represt Fracture 
vm.FigureFracture(prefix,nsteps-1) 

# Save image of a single time step 
vm.ExportSingleFigure(prefix,nsteps-1) 

# Save images of all time steps 
vm.ExportTimeFigure(prefix)    

# Exit
exit()

