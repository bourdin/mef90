#
# run as:  visit -cli -s example_visitmef90.py
# Author corrado.maurini@upmc.fr
# 
from visit import *
import os; import sys; import glob; 
import pymef90;
from visitmef90 import *

# Launch makefigure
MakeFiguresAllSteps()
MakeFiguresLastStep()
        
# Exit
exit()

