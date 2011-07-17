#
# run as:  visit -cli -s example_visitmef90.py
# Author corrado.maurini@upmc.fr
# 

import pymef90;
from visitmef90 import *

#
D = pymef90.Dictreadtxt('00_INFO.txt')
#
ImageOptions = {}
ImageOptions['x0']=0;
ImageOptions['y0']=0;
ImageOptions['lx']=D['lx'];
ImageOptions['ly']=D['ly'];
ImageOptions['res']=1000;
ImageOptions['lxv']=.9;
ImageOptions['lyv']=.9;
ImageOptions['x0v']=0.1;
ImageOptions['y0v']=0.1;


# Launch makefigure
#MakeFiguresAllSteps()
MakeFiguresLastStep(ImageOptions)
        
# Exit
exit()

