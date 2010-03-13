#        Visit Script to save PNG images of a mef90-sieve output
#
#                  How to use
#
# 1) set the environmental variable FILE_PREFIX to set the file name:
#    e.g., for the file set ShrinkingSlab2D-X1-Y2-nc1-lc0.05-h0.025-XXXX.gen
#    "export FILE_PREFIX=ShrinkingSlab2D-X1-Y2-nc1-lc0.05-h0.025"
#    
# 2) run VisIt with Command Line Interface and No Windows:
#    "visit -cli -nowin -s Visit_FracturePNG.py"
#
# Note: the filename is passed as environmental variable. This should be improved
#
#
import os;
import sys
# Set how many images to write 
Nimages=20
#
# Step 0; Files and directories
#
## Get the working directory
#Directory=os.environ.get('WORKDIR')+'/'
Directory=os.environ.get('PWD') + '/'
## Get the prefix for the file name
Prefix=os.environ.get('FILE_PREFIX')
#Prefix='ShrinkingSlab2D-X1-Y1-nc1-lc0.1-h0.025'
## Create directory to store images
ImageDir=Directory+'Images_'+Prefix+'/'
if not os.path.isdir(ImageDir):
  os.mkdir(ImageDir)
#  
# Step 1: Open the database
#
MyDatabase= Prefix+'-*.gen database'
OpenDatabase(MyDatabase)
#
# Step 2: Add plots and set their properties
#
AddPlot('Pseudocolor', 'Fracture')
#AddPlot('Mesh','mesh')
p = PseudocolorAttributes()
# Set the min/max values
p.minFlag = 1
p.maxFlag = 1
p.min=0.0
p.max=1.0
p.legendFlag=0
#p.centering='Nodal'
#p.colorTableName='hot'
SetPlotOptions(p)
# Step 3: Draw the plots
DrawPlots()
# Set the view
#v = GetView2D()
#windowCoords = (-0.776642, 1.70354, -0.425681, 0.626574)
#viewportCoords = (0, 2, 0, 1)
#SetView2D(v)
#
# Step 3: Set Save file parameters 
#
s = SaveWindowAttributes()
s.format = s.PNG
s.fileName = ImageDir+Prefix+'-'
s.saveTiled=0
s.width = 600
s.height = 200
s.screenCapture = 0
s.progressive = 1
s.quality = 80
SetSaveWindowAttributes(s)
#
# Save images of all time steps and add each image filename to a list.
#
names = []
for state in range(1,TimeSliderGetNStates(),TimeSliderGetNStates()/Nimages):
  SetTimeSliderState(state)
  # Save the image
  n = SaveWindow()
  names = names + [n]
print names
exit()
