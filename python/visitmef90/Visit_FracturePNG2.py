###############################################################################
# Visit Script to generate and save PNG images of a mef90-sieve output
#
# How to use:
#  
#  run VisIt IN THE DIRECTORY WHERE THE FILE IS LOCATED with Command Line Interface and No Windows:
#    "visit -cli -nowin -s $MEF_HOME/python/visitmef90/Visit_FracturePNG2.py"
#
# Note: the script identify the database name from the .e2c filename 
#
# Author: Corrado Maurini: cmaurini@gmail.com
###############################################################################

import os; import sys; import glob;
# --------------------------------------------------
# Step 0; Files and directories
# --------------------------------------------------
## Get the working directory
Directory=os.environ.get('PWD') + '/'
## Get the prefix for the file name
#############
energyfiles=glob.glob('*.e2c')
if len(energyfiles)==0:
  print "It seems there are not files to process"
else:
  for energyfilename in energyfiles:
      Prefix=energyfilename.split('.e2c')[0]
      ##  
      ## Step 1: Open the database
      ##
      if os.path.exists(Prefix+'-0001.gen'):
        MyDatabase= Prefix+'-*.gen database'
      else:
        MyDatabase= Prefix+'-0000.gen'
      OpenDatabase(MyDatabase,0)
      # --------------------------------------------------
      # Step 2: Add plots and set their properties
      # --------------------------------------------------
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
      AnnotationAtts = AnnotationAttributes()
      AnnotationAtts.axes2D.visible = 0
      AnnotationAtts.userInfoFlag = 0
      AnnotationAtts.databaseInfoFlag = 0
      AnnotationAtts.legendInfoFlag = 0
      AnnotationAtts.axesArray.visible = 0
      SetAnnotationAttributes(AnnotationAtts)
      # --------------------------------------------------
      # Step 3: Set Save file parameters 
      # --------------------------------------------------
      s = SaveWindowAttributes()
      s.format = s.PNG
      s.saveTiled = 0
      s.family=0
      s.width = 600
      s.height = 200
      s.screenCapture = 0
      s.progressive = 1
      forceMerge = 1
      s.quality = 80
      #
      # Save images of all time steps and add each image filename to a list.
      #
      names = []
      for state in range(0,TimeSliderGetNStates()):
         SetTimeSliderState(state)
         # Save the image
         s.fileName = Directory+ Prefix+'-%04i'%state
         SetSaveWindowAttributes(s)
         n = SaveWindow()
         names = names + [n]
         print "Processing image for frame %s" %state
         #print "Processing image for frame %s" %n
      #print names
exit()
