###############################################################################
# Process all the jobs in all the subdirectoris to extract PNG images with visit
#
# How to use : It assue that the visit script to process the images is the current directory. If not the path must be specified (see the variable visitscript)
#
# Author: Corrado Maurini: cmaurini@gmail.com
###############################################################################
import os; import sys; import glob;
# get root directory
rootdir=os.getcwd()
filelist=os.listdir('.')
visitscript_default=rootdir+'/'+'Visit_FracturePNG2.py'
### Parse input
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--s", dest='visitscript',type='string',default=visitscript_default,help="path to visit script to apply in each directory")
(options, args) = parser.parse_args()
print options
print args
# Start applying the script in all the subdirectories
for DIR in filelist:
    if os.path.isdir(DIR)==True:
        print DIR
        os.chdir(rootdir+'/'+DIR) 
        # Postprocessing with VisIT
        print "----------------------------------------------------------"
        print "    Start postprocessing of files in the directory        "
        print "    %s "%DIR
        print "----------------------------------------------------------"
        command='visit -cli -nowin  -s '+ options.visitscript
        os.system(command)
        print "----------------------------------------------------------"
        print "-   Done postprocessing of files in the current directory "
        print "----------------------------------------------------------"
        os.chdir(rootdir) 
       "----------------------------------------------------------"
print ':::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
print ': Done postprocessing all directories '
print ':::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
   