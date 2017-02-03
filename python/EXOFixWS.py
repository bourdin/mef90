import exodus as exo
import sys
import os.path

f = sys.argv[-1]
if not os.path.isfile(f) or not (len(sys.argv) == 2):
    print "usage: {0} <filename>".format(sys.argv[0])
    exit
else:
    fout = f.split(".gen")[0]+"_fixed.gen"
    if os.path.exists(fout):
        print "{0} exists. Delete and run again".format(fout)
    else:
        exo.exodus(f,"r").copy(fout)
 
