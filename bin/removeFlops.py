#!/usr/bin/env python3
def main():
    import sys
    import shutil
    if not len(sys.argv) == 2:
        print("usage: {} file".format(sys.argv[0]))
        return -1

    shutil.copy2(sys.argv[1],sys.argv[1]+'~')
    file  = open(sys.argv[1],'r')

    output = ''
    for l in file.readlines():
        if not (('petsclogdouble' in l.lower()) or ('flops' in l.lower()) or ('petsclogflops' in l.lower())):
            output+=l
    file.close()

    file = open(sys.argv[1],'w')
    file.write(output)
    file.close()
    
if __name__ == "__main__":
    import sys
    sys.exit(main())

