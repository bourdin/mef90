import fnmatch
import os

from . import MEF90Action, scan, check_mtime
class all(MEF90Action):
    def action(self, args):
        import os.path
        data = scan(args.outdir,args.paths,args.excludedirs,
        args.excludefiles,args.project)
        
        scriptpath = os.path.dirname(__file__)
        energyplotscript = os.path.join(os.getenv('MEF90_DIR'),'bin','plotener.py')
        pngscript = 'visit -cli -nowin -s ' + os.path.join(os.getenv('MEF90_DIR'),'bin','PNGplot.py')
        movtransientscript = 'visit -cli -nowin -s ' + os.path.join(os.getenv('MEF90_DIR'),'bin','MOVTransient.py')
        mov3Dscript = 'visit -cli -nowin -s ' + os.path.join(os.getenv('MEF90_DIR'),'bin','MOV3D.py')
        jobs = data['jobs']
        print 'Project: {0}'.format(data['project'])
    
        for jobid in sorted(jobs.keys()):
            job = jobs[jobid]
            print '{0:_^40}'.format(jobid)
            try:
                prefix = os.path.join(job['path'],job['info']['JOBID'])
                print prefix
            except KeyError:
                try:
                    prefix = os.path.join(job['path'],job['info']['prefix'])
                    print prefix
                except KeyError:
                    try: 
                        prefix = os.path.join(job['path'],job['info']['PBS_JOBID'])
                        print prefix
                    except KeyError: 
                            print "Could not find computation base name. Does 00_INFO.txt contain a JOBID or prefix entry?"
                            break

            if args.plotener: 
                infile = prefix+'.ener'
                outfile = prefix+'_ener.pdf'
                print infile, outfile
                if not check_mtime(infile, outfile):
                    print "Generating energy plot   %s"%outfile
                    cmd = 'python %s %s -o %s'%(energyplotscript,infile,outfile)
                    os.system(cmd)

            if args.plotpng or args.all: 
                infile = prefix+'.ener'
                outfile = prefix+'.png'
                if not check_mtime(infile, outfile):
                    print "Generating png plot %s"%outfile
                    cmd = 'cd %s; %s'%(job['path'],pngscript)
                    #print cmd
                    os.system(cmd)

            if args.movtransient: 
                infile = prefix+'.ener'
                outfile = prefix+'-Transient.avi'
                if not check_mtime(infile, outfile):
                    print "Generating transient movie %s"%outfile
                    cmd = 'cd %s; %s'%(job['path'],movtransientscript)
                    print cmd
                    #os.system(cmd)

            if args.mov3d: 
                infile = prefix+'.ener'
                outfile = prefix+'-3D.avi'
                if not check_mtime(infile, outfile):
                    print "Generating 3D movie %s (if 3D)"%outfile
                    cmd = 'cd %s; %s'%(job['path'],mov3Dscript)
                    print cmd
                    #os.system(cmd)
