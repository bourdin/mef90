#!/usr/bin/env python
import argparse
import sys
import pymef90


def parse(args=None):
    import argparse
    ### Get options from the command line
    parser = argparse.ArgumentParser(description='Remove non-essential_fields from a .gen file.')
    parser.add_argument('inputfile',help='input file')
    parser.add_argument('-o','--outputfile',help='output file',default=None)
    parser.add_argument('--remove',help='list of variables to be deleted',nargs='*',default=[])
    parser.add_argument("--force",action="store_true",default=False,help="Overwrite existing files without prompting")
    return parser.parse_args()

def main():
    import exodus as exo
    import os
    options = parse()


    exoin  = exo.exodus(options.inputfile,mode='r')
    nodalFields = exoin.get_node_variable_names()
    zonalFields = exoin.get_element_variable_names()

    if options.remove == []:
        print('nodal fields: {0}'.format(nodalFields))
        print('zonal fields: {0}'.format(zonalFields))

    keepNodalFields = nodalFields
    for deleteField in options.remove:
        keepNodalFields = [field for field in keepNodalFields if not field.lower().startswith(deleteField.lower())]
    deleteNodalFields = [field for field in nodalFields if not field in keepNodalFields]

    keepZonalFields = zonalFields
    for deleteField in options.remove:
        keepZonalFields = [field for field in keepZonalFields if not field.lower().startswith(deleteField.lower())]
    deleteZonalFields = [field for field in zonalFields if not field in keepZonalFields]


    if len(deleteNodalFields+deleteZonalFields) > 0:
        print('\nDeleting the following fields:')
        for f in deleteNodalFields:
            print('\t{0} (nodal)'.format(f))
        for f in deleteZonalFields:
            print('\t{0} (zonal)'.format(f))

        if options.outputfile == None or options.inputfile == options.outputfile:
            if not options.force:
                print('THIS WILL OVERWRITE INPUT FILE {0}'.format(options.inputfile))
                if not pymef90.confirm("Is this OK?"):
                    print('Exiting without writing')
                    sys.exit()
            outputfile = options.inputfile+'.tmp'
        else:
            if (os.path.exists(options.outputfile)):
                if not options.force:
                    print('THIS WILL OVERWRITE OUTPUT FILE {0}'.format(options.outputfile))
                    if not pymef90.confirm("Is this OK?"):
                        print('Exiting without writing')
                        sys.exit()
                if os.path.exists(options.outputfile):
                    os.remove(options.outputfile)
            outputfile = options.outputfile

        print('inputfile: {0}'.format(options.inputfile))
        print('outputfile: {0}'.format(outputfile))

        print("Copying geometry in {0}".format(outputfile))
        try:
            exoout = exo.copy_mesh(options.inputfile,options.outputfile, array_type='numpy')
        except:
            print('\n\nCopying the background mesh using exodus.copy_mesh failed, trying again using exodus.copy.')
            print('Note that the resulting file may not readable with paraview < 5.8.0 or visit\n\n')
            os.remove(options.outputfile)
            exoin  = exo.exodus(options.inputfile,mode='r')
            exoout = exoin.copy(options.outputfile)

        print("formatting {0}".format(outputfile))
        exoout.set_global_variable_number(0)
        exoout.set_node_variable_number(len(keepNodalFields))
        for i in range(len(keepNodalFields)):
            exoout.put_node_variable_name(keepNodalFields[i],i+1)
        exoout.set_element_variable_number(len(keepZonalFields))
        for i in range(len(keepZonalFields)):
            exoout.put_element_variable_name(keepZonalFields[i],i+1)
        exoout.set_element_variable_truth_table([True] * exoout.numElemBlk.value * len(keepZonalFields))

        step = 1
        for t in exoin.get_times():
            print("processing time step {0}: t={1:.2e}".format(step,t))
            exoout.put_time(step,t)
            for f in keepNodalFields:
                exoout.put_node_variable_values(f,step,exoin.get_node_variable_values(f,step))
            for f in keepZonalFields:
                for cs in exoin.get_elem_blk_ids():
                    exoout.put_element_variable_values(cs,f,step,exoin.get_element_variable_values(cs,f,step))

            step += 1
        exoout.close()

        if options.outputfile == None:
            os.rename(outputfile,options.inputfile)
    else:
        print('Nothing to delete')
    return 0




if __name__ == "__main__":
        sys.exit(main())
