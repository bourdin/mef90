#!/usr/bin/env python
import argparse
import sys
import pymef90
#try:
#    import exomerge3 as exomerge
#except:
#    import exomerge2 as exomerge
import exomerge


def parse(args=None):
    import argparse
    ### Get options from the command line
    parser = argparse.ArgumentParser(description='Remove non-essential_fields from a .gen file.')
    parser.add_argument('inputfile',help='input file')
    parser.add_argument('-o','--outputfile',help='output file',default=None)
    parser.add_argument('--field',help='list of variables to be deleted',nargs='*',default=[])
    parser.add_argument("--force",action="store_true",default=False,help="Overwrite existing files without prompting")
    return parser.parse_args()

def getFields(filename):
    import exodus as exo
    nodalFields = []
    zonalFields = []
    e = exo.exodus(filename) 
    nodalFields = e.get_node_variable_names()
    zonalFields = e.get_element_variable_names()
    e.close()
    return nodalFields,zonalFields

def main():
    options = parse()
    nodalFields,zonalFields = getFields(options.inputfile)

    print('nodal fields: {0}'.format(nodalFields))

    deletedNodalFields = []
    for f in nodalFields[::-1]:
        for ff in options.field:
            if f.lower().startswith(ff.lower()) and f not in deletedNodalFields:
                deletedNodalFields.append(f)

    print('zonal fields: {0}'.format(zonalFields))
    deletedZonalFields = []
    for f in zonalFields[::-1]:
        for ff in options.field:
            if f.lower().startswith(ff.lower()) and f not in deletedZonalFields:
                deletedZonalFields.append(f)

    if len(deletedNodalFields+deletedZonalFields) > 0:
        print('You are about to delete the following fields:')
        for f in deletedNodalFields[::-1]:
            print('{0} (nodal)'.format(f))
        for f in deletedZonalFields[::-1]:
            print('{0} (zonal)'.format(f))
        if not options.force:
            if not pymef90.confirm("Is this OK?"):
                print('Exiting without writing')
                sys.exit()

        if options.outputfile == None:
            outputfile = options.inputfile
        else:
            outputfile = options.outputfile
        
        print('loading model {0}. This may take a while.'.format(options.inputfile))
        model = exomerge.import_model(options.inputfile)

        for f in deletedZonalFields[::-1]:
            model.delete_element_field(f)
            print('\tdeleted {0}'.format(f))
        for f in deletedNodalFields[::-1]:
            model.delete_node_field(f)
            print('\tdeleted {0}'.format(f))

        print('Saving model. This may take a while')
        model.export_model(outputfile)
    else:
        print('Nothing to delete')
    return 0




if __name__ == "__main__":
        sys.exit(main())
