#!/usr/bin/env python
import argparse
import sys
import pymef90
try:
    import exomerge3 as exomerge
except:
    import exomerge2 as exomerge



def parse(args=None):
    import argparse
    ### Get options from the command line
    parser = argparse.ArgumentParser(description='Remove non-essential_fields from a .gen file.')
    parser.add_argument('inputfile',help='input file')
    parser.add_argument('-o','--outputfile',help='output file',default=None)
    parser.add_argument('--field',help='list of variables to be deleted',nargs='*',default=[])
    parser.add_argument("--force",action="store_true",default=False,help="Overwrite existing files without prompting")
    return parser.parse_args()

def main():
    options = parse()
    model = exomerge.import_model(options.inputfile)

    nodalFields = model.get_node_field_names()
    print('nodal fields: {0}'.format(nodalFields))
    deletedNodalFields = []
    for f in nodalFields[::-1]:
        for ff in options.field:
            if f.lower().startswith(ff.lower()):
                model.delete_node_field(f)
                deletedNodalFields.append(f)

    zonalFields = model.get_element_field_names()
    print('zonal fields: {0}'.format(zonalFields))
    deletedZonalFields = []
    for f in zonalFields[::-1]:
        for ff in options.field:
            if f.lower().startswith(ff.lower()):
                model.delete_element_field(f)
                deletedZonalFields.append(f)


    if len(deletedNodalFields+deletedZonalFields) > 0:
        if not options.force:
            print('You are about to delete the following fields:')
            for f in deletedNodalFields[::-1]:
                print('{0} (nodal)'.format(f))
            for f in deletedZonalFields[::-1]:
                print('{0} (zonal)'.format(f))
            if not pymef90.confirm("Is this OK?"):
                print('Exiting without writing')
                sys.exit()

        if options.outputfile == None:
            outputfile = options.inputfile
        else:
            outputfile = options.outputfile
        model.export_model(outputfile)
    else:
        print('Nothing to delete')
    return 0




if __name__ == "__main__":
        sys.exit(main())
