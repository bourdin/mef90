#!/usr/bin/env python
import sys

def parse(args=None):
    import argparse
    ### Get options from the command line
    parser = argparse.ArgumentParser(description='Substitute values from <json> file into input file')
    parser.add_argument('inputfile',type=argparse.FileType('r'),nargs='?',help='Input file',default=sys.stdin)
    parser.add_argument('-o','--outputfile',type=argparse.FileType('w'),help='output file',default=sys.stdout)
    parser.add_argument('-j','--json',type=argparse.FileType('r'),default=None,help='JSON file containig all values to substitute')
    return parser.parse_args()

def main():
    import json
    options = parse()

    D = json.load(options.json)
    options.outputfile.write(options.inputfile.read().format(**D))
    options.outputfile.close()

if __name__ == "__main__":
        sys.exit(main())
