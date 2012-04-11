#!/usr/bin/env python
# encoding: utf-8
#
USAGE = "python run_bar.py  --dim 2 --tmax 2 --epsilon .2 --meshsize .05 --atnum 1 --numproc 1 --numsteps 10 --width .1 --length 1"
"""
Example:
{0}

:author: Paul Sicsic
:date: 2011/01/02
""".format(USAGE)

import argparse
try :
    from launcher import main
    from parsing_options import mef90_parser, geom_bar_parser, base_pymef_parser
#    from pymef.input_output import save_dat
except ImportError as err:
    print "\n Unable to import pymef"
    print "check PYTHONPATH"
    print "ImportError: {0} \n".format(err)
    raise

if __name__ == '__main__':

    geom_parser = geom_bar_parser()
    py_parser = base_pymef_parser()
    pymef_parser = argparse.ArgumentParser(parents=[geom_parser, py_parser])
    pymef_opt, residual = pymef_parser.parse_known_args()
    pymef_opt.type = 'bar'
    print pymef_opt

    generic_parser = mef90_parser()
    mef90_parser = argparse.ArgumentParser(parents=[generic_parser])
    mef90_opt, residual2 = mef90_parser.parse_known_args(residual)
    print mef90_opt

    concrete = {'YOUNG': 1, 'POISSON':.2, 'TOUGHNESS':1 } 
    materials = [concrete]
    set_fix = {'Ux' : 0, 'Uy' : 0, 'Uz':0, 'V' :1}
    set_traction = {'Ux' : 1, 'Uy' : 0, 'Uz':0, 'V' :1}
    node_sets = [set_fix, set_traction]
#    save_dat(vars(pymef_opt), vars(mef90_opt))
    main(pymef_opt, mef90_opt, materials, node_sets)
