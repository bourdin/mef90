#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
cli module 

Classes to add easily command lines

:date: 2011/02/01
:author: Paul Sicsic
"""

import gettext
_ = gettext.gettext

import argparse


def mef90_parser():
    """
    generic options for mef90 computation

    .. seealso:: m_VarFrac_Struc/m_VarFrac_Struc.F90

    :returns: argparse parser 
    """
    parser = argparse.ArgumentParser(add_help=False) 
    comp = parser.add_argument_group('comp', 'group1 description')
    comp.add_argument('-v', '--verbose', default=0, type=int,
                    help='increse verbosity')
    comp.add_argument('--epsilon', default=0.2, type=float, 
                    help='AT regularization parameter')
    comp.add_argument('--atnum', default=1, type=int, 
                    help='type of damage law')
    comp.add_argument('--altminmaxiter', default=1000, 
                help='maximum number of alternate minimization iterations')
#    comp.add_argument('--force', default=False, 
#                help='overwrite existing computation')
    comp.add_argument("--unilateral",  type=int, default=0,
                    help="Unilateral beahvior" 
                      "0: None"
                      "1: Full"
                      "2: Shear ")
    comp.add_argument("--irrev",  type=int, default=2,
                    help="Irreversibility type "
                  "0: No irreversibility"
                  "1: Equality constraint with threshold"
                  "2: Inequality constrains with threshold")
    comp.add_argument("--irrevtol",  type=float, default=10e-2,
                    help="Irreversibility tolerance")
    comp.add_argument("--kepsilon",  type=float, default=1e-6,
                    help="")
    comp.add_argument("--planestrain",  type=bool, default=False,
                    help="Plane strain computation ? ")

#    comp.add_argument("--initw" ,  type=int    , default=6,
#                    help="Iniatialization type for W [default: %default]")
    comp.add_argument("--initv" ,  type=int    , default=0,
                    help="Iniatialization type for V [default: %default]"
       "0: VarFrac_Init_V_PREV : use solution at previous time step (if available) or zero/"
       "1: VarFrac_Init_V_RND : NOT IMPLEMENTED"
       "2: VarFrac_Init_V_SPH : NOT IMPLEMENTED"
       "3: VarFrac_Init_V_CRACKS: NOT IMPLEMENTED"
       "4: VarFrac_Init_V_ONE : Inizialize with V=1"
       "5: VarFrac_Init_V_OSC : Initialize with a oscillatory field")
    comp.add_argument("--InitVLength"   ,  type=float  , default=0,
                    help="[default: %default]")
    comp.add_argument("--extra_cmd"   ,  type=str  , default='', 
                      help="extra cmd arguments")
    comp.add_argument("--saveblk",  type=bool, default=False, 
                help="Energy saved by blocks ? ")
    comp.add_argument("--savestress",  type=bool, default=False, 
                help="Save Stress ? ")
    comp.add_argument("--savestrain",  type=bool, default=False, 
                help="Save Strain ? ")
    comp.add_argument("--altminsaveint",  type=int, default=25, help=" ")

    return parser 

def base_pymef_parser(): 
    """pymef documentation  """
    parser = argparse.ArgumentParser(add_help=False)
    pymef = parser.add_argument_group('pymef', 'group1 description')  
    pymef.add_argument('--tmin', default=0, type=float, 
                    help='minimum value of the loading parameter')
    pymef.add_argument('--tmax', default=1, type=float, 
                    help='maximum value of the loading parameter')
    pymef.add_argument('--numproc', default=1, type=int, 
                    help='cpu count')
    pymef.add_argument('--meshsize', default= 0.05, type=float, 
                    help='mesh size')
    pymef.add_argument('-d', '--dim', default=2, type=int, 
                    help='dimension of the simulation')
    pymef.add_argument('--visit', default=1, 
                help='do postprocessing with visit')
    pymef.add_argument('--prep_mef90',
                    default='$MEF90_DIR/bin/$PETSC_ARCH/PrepVarFracNG', type=str,
                    help='mef90 command for the prep stage')
    pymef.add_argument('--comp_mef90',
                    default='$MEF90_DIR/bin/$PETSC_ARCH/VarFracQS{0}D', type=str,
                    help='mef90 command computation')
    pymef.add_argument('--numsteps', default=10, type=int, 
                       help='Number of time steps of the computation')
    pymef.add_argument('--batch', default=False, type=bool, 
                       help='batch or direct computation')
    pymef.add_argument('--pbsfile', default='sub.sh', type=str)

    return parser  

def geom_parser():
    """
    Basic parser for a geometry with a wifth, length and height 

    :returns: argparse parser 
    """
    parser = argparse.ArgumentParser(add_help=False)
    geom = parser.add_argument_group('geom', 'Geometry arguments') 

    geom.add_argument("--width",  type=float, default=.1, 
                    help="width of the bar")
    geom.add_argument("--length", type=float, default=1, 
                     help="length of the bar")
    geom.add_argument("--height",  type=float, default=.1,  
                    help="width of the bar")

    return parser 


def geom_bar_parser():
    """
    Parsing options fo the geometry of the 2d or 3d bar

    :returns: argparse parser 
    """
    pars_geom = geom_parser()
    parser = argparse.ArgumentParser(add_help=False, parents=[pars_geom])

    return parser

def geom_inclusions_parser():
    """
    Parsing options fo the geometry of the drying computations with inclusions

    :returns: argparse parser 
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--width",  type=float, default=30.,
                        help="width of the bar")
    parser.add_argument("--length", type=float, default=30.,
                        help="length of the bar")
    parser.add_argument("--radius", type=float, default=5.,
                        help="Radius size of the inclusion ")
    parser.add_argument("--sphere_ratio", type=float, default=0.,
                        help="Radius size of the inclusion ")
    parser.add_argument("--nb_inclusions", type=float, default=0.,
                        help="Radius size of the inclusion ")
    parser.add_argument("--inclusions_positions", type=list, default=[],
                        help="list of the positions of the inclusions ")

    return parser

def geom_ring_parser():
    """
    parsing options for the ring and peach problem

    :returns: argparse parser
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--in_radius", type=float, default=3.,
                        help="inner Radius of the ring ")
    parser.add_argument("--out_radius", type=float, default=5.,
                        help="outer radius  of the ring ")

    return parser 

def geom_fiber_parser():
    """
    parsing options for the fiber pull out problem 

    :returns: argparse parser
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--width",  type=float, default=30.,
                        help="width of the bar")
    parser.add_argument("--length", type=float, default=30.,
                        help="length of the bar")
    parser.add_argument("--radius", type=float, default=5.,
                        help="Radius size of the inclusion ")
    parser.add_argument("--height",  type=float, default=.1,
                        help="width of the bar")

    return parser 

def geom_pied_parser():
    """
    Parsing options for the geometry of the Pied Experiment

    Pour Identification de l'Endommagement Diffus

    :returns: argparse parser
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--width",  type=float, default=.1,
                    help="width of the bar")
    parser.add_argument("--length", type=float, default=1,
                    help="length of the bar")
    parser.add_argument("--height",  type=float, default=.1, 
                    help="width of the bar")
    parser.add_argument("--width_alu",  type=float, default=.1,
                    help="width of the bar")
    return parser 

def geom_3pbt_parser():
    """
    Parsing for the geometry of the Three point Bending Test

    :returns: argparse parser
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--width",  type=float, default=1,
                    help="width of the bar")
    parser.add_argument("--length", type=float, default=4,
                    help="length of the bar")
    parser.add_argument("--height",  type=float, default=.1,
                    help="")
    parser.add_argument("--interval",  type=float, default=.5,
                    help="")
    parser.add_argument("--e",  type=float, default=.25,
                    help="")
    parser.add_argument("--pbt_type",  type=str, default='simple',
                        help="type on fomputation : simple, crack, notched")
    parser.add_argument("--width_notch",  type=float, default=.1,
                    help="width of the bar")
    parser.add_argument("--length_notch", type=float, default=.2,
                    help="length of the bar")

    return parser 


def geom_4pbt_parser():
    """
    Parsing for the geometry of the Three point Bending Test

    :returns: argparse parser
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--width",  type=float, default=1,
                    help="width of the bar")
    parser.add_argument("--length", type=float, default=4,
                    help="length of the bar")
    parser.add_argument("--height",  type=float, default=.1,
                    help="width of the bar")
    parser.add_argument("--interval",  type=float, default=.5,
                    help="width of the bar")
    parser.add_argument("--interval_up",  type=float, default=.5,
                    help="width of the bar")
    parser.add_argument("--e",  type=float, default=.25,
                    help="width of the bar")


    return parser 



def geom_sen_parser():
    """
    Parsing the geometry of the Sent test
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--width",  type=float, default=.1,
                    help="width of the bar")
    parser.add_argument("--x1",  type=float, default=.02,
                    help="x1")
    parser.add_argument("--x2",  type=float, default=.18,
                    help="x2")
    parser.add_argument("--x3",  type=float, default=.04,
                    help="x3")
    parser.add_argument("--x4",  type=float, default=.18,
                    help="x4")
    parser.add_argument("--x5",  type=float, default=.02,
                    help="x5")
    parser.add_argument("--width_notch",  type=float, default=.005,
                    help="width of the notch")
    parser.add_argument("--depth_notch",  type=float, default=.02,
                    help="depth of the notch")
    parser.add_argument("--e",  type=float, default=.02,
                    help="width of the NS")
    return parser 





def backtracking_parser():
    """
    Backtracking 

    :returns: argparse parser 
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--dobt"    , action="store_true"   , default=False,  
                                     help="do BT or not? [default: %default]")
    parser.add_argument("--bttol"   , type=float, default="10e-2",  help="")
    parser.add_argument("--btint"   , type=int  , default="10",     help="")
    parser.add_argument("--btscope" , type=int  , default="10000",  help="")

    return parser 

def tao_parser():
    """
    Tao Parser

    :returns: argparse parser
    """
    parser = argparse.ArgumentParser(add_help=False) 
    parser.add_argument("--U_tao"   , type=float, default="10e-2",  help="") 

    return parser  


