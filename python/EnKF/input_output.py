#!/usr/bin/env python
# encoding: utf-8
 
 
"""
Function to save all the parameters to the .dat file or read from in 
:author: Paul Sicsic
:date: 2011/02/09
"""

import ConfigParser
#The ConfigParser module has been renamed to configparser in Python 3.0 
from parsing_options import mef90_parser, base_pymef_parser  



def save_dat(py_opt, mef90_opt, dat_file="options.dat"):
    """ saved the two dictionaries to the file""" 
    config = ConfigParser.RawConfigParser()
    config.add_section("PYTHON_OPTIONS")
    for key, value in py_opt.iteritems() : 
        config.set("PYTHON_OPTIONS", str(key), value)

    config.add_section("MEF90_OPTIONS")
    for key, value in mef90_opt.iteritems() : 
        config.set("MEF90_OPTIONS", str(key), value)


    with open(dat_file, 'wb') as configfile:
             config.write(configfile)


class SetClass( object ):
    pass

def read_dat(input_file):
    """ 
    read from the input file and retruns two dictionaries with the python and mef90 parameters

    :param input_file: Name of the saved file with all the parameters to be loaded  
    :type input_file: string
    """ 
    config = ConfigParser.SafeConfigParser()
    config.read(input_file)

    mef_parser = mef90_parser() 
    mef90_opt = mef_parser.parse_args("")
    py_parser = base_pymef_parser()
    py_opt = py_parser.parse_args("")

    for name_section in config.sections():
        if name_section == "PYTHON_OPTIONS":
            for  n_param in config.options(name_section):
                #print n_param
                if n_param in ('comp_mef90', 'prep_mef90', 'type', 'root_dir',
                              'pbsfile', 'job_id', 'argsfile', 'meshfile',
                               'base_dir', 'pbt_type'):
                    setattr(py_opt, n_param, config.get(name_section, n_param).split(',')[0])  
                elif n_param in ('dim', 'numproc', 'numsteps'): 
                    setattr(py_opt, n_param, config.getint(name_section, n_param))
                elif n_param in ('batch') :
                    setattr(mef90_opt, n_param, config.getboolean(name_section, n_param))
                else : 
                    setattr(py_opt, n_param, config.getfloat(name_section, n_param)) 
        elif name_section == "MEF90_OPTIONS":
            for  n_param in config.options(name_section):
                #print n_param
                if n_param in ('force', 'extra_cmd'):
                    setattr(mef90_opt, n_param, config.get(name_section, n_param).split(',')[0])
                elif n_param in ('atnum', 'verbose', 'altminmaxiter') :
                   setattr(mef90_opt, n_param, config.getint(name_section, n_param))  
                elif n_param in ('savestrain', 'savestress', 'saveblk') :
                    setattr(mef90_opt, n_param, config.getboolean(name_section, n_param))
                else:
                    setattr(mef90_opt, n_param, config.getfloat(name_section, n_param))  

    return py_opt, mef90_opt 

