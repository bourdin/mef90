#!/usr/bin/env python
# encoding: utf-8

""" 
A set of functions to launch a mef90-sieve computation from python cli

:date: 2010/12
"""

import os, sys, uuid 
import subprocess
import fileinput, shutil 
try :
    from input_output import save_dat
#    from pymef.plot_energy import plot_from_path
    from meshes import mesh_selection
#    from pymef.visit.wrapper import wrapper_visit
    from dictionaries import JOB_ID, MESH_FILE, ARGS_FILE
except ImportError as err:
    print "\n Unable to import pymef"
    print "check PYTHONPATH"
    print "ImportError: {0} \n".format(err)
    raise

def infofile(base_dir, options, job_id, mef90_opt):
    """
    print in a file the main parameters of the computation 

    :param base_dir: The directory xhere the computaion is taking place
    :type base_dir: string
    :param options: The dictionnary/namespace used during the current computation
    :param job_id: The job which in going to be executed (directory name usually)
    :type job_id: string 

    :date: 2010
    :author: Corrado Maurini
    """
    info_file = open("%s/00_INFO.txt"% base_dir, 'w')
    info_file.write("JOB_ID            %s\n"%job_id)
    info_file.write("EPSILON           %f\n"%mef90_opt.epsilon)
    info_file.write("MESHSIZE          %f\n"%options.meshsize)
    info_file.write("NUMSTEPS          %s\n"%options.numsteps)
    info_file.write("TMAX              %f\n"%options.tmax)
    info_file.write("ATNUM              %f\n"%mef90_opt.atnum)
    info_file.close()


def debug(options):
    """ 
    Printing information for debugging

    :date: 2010
    :author: Corrado Maurini 

    :param options: Set of options for a mef90 computation to pe printed  
    :type options: namespace
    """
    print("root_dir is %s"% options.options.root_dir)
    print("base_dir is %s"% options.base_dir)
    print("mesh file is %s"% options.meshfile)
    print("argsfile is %s"% options.argsfile)
    print("joufile is %s.jou"% options.joufile)


def replace_expression(file_in, keyword, new_expression):
    """
    :author: Paul Sicsic
    :date: 2010/08/05

    Replace in file_in all the occurence of keyword by new_expression

    :param file_in: a file name
    :param file_in: string
    :param keyword: keyword to be replaced
    :type keyword: string
    :param new_expression: By what the key_word should be replaced 
    :type new_expression: string
    """
#TODO : select case if in is a list of keys are just a key
    for line in fileinput.input(file_in, inplace=1):
        if keyword in line:
            line = line.replace(keyword, new_expression)
        sys.stdout.write(line)

def edit_args(args_file, list_dict): 
    """ 
    Takes a list of dictionaries and for each item replaces the key by it's value

    :author: paul
    :date: 2011/02.09
    """
    for ind, material in enumerate(list_dict) :
        for key, value in material.iteritems() :
            string = "%"+key+"{0:04d}%".format(ind+1) 
            replace_expression(args_file, string, str(value))


def set_args_file(options, materials, node_sets):
    """
    Set specific for the .args file 

    :date: 2011/02/03
    :author: Paul Sicsic 

    :param options: Options for the mef90 computation
    :param material_prop: Material properties for the different blocks to edit the .args file 
    :type material_prop: list, tuple
    :param ns_prop: Node sets properties to est the .args file
    :type ns_prop: list, tuple
    """
    ### generating the .args file
    input_args = os.sep.join((options.root_dir, options.argsfile))
    if not(os.path.isfile(input_args)):
        args_main_dir = os.sep.join((os.path.split( __file__ )[0], 'args'))
        input_args = os.sep.join((args_main_dir, options.argsfile)) 

    try : 
        shutil.copy(input_args, options.base_dir)
    except IOError as err: 
        print("The file {0} is neither in".format(options.argsfile))
        print("The current directory")
        print("Pymef.args directory")
        print("ImportError: {0} \n".format(err))
        raise

    replace_expression(options.argsfile, "%TMIN%", "%f"%options.tmin)
    replace_expression(options.argsfile, "%TMAX%", "%f"%options.tmax)
    replace_expression(options.argsfile, "%NUMSTEPS%",
                       "%i"%options.numsteps)

    edit_args(options.argsfile, materials)
    edit_args(options.argsfile, node_sets)

    my_file = open(options.argsfile, "r") #Opens the file in read-mode
    text = my_file.read()
    if text.find("%") >0 : 
        raise ValueError("The character # is still in the .args file")

def set_pbs_file(options, mef90_opt ):
    """
    edit pbs file to launch a batch comutation

    :date: 2011/03/02
    :author: Paul Sicsic
    """
    ### generating the .args file
    input_pbs = os.sep.join((options.root_dir, options.pbsfile))

    if not(os.path.isfile(input_pbs)):
        pbs_main_dir = os.sep.join((os.path.split( __file__ )[0], 'pbs'))
        input_pbs = os.sep.join((pbs_main_dir, options.pbsfile)) 
    try :
        print input_pbs
        shutil.copy(input_pbs, options.base_dir)
    except IOError as err: 
        print("The file {0} is neither in".format(options.pbsfile))
        print("The current directory")
        print("Pymef.pbs directory")
        print("ImportError: {0} \n".format(err))
        raise

    cmd_prep, cmd_varfrac = commands(options, mef90_opt)

    replace_expression(options.pbsfile, "%JOB_ID%", options.job_id)
    replace_expression(options.pbsfile, "%PREP_CMD%", cmd_prep)
    replace_expression(options.pbsfile, "%VARFRAC_CMD%", cmd_varfrac)
#    replace_expression(options.pbsfile, "%NB_PROC%", "%i"%options.numproc)
    replace_expression(options.pbsfile, "%JOB_PATH%", os.sep.join((options.root_dir, options.job_id)))

    # Ask for the correct ressource
    if options.numproc < 13 :
        replace_expression(options.pbsfile, "%NODES%", '1')
        replace_expression(options.pbsfile, "%NODE_CPU%", "%i"%options.numproc)
    elif options.numproc < 25:
        replace_expression(options.pbsfile, "%NODES%", '2')
        replace_expression(options.pbsfile, "%NODE_CPU%", '12')
    elif options.numproc < 37:
        replace_expression(options.pbsfile, "%NODES%", '3')
        replace_expression(options.pbsfile, "%NODE_CPU%", '12')
    elif options.numproc < 49:
        replace_expression(options.pbsfile, "%NODES%",'4')
        replace_expression(options.pbsfile, "%NODE_CPU%", '12' )
    else:
        raise ValueError("too many proc")


    my_file = open(options.pbsfile, "r") #Opens the file in read-mode
    text = my_file.read()
    if text.find("%") >0 : 
        raise ValueError("The character % is still in the .pbs file")



def runcircle(options, mef90_opt, materials_prop, ns_prop, extra_cmd="", debug=False):
    """ 
    main computation cycle

    :param options: All the options for the mef90 computation 
    :param materials_prop: material properties 
    :type materials_prop: list, tuple 
    :param ns_prop: node stes properties values
    :type ns_prop: list, tuple  
    """

     ### populate a 00_INFO.txt file
    infofile(options.base_dir, options, options.job_id, mef90_opt)

    set_args_file(options, materials_prop, ns_prop) 

    ### run the Prep stage
    os.chdir(options.base_dir)
    save_dat(vars(options), vars(mef90_opt))
    if not(options.batch):
        rundirect(options, mef90_opt)
    elif options.batch : 
        runbatch(options, mef90_opt)

def runbatch(options, mef90_opt):  
    """
    Run directly the computation on the current computer

    :author: Paul Sicsic
    :date: 2011/03/02
    """
    #copy pbs file
    #edit pbs file
    set_pbs_file(options, mef90_opt) 
    # BE CAREFUL WITH THE EN?VIRONNEMENT 
    #launch script
    os.system("qsub {0}".format(options.pbsfile)) 


def commands(options, mef90_opt):
    """
    :date 2011/03/02
    :author: Paul Sicsic
    """

    prep = "mpiexec -n {0} {1} -p {2} -i {3}".format(
                            options.numproc, 
                            options.prep_mef90,
                            options.meshfile.split('.gen')[0], 
                            options.argsfile )

    cmd_args = ''
    dict_opt =  vars(mef90_opt)
    for opt in dict_opt:
        if opt is not 'extra_cmd':
            cmd_args = cmd_args + ' -' + opt +' ' +  str(dict_opt[opt])
        else :
            cmd_args = cmd_args + ' ' + str(dict_opt['extra_cmd'])
#TODO not robust enough

    varfrac = "mpiexec -n {0} {1} -p {2}  {3}".format(
                            options.numproc, 
                            options.comp_mef90.format(options.dim),
                            options.meshfile.split('.gen')[0], 
                            cmd_args )

    varfrac += mef90_opt.extra_cmd

    return prep, varfrac

def rundirect(options, mef90_opt):  
    """
    Run directly the computation on the current computer

    :author: Paul Sicsic
    :date: 2011/03/02
    """
    cmd_prep, cmd_varfrac = commands(options, mef90_opt) 
    print cmd_prep
    print cmd_varfrac

    try:
        retcode = subprocess.call(cmd_prep, shell=True)
        if retcode < 0:
            print >>sys.stderr, "Child was terminated by signal", -retcode
        else:
            print >>sys.stderr, "Child returned", retcode
            raise OSError
    except OSError, e:
            print >>sys.stderr, "Execution failed:", e

    ### Run VarFraq
    try:
        retcode = subprocess.call(cmd_varfrac, shell=True)
        if retcode < 0:
            print >>sys.stderr, "Child was terminated by signal", -retcode
        else:
            print >>sys.stderr, "Child returned", retcode
            raise OSError
    except OSError, e:
            print >>sys.stderr, "Execution failed:", e

def conf_name_files(options, mef90_opt, call, debug=False):
    """ 
    Configure Files and Directories names

    :type options: namespace
    :param call: type of computation is it a test ? 
    :type call: string 
    :type debug: boolean  
    :param options: list of all the options to set a Variational Computation with mef90
    """

    job_uuid4 = '%s'% uuid.uuid4()
    job_uuid4_r = '-JOB_'+job_uuid4[1:8]
    job_par = JOB_ID[options.type].format(options, mef90_opt)
    if call == 'test' : 
        options.job_id = job_par 
    else : 
        options.job_id = job_par + job_uuid4_r 
    options.root_dir = os.getcwd()
    options.base_dir = "%s/%s"% (options.root_dir, options.job_id)
    options.meshfile = MESH_FILE[options.type].format(options)
    options.argsfile = ARGS_FILE[options.type].format(options)
  #  options.joufile = "Bar2D"

    if debug:
        debug(options)

    if os.path.exists(options.base_dir):
        if not options.force:
            print "Destination folder %s exists. Erase or add --force flag"
            return -1
    else:
        os.mkdir(options.base_dir)

def main(options, mef90_opt, materials_prop, ns_prop, call='cline', dict_values='', mesh_function=''):
    """
    - run cycle
    - postprocessing with visit

    :param options: Options for the mef90 computation
    :param material_prop: Material properties for the different blocks to edit the .args file 
    :type material_prop: list, tuple
    :param ns_prop: Node sets properties to est the .args file
    :type ns_prop: list, tuple
    """

    print options
###### If apped the dictionary with the correct names
    if not(dict_values==''):#not robust enough
        try :
            JOB_ID[options.type] =   dict_values['job_id']
            MESH_FILE[options.type] = dict_values['mesh']
            ARGS_FILE[options.type] = dict_values['args']
        except:
            raise IOError("Have you correctly configerd the dict_values dictionary ? ")

## Configure files and directories name ###
    conf_name_files(options, mef90_opt, call)

    if os.path.isfile(MESH_FILE[options.type]):
#the mesh file exists in the current directory - then just copy it
        shutil.copy(MESH_FILE[options.type], options.base_dir)
        os.chdir(options.base_dir)
    else : 
#### Mesh generation through cubit command line####     
        os.chdir(options.base_dir)
        if options.type in ('bar', 'therm', 'PM'): 
            mesh_selection(options)
        else :
            raise ValueError("Type of computation is not supported, modify pymef before launching computation")
            #mesh_function(options)

    runcircle(options, mef90_opt, materials_prop, ns_prop, extra_cmd='')
    print "Done with computation in %s"% options.job_id

### Plot energies ### 
#    plot_from_path(options.base_dir, options) 
    print os.path.dirname(options.base_dir)
    os.chdir(os.path.dirname(options.base_dir))
