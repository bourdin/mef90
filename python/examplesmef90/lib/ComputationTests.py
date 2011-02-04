#!/usr/bin/python2.5
# waiting for 2.6+ support for cubit

# set default options
# set default paths
# create dirtree for each test
import os

class Computation_Class(object):
	def __init__(self, options, TestCtx):
		opts=options.parser.values
		
		# create directory tree
		if os.path.exists(TestCtx.base_dir):
			if not opts.force:
				print "Destination folder %s exists. Erase or add --force flag" %TestCtx.base_dir
				return -1
		else:
			# os.system("mkdir %s"%base_dir)
			cmd='mkdir -p %s'%TestCtx.base_dir
			os.system(cmd)
	
	
	def mesh(self, options, TestCtx):
		opts=options.parser.values
		
		print '**** looking for %s'%(TestCtx.root_dir+TestCtx.meshpath+TestCtx.meshfile)
		if not os.path.exists('%s'%(TestCtx.root_dir+TestCtx.meshpath+TestCtx.meshfile)):
			if opts.debug:
				print "pre-generated mesh %s not found. Generating new mesh "%(TestCtx.root_dir+TestCtx.meshpath+TestCtx.meshfile)
			if not os.path.exists('%s'%(TestCtx.root_dir+TestCtx.meshpath)):
				cmd='mkdir -p %s'%(TestCtx.root_dir+TestCtx.meshpath)
				os.system(cmd)
			self.generatemesh(options, TestCtx)
		else: print '**** Mesh found!'
	
	
	def generatemesh(self, options, TestCtx):
		cubit_cmd="cubit.command"
		
		injoufile_name="%s/Cubit/%s.jou"%(TestCtx.root_dir, TestCtx.joufile)
		outjoufile_name="%s/Cubit/%s_mesh%.3f.jou"%(TestCtx.root_dir, TestCtx.joufile, TestCtx.meshsize)
		if options.debug:
			print "jououtfile is %s"%outjoufile_name
			print "jouinfile is %s"%injoufile_name
		outjoufile=open(outjoufile_name,'w')
		outjoufile.write(injoufile.read())
		outjoufile.close()
		cmd='%s -nojournal  -batch -nographics -input %s'%(cubit_cmd, outjoufile_name)
		os.system(cmd)
	
	def args(self, options, TestCtx):
		"""Populates .args file relating to current computation"""
		opts=options.parser.values
		self.wd.update({'%TMAX%': opts.tmax,
				'%NUMSTEPS%': opts.numsteps})
		# open file
		inargsfile_name="%s"%(TestCtx.root_dir+TestCtx.argspath+TestCtx.argsfile)
		outargsfile_name="%s"%(TestCtx.base_dir+TestCtx.argsfile)
		inargsfile=open(inargsfile_name,'r')
		outargsfile=open(outargsfile_name,'w')
		outargs=inargsfile.read()
		# replace
		for key in self.wd:
		        outargs = outargs.replace(key, str(self.wd[key]))
		# write
		outargsfile.write(outargs)
	
	def prep(self, options, TestCtx):
		"""Runs the prep stage"""
		import sys
		import os
		opts=options.parser.values
		print 'base:' +TestCtx.base_dir
		print 'root:' +TestCtx.root_dir
		print 'cubit:' +TestCtx.root_dir+TestCtx.meshpath
		cmd="cp %s %s"%(TestCtx.root_dir+TestCtx.meshpath+TestCtx.meshfile, TestCtx.base_dir)
		if options.parser.values.debug==True:
			print cmd
		os.system(cmd)
		
		cmd="cd %s && mpiexec -n %i %s -p %s -i %s && cd %s"%(
			TestCtx.base_dir,
			opts.n,
			os.environ.get('MEF90_DIR')+'/'+TestCtx.prep+'/'+os.environ.get('PETSC_ARCH')+'/'+TestCtx.prep,
			TestCtx.meshfile[:-4], # weak! FIX
			TestCtx.argsfile,
			TestCtx.root_dir)
		if options.parser.values.debug==True:
			print cmd
		os.system(cmd)
				
	
	def qsevo(self, options, TestCtx):
		"""Runs the QS code"""
		opts=options.parser.values
		args=''
		for optgroup in options.parser.option_groups: 
			if optgroup.title != 'Custom' and optgroup.title != 'Misc params':
				# print optgroup.title 
				# print 'title prted'
				for optitem in optgroup.option_list:
					optname=optitem.get_opt_string()[2:]
					if optitem.action == 'store_true':
						if getattr(options.parser.values, optname)!='False':
							args += ' -%s'%(optname)
					elif optitem.action == 'store':
						args += ' -%s %s'%(optname, getattr(options.parser.values, optname))
		
		cmd="cd %s && mpiexec -n %i %s -p %s %s && cd %s"%(
			TestCtx.base_dir,
			opts.n,
			os.environ.get('MEF90_DIR')+'/'+TestCtx.qs+'/'+os.environ.get('PETSC_ARCH')+'/'+TestCtx.qs,
			TestCtx.meshfile[:-4], # weak! FIX
			args,
			TestCtx.root_dir)
		if options.parser.values.debug==True:
			print args
			print cmd
		os.system(cmd)
		
	def pproc(self, options, TestCtx):
		"""Runs the post processing stage"""
		
		if opts.debug:
			print "Start postprocessing of %s"%TestCtx.jid
		cmd="cd %s"%TestCtx.base_dir
		os.system(cmd)
		os.system('visit -nowin -cli -s '+visitscript)
		if opts.debug:
			print "Done postprocessing of %s"%TestCtx.jid
		os.system('cd ..')
			
		
class TestOptions_Class(object):
	def __init__(self): 
		from optparse import OptionParser, OptionGroup
			
		self.parser=OptionParser()
		
		group_misc=OptionGroup(self.parser, "Misc params")
		group_misc.add_option("--tmax", dest="tmax", type="float", default=5, help="maximum value of the loading parameter")
		group_misc.add_option("--numsteps", dest="numsteps", type="int", default=11, help="number of time steps")
		group_misc.add_option("--visit", dest="visit", type="int",default=1,help="do postprocessing with visit")
		group_misc.add_option("--force", dest="force", action="store_true", help="overwrite existing computations")
		group_misc.add_option("--meshsize", dest="meshsize", type="float", help="mesh size scale")
		group_misc.add_option("--n", dest="n", type="int", default=4, help="cpu count")
		group_misc.add_option("--debug", dest="debug", action="store_true", default="False",help="Show debug info [default: %default]")
		self.parser.add_option_group(group_misc)
		
		# Computation type
		self.parser.add_option("--comptype"         , dest="comptype"		, type="string"	, default='single', help=" Computation type [default: %default]")
		self.parser.add_option("--pattern"	, dest="pattern", type="string", help="Pattern for computation, i.e. start:step:end")
		
		# Parallel 
		group_parallel=OptionGroup(self.parser, "Parallel and MPI environnement params")
		self.parser.add_option_group(group_parallel)
		
		# Numerics and Solver parameters
		group_numsolv=OptionGroup(self.parser, "Numerics and Solver parameters")
		group_numsolv.add_option("--keps"         , dest="keps"		, type="float"	, default="10e-6", help="[default: %default]")
		group_numsolv.add_option("--epsilon", dest="epsilon",type="float", help="AT regularization parameter")
		group_numsolv.add_option("--integorder"   , dest="integorder"	, type="int"  	, default="2", help="[default: %default]")
		group_numsolv.add_option("--altminmaxiter", dest="altminmaxiter"	, type="int"  	, default="1000", help="maximum number of alternate minimization iterations [default: %default]")
		group_numsolv.add_option("--altmintol"    , dest="altmintol"	, type="float"	, default="10e-4", help="[default: %default]")
		self.parser.add_option_group(group_numsolv)

		
		# Output options
		group_output=OptionGroup(self.parser, "Output options")
		group_output.add_option("--saveblk"	, dest="saveblk"	, action="store_true"	, default="False", 	help="[default: %default]")		
		group_output.add_option("--altminsaveint"	, dest="altminsaveint"	, type="int"		, default="25", 	help="")     
		group_output.add_option("--savestress"	, dest="savestress"	, action="store_true"	, default="False", 	help="[default: %default]")		
		group_output.add_option("--savestrain"	, dest="savestrain"	, action="store_true"	, default="False", 	help="[default: %default]")		
		self.parser.add_option_group(group_output)
		
		
		# Unilateral
		group_unilat=OptionGroup(self.parser, "Unilateral condition")
		group_unilat.add_option("--unilateral"         , dest="unilateral"		, type="int"	, default="0", help="")		
		self.parser.add_option_group(group_unilat)

		# Model
		self.parser.add_option("--atnum"         , dest="atnum"		, type="int"	, default="1", help="AT approximation model [default: %default]"
			"1: (3/8)*(V/epsilon+epilon*grad(V)^2) as approximation of the crack length"
			"2: (1/2)*(V^2/epsilon+epilon*grad(V)^2) as approximation of the crack length")		

		# Backtracking
		group_bt=OptionGroup(self.parser, "BT options")
		group_bt.add_option("--dobt"	, dest="dobt"	, action="store_true"	, default="False", 	help="do BT or not? [default: %default]")		
		group_bt.add_option("--bttol"	, dest="bttol"	, type="float"		, default="10e-2", 	help="")     
		group_bt.add_option("--btint"	, dest="btint"	, type="int"		, default="10", 	help="")     
		group_bt.add_option("--btscope"	, dest="btscope"	, type="int"		, default="10000", 	help="")     
		self.parser.add_option_group(group_bt)
		
		# Irreversibility
		group_irrev=OptionGroup(self.parser, "Irreversibility options")
		group_irrev.add_option("--irrevtype"	, dest="irrevtype"	, type="int"		, default="2"		, 	help="Irreversibility type"
			"0: No irreversibility"                
			"1: Equality constraint with threshold"
			"2: Inequality constrains with threshold")     
		group_irrev.add_option("--irrevtol"	, dest="irrevtol"	, type="float"		, default="10e-2"		, 	help="Irreversibility tolerance")     
		
		self.parser.add_option_group(group_irrev)

		# Initialization

		group_init=OptionGroup(self.parser, "Initialization")
		group_init.add_option("--initw"	, dest="initw"	, type="int"	, default="6", 	help="Iniatialization type for W [default: %default]")
		group_init.add_option("--initv"	, dest="initv"	, type="int"	, default="0", 	help="Iniatialization type for V [default: %default]"
			"0: VarFrac_Init_V_PREV : use solution at previous time step (if available) or zero/"
			"1: VarFrac_Init_V_RND : NOT IMPLEMENTED"                                            
			"2: VarFrac_Init_V_SPH : NOT IMPLEMENTED"                                            
			"3: VarFrac_Init_V_CRACKS: NOT IMPLEMENTED"                                          
			"4: VarFrac_Init_V_ONE : Inizialize with V=1"                                        
			"5: VarFrac_Init_V_OSC : Initialize with a oscillatory field")		
		group_init.add_option("--initvlenght"	, dest="initvlenght"	, type="float"	, default="0", 	help="[default: %default]")		
		self.parser.add_option_group(group_init)

	def checkopts(self):
		"""Checks whether options are consistent or not"""
		# optvals=options.parser.values
		# if optvals.dobt and not opts.a:
		# 	print "Option B requires option A\n"
		# 	parser.print_help()
		# 	exit(-1)

	def update(self, i):
		"""Updates options"""
		
		start=float(self.parser.values.pattern.split(':')[0])
		stop=float(self.parser.values.pattern.split(':')[2])
		nsteps=float(self.parser.values.pattern.split(':')[1])
		inc=(stop-start)/nsteps
		optname=self.parser.values.comptype
		
		if i==0:
			optname=self.parser.values.comptype

			if self.parser.values.comptype!='single':
				if optname=='numproc': 
					self.parser.values.numproc=int(self.parser.values.pattern.split(':')[0])
				else: 
					setattr(self.parser.values, optname, float(self.parser.values.pattern.split(':')[0]))	
		setattr(self.parser.values, optname, getattr(self.parser.values, optname)+i*inc)
		
	
class TestCtx_Class(object):
	def __init__(self, options):
		import uuid
		import os
		import socket
		
		self.jid=str(uuid.uuid4())[1:8]
		
		self.debug=True
		
		self.root_dir=os.getcwd()+'/'
		if options.parser.values.comptype!='single':
			patlist=options.parser.values.pattern.split(':')
			compstr=options.parser.values.comptype+'-%.2f-%.2f-%.2f/'%(float(patlist[0]), float(patlist[1]), float(patlist[2]))
		else:
			compstr=''
		self.base_dir="%s-%s"%(self.root_dir+compstr+self.job, self.jid)+'/'
	
	def __str__(self):
		print "root_dir is %s"%self.root_dir
		print "base_dir is %s"%self.base_dir
		print "mesh file is %s"%self.meshfile
		print "argsfile is %s"%self.argsfile
		print "joufile is %s.jou"%self.joufile
	
