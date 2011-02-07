
#!/usr/bin/python2.5
# encoding: utf-8
import sys
import os
sys.path.append('%s'%(os.getcwd()+'/../lib/'))

from ComputationTests import *

class TractionOptions_Class(TestOptions_Class):
	def __init__(self):
		super(TractionOptions_Class, self).__init__()
		from optparse import OptionParser, OptionGroup

		group_Tractionopts=OptionGroup(self.parser, "Custom")
		group_Tractionopts.add_option("--height"	, dest="height"	, type="float"	, default=".1", 	help="Height of the bar [default: %default]")
		group_Tractionopts.add_option("--young"	, dest="young"	, type="float"	, default="1", 	help="Young modulus [default: %default]")
		group_Tractionopts.add_option("--poisson"	, dest="poisson"	, type="float"	, default="0", 	help="Poisson ratio [default: %default]")
		group_Tractionopts.add_option("--thexp"	, dest="thexp"	, type="float"	, default="1", 	help="Thermal expansion coeff [default: %default]")
		group_Tractionopts.add_option("--gc"	, dest="gc"	, type="float"	, default="1", 	help="Toughness [default: %default]")
		group_Tractionopts.add_option("--theta"	, dest="theta"	, type="float"	, default="1", 	help="Theta scale [default: %default]")

		self.parser.add_option_group(group_Tractionopts)
		self.parser.parse_args()
	
class TractionCtx_Class(TestCtx_Class):
	def __init__(self, options):
		self.meshpath='Cubit/'
		self.argspath='Cubit/'
		self.joupath=''

		self.meshfile="Traction-mesh-%.3e.gen"%(options.parser.values.meshsize)
		self.argsfile='Traction.args'
		self.joufile='Traction.jou'
		self.job='Traction-mesh%.3e' %(options.parser.values.meshsize)
		self.prep='PrepVarFrac'
		self.qs='VarFracQS2D'
		self.qsdir='VarFracQS'
		super(TractionCtx_Class, self).__init__(options)
	

class TractionComputation_Class(Computation_Class):
	
	def generatemesh(self, options, TestCtx):
		opts=options.parser.values

		import cubit
		cubit.init([""])
		cubit.cmd('reset')
		cubit.cmd('set journal off')
		cubit.cmd('set echo off')
		
		#Geometry
		cubit.cmd('create surface rectangle width 1 height %.3f zplane'%opts.height) 

		#Mesh
		cubit.cmd('surface 1 size %.3f scheme TriDelaunay'%opts.meshsize)
		cubit.cmd('mesh surface 1')

		# Blocks 
		cubit.cmd('block 1 surface 1')

		# Nodesets
		cubit.cmd('nodeset 1 curve  2')
		cubit.cmd('nodeset 2 curve  4')

		# mesh generation commands
		cubit.cmd('export mesh "%s" dimension 2 overwrite'%(TestCtx.meshpath+TestCtx.meshfile))
	
	def args(self, options, TestCtx):
		"""docstring for args"""
		opts=options.parser.values
		# word dictionary for custom replacement
		self.wd={'%YOUNG%': opts.young,
			'%THEXP%': opts.thexp,
			'%THETA%': opts.theta,
			'%POISSON%': opts.poisson,
			'%GC%': opts.gc}
		super(TractionComputation_Class, self).args(options, TestCtx)
	
	def pproc(self, options, TestCtx):
		"""Runs the post processing stage"""
		
		print "Start postprocessing of %s"%TestCtx.jid
		cmd="cd %s"%TestCtx.base_dir
		os.system(cmd)
		os.system('visit -nowin -cli -s '+visitscript)
		print "Done postprocessing of %s"%TestCtx.jid
		os.system('cd ..')
	
