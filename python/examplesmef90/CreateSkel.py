#!/usr/bin/env python
# encoding: utf-8
import os
import sys
import re
from optparse import OptionParser


Readme="""
%TESTNAME% Test 
with the variational approach to fracture mechanics: ultifissuration and delamination of thin films

Refs:
Francfort, G. A., & Marigo, J.-J. (1998). Revisiting Brittle Fracture as an Energy Minimization Problem. Journal of the Mechanics and Physics of Solids, 46, 1319-1342.
Bourdin, B., Francfort, G. A., & Marigo, J.-J. (2007). The Variational Approach to Fracture. Journal of Elasticity.
"""

Args="""
"""

Testpy="""
#!/usr/bin/python2.5
# encoding: utf-8
import sys
import os
sys.path.append('%s'%(os.getcwd()+'/../lib/'))

from ComputationTests import *

class %TESTNAME%Options_Class(TestOptions_Class):
	def __init__(self):
		super(FourOptions_Class, self).__init__()
		from optparse import OptionParser, OptionGroup

		group_%TESTNAME%opts=OptionGroup(self.parser, "Custom")

		self.parser.add_option_group(group_%TESTNAME%opts)
		self.parser.parse_args()
	
class %TESTNAME%Ctx_Class(TestCtx_Class):
	def __init__(self, options):
		self.meshpath='Cubit/'
		self.argspath='Cubit/'
		self.joupath=''

		self.meshfile="%TESTNAME%-mesh-%.3e.gen"%(options.parser.values.meshsize)
		self.argsfile='%TESTNAME%.args.txt'
		self.joufile='%TESTNAME%.jou.txt'
		self.job='%TESTNAME%-mesh%.3e' %(options.parser.values.meshsize)
		self.prep='PrepFilm'
		self.qs='VarFilmQS'
		super(FourCtx_Class, self).__init__(self)
	

class %TESTNAME%Computation_Class(Computation_Class):
	
	def generatemesh(self, options, TestCtx):
		opts=options.parser.values

		import cubit
		cubit.init([""])
		cubit.cmd('reset')
		cubit.cmd('set journal off')
		cubit.cmd('set echo off')

		# mesh generation commands
		cubit.cmd('export mesh "%s" dimension 2 overwrite'%(TestCtx.meshpath+TestCtx.meshfile))
	
	def args(self, options, TestCtx):
		\"""docstring for args\"""
		opts=options.parser.values
		# word dictionary for custom replacement
		self.wd={}
		super(%TESTNAME%Computation_Class, self).args(options, TestCtx)
	
	def pproc(self, options, TestCtx):
		\"""Runs the post processing stage\"""
		
		print "Start postprocessing of %s"%TestCtx.jid
		cmd="cd %s"%TestCtx.base_dir
		os.system(cmd)
		os.system('visit -nowin -cli -s '+visitscript)
		print "Done postprocessing of %s"%TestCtx.jid
		os.system('cd ..')
	
"""

Run="""
#!/usr/bin/env python
# encoding: utf-8

import sys
import os
sys.path.append('%s'%(os.getcwd()+'/../lib/'))

from %TESTNAME%Test         import *

def runcircle(options, TestCtx):
	# Init
	%TESTNAME%=%TESTNAME%Computation_Class(options, TestCtx)
	
	# Mesh generation
	
	Four.mesh(options, TestCtx)
	Four.args(options, TestCtx)
	Four.prep(options, TestCtx)
	#Four.qsevo(options, TestCtx)
	
def main():
	# Parse cmdline opts
	cmd="MEF90_DIR=\"/Users/kumiori/Documents/Universita/Dottorato/Code/mef90-sieve\""
	os.system(cmd)
	os.environ.get("CUBIT_HOME")
	options=FourOptions_Class()
	(opts, args) = options.parser.parse_args()
	# print opts
	print options.parser.values
	
	if options.parser.values.comptype!='single':
		nsteps=int(options.parser.values.pattern.split(':')[1])
		
		for step in range(nsteps):
			# (opts, args) = options.parser.parse_args()
			options.update(step)
			TestCtx=%TESTNAME%Ctx_Class(options)
			print "Computation id: %s"%TestCtx.jid
			
			runcircle(options, TestCtx)
			%TESTNAME%.pproc(options, TestCtx)
			
	elif opts.comptype=='single':
		TestCtx=%TESTNAME%Ctx_Class(options)
		runcircle(options, TestCtx)
		print "Computation id: %s"%TestCtx.jid
		
		
	print "Done with computation with id: %s"%TestCtx.jid
	%TESTNAME%.pproc(options, TestCtx)
	
if __name__ == '__main__':
	main()

"""

def main():
	parser=OptionParser()
	parser.add_option("-n", dest="n", type="string"	, help="")
	(opts, args)=parser.parse_args()
	testname=opts.n
	
	os.mkdir(testname)
	os.mkdir(testname+'/Cubit')
	argsf=open("%s/Cubit/%s.args"%(testname, testname),"w")
	testpyf=open("%s/%s.py"%(testname, testname),"w")
	readmef=open("%s/Readme.txt"%(testname),"w")
	runf=open("%s/Run%s.py"%(testname, testname),"w")
	
	argsf.write(re.sub("%TESTNAME%",testname, Args))
	testpyf.write(re.sub("%TESTNAME%",testname, Testpy))
	readmef.write(re.sub("%TESTNAME%",testname, Readme))
	runf.write(re.sub("%TESTNAME%",testname, Run))
	
if __name__ == '__main__':
	main()

