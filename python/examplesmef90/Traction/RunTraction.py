
#!/usr/bin/env python
# encoding: utf-8

import sys
import os
sys.path.append('%s'%(os.getcwd()+'/../lib/'))

from Traction         import *

def runcircle(options, TestCtx):
	# Init
	Traction=TractionComputation_Class(options, TestCtx)
	
	# Mesh generation
	
	Traction.mesh(options, TestCtx)
	Traction.args(options, TestCtx)
	Traction.prep(options, TestCtx)
	Traction.qsevo(options, TestCtx)
	
def main():
	# Parse cmdline opts
	# cmd="MEF90_DIR="/Users/kumiori/Documents/Universita/Dottorato/Code/mef90-sieve""
	# os.system(cmd)
	os.environ.get("CUBIT_HOME")
	options=TractionOptions_Class()
	(opts, args) = options.parser.parse_args()
	# print opts
	print options.parser.values
	
	if options.parser.values.comptype!='single':
		nsteps=int(options.parser.values.pattern.split(':')[1])
		
		for step in range(nsteps):
			# (opts, args) = options.parser.parse_args()
			options.update(step)
			TestCtx=TractionCtx_Class(options)
			print "Computation id: %s"%TestCtx.jid
			
			runcircle(options, TestCtx)
			# Traction.pproc(options, TestCtx)
			
	elif opts.comptype=='single':
		TestCtx=TractionCtx_Class(options)
		runcircle(options, TestCtx)
		print "Computation id: %s"%TestCtx.jid
		
		
	print "Done with computation with id: %s"%TestCtx.jid
	# Traction.pproc(options, TestCtx)
	
if __name__ == '__main__':
	main()

