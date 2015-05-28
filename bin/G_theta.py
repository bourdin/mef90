#visit -cli -nowin -s G_theta.py 
#visit -cli -nowin -np 12 -l srun -s G_theta.py

#! /usr/bin/env python
from visit import *


def Export_data_from_visit(rootdir,prefix,R_int,R_out,X0,Y0,Time,E0,nu,AT_model,k):
	import json
	import os
	import os.path
	import shutil
	import math

	MyDatabase = os.path.join(rootdir+'/'+prefix+'_out.gen')
	print MyDatabase
	OpenDatabase(MyDatabase)

	DefineScalarExpression("X", 'coord(Mesh)[0]-'+str(X0)+'')
	DefineScalarExpression("Y", 'coord(Mesh)[1]-'+str(Y0)+'')
	DefineScalarExpression("R", "sqrt(X*X+Y*Y)")
	DefineScalarExpression("k", ''+k+'')
	DefineScalarExpression("E0", ''+E0+'')
	DefineScalarExpression("v", ''+nu+'')
	print AT_model
	if AT_model=='ATk':
		DefineScalarExpression("E", 'E0*(1-Damage)^2./(  1+(k-1)*( 1 - (1-Damage)^2 )  )')
	elif AT_model=='AT1' or AT_model=='AT2':
		DefineScalarExpression("E", 'E0*(1-Damage)^2.')
	DefineVectorExpression("Displacement_2D", "{Displacement_X,Displacement_Y}")
	DefineTensorExpression("EPS_2D", "{{gradient(Displacement_X)[0],0.5*(gradient(Displacement_X)[1]+gradient(Displacement_Y)[0])},{0.5*(gradient(Displacement_X)[1]+gradient(Displacement_Y)[0]),gradient(Displacement_Y)[1]}} ")
	DefineVectorExpression("Theta_vec", '{((R-'+R_out+')/('+R_int+'-'+R_out+')),0}')
	DefineTensorExpression("GRAD_Theta_vec", "{{gradient(Theta_vec[0])[0],gradient(Theta_vec[0])[1]},{gradient(Theta_vec[1])[0],gradient(Theta_vec[1])[1]}}")
	DefineScalarExpression("G_theta", "trace(STRESS_2D*(GRAD_U*GRAD_Theta_vec))-.5*trace(STRESS_2D*EPS_2D)*divergence(Theta_vec)")
	DefineTensorExpression("GRAD_U", "{{gradient(Displacement_2D[0])[0],gradient(Displacement_2D[0])[1]},{gradient(Displacement_2D[1])[0],gradient(Displacement_2D[1])[1]}}")
	DefineTensorExpression("STRESS_2D", "E/((1+v))*(EPS_2D+v/(1-2*v)*trace(EPS_2D)*{{1,0},{0,1}})")

	SetTimeSliderState(Time)
	AddPlot("Pseudocolor", "G_theta")
	DrawPlots()
	ExportDBAtts = ExportDBAttributes()
	ExportDBAtts.db_type = "ExtrudedVol"
	ExportDBAtts.filename = 'Data_G_theta_'+str(Time)+''
	ExportDBAtts.dirname = ""
	ExportDBAtts.variables = ('G_theta')
	#ExportDBAtts.opts.types = () 
	ExportDatabase(ExportDBAtts)
	DeleteActivePlots()
	CloseDatabase(MyDatabase)



if __name__ == "__main__":
	import sys
	import os.path
	import numpy as np
	import json
	import os
	import os.path
	import shutil
	import math

	rootdir = os.getenv("PWD")

	if os.path.exists(os.path.join(rootdir,'00_INFO.json')):
		json_file = open(os.path.join(rootdir,'00_INFO.json'))
		D = json.load(json_file)
		Geometry = dict(D['Geometry'])
		Mat_prop = dict(D['Mat_prop'])
		Param = dict(D['Param'])
		Param.update()
		json_file.close()

	for Time in range(int(Param['time_min']),int(Param['time_max'])):
		# 
		Export_data_from_visit(str(rootdir),str(Param['prefix']),str(Geometry['R_inn']),str(Geometry['R_out']),float(Geometry['X_center'][Time]),float(Geometry['Y_center']),Time,str(Mat_prop['E0']),str(Mat_prop['nu']),str(Mat_prop['AT_model']),str(Mat_prop['k']))

exit()
