#! /usr/bin/env python


def Data_processing_python(data_name,R_int,R_out,W,X0):
	import numpy as np
	import scipy as sp
	Max_elements=np.genfromtxt(data_name, usecols=(0), skip_header=0, skip_footer=-1 , dtype=float, unpack=True)
	data = np.genfromtxt(data_name, usecols=(0,1), skip_header=1, skip_footer=Max_elements[0] )
	np.matrix(data)
	g = np.genfromtxt(data_name, usecols=(0), skip_header=int(Max_elements[0]+1.) , dtype=float, unpack=True)
	from scipy.interpolate import LinearNDInterpolator
	f_cart=sp.interpolate.LinearNDInterpolator(data, g, rescale=False)
	r = np.linspace(R_int, R_out, 4000)
	o = np.linspace(-W*np.pi,W*np.pi,4000)
	x = np.linspace(X0, X0, 4000)
	f = f_cart(r[:,np.newaxis]*np.cos(o[np.newaxis,:])+x[:,np.newaxis],r[:,np.newaxis]*np.sin(o[np.newaxis,:]))*(r[:,np.newaxis])
	G_theta=sp.trapz(sp.trapz(f, o[np.newaxis,:], axis=1), r, axis=0)
	return G_theta


if __name__ == "__main__":
	import sys
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

	G_theta=[]
	for Time in range(int(Param['time_min']),int(Param['time_max'])):
		G=Data_processing_python('Data_G_theta_'+str(Time)+'.0.0.exvol_var',float(Geometry['R_inn']),float(Geometry['R_out']),float(Geometry['omega']),float(Geometry['X_center'][Time]))
		print 'evaluating G : {0:.3f}'.format(G)
		G_theta.append(G)
	print "G theta   :"
	print G_theta

	#save data un G.txt
	Data_ener=np.loadtxt(rootdir+'/'+str(Param['prefix'])+'.ener')
	Complete = np.zeros(len(Data_ener[:,1])-len(G_theta))
	G_theta.extend(Complete)
	Tabl=np.concatenate(([Data_ener[:,0]],[Data_ener[:,1]],[Data_ener[:,2]],[Data_ener[:,3]],[Data_ener[:,4]],[Data_ener[:,5]],[G_theta]),axis=0).T
	np.savetxt('G.txt', Tabl, fmt='%1.4e',delimiter='\t\t',header='step \t\t\t load \t\t\t elastic energy \t work \t\t\t surface energy \t total energy \t\t G' )
#remove visit data out
#os.system('rm Data_G_theta_*')

exit()