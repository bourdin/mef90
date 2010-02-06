#!/usr/bin/env python 

def plotener(energies):
  import matplotlib.pyplot as plt

  plt.plot(energies[:,1], energies[:,2], 'r', label='Elastic energy')
  plt.plot(energies[:,1], energies[:,4], 'g', label='Surface energy')
  plt.plot(energies[:,1], energies[:,5], 'b', label='Total energy', lw=2)
  plt.grid()
  plt.legend(loc=0)
  plt.xlabel('t')
  plt.ylabel('Energy')
  plt.title('Energies vs normalized time')
  
  return 0

def plotener_withforces(energies):
  import matplotlib.pyplot as plt

  plt.plot(energies[:,1], energies[:,2], 'r', label='Elastic energy')
  plt.plot(energies[:,1], energies[:,3], 'k', label='External Forces')
  plt.plot(energies[:,1], energies[:,4], 'g', label='Surface energy')
  plt.plot(energies[:,1], energies[:,5], 'b', label='Total energy', lw=2)
  plt.grid()
  plt.legend(loc=0)
  plt.xlabel('t')
  plt.ylabel('Energy')
  plt.title('Energies vs normalized time')

  return 0
  
def energetlaststep(fname):
### open file
  f=open(fname)
### Read last line in a string
  lastline = f.readlines()[-1]
  laststep = lastline.rsplit()[0] 
  return(int(laststep))
  
def enersave(fname, energies):
  import numpy as np
  np.savetxt(fname, energies, fmt='%7d  %13.5E %13.5E %13.5E %13.5E %13.5E')
  
def enerfixBT(energies):
  import numpy as np
  
  laststep = energies[-1,0]
  maxstep  = energies[:,0].max()
  
  energiesBT = np.zeros([laststep,energies.shape[1]])
  for i in range(energies.shape[0]):
    step = energies[i,0]
    if step <= laststep:
      energiesBT[step-1,:] = energies[i,:]
  return energiesBT