def plot(energies, showwork=False):
  import matplotlib.pyplot as plt
  plt.plot(energies[:,1], energies[:,2], 'r', label='Elastic energy')
  if showwork:
    plt.plot(energies[:,1], energies[:,3], 'k', label='External Forces')
  plt.plot(energies[:,1], energies[:,4], 'g', label='Surface energy')
  plt.plot(energies[:,1], energies[:,5], 'b', label='Total energy', lw=2)
  plt.grid()
  plt.legend(loc=0)
  plt.xlabel('t')
  plt.ylabel('Energy')
  plt.title('Energies vs normalized time')
  ###
  return 0

def getlaststep(fname):
  ### open file
  f=open(fname)
  ### Read last line in a string
  lastline = f.readlines()[-1]
  laststep = lastline.rsplit()[0] 
  return(int(laststep))

def save(fname, energies):
  import numpy as np
  np.savetxt(fname, energies, fmt='%7d  %13.5E %13.5E %13.5E %13.5E %13.5E')

def fixBT(energies):
  import numpy as np
  ###
  laststep = energies[-1,0]
  maxstep  = energies[:,0].max()
  ###
  energiesBT = np.zeros([laststep,energies.shape[1]])
  for i in range(energies.shape[0]):
    step = energies[i,0]
    if step <= laststep:
      energiesBT[step-1,:] = energies[i,:]
    ###
  return energiesBT

def computeG(t, Eel, l):
  import numpy as np
  ###
  G = -(Eel[1:]/np.power(t[1:],2)-Eel[:-1]/np.power(t[:-1],2)) / (l[1:]-l[:-1]) 
  ###
  return G

def ReadCompositeEnergies(prefix, stepmin=None, stepmax=None):
  import numpy as np
  ###
  toughness   = np.loadtxt(prefix+'.CST', skiprows=1, usecols=[1])
  all_energies = []
  for blk in range(toughness.shape[0]):
    blkfile="%s-%.4i.enerblk" % (prefix, blk+1) 
    all_energies.append(np.loadtxt(blkfile))
  ###
  if stepmin == None:
    tmin = 0
  else:
    tmin = int(stepmin)
  if stepmax == None:
    tmax = all_energies[0].shape[0]
  else:
    tmax = int(stepmax)
  ###
  Eel = np.sum(e[tmin:tmax,2] for e in all_energies)
  l   = np.sum(e[tmin:tmax,4]/k for (e,k) in zip(all_energies, toughness))
  t   = all_energies[0][tmin:tmax,1]
  return Eel, l, t, toughness

