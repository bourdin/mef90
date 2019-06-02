def plot(energies, showwork=False):
  import matplotlib.pyplot as plt
  plt.plot(energies[:,1], energies[:,2], 'r-', label='Elastic energy')
  if showwork:
    plt.plot(energies[:,1], energies[:,3], 'k-', label='External Forces')
  plt.plot(energies[:,1], energies[:,4], 'g-', label='Surface energy')
  plt.plot(energies[:,1], energies[:,5], 'b-', label='Total energy', lw=2)
  plt.grid()
  plt.legend(loc=0)
  plt.xlabel('t')
  plt.ylabel('Energy')
  plt.title('Energies vs normalized time')
  ###
  return 0

def getlaststep2(fname):
  ### open file
  f=open(fname)
  ### Read last line in a string
  lastline = f.readlines()[-1]
  laststep = lastline.rsplit()[0] 
  return(int(laststep))

def getlaststep(fname):
  import numpy as np
  ### open file
  ener = np.readtxt(fname)
  print(ener[-1])
  return ener[-1][0]

def save(fname, energies):
  import numpy as np
  np.savetxt(fname, energies, fmt='%7d  %13.5E %13.5E %13.5E %13.5E %13.5E %13.5E')

def fixBT(energies, laststep=None):
  import numpy as np
  ###
  if laststep == None:
    laststep = int(energies[-1,0])
  else:
    laststep = min(int(energies[-1,0]), int(laststep))
  maxstep  = energies.shape[0]
  ###
  energiesBT = np.zeros([laststep,energies.shape[1]])
  ###
  i = 0
  while True:
    step = energies[i,0]
    energiesBT[step-1,:] = energies[i,:]
    i += 1
    if step == laststep or i >= maxstep:
      break
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

def ReadCompositeEnergiesBT(prefix, laststep=None):
  import numpy as np
  ###
  toughness   = np.loadtxt(prefix+'.CST', skiprows=1, usecols=[1])
  all_energies = []
  for blk in range(toughness.shape[0]):
    blkfile="%s-%.4i.enerblk" % (prefix, blk+1) 
    all_energies.append(np.loadtxt(blkfile))
  ###
  ### Compute last step and larger step
  ###
  if laststep == None:
    laststep = int(all_energies[0][-1,0])
  else:
    lasstep = min(int(laststep), int(all_energies[0][-1,0]))
  maxstep = all_energies[0].shape[0]
  ###
  ### Initialize all_energies_BT
  ###
  all_energiesBT=[]
  for e in all_energies:
    all_energiesBT.append(np.zeros([int(laststep),e.shape[1]]))
  ###
  ### Remove redundant computations
  ###
  i = 0
  while True:
    step = all_energies[0][i,0]
    for (energiesBT, energies) in zip(all_energiesBT, all_energies):
      energiesBT[step-1,:] = energies[i,:]
    i += 1
    if step == int(laststep) or i >= maxstep:
      break
  ###
  ### Extract Elastic Energy, length, time, toughness 
  ##
  Eel = np.sum(e[:,2] for e in all_energiesBT)
  l   = np.sum(e[:,4]/k for (e,k) in zip(all_energiesBT, toughness))
  t   = all_energiesBT[0][:,1]
  return Eel, l, t, toughness

  