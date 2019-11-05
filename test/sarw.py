#!/env/python

import sys,os,glob
import numpy as np
sys.path.append('../build/')
import walkers

# Variable init
nMonomers = 5
bLength   = 10
tol = 9
nTrials = 10
nReps   = 10
rg = np.zeros(nReps)

# Run the walker
#for i in range(nReps):
dists = np.asarray([10,20,30,40,50])
sarw = walkers.Walker()
sarw.chain_growth(dists)
for i in range(nTrials):
    excludeBool = sarw.pivot()
    if excludeBool:
        idx += 1
print(sarw.get_coord())
#print(sarw.get_coord())
#rg[i] = np.sqrt(sarw.get_rg())
#print(i)

# Print results
print("Average Rg = %.3f\nVariance of Rg = %.3f" % (np.mean(rg),np.std(rg)**2))
