#!/env/python

import sys,os,glob
import numpy as np
sys.path.append('../build/')
import walkers

# Variable init
nMonomers = 100
bLength   = 10
tol = bLength/2.0
nTrials = 200
nReps   = 500
rg = np.zeros(nReps)

# Run the walker
sarw = walkers.Walker(nMonomers,bLength)
for i in range(nReps):
    sarw.chain_growth(tol=tol,max_trials=nTrials)
    rg[i] = np.sqrt(sarw.get_rg())

# Print results
print("Average Rg = %.3f\nVariance of Rg = %.3f" % (np.mean(rg),np.std(rg)**2))
