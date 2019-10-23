#!/env/python

import sys,os,glob
import numpy as np
sys.path.append('../build/')
import walkers

# Variable init
nMonomers = 30
bLength   = 10
tol = 9
nTrials = 10
nReps   = 10
rg = np.zeros(nReps)

# Run the walker
for i in range(nReps):
    sarw = walkers.Walker(nMonomers,bLength)
    sarw.chain_growth(tol=tol,max_trials=nTrials)
    #print(sarw.get_coord())
    rg[i] = np.sqrt(sarw.get_rg())
    print(i)

# Print results
print("Average Rg = %.3f\nVariance of Rg = %.3f" % (np.mean(rg),np.std(rg)**2))
