#!/env/python

import sys,os,glob
sys.path.append('../build/')
import walkers

nMonomers = 1000
bLength   = 10
tol = bLength/2.0
nTrials = 200

sarw = walkers.Walker(nMonomers,bLength)
sarw.chain_growth(tol=tol,max_trials=nTrials)
print(sarw.get_coord())
