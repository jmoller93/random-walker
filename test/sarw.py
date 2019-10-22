#!/env/python

import sys,os,glob
sys.path.append('../build/')
import walkers

nMonomers = 100
bLength   = 10

sarw = walkers.Walker(nMonomers,bLength)
sarw.chain_growth()
print(sarw.get_coord())
