#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import glob
import pandas as pd
import sys
import os
from funcs import mkfls_lam, writef, tofits, getrvlist, rvpyraf, plotspec, runlam, dorv

### ----------------------------------------
### Input files/parameters

home = '/Users/ankearentsen/Documents/Werk/code/ferre-master/explore/'
gridfile = '../f_crump4f_Rlamost2.dat'
specfol = '/Users/ankearentsen/Documents/Werk/lamost/spec_martin/'
data = pd.read_csv('{}lamost_martin.csv'.format(specfol))

run0 = 'run_lamost_martin_oldgrid_test'
newgrid = False
rviteration = False
vs = ''

### ----------------------------------------
### set some further parameters based on that & filter LAMOST input file

outfol = 'files/'
outfolferre = '{}ferrespec/'.format(outfol)

os.makedirs(outfol, exist_ok=True)
os.makedirs(outfolferre, exist_ok=True)

rvfile0 = '{}ferre-rvs.csv'.format(outfolferre)
    
# sort LAMOST entries by SNR and keep only the highest
data = data.sort_values(by='combined_snrg', ascending=False)
data = data.drop_duplicates(subset='combined_gaia_source_id', keep='first')


### ----------------------------------------
### CREATE FILES & RUN FERRE

rviteration = False
rvfile = None
run = run0

runlam(data, specfol, outfol, run, gridfile, rviteration, rvfile, newgrid)


#### ----------------------------------------
#### To get RV correction

dorv(run, outfol, outfolferre, rvfile0, vs, newgrid)


#### ----------------------------------------
#### Do it with again with RV correction

rviteration = True
rvfile = rvfile0
run = run0+'_rv2'

runlam(data, specfol, outfol, run, gridfile, rviteration, rvfile, newgrid)


## ----------------------------------------
## ----------------------------------------
## ----------------------------------------
### Plotting


p1 = pd.read_csv('{}ferre_{}.csv'.format(outfol,run0))
plt.scatter(p1.feh_ferre, p1.cfe_ferre, label='no RVcor')
p2 = pd.read_csv('{}ferre_{}.csv'.format(outfol,run))
plt.scatter(p2.feh_ferre, p2.cfe_ferre, label='RVcor')
plt.legend()
plt.xlabel('[Fe/H]')
plt.ylabel('[C/Fe]')
plt.show()


plotspec(['3650671417907498752_132312046'], run0+'_rv2', outfol, vs='', newgrid=False)
