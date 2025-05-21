#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import glob
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib
import sys
font = {'size': 14}
matplotlib.rc('font', **font)

#run  = 'run_lamost_martin_newgrid_rv2'
#run0 = 'run_lamost_martin_oldgrid_rv2'

run  = 'run_lamost_martin_oldgrid'
run0 = run


#vs = 'v2'
#vs = '_ssppcool1'
#vs = '_ssppwarm1'
#vs = '_int3al3'
#vs = '_int3al3_flr'
vs = ''

#erv = '_10'
erv = ''

ster =  '131045813048526976_162805073'

if 'newgrid' in run:
    newgrid = True
else:
    newgrid = False
    
print('newgrid = ', newgrid)



obs_orig = pd.read_csv('{}.frd'.format(run0), sep='\s+', header=None)
err = pd.read_csv('{}{}.err'.format(run0,erv), sep='\s+', header=None)
wl = pd.read_csv('{}.wav'.format(run0), sep='\s+', header=None)
#lsfs = pd.read_csv('{}.lsf'.format(run0), sep='\s+', header=None)
obs = pd.read_csv('{}{}.obs'.format(run,vs), sep='\s+', header=None)
mod = pd.read_csv('{}{}.mdl'.format(run,vs), sep='\s+', header=None)
if newgrid == True:
    par = pd.read_csv('{}{}.opf'.format(run,vs), sep='\s+', header=None, names=['name','teff_ferre','logg_ferre','feh_ferre','cfe_ferre','e_teff_ferre','e_logg_ferre','e_feh_ferre','e_cfe_ferre','zero_ferre','SNR_ferre','chi2_ferre'])
else:
    par = pd.read_csv('{}{}.opf'.format(run,vs), sep='\s+', header=None, names=['name','feh_ferre','cfe_ferre','teff_ferre','logg_ferre','e_feh_ferre','e_cfe_ferre','e_teff_ferre','e_logg_ferre','zero_ferre','SNR_ferre','chi2_ferre'])


print('Number of stars: ', len(par), len(wl))


par.to_csv('ferre_{}{}.csv'.format(run,vs), index=None)


#raise SystemExit()

#sub = par[(par.teff_ferre > 5400) & (par.SNR_ferre > 30)]     #[(par.feh_ferre < -3.9) & (par.SNR_ferre > 25)]

#print('number of stars: {}'.format(len(sub)))

sub = par[(par.name == ster)]

print(len(sub))


for i in [sub.index.values[0]]:
    fig = plt.figure(figsize=(12,5))
    gs1 = gridspec.GridSpec(2,2, height_ratios=[1,1], width_ratios=[3,1])
    ax1 = fig.add_subplot(gs1[0,0])
    ax2 = fig.add_subplot(gs1[1,0])
    ax3 = fig.add_subplot(gs1[0,1])
    ax4 = fig.add_subplot(gs1[1,1])

    ster = sub.name.loc[i]
    print(ster)
    
#    print(wl.iloc[i])

    minwl = 3850
    maxwl = 5750

    ## FIT
    ax1.set_title('teff={}, logg={}, feh={}, cfe={}'.format(int(par.teff_ferre.iloc[i]),round(par.logg_ferre.iloc[i],1), round(par.feh_ferre.iloc[i],1), round(par.cfe_ferre.iloc[i],1)))
    ax1.plot(wl.iloc[i], obs.iloc[i], color='black', label='Normalised observation')
    ax1.plot(wl.iloc[i], mod.iloc[i], color='red', label='Normalised model')
    ax1.set_xlim([minwl,maxwl])
    ax1.get_xaxis().set_visible(False)
    ax1.legend(fontsize=12)
    ax1.set_ylim(0,1.7)

    ax3.plot(wl.iloc[i], obs.iloc[i], color='black', label='Normalised observation')
    ax3.plot(wl.iloc[i], mod.iloc[i], color='red', label='Normalised model')
    ax3.set_xlim([3900,4000])
    ax3.get_xaxis().set_visible(False)
    ax3.get_yaxis().set_visible(False)
    ax3.set_ylim(0,1.7)

    ## ORIG + ERR
    ax2.plot(wl.iloc[i], obs_orig.iloc[i], color='black', label='non-normalised observation')
    ax2.plot(wl.iloc[i], err.iloc[i], color='green', label='noise spectrum')
    ax2.set_ylim(0,np.median(obs_orig.iloc[i])+np.median(obs_orig.iloc[i])*0.5)
    ax2.set_xlim([minwl,maxwl])
    ax2.legend(fontsize=12)
    ax2.set_xlabel('Wavelength (angstrom)')

    ax4.plot(wl.iloc[i], obs_orig.iloc[i], color='black', label='non-normalised observation')
    ax4.plot(wl.iloc[i], err.iloc[i], color='green', label='noise spectrum')
    ax4.set_ylim(0,np.median(obs_orig.iloc[i])+np.median(obs_orig.iloc[i])*0.5)
    ax4.set_xlim([3900,4000])
    ax4.get_yaxis().set_visible(False)

    gs1.tight_layout(fig)
    gs1.update(hspace=0, wspace=0.02)

#    plt.savefig('sub/figures/{}.jpg'.format(ster), dpi=150)

    plt.show()
    plt.close()

