#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import glob
import pandas as pd
import sys

#import warnings
#warnings.simplefilter(action='ignore', category=FutureWarning)

c = 299792.458 # km/s

wh = '_lamost_martin_oldgrid_rv2'
run = 'run{}'.format(wh)

#lsftype = 'sspp'
lsftype = 'ferre'
#lsftype = 'aaint2024'

#ipfversion = '_isologg'
ipfversion = ''

#snrlim = 10
#erv = '_{}'.format(snrlim)

erv = ''

flr = False
rviteration = True
newgrid = False


data = pd.read_csv('/Users/ankearentsen/Documents/Werk/lamost/spec_martin/lamost_martin.csv')
#print(len(data))
data = data.sort_values(by='combined_snrg', ascending=False)
data = data.drop_duplicates(subset='combined_gaia_source_id', keep='first')
#print(len(data))
#sys.exit()

### ------------------------------------------------
if rviteration == True:

    rvs = pd.read_csv( '/Users/ankearentsen/Documents/Werk/lamost/spec_martin/ferrespec/ferre-rvs.csv', header=None, names=['name','rv','rv_err'])
    
    rvs[['source_id', 'obsid']] = rvs['name'].str.split('_', expand=True)
    rvs['source_id'] = rvs.source_id.astype(np.int64)

    g1 = data.merge(rvs, left_on='combined_gaia_source_id', right_on='source_id') #, how='left')

    good = g1
else:
    good = data


if ipfversion == '_isologg':
    isologgs = pd.read_csv('/Users/ankearentsen/Documents/Werk/CEMP/XP/Sarah/lamost_cemp_isologg.csv')
#    isologgs = pd.read_csv('/Users/ankearentsen/Documents/Werk/CEMP/XP/Sarah/lamost_normal_isologg.csv')
    
    g2 = good.merge(isologgs, on='source_id')
    
    
    good = g2

#tryname = 3949843210360996224
#good = good[good.combined_gaia_source_id == tryname]



#### TAKE ONLY 10% good spec FOR TESTING!!!!

#good = good[good.combined_snrg > 20].sample(frac=0.1, random_state=1)
#
#
#nms =  [2650009561959801216]

#good = good[good.source_id.isin(nms)]

print('good len = ', len(good))

#######

special = False

print('Total number of stars: {}'.format(len(good)))

#sys.exit()


ipfs = pd.DataFrame()
frds = pd.DataFrame()
wavs = pd.DataFrame()
errs = pd.DataFrame()
orms = pd.DataFrame()
lsfs = pd.DataFrame()

lsn = []

wls = []
for naam in good['combined_gaia_source_id']:
     
     sub = good[(good['combined_gaia_source_id'] == naam)]
     
     print('{}_{}'.format(naam, sub.combined_obsid.values[0]))
     
          

     file = glob.glob( '/Users/ankearentsen/Documents/Werk/lamost/spec_martin/spec-{}-{}_*{}-*{}.fits'.format(sub.combined_lmjd.values[0], sub.combined_planid.values[0], sub.combined_spid.values[0], sub.combined_fiberid.values[0]))[0]


     with pyfits.open(file) as fits:
#        print(fits[1].header)
        spec = fits[1].data['FLUX'][0]
        ivar = fits[1].data['IVAR'][0]
        wl_vac = fits[1].data['WAVELENGTH'][0]
        orm = fits[1].data['ORMASK'][0]
        head = fits[0].header
        x = np.arange(0,head['NAXIS1'])
        wl_head = 10**(head['CRVAL1'] + x*head['CD1_1'])


     ## Convert vacuum to air wavelengths (FROM VACTOAIR.pro)
     sigma2 = (10000./wl_vac)**2.
     fact = 1. +  5.792105e-2/(238.0185 - sigma2) + 1.67917e-3/(57.362 - sigma2)
     wl_air = wl_vac/fact

#     wl_air = wl_vac

     ## Correct for radial velocity
     if (sub.combined_z.values[0] == -9999) | (sub.combined_z.values[0] > 0.5):
        v = 0
     else:
        v = -sub.combined_z.values[0] * c
     wl = wl_air*np.sqrt((1+v/c)/(1-v/c))

     if rviteration == True:
         v2 = -sub.rv.values[0]
#         print(v2)
         if np.abs(v2) < 200:
            wl = wl*np.sqrt((1+v2/c)/(1-v2/c))
#            print('good rv')


     ## Limit to 3700-5250 ang, roughly 1527 pixels  (1730 pix for 5510 ang)  (3765 is for 8805 ang)
#     i1 = np.searchsorted(wl, 3900)
#     i2 = np.searchsorted(wl, 3900) + 3540
     i1 = np.searchsorted(wl, 3850)
     i2 = np.searchsorted(wl, 3850) + 1700
     spec_ferre = spec[i1:i2]
     ivar_ferre = ivar[i1:i2]
     wl_ferre = wl[i1:i2]
     orm_ferre = orm[i1:i2]

     ## Filter file
     ## 0 = excluded, 1 = included // I need to invert the LAMOST array because FERRE interprets it opposite
     orm_ferre = pd.Series(~orm_ferre.astype(bool)).astype(float)

     if flr == True:
         ## Downweight the rest of the spectrum compared to Ca H&K
         orm_ferre[(orm_ferre == 1) & (wl_ferre > 4000)] = 0.1



#     print(wl_ferre)


     # replace zero values in errors, and set error to high for negative SPEC values
     if special == True:
         ivar_ferre[t] =  1e-10

     ivar_ferre2 = np.where(ivar_ferre==0, 1e-10, ivar_ferre)
     ivar_ferre2 = np.where(spec_ferre<=0, 1e-10, ivar_ferre2)
     err_ferre = 1/np.sqrt(ivar_ferre2)

     # replace negative values in SPEC
     spec_ferre = np.where(spec_ferre<=0, np.nanmedian(spec_ferre), spec_ferre)

#     print(np.where(spec_ferre<0), np.nanmedian(spec_ferre))

     ### RESCALE ERRORS TO TEST
#     SNR = np.median(spec_ferre/err_ferre)
#     if SNR > snrlim:
#        err_ferre = err_ferre * (SNR/snrlim)

#     print('SNR = {}'.format(SNR))


#
#     plt.plot(wl_ferre, spec_ferre/err_ferre)
##     plt.plot(wl_ferre, err_ferre)
#     plt.show()



#     wl_ferre2 = wl2[i1:i2]
#
#     if naam == tryname:
#
#         def getspec(file):
#             fits = pyfits.open(file)
#             spec = fits[0].data
#             head = fits[0].header
#             print(spec)
#
#             x = np.arange(0,head['NAXIS1'])
#             wave = head['CRVAL1'] + x*head['CDELT1']
#
#             return spec, wave
#
#         spec_syn, wl_syn = getspec('/Users/ankearentsen/Documents/Werk/CEMP/XP/lamost_allLucey/ferrespec/spec_{}_ferre-syn.fits'.format(tryname))
#
#         spec_fit, wl_fit = getspec('/Users/ankearentsen/Documents/Werk/CEMP/XP/lamost_allLucey/ferrespec/spec_{}_ferre-obs.fits'.format(tryname))
#
#         plt.figure(figsize=(15,4))
#         plt.plot(wl_ferre, spec_ferre/np.nanmedian(spec_ferre))
##         plt.plot(wl_ferre2, spec_ferre/np.nanmedian(spec_ferre))
#         plt.plot(wl_fit, spec_fit/np.nanmedian(spec_fit), color='magenta')
#         plt.plot(wl_syn, spec_syn/np.nanmedian(spec_syn), color='black')
#         plt.show()

     if special == True:
         t = np.where(spec_ferre > 500)
         spec_ferre[t] = 350

#     print(wl_ferre)


    ## ipf file: name + param
     if ipfversion == '_isologg':
        logg_i = sub.logg.values[0]
     else:
        logg_i = 2.0

     par = pd.Series({'name': '{}_{}'.format(naam, sub.combined_obsid.values[0]), 'feh': -1.0, 'cfe': 0, 'teff': 4800, 'logg': logg_i})

     ipfs = pd.concat([ipfs, pd.DataFrame(par).T], ignore_index=True)

     ## frd file: spectrum
     frds = pd.concat([frds, pd.DataFrame(pd.Series(spec_ferre.byteswap().newbyteorder())).T], ignore_index=True)

     ## wav file: wavelengths
     wavs = pd.concat([wavs, pd.DataFrame(pd.Series(wl_ferre)).T], ignore_index=True)

     # err file: errors
     errs = pd.concat([errs, pd.DataFrame(pd.Series(err_ferre)).T], ignore_index=True)


#     print(np.median(spec_ferre[wl_ferre < 3800]/errs[wl_ferre < 3800]))


     if special == True:
         plt.plot(wl_ferre,spec_ferre)
         plt.plot(wl_ferre,err_ferre)
         plt.show()

     ## flr file: filter
     orms = pd.concat([orms, pd.DataFrame(orm_ferre).T], ignore_index=True)


ipfs = ipfs[['name','feh','cfe','teff','logg']]


ipfs.to_csv('{}.ipf'.format(run), header=None, index=None, sep=' ')
frds.to_csv('{}.frd'.format(run), header=None, index=None, sep=' ')
wavs.to_csv('{}.wav'.format(run), header=None, index=None, sep=' ')
errs.to_csv('{}{}.err'.format(run,erv), header=None, index=None, sep=' ')
orms.to_csv('{}.flr'.format(run), header=None, index=None, sep=' ')

print('Total number of stars: {}'.format(len(good)))

### -------------------------------------
### for Model spec: calculate LSF
### -------------------------------------
#
#
#if lsftype == 'ferre':
#    Npix = 32968
#    x = np.arange(0,Npix)
#    wave = 10**(3.55022035572257 + x*1.36436155109925e-05)
#    filesub = ''
#    R1 = 10000
#elif lsftype == 'sspp':
#    Npix = 28001
#    x = np.arange(0,Npix)
#    wave = 3000.0 + x*0.25
#    filesub = '_{}'.format(lsftype)
#    R1 = 10000
#elif lsftype == 'aaint2024':
#    Npix = 340467
#    x = np.arange(0,Npix)
#    wave = 3000.00011 + x*0.020559999999932188
#    filesub = '_{}'.format(lsftype)
#    R1 = 100000
#
##R1 = 10000
#R2 = 1800       #LAMOST R
#
#fwhm = []
#for i in np.arange(0,Npix-1):
#    # FWHM in Angstrom
#    fwhm1 = wave[i]/R1
#    fwhm2 = R1/R2 * fwhm1         #angstrom
#    dwl = (wave[i+1] - wave[i])   #angstrom / pixel
#
#    fwhm2_pix = fwhm2/dwl         #pixel
#    fwhm.append(fwhm2_pix)
#
## add one more value for the final pixel (same as the one before)
#fwhm.append(fwhm2_pix)
#
## for log wavelength, the FWHM is constant, so adopt a wl single value
#if lsftype == 'ferre':
#    lsf = fwhm[0]
#elif lsftype == 'sspp':
#    lsf = fwhm
#elif lsftype == 'aaint2024':
#    lsf = fwhm
#
#print(fwhm[0], fwhm[-1])
##print(lsf)
##
#with open ('{}{}.lsf'.format(run,filesub), 'w') as f:
#    f.write('{}'.format(lsf))
