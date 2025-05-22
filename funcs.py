#!/usr/bin/env python

import numpy as np
import astropy.io.fits as pyfits
import glob
import pandas as pd
import os.path
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
import sys

c = 299792.458 # km/s

### Inspired by APS pipeline
def mknml(outfol, ipf_file, frd_file, wav_file, err_file, flr_file, root, outroot, sf1, sf2=None, ngrid=1):

    nml={}
    ndim = 4
    nml['NDIM']=ndim
    nml['NOV']=ndim
    nml['INDV']=' '.join(map(str,np.arange(ndim)+1))
    nml['SYNTHFILE(1)'] = "'"+sf1+"'"
    if ngrid == 2:
        nml['SYNTHFILE(2)'] = "'"+sf2+"'"
    nml['PFILE'] = "'"+outfol+ipf_file+"'"
    nml['FFILE'] = "'"+outfol+frd_file+"'"
    nml['WFILE'] = "'"+outfol+wav_file+"'"
    nml['ERFILE'] = "'"+outfol+err_file+"'"
    nml['FILTERFILE'] = "'"+outfol+flr_file+"'"
    nml['SFFILE'] = "'"+outfol+outroot+".obs"+"'"
    nml['OPFILE'] = "'"+outfol+outroot+".opf"+"'"
    nml['OFFILE'] = "'"+outfol+outroot+".mdl"+"'"
    nml['ERRBAR']=1
    nml['ALGOR']=3
    nml['INTER']=3
    nml['WINTER']=2
    nml['CONT']=3
    nml['NCONT']=30
    nml['NTHREADS']=5

    return nml


### From APS pipeline
def writenml(nml,nmlfile='input.nml',path=None):
    if path is None: path='./'
    f=open(os.path.join(path,nmlfile),'w')
    f.write('&LISTA\n')
    for item in nml.keys():
        f.write(str(item))
        f.write("=")
        f.write(str(nml[item]))
        f.write("\n")
    f.write(" /\n")
    f.close()
    print('Wrote {}'.format(nmlfile))
    return None


def mkfls_lam(data, specfol, outfol, run, gridfile, rviteration=False, rvfile=None, newgrid=False):
    '''
    Parameters
    ----------
    data : pandas dataframe
          summary file from LAMOST archive for the downloaded spectra,
          incl. combined_gaia_source_id
    specfol: str
          folder where the LAMOST spectra are stored
    outfol: str
          folder where to save the FERRE files
    run: str
          name for this FERRE run, used to refer to a unique analysis
    rviteration: bool
          when True a radial velocity correction will be used on top of
          the default LAMOST combined_z value, from rvfile
    newgrid: bool
          set to True if using a grid with order [teff, logg, feh, cfe]
          rather than the default [feh, cfe, teff, logg]
    rvfile: str
          path to a file with radial velocity corrections, with columns
          name, rv and rv_err
    
    Returns
    -------
    wavelength: numpy array
       corresponding air wavelengths in AA (keeps vacuum for lambda<=2000. A)
       
    '''
    

    ### ------------------------------------------------
    if rviteration == True:

        rvs = pd.read_csv(rvfile, header=None, names=['name','rv','rv_err'])
        
        a = rvs.name.str.split('/', expand=True).iloc[:, -1].str.split('_', expand=True)
        source_ids = a.iloc[:,-2]
        obsids = a.iloc[:,-1]
                
        rvs['obsid'] = obsids.astype(np.int64)
        rvs['source_id'] = source_ids.astype(np.int64)

        g1 = data.merge(rvs, left_on='combined_gaia_source_id', right_on='source_id')

        good = g1
    else:
        good = data

    #######

    print('Total number of stars: {}'.format(len(good)))


    ipfs = pd.DataFrame()
    frds = pd.DataFrame()
    wavs = pd.DataFrame()
    errs = pd.DataFrame()
    orms = pd.DataFrame()
    lsfs = pd.DataFrame()

    lsn = []

    wls = []
    m = 1
    for naam in good['combined_gaia_source_id']:
         
         sub = good[(good['combined_gaia_source_id'] == naam)]
         
         print('{} {}_{}'.format(m, naam, sub.combined_obsid.values[0]))
         
         file = glob.glob( '{}spec-{}-{}_*{}-*{}.fits'.format(specfol, sub.combined_lmjd.values[0], sub.combined_planid.values[0], sub.combined_spid.values[0], sub.combined_fiberid.values[0]))[0]


         with pyfits.open(file) as fits:
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


         ## Limit wavelength range, make sure they're all the same length
         ## 1700 pixels results in ~3850-5690 A
         i1 = np.searchsorted(wl, 3850)
         i2 = np.searchsorted(wl, 3850) + 1700
         spec_ferre = spec[i1:i2]
         ivar_ferre = ivar[i1:i2]
         wl_ferre = wl[i1:i2]
         orm_ferre = orm[i1:i2]

         ## Filter file
         ## 0 = excluded, 1 = included // I need to invert the LAMOST array because FERRE interprets it opposite
         orm_ferre = pd.Series(~orm_ferre.astype(bool)).astype(float)

         ivar_ferre2 = np.where(ivar_ferre==0, 1e-10, ivar_ferre)
         ivar_ferre2 = np.where(spec_ferre<=0, 1e-10, ivar_ferre2)
         err_ferre = 1/np.sqrt(ivar_ferre2)

         # replace negative values in SPEC
         spec_ferre = np.where(spec_ferre<=0, np.nanmedian(spec_ferre), spec_ferre)

        ## ipf file: name + param
         par = pd.Series({'name': '{}_{}'.format(naam, sub.combined_obsid.values[0]), 'feh': -1.0, 'cfe': 0, 'teff': 4800, 'logg': 2.0})

         ipfs = pd.concat([ipfs, pd.DataFrame(par).T], ignore_index=True)

         ## frd file: spectrum
         frds = pd.concat([frds, pd.DataFrame(pd.Series(spec_ferre.byteswap().newbyteorder())).T], ignore_index=True)

         ## wav file: wavelengths
         wavs = pd.concat([wavs, pd.DataFrame(pd.Series(wl_ferre)).T], ignore_index=True)

         # err file: errors
         errs = pd.concat([errs, pd.DataFrame(pd.Series(err_ferre)).T], ignore_index=True)

         ## flr file: filter
         orms = pd.concat([orms, pd.DataFrame(orm_ferre).T], ignore_index=True)
         
         m+=1


    ipfs = ipfs[['name','feh','cfe','teff','logg']]


    ipfs.to_csv('{}{}.ipf'.format(outfol,run), header=None, index=None, sep=' ')
    frds.to_csv('{}{}.frd'.format(outfol,run), header=None, index=None, sep=' ')
    wavs.to_csv('{}{}.wav'.format(outfol,run), header=None, index=None, sep=' ')
    errs.to_csv('{}{}.err'.format(outfol,run), header=None, index=None, sep=' ')
    orms.to_csv('{}{}.flr'.format(outfol,run), header=None, index=None, sep=' ')

    print('Total number of stars in FERRE spec files: {}'.format(len(good)))
    
    root = '{}'.format(run)
    ipf_file = '{}.ipf'.format(root)
    frd_file = '{}.frd'.format(root)
    wav_file = '{}.wav'.format(root)
    err_file = '{}.err'.format(root)
    orm_file = '{}.flr'.format(root)

    outroot = root
    
    nml = mknml(outfol, ipf_file, frd_file, wav_file, err_file, orm_file, root, outroot, ngrid=1, sf1=gridfile)
    writenml(nml,nmlfile='input.nml_{}'.format(outroot), path=outfol)

    return None

def writef(run, outfol, vs='', newgrid=False):

    if newgrid == True:
        par = pd.read_csv('{}{}{}.opf'.format(outfol,run,vs), sep='\s+', header=None, names=['name','teff_ferre','logg_ferre','feh_ferre','cfe_ferre','e_teff_ferre','e_logg_ferre','e_feh_ferre','e_cfe_ferre','zero_ferre','SNR_ferre','chi2_ferre'])
    else:
        par = pd.read_csv('{}{}{}.opf'.format(outfol,run,vs), sep='\s+', header=None, names=['name','feh_ferre','cfe_ferre','teff_ferre','logg_ferre','e_feh_ferre','e_cfe_ferre','e_teff_ferre','e_logg_ferre','zero_ferre','SNR_ferre','chi2_ferre'])

    par.to_csv('{}ferre_{}{}.csv'.format(outfol, run,vs), index=None)
    
    print('wrote {}ferre_{}{}.csv'.format(outfol, run,vs))
    
    return None


def tofits(run, outfol, outfolferre, vs='', newgrid=False):

    obs_orig = pd.read_csv('{}{}.frd'.format(outfol,run), sep='\s+', header=None)
    err = pd.read_csv('{}{}.err'.format(outfol,run), sep='\s+', header=None)
    wl = pd.read_csv('{}{}.wav'.format(outfol,run), sep='\s+', header=None)
    obs = pd.read_csv('{}{}{}.obs'.format(outfol,run,vs), sep='\s+', header=None)
    mod0 = pd.read_csv('{}{}{}.mdl'.format(outfol,run,vs), sep='\s+', header=None)
    if newgrid == False:
        par = pd.read_csv('{}{}{}.opf'.format(outfol,run,vs), sep='\s+', header=None, names=['name','feh_ferre','cfe_ferre','teff_ferre','logg_ferre','e_feh_ferre','e_cfe_ferre','e_teff_ferre','e_logg_ferre','zero_ferre','SNR_ferre','chi2_ferre'])
    else:
        par = pd.read_csv('{}{}{}.opf'.format(outfol,run,vs), sep='\s+', header=None, names=['name','teff_ferre','logg_ferre','feh_ferre','cfe_ferre','e_teff_ferre','e_logg_ferre','e_feh_ferre','e_cfe_ferre','zero_ferre','SNR_ferre','chi2_ferre'])

    mod1 = mod0.replace(to_replace='Infinity', value=0)
    mod = mod1.apply(pd.to_numeric)

    pad = outfolferre
    
    os.makedirs(pad, exist_ok=True)

    for k in par.index.values:
        name = par.name[k]
        print(name)

        fl_obs = obs.iloc[k]
        fl_syn = mod.iloc[k]
        wl_both =wl.iloc[k]

        wl_goal = np.arange(wl_both[0], wl_both[len(wl_both)-1]+0.000001, (wl_both[len(wl_both)-1] - wl_both[0])/(len(wl_both)))

        fl_obs_int = pd.Series(np.interp(wl_goal, wl_both, fl_obs))
        fl_syn_int = pd.Series(np.interp(wl_goal, wl_both, fl_syn))

        ### OBS
        hdul = pyfits.PrimaryHDU(fl_obs_int)

        head = hdul.header
        head['CRPIX1'] = 1
        head['CRVAL1'] = wl_goal[0]
        head['CDELT1'] = wl_goal[1]-wl_goal[0]

        hdul.writeto('{}spec_{}_ferre-obs.fits'.format(pad,name), overwrite=True)

        ### SYN
        hdul = pyfits.PrimaryHDU(fl_syn_int)

        head = hdul.header
        head['CRPIX1'] = 1
        head['CRVAL1'] = wl_goal[0]
        head['CDELT1'] = wl_goal[1]-wl_goal[0]

        hdul.writeto('{}spec_{}_ferre-syn.fits'.format(pad,name), overwrite=True)
        
    print('Created fits files in {}'.format(outfolferre))
    
    par[['name']].to_csv('{}starnames.txt'.format(outfolferre), index=None, header=None)
    
    return None


def getrvlist(outfolferre, rvfile):

    files = glob.glob('{}rvout_ferrespec/*3900-5500.txt'.format(outfolferre))
    outrv = rvfile

    list = pd.DataFrame()
    bestlist = pd.DataFrame()

    b=0
    for file in files:
        try:
            data = pd.read_csv(file,sep='\s+', comment="#", header=None)
        except:
            print('ERROR')
            b+=1
            continue

        name = data[0].iloc[0][:-15]
        print(name)
        data['name'] = name[5:]

        bestlist = pd.concat([bestlist, data], ignore_index=True)

    bestlist = bestlist.sort_values(by='name')

    bestlist = bestlist[bestlist[10] != "INDEF"]

    bestlist.to_csv(outrv,index=None, columns=['name',10,12], header=None)

    print('made {} file'.format(outrv))
    
    return None


def rvpyraf(outfolferre):

    from pyraf import iraf

    # Load necessary IRAF packages
    iraf.rv()
    
    os.system("mkdir {}rvout_ferrespec".format(outfolferre))

    # Read starnames from file
    with open("{}starnames.txt".format(outfolferre), "r") as f:
        for line in f:
            s1 = line.strip()
            print(s1)
            
            infile    = "{}spec_".format(outfolferre) + s1 + "_ferre-obs.fits"
            outfile   = "{}spec_".format(outfolferre) + s1 + "_ferre-syn.fits"
            rvoutfile = "{}rvout_ferrespec/".format(outfolferre) + s1 + "_3900-5500"

            iraf.fxcor(
                objects=infile,
                template=outfile,
                apertures="*",
                cursor="",
                continuum="both",
                filter="none",
                rebin="object",
                pixcorr="no",
                osample="3900-5500",
                rsample="3900-5500",
                apodize=0.1,
                function="gaussian",
                width="INDEF",
                height=0.0,
                peak="no",
                minwidth=3.0,
                maxwidth=21.0,
                weights=1.0,
                background="INDEF",
                window=200.0,
                wincenter="INDEF",
                output=rvoutfile,
                verbose="txtonly",
                imupdate="no",
                graphics="stdgraph",
                interactive="no",
                autowrite="yes",
                autodraw="yes",
                ccftype="text",
                continpars="",
                filtpars="",
                keywpars=""
            )

def plotspec(starnames, run, outfol, vs='', newgrid=False):

    font = {'size': 14}
    matplotlib.rc('font', **font)

    obs_orig = pd.read_csv('{}{}.frd'.format(outfol,run), sep='\s+', header=None)
    err = pd.read_csv('{}{}.err'.format(outfol,run), sep='\s+', header=None)
    wl = pd.read_csv('{}{}.wav'.format(outfol,run), sep='\s+', header=None)
    obs = pd.read_csv('{}{}{}.obs'.format(outfol,run,vs), sep='\s+', header=None)
    mod = pd.read_csv('{}{}{}.mdl'.format(outfol,run,vs), sep='\s+', header=None)
    if newgrid == False:
        par = pd.read_csv('{}{}{}.opf'.format(outfol,run,vs), sep='\s+', header=None, names=['name','feh_ferre','cfe_ferre','teff_ferre','logg_ferre','e_feh_ferre','e_cfe_ferre','e_teff_ferre','e_logg_ferre','zero_ferre','SNR_ferre','chi2_ferre'])
    else:
        par = pd.read_csv('{}{}{}.opf'.format(outfol,run,vs), sep='\s+', header=None, names=['name','teff_ferre','logg_ferre','feh_ferre','cfe_ferre','e_teff_ferre','e_logg_ferre','e_feh_ferre','e_cfe_ferre','zero_ferre','SNR_ferre','chi2_ferre'])


    sub = par[(par.name.isin(starnames))]


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
        ax2.plot(wl.iloc[i], err.iloc[i], color='deepskyblue', label='noise spectrum')
        ax2.set_ylim(0,np.median(obs_orig.iloc[i])+np.median(obs_orig.iloc[i])*0.5)
        ax2.set_xlim([minwl,maxwl])
        ax2.legend(fontsize=12)
        ax2.set_xlabel('Wavelength (angstrom)')

        ax4.plot(wl.iloc[i], obs_orig.iloc[i], color='black', label='non-normalised observation')
        ax4.plot(wl.iloc[i], err.iloc[i], color='deepskyblue', label='noise spectrum')
        ax4.set_ylim(0,np.median(obs_orig.iloc[i])+np.median(obs_orig.iloc[i])*0.5)
        ax4.set_xlim([3900,4000])
        ax4.get_yaxis().set_visible(False)

        gs1.tight_layout(fig)
        gs1.update(hspace=0, wspace=0.02)

    #    plt.savefig('figures/{}.jpg'.format(ster), dpi=150)

        plt.show()

