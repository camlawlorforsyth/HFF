
import os
import numpy as np
from distutils.util import strtobool

import astropy.constants as const
from astropy.io import fits
from astropy.table import join, Table
import astropy.units as u

from core import save_cutout

def combine_with_issues(cluster) :
    '''
    Join the final science objects with the issues file, saving the result.
    
    Parameters
    ----------
    cluster : string
        Operate on the files of this cluster.
    
    Returns
    -------
    None.
    
    '''
    
    sci_objs = Table.read('{}/{}_sample.fits'.format(cluster, cluster))
    issues = Table.read('{}/{}_issues.csv'.format(cluster, cluster))
    
    combined = join(sci_objs, issues, keys='id')
    
    # ensure that only galaxies with 'f160w_use' are included
    f160w_use = np.array([strtobool(string.lower()) for
                          string in combined['f160w_use']])
    combined = combined[f160w_use > 0]
    
    combined.write('{}/{}_sample-with-use-cols.fits'.format(cluster, cluster))
    
    return

def check_number_of_bins(cluster) :
    # Check for galaxies that don't have a photometry file, and therefore have
    # no bins. Also check for any galaxies that only have a single bin.
    
    table = Table.read('{}/{}_sample-with-use-cols.fits'.format(
        cluster, cluster))
    
    bins = []
    for ID in table['id'] :
        phot_file = '{}/photometry/{}_ID_{}_photometry.fits'.format(
            cluster, cluster, ID)
        if os.path.exists(phot_file) :
            phot = Table.read(phot_file)
            bins.append(len(phot))
        else :
            bins.append(0)
    
    table['bins'] = bins
    table = table[table['bins'] > 1]
    table.write('{}/{}_sample-with-use-cols-and-bins.fits'.format(
        cluster, cluster))
    
    return

def determine_rms(segPath, files) :
    '''
    Estimate the RMS noise level in each flter image by using the
    segmentation map as a basic mask. Mask all the pixels that are
    associated with an object and determine the noise level from the
    remaining background pixels that have significant values.
    
    Returns
    -------
    rmses : list
        List of rms values, one for each filter.
    
    '''
    
    with fits.open(segPath) as hdu :
        segMap = hdu[0].data
    
    rmses = []
    for file in files :
        with fits.open(file) as hdu :
            sci = hdu[0].data
            exptime = hdu[0].header['EXPTIME']
        
        threshold = 0.001/exptime
        
        bkg = sci[segMap == 0].flatten()
    
        neg = -1*bkg[bkg < 0]
        neg = neg[neg > threshold]
        neg = -1*neg
        
        pos = bkg[bkg >= 0]
        pos = pos[pos > threshold]
        
        background = np.concatenate((neg, pos), axis=None)
        
        # mask out large values that will influence the rms
        background[background > 10] = np.nan
        background[background < -10] = np.nan
        
        rms = np.sqrt(np.nanmean(np.square(background)))
        rmses.append(rms)
    
    return rmses

def save_cutouts(cluster, sample_path, filters, rms, files, segPath, bCGsegPath,
                 models, redshift_type='z_spec') :
    '''
    Save the cutouts for all science objects in a given filter, as well as
    noise cutouts and segmentation map cutouts. Then move to the next filter,
    saving those cutouts. If there are N filters and M final science objects,
    this will result in 2*(N+1)*M cutout images.
    
    Returns
    -------
    None.
    
    '''
    
    extent = 5 # maximum radius of image will be 5 times R_e
    
    sci_objs = Table.read(sample_path)
    xx, yy = sci_objs['x'], sci_objs['y']
    redshift = sci_objs[redshift_type]
    ID, r_e = sci_objs['id'], sci_objs['flux_radius']*u.pixel
    sma, smb = sci_objs['a_image'], sci_objs['b_image']
    pa = sci_objs['theta_J2000']
    
    dictionary = {}
    for i in range(len(filters)) :
        dictionary[filters[i]] = {'rms':rms[i], 'file':files[i],
                                  'model':models[i]}
    
    os.makedirs('{}/cutouts'.format(cluster), exist_ok=True) # ensure the
        # output directory for the cutouts is available
    
    for filt in dictionary :
        with fits.open(dictionary.get(filt).get('file')) as hdu :
            hdr = hdu[0].header
            exposure = hdr['EXPTIME']*u.s
            # photfnu = hdr['PHOTFNU'] # the ACS images don't have PHOTFNU
            photflam = hdr['PHOTFLAM']*u.erg/(u.cm**2)/u.AA/u.electron
            photplam = hdr['PHOTPLAM']*u.AA
            photfnu = ((photplam**2)*
                       photflam/const.c).to(u.Jy*u.s/u.electron)
            try :
                scale = hdr['D001SCAL']*u.arcsec/u.pixel
            except KeyError :
                scale = 0.06*u.arcsec/u.pixel
            science = hdu[0].data
        
        with fits.open(dictionary.get(filt).get('model')) as hdum :
            bcg_model = hdum[0].data
        
        for i in range(len(sci_objs)) :
            # add the bCG model into the bCG-subtracted images for the bCGs
            if ID[i] < 20000 :
                science = science
            else :
                science = science + bcg_model
            
            # calculate the noise for the filter
            rms = dictionary.get(filt).get('rms')
            noise = np.sqrt(science/(exposure.value) + rms**2)
            
            # science file
            outfile = '{}/cutouts/{}_ID_{}_{}.fits'.format(
                cluster, cluster, str(ID[i]), filt)
            save_cutout(xx[i], yy[i], extent*r_e[i],
                        science, outfile, exposure.value,
                        photfnu.value, scale.value, rms, r_e.value[i],
                        redshift[i], sma[i], smb[i], pa[i])
            
            # noise file
            noise_outfile = '{}/cutouts/{}_ID_{}_{}_noise.fits'.format(
                cluster, cluster, str(ID[i]), filt)
            save_cutout(xx[i], yy[i], extent*r_e[i],
                        noise, noise_outfile, exposure.value,
                        photfnu.value, scale.value, rms, r_e.value[i],
                        redshift[i], sma[i], smb[i], pa[i])
    
    # segmentation map
    for i in range(len(sci_objs)) :
        with fits.open(segPath) as hdu :
            segMap = hdu[0].data
        segmap_outfile = '{}/cutouts/{}_ID_{}_segmap.fits'.format(
            cluster, cluster, str(ID[i]))
        save_cutout(xx[i], yy[i], extent*r_e[i],
                    segMap, segmap_outfile, -999,
                    -999, scale.value, -999, r_e.value[i],
                    redshift[i], sma[i], smb[i], pa[i])
    
    # bCG segmentation map
    for i in range(len(sci_objs)) :
        with fits.open(bCGsegPath) as hdu :
            bCGsegMap = hdu[0].data
        bCG_segmap_outfile = '{}/cutouts/{}_ID_{}_segmap-bCG.fits'.format(
            cluster, cluster, str(ID[i]))
        save_cutout(xx[i], yy[i], extent*r_e[i],
                    bCGsegMap, bCG_segmap_outfile, -999,
                    -999, scale.value, -999, r_e.value[i],
                    redshift[i], sma[i], smb[i], pa[i])
    
    return

def save_filterset(cluster, filters) :
    '''
    Save a file containing the filters for the cluster, one per line. This
    will be used when fitting the photometry in prospector.
    
    Returns
    -------
    None.
    
    '''
    
    filterset_file = '{}/{}_filters.txt'.format(cluster, cluster)
    with open(filterset_file, 'a+') as file :
        for filt in filters :
            file.write('hff_{}\n'.format(filt))
    
    return
