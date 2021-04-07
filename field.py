
import os
import errno
import numpy as np

import astropy.constants as const
from astropy.io import fits
from astropy.table import Table
import astropy.units as u
from astropy.wcs import WCS

from core import save_cutout

def determine_rms(segPath, files) :
    '''
    Estimate the RMS noise level in each flter image by using the
    segmentation map as a basic mask. Mask all the pixels that are
    associated with an object and determine the noise level from the
    remaining background pixels.
    
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
        
        mask = np.where(segMap == 0)
        masked_sci = sci[mask]
        masked_sci[masked_sci > 10] = np.nan
        rms = np.sqrt(np.nanmean(np.square(masked_sci)))
        rmses.append(rms)
    
    return rmses

def save_cutouts(cluster, final_objs_path, filters, rms, files, segPath) :
    '''
    Save the cutouts for all science objects in a given filter, as well as
    noise cutouts. Then move to the next filter, saving those cutouts. If
    there are N filters and M final science objects, this will result in
    2*N*M cutout images.
    
    Returns
    -------
    None.
    
    '''
    
    extent = 5 # maximum radius of image will be 5 times R_e
    
    sci_objs = Table.read(final_objs_path)
    sci_objs = sci_objs[np.where(sci_objs['id'] == 5830)] # for testing
    ID, r_e = sci_objs['id'], sci_objs['flux_radius']*u.pixel
    ra, dec, redshift = sci_objs['ra'], sci_objs['dec'], sci_objs['z']
    sma, smb = sci_objs['a_image'], sci_objs['b_image']
    pa = sci_objs['theta_J2000']
    
    dictionary = {}
    for i in range(len(filters)) :
        dictionary[filters[i]] = {'rms':rms[i], 'file':files[i]}
    
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
            scale = hdr['D001SCAL']*u.arcsec/u.pixel
            wcs = WCS(hdr)
            science = hdu[0].data
            
            for galaxy in sci_objs :
                rms = dictionary.get(filt).get('rms')
                noise = np.sqrt(science/(exposure.value) + rms**2)
                id_string = (str(ID.data)).strip('[]')
                
                # science file
                outfile = '{}/cutouts/{}_ID_{}_{}.fits'.format(
                    cluster, cluster, id_string, filt)
                save_cutout(ra, dec, extent*r_e*scale, science, wcs,
                            outfile, exposure.value, photfnu.value,
                            scale.value, rms, r_e.value[0], redshift[0],
                            sma[0], smb[0], pa[0])
                
                # noise file
                noise_outfile = '{}/cutouts/{}_ID_{}_{}_noise.fits'.format(
                    cluster, cluster, id_string, filt)
                save_cutout(ra, dec, extent*r_e*scale, noise, wcs,
                            noise_outfile, exposure.value, photfnu.value,
                            scale.value, rms, r_e.value[0], redshift[0],
                            sma[0], smb[0], pa[0])
                
                # segmentation map
                with fits.open(segPath) as hdu :
                    segMap = hdu[0].data
                segmap_outfile = '{}/cutouts/{}_ID_{}_segmap.fits'.format(
                    cluster, cluster, id_string)
                try :
                    save_cutout(ra, dec, extent*r_e*scale, segMap, wcs,
                                segmap_outfile, -999, -999, scale.value,
                                -999, r_e.value[0], redshift[0],
                                sma[0], smb[0], pa[0])
                except OSError as error :
                    if error.errno == errno.EEXIST :
                        pass
    
    return

def save_filterset(cluster, filters) :
    '''
    Save a file containing the filters for the cluster, one per line. This
    will be used when fitting the photometry in prospector.
    
    Returns
    -------
    None.
    
    '''
    
    filterset_file = '{}/filters.txt'.format(cluster)
    with open(filterset_file, 'a+') as file :
        for filt in filters :
            file.write('hff_{}\n'.format(filt))
        
        return
