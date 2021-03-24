
import os
import glob
import numpy as np

from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
import astropy.units as u

from core import open_cutout

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def determine_fluxes(binDir, cutoutDir, outDir, filters) :
    '''
    Determine the flux in every Voronoi bin for a given object for a given
    filter. Then move to the next subsequent filter and determine the fluxes
    in the corresponding Voronoi bins for that filter. Repeat for all filters.
    Create a table which includes all determined fluxes. Save table to file
    for subsequent use with Prospector.
    
    Parameters
    ----------
    vorbinDir : TYPE
        DESCRIPTION.
    cutoutDir : TYPE
        DESCRIPTION.
    outDir : TYPE
        DESCRIPTION.
    filters : TYPE
        DESCRIPTION.
    
    Returns
    -------
    None.
    
    '''
    
    cluster = binDir[:-6] # get the name of the cluster from the input directory
    
    bin_paths = '{}{}_ID_*_annuli.npz'.format(binDir, cluster)
    bin_file_list = glob.glob(bin_paths) # get all bin numpy files
    bin_files, IDs = [], [] # get a list of IDs corresponding to each galaxy
    for file in bin_file_list :
        file = file.replace(os.sep, '/')
        bin_files.append(file)
        IDs.append(file.split('_')[2])
    
    for i in range(len(bin_files)) :
        outfile = '{}{}_ID_{}_photometry.fits'.format(outDir, cluster, IDs[i])
        
        bin_data = np.load(bin_files[i])
        bins_image = bin_data['image']
        sma, smb = bin_data['sma'], bin_data['smb']
        flux, err = bin_data['flux'], bin_data['err']
        nPixels = bin_data['nPixels']
        
        numBins = np.nanmax(bins_image) + 1 # accounts for python 0-index
        
        photometry = Table()
        photometry['bin'] = range(int(numBins))
        photometry['sma'], photometry['smb'] = sma, smb
        photometry['flux'], photometry['err'] = flux, err
        photometry['SN'], photometry['nPixels'] = flux/err, nPixels
        
        for filt in filters :
            sci_file = '{}{}_ID_{}_{}.fits'.format(cutoutDir, cluster, IDs[i],
                                                   filt)
            noise_file = '{}{}_ID_{}_{}_noise.fits'.format(cutoutDir, cluster,
                                                           IDs[i], filt)
            segmap_file = '{}{}_ID_{}_segmap.fits'.format(cutoutDir, cluster,
                                                          IDs[i])
            
            (sci, dim, photfnu, r_e,
             redshift, sma, smb, pa) = open_cutout(sci_file)
            noise, _, _, _, _, _, _, _ = open_cutout(noise_file)
            segMap, _, _, _, _, _, _, _ = open_cutout(segmap_file)
            
            # make a copy of the science image based on the segmentation map
            ID = int(IDs[i])
            new_sci = sci.copy()
            new_sci[(segMap > 0) & (segMap != ID)] = 0 # don't mask out the sky
            
            # and for the noise image as well
            new_noise = noise.copy()
            new_noise[(segMap > 0) & (segMap != ID)] = 0
            
            lumDist = cosmo.luminosity_distance(redshift)
            
            fluxes, uncerts, R_e, redshifts, lumDists = [], [], [], [], []
            for val in range(int(numBins)) :
                mask = np.where(bins_image == val)
                
                masked_sci = new_sci[mask]
                flux = photfnu*np.nansum(masked_sci)
                fluxes.append(flux)
                
                masked_noise = new_noise[mask]
                uncert = photfnu*np.nansum(masked_noise)
                uncerts.append(uncert)
                
                R_e.append(r_e)
                redshifts.append(redshift)
                lumDists.append(lumDist.value)
            
            photometry[filt + '_flux'] = fluxes*u.Jy
            photometry[filt + '_err'] = uncerts*u.Jy
        
        photometry['R_e'] = R_e*u.pix
        photometry['z'] = redshifts
        photometry['lumDist'] = lumDists*u.Mpc
        
        photometry.write(outfile)
    
    return
