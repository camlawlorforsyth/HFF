
import os
import glob
import numpy as np

from astropy.table import Table
import astropy.units as u

from core import open_cutout

def determine_fluxes(vorbinDir, cutoutDir, outDir, filters) :
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
    cluster = vorbinDir[:-9] # get the name of the cluster from the input directory
    R_e_percent = 0.1
    
    vorbin_paths = '{}{}_ID_*_vorbins.npz'.format(vorbinDir, cluster)
    vorbin_file_list = glob.glob(vorbin_paths) # get all vorbin numpy files
    vorbin_files, IDs = [], [] # get a list of IDs corresponding to each galaxy
    for file in vorbin_file_list :
        file = file.replace(os.sep, '/')
        vorbin_files.append(file)
        IDs.append(file.split('_')[2])
    
    for i in range(len(vorbin_files)) :
        outfile = '{}{}_ID_{}_photometry.fits'.format(outDir, cluster, IDs[i])
        vorbin_data = np.load(vorbin_files[i])
        bins_image = vorbin_data['image']
        xBar, yBar = vorbin_data['xbar'], vorbin_data['ybar']
        SN, nPixels = vorbin_data['SN'], vorbin_data['nPixels']
        numBins = np.max(bins_image) + 1 # accounts for python 0-index
        
        photometry = Table()
        photometry['vorbin'] = range(numBins)
        photometry['xBar'], photometry['yBar'] = xBar, yBar
        photometry['SN'], photometry['nPixels'] = SN, nPixels
        
        for filt in filters :
            sci_file = '{}{}_ID_{}_{}.fits'.format(cutoutDir, cluster, IDs[i],
                                                   filt)
            noise_file = '{}{}_ID_{}_{}_noise.fits'.format(cutoutDir, cluster,
                                                           IDs[i], filt)
            sci, dim, photfnu, r_e = open_cutout(sci_file)
            noise, _, _, _ = open_cutout(noise_file)
            fluxes, uncerts, R_e = [], [], []
            for val in range(numBins) :
                mask = np.where(bins_image == val)
                masked_sci = sci[mask]
                flux = photfnu*np.nansum(masked_sci)
                fluxes.append(flux)
                masked_noise = noise[mask]
                uncert = photfnu*np.nansum(masked_noise)
                uncerts.append(uncert)
                R_e.append(r_e)
            photometry[filt + '_flux'] = fluxes*u.Jy
            photometry[filt + '_err'] = uncerts*u.Jy
        photometry['R_e'] = R_e
        photometry['use'] = np.sqrt(nPixels/np.pi) < R_e_percent*np.array(R_e)
        photometry.write(outfile)
    
    return
