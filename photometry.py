
import numpy as np

from astropy.table import Table
import astropy.units as u

from core import open_cutout

def determine_flux(inDir) :
    '''
    Determine the flux in every Voronoi bin for a given object for a given
    filter. Then move to the next subsequent filter and determine the fluxes
    in the corresponding Voronoi bins for that filter. Repeat for all filters.
    Create a table which includes all determined fluxes. Save table to file
    for subsequent use with Prospector.
    
    Returns
    -------
    None.
    
    '''
    
    import numpy.ma as ma
    
    cluster = 'abell2744'
    inDir = 'tests_using_A2744/cutouts/'
    infile = 'abell2744_ID_5830_vorbins.npz'
    
    filters = ['f275w', 'f336w', 'f435w', 'f606w', 'f814w', 'f105w', 'f125w',
               'f140w', 'f160w']
    
    vorbin_data = np.load(inDir + infile)
    bins_image = vorbin_data['image']
    xbar = vorbin_data['xbar']
    ybar = vorbin_data['ybar']
    
    numBins = np.max(bins_image) + 1
    
    base_string = infile[:-11]
    table = Table()
    table['vorbin'] = range(numBins)
    for filt in filters :
        data, dim, photfnu = open_cutout(inDir + base_string + filt + '.fits')
        fluxes = []
        for val in range(numBins) :
            mask = np.where(bins_image == val, 0, 1)
            masked_data = ma.masked_array(data, mask)
            flux = photfnu*np.sum(masked_data)
            fluxes.append(flux)
        table[filt + '_flux'] = fluxes*u.Jy
    table.write(inDir + base_string + 'photometry.fits')
    
    return
