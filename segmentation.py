
import numpy as np

from astropy.io import fits

import initial

def add_gaussian_background(infile, outfile, rms) :
    
    with fits.open(infile) as hdu :
        hdr = hdu[0].header
        data = hdu[0].data
    
    background = np.random.normal(loc=0.0, scale=rms, size=data.shape)
    
    hdu = fits.PrimaryHDU(data + background)
    hdu.header = hdr
    hdu.writeto(outfile)
    
    return

def sextractor_prep(cluster) :
    
    # CLUSTERS
    
    if cluster == 'a370' :
        _, _, _, _, _, _, _, _, _, models, RMS = initial.a370_params()
    
    if cluster == 'a1063' :
        _, _, _, _, _, _, _, _, _, models, RMS = initial.a1063_params()
    
    if cluster == 'a2744' :
        _, _, _, _, _, _, _, _, _, models, RMS = initial.a2744_params()
    
    if cluster == 'm416' :
        _, _, _, _, _, _, _, _, _, models, RMS = initial.m416_params()
    
    if cluster == 'm717' :
        _, _, _, _, _, _, _, _, _, models, RMS = initial.m717_params()
    
    if cluster == 'm1149' :
        _, _, _, _, _, _, _, _, _, models, RMS = initial.m1149_params()
    
    # PARALLEL FIELDS
    
    if cluster == 'a370par' :
        _, _, _, _, _, _, _, _, _, models, RMS = initial.a370par_params()
    
    if cluster == 'a1063par' :
        _, _, _, _, _, _, _, _, _, models, RMS = initial.a1063par_params()
    
    if cluster == 'a2744par' :
        _, _, _, _, _, _, _, _, _, models, RMS = initial.a2744par_params()
    
    if cluster == 'm416par' :
        _, _, _, _, _, _, _, _, _, models, RMS = initial.m416par_params()
    
    if cluster == 'm717par' :
        _, _, _, _, _, _, _, _, _, models, RMS = initial.m717par_params()
    
    if cluster == 'm1149par' :
        _, _, _, _, _, _, _, _, _, models, RMS = initial.m1149par_params()
    
    outfile = 'sextractor/images/{}_f160w_model_sextractor.fits'.format(cluster)
    add_gaussian_background(models[-1], outfile, RMS[-1])
    
    return
