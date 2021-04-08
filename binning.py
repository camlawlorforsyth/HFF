
import os
import glob
import numpy as np

from core import open_cutout

def bin_all(cluster):
    '''
    Bin all the F160W cutout images using elliptical annuli. Save the image
    of the annuli bins for subsequent masking, as well as other useful
    quantities.
    
    Parameters
    ----------
    cluster : string
        Operate on the files of this cluster.
    
    Returns
    -------
    None.
    
    '''
    
    inDir = '{}/cutouts'.format(cluster)
    outDir = '{}/bins'.format(cluster)
    
    f160w_paths = '{}/{}_ID_*_f160w.fits'.format(inDir, cluster)
    f160w_file_list = glob.glob(f160w_paths) # get all F160W images
    f160w_files, IDs = [], [] # get a list of IDs corresponding to each galaxy
    for file in f160w_file_list :
        file = file.replace(os.sep, '/')
        f160w_files.append(file)
        IDs.append(file.split('_')[2])
    
    noise_paths = '{}/{}_ID_*_f160w_noise.fits'.format(inDir, cluster)
    noise_file_list = glob.glob(noise_paths) # get all F160W noise images
    noise_files = []
    for file in noise_file_list :
        file = file.replace(os.sep, '/')
        noise_files.append(file)
    
    segmap_paths = '{}/{}_ID_*_segmap.fits'.format(inDir, cluster)
    segmap_file_list = glob.glob(segmap_paths) # get all F160W noise images
    segmap_files = []
    for file in segmap_file_list :
        file = file.replace(os.sep, '/')
        segmap_files.append(file)
    
    os.makedirs('{}'.format(outDir), exist_ok=True) # ensure the output
        # directory for the bins is available
    
    for i in range(len(f160w_files)) :
        outfile = '{}/{}_ID_{}_annuli.npz'.format(outDir, cluster, IDs[i])
        
        (annuli_map, smas, smbs, fluxes, errs,
         nPixels, widths, pas) = annuli_bins(f160w_files[i], noise_files[i],
                                             segmap_files[i], int(IDs[i]))
        
        np.savez(outfile, image=annuli_map, sma=smas, smb=smbs,
                 flux=fluxes, err=errs, nPixels=nPixels, width=widths, pa=pas)
    
    return

def annuli_bins(f160w_file, noise_file, segmap_file, ID) :
    '''
    Determine the annuli to use for subsequent analysis.
    
    Parameters
    ----------
    f160w_file : string
        The F160W science file to bin.
    noise_file : string
        The noise file used in determining the bins.
    segmap_file : string
        Segmentation map file used to mask the above files.
    ID : int
        ID of the galaxy that is being binned.
    
    Returns
    -------
    annuli_map : numpy.ndarray
        DESCRIPTION.
    smas : numpy.ndarray
        The semi-major axis of the inner edge of each annulus.
    smbs : numpy.ndarray
        The semi-minor axis of the inner edge of each annulus.
    fluxes : list
        Flux contained in each annulus.
    errs : list
        Uncertainty on the flux contained in each annulus.
    nPixels_list : list
        Number of pixels contained in each annulus.
    widths : numpy.ndarray
        Widths of each annulus.
    pas : numpy.ndarray
        Array of position angles for the annuli.
    
    '''
    
    IR_pixel_scale = 0.128 # arcsec/pixel
    image_pixel_scale = 0.06 # arcsec/pixel
    targetSN = 30 # the target signal-to-noise ratio to use
    rin = 0 # the starting value for the semi-major axis of the ellipse
    
    (sci, dim, photnu, r_e,
     redshift, sma, smb, pa_deg) = open_cutout(f160w_file)
    noise, _, _, _, _, _, _, _ = open_cutout(noise_file)
    segMap, _, _, _, _, _, _, _ = open_cutout(segmap_file)
    
    eta = 1 - smb/sma # the ellipticity of the ellipse
    
    # determine the center of the ellipse, based on the size of the cutout
    xx, yy = np.indices(dim)
    x0, y0 = np.median(xx), np.median(yy)
    
    # make a copy of the science image based on the segmentation map
    new_sci = sci.copy()
    new_sci[(segMap > 0) & (segMap != ID)] = 0 # don't mask out the sky
    
    # and for the noise image as well
    new_noise = noise.copy()
    new_noise[(segMap > 0) & (segMap != ID)] = 0
    
    dr = max(IR_pixel_scale/image_pixel_scale, 0.1*r_e)
    pa = np.pi*pa_deg/180
    
    rins, fluxes, errs, widths, annuli, nPixels_list = [], [], [], [], [], []
    while rin < 5*r_e :
        (flux, err, rnew, width,
         annulus, nPixels) = compute_annuli(new_sci, new_noise, dim, (x0, y0),
                                            rin, dr, eta, pa, targetSN)
        if rnew < 5*r_e :
            fluxes.append(flux)
            errs.append(err)
            rins.append(rnew)
            widths.append(width)
            annuli.append(annulus)
            nPixels_list.append(nPixels)
        rin = rnew
    
    # create the annuli map for subsequent determination of photometric fluxes
    annuli_map = np.zeros(dim)
    for i in range(len(annuli)) :
        annuli_map += (i+1)*annuli[i]
    annuli_map[annuli_map == 0] = np.nan
    annuli_map -= 1
    
    smas = np.array(rins) # the semi-major axes of the inner annuli
    smbs = (1-eta)*smas # the semi-minor axes of the inner annuli
    widths = np.array(widths)
    pas = np.array([pa_deg]*len(smas))
    
    return annuli_map, smas, smbs, fluxes, errs, nPixels_list, widths, pas

def compute_annuli(sci, noise, dim, xy, rin, dr, eta, pa, targetSN) :
    '''
    
    
    Parameters
    ----------
    sci : numpy.ndarray
        Science data to bin.
    noise : numpy.ndarray
        Corresponding noise data used in determining bins.
    dim : tuple
        Shape of the science/noise data.
    xy : tuple
        Coordinates of the center of the annulus.
    rin : float
        Inner semi-major axis of the annulus.
    dr : float
        Width of the annulus.
    eta : float
        Ellipticity of the annulus.
    pa : float
        Position angle of the annulus in radians, measured anti-clockwise.
    targetSN : float
        Target signal-to-noise ratio to use for binning.
    
    Returns
    -------
    flux : float
        Flux in the annulus.
    err : float
        Total uncertainty, added in quadrature, within the annulus.
    rnew : float
        The new inner semi-major axis to use for the next annulus.
    dr : float
        The width of the annulus.
    annulus : numpy.ndarray
        Array which corresponds to a given annulus.
    nPixels : int
        Number of pixels in the annnulus.
    
    '''
    
    annulus, nPixels = elliptical_annulus(dim, xy, rin, dr, eta, pa)
    flux, err = compute_SNR(sci, noise, annulus)
    if flux/err < targetSN :
        try :
            return compute_annuli(sci, noise, dim, xy, rin, dr+0.01, eta, pa,
                                  targetSN)
        except RecursionError :
            return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    else :
        return flux, err, rin+dr, dr, annulus, nPixels

def compute_SNR(sci, noise, annulus) :
    '''
    
    
    Parameters
    ----------
    sci : numpy.ndarray
        Science data to bin.
    noise : numpy.ndarray
        Corresponding noise data used in determining bins.
    annulus : numpy.ndarray
        Array which corresponds to a given annulus.
    
    Returns
    -------
    flux : float
        Flux in the annulus.
    err : float
        Total uncertainty, added in quadrature, within the annulus.
    
    '''
    
    sci, noise = sci.copy(), noise.copy()
    sci[~annulus], noise[~annulus] = 0, 0
    flux, err = np.sum(sci), np.sqrt(np.sum(np.square(noise)))
    
    return flux, err

def elliptical_annulus(dim, xy, rin, dr, eta, pa) :
    '''
    
    
    Parameters
    ----------
    dim : tuple
        Shape of the science/noise data.
    xy : tuple
        Coordinates of the center of the annulus.
    rin : float
        Inner semi-major axis of the annulus.
    dr : float
        Width of the annulus.
    eta : float
        Ellipticity of the annulus.
    pa : float
        Position angle of the annulus in radians, measured anti-clockwise.
    
    Returns
    -------
    annulus : numpy.ndarray
        Array which corresponds to a given annulus.
    nPixels : int
        Number of pixels in the annnulus.
    
    '''
    
    inner = elliptical_mask(dim, xy, rin, eta, pa)
    outer = elliptical_mask(dim, xy, rin+dr, eta, pa)
    annulus = np.bitwise_xor(inner, outer) # this is the region of overlap
    
    nPixels = np.sum(annulus)
    
    return annulus, nPixels

def elliptical_mask(dim, xy, rin, eta, pa) :
    '''
    
    
    Parameters
    ----------
    dim : tuple
        Shape of the science/noise data.
    xy : tuple
        Coordinates of the center of the annulus.
    rin : float
        Inner semi-major axis of the annulus.
    eta : float
        Ellipticity of the annulus.
    pa : float
        Position angle of the annulus in radians, measured anti-clockwise.
    
    Returns
    -------
    mask : numpy.ndarray
        Array which corresponds to a given annulus.
    
    '''
    
    YY, XX = np.ogrid[:dim[0], :dim[1]]
    
    left_num = np.square((XX-xy[0])*np.cos(pa) + (YY-xy[1])*np.sin(pa)) 
    left_denom = np.square(rin)
    
    right_num = np.square(-(XX-xy[0])*np.sin(pa) + (YY-xy[1])*np.cos(pa))
    right_denom = np.square((1-eta)*rin)
    
    ellipse = left_num/left_denom + right_num/right_denom
    mask = ellipse <= 1
    
    return mask
