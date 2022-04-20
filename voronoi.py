
import os
import numpy as np

from astropy.table import Table
from vorbin.voronoi_2d_binning import voronoi_2d_binning

import plotting as plt
from core import open_cutout

def vorbin_data(signal, noise, dim, targetSN, display=False, printraw=False) :
    '''
    Use the vorbin package to bin the data using Voronoi tesselations.
    
    Parameters
    ----------
    signal : numpy.ndarray
        The data to adaptively bin.
    noise : numpy.ndarray
        The noise corresponding to the signal.
    dim : tuple
        The dimension of the input arrays.
    targetSN : float
        The target signal-to-noise ratio to use for binning.
    display : bool, optional
        Boolean to plot the resulting bin image. The default is False.
    printraw : bool, optional
        Boolean to print the resulting bin image. The default is False.
    
    Returns
    -------
    xs : numpy.ndarray
        The x-coordinate of the pixels to bin.
    ys : numpy.ndarray
        The y-coordinate of the pixels to bin.
    binNum : numpy.ndarray
        The bin number of each pixel in the image.
    xbar : numpy.ndarray
        The x-coordinates of the bins luminosity weighted centroids.
    ybar : numpy.ndarray
        The y-coordinates of the bins luminosity weighted centroids.
    sn : numpy.ndarray
        The final signal-to-noise ratio of each bin.
    npixels : numpy.ndarray
        The number of pixels contained in each bin.
    
    '''
    
    signal = signal.flatten() # vorbin requires 1D arrays
    noise = noise.flatten()
    
    ys, xs = np.mgrid[:dim[0], :dim[1]] # position arrays for each pixel,
    xs = xs.flatten() - dim[1]/2        # relative to the center of the image
    ys = ys.flatten() - dim[0]/2
    
    (binNum, xNode, yNode, xbar, ybar, sn, npixels,
     scale) = voronoi_2d_binning(xs, ys, signal, noise, targetSN, pixelsize=1,
                                 plot=False, quiet=True)
    
    if display :
        bins = binNum.reshape(dim)
        plt.display_image_simple(bins)
    
    if printraw :
        bins = binNum.reshape(dim)
        print(np.flipud(bins))
    
    return xs, ys, binNum, xbar, ybar, sn, npixels

def vorbin_all(cluster):
    '''
    Open the F160W cutout image and compute the Voronoi bins. Save the
    image of the Voronoi bins for subsequent masking, as well as the locations
    of the bins's luminosity weighted centroids.
    
    Parameters
    ----------
    cluster : string
        Operate on the files of this cluster.
    
    Returns
    -------
    None.
    
    '''
    
    os.makedirs('{}/vorbins'.format(cluster), exist_ok=True) # ensure the output
        # directory for the vorbins is available
    
    table = Table.read('{}/{}_sample-with-use-cols.fits'.format(cluster, cluster))
    IDs = table['ID']    
    
    for ID in IDs :
        (sci, dim, photnu, r_e,
         redshift, sma, smb, pa) = open_cutout(
             '{}/cutouts/{}_ID_{}_f160w.fits'.format(cluster, cluster, ID))
        noise = open_cutout(
            '{}/cutouts/{}_ID_{}_f160w_noise.fits'.format(cluster, cluster, ID),
            simple=True)
        
        # target a signal-to-noise ratio of 400
        xs, ys, binNum, xBar, yBar, SN, nPixels = vorbin_data(sci, noise,
                                                              dim, 400)
        bins_image = binNum.reshape(dim) # reshape the bins into an image
        
        outfile = '{}/vorbins/{}_ID_{}_vorbins.npz'.format(cluster, cluster, ID)
        np.savez(outfile, image=bins_image, x=xs, y=ys, binNum=binNum,
                 xbar=xBar, ybar=yBar, SN=SN, nPixels=nPixels)
    
    return
