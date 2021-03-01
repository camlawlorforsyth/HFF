
import os
import glob
import numpy as np

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

def vorbin_all(inDir, outDir):
    '''
    Open the F160W cutout image and compute the Voronoi bins. Save the
    image of the Voronoi bins for subsequent masking, as well as the locations
    of the bins's luminosity weighted centroids.
    
    Parameters
    ----------
    inDir : string
        The input directory to access files from.
    outDir : string
        The output directory to save files to.
    
    Returns
    -------
    None.
    
    '''
    
    cluster = inDir[:-9] # get the name of the cluster from the input directory
    SNR = 120 # the target signal-to-noise ratio to use
    
    f160w_paths = '{}{}_ID_*_f160w.fits'.format(inDir, cluster)
    f160w_file_list = glob.glob(f160w_paths) # get all F160W images
    f160w_files, IDs = [], [] # get a list of IDs corresponding to each galaxy
    for file in f160w_file_list :
        file = file.replace(os.sep, '/')
        f160w_files.append(file)
        IDs.append(file.split('_')[2])
    
    noise_paths = '{}{}_ID_*_f160w_noise.fits'.format(inDir, cluster)
    noise_file_list = glob.glob(noise_paths) # get all F160W noise images
    noise_files = []
    for file in noise_file_list :
        file = file.replace(os.sep, '/')
        noise_files.append(file)
    
    for i in range(len(f160w_files)) :
        outfile = '{}{}_ID_{}_vorbins.npz'.format(outDir, cluster, IDs[i])
        sci, dim, photnu = open_cutout(f160w_files[i])
        noise, _, _ = open_cutout(noise_files[i])
        
        xs, ys, binNum, xBar, yBar, SN, nPixels = vorbin_data(sci, noise,
                                                              dim, SNR)
        bins_image = binNum.reshape(dim) # reshape the bins into an image
        
        np.savez(outfile, image=bins_image, x=xs, y=ys, binNum=binNum,
                  xbar=xBar, ybar=yBar, SN=SN, nPixels=nPixels)
    
    return
