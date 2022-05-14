
import os
import glob
import numpy as np

from astropy.io import fits

import core
import plotting

def sample_randomly(cluster, filt, dims) :
    
    if cluster[0] == 'a' :
        name = 'abell' + cluster[1:]
    if cluster[0] == 'm' :
        name = 'macs0' + cluster[1:]
    
    if name == 'macs01149' :
        name = 'macs1149'
    
    path = ('misc/{}clu_misc/images/psf_matched/'.format(name) +
            '{}_bcgs_out_{}_psf_bkg_drz.fits.gz'.format(name, filt))
    
    segPath = ('misc/{}clu_misc/photometry/'.format(name) +
               '{}_bcgs_out_detection_seg.fits.gz'.format(name))
    
    with fits.open(segPath) as hdu :
        segMap = hdu[0].data
    
    with fits.open(path) as hdu :
        sci = hdu[0].data
        # full_dims = data.shape
        exptime = hdu[0].header['EXPTIME']
    
    threshold = 0.001/exptime
    
    bkg = sci[segMap == 0].flatten()
    
    neg = -1*bkg[bkg < 0]
    neg = neg[neg > threshold]
    neg = -1*neg
    
    pos = bkg[bkg >= 0]
    pos = pos[pos > threshold]
    
    background = np.concatenate((neg, pos), axis=None)
    rms = np.sqrt(np.nanmean(np.square(background)))
    rms_in_electrons = exptime*rms
    print(rms_in_electrons)
    
    '''
    random_signals, random_rmses = [], []
    for i in range(len(dims)) :
        random_x, random_y = np.int_(full_dims*np.random.rand(2))
        
        dim = dims[i]
        
        # try :
        random_x_end = random_x + dim
        random_y_end = random_y + dim
        
        excess_x = full_dims[0] - random_x_end
        excess_y = full_dims[1] - random_y_end
        
        if excess_x < 0 :
            random_x = random_x + excess_x - 1
        if excess_y < 0 :
            random_y = random_y + excess_y - 1
        
        random_cutout = data[random_x:random_x_end, random_y:random_y_end]
        median = np.median(random_cutout)
        rms = np.sqrt(np.nanmean(np.square(random_cutout)))
        
        random_signals.append(median)
        random_rmses.append(rms)
    
    random_signals = np.array(random_signals)
    random_rmses = np.array(random_rmses)
    '''
    
    return

def compare_rms(cluster, filt) :
    
    cutout_files = glob.glob('{}/cutouts/{}_ID_*_{}.fits'.format(cluster,
                                                                 cluster,
                                                                 filt))
    
    cutouts = []
    for file in cutout_files :
        file = file.replace(os.sep, '/') # compatibility for Windows
        cutouts.append(file)
    
    median_signals, rmses, dims = [], [], []
    for cutout in cutouts :
        data, exptime = core.open_cutout(cutout, exp=True)
        dim = data.shape[0]
        dims.append(dim)
        
        electrons_image = exptime*data
        
        median = np.median(electrons_image)
        rms = np.sqrt(np.nanmean(np.square(electrons_image)))
        
        median_signals.append(median)
        rmses.append(rms)
    
    median_signals = np.array(median_signals)
    rmses = np.array(rmses)
    dims = np.array(dims)
    
    neg_mask = (median_signals < 0)
    neg_signals = np.abs(median_signals[neg_mask])
    neg_rmses = rmses[neg_mask]
    
    random_signals, random_rmses = sample_randomly(cluster, filt, dims)
    neg_mask = (random_signals < 0)
    neg_random_signals = np.abs(random_signals[neg_mask])
    neg_random_rmses = random_rmses[neg_mask]
    
    xs = [median_signals, neg_signals, random_signals, neg_random_signals]
    ys = [rmses, neg_rmses, random_rmses, neg_random_rmses]
    labels = ['signal', 'abs(signal)', 'random signal', 'abs(random signal)']
    colors = ['k', 'r', 'b', 'g']
    
    # plotting.plot_simple_multi(xs, ys, labels, colors,
    #                             xlabel='median signal', ylabel='RMS')
    
    return

# sample_randomly('a2744', 'f275w', [0])
# compare_rms('a2744', 'f275w')

with fits.open('a2744/cutouts/a2744_ID_4369_f275w.fits') as hdu :
    sci = hdu[0].data
    exptime = hdu[0].header['EXPTIME']
    photfnu = hdu[0].header['PHOTFNU']

with fits.open('a2744/cutouts/a2744_ID_4369_f275w_noise.fits') as hdu :
    noise = hdu[0].data

bin_data = np.load('a2744/bins/a2744_ID_4369_annuli.npz')
bin_image = bin_data['image']

central_bin = np.int64(bin_image == 0)*sci
total = np.sum(central_bin) # electron/s
janskys = photfnu*total # janskys
electrons = total*exptime
print(janskys)

central_noise = np.int64(bin_image == 0)*noise
total_noise = np.sqrt(np.sum(np.square(central_noise))) # electron/s
janskys_noise = photfnu*total_noise # janskys
electrons_noise = total_noise*exptime
print(janskys_noise)

plotting.display_image_simple(central_bin*exptime, norm=None)

