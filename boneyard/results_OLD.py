
import os
import numpy as np

from astropy.table import Table
from matplotlib import cm
from matplotlib.colors import LogNorm
import prospect.io.read_results as reader

from core import open_cutout
import plotting as plt

# TESTS #

(sci, dim, photfnu, R_e,
 redshift, sma, smb, pa) = open_cutout('a2744/cutouts/a2744_ID_5830_f160w.fits')
noise, _, _, _, _, _, _, _ = open_cutout('a2744/cutouts/a2744_ID_5830_f160w_noise.fits')

# plt.display_image_simple(sci, bad='black')
# plt.display_image_simple(noise, bad='black', vmax=0.001)
# plt.display_image_simple(sci/noise, bad='black', vmax=450)

vorbins = np.load('a2744/vorbins/a2744_ID_5830_vorbins.npz')
bins_image = vorbins['image']
plt.display_image_simple(bins_image[25:80, 25:80], bad='black', cmap=cm.prism,
                          norm=None)

# PLOT RESULTS FROM FITTING #
cluster = 'a2744'
inDir = cluster + '/results/'
vorbinDir = cluster + '/vorbins/'
galaxy = 5830

results_file = '{}{}_ID_{}_results.fits'.format(inDir, cluster, str(galaxy))
results = Table.read(results_file)

vorbin_file = '{}{}_ID_{}_vorbins.npz'.format(vorbinDir, cluster, str(galaxy))
bins_image = np.load(vorbin_file)['image']
numBins = np.max(bins_image) + 1 # accounts for python 0-index

mass_map = np.ones(bins_image.shape)
for val in range(numBins) :
    nPixels = results['nPixels'][val]
    mask = np.where(bins_image == val)
    mass_map[mask] = np.log10(np.power(10, results['logmass'][val])/nPixels)

plt.display_image_simple(mass_map[25:80, 25:80], bad='white', norm=None,
                         cmap=cm.viridis, vmin=7.2, vmax=8.8)

# INTERACT WITH RESULTS
cluster = 'a2744'
inDir = cluster + '/h5/'
outDir = cluster + '/results/'
galaxies = os.listdir(inDir) # get a list of all subdirectories corresponding
                             # to galaxies for the given cluster

for galaxy in galaxies :
    
    photDir = cluster + '/photometry/'
    photometry = '{}{}_ID_{}_photometry.fits'.format(photDir, cluster,
                                                     str(galaxy))
    
    
    
    table = Table.read(photometry)
    vorbins = table['vorbin']
    
    logmasses = []
    for vorbin in vorbins :
        if table['use'][vorbin] :
            h5_file = '{}{}/{}_ID_{}_BIN_{}_dynesty.h5'.format(inDir,
                                                               str(galaxy),
                                                               cluster,
                                                               str(galaxy),
                                                               str(vorbin))
            result, obs, model = reader.results_from(h5_file, dangerous=False)
            imax = np.argmax(result['lnprobability'])
            logmass = result['chain'][imax, :][2]
            logmasses.append(logmass)
        else :
            logmasses.append(np.nan)
    
    table['logmass'] = np.array(logmasses)
    
    final = '{}{}_ID_{}_results.fits'.format(outDir, cluster, str(galaxy))
    table.write(final)
    
    '''
    import glob
    
    h5_files = '{}{}/{}_ID_{}_BIN_*_dynesty.h5'.format(inDir, str(galaxy),
                                                       cluster, str(galaxy))
    # get a list of HDF5 files containing fit results for all vorbins for
    # a given galaxy, as denoted by vorbin number
    results = glob.glob(h5_files)
    # loop over all the results files in the directory
    for file in results :
        result, obs, model = reader.results_from(file, dangerous=False)
        # print(result.keys())
        # print()
        imax = np.argmax(result['lnprobability'])
        # theta_max = result['chain'][imax, :]
        # cornerfig = reader.subcorner(result) #, start=0, thin=1), #truths=theta_best,
                                      # fig=subplots(5,5,figsize=(27,27))[0])
        logmass = result['chain'][imax, :][2]
    '''

# USE DICTIONARY TO SAVE VALUES
# then use keys to get values when creating image        

'''
results = '{}/h5/{}_ID_*_.fits'.format(cluster, cluster)

# get a list of HDF5 files containing the results for all vorbins for a given
# galaxy, as denoted by ID and vorbin number
# a given galaxy, as denoted by ID
tables = glob.glob(photometries)

# loop over all the fits files in the directory
for file in tables :
    ID = file.split('_')[2] # the galaxy ID to fit the vorbins for
    table = Table.read(file)
    vorbins = table['vorbin'] # get a list of vorbin values
    
'''

def isophote_data() :
    
    from photutils.isophote import EllipseGeometry
    from photutils import EllipticalAperture
    from photutils.isophote import Ellipse
    
    test = 'a2744/cutouts/a2744_ID_5830_f160w.fits'
    (data, dim, photfnu, R_e,
     redshift, sma, smb, pa) = open_cutout(test)
    (noise, _, _, _,
     _, _, _, _) = open_cutout('a2744/cutouts/a2744_ID_5830_f160w_noise.fits')
    plt.display_image_simple(data, cmap=cm.viridis)
    
    xlen, ylen = dim[1], dim[0]
    
    geometry = EllipseGeometry(x0=xlen/2, y0=ylen/2, sma=20, eps=0.5,
                               pa=70*np.pi/180)
    aper = EllipticalAperture((geometry.x0, geometry.y0), geometry.sma,
                               geometry.sma*(1 - geometry.eps), geometry.pa)
    
    # pyp.imshow(data, origin='lower')
    # aper.plot(color='white')
    
    ellipse = Ellipse(data, geometry)
    isolist = ellipse.fit_image()
    
    # print(isolist.tflux_e)
    
    isophotes = True
    if isophotes :
        plt.display_isophotes(data, isolist, cmap=cm.viridis)
        
        # from photutils.isophote import build_ellipse_model
        # model = build_ellipse_model(data.shape, isolist)
        # residual = data - model
        # plt.display_image_simple(residual, norm=None)
    
    annuli = True
    if annuli :
        from photutils import EllipticalAnnulus
        from photutils import aperture_photometry
        
        center = (isolist[0].x0, isolist[0].y0)
        # print(center)
        
        last = np.where(isolist.stop_code == 0)[0][-1]
        isolist = isolist[last]
        
        pa = isolist.pa
        # print(pa*180/np.pi)
        
        a_outs = np.arange(1e-5, isolist.sma, isolist.sma/11)
        b_outs = a_outs*(1-isolist.eps)
        
        for i in range(len(a_outs) - 1) :
            a_in = a_outs[i]
            a_out = a_outs[i+1]
            b_in = b_outs[i]
            b_out = b_outs[i+1]
            # print(a_in, a_out, b_in, b_out)
            aper = EllipticalAnnulus(center, a_in, a_out, b_out, b_in=b_in,
                                     theta=isolist.pa)
            phot_table = aperture_photometry(data, aper, error=noise)
            flux = phot_table['aperture_sum'][0]
            flux_err = phot_table['aperture_sum_err'][0]
            # print(flux)
            # print(flux_err)
            
            annulus_mask = aper.to_mask()
            
            annulus_data = annulus_mask.multiply(data)
            # plt.display_image_simple(annulus_data, cmap=cm.viridis, norm=None)
            # print(np.sum(annulus_data))
            
            err_table = aperture_photometry(noise, aper)
            flux_err_alt = err_table['aperture_sum'][0]
            # print(flux_err)
            
            err_data = annulus_mask.multiply(noise)
            # print(np.sum(err_data))
            
            # print(flux/flux_err)
            print(flux/flux_err_alt)
            # print()
            
    return

def display_isophotes(data, isolist, bad='black', cbar_label='', cmap=cm.gray, 
                      norm=LogNorm(), vmin=None, vmax=None) :
    
    global currentFig
    fig = plt.figure(currentFig)
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    cmap = cmap
    cmap.set_bad(bad, 1)
    
    frame = ax.imshow(data, origin='lower', cmap=cmap, norm=norm,
                      vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(frame)
    cbar.set_label(cbar_label)
    
    smas = isolist.sma
    for i in range(len(isolist)) :
        if (isolist.stop_code[i] == 0) and (isolist.nflag[i] == 0) :
            iso = isolist.get_closest(isolist.sma[i])
            x, y, = iso.sampled_coordinates()
            ax.plot(x, y, color='white')
    
    plt.tight_layout()
    plt.show()
    
    return
