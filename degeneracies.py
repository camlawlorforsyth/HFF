
import os
# import glob
import numpy as np

from matplotlib import cm

from astropy.table import Table
import astropy.units as u
from astropy.visualization import make_lupton_rgb
import prospect.io.read_results as reader
from prospect.plotting.sfh import parametric_mwa as mwa_calc
from prospect.plotting.utils import sample_posterior
from prospect.sources.constants import cosmo

from core import open_cutout
import checks
import plotting as plt

hst_pixelscale = u.pixel_scale(0.06*u.arcsec/u.pixel)

def save_age_gradients() :
    
    sample = Table.read('output/tables/nbCGs_with-Shipley-mass_currentlyFitted.fits')
    
    all_radii = []
    all_mwa = []
    for cluster, ID in zip(sample['cluster'], sample['ID']) :
        
        table = Table.read('{}/photometry/{}_ID_{}_photometry.fits'.format(
            cluster, cluster, ID))
        bins = table['bin'] # get a list of bin values
        sma, smb = table['sma'], table['smb']
        R_e, width = table['R_e'], table['width']
        
        radius_array = list((sma - width)*np.sqrt(smb/sma)/R_e)
        
        radii, mwas = determine_age_gradients(cluster, ID, bins, radius_array,
                                              version='_gradient')
        
        all_radii.append(radii)
        all_mwa.append(mwas)
    
    all_radii = np.concatenate(all_radii).ravel()
    all_mwa = np.concatenate(all_mwa).ravel()
    
    test_table = Table([all_radii, all_mwa], names=('radius', 'MWA'))
    test_table.write('output/tables/radii_and_MWA.fits')
    
    return

def plot_age_gradients() :
    
    sample = Table.read('output/tables/radii_and_MWA.fits')
    radii = sample['radius']
    MWA = sample['MWA']
    
    edges = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3]
    
    los, medians, his = [], [], []
    centers, nSamples = [], []
    for first, second in zip(edges, edges[1:]) :
        
        centers.append(np.mean([first, second]))
        
        mask = (radii >= first) & (radii < second)
        
        MWA_in_bin = MWA[mask]
        nSamples.append(len(MWA_in_bin))
        
        lo, median, hi = np.percentile(MWA_in_bin, [16, 50, 84])
        
        los.append(lo)
        medians.append(median)
        his.append(hi)
    
    centers = np.array(centers)
    los, medians, his = np.array(los), np.array(medians), np.array(his)
    
    plt.plot_simple_multi([centers, centers, centers], [los, medians, his],
                          ['', '', ''], ['k', 'k', 'k'], ['', '', ''], ['--', '-', '--'],
                          xlabel=r'Radius ($R_{\rm e}$)', ylabel='MWA (Gyr)',
                          xmin=0, xmax=3, scale='linear')
    
    return

def determine_age_gradients(cluster, ID, bins, radius_array, version='') :
    
    mwas = []
    radii = []
    for binNum, radius in zip(bins, radius_array) :
        infile = '{}/h5/{}_ID_{}_bin_{}{}.h5'.format(
            cluster, cluster, ID, binNum, version)
        result, obs, _ = reader.results_from(infile, dangerous=True)
        
        samples = sample_posterior(result['chain'], weights=result['weights'])
        
        mwas.append(mwa_calc(samples[:, 5], samples[:, 4], power=1))
        
        radii.append(np.full(len(samples), radius))
    
    return np.concatenate(radii).ravel(), np.concatenate(mwas).ravel()






def comparison(cluster, ID, loc=0) :
    
    table = Table.read('{}/photometry/{}_ID_{}_photometry.fits'.format(
        cluster, cluster, ID))
    bins = table['bin'] # get a list of bin values
    sma, smb = table['sma'], table['smb']
    R_e, width = table['R_e'], table['width']
    
    # xs = (sma - 0.5*width)/R_e
    # xerrs = 0.5*width/R_e
    
    xs = (sma - width)*np.sqrt(smb/sma)/R_e
    # xerrs_lo, xerrs_hi = np.zeros(len(xs)), width*np.sqrt(smb/sma)/R_e
    
    angSize = (R_e[0]*u.pix).to(u.arcsec, hst_pixelscale)
    physSize = angSize/cosmo.arcsec_per_kpc_comoving(table['z'][0])
    plot_title = ('{} ID {}'.format(cluster, ID) +
                  r' [$R_{\rm e}$' +
                  ' = {:.3f} = {:.3f} at z = {}]'.format(angSize,
                                                         physSize.to(u.kpc),
                                                         table['z'][0]))
    
    (Z_flat_median, Z_flat_lo, Z_flat_hi,
     dust_flat_median, dust_flat_lo, dust_flat_hi,
     mwa_flat_median, mwa_flat_lo, mwa_flat_hi) = determine_lines(
         cluster, ID, bins, version='_flat')
    
    (Z_grad_median, Z_grad_lo, Z_grad_hi,
     dust_grad_median, dust_grad_lo, dust_grad_hi,
     mwa_grad_median, mwa_grad_lo, mwa_grad_hi) = determine_lines(
         cluster, ID, bins, version='_gradient')
    
    xs = [xs, xs+0.05, xs, xs+0.05]
    # xs = [xs, xs, xs, xs]
    ys = [Z_flat_median, Z_grad_median, mwa_flat_median, mwa_grad_median]
    lo = [Z_flat_lo, Z_grad_lo, mwa_flat_lo, mwa_grad_lo]
    hi = [Z_flat_hi, Z_grad_hi, mwa_flat_hi, mwa_grad_hi]
    labels = [r'$\nabla = 0$', r'$\nabla = -0.1$',
              r'$\nabla = 0$', r'$\nabla = -0.1$']
    colors = ['red', 'darkred', 'grey', 'darkgrey']
    styles = ['-', '--', '-', '--']
    
    plt.plot_degens_small(xs, ys, lo, hi, labels, colors, styles,
                          xlabel=r'Radius ($R_{\rm e}$)',
                          ylabels=[r'$\log(Z/Z_{\odot})$',
                                   r'Dust ($\hat{\tau}_{\lambda, 2}$)',
                                   'MWA (Gyr)'], 
                          title=plot_title,            
                          xmin=np.min(xs) - 0.05,
                          xmax=np.max(xs) + 0.05, loc=loc,
                          outfile='{}_ID_{}.pdf'.format(cluster, ID),
                          save=True)
    
    return

def degenerate_results(cluster, ID, loc=0, save=False, version='') :
    
    if save :
        os.makedirs('{}/images_degen_plots'.format(cluster), # ensure the output
                    exist_ok=True) # directory for the figures is available
    
    '''
    # get a list of fits files containing photometric data for all bins for
    # a given galaxy, as denoted by ID
    photometries = '{}/photometry/{}_ID_*_photometry.fits'.format(cluster,
                                                                  cluster)
    phot_files = glob.glob(photometries)
    
    phots = []
    for file in phot_files :
        file = file.replace(os.sep, '/') # compatibility for Windows
        phots.append(file)
    
    # loop over all the fits files in the directory
    for file in phots :
        ID = file.split('_')[2] # the galaxy ID to fit the bins for
    '''
    
    f435w = open_cutout('{}/cutouts/{}_ID_{}_f435w.fits'.format(
        cluster, cluster, ID), phot=True)
    f606w = open_cutout('{}/cutouts/{}_ID_{}_f606w.fits'.format(
        cluster, cluster, ID), phot=True)
    f814w = open_cutout('{}/cutouts/{}_ID_{}_f814w.fits'.format(
        cluster, cluster, ID), phot=True)
    f105w = open_cutout('{}/cutouts/{}_ID_{}_f105w.fits'.format(
        cluster, cluster, ID), phot=True)
    f125w = open_cutout('{}/cutouts/{}_ID_{}_f125w.fits'.format(
        cluster, cluster, ID), phot=True)
    f140w = open_cutout('{}/cutouts/{}_ID_{}_f140w.fits'.format(
        cluster, cluster, ID), phot=True)
    f160w = open_cutout('{}/cutouts/{}_ID_{}_f160w.fits'.format(
        cluster, cluster, ID), phot=True)
    
    rgb = make_lupton_rgb(7e7*(f125w + f140w + f160w),
                          7e7*(f814w + f105w),
                          7e7*(f435w + f606w))
    
    file = '{}/photometry/{}_ID_{}_photometry.fits'.format(cluster, cluster, ID)
    
    table = Table.read(file)
    bins = table['bin'] # get a list of bin values
    sma, smb = table['sma'], table['smb']
    R_e, width = table['R_e'], table['width']
    
    # xs = (sma - 0.5*width)/R_e
    # xerrs = 0.5*width/R_e
    
    xs = (sma - width)*np.sqrt(smb/sma)/R_e    
    xerrs_lo, xerrs_hi = np.zeros(len(xs)), width*np.sqrt(smb/sma)/R_e
    
    angSize = (R_e[0]*u.pix).to(u.arcsec, hst_pixelscale)
    physSize = angSize/cosmo.arcsec_per_kpc_comoving(table['z'][0])
    plot_title = ('{} ID {}'.format(cluster, ID) +
                  r' [$R_{\rm e}$' +
                  ' = {:.3f} = {:.3f} at z = {}]'.format(angSize,
                                                         physSize.to(u.kpc),
                                                         table['z'][0]))
    
    FUV_mag = table['FUV_mag'][0]
    V_mag = table['V_mag'][0]
    J_mag = table['J_mag'][0]
    
    (metal_median, metal_lo, metal_hi,
     dust_median, dust_lo, dust_hi,
     mwa_median, mwa_lo, mwa_hi,
     metal_best, dust_best, mwa_best) = determine_lines(cluster, ID, bins,
                                                        best=True,
                                                        version=version)
    
    metal_slope, metal_int = np.polyfit(xs, metal_median, 1)
    dust_slope, dust_int = np.polyfit(xs, dust_median, 1)
    mwa_slope, mwa_int = np.polyfit(xs, mwa_median, 1)
    # print(mwa_slope)
    
    outfile = '{}/images_degen_plots/{}_ID_{}_new.pdf'.format(cluster, cluster, ID)
    q_xi, q_yi, q_z, sf_xi, sf_yi, sf_z = checks.load_FUVVJ_contours()
    
    plt.plot_degeneracies([xs, xs, xs],
                          [xerrs_lo, xerrs_lo, xerrs_lo],
                          [xerrs_hi, xerrs_hi, xerrs_hi],
                          [metal_median, dust_median, mwa_median],
                          [metal_lo, dust_lo, mwa_lo],
                          [metal_hi, dust_hi, mwa_hi],
                          [metal_best, dust_best, mwa_best],
                          [metal_slope, dust_slope, mwa_slope],
                          [metal_int, dust_int, mwa_int],
                          q_xi, q_yi, q_z, sf_xi, sf_yi, sf_z,
                          V_mag-J_mag, FUV_mag-V_mag, rgb,
                          labels=['best fit/MAP', 'linear fit',
                                  r'median$\pm 1\sigma$'],
                          xlabel=r'Radius ($R_{\rm e}$)',
                          list_of_ylabels=[r'$\log(Z/Z_{\odot})$',
                                           r'Dust ($\hat{\tau}_{\lambda, 2}$)',
                                           'MWA (Gyr)'], title=plot_title,
                          xmin=-0.05, xmax=(xs[-1] + xerrs_hi[-1] + 0.05),
                          outfile=outfile, save=False, loc=loc)
    
    return

def determine_lines(cluster, ID, bins, best=False, version='') :
    
    metals_16, metals_50, metals_84 = [], [], []
    dusts_16, dusts_50, dusts_84 = [], [], []
    mwas_16, mwas_50, mwas_84 = [], [], []
    
    metals_best, dusts_best, mwas_best = [], [], []
    
    for binNum in bins : # loop over all the bins in the table
        infile = '{}/h5/{}_ID_{}_bin_{}{}.h5'.format(
            cluster, cluster, ID, binNum, version)
        result, obs, _ = reader.results_from(infile, dangerous=True)
        
        zred, mass, logzsol, dust2, tage, tau = result['bestfit']['parameter']
        
        dusts_best.append(dust2)            
        metals_best.append(logzsol)
        mwas_best.append(mwa_calc(tau, tage, power=1))
        
        # now pull samples to compute 1sigma errors
        samples = sample_posterior(result['chain'],
                                   weights=result['weights'])
        
        metal_16, metal_50, metal_84 = np.percentile(samples[:, 2],
                                                     [16, 50, 84])
        metals_16.append(metal_16)
        metals_50.append(metal_50)
        metals_84.append(metal_84)
        
        dust_16, dust_50, dust_84 = np.percentile(samples[:, 3],
                                                  [16, 50, 84])
        dusts_16.append(dust_16)
        dusts_50.append(dust_50)
        dusts_84.append(dust_84)
        
        mwa_16, mwa_50, mwa_84 = np.percentile(mwa_calc(samples[:, 5],
                                                        samples[:, 4],
                                                        power=1),
                                               [16, 50, 84])
        mwas_16.append(mwa_16)
        mwas_50.append(mwa_50)
        mwas_84.append(mwa_84)
    
    metal_lo = np.abs(np.array(metals_50) - np.array(metals_16))
    metal_hi = np.abs(np.array(metals_84) - np.array(metals_50))
    
    dust_lo = np.abs(np.array(dusts_50) - np.array(dusts_16))
    dust_hi = np.abs(np.array(dusts_84) - np.array(dusts_50))
    
    mwa_lo = np.abs(np.array(mwas_50) - np.array(mwas_16))
    mwa_hi = np.abs(np.array(mwas_84) - np.array(mwas_50))
    
    metal_median = np.array(metals_50)
    dust_median = np.array(dusts_50)
    mwa_median = np.array(mwas_50)
    
    metal_best = np.array(metals_best)
    dust_best = np.array(dusts_best)
    mwa_best = np.array(mwas_best)
    
    if best :
        return (metal_median, metal_lo, metal_hi,
                dust_median, dust_lo, dust_hi,
                mwa_median, mwa_lo, mwa_hi,
                metal_best, dust_best, mwa_best)
    else :
        return (metal_median, metal_lo, metal_hi,
                dust_median, dust_lo, dust_hi,
                mwa_median, mwa_lo, mwa_hi)

# degenerate_results('a370', 3337, loc=1)
# degenerate_results('a1063', 1366, loc=3)
# degenerate_results('a1063', 2455, loc=3)
# degenerate_results('a1063', 3550, loc=1)
# degenerate_results('a1063', 4823, loc=1)
# degenerate_results('a2744', 3859, loc=3)
# degenerate_results('a2744', 3964, loc=3)
# degenerate_results('a2744', 4173, loc=3)
# degenerate_results('a2744', 4369, loc=1)
# degenerate_results('a2744', 4765, loc=1)
# degenerate_results('a2744', 4862, loc=3)
# degenerate_results('a2744', 7427, loc=4)
# degenerate_results('m416', 5997, loc=2)
# degenerate_results('m416', 6255, loc=1)
# degenerate_results('m717', 861, loc=1)
# degenerate_results('m1149', 1967, loc=3)
# degenerate_results('m1149', 2403, loc=2)
# degenerate_results('m1149', 3531, loc=4)
# degenerate_results('m1149', 4246, loc=1)
# degenerate_results('m1149', 5095, loc=2)

# degenerate_results('a1063', 4823, loc=1)
# degenerate_results('a1063', 5156, loc=4)
# degenerate_results('a1063', 5771, loc=4)
# degenerate_results('a2744', 4369, loc=1)
# degenerate_results('m416', 2452, loc=2)
# degenerate_results('m416', 3638, loc=1)
# degenerate_results('m416', 4030, loc=1)
# degenerate_results('m416', 5607, loc=1)
# degenerate_results('m717', 2786, loc=1)
# degenerate_results('m1149', 4322, loc=1)

# comparison('a1063', 4823, loc=[0.6, 0.8, 0.3, 0.2]) # the ten logM~10 galaxies
# comparison('a1063', 5156, loc=4)
# comparison('a1063', 5771, loc=2)
# comparison('a2744', 4369, loc=3)
# comparison('m416', 2452, loc=2)
# comparison('m416', 3638, loc=[0.6, 0.8, 0.3, 0.2])
# comparison('m416', 4030, loc=[0.65, 0.8, 0.3, 0.2])
# comparison('m416', 5607, loc=3)
# comparison('m717', 2786, loc=1)
# comparison('m1149', 4322, loc=4)

# comparison('a370', 3337, loc=1) # the twenty "representative" galaxies
# comparison('a1063', 1366, loc=1)
# comparison('a1063', 2455, loc=3)
# comparison('a1063', 3550, loc=1)
# comparison('a1063', 4823, loc=[0.6, 0.8, 0.3, 0.2])
# comparison('a2744', 3859, loc=3)
# comparison('a2744', 3964, loc=3)
# comparison('a2744', 4173, loc=[0.65, 0.8, 0.3, 0.2])
# comparison('a2744', 4369, loc=3)
# comparison('a2744', 4765, loc=1)
# comparison('a2744', 4862, loc=3)
# comparison('a2744', 7427, loc=4)
# comparison('m416', 5997, loc=[0.6, 0.8, 0.3, 0.2])
# comparison('m416', 6255, loc=1)
# comparison('m717', 861, loc=[0.763, 0.8, 0.2, 0.2])
# comparison('m1149', 1967, loc=3)
# comparison('m1149', 2403, loc=1)
# comparison('m1149', 3531, loc=4)
# comparison('m1149', 4246, loc=1)
# comparison('m1149', 5095, loc=[0.6, 0.8, 0.3, 0.2])

# new before committee meeting
# degenerate_results('a370', 3337, loc=1, version='_gradient') # bad
# degenerate_results('a1063', 1366, loc=3, version='_gradient') # decent
# degenerate_results('a1063', 2455, loc=3, version='_gradient') # bad
# degenerate_results('a1063', 3550, loc=1, version='_gradient') # bad
# degenerate_results('a1063', 4823, loc=1, version='_gradient') # decent
# degenerate_results('a2744', 3859, loc=3, version='_gradient') # saved
# degenerate_results('a2744', 3964, loc=3, version='_gradient') # saved
# degenerate_results('a2744', 4173, loc=3, version='_gradient') # saved
# degenerate_results('a2744', 4369, loc=1, version='_gradient') # saved
# degenerate_results('a2744', 4765, loc=1, version='_gradient') # saved
# degenerate_results('a2744', 4862, loc=3, version='_gradient') # saved
# degenerate_results('a2744', 7427, loc=4, version='_gradient') # saved
# degenerate_results('m416', 5997, loc=2, version='_gradient') # saved
# degenerate_results('m416', 6255, loc=1, version='_gradient') # saved
# degenerate_results('m717', 861, loc=1, version='_gradient') # saved
# degenerate_results('m1149', 1967, loc=3, version='_gradient') # saved
# degenerate_results('m1149', 2403, loc=2, version='_gradient') # bad
# degenerate_results('m1149', 3531, loc=4, version='_gradient') # bad
# degenerate_results('m1149', 4246, loc=1, version='_gradient') # saved
# degenerate_results('m1149', 5095, loc=2, version='_gradient') # bad

# degenerate_results('a1063', 4823, loc=1, version='_gradient') # bad
# degenerate_results('a1063', 5156, loc=4, version='_gradient') # decent
# degenerate_results('a1063', 5771, loc=4, version='_gradient') # bad
# degenerate_results('a2744', 4369, loc=1, version='_gradient') # decent
# degenerate_results('m416', 2452, loc=2, version='_gradient') # decent
# degenerate_results('m416', 3638, loc=1, version='_gradient') # decent
# degenerate_results('m416', 4030, loc=1, version='_gradient') # bad
# degenerate_results('m416', 5607, loc=1, version='_gradient') # saved (main)
# degenerate_results('m717', 2786, loc=1, version='_gradient') # decent
# degenerate_results('m1149', 4322, loc=1, version='_gradient') # bad

# plot_age_gradients()
