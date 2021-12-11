
import os
# import glob
import numpy as np

from astropy.table import Table
import astropy.units as u
import prospect.io.read_results as reader
from prospect.plotting.sfh import parametric_mwa as mwa_calc
from prospect.plotting.utils import sample_posterior
from prospect.sources.constants import cosmo

import checks
import plotting as plt

hst_pixelscale = u.pixel_scale(0.06*u.arcsec/u.pixel)

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
    
    file = '{}/photometry/{}_ID_{}_photometry.fits'.format(cluster, cluster, ID)
    
    table = Table.read(file)
    bins = table['bin'] # get a list of bin values
    R_e = table['R_e']
    xs = (table['sma'] - table['width']/2)/R_e
    x_errs = (table['width']/2)/R_e
    
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
    
    metals_16, metals_50, metals_84, metals_best = [], [], [], []
    dusts_16, dusts_50, dusts_84, dusts_best = [], [], [], []
    mwas_16, mwas_50, mwas_84, mwas_best = [], [], [], []
    for binNum in bins : # loop over all the bins in the table
        infile = '{}/h5/{}_ID_{}_bin_{}.h5'.format(cluster, cluster, ID, binNum)
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
    
    metal_slope, metal_int = np.polyfit(xs, metals_50, 1)
    dust_slope, dust_int = np.polyfit(xs, dusts_50, 1)
    mwa_slope, mwa_int = np.polyfit(xs, mwas_50, 1)
    print(mwa_slope)
    
    outfile = '{}/images_degen_plots/{}_ID_{}_new.pdf'.format(cluster, cluster, ID)
    q_xi, q_yi, q_z, sf_xi, sf_yi, sf_z = checks.load_FUVVJ_contours()
    '''
    plt.plot_degeneracies([xs, xs, xs], [x_errs, x_errs, x_errs],
                          [x_errs, x_errs, x_errs],
                          [metals_50, dusts_50, mwas_50],
                          [metal_lo, dust_lo, mwa_lo],
                          [metal_hi, dust_hi, mwa_hi],
                          [metals_best, dusts_best, mwas_best],
                          [metal_slope, dust_slope, mwa_slope],
                          [metal_int, dust_int, mwa_int],
                          q_xi, q_yi, q_z, sf_xi, sf_yi, sf_z,
                          V_mag-J_mag, FUV_mag-V_mag,
                          labels=['best fit/MAP', 'linear fit',
                                  r'median$\pm 1\sigma$'],
                          xlabel=r'Center of Annulus ($R_{\rm e}$)',
                          list_of_ylabels=[r'$\log(Z/Z_{\odot})$',
                                           r'Dust ($\hat{\tau}_{\lambda, 2}$)',
                                           'MWA (Gyr)'], title=plot_title,
                          xmin=0, xmax=(xs[-1]+x_errs[-1]),
                          outfile=outfile, save=False, loc=loc)
    '''
    return

degenerate_results('a370', 3337, loc=1)
degenerate_results('a1063', 1366, loc=3)
degenerate_results('a1063', 2455, loc=3)
degenerate_results('a1063', 3550, loc=1)
degenerate_results('a1063', 4823, loc=1)
degenerate_results('a2744', 3859, loc=3)
degenerate_results('a2744', 3964, loc=3)
degenerate_results('a2744', 4173, loc=3)
degenerate_results('a2744', 4369, loc=1)
degenerate_results('a2744', 4765, loc=1)
degenerate_results('a2744', 4862, loc=3)
degenerate_results('a2744', 7427, loc=4)
degenerate_results('m416', 5997, loc=2)
degenerate_results('m416', 6255, loc=1)
degenerate_results('m717', 861, loc=1)
degenerate_results('m1149', 1967, loc=3)
degenerate_results('m1149', 2403, loc=2)
degenerate_results('m1149', 3531, loc=4)
degenerate_results('m1149', 4246, loc=1)
degenerate_results('m1149', 5095, loc=2)
