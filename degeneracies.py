
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

def degenerate_results(cluster, save=False) :
    
    if save :
        os.makedirs('{}/pngs_degen'.format(cluster), # ensure the output
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
    '''
    
    if cluster == 'a370' :
        phots = ['a370/photometry/a370_ID_3337_photometry.fits']
    
    if cluster == 'a1063' :
        phots = ['a1063/photometry/a1063_ID_1366_photometry.fits',
                 'a1063/photometry/a1063_ID_2455_photometry.fits',
                 'a1063/photometry/a1063_ID_3550_photometry.fits',
                 'a1063/photometry/a1063_ID_4823_photometry.fits']
    
    if cluster == 'a2744' :
        phots = ['a2744/photometry/a2744_ID_3859_photometry.fits',
                 'a2744/photometry/a2744_ID_3964_photometry.fits',
                 'a2744/photometry/a2744_ID_4173_photometry.fits',
                 'a2744/photometry/a2744_ID_4369_photometry.fits',
                 'a2744/photometry/a2744_ID_4765_photometry.fits',
                 'a2744/photometry/a2744_ID_4862_photometry.fits',
                 'a2744/photometry/a2744_ID_7427_photometry.fits']
    
    if cluster == 'm416' :
        phots = ['m416/photometry/m416_ID_5997_photometry.fits',
                 'm416/photometry/m416_ID_6255_photometry.fits']
    
    if cluster == 'm717' :
        phots = ['m717/photometry/m717_ID_861_photometry.fits']
    
    if cluster == 'm1149' :
        phots = ['m1149/photometry/m1149_ID_1967_photometry.fits',
                 'm1149/photometry/m1149_ID_2403_photometry.fits',
                 'm1149/photometry/m1149_ID_3531_photometry.fits',
                 'm1149/photometry/m1149_ID_4246_photometry.fits',
                 'm1149/photometry/m1149_ID_5095_photometry.fits']
    
    # loop over all the fits files in the directory
    for file in phots :
        ID = file.split('_')[2] # the galaxy ID to fit the bins for
        
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
        
        dusts, dusts_16, dusts_84, dust_best = [], [], [], []
        metals, metals_16, metals_84, metal_best = [], [], [], []
        mwas, mwa_16, mwa_84, mwa_best = [], [], [], []
        for binNum in bins : # loop over all the bins in the table
            infile = '{}/h5/ID_{}_bin_{}.h5'.format(cluster, ID, binNum)
            result, obs, _ = reader.results_from(infile, dangerous=True)
            zred, mass, logzsol, dust2, tage, tau = result['bestfit']['parameter']
            
            dust_best.append(dust2)            
            metal_best.append(logzsol)
            mwa_best.append(mwa_calc(tau, tage, power=1))
            
            # now pull samples to compute 1sigma errors
            samples = sample_posterior(result['chain'],
                                       weights=result['weights'])
            
            sample_dusts = samples[:, 3]
            sample_metals = samples[:, 2]
            sample_mwas = mwa_calc(samples[:, 5], samples[:, 4], power=1)
            
            dusts_16.append(np.percentile(sample_dusts, 16))
            dusts.append(np.percentile(sample_dusts, 50))
            dusts_84.append(np.percentile(sample_dusts, 84))
            
            metals_16.append(np.percentile(sample_metals, 16))
            metals.append(np.percentile(sample_metals, 50))
            metals_84.append(np.percentile(sample_metals, 84))
            
            mwa_16.append(np.percentile(sample_mwas, 16))
            mwas.append(np.percentile(sample_mwas, 50))
            mwa_84.append(np.percentile(sample_mwas, 84))
        
        dust_lo = np.abs(dusts - np.array(dusts_16))
        dust_hi = np.abs(np.array(dusts_84) - dusts)
        
        metal_lo = np.abs(metals - np.array(metals_16))
        metal_hi = np.abs(np.array(metals_84) - metals)
        
        mwa_lo = np.abs(mwas - np.array(mwa_16))
        mwa_hi = np.abs(np.array(mwa_84) - mwas)
        
        outfile = '{}/pngs_degen/{}_ID_{}.pdf'.format(cluster, cluster, ID)
        q_xi, q_yi, q_z, sf_xi, sf_yi, sf_z = checks.load_FUVVJ_contours()
        plt.plot_degeneracies([xs, xs, xs], [x_errs, x_errs, x_errs],
                              [x_errs, x_errs, x_errs], [dusts, metals, mwas],
                              [dust_lo, metal_lo, mwa_lo],
                              [dust_hi, metal_hi, mwa_hi],
                              [dust_best, metal_best, mwa_best],
                              q_xi, q_yi, q_z, sf_xi, sf_yi, sf_z,
                              V_mag-J_mag, FUV_mag-V_mag,
                              labels=['best fit/MAP', r'median$\pm 1\sigma$'],
                              xlabel=r'Center of Annulus ($R_{\rm e}$)',
                              list_of_ylabels=[r'Dust ($\tau_{\lambda, 5500 {\rm \AA}}$)',
                                               r'$\log(Z/Z_{\odot})$',
                                               'MWA (Gyr)'], title=plot_title,
                              xmin=0, xmax=(xs[-1]+x_errs[-1]),
                              outfile=outfile, save=save, loc='lower right')
    
    return
