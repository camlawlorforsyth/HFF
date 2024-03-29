
import os
import glob
import numpy as np

from astropy.table import Table
import matplotlib.pyplot as plt

import prospect.io.read_results as reader
from prospect.plotting.sfh import parametric_mwa as mwa_calc
from prospect.plotting.sfh import parametric_mwa_numerical as mwa_calc_alt
from prospect.plotting.utils import sample_posterior

from prospect.sources.constants import cosmo

currentFig = 1

def plot_mwa_metal(list_of_xs,
                   list_of_ys,
                   xerr=None,
                   yerr_lo=None, yerr_hi=None,
                   label1='',
                   xs2=None, ys2=None, y2_lo=None, y2_hi=None, label2='',
                   xlabel=None, ylabel=None, title=None,
                   xmin=None, xmax=None, ymin=None, ymax=None,
                   figsizewidth=9, figsizeheight=6) :
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(figsizewidth, figsizeheight))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    ax.plot(list_of_xs[0], list_of_ys[0], 'k-', zorder=1)
    ax.errorbar(list_of_xs[0], list_of_ys[0],
                xerr=xerr,
                yerr=(yerr_lo, yerr_hi),
                fmt='k.', label=label1, zorder=2, elinewidth=1)
    # ax.set_ylabel(ylabel, fontsize=15, color=colors[0])
    
    secax = ax.twinx()
    secax.set_ylabel('Metallicity', fontsize=15)
    secax.set_yticks([-3, -2, -1, 0, 1])
    
    secax.plot(xs2, ys2, 'r-', zorder=1)
    secax.errorbar(xs2, ys2, xerr=xerr, yerr=(abs(ys2-y2_lo), abs(y2_hi-ys2)),
                   fmt='r.', label=label2, zorder=2, elinewidth=1)
    
    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_title(title, fontsize=18)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    # ax.legend(facecolor='whitesmoke', framealpha=1, fontsize=15)
    
    plt.tight_layout()
    plt.show()
    
    return

# clusters = ['a370', 'a1063', 'a2744', 'm416', 'm717', 'm1149']
clusters = ['a370'] #, 'a1063']

for cluster in clusters :
    
    if cluster == 'a370' :
        phots = ['a370/photometry/a370_ID_3337_photometry.fits']
    
    if cluster == 'a1063' :
        phots = ['a1063/photometry/a1063_ID_1366_photometry.fits',
                 'a1063/photometry/a1063_ID_2455_photometry.fits',
                 'a1063/photometry/a1063_ID_3550_photometry.fits']
    
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
    
    # loop over all the fits files in the directory
    for file in phots :
        ID = file.split('_')[2] # the galaxy ID to fit the bins for
        
        table = Table.read(file)
        bins = table['bin'] # get a list of bin values
        xs = (table['sma'] - table['width']/2)/table['R_e']
        xs_err = (table['width']/2)/table['R_e']
        
        mwas, mwa_16, mwa_84 = [], [], []
        metals, metals_16, metals_84 = [], [], []
        for binNum in bins : # loop over all the bins in the table
            infile = '{}/h5/ID_{}_bin_{}.h5'.format(cluster, ID, binNum)
            result, obs, _ = reader.results_from(infile, dangerous=True)
            zred, mass, logzsol, dust2, tage, tau = result['bestfit']['parameter']
            
            metals.append(logzsol)
            mwa = mwa_calc(tau, tage, power=1)
            mwas.append(mwa)
            
            # print(mwa, mwa_calc_alt(tau, tage, power=1))
            # print(cosmo.age(zred).value, tage, tau, mwa)
            
            samples = sample_posterior(result['chain'],
                                       weights=result['weights'])
            
            sample_metals = samples[:, 2]
            sample_mwas = mwa_calc(samples[:, 5], samples[:, 4], power=1)
            
            mwa_16.append(np.percentile(sample_mwas, 16))
            mwa_84.append(np.percentile(sample_mwas, 84))
            
            metals_16.append(np.percentile(sample_metals, 16))
            metals_84.append(np.percentile(sample_metals, 84))
        
        mwa_lo = np.abs(mwas - np.array(mwa_16))
        mwa_hi = np.abs(np.array(mwa_84) - mwas)
        
        metals_lo = np.abs(metals - np.array(metals_16))
        metals_hi = np.abs(np.array(metals_84) - metals)
        
        print(mwa)
        print(mwa_lo)
        
        
        plot_mwa_metal([xs], [mwas], xerr=xs_err, yerr_lo=mwa_lo, yerr_hi=mwa_hi,
                       label1='MWA',
                       xs2=xs, ys2=metals,
                       y2_lo=np.array(metals_lo),
                       y2_hi=np.array(metals_hi),
                       label2='Metallicity',
                       xlabel=r'Center of Annulus ($R_{\rm e}$)',
                       ylabel='Value',
                       title='{} ID {}'.format(cluster, ID))
        
        
        '''
        plot_simple(xs, mwas,
                    yerr_lo=np.array(mwa_lo),
                    yerr_hi=np.array(mwa_hi),
                    xlabel=r'Semi Major Axis ($R_{\rm e}$)',
                    ylabel='Mass-Weighted Age (Gyr)',
                    # ylabel='Metallicity',
                    title='{} ID {}'.format(cluster, ID))
        '''
