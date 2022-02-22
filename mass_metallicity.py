
import os
import numpy as np

from astropy.io import fits
from astropy.table import Table
import prospect.io.read_results as reader
from prospect.plotting.utils import sample_posterior
from scipy.optimize import curve_fit

import core
import plotting as plt

def compare_masses() :
    
    shipley = Table.read('output/tables/nbCGs_with-Shipley-mass.fits')
    fitted = Table.read(
        'output/tables/nbCGs_integrated-logM-logZ-from-fitting.fits')
    
    fittedMass, shipleyMass = fitted['logM'], shipley['logM']
    
    mask = ~np.isnan(fittedMass)
    shipleyMass, fittedMass = shipleyMass[mask], fittedMass[mask]
    
    popt, pcov = curve_fit(core.linear, shipleyMass, fittedMass - shipleyMass)
        # p0=[0.4, 8, 0.5] # initial parameters for `core.exp` curve
    
    xs = np.linspace(8.3, 11.5, 1000)
    ys = core.linear(xs, *popt)
    
    # compare the masses from the Shipley table to the integrated fitted ones
    plt.plot_simple_multi([shipleyMass, xs], [fittedMass - shipleyMass, ys],
                          ['', 'fit'], ['k', 'r'], ['o', ''], ['', '-'],
                          xlabel=r'$\log(M/M_{\odot})_{\rm DeepSpace}$',
                          ylabel=(r'$\log(M/M_{\odot})_{\rm fit} - ' +
                                  r'\log(M/M_{\odot})_{\rm DeepSpace}$'),
                          scale='linear')#,
                          # xmin=8, xmax=11.4, ymin=-0.4, ymax=0.7)
    
    # compare the mass from the Shipley table to the integrated fitted ones
    # including a correction based on the exponential fit above
    y_correction = core.linear(shipleyMass, *popt)
    plt.plot_simple_multi([shipleyMass + y_correction, xs], [fittedMass, xs],
                          ['', 'equality'], ['k', 'r'], ['o', ''], ['', '--'],
                          xlabel=r'$\log(M/M_{\odot})_{\rm DeepSpace, corrected}$',
                          ylabel=r'$\log(M/M_{\odot})_{\rm fit}$',
                          scale='linear',
                          xmin=8.3, xmax=11.5, ymin=8.3, ymax=11.5)
    
    # compare the histograms of the masses from the Shipley table to the
    # integrated fitted ones
    plt.histogram_multi([fittedMass, shipleyMass], r'$\log(M/M_{\odot})$',
                        [20, 20], colors=['k', 'r'],
                        labels=['fit', 'DeepSpace'], styles=['-', '-'],
                        histtype='step', loc=2)
    
    # investigate the histogram of the integrated fitted metallicities
    plt.histogram(fitted['avglogZ'],
                  r'$\langle \log(Z/Z_{\odot}) \rangle_{L}$',
                  bins=20, histtype='step')
    
    return

def determine_expected_integrated_logZ() :
    
    infile = 'output/tables/Gallazzi_SDSS_mass-metallicity_with_SFR.fits'    
    sdss = Table.read(infile) # uses the metallicities from Gallazzi+2005
        # based on SDSS DR4, but updated total mass estimates from fits to DR7
        # photometry, along with SFRs from DR7
    
    # mask the data based on reasonable stellar masses
    sdss = sdss[sdss['MASS_P50'] > 8]
    
    # mask the table based on reasonable solar metallicities
    sdss = sdss[sdss['METALLICITY_P50'] - np.log10(0.02) > -2]
    
    # mask the table based on SFR
    sdss = sdss[sdss['SFR_P50'] > -99]
    
    # mask the table based on sSFR (ie. to select quiescent galaxies)
    sdss = sdss[sdss['SFR_P50'] - sdss['MASS_P50'] < -11]
    
    # now select only galaxies with sufficient signal to noise in their pixels
    # table = table[table['SN_MEDIAN'] > 3]
    
    # populate values from the table
    mass = sdss['MASS_P50']
    metal = sdss['METALLICITY_P50'] - np.log10(0.02) # convert to solar values
    
    # define the edges for the mass bins
    edges = [8, 8.2, 8.4, 8.6, 8.8, 9, 9.2, 9.4, 9.6, 9.8,
             10, 10.2, 10.4, 10.6, 10.8, 11, 11.2, 11.4, 11.6, 11.8, 12]
    
    tinies, los, medians, his, larges = [], [], [], [], []
    centers, nGals = [], []
    for first, second in zip(edges, edges[1:]) :
        
        centers.append(np.mean([first, second]))
        
        mass_mask = (mass >= first) & (mass < second)
        
        metallicities_in_bin = metal[mass_mask]
        nGals.append(len(metallicities_in_bin))
        
        tiny, lo, median, hi, large = np.percentile(metallicities_in_bin,
                                                    [2.5, 16, 50, 84, 97.5])
        
        tinies.append(tiny)
        los.append(lo)
        medians.append(median)
        his.append(hi)
        larges.append(large)
    
    centers = np.array(centers)
    los, medians, his = np.array(los), np.array(medians), np.array(his)
    
    popt, pcov = curve_fit(core.tanh, centers[centers > 8.7],
                           medians[centers > 8.7],
                           p0=[-0.452, 0.572, 9.66, 1.04]) # initial parameters
        # come from Panter+2008, MNRAS, 391, 1117
    
    fitx = np.linspace(8, 12, 4001)
    fity = core.tanh(fitx, *popt)
    
    HFF = Table.read(
        'output/tables/nbCGs_integrated-logM-logZ-from-fitting.fits')
    
    # plot the Gallazzi relation for quiescent galaxies, and overplot
    # the integrated values from the zeroth fitting run, ie. integrated values
    # from using the `params_earlyJan2022.py` file
    plt.histogram_2d(mass, metal, HFF['logM'], HFF['avglogZ'],
                     [centers, centers, centers, centers, centers],
                     [tinies, los, medians, his, larges],
                     [fitx], [fity], ['', 'fit', '', 'HFF'],
                     [':', '--', '-', '--', ':'],
                     xlabel=r'$\log(M/M_{\odot})$',
                     ylabel=r'$\log(Z/Z_{\odot})$',
                     xmin=8, xmax=12, ymin=-1.6, ymax=0.4)
    
    # now fit the +/- 1 sigma lines using galaxies with logM > 8.7
    popt_lo, pcov_16 = curve_fit(core.tanh, centers[centers > 8.7],
                                 los[centers > 8.7],
                                 p0=[-0.53348847, 0.46418104,
                                      9.72088385, 0.93802711])
    fity_lo = core.tanh(fitx, *popt_lo)
    
    popt_hi, pcov_hi = curve_fit(core.tanh, centers[centers > 8.7],
                                 his[centers > 8.7],
                                 p0=[0.02577337, 0.24113023,
                                     9.9563868,  1.34806179])
    fity_hi = core.tanh(fitx, *popt_hi)
    
    # we'll use the maximum offset as the offset on either side
    sigma = np.maximum(fity_hi - fity, fity - fity_lo)
    
    # now plot the derived relations with the fitted relations
    plt.plot_simple_multi([centers, centers, centers, centers, centers,
                           fitx, fitx, fitx, fitx, fitx],
                          [tinies, los, medians, his, larges,
                           fity_lo, fity, fity_hi, fity - sigma, fity + sigma],
                          ['', '', '', '', '', '', r'${\rm fit}_{16/50/84th}$',
                           '', r'${\rm fit}_{50th} \pm \sigma$', ''],
                          ['k', 'k', 'k', 'k', 'k',
                           'r', 'r', 'r', 'b', 'b'],
                          ['', '', '', '', '',
                           '', '', '', '', ''],
                          [':', '--', '-', '--', ':',
                           '-', '-', '-', '--', '--'],
                          xlabel=r'$\log(M/M_{\odot})$',
                          ylabel=r'$\log(Z/Z_{\odot})$',
                          xmin=8, xmax=12, ymin=-1.6, ymax=0.4, scale='linear',
                          loc=4)
    
    # now use the fitted relations above to derive expected values for the HFF
    shipley = Table.read('output/tables/nbCGs_with-Shipley-mass.fits')
    
    expected = core.tanh(shipley['logM'], *popt)
    expected_lo = core.tanh(shipley['logM'], *popt_lo)
    expected_hi = core.tanh(shipley['logM'], *popt_hi)
    
    shipley['GallazzilogZ'] = expected
    shipley['sigma'] = np.maximum(expected_hi-expected, expected-expected_lo)
    # shipley.write('output/tables/nbCGs_GallazzilogZ_from-Shipley-mass.fits')
    
    return

def estimate_gradient_from_literature() :
    
    califa_logM_lo = np.array([10.6, 10.5, 10.3])
    califa_logM_hi = np.array([11.8, 11.9, 11.9])
    
    califa_result = np.array([-0.1, -0.248, -0.2])
    califa_extent = np.array([1, 1, 2])
    
    sami_logM_lo = np.array([9.6, 9.5])
    sami_logM_hi = np.array([11.7, 11.7])
    
    sami_result = np.array([-0.31, -0.275])
    sami_extent = np.array([5, 2])
    
    manga_logM_lo = np.array([9, 9, 8.4, 9.1, 9.9, 9.9,
                              10.9, 9.5, 8.8, 9, 9.4])
    manga_logM_hi = np.array([11.9, 11.9, 11.9, 13.1, 10.8, 10.8,
                              12, 11.9, 11.3, 11.8, 12])
    
    manga_result = np.array([-0.12, -0.11, -0.09, -0.102, -0.14, -0.104,
                             -0.202, -0.112, -0.106, -0.092, -0.18])
    manga_extent = np.array([1.5, 1.5, 2, 1, 1, 1,
                             1, 1.5, 1.5, 1, 1])
    
    logM = np.concatenate([califa_logM_lo, sami_logM_lo, manga_logM_lo])
    results = np.concatenate([califa_result, sami_result, manga_result])
    
    xhi = np.concatenate([califa_logM_hi, sami_logM_hi, manga_logM_hi]) - logM
    xerr = np.array([np.zeros(xhi.shape), xhi])
    
    extents = np.concatenate([califa_extent, sami_extent, manga_extent])
    
    plt.plot_scatter_err(logM, results, xerr, extents, 'o',
                         cbar_label=r'$R_{\rm e}$ extent',
                         xlabel=r'Mass Range $(\log M/M_{\odot})$',
                         ylabel=r'$\nabla \log(Z/Z_{\odot})$')
    
    plt.plot_scatter(extents, results, logM, '', 'o',
                     cbar_label=r'Lowest Mass $(\log M/M_{\odot})$',
                     xlabel=r'$R_{\rm e}$ extent',
                     ylabel=r'$\nabla \log(Z/Z_{\odot})$')
    
    return

def luminosity_weighted_average_metallicity() :
    
    HFF = Table.read('output/tables/nbCGs.fits')
    
    masses, metallicities = [], []
    for cluster, ID in zip(HFF['cluster'], HFF['ID']) :
        if (cluster == 'm717') & (ID == 3692) : # fitting issue for bin_4
            masses.append(np.nan)
            metallicities.append(np.nan)    
        else :
            mass_image = core.open_image(cluster, ID, 'logM')
            masses.append(np.log10(np.nansum(mass_image)))
            
            fluxfrac = core.open_image(cluster, ID, 'fluxfrac')
            metal_image = core.open_image(cluster, ID, 'logZ')
            metallicities.append(np.nansum(fluxfrac*metal_image))
    
    HFF['logM'] = masses
    HFF['avglogZ'] = metallicities
    HFF.write('output/tables/nbCGs_integrated-logM-logZ-from-fitting_NEW.fits')
    
    return

def metallicity_normalization_estimation() :
    
    HFF = Table.read(
        'output/tables/nbCGs_GallazzilogZ_from-Shipley-mass.fits')
    
    vals = []
    for cluster, ID in zip(HFF['cluster'], HFF['ID']) :
        mask = (HFF['cluster'] == cluster) & (HFF['ID'] == ID)
        vals.append(metallicity_normalization_helper(
            cluster, ID, HFF['GallazzilogZ'][mask][0]))
    
    plt.plot_scatter(HFF['logM'], HFF['GallazzilogZ'] - np.array(vals),
                     HFF['bins'], 'HFF', 'o', cbar_label='Number of Annuli',
                     xlabel=r'$\log(M/M_{\odot})$',
                     ylabel=r'$\Delta \langle \log(Z/Z_{\odot}) \rangle$',
                     xmin=7.9, xmax=11.5, loc=2)
    
    offset = HFF['GallazzilogZ'] - np.array(vals)
    
    finals = []
    for cluster, ID in zip(HFF['cluster'], HFF['ID']) :
        mask = (HFF['cluster'] == cluster) & (HFF['ID'] == ID)
        finals.append(metallicity_normalization_helper(
            cluster, ID, HFF['GallazzilogZ'][mask][0] + offset[mask]))
    
    plt.plot_scatter(HFF['logM'], HFF['GallazzilogZ'] - np.array(finals),
                     HFF['bins'], 'HFF', 'o', cbar_label='Number of Annuli',
                     xlabel=r'$\log(M/M_{\odot})$',
                     ylabel=r'$\Delta \langle \log(Z/Z_{\odot}) \rangle$',
                     xmin=7.9, xmax=11.5, loc=2)
    
    HFF['centrallogZ'] = HFF['GallazzilogZ'] + offset
    # HFF.write('output/tables/nbCGs_GallazzilogZ_from-Shipley-mass_final.fits')
    
    return

def metallicity_normalization_helper(cluster, ID, central_value) :
    
    metal_image = core.open_image(cluster, ID, 'bins')
    fluxfrac = core.open_image(cluster, ID, 'fluxfrac')
    
    table = Table.read('{}/photometry/{}_ID_{}_photometry.fits'.format(
        cluster, cluster, ID))
    radii = np.array((table['sma'] - table['width'])*np.sqrt(
        table['smb']/table['sma'])/table['R_e'])
    
    for binNum, radius in zip(table['bin'], radii) :
        metal_image[metal_image == binNum] = central_value - 0.1*radius
    
    return np.nansum(fluxfrac*metal_image)

def save_mass_metallicity_images() :
    
    HFF = Table.read('output/tables/nbCGs.fits')
    
    for cluster, ID in zip(HFF['cluster'], HFF['ID']) :
        # ensure the output directories for the images are available
        os.makedirs('{}/logM_images'.format(cluster), exist_ok=True)
        os.makedirs('{}/logZ_images'.format(cluster), exist_ok=True)
        
        bins_image = core.open_image(cluster, ID, 'bins')
        bins = np.sort(np.unique(bins_image[~np.isnan(bins_image)])).astype(int)
        
        mass = bins_image.copy()
        metal = bins_image.copy()
        
        if not (cluster == 'm717') & (ID == 3692) : # fitting issue for bin_4
            for binNum in bins : # loop over all the bins
                result, obs, _ = reader.results_from(
                    '{}/h5/{}_ID_{}_bin_{}.h5'.format(
                        cluster, cluster, ID, binNum), dangerous=True)
                samples = sample_posterior(result['chain'],
                                           weights=result['weights'])
                
                mass[mass == binNum] = np.percentile(samples[:, 1], 50)
                metal[metal == binNum] = np.percentile(samples[:, 2], 50)
            
            hdu = fits.PrimaryHDU(mass)
            hdu.writeto('{}/logM_images/{}_ID_{}_logM.fits'.format(
                cluster, cluster, ID))
            
            hdu = fits.PrimaryHDU(metal)
            hdu.writeto('{}/logZ_images/{}_ID_{}_logZ.fits'.format(
                cluster, cluster, ID))
    
    return
