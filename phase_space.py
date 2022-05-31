
import numpy as np

import astropy.constants as c
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
import astropy.units as u
from matplotlib import cm
from scipy.optimize import curve_fit

import core
import plotting as plt

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def color_err_from_flux(flux1, flux1_err, flux2, flux2_err) :
    return np.array(2.5/np.log(10)*np.sqrt(
        np.square(flux1_err/flux1) + np.square(flux2_err/flux2)))

def color_from_flux(flux1, flux2) :
    return flux_to_mag(flux1) - flux_to_mag(flux2)

def compute_color_gradients(catalog) :
    
    dictionary = {'a370':{'lo':'f475w', 'hi':'f625w'},
                  'a1063':{'lo':'f435w', 'hi':'f625w'},
                  'a2744':{'lo':'f435w', 'hi':'f606w'}, # 'hi':'f625w' would be ideal,
                  'm416':{'lo':'f475w', 'hi':'f625w'},  # but is unavailable
                  'm717':{'lo':'f555w', 'hi':'f775w'},
                  'm1149':{'lo':'f555w', 'hi':'f775w'},
                  'a370par':{'lo':'f435w', 'hi':'f814w'},
                  'a1063par':{'lo':'f435w', 'hi':'f814w'},
                  'a2744par':{'lo':'f435w', 'hi':'f814w'},
                  'm416par':{'lo':'f435w', 'hi':'f775w'},
                  'm717par':{'lo':'f435w', 'hi':'f814w'},
                  'm1149par':{'lo':'f435w', 'hi':'f814w'}}
    
    gradients = []
    for cluster, ID in zip(catalog['cluster'], catalog['id']) :
        
        phot = Table.read('{}/photometry/{}_ID_{}_photometry.fits'.format(
            cluster, cluster, ID))
        sma, smb = phot['sma'], phot['smb']
        R_e, width = phot['R_e'], phot['width']
        
        lo = phot['{}_flux'.format(dictionary[cluster]['lo'])]
        hi = phot['{}_flux'.format(dictionary[cluster]['hi'])]
        
        lo_use = phot['{}_use'.format(dictionary[cluster]['lo'])][0]
        hi_use = phot['{}_use'.format(dictionary[cluster]['hi'])][0]
        
        if lo_use and hi_use :
            check = color_from_flux(lo, hi) # calculate colors to check values
            mask = (~np.isnan(check)) & (~np.isinf(check)) # mask based on valid entries
            
            # mask quantities based on valid entries
            color = color_from_flux(lo[mask], hi[mask])
            xs = ((sma - width)*np.sqrt(smb/sma)/R_e)[mask]
            lo_err = phot['{}_err'.format(dictionary[cluster]['lo'])][mask]
            hi_err = phot['{}_err'.format(dictionary[cluster]['hi'])][mask]
            
            if len(color) == 0 :
                value = np.nan
            else :
                if len(color) == 2 :
                    value = float( (color[1] - color[0])/(xs[1] - xs[0]) )
                else :
                    color_err = color_err_from_flux(lo[mask], lo_err, hi[mask], hi_err)
                    popt, pcov = curve_fit(core.linear, xs, color, sigma=color_err)
                    value = popt[0]
            
            gradients.append(value)
        else :
            gradients.append(np.nan)
    
    catalog['color_grad'] = gradients
    
    return catalog

def compute_projected_distances(catalog) :
    
    dictionary = {'a370':{'z_clus':0.375, 'sigma':1170*u.km/u.s, 'R200':2.66*u.Mpc},
                  'a370par':{'z_clus':0.375, 'sigma':1170*u.km/u.s, 'R200':2.66*u.Mpc},
                  'a1063':{'z_clus':0.348, 'sigma':1840*u.km/u.s, 'R200':2.38*u.Mpc},
                  'a1063par':{'z_clus':0.348, 'sigma':1840*u.km/u.s, 'R200':2.38*u.Mpc},
                  'a2744':{'z_clus':0.308, 'sigma':1497*u.km/u.s, 'R200':2.35*u.Mpc},
                  'a2744par':{'z_clus':0.308, 'sigma':1497*u.km/u.s, 'R200':2.35*u.Mpc},
                  'm416':{'z_clus':0.396, 'sigma':955*u.km/u.s, 'R200':1.88*u.Mpc},
                  'm416par':{'z_clus':0.396, 'sigma':955*u.km/u.s, 'R200':1.88*u.Mpc},
                  'm717':{'z_clus':0.545, 'sigma':1660*u.km/u.s, 'R200':2.36*u.Mpc},
                  'm717par':{'z_clus':0.545, 'sigma':1660*u.km/u.s, 'R200':2.36*u.Mpc},
                  'm1149':{'z_clus':0.543, 'sigma':1840*u.km/u.s, 'R200':2.35*u.Mpc},
                  'm1149par':{'z_clus':0.543, 'sigma':1840*u.km/u.s, 'R200':2.35*u.Mpc}}
    
    all_z_clus, all_sigma, all_R200 = [], [], []
    for cluster in catalog['cluster'] :
        all_z_clus.append(dictionary[cluster]['z_clus'])
        all_sigma.append(dictionary[cluster]['sigma'])
        all_R200.append(dictionary[cluster]['R200'])
    
    catalog['z_clus'] = all_z_clus
    catalog['sigma'] = all_sigma
    catalog['R200'] = all_R200
    
    clusters = ['a370', 'a370par',
                'a1063', 'a1063par',
                'a2744', 'a2744par',
                'm416', 'm416par',
                'm717', 'm717par',
                'm1149', 'm1149par']
    coords = [('02h39m52.9s', '-01d34m36.5s'), ('02h39m52.9s', '-01d34m36.5s'),
              ('22h48m44.4s', '-44d31m48.5s'), ('22h48m44.4s', '-44d31m48.5s'),
              ('00h14m21.2s', '-30d23m50.1s'), ('00h14m21.2s', '-30d23m50.1s'),
              ('04h16m08.9s', '-24d04m28.7s'), ('04h16m08.9s', '-24d04m28.7s'),
              ('07h17m34.0s', '+37d44m49.0s'), ('07h17m34.0s', '+37d44m49.0s'),
              ('11h49m36.3s', '+22d23m58.1s'), ('11h49m36.3s', '+22d23m58.1s')]
    
    proj_dists = []
    for cluster, coord in zip(clusters, coords) :
        center = SkyCoord(coord[0], coord[1])
        
        cluster_specific_cat = catalog
        cat = cluster_specific_cat[catalog['cluster'] == cluster] # mask
        
        scales = cosmo.kpc_comoving_per_arcmin(cat['z_clus'])
        
        positions = SkyCoord(cat['ra']*u.deg, cat['dec']*u.deg)
        
        seps = center.separation(positions).arcminute*u.arcmin
        
        proj_dist = (scales*seps).to(u.Mpc)
        proj_dists.append(proj_dist)
    
    catalog['R'] = np.concatenate(proj_dists).ravel()
    
    return catalog

def flux_to_mag(flux) :
    return -2.5*np.log10(flux/(3631*u.Jy))

def phase_space_plot() :
    
    table = Table.read('output/tables/sample_phase-space_with-color-gradients.fits')
    colors = table['color_grad']
    xs = table['R']/table['R200']
    ys = np.absolute(c.c.to(u.km/u.s)*(table['z_best'] -
                                       table['z_clus'])/table['sigma'])
    
    # check the distribution for the color gradients
    # plt.histogram(table['color_grad'], 'color gradients', bins=20, histtype='step')
    
    # change the range for the colorbar
    # vmin, vmax = -0.4, 0.4
    # colors[colors < vmin] = vmin
    # colors[colors > vmax] = vmax
    
    masks = [(table['pop'] == 'Q'), # 617 quiescent galaxies
             (table['pop'] == 'SF'), # 160 star forming galaxies
             (table['env'] == 'cluster'), # 672 cluster galaxies
             (table['env'] == 'field')] # 105 field galaxies
    for mask in masks :
        temp = colors[mask]
        temp = temp[temp > -99]
        # print(len(temp))
        lo, med, hi = np.percentile(temp, [2.5, 50, 97.5])
        # print(med)
        
        # colors[colors < lo] = lo
        # colors[colors > hi] = hi
        '''
        plt.plot_scatter(xs[mask], ys[mask], colors[mask], '', 'o', cmap=cm.bwr,
                         xlabel=r'$R/R_{200}$', ylabel=r'$\Delta v/\sigma$',
                         cbar_label='color gradient', scale='log',
                         xmin=0.001, xmax=10, ymin=0.01, ymax=100,
                         figsizewidth=16, figsizeheight=9)
        '''
    
    xx = [xs[(table['env'] == 'cluster') & (table['pop'] == 'Q')],
          xs[(table['env'] == 'cluster') & (table['pop'] == 'SF')],
          xs[(table['env'] == 'field') & (table['pop'] == 'Q')],
          xs[(table['env'] == 'field') & (table['pop'] == 'SF')]]
    
    yy = [ys[(table['env'] == 'cluster') & (table['pop'] == 'Q')],
          ys[(table['env'] == 'cluster') & (table['pop'] == 'SF')],
          ys[(table['env'] == 'field') & (table['pop'] == 'Q')],
          ys[(table['env'] == 'field') & (table['pop'] == 'SF')]]
    
    plt.plot_simple_multi(xx, yy,
                          ['Cluster QGs', 'Cluster SFGs', 'Field QGs', 'Field SFGs'],
                          ['darkred', 'darkblue', 'r', 'dodgerblue'],
                          ['s', 's', 'o', 'o'], ['', '', '', ''],
                          alphas=[0.3, 0.6, 0.5, 0.5],
                          xlabel=r'$R/R_{200}$', ylabel=r'$\Delta v/\sigma$',
                          figsizewidth=16, figsizeheight=9)
    
    # plt.plot_colorcolor_multi(xs, ys,
    #                           ,
    #                           [cluster_q_len, cluster_sf_len,
    #                            field_q_len, field_sf_len],
    #                           ,
    #                           , [19, 19, 15, 15],
    #                           ,
    #                           [q_xi, sf_xi], [q_yi, sf_yi], [q_z, sf_z],
    #                           version='FUVVJ',
    #                           xlabel=r'V$-$J', ylabel=r'FUV$-$V',
    #                           xmin=0, xmax=2.1, ymin=0, ymax=8.4, loc=2)
    
    
    return

def phase_space_prep() :
    
    catalog = Table.read('output/tables/sample_final.fits')
    
    catalog = compute_projected_distances(catalog)
    
    catalog = compute_color_gradients(catalog)
    
    catalog.write('output/tables/sample_phase-space_with-color-gradients.fits')
    
    return
