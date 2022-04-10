
import numpy as np

# import astropy.constants as c
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u
# from matplotlib import cm
from prospect.sources.constants import cosmo

import plotting as plt

def SFMS() :
    # plot the star forming main sequence
    
    sdss_infile = 'catalogs/SDSS_DR4_from_Gallazzi/galinfo_totlgm_totsfr_totssfr_dr7_v5_2b.fits'
    sdss = Table.read(sdss_infile)
    
    cosmos = Table.read('catalogs/COSMOS/3D-HST_COSMOS.fits')
    
    sdss_mask = (sdss['lmass'] >= 8) & (sdss['lsfr'] >= -6) & (sdss['lssfr'] >= -15)
    sdss = sdss[sdss_mask]
    # sdss = sdss[(sdss['Z'] >= 0.25) & (sdss['Z'] <= 0.6)]
    
    clus = Table.read('output/tables/all_clusters.fits')
    par = Table.read('output/tables/all_parallels.fits')
    
    xx = np.linspace(8, 11.4, 1000)
    # yy = -7.64 + 0.76*xx # fit from Renzini+Peng 2015, defines the ridge line?
    # yy = xx -9.83 - 0.35*(xx - 10) # fit from Salim+ 2007
    # yy = np.log10(0.87) - 0.77*10 + 0.77*xx # fit from Elbaz+ 2007, eq. 5
    yy = xx - 10.0 - 0.1*(xx - 10.0) # fit from Peng+ 2010
    
    xs = [sdss['lmass'], clus['lmass'], par['lmass'], cosmos['lmass'], xx]
    ys = [sdss['lsfr'], clus['lsfr'], par['lsfr'], cosmos['lsfr'], yy]
    
    plt.plot_scatter_multi(xs, ys, [sdss['lssfr'], 'k', 'grey', 'm'],
                           ['SDSS DR7',
                            'HFF cluster ({})'.format(len(clus)),
                            'HFF field ({})'.format(len(par)),
                            'COSMOS ({})'.format(len(cosmos)),
                            'Peng+2010'],
                           ['.', 'o', 's', '^'], [0.4, 1, 1, 1],
                           cbar_label=r'sSFR (${\rm yr}^{-1}$)',
                           xlabel=r'$\log(M/M_{\odot})$',
                           ylabel=r'$\log({\rm SFR}/M_{\odot}~{\rm yr}^{-1})$',
                           xmin=8, xmax=11.4, ymin=-4, ymax=3, loc=0)
    
    return

def cluster_centric_radius() :
    
    clusters = ['a370', 'a1063', 'a2744', 'm416', 'm717', 'm1149']
    coords = [('02h39m52.9s', '-01d34m36.5s'),
              ('22h48m44.4s', '-44d31m48.5s'),
              ('00h14m21.2s', '-30d23m50.1s'),
              ('04h16m08.9s', '-24d04m28.7s'),
              ('07h17m34.0s', '+37d44m49.0s'),
              ('11h49m36.3s', '+22d23m58.1s')]
    
    clus = Table.read('output/tables/all_clusters.fits')
    
    # d2d = []
    d3d = []
    for cluster, coord in zip(clusters, coords) :
        center = SkyCoord(coord[0], coord[1])
        
        new = clus
        new = new[clus['cluster'] == cluster]
        
        scales = cosmo.kpc_comoving_per_arcmin(new['z_spec'])
        
        cat = SkyCoord(new['ra']*u.deg, new['dec']*u.deg)
        
        seps = center.separation(cat).arcminute*u.arcmin
        
        dist = scales*seps
        d3d.append(dist)
    
    d3d = np.concatenate(d3d).ravel()
    
    print(d3d)
    
    # par = Table.read('output/tables/all_parallels.fits')
    
    return

# cluster_centric_radius()

SFMS()
