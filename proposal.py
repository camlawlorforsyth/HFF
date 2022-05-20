
import numpy as np

from astropy.table import Table

import plotting as plt

def SFMS() :
    # plot the star forming main sequence
    
    sdss_infile = 'catalogs/SDSS_DR4_from_Gallazzi/galinfo_totlgm_totsfr_totssfr_dr7_v5_2b.fits'
    sdss = Table.read(sdss_infile)
    
    cosmos = Table.read('catalogs/COSMOS/3D-HST_COSMOS.fits')
    
    sdss_mask = (sdss['lmass'] >= 8) & (sdss['lsfr'] >= -6) & (sdss['lssfr'] >= -15)
    sdss = sdss[sdss_mask]
    # sdss = sdss[(sdss['Z'] >= 0.25) & (sdss['Z'] <= 0.6)]
    
    HFF = Table.read('output/tables/sample_final.fits')
    clus = HFF[HFF['env'] == 'cluster']
    field = HFF[HFF['env'] == 'field']
    
    xx = np.linspace(8, 11.4, 1000)
    # yy = -7.64 + 0.76*xx # fit from Renzini+Peng 2015, defines the ridge line?
    # yy = xx -9.83 - 0.35*(xx - 10) # fit from Salim+ 2007
    # yy = np.log10(0.87) - 0.77*10 + 0.77*xx # fit from Elbaz+ 2007, eq. 5
    yy = xx - 10.0 - 0.1*(xx - 10.0) # fit from Peng+ 2010
    
    xs = [sdss['lmass'], clus['lmass'], field['lmass'], cosmos['lmass'], xx]
    ys = [sdss['lsfr'], clus['lsfr'], field['lsfr'], cosmos['lsfr'], yy]
    
    plt.plot_scatter_multi(xs, ys, [sdss['lssfr'], 'k', 'grey', 'm'],
                           ['SDSS DR7',
                            'HFF cluster ({})'.format(len(clus)),
                            'HFF field ({})'.format(len(field)),
                            'COSMOS ({})'.format(len(cosmos)),
                            'Peng+2010'],
                           ['.', 'o', 's', '^'], [0.4, 1, 1, 1],
                           cbar_label=r'sSFR (${\rm yr}^{-1}$)',
                           xlabel=r'$\log(M/M_{\odot})$',
                           ylabel=r'$\log({\rm SFR}/M_{\odot}~{\rm yr}^{-1})$',
                           xmin=8, xmax=11.4, ymin=-4, ymax=3, loc=0,
                           figsizewidth=14, figsizeheight=10)
    
    return
