
import numpy as np

from astropy.table import Table, join

import plotting as plt

def prep(save=False) :
    
    info = Table.read('catalogs/SDSS_DR4_from_Gallazzi/gal_info_dr4_v5_1b.fit.gz')
    
    masses = Table.read('catalogs/SDSS_DR4_from_Gallazzi/all_stat_mstar.dat.gz',
                        format='ascii.no_header')
    masses.rename_column('col1', 'PLATEID')
    masses.rename_column('col2', 'MJD')
    masses.rename_column('col3', 'FIBERID')
    masses.rename_column('col4', 'M_2p5')
    masses.rename_column('col5', 'M_16')
    masses.rename_column('col6', 'M_50')
    masses.rename_column('col7', 'M_84')
    masses.rename_column('col8', 'M_97p5')
    del masses['col9']
    del masses['col10']
    
    metallicities = Table.read('catalogs/SDSS_DR4_from_Gallazzi/all_stat_z_log.dat.gz',
                               format='ascii.no_header')
    metallicities.rename_column('col1', 'PLATEID')
    metallicities.rename_column('col2', 'MJD')
    metallicities.rename_column('col3', 'FIBERID')
    metallicities.rename_column('col4', 'Z_2p5')
    metallicities.rename_column('col5', 'Z_16')
    metallicities.rename_column('col6', 'Z_50')
    metallicities.rename_column('col7', 'Z_84')
    metallicities.rename_column('col8', 'Z_97p5')
    del metallicities['col9']
    del metallicities['col10']
    
    # convert the metallicities into solar values
    metallicities['Z_2p5'] = metallicities['Z_2p5'] - np.log10(0.02)
    metallicities['Z_16'] = metallicities['Z_16'] - np.log10(0.02)
    metallicities['Z_50'] = metallicities['Z_50'] - np.log10(0.02)
    metallicities['Z_84'] = metallicities['Z_84'] - np.log10(0.02)
    metallicities['Z_97p5'] = metallicities['Z_97p5'] - np.log10(0.02)
    
    derived = join(masses, metallicities, keys=['PLATEID', 'MJD', 'FIBERID'])
    final = join(info, derived, keys=['PLATEID', 'MJD', 'FIBERID'])
    del final['PHOTOID']
    del final['RA']
    del final['DEC']
    del final['PLUG_MAG']
    del final['PRIMTARGET']
    del final['SECTARGET']
    del final['TARGETTYPE']
    del final['SPECTROTYPE']
    del final['SUBCLASS']
    del final['Z_WARNING']
    del final['V_DISP']
    del final['V_DISP_ERR']
    del final['E_BV_SFD']
    del final['ZTWEAK']
    del final['ZTWEAK_ERR']
    del final['SPECTRO_MAG']
    del final['KCOR_MAG']
    del final['KCOR_MODEL_MAG']
    del final['RELEASE']
    final = final[final['Z_50'] > -2]
    final = final[final['SN_MEDIAN'] > 20]
    
    if save :
        final.write('output/MZR.fits')
    
    return

def determine_medians(infile='output/MZR.fits') :
    
    table = Table.read(infile)
    table = table[table['M_50'] > 8]
    # table.sort('m50')
    
    mass = table['M_50']
    metal = table['Z_50']
    
    edges = np.arange(8, 12.2, 0.2)
    
    los, medians, his, centers, counts = [], [], [], [], []
    for first, second in zip(edges, edges[1:]) :
        mass_mask = (mass >= first) & (mass < second)
        lo, median, hi = np.percentile(metal[mass_mask], [16, 50, 84])
        
        los.append(lo)
        medians.append(median)
        his.append(hi)
        
        centers.append(np.mean([first, second]))
        
        counts.append(len(metal[mass_mask]))
    
    counts = np.log10(counts)
    plt.histogram_2d(mass, metal, [20, 20], centers, [los, medians, his, counts],
                     xlabel=r'$\log(M/M_{\odot})$', ylabel=r'$\log(Z/Z_{\odot})$',
                     xmin=8, xmax=12, ymin=-1.6, ymax=0.4)
    
    return

determine_medians()

