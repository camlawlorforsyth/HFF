
import os
import glob
import numpy as np

from astropy.io import fits
from astropy.table import Table, join
import prospect.io.read_results as reader
from prospect.plotting.utils import sample_posterior
from scipy.optimize import curve_fit

from core import open_cutout
import plotting as plt

def prep(save=False) :
    
    info = Table.read('catalogs/SDSS_DR4_from_Gallazzi/gal_info_dr4_v5_1b.fit.gz')
    
    masses = Table.read('catalogs/SDSS_DR4_from_Gallazzi/all_stat_mstar.dat.gz',
                        format='ascii.no_header')
    masses.rename_column('col1', 'PLATEID')
    masses.rename_column('col2', 'MJD')
    masses.rename_column('col3', 'FIBERID')
    masses.rename_column('col4', 'MASS_P2P5')
    masses.rename_column('col5', 'MASS_P16')
    masses.rename_column('col6', 'MASS_P50')
    masses.rename_column('col7', 'MASS_P84')
    masses.rename_column('col8', 'MASS_P97P5')
    del masses['col9']
    del masses['col10']
    
    metallicities = Table.read('catalogs/SDSS_DR4_from_Gallazzi/all_stat_z_log.dat.gz',
                               format='ascii.no_header')
    metallicities.rename_column('col1', 'PLATEID')
    metallicities.rename_column('col2', 'MJD')
    metallicities.rename_column('col3', 'FIBERID')
    metallicities.rename_column('col4', 'METALLICITY_P2P5')
    metallicities.rename_column('col5', 'METALLICITY_P16')
    metallicities.rename_column('col6', 'METALLICITY_P50')
    metallicities.rename_column('col7', 'METALLICITY_P84')
    metallicities.rename_column('col8', 'METALLICITY_P97P5')
    del metallicities['col9']
    del metallicities['col10']
    
    # convert the metallicities into solar values
    # metallicities['METALLICITY_P2P5'] = metallicities['METALLICITY_P2P5'] - np.log10(0.02)
    # metallicities['METALLICITY_P16'] = metallicities['METALLICITY_P16'] - np.log10(0.02)
    # metallicities['METALLICITY_P50'] = metallicities['METALLICITY_P50'] - np.log10(0.02)
    # metallicities['METALLICITY_P84'] = metallicities['METALLICITY_P84'] - np.log10(0.02)
    # metallicities['METALLICITY_P97P5'] = metallicities['METALLICITY_P97P5'] - np.log10(0.02)
    
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
    del final['Z']
    del final['Z_ERR']
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
    # final = final[final['METALLICITY_P50'] > -2]
    # final = final[final['SN_MEDIAN'] > 20]
    
    if save :
        final.write('output/MZR_alt.fits')
    
    return

def tanh(xx, aa, bb, cc, dd) :
    return aa + bb*np.tanh((xx - cc)/dd)

def determine_medians(infile='output/MZR_with_SFR.fits') :
    
    table = Table.read(infile) # MZR uses the metallicities from Gallazzi+2005
        # based on DR4, but updated total mass estimates from fits to DR7
        # photometry
    
    # mask the data based on reasonable stellar masses
    table = table[table['MASS_P50'] > 8]
    
    # mask the table based on reasonable metallicities
    table = table[table['METALLICITY_P50'] - np.log10(0.02) > -2]
    
    # mask the table based on SFR
    table = table[table['SFR_P50'] > -99]
    
    # mask the table based on sSFR
    table = table[table['SFR_P50'] - table['MASS_P50'] < -11]
    
    # now select only galaxies with sufficient signal to noise in their pixels
    # table = table[table['SN_MEDIAN'] > 3]
    
    # populate values from the table
    mass = table['MASS_P50']
    metal = table['METALLICITY_P50'] - np.log10(0.02) # convert to solar values
    
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
    
    popt, pcov = curve_fit(tanh, centers, medians,
                           p0=[-0.452, 0.572, 9.66, 1.04]) # initial parameters
        # come from Panter+2008, MNRAS, 391, 1117
    
    # y_model = tanh(centers, *popt)
    
    xs = np.linspace(8, 12, 1000)
    ys = tanh(xs, *popt)
    
    if not os.path.isfile('output/HFF_integrated_logM_and_logZ.fits') :
        all_masses_and_metallicities()
        HFF = {}
        HFF['logM'], HFF['avglogZ'] = [0], [0]
    else :
        HFF = Table.read('output/HFF_integrated_logM_and_logZ.fits')
    
    plt.histogram_2d(mass, metal, [20, 20], centers,
                     [tinies, los, medians, his, larges], xs, ys,
                     HFF['logM'], HFF['avglogZ'], ['fit', 'HFF'],
                     xlabel=r'$\log(M/M_{\odot})$', ylabel=r'$\log(Z/Z_{\odot})$',
                     xmin=8, xmax=12, ymin=-1.6, ymax=0.4)
    
    return

def all_masses_and_metallicities() :
    
    clusters = ['a370', 'a1063', 'a2744', 'm416', 'm717', 'm1149']
    
    all_clus, all_IDs, all_masses, all_metallicities = [], [], [], []
    for cluster in clusters :
        (clusterID, IDs, masses,
         metallicities) = luminosity_weighted_average_metallicity(cluster)
        
        for clus in clusterID :
            all_clus.append(clus)
        
        for ID in IDs :
            all_IDs.append(ID)
        
        for mass in masses :
            all_masses.append(mass)
        
        for metallicity in metallicities :
            all_metallicities.append(metallicity)
    
    table = Table([all_clus, all_IDs, all_masses, all_metallicities],
                  names=('cluster', 'ID', 'logM', 'avglogZ'))
    table.write('output/HFF_integrated_logM_and_logZ.fits')
    
    return 

def luminosity_weighted_average_metallicity(cluster) :
    
    paths = '{}/spatial_result_images/{}_ID_*_logZ.fits'.format(cluster, cluster)
    images = glob.glob(paths)
    
    clusters, IDs, masses, metallicities = [], [], [], []
    for image in images :
        file = image.replace(os.sep, '/') # compatibility for Windows
        ID = int(file.split('_')[4]) # the galaxy ID
        
        with fits.open(file) as hdu :
            metal_image = hdu[0].data
        
        # get F160W cutout for luminosity-weighting
        flux = open_cutout('{}/cutouts/{}_ID_{}_f160w.fits'.format(cluster, cluster, ID),
                           simple=True)
        segmap = open_cutout('{}/cutouts/{}_ID_{}_segmap.fits'.format(cluster, cluster, ID),
                             simple=True)
        
        new_sci = flux.copy()
        new_sci[(segmap >= 0) & (segmap != ID)] = 0 # mask out the background
            # and pixels associated with other galaxies
        
        flux_fraction = new_sci/np.sum(new_sci)
        
        average_metallicity = np.nansum(flux_fraction*metal_image)
        metallicities.append(average_metallicity)
        
        with fits.open(file.replace('Z', 'M')) as hdu :
            mass_image = hdu[0].data
        
        total_mass = np.log10(np.nansum(mass_image))
        masses.append(total_mass)
        
        clusters.append(cluster)
        IDs.append(ID)
    
    return clusters, IDs, masses, metallicities

def construct_images(cluster, disallowed_IDs, selection) :
    
    if selection == 'mass' :
        index, string = 1, 'logM'
    if selection == 'metallicity' :
        index, string = 2, 'logZ'
    
    os.makedirs('{}/spatial_result_images'.format(cluster), exist_ok=True)
        # ensure the output directory for the metallicity images is available
    
    phot_paths = '{}/photometry/{}_ID_*_photometry.fits'.format(cluster,
                                                                cluster)
    phots = glob.glob(phot_paths)
    
    for file in phots :
        file = file.replace(os.sep, '/') # compatibility for Windows
        ID = int(file.split('_')[2]) # the galaxy ID
        
        if (ID not in disallowed_IDs) :
            table = Table.read(file)
            
            bin_data = np.load('{}/bins/{}_ID_{}_annuli.npz'.format(cluster, cluster, ID))
            
            bins_image = bin_data['image']
            
            bins = table['bin'] # get a list of bin values
            for binNum in bins : # loop over all the bins in the table
                infile = '{}/h5/{}_ID_{}_bin_{}.h5'.format(cluster, cluster,
                                                           ID, binNum)
                
                result, obs, _ = reader.results_from(infile, dangerous=True)
                samples = sample_posterior(result['chain'],
                                           weights=result['weights'])
                
                median = np.percentile(samples[:, index], 50)
                
                bins_image[bins_image == binNum] = median
            
            hdu = fits.PrimaryHDU(bins_image)
            outfile = '{}/spatial_result_images/{}_ID_{}_{}.fits'.format(cluster, cluster,
                                                                         ID, string)
            hdu.writeto(outfile)
    
    return

# construct_images('a370', [3826, 4094, 4966, 5069], 'mass')
# construct_images('a370', [3826, 4094, 4966, 5069], 'metallicity')
# construct_images('a1063', [], 'mass')
# construct_images('a1063', [], 'metallicity')
# construct_images('a2744', [6727], 'mass')
# construct_images('a2744', [6727], 'metallicity')
# construct_images('m416', [3594, 3638, 5413], 'mass')
# construct_images('m416', [3594, 3638, 5413], 'metallicity')
# construct_images('m717', [1358, 2774, 3692, 3731, 4784, 4814, 5216], 'mass')
# construct_images('m717', [1358, 2774, 3692, 3731, 4784, 4814, 5216], 'metallicity')
# construct_images('m1149', [2436, 2978, 3751, 5808], 'mass')
# construct_images('m1149', [2436, 2978, 3751, 5808], 'metallicity')

# determine_medians()
