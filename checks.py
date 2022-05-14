
import os
import glob
import numpy as np

from astropy.table import Table, vstack, join
from scipy.stats import kde
from sedpy.observate import load_filters
import statmorph

import core
import plotting as plt

def check_all_bins(plot=True) :
    
    table = Table.read('output/number_of_annuli.csv')
    bins = table['bins']
    nbCG_QG_bins = bins[table['ID'] < 20000]
    bCG_bins = bins[table['ID'] > 20000]
    
    total_annuli = np.sum(nbCG_QG_bins) + np.sum(bCG_bins)
    print(total_annuli)
    
    if plot :
        numBins1 = int(np.ceil(np.sqrt(len(nbCG_QG_bins))))
        # numBins2 = int(np.ceil(np.sqrt(len(bCG_bins))))
        
        plt.histogram_multi([nbCG_QG_bins], 'Number of Annuli',
                            bins=[numBins1], colors=['k'],
                            labels=['non-bCG QGs'], styles=['-'])
    
    return

def check_asymmetries() :
    
    mass_table = Table.read('output/tables/nbCGs_with-Shipley-mass.fits')
    asymm_table = Table.read('output/tables/nbCGs_asymmetries.fits')
    
    table = join(mass_table, asymm_table, keys=['cluster', 'ID'])
    # table.sort('asymmetry')
    # table.pprint(max_lines=-1, max_width=-1)
    
    plt.plot_simple_multi([table['logM']], [table['asymmetry']],
                          [''], ['k'], ['o'], [''],
                          xlabel=r'$\log(M/M_{\odot})$', ylabel='Asymmetry',
                          scale='linear', xmin=8, xmax=11.5, ymin=-1, ymax=1)
    
    return

def check_distributions() :
    
    all_clusters = concatenate_all()
    all_clusters = all_clusters[all_clusters['pop'] == 'Q']
    
    bins = int(np.round(np.sqrt(len(all_clusters))))
    
    redshift = all_clusters['z_spec']
    mstar = all_clusters['lmass']
    r_e = all_clusters['flux_radius']
    
    plt.histogram(redshift, r'$z_{\rm spec}$', bins=bins, histtype='step',
                  vlines=[0.308, 0.375, 0.396, 0.545, 0.543, 0.348],
                  colors=['purple', 'g', 'gold', 'r', 'orange', 'b'],
                  labels=['A2744', 'A370', 'M416', 'M717', 'M1149', 'AS1063'])
    
    plt.histogram(mstar, r'$\log(M_{*}/M_{\odot})$', bins=bins,
                  histtype='step')
    
    plt.histogram(np.log10(r_e), r'$\log(R_{\rm e}/{\rm pix})$',
                  bins=bins, histtype='step')
    
    return

def check_fit_progress() :
    
    clusters = ['a370', 'a1063', 'a2744', 'm416', 'm717', 'm1149']
    
    for cluster in clusters :
        phot_paths = '{}/photometry/{}_ID_*_photometry.fits'.format(cluster,
                                                                    cluster)
        phots = glob.glob(phot_paths)
        
        for file in phots :
            file = file.replace(os.sep, '/') # compatibility for Windows
            ID = int(file.split('_')[2]) # the galaxy ID to fit the bins for
            
            table = Table.read(file)
            bins = table['bin'] # get a list of bin values
            for binNum in bins : # loop over all the bins in the table
                h5_file = '{}/h5/{}_ID_{}_bin_{}.h5'.format(cluster, cluster,
                                                            ID, binNum)
                
                if not os.path.isfile(h5_file) :
                    print('Missing {}'.format(h5_file))
    
    return

def check_for_missing_issues() :
    
    from astropy.table import Table, setdiff
    
    for cluster in ['a370', 'a1063', 'a2744', 'm416', 'm717', 'm1149',
                    'a370par', 'a1063par', 'a2744par', 'm416par', 'm717par',
                    'm1149par'] :
        
        old = Table.read('{}/{}_issues.csv'.format(cluster, cluster))
        old = Table([old['id'], old['pop']], names=('id', 'pop'))
        
        new = Table.read('{}/{}_sample.fits'.format(cluster, cluster))
        new = Table([new['id'], new['pop']], names=('id', 'pop'))
        
        unique = setdiff(new, old, keys=['id', 'pop'])
        print(cluster, len(unique))
        if len(unique > 0) :
            unique.pprint(max_lines=-1)
        print()
    
    return

def check_radial_extents() :
    
    HFF = Table.read('output/tables/nbCGs.fits')
    
    extents = []
    for cluster, ID in zip(HFF['cluster'], HFF['ID']) :
        table = Table.read('{}/photometry/{}_ID_{}_photometry.fits'.format(
            cluster, cluster, ID))
        
        sma, smb = table['sma'], table['smb']
        R_e, width = table['R_e'], table['width']
        
        xs = (sma - width)*np.sqrt(smb/sma)/R_e
        # xerrs_hi = width*np.sqrt(smb/sma)/R_e
        
        extents.append(xs[-1])
    
    HFF['max_Re'] = extents
    HFF.write('output/tables/nbCGs_radial-extent.fits')
    
    return

def check_SNR(filt) :
    
    clusters = ['a370', 'a1063', 'a2744', 'm416', 'm717', 'm1149']
    
    SNRs_lo, SNRs_med, SNRs_hi = [], [], []
    for cluster in clusters :
        
        inDir = '{}/cutouts'.format(cluster)
        
        cluster_table = Table.read('{}/{}_sample.fits'.format(
            cluster, cluster))
        
        mask = (cluster_table['id'] < 20000) & (cluster_table['pop'] == 'Q')
        cluster_table = cluster_table[mask]
        
        lo_mask = cluster_table['lmass'] < 8.78
        med_mask = ((cluster_table['lmass'] >= 9.52) &
                    (cluster_table['lmass'] <= 9.76))
        hi_mask = cluster_table['lmass'] >= 10.5
        
        lo_table = cluster_table.copy()
        med_table = cluster_table.copy()
        hi_table = cluster_table.copy()
        
        lo_table = lo_table[lo_mask]
        med_table = med_table[med_mask]
        hi_table = hi_table[hi_mask]
        
        for ID in lo_table['id'] :
            sci_file = '{}/{}_ID_{}_{}.fits'.format(inDir, cluster, ID, filt)
            noise_file = '{}/{}_ID_{}_{}_noise.fits'.format(inDir, cluster,
                                                            ID, filt)
            segmap_file = '{}/{}_ID_{}_segmap.fits'.format(inDir, cluster, ID)
            
            try :
                sci = core.open_cutout(sci_file, simple=True)
                noise = core.open_cutout(noise_file, simple=True)
                segmap = core.open_cutout(segmap_file, simple=True)
                
                SNR = sci/noise
                SNR[(segmap != ID)] = np.nan
                SNR_flat = SNR.flatten()
                
                SNR_flat = SNR_flat[~np.isnan(SNR_flat)] # remove nans
                
                for val in SNR_flat :
                    SNRs_lo.append(val)
            except :
                pass
        
        for ID in med_table['id'] :
            sci_file = '{}/{}_ID_{}_{}.fits'.format(inDir, cluster, ID, filt)
            noise_file = '{}/{}_ID_{}_{}_noise.fits'.format(inDir, cluster,
                                                            ID, filt)
            segmap_file = '{}/{}_ID_{}_segmap.fits'.format(inDir, cluster, ID)
            
            try :
                sci = core.open_cutout(sci_file, simple=True)
                noise = core.open_cutout(noise_file, simple=True)
                segmap = core.open_cutout(segmap_file, simple=True)
                
                SNR = sci/noise
                SNR[(segmap != ID)] = np.nan
                SNR_flat = SNR.flatten()
                
                SNR_flat = SNR_flat[~np.isnan(SNR_flat)]
                
                for val in SNR_flat :
                    SNRs_med.append(val)
            except :
                pass
        
        for ID in hi_table['id'] :
            sci_file = '{}/{}_ID_{}_{}.fits'.format(inDir, cluster, ID, filt)
            noise_file = '{}/{}_ID_{}_{}_noise.fits'.format(inDir, cluster,
                                                            ID, filt)
            segmap_file = '{}/{}_ID_{}_segmap.fits'.format(inDir, cluster, ID)
            
            try :
                sci = core.open_cutout(sci_file, simple=True)
                noise = core.open_cutout(noise_file, simple=True)
                segmap = core.open_cutout(segmap_file, simple=True)
                
                SNR = sci/noise
                SNR[(segmap != ID)] = np.nan
                SNR_flat = SNR.flatten()
                
                SNR_flat = SNR_flat[~np.isnan(SNR_flat)]
                
                for val in SNR_flat :
                    SNRs_hi.append(val)
            except :
                pass
    
    lo_label = r'$\log(M_{*}/M_\odot) < 8.78$'
    med_label = r'$9.52 \leq \log(M_{*}/M_\odot) \leq 9.76$'
    hi_label = r'$\log(M_{*}/M_\odot) \geq 10.5$'
    
    plt.histogram_multi([SNRs_lo, SNRs_med, SNRs_hi],
                        '{} SNR'.format(filt),
                        bins=[70, 70, 70],
                        log=True,
                        histtype='step',
                        colors=['k', 'g', 'm'],
                        labels=[lo_label, med_label, hi_label],
                        styles=[':', '--', '-'],
                        xmin=0.1, xmax=1300, ymin=1, ymax=2e5)
    
    return

def check_total_flux(filt) :
    
    all_clusters = concatenate_all()
    all_clusters = all_clusters[(all_clusters['pop'] == 'Q') &
                                (all_clusters['id'] < 20000)]
    all_clusters.remove_row(92) # A1063 ID 4746 - no photometry file
    all_clusters.remove_row(109) # A1063 ID 5638 - no photometry file
    
    IDs = all_clusters['id']
    clusters = all_clusters['cluster']
    cat_filt = filt.upper()
    f_cat = all_clusters['f_{}'.format(cat_filt)]
    e_cat = all_clusters['e_{}'.format(cat_filt)]
    
    f_calc, e_calc = [], []
    for i in range(len(all_clusters)) :
        phot_path = '{}/photometry/{}_ID_{}_photometry.fits'.format(
                clusters[i], clusters[i], IDs[i])
        table = Table.read(phot_path)
        flux_col = table['{}_flux'.format(filt)]
        flux = np.sum(flux_col)
        err_col = table['{}_err'.format(filt)]
        err = np.sqrt(np.sum(np.square(err_col)))
        f_calc.append(flux)
        e_calc.append(err)
    
    m_cat = -2.5*np.log10(f_cat) + 25
    m_calc = -2.5*np.log10(f_calc) + 8.9
    
    delta_m_cat = 2.5/np.log(10)*np.abs(e_cat/f_cat)
    delta_m_calc = 2.5/np.log(10)*np.abs(np.array(e_calc)/np.array(f_calc))
    
    plt.plot_simple(m_cat, m_calc, delta_m_calc, xerr=delta_m_cat,
                    xlabel=r'$m_{\rm AB_{\rm cat}}$',
                    ylabel=r'$m_{\rm AB_{\rm calc}}$',
                    xmin=25, xmax=17.8, ymin=26, ymax=17.8)
    
    return

def concatenate_all(both=False, save=False) :
    
    clusters = ['a370', 'a1063', 'a2744', 'm416', 'm717', 'm1149']
    all_clusters_list = []
    
    for cluster in clusters :
        # infile = '{}/{}_sample.fits'.format(cluster, cluster)
        infile = '{}/{}_non-bCG_QGs.fits'.format(cluster, cluster)
        cluster_table = Table.read(infile)
        cluster_table['cluster'] = cluster
        all_clusters_list.append(cluster_table)
    all_clusters = vstack(all_clusters_list)
    if save :
        all_clusters.write('all_clusters.fits')
    
    if both :
        all_parallels_list = []
        for cluster in clusters :
            parallel_table = Table.read(
                '{}par/{}par_sample.fits'.format(cluster, cluster))
            parallel_table['field'] = '{}par'.format(cluster)
            all_parallels_list.append(parallel_table)
        all_parallels = vstack(all_parallels_list)
        if save :
            all_parallels.write('all_parallels.fits')
        return all_clusters, all_parallels
    else :
        return all_clusters

def load_FUVVJ_contours(infile='output/FUVVJ_contours.npz') :
    
    with np.load(infile) as data :
        q_xi, q_yi, q_z = data['q_xi'], data['q_yi'], data['q_z']
        sf_xi, sf_yi, sf_z = data['sf_xi'], data['sf_yi'], data['sf_z']
    
    return q_xi, q_yi, q_z, sf_xi, sf_yi, sf_z

def plot_all_FUVVJ(plot=True, save_contours=False) :
    
    all_clusters, all_parallels = concatenate_all(both=True)
    
    cluster_q = all_clusters[all_clusters['pop'] == 'Q']
    cluster_sf = all_clusters[all_clusters['pop'] == 'SF']
    
    field_q = all_parallels[all_parallels['pop'] == 'Q']
    field_sf = all_parallels[all_parallels['pop'] == 'SF']
    
    cluster_q_len, cluster_sf_len = len(cluster_q), len(cluster_sf)
    field_q_len, field_sf_len = len(field_q), len(field_sf)
    
    cluster_q_x = cluster_q['M_AB_V'] - cluster_q['M_AB_J']
    cluster_q_y = cluster_q['M_AB_FUV'] - cluster_q['M_AB_V']
    cluster_sf_x = cluster_sf['M_AB_V'] - cluster_sf['M_AB_J']
    cluster_sf_y = cluster_sf['M_AB_FUV'] - cluster_sf['M_AB_V']
    
    field_q_x = field_q['M_AB_V'] - field_q['M_AB_J']
    field_q_y = field_q['M_AB_FUV'] - field_q['M_AB_V']
    field_sf_x = field_sf['M_AB_V'] - field_sf['M_AB_J']
    field_sf_y = field_sf['M_AB_FUV'] - field_sf['M_AB_V']
    
    xs = [cluster_q_x, cluster_sf_x, field_q_x, field_sf_x]
    ys = [cluster_q_y, cluster_sf_y, field_q_y, field_sf_y]
    
    q_x = list(cluster_q_x) + list(field_q_x)
    q_y = list(cluster_q_y) + list(field_q_y)
    sf_x = list(cluster_sf_x) + list(field_sf_x)
    sf_y = list(cluster_sf_y) + list(field_sf_y)
    
    # use Gaussian kernel density estimate
    nbins = 100
    q_data = np.vstack([q_x, q_y])
    q_k = kde.gaussian_kde(q_data)
    q_xi, q_yi = np.mgrid[min(q_x):max(q_x):nbins*1j,
                          min(q_y):max(q_y):nbins*1j]
    q_zi = q_k(np.vstack([q_xi.flatten(), q_yi.flatten()]))
    q_z = q_zi.reshape(q_xi.shape)
    
    sf_data = np.vstack([sf_x, sf_y])
    sf_k = kde.gaussian_kde(sf_data)
    sf_xi, sf_yi = np.mgrid[min(sf_x):max(sf_x):nbins*1j,
                            min(sf_y):max(sf_y):nbins*1j]
    sf_zi = sf_k(np.vstack([sf_xi.flatten(), sf_yi.flatten()]))
    sf_z = sf_zi.reshape(sf_xi.shape)
    
    if plot :
        plt.plot_colorcolor_multi(xs, ys,
                                  ['Cluster QGs', 'Cluster SFGs',
                                   'Field QGs', 'Field SFGs'],
                                  [cluster_q_len, cluster_sf_len,
                                   field_q_len, field_sf_len],
                                  ['darkred', 'darkblue', 'r', 'dodgerblue'],
                                  ['s', 's', 'o', 'o'], [19, 19, 15, 15],
                                  [0.6, 0.6, 0.5, 0.5],
                                  [q_xi, sf_xi], [q_yi, sf_yi], [q_z, sf_z],
                                  version='FUVVJ',
                                  xlabel=r'V$-$J', ylabel=r'FUV$-$V',
                                  xmin=0, xmax=2.1, ymin=0, ymax=8.4, loc=2)
    
    if save_contours :
        np.savez('FUVVJ_contours.npz', q_xi=q_xi, q_yi=q_yi, q_z=q_z,
                 sf_xi=sf_xi, sf_yi=sf_yi, sf_z=sf_z)
    
    return

def plot_all_UVJ() :
    
    all_clusters, all_parallels = concatenate_all(both=True)
    
    cluster_q = all_clusters[all_clusters['pop'] == 'Q']
    cluster_sf = all_clusters[all_clusters['pop'] == 'SF']
    
    field_q = all_parallels[all_parallels['pop'] == 'Q']
    field_sf = all_parallels[all_parallels['pop'] == 'SF']
    
    cluster_q_len, cluster_sf_len = len(cluster_q), len(cluster_sf)
    field_q_len, field_sf_len = len(field_q), len(field_sf)
    
    cluster_q_x = cluster_q['M_AB_V'] - cluster_q['M_AB_J']
    cluster_q_y = cluster_q['M_AB_U'] - cluster_q['M_AB_V']
    cluster_sf_x = cluster_sf['M_AB_V'] - cluster_sf['M_AB_J']
    cluster_sf_y = cluster_sf['M_AB_U'] - cluster_sf['M_AB_V']
    
    field_q_x = field_q['M_AB_V'] - field_q['M_AB_J']
    field_q_y = field_q['M_AB_U'] - field_q['M_AB_V']
    field_sf_x = field_sf['M_AB_V'] - field_sf['M_AB_J']
    field_sf_y = field_sf['M_AB_U'] - field_sf['M_AB_V']
    
    xs = [cluster_q_x, cluster_sf_x, field_q_x, field_sf_x]
    ys = [cluster_q_y, cluster_sf_y, field_q_y, field_sf_y]
    
    q_x = list(cluster_q_x) + list(field_q_x)
    q_y = list(cluster_q_y) + list(field_q_y)
    sf_x = list(cluster_sf_x) + list(field_sf_x)
    sf_y = list(cluster_sf_y) + list(field_sf_y)
    
    # use Gaussian kernel density estimate
    nbins = 100
    q_data = np.vstack([q_x, q_y])
    q_k = kde.gaussian_kde(q_data)
    q_xi, q_yi = np.mgrid[min(q_x):max(q_x):nbins*1j,
                          min(q_y):max(q_y):nbins*1j]
    q_zi = q_k(np.vstack([q_xi.flatten(), q_yi.flatten()]))
    q_z = q_zi.reshape(q_xi.shape)
    
    sf_data = np.vstack([sf_x, sf_y])
    sf_k = kde.gaussian_kde(sf_data)
    sf_xi, sf_yi = np.mgrid[min(sf_x):max(sf_x):nbins*1j,
                            min(sf_y):max(sf_y):nbins*1j]
    sf_zi = sf_k(np.vstack([sf_xi.flatten(), sf_yi.flatten()]))
    sf_z = sf_zi.reshape(sf_xi.shape)
    
    plt.plot_colorcolor_multi(xs, ys,
                              ['Cluster QGs', 'Cluster SFGs',
                               'Field QGs', 'Field SFGs'],
                              [cluster_q_len, cluster_sf_len,
                               field_q_len, field_sf_len],
                              ['darkred', 'darkblue', 'r', 'dodgerblue'],
                              ['s', 's', 'o', 'o'], [19, 19, 15, 15],
                              [0.6, 0.6, 0.5, 0.5],
                              [q_xi, sf_xi], [q_yi, sf_yi], [q_z, sf_z],
                              version='UVJ',
                              xlabel=r'V$-$J', ylabel=r'U$-$V',
                              xmin=0, xmax=2.1, ymin=0, ymax=2.4, loc=4)
    
    return

def plot_bins(cluster, ID) :
    
    from matplotlib import cm
    
    bins_image = core.open_image(cluster, ID, 'bins')
    
    plt.display_image_simple(bins_image, cmap=cm.gist_rainbow, norm=None)
    
    return

def plot_parallel_objects_all() :
    
    all_clusters, all_parallels = concatenate_all(both=True)
    
    xs = [all_parallels['lmass']]
    ys = [all_parallels['z']]
    
    xx = list(all_parallels['lmass'])
    yy = list(all_parallels['z'])
    
    nbins = 100
    data = np.vstack([xx, yy])
    kk = kde.gaussian_kde(data)
    xi, yi = np.mgrid[min(xx):max(xx):nbins*1j,
                      min(yy):max(yy):nbins*1j]
    zi = kk(np.vstack([xi.flatten(), yi.flatten()]))
    zz = zi.reshape(xi.shape)
    
    plt.plot_colorcolor_multi(xs, ys, ['Parallel Finals'],
                              [len(all_parallels['lmass'])], ['cyan'], ['o'],
                              [19], [0.6], [xi], [yi], [zz], version='none',
                              xlabel=r'$\log(M_{*}/M_{\odot})$',
                              ylabel=r'$z$',
                              xmin=2, xmax=14, ymin=0, ymax=0.9)
    
    # this isn't working correctly:
    # plt.plot_objects(xs, ys, 0.3855, ['final'], ['o'], ['cyan'],
    #                  redshift_tol_lo=0.1275, redshift_tol_hi=0.2095,
    #                  xlabel=r'$\log(M_{*}/M_{\odot})$',
    #                  ylabel=r'$z$',
    #                  xmin=2, xmax=14, ymin=0, ymax=0.9)
    
    return

def save_all_bins() :
    
    clusters = ['a370', 'a1063', 'a2744', 'm416', 'm717', 'm1149']
    
    all_clusters, all_IDs, all_lengths = [], [], []
    for cluster in clusters :
        phot_paths = '{}/photometry/{}_ID_*_photometry.fits'.format(cluster,
                                                                    cluster)
        phot_file_list = glob.glob(phot_paths)
        IDs = []
        for file in phot_file_list :
            file = file.replace(os.sep, '/')
            IDs.append(file.split('_')[2])
        
        IDs = np.array(IDs)
        IDs = IDs.astype(int)
        
        for i in range(len(IDs)) :
            table = Table.read(
                '{}/photometry/{}_ID_{}_photometry.fits'.format(cluster,
                                                                cluster,
                                                                IDs[i]))
            length = len(table)
            
            all_clusters.append(cluster)
            all_IDs.append(IDs[i])
            all_lengths.append(length)
    
    master_table = Table([all_clusters, all_IDs, all_lengths],
                         names=('cluster', 'ID', 'bins'))
    master_table.write('output/number_of_annuli.csv')
    
    return

def save_asymmetries(filt='f160w') :
    
    HFF = Table.read('output/tables/nbCGs.fits')
    
    asymmetries = []
    for cluster, ID in zip(HFF['cluster'], HFF['ID']) :
        sci = core.open_cutout('{}/cutouts/{}_ID_{}_{}.fits'.format(
            cluster, cluster, ID, filt), simple=True)
        
        err = core.open_cutout('{}/cutouts/{}_ID_{}_{}_noise.fits'.format(
            cluster, cluster, ID, filt), simple=True)
        
        segmap = core.open_cutout('{}/cutouts/{}_ID_{}_segmap.fits'.format(
            cluster, cluster, ID), simple=True)
        
        morph = statmorph.source_morphology(sci, segmap, weightmap=err,
                                            label=ID)
        asymmetries.append(morph.asymmetry)
    
    HFF['asymmetry'] = asymmetries
    HFF.write('output/tables/nbCGs_asymmetries_NEW.fits')
    
    return

def save_bkgshists(cluster, filters, population) :
    
    num_filts = len(filters)
    if num_filts == 9 :
        nrows, ncols = 3, 3
    elif num_filts == 12 :
        nrows, ncols = 3, 4
    elif num_filts == 16 :
        nrows, ncols = 4, 4
    else :
        nrows, ncols = 4, 5
    
    outDir = '{}/images_bkg_dists'.format(cluster)
    os.makedirs(outDir, exist_ok=True) # ensure the output direc. is available
    
    # now loop here over all the IDs that are available
    # open the table of all the objects 
    clustable = Table.read('{}/{}_sample.fits'.format(cluster, cluster))
    
    # mask the table based on the population of interest
    clustable = clustable[clustable['pop'] == population]
    
    IDs = list(clustable['id'])
    
    for i in range(len(IDs)) :
        outfile = '{}/images_bkg_dists/{}_ID_{}.png'.format(cluster, cluster,
                                                            IDs[i])
        
        segPath = '{}/cutouts/{}_ID_{}_segmap.fits'.format(cluster, cluster,
                                                           IDs[i])
        segMap = core.open_cutout(segPath, simple=True)
        
        cutout_segmapped_data = []
        medians = []
        bins = []
        for filt in filters :
            infile = '{}/cutouts/{}_ID_{}_{}.fits'.format(cluster, cluster,
                                                          IDs[i], filt)
            data = core.open_cutout(infile, simple=True)
            
            if IDs[i] > 20000:
                mask = (segMap > 0)
            else :
                mask = (segMap > 0) & (segMap != IDs[i])
            segmapped_data = data.copy()
            segmapped_data[mask] = np.nan
            
            segmapped_data = segmapped_data.flatten()
            non_nan_data = segmapped_data[~np.isnan(segmapped_data)]
            
            median = np.median(non_nan_data)
            num_bins = int(np.ceil(np.cbrt(len(non_nan_data)))) # Rice rule
                # see https://en.wikipedia.org/wiki/Histogram#Rice_Rule
            
            cutout_segmapped_data.append(non_nan_data)
            medians.append(median)
            bins.append(num_bins)
        
        plt.display_hists(cutout_segmapped_data, nrows, ncols, filters,
                          medians, bins, outfile, save=True)
    
    return

def save_pngs(cluster, filters) :
    
    num_filts = len(filters)
    if num_filts == 7 :
        nrows, ncols = 3, 3
    elif num_filts == 9 :
        nrows, ncols = 3, 3
    elif num_filts == 12 :
        nrows, ncols = 3, 4
    elif num_filts == 16 :
        nrows, ncols = 4, 4
    else :
        nrows, ncols = 4, 5
    
    # outDir = '{}/pngs'.format(cluster)
    outDirAlt = '{}/images_cutouts_segmapped'.format(cluster)
    
    # os.makedirs('{}'.format(outDir), exist_ok=True) # ensure the output
        # directory for the pngs is available
    os.makedirs('{}'.format(outDirAlt), exist_ok=True)
    
    # open the table of all the objects 
    table = Table.read('{}/{}_sample.fits'.format(cluster, cluster))
    
    # determine the band flags that are in the catalog
    band_flags = [string for string in table.colnames if 'flag_F' in string]
    
    # use only those columns to create a new subtable
    flag_table = Table([table[band_flag] for band_flag in band_flags],
                       names=tuple(band_flags))
    
    # create a list of the flags per band for each galaxy
    flags = []
    for row in flag_table.iterrows() :
        flags.append(list(row))
    
    # next two lines are temporary - May 9th, 2022
    IDs = table['id']
    IDs = list(IDs[IDs > 20000])
    
    # IDs = list(table['id'])
    for i in range(len(IDs)) :
        segPath = '{}/cutouts/{}_ID_{}_segmap.fits'.format(cluster, cluster,
                                                           IDs[i])
        segMap = core.open_cutout(segPath, simple=True)
        
        bCGsegPath = '{}/cutouts/{}_ID_{}_segmap-bCG.fits'.format(
            cluster, cluster, IDs[i])
        bCGsegMap = core.open_cutout(bCGsegPath, simple=True)
        
        # outfile = outDir + '/{}_ID_{}.png'.format(cluster, IDs[i])
        outfile_segmapped = outDirAlt + '/{}_ID_{}.png'.format(cluster, IDs[i])
        
        # cutout_data = []
        cutout_segmapped = []
        for filt in filters :
            infile = '{}/cutouts/{}_ID_{}_{}.fits'.format(cluster, cluster,
                                                          IDs[i], filt)
            data = core.open_cutout(infile, simple=True)
            # cutout_data.append(data)
            
            if IDs[i] > 20000 :
                mask = (segMap > 0) | ((bCGsegMap > 0) & (bCGsegMap != IDs[i]))
            else :
                mask = (segMap > 0) & (segMap != IDs[i])
            segmapped_data = data.copy()
            segmapped_data[mask] = 0
            cutout_segmapped.append(segmapped_data)
        
        # plt.display_cutouts(cutout_data, nrows, ncols, filters, flags[i],
        #                     outfile, save=True)
        plt.display_cutouts(cutout_segmapped, nrows, ncols, filters, flags[i],
                            outfile_segmapped, save=True)
    
    return

def save_sbps(cluster, subpop) :
    
    outDir = '{}/images_sbps'.format(cluster)
    os.makedirs(outDir, exist_ok=True) # ensure the output direc. is available
    
    # open the table of all the objects
    if subpop == 'non-bCG' :
        clustable = Table.read('{}/{}_non-bCG_QGs.fits'.format(cluster,
                                                               cluster))
    
    # clustable = Table.read('{}/{}_sample.fits'.format(cluster, cluster))
    # mask the table based on the population of interest
    # clustable = clustable[clustable['pop'] == population]
    
    IDs = list(clustable['id'])
    
    for i in range(len(IDs)) :
        infile = '{}/photometry/{}_ID_{}_photometry.fits'.format(cluster,
                                                                 cluster,
                                                                 IDs[i])
        outfile = '{}/images_sbps/{}_ID_{}.png'.format(cluster, cluster, IDs[i])
        
        try :
            table = Table.read(infile)
            flux_columns = [col for col in table.colnames if col.endswith('_flux')]
            nPix_columns = [col.replace('flux', 'nPix') for col in flux_columns]
            filternames = [col.replace('_flux', '') for col in flux_columns]
            
            if len(filternames) == 9 :
                colors = ['violet',        # F275W
                          'mediumorchid',  # F336W
                          'darkslateblue', # F435W
                          'deepskyblue',   # F606W
                          'lime',          # F814W
                          'yellow',        # F105W
                          'darkorange',    # F125W
                          'orangered',     # F140W
                          'red']           # F160W
            if len(filternames) == 12 :
                colors = ['violet',        # F275W
                          'mediumorchid',  # F336W
                          'darkslateblue', # F435W
                          'royalblue',     # F475W
                          'deepskyblue',   # F606W
                          'cyan',          # F625W
                          'lime',          # F814W
                          'yellow',        # F105W
                          'gold',          # F110W
                          'darkorange',    # F125W
                          'orangered',     # F140W
                          'red']           # F160W
            if len(filternames) == 16 :
                colors = ['hotpink',       # F225W
                          'violet',        # F275W
                          'mediumorchid',  # F336W
                          'darkviolet',    # F390W
                          'darkslateblue', # F435W
                          'royalblue',     # F475W
                          'deepskyblue',   # F606W
                          'cyan',          # F625W
                          'springgreen',   # F775W
                          'lime',          # F814W
                          'greenyellow',   # F850LP
                          'yellow',        # F105W
                          'gold',          # F110W
                          'darkorange',    # F125W
                          'orangered',     # F140W
                          'red']           # F160W
            if len(filternames) == 17 :
                colors = ['hotpink',       # F225W
                          'violet',        # F275W
                          'mediumorchid',  # F336W
                          'darkviolet',    # F390W
                          'darkslateblue', # F435W
                          'royalblue',     # F475W
                          'dodgerblue',    # F555W
                          'deepskyblue',   # F606W
                          'cyan',          # F625W
                          'springgreen',   # F775W
                          'lime',          # F814W
                          'greenyellow',   # F850LP
                          'yellow',        # F105W
                          'gold',          # F110W
                          'darkorange',    # F125W
                          'orangered',     # F140W
                          'red']           # F160W
            
            xs = table['sma']/table['R_e']
            x_max = max(xs) + 1
            
            ys = []
            for j in range(len(flux_columns)) :
                maggies_per_pix = list(table[flux_columns[j]]/3631/table[nPix_columns[j]])
                ys.append(maggies_per_pix)
            
            plt.plot_sed(xs, ys, filternames, colors, outfile,
                         xlabel=r'Radius ($R_{\rm e}$)',
                         ylabel=r'Spatial Flux Density (maggies pixel$^{-1}$)',
                         xmax=x_max, save=True)
        
        except FileNotFoundError : # passes over A1063 ID 4746 and ID 5638
            pass
    
    return

def save_seds(cluster, subpop) :
    
    outDir = '{}/images_seds'.format(cluster)
    os.makedirs(outDir, exist_ok=True) # ensure the output direc. is available
    
    # open the table of all the objects
    if subpop == 'non-bCG' :
        clustable = Table.read('{}/{}_non-bCG_QGs.fits'.format(cluster,
                                                               cluster))
    
    # clustable = Table.read('{}/{}_sample.fits'.format(cluster, cluster))
    # mask the table based on the population of interest
    # clustable = clustable[clustable['pop'] == population]
    
    IDs = list(clustable['id'])
    
    for i in range(len(IDs)) :
        infile = '{}/photometry/{}_ID_{}_photometry.fits'.format(cluster,
                                                                 cluster,
                                                                 IDs[i])
        outfile = '{}/images_seds/{}_ID_{}.png'.format(cluster, cluster, IDs[i])
        
        try :
            table = Table.read(infile)
            flux_columns = [col for col in table.colnames if col.endswith('_flux')]
            nPix_columns = [col.replace('flux', 'nPix') for col in flux_columns]
            filternames = ['hff_' + col.replace('_flux', '') for col in flux_columns]
            
            colors = ['red', 'orangered', 'darkorange', 'gold', 'yellow',
                      'greenyellow', 'lime', 'springgreen', 'cyan',
                      'deepskyblue', 'dodgerblue', 'royalblue',
                      'darkslateblue', 'darkviolet', 'mediumorchid', 'violet',
                      'hotpink', 'red', 'orangered', 'darkorange', 'gold',
                      'yellow', 'greenyellow', 'lime', 'springgreen', 'cyan',
                      'deepskyblue', 'dodgerblue', 'royalblue',
                      'darkslateblue', 'darkviolet', 'mediumorchid', 'violet',
                      'hotpink', 'red', 'orangered', 'darkorange', 'gold',
                      'yellow', 'greenyellow', 'lime', 'springgreen', 'cyan',
                      'deepskyblue', 'dodgerblue', 'royalblue',
                      'darkslateblue', 'darkviolet', 'mediumorchid', 'violet',
                      'hotpink']
            
            ys = []
            labels = []
            for j in range(len(table)) :
                fluxes = np.array(list(table[flux_columns][j]))/3631 # in maggies
                nPixels = np.array(list(table[nPix_columns][j]))
                
                maggies_per_pix = fluxes/nPixels
                ys.append(maggies_per_pix)
                labels.append('bin {}'.format(j))
            
            filters = load_filters(filternames)
            xs = np.array([f.wave_effective for f in filters])
            
            x_min, x_max = np.min(xs)*0.8, np.max(xs)/0.8
            
            yflat = np.array(ys).flatten()
            ytemp = yflat[~np.isnan(yflat)]
            y_min, y_max = np.max([np.min(ytemp), 1e-15])*0.8, np.max(ytemp)/0.4
            
            plt.plot_sed(xs, ys, labels, colors, outfile,
                          xlabel=r'Wavelength ($\rm \AA$)',
                          ylabel=r'Spatial Flux Density (maggies pixel$^{-1}$)',
                          xmin=x_min, xmax=x_max, ymin=y_min, ymax=y_max, save=True)
        
        except FileNotFoundError : # passes over A1063 ID 4746 and ID 5638
            pass
    
    return
