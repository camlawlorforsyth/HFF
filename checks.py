
import os
import numpy as np

from astropy.table import Table, vstack
from scipy.stats import kde

import plotting as plt

def check_bins(cluster, plot=True) :
    
    photDir = '{}/photometry'.format(cluster)
    files = os.listdir(photDir)
    
    lengths = []
    for file in files :
        table = Table.read('{}/{}'.format(photDir, file))
        length = len(table)
        lengths.append(length)
    
    if plot :
        numBins = int(np.ceil(1.3*np.sqrt(len(lengths))))
        plt.histogram(lengths, 'Number of Annuli per Galaxy', title=cluster,
                      bins=numBins)
    
    return

def concatenate_all() :
    
    clusters = ['a370', 'a1063', 'a2744', 'm416', 'm717', 'm1149']
    all_clusters_list = []
    all_parallels_list = []
    for cluster in clusters :
        cluster_table = Table.read('{}/{}_final_objects.fits'.format(
            cluster, cluster))
        cluster_table['cluster'] = cluster
        all_clusters_list.append(cluster_table)
        parallel_table = Table.read('{}par/{}par_final_objects.fits'.format(
            cluster, cluster))
        parallel_table['field'] = '{}par'.format(cluster)
        all_parallels_list.append(parallel_table)
    
    all_clusters = vstack(all_clusters_list)
    all_parallels = vstack(all_parallels_list)
    
    return all_clusters, all_parallels

def plot_all_FUVVJ() :
    
    all_clusters, all_parallels = concatenate_all()
    
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
                              xlabel=r'$V - J$', ylabel=r'$FUV - V$',
                              xmin=0, xmax=2.1, ymin=0, ymax=8.4, loc=2)

def plot_all_UVJ() :
    
    all_clusters, all_parallels = concatenate_all()
    
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
                              xlabel=r'$V - J$', ylabel=r'$U - V$',
                              xmin=0, xmax=2.1, ymin=0, ymax=2.4, loc=4)
    
    return


def plot_parallel_objects_all() :
    
    all_clusters, all_parallels = concatenate_all()
    
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
