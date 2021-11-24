
import os
import numpy as np

import prospect.io.read_results as reader
from prospect.plotting.sfh import parametric_mwa as mwa_calc
from prospect.plotting.utils import sample_posterior

import plotting as plt

''' # to review and delete - believed to be mostly redundant
def show_results(cluster, ID, binNum, results_type='dynesty', save=False) :
    
    result, obs, model, sps = get_results(cluster, ID, binNum,
                                          results_type=results_type)
    
    # plot SED with only photometry
    # plt.plot_sed_alone(obs)
    
    # plot SED with initial theta
    # theta = model.theta.copy()
    # initial_spec, initial_phot, _ = model.sed(theta, obs=obs, sps=sps)
    # plt.plot_sed_and_model(obs, model, initial_phot, initial_spec, sps)
    
    # investigate the parameter traces
    # plt.plot_chains(result, results_type=results_type)
    
    return
'''

def get_results(cluster, ID, binNum, results_type='dynesty', version='') :
    # convenience function
    
    if (cluster == 'a370') and (ID == 3337) and (binNum == 0) :
        infile = '{}/h5/ID_{}_bin_{}{}.h5'.format(cluster, ID, binNum, version)
    elif (cluster == 'a1063') and (ID == 1366) and (binNum == 0) :
        infile = '{}/h5/ID_{}_bin_{}{}.h5'.format(cluster, ID, binNum, version)
    elif (cluster == 'a2744') and (ID == 3859) and (binNum == 0) :
        infile = '{}/h5/ID_{}_bin_{}{}.h5'.format(cluster, ID, binNum, version)
    elif (cluster == 'a2744') and (ID == 4369) and (binNum == 0) :
        infile = '{}/h5/ID_{}_bin_{}{}.h5'.format(cluster, ID, binNum, version)
    elif (cluster == 'a2744') and (ID == 7427) and (binNum == 0) :
        infile = '{}/h5/ID_{}_bin_{}{}.h5'.format(cluster, ID, binNum, version)
    elif (cluster == 'm416') and (ID == 5997) and (binNum == 0) :
        infile = '{}/h5/ID_{}_bin_{}{}.h5'.format(cluster, ID, binNum, version)
    elif (cluster == 'm1149') and (ID == 1967) and (binNum == 0) :
        infile = '{}/h5/ID_{}_bin_{}{}.h5'.format(cluster, ID, binNum, version)
    elif (cluster == 'm1149') and (ID == 1967) and (binNum == 1) :
        infile = '{}/h5/ID_{}_bin_{}{}.h5'.format(cluster, ID, binNum, version)
    else :
        infile = '{}/h5/ID_{}_bin_{}.h5'.format(cluster, ID, binNum)
    
    result, obs, _ = reader.results_from(infile, dangerous=True)
    model = reader.get_model(result)
    sps = reader.get_sps(result)
    
    # t_hr = result['sampling_duration']/3600
    # print('The sampling took {:.2f} hours'.format(t_hr))
    
    return result, obs, model, sps

def results_corner(cluster, ID, binNum, results_type='dynesty', save=False,
                   version='', full=False) :
    
    # create corner plots
    result, obs, model, sps = get_results(cluster, ID, binNum,
                                          results_type=results_type,
                                          version=version)
    
    if save :
        os.makedirs('{}/pngs_corner'.format(cluster), # ensure the output
                    exist_ok=True) # directory for the figures is available
    outfile = '{}/pngs_corner/{}_ID_{}_bin_{}{}.pdf'.format(cluster, cluster,
                                                          ID, binNum, version)
    
    samples = sample_posterior(result['chain'], weights=result['weights'],
                               nsample=len(result['chain']))
    mwas = mwa_calc(samples[:, 5], samples[:, 4], power=1)
    
    # logify the mass and tau parameters
    samples[:, 1] = np.log10(samples[:, 1])
    samples[:, 5] = np.log10(samples[:, 5])
    
    samples = np.c_[samples, mwas]
    
    z_lo, z_hi = np.percentile(samples[:, 0], [0.15, 99.85])
    M_lo, M_hi = np.percentile(samples[:, 1], [0.15, 99.85])
    Z_lo, Z_hi = np.percentile(samples[:, 2], [0.15, 99.85])
    d_lo, d_hi = np.percentile(samples[:, 3], [0.15, 99.85])
    t_lo, t_hi = np.percentile(samples[:, 4], [0.15, 99.85])
    T_lo, T_hi = np.percentile(samples[:, 5], [0.15, 99.85])
    MWA_lo, MWA_hi = np.percentile(samples[:, 6], [0.15, 99.85])
    
    if full :
        ranges = [(z_lo, z_hi), (M_lo, M_hi), (Z_lo, Z_hi),
                  (d_lo, d_hi), (t_lo, t_hi), (T_lo, T_hi), (MWA_lo, MWA_hi)]
        labels = [r'$z_{\rm red}$', r'$\log(M_{*})$', r'$\log(Z_{*})$',
                  r'$\hat{\tau}_{\lambda, 2}$', r'$T_{0}$', r'$\log(\tau)$', 'MWA']
        plt.plot_corner(samples, labels, len(labels), result, ranges=ranges,
                        outfile=outfile, save=save)
    else :
        sub_samples = np.delete(samples, [4, 5], 1)
        sub_labels = [r'$z_{\rm red}$', r'$\log(M_{*})$', r'$\log(Z_{*})$',
                      r'$\hat{\tau}_{\lambda, 2}$', 'MWA']
        sub_ranges = [(z_lo, z_hi), (M_lo, M_hi), (Z_lo, Z_hi),
                      (d_lo, d_hi), (MWA_lo, MWA_hi)]
        plt.plot_corner(sub_samples, sub_labels, len(sub_labels), result,
                        ranges=sub_ranges, outfile=outfile, save=save)
    
    return

def get_map_models(result, obs, model, sps, cluster, results_type='dynesty',
                   load=False, mapfile=None, save=False) :
    '''
    Get the set of parameters of the sample with the highest posterior
    probability, as discussed here: https://github.com/bd-j/prospector/issues/215
    
    Note that in the prospector documentation, `theta_best` is the same as
    `theta_max`, as is `theta_map`. For all cases, the result is also stored in
    `result['bestfit']['parameter']`.
    
    Parameters
    ----------
    result : TYPE
        DESCRIPTION.
    obs : TYPE
        DESCRIPTION.
    model : TYPE
        DESCRIPTION.
    sps : TYPE
        DESCRIPTION.
    cluster : TYPE
        DESCRIPTION.
    results_type : TYPE, optional
        DESCRIPTION. The default is 'dynesty'.
    load : TYPE, optional
        DESCRIPTION. The default is False.
    mapfile : TYPE, optional
        DESCRIPTION. The default is None.
    save : TYPE, optional
        DESCRIPTION. The default is False.
    
    Returns
    -------
    map_spec : TYPE
        DESCRIPTION.
    map_phot : TYPE
        DESCRIPTION.
    
    '''
    if save :
        os.makedirs('{}/map_arrays'.format(cluster), # ensure the output
                    exist_ok=True) # directory for the arrays is available
        
        # find the MAP parameter set
        imax = np.argmax(result['lnprobability'])
        if results_type == 'emcee' :
            i, j = np.unravel_index(imax, result['lnprobability'].shape)
            theta_max = result['chain'][i, j, :].copy()
        else :
            theta_max = result['chain'][imax, :]
        
        map_spec, map_phot, _ = model.mean_model(theta_max, obs, sps=sps)
        
        np.savez(mapfile, map_spec=map_spec, map_phot=map_phot)
        
        print('Values saved to arrays.')
        
        return None, None
    
    if load :
        mapfilenpz = np.load(mapfile)
        map_spec, map_phot = mapfilenpz['map_spec'], mapfilenpz['map_phot']
        
        return map_spec, map_phot

def drop_band(waves, fluxes, e_fluxes, mask, map_phot, band) :
    return waves[band:], fluxes[band:], e_fluxes[band:], mask[band:], map_phot[band:]

def results_sed(cluster, ID, binNum, results_type='dynesty', version='',
                show=False) :
    # plot SED and residuals
    
    # get results
    result, obs, model, sps = get_results(cluster, ID, binNum,
                                          results_type=results_type,
                                          version=version)
    
    mapfile = '{}/map_arrays/{}_ID_{}_bin_{}{}.npz'.format(cluster, cluster, ID,
                                                           binNum, version)
    map_spec, map_phot = get_map_models(result, obs, model, sps, cluster,
                                        mapfile=mapfile, save=False, load=True)
    
    model_waves = sps.wavelengths*(1.0 + model.params.get('zred', 0.0))
    waves, mask = obs['phot_wave'], obs['phot_mask']
    fluxes, e_fluxes = obs['maggies'], obs['maggies_unc']
    
    # drop certain bands for figuring out what issue drives high chisq values
    waves, fluxes, e_fluxes, mask, map_phot = drop_band(waves, fluxes, e_fluxes,
                                                        mask, map_phot, 0)
    
    chisq = np.sum(np.square((fluxes - map_phot)/e_fluxes))
    df = len(fluxes[mask]) - len(result['theta_labels'])
    outfile = '{}_ID_{}_bin_{}_new'.format(cluster, ID, binNum)
    if show :
        plt.plot_sed_from_fit(waves, fluxes, e_fluxes, mask, map_spec, map_phot,
                              model_waves, chisq=chisq/df, title=outfile,
                              save=False)
    
    # print('{:.2f}'.format(chisq/df))
    # print(fluxes[0]/e_fluxes[0])
    
    return

# results_sed('a370', 3337, 0, version='_v2', show=True)
# results_sed('a1063', 1366, 0, version='_v2', show=True)
# results_sed('a2744', 3859, 0, version='_v2', show=True) # F275W drives high chisq
results_sed('a2744', 4369, 0, version='_v2', show=True) # F275W (and minorly F336W) drive high chisq
# results_sed('a2744', 7427, 0, version='_v2', show=True)
# results_sed('m416', 5997, 0, version='_v2', show=True)
# results_sed('m1149', 1967, 0, version='_v2', show=True)
# results_sed('m1149', 1967, 1, version='_v2', show=True) # F225W drives high chisq

'''
results_sed('a370', 3337, 0, version='_v2')
# results_sed('a370', 3337, 1)
# results_sed('a370', 3337, 2)
# results_sed('a370', 3337, 3)
# results_sed('a370', 3337, 4)
# results_sed('a370', 3337, 5)
# results_sed('a370', 3337, 6)

results_sed('a1063', 1366, 0, version='_v2')
# results_sed('a1063', 1366, 1)
# results_sed('a1063', 1366, 2)
# results_sed('a1063', 1366, 3)
# results_sed('a1063', 1366, 4)
# results_sed('a1063', 1366, 5)

# results_sed('a1063', 2455, 0)
# results_sed('a1063', 2455, 1)
# results_sed('a1063', 2455, 2)
# results_sed('a1063', 2455, 3)
# results_sed('a1063', 2455, 4)
# results_sed('a1063', 2455, 5)
# results_sed('a1063', 2455, 6)
# results_sed('a1063', 2455, 7)
# results_sed('a1063', 2455, 8)

# results_sed('a1063', 3550, 0)
# results_sed('a1063', 3550, 1)
# results_sed('a1063', 3550, 2)
# results_sed('a1063', 3550, 3)
# results_sed('a1063', 3550, 4)

# results_sed('a1063', 4823, 0)
# results_sed('a1063', 4823, 1)
# results_sed('a1063', 4823, 2)
# results_sed('a1063', 4823, 3)
# results_sed('a1063', 4823, 4)
# results_sed('a1063', 4823, 5)

results_sed('a2744', 3859, 0, version='_v2')
# results_sed('a2744', 3859, 1)
# results_sed('a2744', 3859, 2)
# results_sed('a2744', 3859, 3)
# results_sed('a2744', 3859, 4)
# results_sed('a2744', 3859, 5)
# results_sed('a2744', 3859, 6)

# results_sed('a2744', 3964, 0)
# results_sed('a2744', 3964, 1)
# results_sed('a2744', 3964, 2)
# results_sed('a2744', 3964, 3)
# results_sed('a2744', 3964, 4)
# results_sed('a2744', 3964, 5)
# results_sed('a2744', 3964, 6)
# results_sed('a2744', 3964, 7)
# results_sed('a2744', 3964, 8)
# results_sed('a2744', 3964, 9)

# results_sed('a2744', 4173, 0)
# results_sed('a2744', 4173, 1)
# results_sed('a2744', 4173, 2)
# results_sed('a2744', 4173, 3)
# results_sed('a2744', 4173, 4)
# results_sed('a2744', 4173, 5)
# results_sed('a2744', 4173, 6)
# results_sed('a2744', 4173, 7)
# results_sed('a2744', 4173, 8)
# results_sed('a2744', 4173, 9)

results_sed('a2744', 4369, 0, version='_v2')
# results_sed('a2744', 4369, 1)
# results_sed('a2744', 4369, 2)
# results_sed('a2744', 4369, 3)
# results_sed('a2744', 4369, 4)
# results_sed('a2744', 4369, 5)
# results_sed('a2744', 4369, 6)
# results_sed('a2744', 4369, 7)
# results_sed('a2744', 4369, 8)

# results_sed('a2744', 4765, 0)
# results_sed('a2744', 4765, 1)
# results_sed('a2744', 4765, 2)
# results_sed('a2744', 4765, 3)
# results_sed('a2744', 4765, 4)
# results_sed('a2744', 4765, 5)
# results_sed('a2744', 4765, 6)
# results_sed('a2744', 4765, 7)
# results_sed('a2744', 4765, 8)
# results_sed('a2744', 4765, 9)
# results_sed('a2744', 4765, 10)
# results_sed('a2744', 4765, 11)
# results_sed('a2744', 4765, 12)
# results_sed('a2744', 4765, 13)

# results_sed('a2744', 4862, 0)
# results_sed('a2744', 4862, 1)
# results_sed('a2744', 4862, 2)
# results_sed('a2744', 4862, 3)
# results_sed('a2744', 4862, 4)
# results_sed('a2744', 4862, 5)
# results_sed('a2744', 4862, 6)
# results_sed('a2744', 4862, 7)
# results_sed('a2744', 4862, 8)
# results_sed('a2744', 4862, 9)
# results_sed('a2744', 4862, 10)

results_sed('a2744', 7427, 0, version='_v2')
# results_sed('a2744', 7427, 1)
# results_sed('a2744', 7427, 2)
# results_sed('a2744', 7427, 3)
# results_sed('a2744', 7427, 4)
# results_sed('a2744', 7427, 5)
# results_sed('a2744', 7427, 6)
# results_sed('a2744', 7427, 7)
# results_sed('a2744', 7427, 8)
# results_sed('a2744', 7427, 9)

results_sed('m416', 5997, 0, version='_v2')
# results_sed('m416', 5997, 1)
# results_sed('m416', 5997, 2)
# results_sed('m416', 5997, 3)
# results_sed('m416', 5997, 4)
# results_sed('m416', 5997, 5)
# results_sed('m416', 5997, 6)
# results_sed('m416', 5997, 7)
# results_sed('m416', 5997, 8)
# results_sed('m416', 5997, 9)
# results_sed('m416', 5997, 10)

# results_sed('m416', 6255, 0)
# results_sed('m416', 6255, 1)
# results_sed('m416', 6255, 2)
# results_sed('m416', 6255, 3)
# results_sed('m416', 6255, 4)
# results_sed('m416', 6255, 5)
# results_sed('m416', 6255, 6)
# results_sed('m416', 6255, 7)
# results_sed('m416', 6255, 8)
# results_sed('m416', 6255, 9)
# results_sed('m416', 6255, 10)

# results_sed('m717', 861, 0)
# results_sed('m717', 861, 1)
# results_sed('m717', 861, 2)
# results_sed('m717', 861, 3)
# results_sed('m717', 861, 4)
# results_sed('m717', 861, 5)
# results_sed('m717', 861, 6)
# results_sed('m717', 861, 7)

results_sed('m1149', 1967, 0, version='_v2')
results_sed('m1149', 1967, 1, version='_v2')
# results_sed('m1149', 1967, 2)
# results_sed('m1149', 1967, 3)
# results_sed('m1149', 1967, 4)
# results_sed('m1149', 1967, 5)
# results_sed('m1149', 1967, 6)
# results_sed('m1149', 1967, 7)

# results_sed('m1149', 2403, 0)
# results_sed('m1149', 2403, 1)
# results_sed('m1149', 2403, 2)
# results_sed('m1149', 2403, 3)
# results_sed('m1149', 2403, 4)
# results_sed('m1149', 2403, 5)

# results_sed('m1149', 3531, 0)
# results_sed('m1149', 3531, 1)
# results_sed('m1149', 3531, 2)
# results_sed('m1149', 3531, 3)
# results_sed('m1149', 3531, 4)

# results_sed('m1149', 4246, 0)
# results_sed('m1149', 4246, 1)
# results_sed('m1149', 4246, 2)
# results_sed('m1149', 4246, 3)
# results_sed('m1149', 4246, 4)
# results_sed('m1149', 4246, 5)
# results_sed('m1149', 4246, 6)
# results_sed('m1149', 4246, 7)

# results_sed('m1149', 5095, 0)
# results_sed('m1149', 5095, 1)
# results_sed('m1149', 5095, 2)
# results_sed('m1149', 5095, 3)
# results_sed('m1149', 5095, 4)
'''