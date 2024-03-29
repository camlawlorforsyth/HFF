
import os
from os.path import exists
import glob
import numpy as np

from astropy.table import Table, vstack
import prospect.io.read_results as reader
from prospect.plotting.sfh import parametric_mwa as mwa_calc
from prospect.plotting.utils import sample_posterior
from scipy.special import gamma, gammainc

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
    
    infile = '{}/h5/{}_ID_{}_bin_{}{}.h5'.format(cluster, cluster, ID, binNum,
                                                 version)
    
    result, obs, _ = reader.results_from(infile, dangerous=True)
    model = reader.get_model(result)
    sps = reader.get_sps(result)
    
    # t_hr = result['sampling_duration']/3600
    # print('The sampling took {:.2f} hours'.format(t_hr))
    
    return result, obs, model, sps

def get_results_integrated(cluster, ID) :
    
    infile = '{}/h5/{}_ID_{}_integrated.h5'.format(cluster, cluster, ID)
    
    result, obs, _ = reader.results_from(infile, dangerous=True)
    model = reader.get_model(result)
    sps = reader.get_sps(result)
    
    return result, obs, model, sps

def results_corner(cluster, ID, binNum, full=False, results_type='dynesty',
                   save=False, version='') :
    # create corner plots
    
    result, obs, model, sps = get_results(cluster, ID, binNum,
                                          results_type=results_type,
                                          version=version)
    
    if save :
        os.makedirs('{}/images_corner_plots'.format(cluster), # ensure the output
                    exist_ok=True) # directory for the figures is available
    outfile = '{}/images_corner_plots/{}_ID_{}_bin_{}{}.pdf'.format(cluster,
                                                                    cluster,
                                                                    ID, binNum,
                                                                    version)
    
    samples = sample_posterior(result['chain'], weights=result['weights'])
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

def results_corner_integrated(cluster, ID) :
    
    result, obs, model, sps = get_results_integrated(cluster, ID)
    
    samples = sample_posterior(result['chain'], weights=result['weights'])
    
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
    d_ratio_lo, d_ratio_hi = np.percentile(samples[:, 6], [0.15, 99.85])
    d_index_lo, d_index_hi = np.percentile(samples[:, 7], [0.15, 99.85])
    MWA_lo, MWA_hi = np.percentile(samples[:, 8], [0.15, 99.85])
    
    sub_samples = np.delete(samples, [4, 5], 1)
    sub_labels = [r'$z_{\rm red}$', r'$\log(M_{*})$', r'$\log(Z_{*})$',
                  r'$\hat{\tau}_{\lambda, 2}$', 'ratio', 'index', 'MWA']
    sub_ranges = [(z_lo, z_hi), (M_lo, M_hi), (Z_lo, Z_hi),
                  (d_lo, d_hi), (d_ratio_lo, d_ratio_hi),
                  (d_index_lo, d_index_hi), (MWA_lo, MWA_hi)]
    plt.plot_corner(sub_samples, sub_labels, len(sub_labels), result,
                    ranges=sub_ranges, save=False)
    
    return

def results_sed(cluster, ID, binNum, results_type='dynesty', save=False,
                version='',) :
    # plot SED and residuals
    
    result, obs, model, sps = get_results(cluster, ID, binNum,
                                          results_type=results_type,
                                          version=version)
    
    mapfile = '{}/map_arrays/{}_ID_{}_bin_{}{}.npz'.format(cluster, cluster, ID,
                                                           binNum, version)
    mapfilenpz = np.load(mapfile)
    map_spec, map_phot = mapfilenpz['map_spec'], mapfilenpz['map_phot']
    
    model_waves = sps.wavelengths*(1.0 + model.params.get('zred', 0.0))
    waves, mask = obs['phot_wave'], obs['phot_mask']
    fluxes, e_fluxes = obs['maggies'], obs['maggies_unc']
    
    chisq = np.sum(np.square((fluxes - map_phot)/e_fluxes))
    df = len(fluxes[mask]) - len(result['theta_labels'])
    
    if save :
        os.makedirs('{}/images_fitted_seds'.format(cluster), # ensure the output
                    exist_ok=True) # directory for the figures is available
    outfile = '{}/images_fitted_seds/{}_ID_{}_bin_{}{}.pdf'.format(cluster,
                                                                   cluster,
                                                                   ID, binNum,
                                                                   version)
    
    plt.plot_sed_from_fit(waves, fluxes, e_fluxes, mask, map_spec, map_phot,
                          model_waves, #chisq=chisq/df,
                          title=outfile,
                          save=save)
    
    # print('{:.2f}'.format(chisq/df))
    # print(fluxes[0]/e_fluxes[0])
    
    return

def save_map_models(cluster, ID, binNum, mapfile=None, results_type='dynesty',
                    verbose=False, version='') :
    '''
    Save the set of parameters of the sample with the highest posterior
    probability, as discussed here:
        https://github.com/bd-j/prospector/issues/215
    
    Note that in the prospector documentation, `theta_best` is the same as
    `theta_max`, as is `theta_map`. For all cases, the result is also stored in
    `result['bestfit']['parameter']`.
    '''
    
    result, obs, model, sps = get_results(cluster, ID, binNum,
                                          results_type=results_type,
                                          version=version)
    
    # find the MAP parameter set
    imax = np.argmax(result['lnprobability'])
    if results_type == 'emcee' :
        i, j = np.unravel_index(imax, result['lnprobability'].shape)
        theta_max = result['chain'][i, j, :].copy()
    else :
        theta_max = result['chain'][imax, :]
    
    map_spec, map_phot, _ = model.mean_model(theta_max, obs, sps=sps)
    
    # now save the MAP parameter set to file for future use
    os.makedirs('{}/map_arrays'.format(cluster), # ensure the output
                exist_ok=True) # directory for the arrays is available
    
    mapfile = '{}/map_arrays/{}_ID_{}_bin_{}{}.npz'.format(cluster, cluster, ID,
                                                           binNum, version)
    np.savez(mapfile, map_spec=map_spec, map_phot=map_phot)
    
    if verbose :
        print('Values saved to arrays.')
    
    return

def test(cluster, ID, version='') :
    
    import astropy.units as u
    
    from prospect.models.transforms import logsfr_ratios_to_sfrs
    from prospect.sources.constants import cosmo
    
    result, obs, model, sps = get_results(cluster, ID, 0, version=version)
    # print(result['theta_labels'])
    zred = result['model_params'][0]['init']
    max_age = cosmo.age(zred).value
    # agebins = result['model_params'][8]['init']
    
    agelims = np.concatenate([np.log10([1e-9, 0.03, 0.1, 0.5]),
                              np.linspace(0, np.log10(0.95*max_age), 6),
                              np.log10([max_age])]) + 9
    agebins = np.array([agelims[:-1], agelims[1:]]).T
    
    bestfit = result['bestfit']
    
    # zred_fit = bestfit['parameter'][0]
    logmass = bestfit['parameter'][3]
    logsfr_ratios = bestfit['parameter'][4:]
    
    logmass = np.full(logsfr_ratios.shape, logmass)
    
    from prospect.plotting.sfh import nonpar_mwa
    
    mwa = nonpar_mwa(logmass, logsfr_ratios, agebins)
    print(mwa)
    
    '''
    # Cam's attempt at MWA calculation
    agelims = np.power(10, agelims)
    agelims[0] = 0
    
    sfrs = logsfr_ratios_to_sfrs(logmass=logmass, logsfr_ratios=logsfr_ratios,
                                 agebins=agebins)
    ts = np.linspace(0, agelims[-1], int(1e3))
    sfr_vec = []
    for i in range(10) :
        if i == 9 :
            mask = (ts >= agelims[i]) & (ts < agelims[i+1] + 0.1)
        else :
            mask = (ts >= agelims[i]) & (ts < agelims[i+1])
        
        length = np.sum(mask)
        
        sfr_val = [sfrs[i]]*length
        
        sfr_vec.append(sfr_val)
    
    sfr_vec = np.concatenate(sfr_vec).ravel()
    
    ts = (ts*u.yr).to(u.Gyr).value
    
    plt.plot_simple_dumb(ts, sfr_vec, xlabel='Lookback Time (Gyr)',
                         ylabel=r'SFR $(M_{\odot}~{\rm yr}^{-1})$')
    MWA = np.trapz(ts*sfr_vec, ts)/np.trapz(sfr_vec, ts)
    # print(MWA)
    '''
    
    
    
    
    
    # samples = sample_posterior(result['chain'], weights=result['weights'])
    
    # labels = [r'$z_{\rm red}$', r'$\log(Z_{*})$', r'$\hat{\tau}_{\lambda, 2}$',
    #           r'$\log(M_{*})$',
    #           'ratio1', 'ratio2', 'ratio3', 'ratio4', 'ratio5',
    #           'ratio6', 'ratio7', 'ratio8', 'ratio9'] # 'MWA'
    
    # z_lo, z_hi = np.percentile(samples[:, 0], [0.15, 99.85])
    # Z_lo, Z_hi = np.percentile(samples[:, 1], [0.15, 99.85])
    # d_lo, d_hi = np.percentile(samples[:, 2], [0.15, 99.85])
    # M_lo, M_hi = np.percentile(samples[:, 3], [0.15, 99.85])
    # rat1_lo, rat1_hi = np.percentile(samples[:, 4], [2.5, 97.5])
    # rat2_lo, rat2_hi = np.percentile(samples[:, 5], [2.5, 97.5])
    # rat3_lo, rat3_hi = np.percentile(samples[:, 6], [2.5, 97.5])
    # rat4_lo, rat4_hi = np.percentile(samples[:, 7], [2.5, 97.5])
    # rat5_lo, rat5_hi = np.percentile(samples[:, 8], [2.5, 97.5])
    # rat6_lo, rat6_hi = np.percentile(samples[:, 9], [2.5, 97.5])
    # rat7_lo, rat7_hi = np.percentile(samples[:, 10], [2.5, 97.5])
    # rat8_lo, rat8_hi = np.percentile(samples[:, 11], [2.5, 97.5])
    # rat9_lo, rat9_hi = np.percentile(samples[:, 12], [2.5, 97.5])
    
    # ranges = [(z_lo, z_hi), (Z_lo, Z_hi), (d_lo, d_hi), (M_lo, M_hi), 
    #           (rat1_lo, rat1_hi), (rat2_lo, rat2_hi), (rat3_lo, rat3_hi),
    #           (rat4_lo, rat4_hi), (rat5_lo, rat5_hi), (rat6_lo, rat6_hi),
    #           (rat7_lo, rat7_hi), (rat8_lo, rat8_hi), (rat9_lo, rat9_hi)]
    
    # logify the mass and tau parameters
    # samples[:, 1] = np.log10(samples[:, 1])
    # samples[:, 5] = np.log10(samples[:, 5])
    
    # plt.plot_corner(samples, labels, len(labels), result, ranges=ranges,
    #                 outfile='output/nonparametric_corner_test.png', save=True)
    
    return

# results_corner('a370', 3337, 0)
# results_corner('a1063', 1366, 0)
# results_corner('a2744', 3859, 0)
# results_corner('a2744', 4369, 0)
# results_corner('a2744', 7427, 0)
# results_corner('m416', 5997, 0)
# results_corner('m1149', 1967, 0)
# results_corner('m1149', 1967, 1)

# results_sed('a370', 3337, 0)
# results_sed('a1063', 1366, 0)
# results_sed('a2744', 3859, 0) # F275W drives high chisq
# results_sed('a2744', 4369, 0) # F275W drives high chisq
# results_sed('a2744', 7427, 0)
# results_sed('m416', 5997, 0)
# results_sed('m1149', 1967, 0)
# results_sed('m1149', 1967, 1) # F225W drives high chisq

'''
# results_sed('a370', 3337, 0)
# results_sed('a370', 3337, 1)
# results_sed('a370', 3337, 2)
# results_sed('a370', 3337, 3)
# results_sed('a370', 3337, 4)
# results_sed('a370', 3337, 5)

# results_sed('a1063', 1366, 0)
# results_sed('a1063', 1366, 1)
# results_sed('a1063', 1366, 2)
# results_sed('a1063', 1366, 3)

# results_sed('a1063', 2455, 0)
# results_sed('a1063', 2455, 1)
# results_sed('a1063', 2455, 2)
# results_sed('a1063', 2455, 3)
# results_sed('a1063', 2455, 4)
# results_sed('a1063', 2455, 5)
# results_sed('a1063', 2455, 6)

# results_sed('a1063', 3550, 0)
# results_sed('a1063', 3550, 1)
# results_sed('a1063', 3550, 2)

# results_sed('a1063', 4823, 0)
# results_sed('a1063', 4823, 1)
# results_sed('a1063', 4823, 2)
# results_sed('a1063', 4823, 3)

# results_sed('a2744', 3859, 0)
# results_sed('a2744', 3859, 1)
# results_sed('a2744', 3859, 2)
# results_sed('a2744', 3859, 3)
# results_sed('a2744', 3859, 4)

# results_sed('a2744', 3964, 0)
# results_sed('a2744', 3964, 1)
# results_sed('a2744', 3964, 2)
# results_sed('a2744', 3964, 3)
# results_sed('a2744', 3964, 4)
# results_sed('a2744', 3964, 5)
# results_sed('a2744', 3964, 6)

# results_sed('a2744', 4173, 0)
# results_sed('a2744', 4173, 1)
# results_sed('a2744', 4173, 2)
# results_sed('a2744', 4173, 3)
# results_sed('a2744', 4173, 4)
# results_sed('a2744', 4173, 5)

# results_sed('a2744', 4369, 0)
# results_sed('a2744', 4369, 1)
# results_sed('a2744', 4369, 2)
# results_sed('a2744', 4369, 3)
# results_sed('a2744', 4369, 4)
# results_sed('a2744', 4369, 5)

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

# results_sed('a2744', 4862, 0)
# results_sed('a2744', 4862, 1)
# results_sed('a2744', 4862, 2)
# results_sed('a2744', 4862, 3)
# results_sed('a2744', 4862, 4)
# results_sed('a2744', 4862, 5)
# results_sed('a2744', 4862, 6)
# results_sed('a2744', 4862, 7)

# results_sed('a2744', 7427, 0)
# results_sed('a2744', 7427, 1)
# results_sed('a2744', 7427, 2)
# results_sed('a2744', 7427, 3)
# results_sed('a2744', 7427, 4)
# results_sed('a2744', 7427, 5)
# results_sed('a2744', 7427, 6)

# results_sed('m416', 5997, 0)
# results_sed('m416', 5997, 1)
# results_sed('m416', 5997, 2)
# results_sed('m416', 5997, 3)
# results_sed('m416', 5997, 4)
# results_sed('m416', 5997, 5)
# results_sed('m416', 5997, 6)
# results_sed('m416', 5997, 7)

# results_sed('m416', 6255, 0)
# results_sed('m416', 6255, 1)
# results_sed('m416', 6255, 2)
# results_sed('m416', 6255, 3)
# results_sed('m416', 6255, 4)
# results_sed('m416', 6255, 5)
# results_sed('m416', 6255, 6)
# results_sed('m416', 6255, 7)

# results_sed('m717', 861, 0)
# results_sed('m717', 861, 1)
# results_sed('m717', 861, 2)
# results_sed('m717', 861, 3)
# results_sed('m717', 861, 4)
# results_sed('m717', 861, 5)

# results_sed('m1149', 1967, 0)
# results_sed('m1149', 1967, 1)
# results_sed('m1149', 1967, 2)
# results_sed('m1149', 1967, 3)
# results_sed('m1149', 1967, 4)
# results_sed('m1149', 1967, 5)

# results_sed('m1149', 2403, 0)
# results_sed('m1149', 2403, 1)
# results_sed('m1149', 2403, 2)
# results_sed('m1149', 2403, 3)

# results_sed('m1149', 3531, 0)
# results_sed('m1149', 3531, 1)
# results_sed('m1149', 3531, 2)
# results_sed('m1149', 3531, 3)

# results_sed('m1149', 4246, 0)
# results_sed('m1149', 4246, 1)
# results_sed('m1149', 4246, 2)
# results_sed('m1149', 4246, 3)
# results_sed('m1149', 4246, 4)
# results_sed('m1149', 4246, 5)

# results_sed('m1149', 5095, 0)
# results_sed('m1149', 5095, 1)
# results_sed('m1149', 5095, 2)
'''

# results_sed('m416', 2876, 0)

# results for the flat gradient
# results_corner('a1063', 4823, 0)
# results_corner('a1063', 4823, 1)
# results_corner('a1063', 4823, 2)
# results_corner('a1063', 4823, 3)

# results_corner('a1063', 5156, 0)
# results_corner('a1063', 5156, 1)
# results_corner('a1063', 5156, 2)
# results_corner('a1063', 5156, 3)
# results_corner('a1063', 5156, 4)

# results_sed('m416', 5607, 0, version='_flat')
# results_sed('m416', 5607, 0, version='_gradient')
# results_corner('m416', 5607, 0)

# test('a2744', 3964, version='_nonpara')
# results_corner('a2744', 3964, 1, version='_flat')
# results_corner('a2744', 3964, 1, version='_gradient')

# results_corner_integrated('a370', 1328)
# results_corner_integrated('a370', 2590)
# results_corner_integrated('a370', 3089)
# results_corner_integrated('a370', 3767)
# results_corner_integrated('a370', 4259)
# results_corner_integrated('a370', 5177)

def results_integrated(result) :
    
    # get the total number of pixels involved in the fit, as the fits
    # were completed on a "per pixel" basis
    nPixels = result['obs']['nPixels'][-1] # use the value for the F160W filter
    
    # get the samples
    samples = sample_posterior(result['chain'], weights=result['weights'])
    
    # calculate the MWAs, and add those to the samples
    mwas = mwa_calc(samples[:, 5], samples[:, 4], power=1)
    samples = np.c_[samples, mwas]
    
    # logify the mass and tau parameters
    samples[:, 1] = np.log10(samples[:, 1]) + np.log10(nPixels)
    samples[:, 5] = np.log10(samples[:, 5])
    
    z_16, z_50, z_84 = np.percentile(samples[:, 0], [16, 50, 84])
    logM_16, logM_50, logM_84 = np.percentile(samples[:, 1], [16, 50, 84])
    logZ_16, logZ_50, logZ_84 = np.percentile(samples[:, 2], [16, 50, 84])
    Av_16, Av_50, Av_84 = np.percentile(samples[:, 3], [16, 50, 84])
    tage_16, tage_50, tage_84 = np.percentile(samples[:, 4], [16, 50, 84])
    logT_16, logT_50, logT_84 = np.percentile(samples[:, 5], [16, 50, 84])
    ratio_16, ratio_50, ratio_84 = np.percentile(samples[:, 6], [16, 50, 84])
    index_16, index_50, index_84 = np.percentile(samples[:, 7], [16, 50, 84])
    MWA_16, MWA_50, MWA_84 = np.percentile(samples[:, 8], [16, 50, 84])
    
    return (z_16, z_50, z_84, logM_16, logM_50, logM_84,
            logZ_16, logZ_50, logZ_84, Av_16, Av_50, Av_84,
            tage_16, tage_50, tage_84, logT_16, logT_50, logT_84,
            ratio_16, ratio_50, ratio_84, index_16, index_50, index_84,
            MWA_16, MWA_50, MWA_84)

def build_summary_results() :
    
    outfile = 'output/tables/integrated_results_summary.fits'
    
    table = Table.read('output/tables/sample_final.fits')
    
    # the values that will be in the final table
    clus, ids = [], []
    logM_16s, logM_50s, logM_84s = [], [], []
    z_16s, z_50s, z_84s = [], [], []
    logM_16s, logM_50s, logM_84s = [], [], []
    logZ_16s, logZ_50s, logZ_84s = [], [], []
    Av_16s, Av_50s, Av_84s = [], [], []
    tage_16s, tage_50s, tage_84s = [], [], []
    logT_16s, logT_50s, logT_84s = [], [], []
    ratio_16s, ratio_50s, ratio_84s = [], [], []
    index_16s, index_50s, index_84s = [], [], []
    MWA_16s, MWA_50s, MWA_84s = [], [], []
    
    for cluster, ID in zip(table['cluster'], table['id']) :
        # get the results
        result, obs, model, sps = get_results_integrated(cluster, ID)
        
        (z_16, z_50, z_84, logM_16, logM_50, logM_84,
         logZ_16, logZ_50, logZ_84, Av_16, Av_50, Av_84,
         tage_16, tage_50, tage_84, logT_16, logT_50, logT_84,
         ratio_16, ratio_50, ratio_84, index_16, index_50, index_84,
         MWA_16, MWA_50, MWA_84) = results_integrated(result)
        
        # append the calculated values
        clus.append(cluster)
        ids.append(ID)
        
        z_16s.append(z_16)
        z_50s.append(z_50)
        z_84s.append(z_84)
        
        logM_16s.append(logM_16)
        logM_50s.append(logM_50)
        logM_84s.append(logM_84)
        
        logZ_16s.append(logZ_16)
        logZ_50s.append(logZ_50)
        logZ_84s.append(logZ_84)
        
        Av_16s.append(Av_16)
        Av_50s.append(Av_50)
        Av_84s.append(Av_84)
        
        tage_16s.append(tage_16)
        tage_50s.append(tage_50)
        tage_84s.append(tage_84)
        
        logT_16s.append(logT_16)
        logT_50s.append(logT_50)
        logT_84s.append(logT_84)
        
        ratio_16s.append(ratio_16)
        ratio_50s.append(ratio_50)
        ratio_84s.append(ratio_84)
        
        index_16s.append(index_16)
        index_50s.append(index_50)
        index_84s.append(index_84)
        
        MWA_16s.append(MWA_16)
        MWA_50s.append(MWA_50)
        MWA_84s.append(MWA_84)
    
    table = Table([clus, ids,
                   z_16s, z_50s, z_84s,
                   logM_16s, logM_50s, logM_84s,
                   logZ_16s, logZ_50s, logZ_84s,
                   Av_16s, Av_50s, Av_84s,
                   tage_16s, tage_50s, tage_84s,
                   logT_16s, logT_50s, logT_84s,
                   ratio_16s, ratio_50s, ratio_84s,
                   index_16s, index_50s, index_84s,
                   MWA_16s, MWA_50s, MWA_84s],
                  names=('cluster', 'ID',
                         'z_16', 'z_50', 'z_84',
                         'logM_16', 'logM_50', 'logM_84',
                         'logZ_16', 'logZ_50', 'logZ_84',
                         'Av_16', 'Av_50', 'Av_84',
                         'tage_16', 'tage_50', 'tage_84',
                         'logT_16', 'logT_50', 'logT_84',
                         'ratio_16', 'ratio_50', 'ratio_84',
                         'index_16', 'index_50', 'index_84',
                         'MWA_16', 'MWA_50', 'MWA_84'))
    table.write(outfile)
    
    return

def append_summary_results(outfile) :
    
    table = Table.read(outfile)
    
    table_vals = [clus + '_' + ID for clus, ID in zip(table['cluster'], table['ID'])]
    fitted_vals = [clus + '_' + ID for clus, ID in zip(fitted_clusters, fitted_IDs)]
    
    new_table = Table(table.columns)
    
    for fitted_val in fitted_vals :
        if fitted_val not in table_vals :
            new_clus, new_ID = fitted_val.split('_')
            
            # get the results
            result, _, _, _ = get_results_integrated(new_clus, new_ID)
            
            # add those results as a new row
            new_row = [new_clus, new_ID] + list(results_integrated(result))
            new_table.add_row(new_row)
    
    if len(new_table) > 0 :
        table = vstack([table, new_table])
        table.write(outfile, overwrite=True)
    
    return

def plot_summary_results() :
    
    table = Table.read('output/tables/integrated_results_summary.fits')
    
    logM = table['logM_50']
    ind_16, ind_50, ind_84 = table['index_16'], table['index_50'], table['index_84']
    
    # slope, intercept = np.polyfit(logM, ind_50, 1)
    # print(slope, intercept)
    
    plt.plot_scatter_err_asymm(logM, ind_50,
                                # np.zeros_like(logM),
                                # np.zeros_like(logM),
                               ind_50 - ind_16,
                               ind_84 - ind_50,
                               'k', 'o',
                               xlabel=r'$\log(M_{*}/M_{\odot})$',
                               ylabel='Calzetti offset index',
                               xmin=8, xmax=13, ymin=-2.1, ymax=0.5,
                               figsizewidth=12, figsizeheight=9)
    
    return

def SFMS() :
    
    table = Table.read('output/tables/integrated_results_summary.fits')
    
    logmass = table['logM_50']
    mass = np.power(10, logmass)
    tage = table['tage_50']
    tau = np.power(10, table['logT_50'])
    
    sfr = calculate_sfr(mass, tage, tau)
    
    norm = np.array( tau*gamma(2)*gammainc(2, tage/tau) )
    
    # alt = []
    # for t, T in zip(tage, tau) :
    #     times = np.linspace(0, t, 1000)
    #     alt.append(np.trapz(sfr_calc(times, T), times))
    # alt = np.array(alt)
    
    
    
    # print(list(np.sort(gammainc(2, tage/tau))))
    
    # print(list(np.sort(norm)))
    
    # table['sfr'] = sfr
    # table.write('output/tables/integrated_results_summary_with_SFR.fits')
    
    
    # plt.plot_scatter(logmass, sfr, 'k', '', 'o',
    #                  xmin=8, xmax=13, ymin=-4, ymax=4)
    
    
    
    
    return

def sfr_calc(times, tau) :
    return (times/tau)*np.exp(-times/tau)

def calculate_sfr(mass, tage, tau) :
    return mass*tage/np.square(tau)*np.exp(-tage/tau)/(gamma(2)*gammainc(2, tage/tau))*1e-9



SFMS()

