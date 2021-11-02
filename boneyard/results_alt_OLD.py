
import numpy as np
import matplotlib.pyplot as plt
import prospect.io.read_results as reader

results_type = 'emcee' # 'emcee' or 'dynesty'
verbose = True
plots = True

# get the results and observation dictionaries
result, obs, _ = reader.results_from('parametric/m416_ID_2876_BIN_10_TT1_mcmc.h5',
                                     dangerous=True)

# also get the sps and model objects. This works if using a parameter file with
# `build_*` methods
sps = reader.get_sps(result)
model = reader.get_model(result)

# print the contents of the `results` dictionary
if verbose :
    print(result.keys())

# STEP 1: investigate the parameter traces
if plots :
    if results_type == 'emcee' :
        chosen = np.random.choice(result['run_params']['nwalkers'], size=10,
                                  replace=False)
        tracefig = reader.traceplot(result, figsize=(16, 9), chains=chosen)
    else :
        tracefig = reader.traceplot(result, figsize=(16, 9))

# STEP 2: create corner plots
imax = np.argmax(result['lnprobability'])

if results_type == 'emcee' :
    i, j = np.unravel_index(imax, result['lnprobability'].shape)
    theta_max = result['chain'][i, j, :].copy()
    thin = 5
else :
    theta_max = result['chain'][imax, :]
    thin = 1

if verbose :
    print('MAP value: {}'.format(theta_max))

if plots :
    cornerfig = reader.subcorner(result, start=0, thin=thin,
                                 fig=plt.subplots(4,4, figsize=(9, 9))[0])
'''
# STEP 3: plot SED and residuals
randint = np.random.randint
if results_type == 'emcee' :
    nwalkers, niter = result['run_params']['nwalkers'], result['run_params']['niter']
    theta = result['chain'][randint(nwalkers), randint(niter)]
else:
    theta = result['chain'][randint(len(result['chain']))]

# get helpful info for plotting
a = 1.0 + model.params.get('zred', 0.0)
wspec = sps.wavelengths
obs['phot_wave'] = np.array([f.wave_effective for f in obs['filters']])
wphot = obs['phot_wave']

xmin, xmax = np.min(wphot)*0.8, np.max(wphot)/0.8
ymin, ymax = obs['maggies'].min()*0.8, obs['maggies'].max()/0.4

# generate models
mspec, mphot, mextra = model.mean_model(theta, obs, sps=sps)
mspec_map, mphot_map, _ = model.mean_model(theta_max, obs, sps=sps)

# make plot of data and model
if plots :
    plt.figure(figsize=(16, 9))
    
    plt.loglog(wspec, mspec, label='Model spectrum (random draw)',
               lw=0.7, color='navy', alpha=0.7)
    plt.loglog(wspec, mspec_map, label='Model spectrum (MAP)',
               lw=0.7, color='green', alpha=0.7)
    plt.errorbar(wphot, mphot, label='Model photometry (random draw)',
                 marker='s', markersize=10, alpha=0.8, ls='', lw=3, 
                 markerfacecolor='none', markeredgecolor='blue', 
                 markeredgewidth=3)
    plt.errorbar(wphot, mphot_map, label='Model photometry (MAP)',
                 marker='s', markersize=10, alpha=0.8, ls='', lw=3, 
                 markerfacecolor='none', markeredgecolor='green', 
                 markeredgewidth=3)
    plt.errorbar(wphot, obs['maggies'], yerr=obs['maggies_unc'], 
                 label='Observed photometry', ecolor='red', 
                 marker='o', markersize=10, ls='', lw=3, alpha=0.8, 
                 markerfacecolor='none', markeredgecolor='red', 
                 markeredgewidth=3)
    
    # plot transmission curves
    for f in obs['filters']:
        w, t = f.wavelength.copy(), f.transmission.copy()
        t = t / t.max()
        t = 10**(0.2*(np.log10(ymax/ymin)))*t * ymin
        plt.loglog(w, t, lw=3, color='gray', alpha=0.7)
    
    plt.xlabel(r'Wavelength ($\rm \AA$)')
    plt.ylabel('Flux Density (maggies)')
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.legend(loc='upper right', fontsize=20)
    plt.tight_layout()
'''