
import os
import glob
import numpy as np

from astropy.table import Table
import astropy.units as u
import prospect.io.read_results as reader
from prospect.plotting.utils import sample_posterior

import plotting as plt

def extract_posteriors() :
    
    logMs, logZs, dusts, t0s, logTaus = [], [], [], [], []
    clusters = ['a370', 'a1063', 'a2744', 'm416', 'm717', 'm1149']
    for cluster in clusters :
        
        h5_files = glob.glob('{}/h5/*.h5'.format(cluster))
        
        h5s = []
        for file in h5_files :
            file = file.replace(os.sep, '/') # compatibility for Windows
            h5s.append(file)
        
        for file in h5s :
            result, obs, _ = reader.results_from(file, dangerous=True)
            
            samples = sample_posterior(result['chain'],
                                       weights=result['weights'],
                                       nsample=len(result['chain']))
            
            samples[:, 1] = np.log10(samples[:, 1])
            samples[:, 5] = np.log10(samples[:, 5])
            
            logMs.append(list(samples[:, 1]))
            logZs.append(list(samples[:, 2]))
            dusts.append(list(samples[:, 3]))
            t0s.append(list(samples[:, 4]))
            logTaus.append(list(samples[:, 5]))
    
    logM = flatten_posteriors(logMs)
    logZ = flatten_posteriors(logZs)
    dust = flatten_posteriors(dusts)
    t0 = flatten_posteriors(t0s)
    logTau = flatten_posteriors(logTaus)
    
    # np.savez('subsample_posteriors.npz', logM=logM, logZ=logZ, dust=dust,
    #          t0=t0, logTau=tau)
    
    solMetal = u.def_unit(['solMetal', 'Z_sun', 'Zsun'], prefixes=False,
                          format={'latex':r'Z_{\odot}', 'unicode': 'Z\N{SUN}'})
    
    table = Table([logM/u.solMass, logZ/solMetal, dust, t0*u.Gyr, logTau/u.Gyr],
                  names=('logM', 'logZ', 'dust', 't0', 'logTau'))
    table.write('boneyard/subsample_posteriors.fits')
    
    return

def flatten_posteriors(array) :
    
    length = max(map(len, array)) # determine maximum length
    
    padded_row = [row + [None]*(length - len(row)) for row in array]
    final = (np.array(padded_row, dtype=object)).flatten()
    
    return np.float64(final[final != None])

def plot_posteriors(param, minmax=False, percentiles=False) :
    
    table = Table.read('boneyard/subsample_posteriors.fits')
    val = table[param]
    
    plt.histogram(val, '{}'.format(param), bins=50, histtype='step')
    
    if percentiles :
        print(np.percentile(val, [0.15, 50, 99.85]))
    
    if minmax :
        print(np.min(val), np.percentile(val, 50), np.max(val))
    
    return

# plot_posteriors('logM')
# plot_posteriors('logZ')
# plot_posteriors('dust')
# plot_posteriors('t0')
# plot_posteriors('logTau')
