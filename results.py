
import os
import numpy as np
from astropy.table import Table

import prospect.io.read_results as reader

import plotting as plt

def show_results(cluster, ID, binNum, results_type='dynesty') :
    infile = '{}/h5/ID_{}_bin_{}.h5'.format(cluster, ID, binNum)

    result, obs, _ = reader.results_from(infile, dangerous=True)
    model = reader.get_model(result)
    sps = reader.get_sps(result)

    # print('The sampling took {:.2f} hours'.format(result['sampling_duration']/3600))
    
    # plot SED with only photometry
    # plt.plot_sed_alone(obs)
    
    # plot SED with initial theta
    # theta = model.theta.copy()
    # initial_spec, initial_phot, _ = model.sed(theta, obs=obs, sps=sps)
    # plt.plot_sed_and_model(obs, model, initial_phot, initial_spec, sps)
    
    # investigate the parameter traces
    # plt.plot_chains(result, results_type=results_type)
    
    # create corner plots
    # plt.plot_corner(model, result, result['bestfit']['parameter'])
    
    # plot SED and residuals
    '''
    if os.path.isfile('one_sigma_arrays.npz') :
        npzfile = np.load('one_sigma_arrays.npz')
        sixteen = npzfile['sixteen']
        eightyfour = npzfile['eightyfour']
    else :
        if os.path.isfile('table_of_sampled_specs.fits') :
            sixteen = []
            eightyfour = []
            table = Table.read('table_of_sampled_specs.fits')
            for i in range(len(table)) :
                row = list(table[i].as_void())
                lo = np.percentile(row, 16)
                hi = np.percentile(row, 84)
                sixteen.append(lo)
                eightyfour.append(hi)
        
            sixteen = np.array(sixteen)
            eightyfour = np.array(eightyfour)
            
            np.savez('one_sigma_arrays.npz', sixteen=sixteen, eightyfour=eightyfour)
        else :
            samples = pplots.utils.sample_posterior(result['chain'],
                                                    weights=result['weights'],
                                                    nsample=995)
            table = Table()
            for i in range(len(samples)) :
                spec, phot, _ = model.predict(samples[i], obs=obs, sps=sps)
                table['spec_{}'.format(i)] = spec
            table.write('table_of_sampled_specs.fits')
    '''
    plt.plot_sed_from_fit(obs, model, result, sps, infile)
    
    return

show_results('a370', 3337, 0)
