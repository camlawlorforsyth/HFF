
import numpy as np
import prospect.io.read_results as reader

import plotting as plt

def show_results(cluster, ID, binNum, results_type='dynesty') :
    infile = '{}/h5/ID_{}_bin_{}.h5'.format(cluster, ID, binNum)
    
    result, obs, _ = reader.results_from(infile, dangerous=True)
    model = reader.get_model(result)
    sps = reader.get_sps(result)
    
    '''
    set of parameters of the sample with the highest posterior probability,
    with discussion here: https://github.com/bd-j/prospector/issues/215
    
    theta_best == theta_max == theta_map in the documentation, where
    theta_max = result['chain'][np.argmax(result['lnprobability']), :]
    
    also stored in result['bestfit']['parameter']
    '''
    
    # t_hr = result['sampling_duration']/3600
    # print('The sampling took {:.2f} hours'.format(t_hr))
    
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
    # plt.plot_sed_from_fit(obs, model, result, sps, infile)
    
    return

# show_results('a370', 3337, 0) # good example of metallicity, dust2, tage values jumping around
# show_results('a370', 3337, 1)
# show_results('a370', 3337, 2)
# show_results('a370', 3337, 3)
# show_results('a370', 3337, 4)
# show_results('a370', 3337, 5)
# show_results('a370', 3337, 6)

# show_results('a1063', 1366, 0)
# show_results('a1063', 1366, 1)
# show_results('a1063', 1366, 2)
# show_results('a1063', 1366, 3)
# show_results('a1063', 1366, 4)
# show_results('a1063', 1366, 5)

# show_results('a1063', 2455, 0) # good example of tage and tau values jumping around
# show_results('a1063', 2455, 1)
# show_results('a1063', 2455, 2)
# show_results('a1063', 2455, 3)
# show_results('a1063', 2455, 4)
# show_results('a1063', 2455, 5)
# show_results('a1063', 2455, 6)
# show_results('a1063', 2455, 7)
# show_results('a1063', 2455, 8)

# show_results('a1063', 3550, 0)
# show_results('a1063', 3550, 1)
# show_results('a1063', 3550, 2)
# show_results('a1063', 3550, 3)
# show_results('a1063', 3550, 4)

# show_results('a1063', 4823, 0)
# show_results('a1063', 4823, 1)
# show_results('a1063', 4823, 2)
# show_results('a1063', 4823, 3)
# show_results('a1063', 4823, 4)
# show_results('a1063', 4823, 5)

# show_results('a2744', 3859, 0)
# show_results('a2744', 3859, 1)
# show_results('a2744', 3859, 2)
# show_results('a2744', 3859, 3)
# show_results('a2744', 3859, 4)
# show_results('a2744', 3859, 5)
# show_results('a2744', 3859, 6)

# show_results('a2744', 3964, 0)
# show_results('a2744', 3964, 1)
# show_results('a2744', 3964, 2)
# show_results('a2744', 3964, 3)
# show_results('a2744', 3964, 4)
# show_results('a2744', 3964, 5)
# show_results('a2744', 3964, 6)
# show_results('a2744', 3964, 7)
# show_results('a2744', 3964, 8)
# show_results('a2744', 3964, 9)

# show_results('a2744', 4173, 0)
# show_results('a2744', 4173, 1)
# show_results('a2744', 4173, 2)
# show_results('a2744', 4173, 3)
# show_results('a2744', 4173, 4)
# show_results('a2744', 4173, 5)
# show_results('a2744', 4173, 6)
# show_results('a2744', 4173, 7)
# show_results('a2744', 4173, 8)
# show_results('a2744', 4173, 9)

# show_results('a2744', 4369, 0)
# show_results('a2744', 4369, 1)
# show_results('a2744', 4369, 2)
# show_results('a2744', 4369, 3)
# show_results('a2744', 4369, 4)
# show_results('a2744', 4369, 5)
# show_results('a2744', 4369, 6)
# show_results('a2744', 4369, 7)
# show_results('a2744', 4369, 8)

# show_results('a2744', 4765, 0)
# show_results('a2744', 4765, 1)
# show_results('a2744', 4765, 2)
# show_results('a2744', 4765, 3)
# show_results('a2744', 4765, 4)
# show_results('a2744', 4765, 5)
# show_results('a2744', 4765, 6)
# show_results('a2744', 4765, 7)
# show_results('a2744', 4765, 8)
# show_results('a2744', 4765, 9)
# show_results('a2744', 4765, 10)
# show_results('a2744', 4765, 11)
# show_results('a2744', 4765, 12)
# show_results('a2744', 4765, 13)

# show_results('a2744', 4862, 0)
# show_results('a2744', 4862, 1)
# show_results('a2744', 4862, 2)
# show_results('a2744', 4862, 3)
# show_results('a2744', 4862, 4)
# show_results('a2744', 4862, 5)
# show_results('a2744', 4862, 6)
# show_results('a2744', 4862, 7)
# show_results('a2744', 4862, 8)
# show_results('a2744', 4862, 9)
# show_results('a2744', 4862, 10)

# show_results('a2744', 7427, 0)
# show_results('a2744', 7427, 1)
# show_results('a2744', 7427, 2)
# show_results('a2744', 7427, 3)
# show_results('a2744', 7427, 4)
# show_results('a2744', 7427, 5)
# show_results('a2744', 7427, 6)
# show_results('a2744', 7427, 7)
# show_results('a2744', 7427, 8)
# show_results('a2744', 7427, 9)

# show_results('m416', 5997, 0)
# show_results('m416', 5997, 1)
# show_results('m416', 5997, 2)
# show_results('m416', 5997, 3)
# show_results('m416', 5997, 4)
# show_results('m416', 5997, 5)
# show_results('m416', 5997, 6)
# show_results('m416', 5997, 7)
# show_results('m416', 5997, 8)
# show_results('m416', 5997, 9)
# show_results('m416', 5997, 10)

# show_results('m416', 6255, 0)
# show_results('m416', 6255, 1)
# show_results('m416', 6255, 2)
# show_results('m416', 6255, 3)
# show_results('m416', 6255, 4)
# show_results('m416', 6255, 5)
# show_results('m416', 6255, 6)
# show_results('m416', 6255, 7)
# show_results('m416', 6255, 8)
# show_results('m416', 6255, 9)
# show_results('m416', 6255, 10)

# show_results('m717', 861, 0)
# show_results('m717', 861, 1)
# show_results('m717', 861, 2)
# show_results('m717', 861, 3)
# show_results('m717', 861, 4)
# show_results('m717', 861, 5)
# show_results('m717', 861, 6)
# show_results('m717', 861, 7)

# show_results('m1149', 1967, 0)
# show_results('m1149', 1967, 1)
# show_results('m1149', 1967, 2)
# show_results('m1149', 1967, 3)
# show_results('m1149', 1967, 4)
# show_results('m1149', 1967, 5)
# show_results('m1149', 1967, 6)
# show_results('m1149', 1967, 7)

# show_results('m1149', 2403, 0)
# show_results('m1149', 2403, 1)
# show_results('m1149', 2403, 2)
# show_results('m1149', 2403, 3)
# show_results('m1149', 2403, 4)
# show_results('m1149', 2403, 5)

# show_results('m1149', 3531, 0)
# show_results('m1149', 3531, 1)
# show_results('m1149', 3531, 2)
# show_results('m1149', 3531, 3)
# show_results('m1149', 3531, 4)

# show_results('m1149', 4246, 0)
# show_results('m1149', 4246, 1)
# show_results('m1149', 4246, 2)
# show_results('m1149', 4246, 3)
# show_results('m1149', 4246, 4)
# show_results('m1149', 4246, 5)
# show_results('m1149', 4246, 6)
# show_results('m1149', 4246, 7)

# show_results('m1149', 5095, 0)
# show_results('m1149', 5095, 1)
# show_results('m1149', 5095, 2)
# show_results('m1149', 5095, 3)
# show_results('m1149', 5095, 4)
