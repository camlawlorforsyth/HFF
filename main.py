
import os

import binning
import core
import checks
import field
import initial
import photometry
import plotting as plt
import resave
import voronoi

def premain(save_catalogs=False, save_filters=False, plot_curves=False) :
    '''
    Operations to complete before the main analysis operations.
    
    Parameters
    ----------
    save_catalogs : bool, optional
        Flag to resave catalogs in fits formats. The default is False.
    save_filters : bool, optional
        Flag to save HFF filters in sedpy format. The default is False.
    plot_curves : bool, optional
        Flag to plot HFF filterset transmission curves. The default is False.
    
    Returns
    -------
    None.
    
    '''
    
    if save_catalogs :
        resave.resave_all()
    
    if save_filters :
        core.build_filter_table(save_par=True)
    
    if plot_curves :
        plt.hst_transmission_curves(core.build_filter_table())
    
    return

def main(cluster, calculate_rms=False, verbose=False, vorbin=False) :
    '''
    The main (preparatory) analysis operations. Determine the final science
    objects for each cluster and parallel field, save all relevant cutouts,
    determine bins (annuli) for each galaxy, and determine fluxes for each bin,
    saving the output into a Prospector-compatible file for subsequent use.
    
    Parameters
    ----------
    cluster : string
        The cluster to operate on.
    calculate_rms : bool, optional
        Flag to calculate the RMS or use existing values. The default is False.
    verbose : bool, optional
        Flag to have messages printed after milestones. The default is False.
    vorbin : bool, optional
        Flag to set binning method. The default is False.
    
    Returns
    -------
    None.
    
    '''
    
    # CLUSTERS
    
    if cluster == 'a370' :
        (redshift, sigma, delta_z_lo, delta_z_hi, redshift_type, path1, path2,
         filters, segPath, bCGsegPath, files, models, RMS) = initial.a370_params()
    
    if cluster == 'a1063' :
        (redshift, sigma, delta_z_lo, delta_z_hi, redshift_type, path1, path2,
         filters, segPath, bCGsegPath, files, models, RMS) = initial.a1063_params()
    
    if cluster == 'a2744' :
        (redshift, sigma, delta_z_lo, delta_z_hi, redshift_type, path1, path2,
         filters, segPath, bCGsegPath, files, models, RMS) = initial.a2744_params()
    
    if cluster == 'm416' :
        (redshift, sigma, delta_z_lo, delta_z_hi, redshift_type, path1, path2,
         filters, segPath, bCGsegPath, files, models, RMS) = initial.m416_params()
    
    if cluster == 'm717' :
        (redshift, sigma, delta_z_lo, delta_z_hi, redshift_type, path1, path2,
         filters, segPath, bCGsegPath, files, models, RMS) = initial.m717_params()
    
    if cluster == 'm1149' :
        (redshift, sigma, delta_z_lo, delta_z_hi, redshift_type, path1, path2,
         filters, segPath, bCGsegPath, files, models, RMS) = initial.m1149_params()
    
    # PARALLEL FIELDS
    
    if cluster == 'a370par' :
        (redshift, sigma, delta_z_lo, delta_z_hi, redshift_type, path1, path2,
         filters, segPath, bCGsegPath, files, models, RMS) = initial.a370par_params()
    
    if cluster == 'a1063par' :
        (redshift, sigma, delta_z_lo, delta_z_hi, redshift_type, path1, path2,
         filters, segPath, bCGsegPath, files, models, RMS) = initial.a1063par_params()
    
    if cluster == 'a2744par' :
        (redshift, sigma, delta_z_lo, delta_z_hi, redshift_type, path1, path2,
         filters, segPath, bCGsegPath, files, models, RMS) = initial.a2744par_params()
    
    if cluster == 'm416par' :
        (redshift, sigma, delta_z_lo, delta_z_hi, redshift_type, path1, path2,
         filters, segPath, bCGsegPath, files, models, RMS) = initial.m416par_params()
    
    if cluster == 'm717par' :
        (redshift, sigma, delta_z_lo, delta_z_hi, redshift_type, path1, path2,
         filters, segPath, bCGsegPath, files, models, RMS) = initial.m717par_params()
    
    if cluster == 'm1149par' :
        (redshift, sigma, delta_z_lo, delta_z_hi, redshift_type, path1, path2,
         filters, segPath, bCGsegPath, files, models, RMS) = initial.m1149par_params()
    
    # START
    
    # ensure the output directory is available
    os.makedirs('{}'.format(cluster), exist_ok=True) 
    if verbose :
        print('Created: Cluster directory.')
    
    # determine the final science objects, based on flags in the catalogs
    core.determine_final_sample(cluster, redshift, sigma, path1, path2,
                                redshift_tol_lo=delta_z_lo,
                                redshift_tol_hi=delta_z_hi,
                                redshift_type='z_spec',
                                plot_all=False, plot_uvj=False,
                                write_final_objs=True,
                                write_regions=True, selection='FUVVJ',
                                verbose=verbose)
    if verbose :
        print('Determined and saved: Sample and population demographics to file.')
    
    # determine the rms values for each filter
    if calculate_rms :
        rms = field.determine_rms(segPath, files)
        print(rms)
    else :
        rms = RMS
    if verbose :
        print('Determined: Background RMS.')
    
    # save all the cutouts for all the science objects
    sample_path = '{}/{}_sample.fits'.format(cluster, cluster)
    field.save_cutouts(cluster, sample_path, filters, rms, files, segPath,
                       bCGsegPath, models, redshift_type=redshift_type)
    if verbose :
        print('Saved: All science, noise, and segmentation map cutouts.')
    
    # save pngs of the cutouts for visual inspection
    checks.save_pngs(cluster, filters)
    if verbose :
        print('Saved: Images of cutouts for visual inspection.')
    
    # now visually inspect all the cutouts to ensure the final objects are
    # reasonable and have no glaring issues, saving issues to
    # `cluster`_issues.csv. Then combine the issues file with the objects file
    field.combine_with_issues(cluster)
    if verbose :
        print('Complete: Issues combined with science objects.')
    
    # now bin the cutouts for each galaxy and save the resulting files
    if vorbin :
        voronoi.vorbin_all(cluster)
    else :
        binning.bin_all(cluster)
    if verbose :
        print('Saved: All bins to numpy array files.')
    
    # then determine the flux for every bin for each galaxy, saving the
    # photometry for each galaxy into a separate file
    photometry.determine_fluxes(cluster, filters)
    if verbose :
        print('Saved: Photometries to fits files.')
    
    # finally, check for galaxies that don't have any bins or only a single bin,
    # saving an updated sample file for galaxies that have at least two bins
    field.check_number_of_bins(cluster)
    if verbose :
        print('Saved: Updated sample file.')
    
    # END
    
    '''
    Note: the sample as determined above from the `determine_final_sample`
    function returns: 1) 668 quiescent galaxies
                      2) 179 star forming galaxies
                      3) 733 cluster galaxies
                      4) 114 field galaxies
    for a total of 847 galaxies.
    
    Of these, 32 should be removed in the step performed by the
    `combine_with_issues` function, as these galaxies all have incomplete or
    poor F160W coverage, and so reliable bins/annuli cannot be determined.
    These galaxies are:
    a370 ID 20001
    a1063 ID 20005
    a1063 ID 20016
    a1063 ID 20023
    a1063 ID 20024
    a1063 ID 20025
    a1063 ID 20056
    a1063 ID 20058
    a1063 ID 20059
    a1063 ID 20060 (not detected in segmentation map)
    a2744 ID 20002
    a2744 ID 20005
    a2744 ID 20008
    a2744 ID 20030
    a2744 ID 20036
    m416 ID 20001
    m416 ID 20002
    m416 ID 20003
    m416 ID 20024
    m717 ID 905
    m717 ID 3896
    m717 ID 4819
    m717 ID 20004
    m717 ID 20005
    m1149 ID 212
    m1149 ID 1678
    m1149 ID 2737
    m1149 ID 20002
    m1149 ID 20003
    m1149 ID 20004
    a2744par ID 4502
    m1149par ID 5208
    
    Following this, 19 should be further removed in the step performed by the
    `photometry.determine_fluxes` function, as these galaxies are all too dim
    to have even a single bin/annulus determined (but note that they will have
    *_annuli.npz files, but will not have information populated into the
    bins_image arrays included in those files). These galaxies are:
    a370 ID 2205
    a1063 ID 3030
    a1063 ID 3089
    a1063 ID 3426
    a1063 ID 4556
    a1063 ID 4746
    a1063 ID 4751
    a1063 ID 4755
    a1063 ID 4761
    a1063 ID 5090
    a1063 ID 5160
    a1063 ID 5188
    a1063 ID 5397
    a1063 ID 5638
    a1063 ID 5705
    a1063 ID 5806
    a1063 ID 6178
    a2744 ID 4358
    m1149 ID 4316
    
    Lastly, 19 should be finally removed manually, as these galaxies all have
    only a single bin/annulus. These galaxies are:
    a1063 ID 2601
    a1063 ID 2959
    a1063 ID 2986
    a1063 ID 3221
    a1063 ID 3312
    a1063 ID 3486
    a1063 ID 4213
    a1063 ID 4377
    a1063 ID 4475
    a1063 ID 5083
    a1063 ID 6093
    a2744 ID 2226
    a2744 ID 3812
    a2744 ID 3892
    a2744 ID 4655
    a2744 ID 4768
    a2744 ID 7384
    m1149 ID 4247
    m416par ID 4628
    
    This results in: 1) 617 quiescent galaxies
                     2) 160 star forming galaxies
                     3) 672 cluster galaxies
                     4) 105 field galaxies
    for a total of 777 galaxies in the final sample.
    '''
    
    return

def postmain(total_FUVVJ=False, total_UVJ=False, parallel_objects=False,
             dists=False, bins=False) :
    '''
    Operations to complete after the main analysis operations.
    
    Parameters
    ----------
    total_FUVVJ : bool, optional
        Flag to plot the FUV-V-J diagram for all objects. The default is False.
    total_UVJ : bool, optional
        Flag to plot the UVJ diagram for all objects. The default is False.
    parallel_objects : bool, optional
        Flag to plot all parallel objects simultaneously. The default is False.
    dists : bool, optional
        Flag to plot redshift and logmass distributions. The default is False.
    bins : bool, optional
        Flag to plot distributions of the number of bins. The default is False.
    
    Returns
    -------
    None.
    
    '''
    
    if concat :
        checks.concatenate()
    
    if total_FUVVJ :
        checks.plot_all_FUVVJ()
    
    if total_UVJ :
        checks.plot_all_UVJ()
    
    if parallel_objects :
        checks.plot_parallel_objects_all()
    
    if dists :
        checks.check_distributions()
    
    if bins :
        checks.check_all_bins()
    
    # plot and save the surface brightness profiles for every filter, using
    # the photometry file saved above
    checks.save_sbps(cluster, population)
    if verbose :
        print('Saved: Surface brightness profiles to file.')
    
    # plot and save the SEDs for every annulus/radial bin, using the
    # photometry file saved above
    checks.save_seds(cluster, population)
    if verbose :
        print('Saved: Spectral energy distributions to file.')
    
    # plot and save histograms of the background pixels to verify that the
    # background subtraction was completed correctly
    checks.save_bkgshists(cluster, filters, population)
    if verbose :
        print('Saved: Background pixel histograms to file.')
    
    return

import warnings
warnings.filterwarnings('ignore')

# main('a370')
# main('a1063')
# main('a2744')
# main('m416')
# main('m717')
# main('m1149')

# main('a370par') # 0 galaxies in the sample for this parallel field
# main('a1063par') # 0 galaxies in the sample for this parallel field
# main('a2744par')
# main('m416par')
# main('m717par')
# main('m1149par')

# import numpy as np
# nfilters = np.array([12, 16, 9, 16, 17, 17])
# ngals = np.array([74, 103, 139, 92, 110, 129]) # quiescent galaxies
# ncutouts = (2*nfilters + 1)*ngals # noise and science, and one segmap

# cc = 299792.458
# zs = np.array([0.375, 0.348, 0.308, 0.396, 0.545, 0.543])
# print(zs)
# sigmas = np.array([1170, 1840, 1497, 955, 1660, 1840])
# deltas = (3*sigmas/cc)*(1+zs)
# print(deltas)
# Ha = 6562.8
# print(Ha*(1+zs))

# from astropy.cosmology import FlatLambdaCDM
# cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
# hi = cosmo.lookback_time(np.inf)
# lims = np.array([0, 0.03, 0.1, 0.5, 1, 3.3587269499417567,
#                  5.717453899883513, 8.07618084982527, 10.434907799767027,
#                  12.793634749708783, 13.466983947061877])*1e9
# agelims = np.log10(lims)
# agelims[0] = 0
# agebins = np.array([agelims[:-1], agelims[1:]]).T
