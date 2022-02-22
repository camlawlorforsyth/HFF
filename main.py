
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
        (redshift, delta_z_lo, delta_z_hi, force_spec, path1, path2,
         filters, segPath, files, models, RMS) = initial.a370_params()
    
    if cluster == 'a1063' :
        (redshift, delta_z_lo, delta_z_hi, force_spec, path1, path2,
         filters, segPath, files, models, RMS) = initial.a1063_params()
    
    if cluster == 'a2744' :
        (redshift, delta_z_lo, delta_z_hi, force_spec, path1, path2,
         filters, segPath, files, models, RMS) = initial.a2744_params()
    
    if cluster == 'm416' :
        (redshift, delta_z_lo, delta_z_hi, force_spec, path1, path2,
         filters, segPath, files, models, RMS) = initial.m416_params()
    
    if cluster == 'm717' :
        (redshift, delta_z_lo, delta_z_hi, force_spec, path1, path2,
         filters, segPath, files, models, RMS) = initial.m717_params()
    
    if cluster == 'm1149' :
        (redshift, delta_z_lo, delta_z_hi, force_spec, path1, path2,
         filters, segPath, files, models, RMS) = initial.m1149_params()
    
    # PARALLEL FIELDS
    
    if cluster == 'a370par' :
        (redshift, delta_z_lo, delta_z_hi, force_spec, path1, path2,
         filters, segPath, files, models, RMS) = initial.a370par_params()
    
    if cluster == 'a1063par' :
        (redshift, delta_z_lo, delta_z_hi, force_spec, path1, path2,
         filters, segPath, files, models, RMS) = initial.a1063par_params()
    
    if cluster == 'a2744par' :
        (redshift, delta_z_lo, delta_z_hi, force_spec, path1, path2,
         filters, segPath, files, models, RMS) = initial.a2744par_params()
    
    if cluster == 'm416par' :
        (redshift, delta_z_lo, delta_z_hi, force_spec, path1, path2,
         filters, segPath, files, models, RMS) = initial.m416par_params()
    
    if cluster == 'm717par' :
        (redshift, delta_z_lo, delta_z_hi, force_spec, path1, path2,
         filters, segPath, files, models, RMS) = initial.m717par_params()
    
    if cluster == 'm1149par' :
        (redshift, delta_z_lo, delta_z_hi, force_spec, path1, path2,
         filters, segPath, files, models, RMS) = initial.m1149par_params()
    
    # START
    
    population = 'Q' # only consider quiescent galaxies
    subpop = 'non-bCG' # only consider the non bCG quiescent galaxies
    
    # ensure the output directory is available
    os.makedirs('{}'.format(cluster), exist_ok=True) 
    if verbose :
        print('Created: Cluster directory.')
    
    # determine the final science objects, based on flags in the catalogs
    core.determine_finalObjs_w_color(cluster, redshift, path1, path2,
                                     redshift_tol_lo=delta_z_lo,
                                     redshift_tol_hi=delta_z_hi,
                                     z_spec=force_spec,
                                     plot_all=False, plot_uvj=False,
                                     write_final_objs=True,
                                     write_regions=True, selection='FUVVJ',
                                     verbose=verbose)
    
    # save a file containing information about the filterset for the cluster
    field.save_filterset(cluster, filters)
    if verbose :
        print('Saved: Filterset to file.')
    
    # determine the rms values for each filter
    if calculate_rms :
        rms = field.determine_rms(segPath, files)
        print(rms)
    else :
        rms = RMS
    if verbose :
        print('Determined: Background RMS.')
    
    # save all the cutouts for all the science objects
    final_objs_path = '{}/{}_final_objects.fits'.format(cluster, cluster)
    field.save_cutouts(cluster, final_objs_path, filters, rms, files, segPath,
                       population, models, z_spec=force_spec)
    if verbose :
        print('Saved: All science, noise, and segmentation map cutouts.')
    
    # save pngs of the cutouts for visual inspection
    checks.save_pngs(cluster, filters, population)
    if verbose :
        print('Saved: Images of cutouts for visial inspection.')
    
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
    photometry.determine_fluxes(cluster, filters, subpop)
    if verbose :
        print('Saved: Photometries to fits files.')
    
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
    
    # END
    if verbose :
        print('\n{} complete.'.format(cluster))
    
    # then on linux, move to the `hff/` directory, and run
    # python prospector.py
    
    '''
    Note: the sample as determined above from the `determine_finalObjs_w_color`
    function returns 421 non-bCG quiescent cluster galaxies.
    
    Of these, 6 should be removed in the step performed by the
    `combine_with_issues` function, as these galaxies all have incomplete
    F160W coverage, and so reliable bins/annuli cannot be determined. These
    galaxies are:
    M717 ID 3503 (14 bins from erroneously including it previously)
    M717 ID 3896 (4 bins)
    M717 ID 4819 (14 bins)
    M717 ID 5357 (14 bins)
    M1149 ID 1678 (8 bins)
    M1149 ID 2737 (10 bins)
    
    Following this, 12 should be further removed in the step performed by the
    `binning.bin_all` function, as these galaxies are all too dim to have even
    a single bin/annulus determined. These galaxies are:
    A1063 ID 3089
    A1063 ID 3426
    A1063 ID 4556
    A1063 ID 4746
    A1063 ID 4751
    A1063 ID 4755
    A1063 ID 5090
    A1063 ID 5188
    A1063 ID 5397
    A1063 ID 5638
    A1063 ID 5806
    A2744 ID 4358
    
    Lastly, 10 should be finally removed manually, as these galaxies all have
    only a single bin/annulus. These galaxies are:
    A1063 ID 2601
    A1063 ID 2959
    A1063 ID 3221
    A1063 ID 3312
    A1063 ID 3486
    A1063 ID 4377
    A1063 ID 4475
    A1063 ID 5083
    A2744 ID 3892
    A2744 ID 4655
    
    This results in 393 galaxies that should be in the final sample for the
    non-bCG quiescent cluster galaxies.
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
    
    return

import warnings
warnings.filterwarnings('ignore')

# main('a370')
# main('a1063')
# main('a2744')
# main('m416')
# main('m717')
# main('m1149')

# main('a370par')
# main('a1063par')
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

# import numpy as np
# factors = [2.181495551643033, 2.1883236657000866, 1.7261923253729665,
#            1.4134601654681689, 1.7324986472791672, 1.7114998726238,
#            1.5857568889472982, 2.6483300822937696, 2.397462563618582,
#            2.8059658333596764, 2.18455305916632, 2.077538303832598,
#            1.9652980170768732, 1.9316625789196562, 1.9315518515817243,
#            1.8138795439082092, 1.392829120841449, 1.4227200431836113,
#            1.3915003268933417, 1.4241272366253315, 1.421152644222854,
#            1.3549649789671274, 1.4192680810266645, 2.170888909291974,
#            2.45585483722647, 2.3614106449535206, 2.1700191051298745,
#            2.219562701601038, 2.1802716636763524, 2.1798622548063715,
#            1.386939736638916, 1.3962530131418174, 1.3956817157033081,
#            2.421668950737322, 2.8516138535540807, 2.7287068061679354,
#            2.713726292101576, 1.9844135844148416, 1.9386668270745775,
#            1.9389045227865715, 2.199571064145755, 1.4193668752069424,
#            1.42821163096129, 1.478674235780314, 1.4284505218185852,
#            1.4356259731907939, 1.4058119574285082, 1.428014667786318,
#            2.214759619956099, 2.535955887224781, 2.239966316587822,
#            2.2560110851031854, 2.217309670710604, 1.9680403862972147,
#            1.8939748331238828, 1.8974302895914317, 2.1879447637719953,
#            1.399918103436626, 1.4409934384603484, 1.536614544618529,
#            1.3661944057852728, 1.4320088314942003, 1.4346914603552603,
#            1.3219041910065046, 1.4388740019642372, 2.1262681548096305,
#            2.569061713530438, 2.2079024297972043, 2.12912125366805,
#            2.124302015834728, 2.129775821991798, 2.0900351323184663,
#            2.0915760347156476, 2.358438881906315, 1.663681474153013,
#            1.531419851673137, 1.6484013034043448, 1.664127294509825,
#            1.538556155560637, 1.536506163681755, 1.6299226536398579,
#            1.5359034863153715, 2.21734593712219, 2.713042736892904,
#            2.0011043985983497, 2.6713350943262917, 1.9946837285777617]
# print(len(factors), np.mean(factors), np.median(factors))

# import numpy as np
# pfactors = [1.6156726833307276, 1.6162379269626066, 1.6135370037225696,
#             2.539293197916151, 2.5454421025138343, 2.5460908723927083,
#             2.5355020590645023, 1.5958949344690876, 1.595828787314081,
#             1.6086396207497886, 2.5200011299491325, 2.1824981528326868,
#             2.5210049093384512, 2.1790724593648436, 1.6052631295505537,
#             1.5605351037936066, 1.6198546113535777, 2.58298961785788,
#             2.5698491616872534, 2.5722946900624084, 2.570783754535875,
#             1.6291833844736454, 1.6302878438361792, 1.317603582205626,
#             1.6274413374833114, 1.3176681380989899, 2.600486169697143,
#             2.590394150479253, 2.5893631111878586, 2.586784746916275,
#             1.5550798933230554, 1.5552989059082531, 1.5816774306856214,
#             2.6480808209830102, 2.289343226248264, 2.6508573986473563,
#             2.274597884284631, 1.6136376564731885, 1.4570434258537033,
#             1.459570978879453, 2.5368405540893133, 2.539279938183213,
#             2.5389961993265047, 2.4505216964809966]
# print(len(pfactors), np.mean(pfactors), np.median(pfactors))
