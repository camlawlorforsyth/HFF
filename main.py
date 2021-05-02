
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
    
    os.makedirs('{}'.format(cluster), exist_ok=True) # ensure the output
        # directory is available
    if verbose :
        print('Created cluster directory.')
    
    # determine the final science objects, based on flags in the catalogs
    core.determine_finalObjs_w_color(cluster, redshift, path1, path2,
                                     redshift_tol_lo=delta_z_lo,
                                     redshift_tol_hi=delta_z_hi,
                                     z_spec=force_spec,
                                     plot_all=False, plot_uvj=False,
                                     write_final_objs=False,
                                     write_regions=False, selection='FUVVJ')
    if verbose :
        print('Determined final objects and saved to fits file.')
    
    # save a file containing information about the filterset for the cluster
    field.save_filterset(cluster, filters)
    if verbose :
        print('Saved filterset to file.')
    
    # determine the rms values for each filter
    if calculate_rms :
        rms = field.determine_rms(segPath, files)
        print(rms)
    else :
        rms = RMS
    if verbose :
        print('RMS determination complete.')
    
    # save all the cutouts for all the science objects
    final_objs_path = '{}/{}_final_objects.fits'.format(cluster, cluster)
    selection = 'Q' # only consider quiescent galaxies
    field.save_cutouts(cluster, final_objs_path, filters, rms, files, segPath,
                       selection, models, z_spec=force_spec)
    if verbose :
        print('Saved all science, noise, and segmentation map cutouts.')
    
    # now bin the cutouts for each galaxy and save the resulting files
    if vorbin :
        voronoi.vorbin_all(cluster)
    else :
        binning.bin_all(cluster)
    if verbose :
        print('Determined all bins and saved numpy arrays to files.')
    
    # then determine the flux for every bin for each galaxy, saving the
    # photometry for each galaxy into a separate file
    photometry.determine_fluxes(cluster, filters)
    if verbose :
        print('Saved photometries to fits files.')
    
    # check how many annuli there are per galaxy
    checks.check_bins(cluster)
    if verbose :
        print('Plotted annuli histogram.')
    
    # END
    
    # then on linux, move to the `hff/` directory, and run
    # python prospector.py
    
    return

def postmain(total_FUVVJ=False, total_UVJ=False, parallel_objects=False) :
    
    if total_FUVVJ :
        checks.plot_all_FUVVJ()
    
    if total_UVJ :
        checks.plot_all_UVJ()
    
    if parallel_objects :
        checks.plot_all_parallel_objects_all()
    
    return
