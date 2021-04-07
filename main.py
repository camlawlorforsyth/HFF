
import binning
import core
import field
import photometry
import plotting as plt
import resave
import voronoi

def premain(save_catalogs=False, save_filters=False, plot_curves=False) :
    '''
    
    
    Parameters
    ----------
    save_catalogs : TYPE, optional
        DESCRIPTION. The default is False.
    save_filters : TYPE, optional
        DESCRIPTION. The default is False.
    plot_curves : TYPE, optional
        DESCRIPTION. The default is False.
    
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

def main(cluster, calculate_rms=False, binning_scheme='annuli') :
    '''
    
    
    Parameters
    ----------
    cluster : TYPE
        DESCRIPTION.
    calculate_rms : bool
        Flag to set if the rms values should be calculated. Default is False.
    binning_scheme : TYPE, optional
        DESCRIPTION. The default is 'annuli'.
    
    Returns
    -------
    None.
    
    '''
    
    U_filtnum, V_filtnum, J_filtnum = 153, 155, 161
    
    if cluster == 'a370' :
        cluster_redshift = 0.375
        first_path, second_path = 'abell370clu_catalogs', 'abell370clu_v3.9'
        misc = 'misc/abell370clu_misc/'
        
        filters = ['f275w', 'f336w', # UVIS
                   'f435w', 'f475w', 'f606w', 'f625w', 'f814w', # ACS
                   'f105w', 'f110w', 'f125w', 'f140w', 'f160w'] # IR
        
        segPath = misc + 'photometry/abell370_bcgs_out_detection_seg.fits.gz'
        
        psf_matched = misc + 'images/psf_matched/'
        bcgs_out = misc + 'images/bcgs_out/'
        paths = [psf_matched + 'abell370_bcgs_out_f275w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell370_bcgs_out_f336w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell370_bcgs_out_f435w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell370_bcgs_out_f606w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell370_bcgs_out_f814w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell370_bcgs_out_f105w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell370_bcgs_out_f125w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell370_bcgs_out_f140w_psf_bkg_drz.fits.gz',
                 bcgs_out + 'abell370_bcgs_out_f160w_bkg_drz_masked.fits.gz']
    
    if cluster == 'a1063' :
        cluster_redshift = 0.348
        first_path, second_path = 'abell1063clu_catalogs', 'abell1063clu_v3.9'
        misc = 'misc/abell1063clu_misc/'
        
        filters = ['f225w', 'f275w', 'f336w','f390w', # UVIS
                   'f435w', 'f475w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp', # ACS
                   'f105w', 'f110w', 'f125w', 'f140w', 'f160w'] # IR
        
        segPath = misc + 'photometry/abell1063_bcgs_out_detection_seg.fits.gz'
        
        psf_matched = misc + 'images/psf_matched/'
        bcgs_out = misc + 'images/bcgs_out/'
        paths = [psf_matched + 'abell1063_bcgs_out_f275w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell1063_bcgs_out_f336w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell1063_bcgs_out_f435w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell1063_bcgs_out_f606w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell1063_bcgs_out_f814w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell1063_bcgs_out_f105w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell1063_bcgs_out_f125w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell1063_bcgs_out_f140w_psf_bkg_drz.fits.gz',
                 bcgs_out + 'abell1063_bcgs_out_f160w_bkg_drz_masked.fits.gz']
    
    if cluster == 'a2744' :
        cluster_redshift = 0.308
        first_path, second_path = 'abell2744clu_catalogs', 'abell2744clu_v3.9'
        misc = 'misc/abell2744clu_misc/'
        
        filters = ['f275w', 'f336w', # UVIS
                   'f435w', 'f606w', 'f814w', # ACS
                   'f105w', 'f125w', 'f140w', 'f160w'] # IR
        
        segPath = misc + 'photometry/abell2744_bcgs_out_detection_seg.fits.gz'
        
        psf_matched = misc + 'images/psf_matched/'
        bcgs_out = misc + 'images/bcgs_out/'
        paths = [psf_matched + 'abell2744_bcgs_out_f275w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell2744_bcgs_out_f336w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell2744_bcgs_out_f435w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell2744_bcgs_out_f606w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell2744_bcgs_out_f814w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell2744_bcgs_out_f105w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell2744_bcgs_out_f125w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell2744_bcgs_out_f140w_psf_bkg_drz.fits.gz',
                 bcgs_out + 'abell2744_bcgs_out_f160w_bkg_drz_masked.fits.gz']
        
        rms = [0.00017821459857111745, 0.00014789605532426083,
               0.0007145838788770763,  0.001578221285547029,
               0.0018031112928780836,  0.0293976362500103,
               0.002737900021748479,   0.016844529370723973,
               0.0037151459]
    
    if cluster == 'm416' :
        cluster_redshift = 0.392
        first_path, second_path = 'macs0416clu_catalogs', 'macs0416clu_v3.9'
        misc = 'misc/macs0416clu_misc/'
        
        filters = ['f225w', 'f275w', 'f336w','f390w', # UVIS
                   'f435w', 'f475w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp', # ACS
                   'f105w', 'f110w', 'f125w', 'f140w', 'f160w'] # IR
        
        segPath = misc + 'photometry/macs0416_bcgs_out_detection_seg.fits.gz'
        
        psf_matched = misc + 'images/psf_matched/'
        bcgs_out = misc + 'images/bcgs_out/'
        paths = [psf_matched + 'macs0416_bcgs_out_f275w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0416_bcgs_out_f336w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0416_bcgs_out_f435w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0416_bcgs_out_f606w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0416_bcgs_out_f814w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0416_bcgs_out_f105w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0416_bcgs_out_f125w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0416_bcgs_out_f140w_psf_bkg_drz.fits.gz',
                 bcgs_out + 'macs0416_bcgs_out_f160w_bkg_drz_masked.fits.gz']
    
    if cluster == 'm717' :
        cluster_redshift = 0.5458
        first_path, second_path = 'macs0717clu_catalogs', 'macs0717clu_v3.9'
        misc = 'misc/macs0717clu_misc/'
        
        filters = ['f225w', 'f275w', 'f336w','f390w', # UVIS
                   'f435w', 'f475w', 'f555w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp', # ACS
                   'f105w', 'f110w', 'f125w', 'f140w', 'f160w'] # IR
        
        segPath = misc + 'photometry/macs0717_bcgs_out_detection_seg.fits.gz'
        
        psf_matched = misc + 'images/psf_matched/'
        bcgs_out = misc + 'images/bcgs_out/'
        paths = [psf_matched + 'macs0717_bcgs_out_f275w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0717_bcgs_out_f336w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0717_bcgs_out_f435w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0717_bcgs_out_f606w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0717_bcgs_out_f814w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0717_bcgs_out_f105w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0717_bcgs_out_f125w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0717_bcgs_out_f140w_psf_bkg_drz.fits.gz',
                 bcgs_out + 'macs0717_bcgs_out_f160w_bkg_drz_masked.fits.gz']
    
    if cluster == 'm1149' :
        cluster_redshift = 0.543
        first_path, second_path = 'macs1149clu_catalogs', 'macs1149clu_v3.9'
        misc = 'misc/macs1149clu_misc/'
        
        filters = ['f225w', 'f275w', 'f336w','f390w', # UVIS
                   'f435w', 'f475w', 'f555w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp', # ACS
                   'f105w', 'f110w', 'f125w', 'f140w', 'f160w'] # IR
        
        segPath = misc + 'photometry/macs1149_bcgs_out_detection_seg.fits.gz'
        
        psf_matched = misc + 'images/psf_matched/'
        bcgs_out = misc + 'images/bcgs_out/'
        paths = [psf_matched + 'macs1149_bcgs_out_f275w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs1149_bcgs_out_f336w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs1149_bcgs_out_f435w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs1149_bcgs_out_f606w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs1149_bcgs_out_f814w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs1149_bcgs_out_f105w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs1149_bcgs_out_f125w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs1149_bcgs_out_f140w_psf_bkg_drz.fits.gz',
                 bcgs_out + 'macs1149_bcgs_out_f160w_bkg_drz_masked.fits.gz']
    
    # determine the final science objects, based on flags in the catalogs
    core.determine_finalObjs_w_UVJ(cluster, 'id', cluster_redshift,
                                   first_path, second_path,
                                   U_filtnum, V_filtnum, J_filtnum,
                                   plot_all=True, plot_uvj=True,
                                   write_final_objs=True, write_regions=True)
    
    # save a file containing information about the filterset for the cluster
    field.save_filterset(cluster, filters)
    
    # determine the rms values for each filter for Abell 2744
    if calculate_rms :
        rms = field.determine_rms(segPath, paths)
    else :
        rms = rms
    
    # save all the cutouts for all the science objects
    final_objs_path = '{}/{}_final_objects.fits'.format(cluster, cluster)
    field.save_cutouts(cluster, final_objs_path, filters, rms, paths, segPath)
    
    # now vorbin the cutouts for each galaxy and save the resulting files
    if binning_scheme == 'annuli' :
        binning.bin_all(cluster)
    else :
        voronoi.vorbin_all(cluster)
    
    # then determine the flux for every vorbin for each galaxy, saving the
    # photometry for each galaxy into a separate file
    photometry.determine_fluxes(cluster, filters)
    
    # then on linux, move to the `hff/` directory, and run
    # python prospector.py
    
    return
