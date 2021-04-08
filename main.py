
import os

import binning
import core
import checks
import field
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
                 psf_matched + 'abell370_bcgs_out_f475w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell370_bcgs_out_f606w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell370_bcgs_out_f625w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell370_bcgs_out_f814w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell370_bcgs_out_f105w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell370_bcgs_out_f110w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell370_bcgs_out_f125w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell370_bcgs_out_f140w_psf_bkg_drz.fits.gz',
                 bcgs_out + 'abell370_bcgs_out_f160w_bkg_drz_masked.fits.gz']
        
        RMS = [0.006150097033513798,   0.00018895763834042654,
               0.00030218235143297724, 0.0019826714331752926,
               0.0009539727215735511,  0.002089619907261409,
               0.0018212787251443347,  0.05877445329889749,
               0.0070035497120136065,  0.0031736748573591906,
               0.006013218785180418,   0.0041546267]
    
    if cluster == 'a1063' :
        cluster_redshift = 0.348
        first_path, second_path = 'abell1063clu_catalogs', 'abell1063clu_v3.9'
        misc = 'misc/abell1063clu_misc/'
        
        filters = ['f225w', 'f275w', 'f336w', 'f390w', # UVIS
                   'f435w', 'f475w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp', # ACS
                   'f105w', 'f110w', 'f125w', 'f140w', 'f160w'] # IR
        
        segPath = misc + 'photometry/abell1063_bcgs_out_detection_seg.fits.gz'
        
        psf_matched = misc + 'images/psf_matched/'
        bcgs_out = misc + 'images/bcgs_out/'
        paths = [psf_matched + 'abell1063_bcgs_out_f225w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell1063_bcgs_out_f275w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell1063_bcgs_out_f336w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell1063_bcgs_out_f390w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell1063_bcgs_out_f435w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell1063_bcgs_out_f475w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell1063_bcgs_out_f606w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell1063_bcgs_out_f625w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell1063_bcgs_out_f775w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell1063_bcgs_out_f814w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell1063_bcgs_out_f850lp_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell1063_bcgs_out_f105w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell1063_bcgs_out_f110w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell1063_bcgs_out_f125w_psf_bkg_drz.fits.gz',
                 psf_matched + 'abell1063_bcgs_out_f140w_psf_bkg_drz.fits.gz',
                 bcgs_out + 'abell1063_bcgs_out_f160w_bkg_drz_masked.fits.gz']
        
        RMS = [0.001009132276649849,   0.00019468210091193593,
               0.00021583139841912794, 0.002584890399809958,
               0.0009133420729312283,  0.007389349023214404,
               0.0016408712958101048,  0.00592355544506508,
               0.006542297661893592,   0.0016884383764566345,
               0.002285859571050604,   0.0033306179274798955,
               0.004810575815582492,   0.005576161652329551,
               0.0064130011504368615,  0.010719087]
    
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
        
        RMS = [0.00017821459857111745, 0.00014789605532426083,
               0.0007145838788770763,  0.001578221285547029,
               0.0018031112928780836,  0.0293976362500103,
               0.002737900021748479,   0.016844529370723973,
               0.0037151459]
    
    if cluster == 'm416' :
        cluster_redshift = 0.396
        first_path, second_path = 'macs0416clu_catalogs', 'macs0416clu_v3.9'
        misc = 'misc/macs0416clu_misc/'
        
        filters = ['f225w', 'f275w', 'f336w', 'f390w', # UVIS
                   'f435w', 'f475w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp', # ACS
                   'f105w', 'f110w', 'f125w', 'f140w', 'f160w'] # IR
        
        segPath = misc + 'photometry/macs0416_bcgs_out_detection_seg.fits.gz'
        
        psf_matched = misc + 'images/psf_matched/'
        bcgs_out = misc + 'images/bcgs_out/'
        paths = [psf_matched + 'macs0416_bcgs_out_f225w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0416_bcgs_out_f275w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0416_bcgs_out_f336w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0416_bcgs_out_f390w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0416_bcgs_out_f435w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0416_bcgs_out_f475w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0416_bcgs_out_f606w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0416_bcgs_out_f625w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0416_bcgs_out_f775w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0416_bcgs_out_f814w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0416_bcgs_out_f850lp_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0416_bcgs_out_f105w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0416_bcgs_out_f110w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0416_bcgs_out_f125w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0416_bcgs_out_f140w_psf_bkg_drz.fits.gz',
                 bcgs_out + 'macs0416_bcgs_out_f160w_bkg_drz_masked.fits.gz']
        
        RMS = [0.0005152721958643792, 0.0015232390831579111,
               0.008695907525215575,  0.0007815615517791924,
               0.0007717699865627695, 0.0020688076710073783,
               0.0025887395507373677, 0.002856857403925015,
               0.003277634974401299,  0.0016246667309102309,
               0.003120975668358854,  0.0034987608992435174,
               0.005269143193058465,  0.003089861725271567,
               0.00485683000893487,   0.0033343788]
    
    if cluster == 'm717' :
        cluster_redshift = 0.545
        first_path, second_path = 'macs0717clu_catalogs', 'macs0717clu_v3.9'
        misc = 'misc/macs0717clu_misc/'
        
        filters = ['f225w', 'f275w', 'f336w', 'f390w', # UVIS
                   'f435w', 'f475w', 'f555w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp', # ACS
                   'f105w', 'f110w', 'f125w', 'f140w', 'f160w'] # IR
        
        segPath = misc + 'photometry/macs0717_bcgs_out_detection_seg.fits.gz'
        
        psf_matched = misc + 'images/psf_matched/'
        bcgs_out = misc + 'images/bcgs_out/'
        paths = [psf_matched + 'macs0717_bcgs_out_f225w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0717_bcgs_out_f275w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0717_bcgs_out_f336w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0717_bcgs_out_f390w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0717_bcgs_out_f435w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0717_bcgs_out_f475w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0717_bcgs_out_f555w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0717_bcgs_out_f606w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0717_bcgs_out_f625w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0717_bcgs_out_f775w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0717_bcgs_out_f814w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0717_bcgs_out_f850lp_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0717_bcgs_out_f105w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0717_bcgs_out_f110w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0717_bcgs_out_f125w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs0717_bcgs_out_f140w_psf_bkg_drz.fits.gz',
                 bcgs_out + 'macs0717_bcgs_out_f160w_bkg_drz_masked.fits.gz']
        
        RMS = [0.00046611518484917884, 0.00017143037839937033,
               0.00018427351205806768, 0.0006596745694923388,
               0.005280505220527691,   0.007807015679126884,
               0.0009968294917573966,  0.017375741800531513,
               0.011421332862687431,   0.01133954743386766,
               0.001586713150211549,   0.009948818178895768,
               0.002251223101431824,   0.004196418715575267,
               0.0023549107302774687,  0.0034093925114306325,
               0.0027945212]
    
    if cluster == 'm1149' :
        cluster_redshift = 0.543
        first_path, second_path = 'macs1149clu_catalogs', 'macs1149clu_v3.9'
        misc = 'misc/macs1149clu_misc/'
        
        filters = ['f225w', 'f275w', 'f336w', 'f390w', # UVIS
                   'f435w', 'f475w', 'f555w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp', # ACS
                   'f105w', 'f110w', 'f125w', 'f140w', 'f160w'] # IR
        
        segPath = misc + 'photometry/macs1149_bcgs_out_detection_seg.fits.gz'
        
        psf_matched = misc + 'images/psf_matched/'
        bcgs_out = misc + 'images/bcgs_out/'
        paths = [psf_matched + 'macs1149_bcgs_out_f225w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs1149_bcgs_out_f275w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs1149_bcgs_out_f336w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs1149_bcgs_out_f390w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs1149_bcgs_out_f435w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs1149_bcgs_out_f475w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs1149_bcgs_out_f555w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs1149_bcgs_out_f606w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs1149_bcgs_out_f625w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs1149_bcgs_out_f775w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs1149_bcgs_out_f814w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs1149_bcgs_out_f850lp_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs1149_bcgs_out_f105w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs1149_bcgs_out_f110w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs1149_bcgs_out_f125w_psf_bkg_drz.fits.gz',
                 psf_matched + 'macs1149_bcgs_out_f140w_psf_bkg_drz.fits.gz',
                 bcgs_out + 'macs1149_bcgs_out_f160w_bkg_drz_masked.fits.gz']
        
        RMS = [0.0006290553961319055,  0.00018143297988564388,
               0.00018041860975653422, 0.0006786344022559748,
               0.00030392911762606366, 0.010212617906354983,
               0.004172468449132965,   0.0010194241574960708,
               0.011533475282206431,   0.012840495729004253,
               0.0013911904180016364,  0.010477723356620857,
               0.0025752394829416637,  0.006466929453150982,
               0.0038990153402697777,  0.003257393011206763,
               0.0040224046]
    
    # START
    
    os.makedirs('{}'.format(cluster), exist_ok=True) # ensure the output
        # directory is available
    if verbose :
        print('Created cluster directory.')
    
    # determine the final science objects, based on flags in the catalogs
    U_filtnum, V_filtnum, J_filtnum = 153, 155, 161
    core.determine_finalObjs_w_UVJ(cluster, 'id', cluster_redshift,
                                   first_path, second_path,
                                   U_filtnum, V_filtnum, J_filtnum,
                                   plot_all=True, plot_uvj=True,
                                   write_final_objs=True, write_regions=True)
    if verbose :
        print('Determined final objects and saved to fits file.')
    
    # save a file containing information about the filterset for the cluster
    field.save_filterset(cluster, filters)
    if verbose :
        print('Saved filterset to file.')
    
    # determine the rms values for each filter
    if calculate_rms :
        rms = field.determine_rms(segPath, paths)
        print(rms)
    else :
        rms = RMS
    if verbose :
        print('RMS determination complete.')
    
    # save all the cutouts for all the science objects
    final_objs_path = '{}/{}_final_objects.fits'.format(cluster, cluster)
    selection = 'Q' # only consider quiescent galaxies
    field.save_cutouts(cluster, final_objs_path, filters, rms, paths, segPath,
                       selection)
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
    
    # then on linux, move to the `hff/` directory, and run
    # python prospector.py
    
    return
