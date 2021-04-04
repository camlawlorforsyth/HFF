
import binning
from cluster import Cluster
import core
import plotting as plt
import photometry
import voronoi

def premain(save_filters=False, plot_curves=False) :
    
    if save_filters :
        core.build_filter_table(save_par=True)
    
    if plot_curves :
        plt.hst_transmission_curves(core.build_filter_table())
    
    return

def main(binning_scheme='annuli') :
    
    # determine the final science objects for Abell 2744, based on flags
    # in the catalogs
    core.determine_finalObjs_w_UVJ('abell_2744', 'id', 0.308,
                                   'abell2744clu_catalogs',
                                   'abell2744clu_v3.9', 153, 155, 161)
    
    # create a Cluster instance for Abell 2744
    a2744 = initialize('a2744')
    
    # determine the rms values for each filter for Abell 2744
    rms = a2744.determine_rms()
    print(rms)
    
    # save all the cutouts for all the science objects for Abell 2744
    a2744.save_cutouts()
    
    # now vorbin the cutouts for each galaxy and save the resulting files
    if binning_scheme == 'annuli' :
        binning.bin_all('a2744/cutouts/', 'a2744/bins/')
    else :
        voronoi.vorbin_all('a2744/cutouts/', 'a2744/vorbins/')
    
    # then determine the flux for every vorbin for each galaxy, saving the
    # photometry for each galaxy into a separate file
    photometry.determine_fluxes('a2744/bins/', 'a2744/cutouts/',
                                'a2744/photometry/', a2744.filters)
    
    # then on linux, move to the `hff/` directory, and run
    # python prospector.py
    
    return

def initialize(cluster) :
    
    if cluster == 'a2744' :
        
        filters = ['f275w', 'f336w', # UVIS
                   'f435w', 'f606w', 'f814w', # ACS
                   'f105w', 'f125w', 'f140w', 'f160w'] # IR
        
        misc = 'misc/abell2744clu_misc/'
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
        
        segPath = misc + 'photometry/abell2744_bcgs_out_detection_seg.fits.gz'
        
        finalSciObjs = cluster + '/' + cluster + '_final_objects.fits'
        
        rms = [0.00017821459857111745, 0.00014789605532426083,
               0.0007145838788770763,  0.001578221285547029,
               0.0018031112928780836,  0.0293976362500103,
               0.002737900021748479,   0.016844529370723973,
               0.0037151459]
        
        a2744 = Cluster(cluster, filters, paths, segPath, finalSciObjs, rms)
        
        return a2744
    else :
        return None

# determine_finalObjs_w_UVJ('abell_2744', 'id', 0.308, 'abell2744clu_catalogs',
#                           'abell2744clu_v3.9', 153, 155, 161)
# determine_finalObjs_w_UVJ('abell_370', 'id', 0.375, 'abell370clu_catalogs',
#                           'abell370clu_v3.9', 153, 155, 161)
# determine_finalObjs_w_UVJ('abell_S1063', 'id', 0.348, 'abell1063clu_catalogs',
#                           'abell1063clu_v3.9', 153, 155, 161)
# determine_finalObjs_w_UVJ('macs_0416', 'id', 0.392, 'macs0416clu_catalogs',
#                           'macs0416clu_v3.9', 153, 155, 161)
# determine_finalObjs_w_UVJ('macs_0717', 'id', 0.5458, 'macs0717clu_catalogs',
#                           'macs0717clu_v3.9', 153, 155, 161)
# determine_finalObjs_w_UVJ('macs_1149', 'id', 0.543, 'macs1149clu_catalogs',
#                           'macs1149clu_v3.9', 153, 155, 161)

# still necessary?
# determine_final_objects('abell_370', 'id', 0.375)
# determine_final_objects('abell_S1063', 'id', 0.348)
# determine_final_objects('abell_2744', 'id', 0.308, plot=True)
# determine_final_objects('macs_0416', 'id', 0.392)
# determine_final_objects('macs_0717', 'id', 0.5458)
# determine_final_objects('macs_1149', 'id', 0.543)

# still necessary?
# combine_UVJ('abell2744clu_catalogs', 'abell2744clu_v3.9', 153, 155, 161)
