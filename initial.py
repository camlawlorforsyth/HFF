
import astropy.constants as const
import astropy.units as u

speed_of_light = const.c.to(u.km/u.s)

def a370_params() :
    
    cluster_redshift, sigma = 0.375, 1170*u.km/u.s # sigma from Dressler+ 1999
    delta_z = (3*sigma/speed_of_light)*(1 + cluster_redshift)
    redshift_type = 'z_spec'
    first_path, second_path = 'abell370clu_catalogs', 'abell370clu_v3.9'
    misc = 'misc/abell370clu_misc/'
    
    filters = ['f275w', 'f336w', # UVIS
               'f435w', 'f475w', 'f606w', 'f625w', 'f814w', # ACS
               'f105w', 'f110w', 'f125w', 'f140w', 'f160w'] # IR
    
    segPath = misc + 'photometry/abell370_bcgs_out_detection_seg.fits.gz'
    
    psf_matched = misc + 'images/psf_matched/'
    bcgs_out = misc + 'images/bcgs_out/'
    files = [psf_matched + 'abell370_bcgs_out_f275w_psf_bkg_drz.fits.gz',
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
    
    models_dir = misc + 'images/bcgs_models/'
    models = [models_dir + 'abell370_f275w_bcgs_model.fits.gz',
              models_dir + 'abell370_f336w_bcgs_model.fits.gz',
              models_dir + 'abell370_f435w_bcgs_model.fits.gz',
              models_dir + 'abell370_f475w_bcgs_model.fits.gz',
              models_dir + 'abell370_f606w_bcgs_model.fits.gz',
              models_dir + 'abell370_f625w_bcgs_model.fits.gz',
              models_dir + 'abell370_f814w_bcgs_model.fits.gz',
              models_dir + 'abell370_f105w_bcgs_model.fits.gz',
              models_dir + 'abell370_f110w_bcgs_model.fits.gz',
              models_dir + 'abell370_f125w_bcgs_model.fits.gz',
              models_dir + 'abell370_f140w_bcgs_model.fits.gz',
              models_dir + 'abell370_f160w_bcgs_model.fits.gz',]
    
    RMS = [0.013416409320783362,  0.0004135004717951535,
           0.000521624855906762,  0.0028024270920049606,
           0.001652756449667403,  0.003576384205110058,
           0.0028881052850907815, 0.0212314134774129,
           0.016790748246994322,  0.008905223215942533,
           0.013136195492602264,  0.008631396]
    
    return (cluster_redshift, delta_z, delta_z, redshift_type,
            first_path, second_path, filters, segPath, files, models, RMS)

def a1063_params() :
    
    cluster_redshift, sigma = 0.348, 1840*u.km/u.s # sigma from Lotz+ 2017
    delta_z = (3*sigma/speed_of_light)*(1 + cluster_redshift)
    redshift_type = 'z_spec'
    first_path, second_path = 'abell1063clu_catalogs', 'abell1063clu_v3.9'
    misc = 'misc/abell1063clu_misc/'
    
    filters = ['f225w', 'f275w', 'f336w', 'f390w', # UVIS
               'f435w', 'f475w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp', # ACS
               'f105w', 'f110w', 'f125w', 'f140w', 'f160w'] # IR
    
    segPath = misc + 'photometry/abell1063_bcgs_out_detection_seg.fits.gz'
    
    psf_matched = misc + 'images/psf_matched/'
    bcgs_out = misc + 'images/bcgs_out/'
    files = [psf_matched + 'abell1063_bcgs_out_f225w_psf_bkg_drz.fits.gz',
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
    
    models_dir = misc + 'images/bcgs_models/'
    models = ['models/' + 'abell1063_f225w_bcgs_model.fits.gz', # user gen.
              models_dir + 'abell1063_f275w_bcgs_model.fits.gz',
              models_dir + 'abell1063_f336w_bcgs_model.fits.gz',
              models_dir + 'abell1063_f390w_bcgs_model.fits.gz',
              models_dir + 'abell1063_f435w_bcgs_model.fits.gz',
              models_dir + 'abell1063_f475w_bcgs_model.fits.gz',
              models_dir + 'abell1063_f606w_bcgs_model.fits.gz',
              models_dir + 'abell1063_f625w_bcgs_model.fits.gz',
              models_dir + 'abell1063_f775w_bcgs_model.fits.gz',
              models_dir + 'abell1063_f814w_bcgs_model.fits.gz',
              models_dir + 'abell1063_f850lp_bcgs_model.fits.gz',
              models_dir + 'abell1063_f105w_bcgs_model.fits.gz',
              models_dir + 'abell1063_f110w_bcgs_model.fits.gz',
              models_dir + 'abell1063_f125w_bcgs_model.fits.gz',
              models_dir + 'abell1063_f140w_bcgs_model.fits.gz',
              models_dir + 'abell1063_f160w_bcgs_model.fits.gz']
    
    RMS = [0.001983245662268219,  0.0003760601291170469,
           0.0004168895372459394, 0.004688679819459995,
           0.0012721294364683093, 0.006964068399754354,
           0.002283272944509662,  0.007449095046361041,
           0.00879840020557511,   0.0022877748692428546,
           0.003244247526901425,  0.007230401519855124,
           0.011814075886542935,  0.013167607483792613,
           0.013916335017667853,  0.023791686]
    
    return (cluster_redshift, delta_z, delta_z, redshift_type,
            first_path, second_path, filters, segPath, files, models, RMS)

def a2744_params() :
    
    cluster_redshift, sigma = 0.308, 1497*u.km/u.s # sigma from Owers+ 2011
    delta_z = (3*sigma/speed_of_light)*(1 + cluster_redshift)
    redshift_type = 'z_spec'
    first_path, second_path = 'abell2744clu_catalogs', 'abell2744clu_v3.9'
    misc = 'misc/abell2744clu_misc/'
    
    filters = ['f275w', 'f336w', # UVIS
               'f435w', 'f606w', 'f814w', # ACS
               'f105w', 'f125w', 'f140w', 'f160w'] # IR
    
    segPath = misc + 'photometry/abell2744_bcgs_out_detection_seg.fits.gz'
    
    psf_matched = misc + 'images/psf_matched/'
    bcgs_out = misc + 'images/bcgs_out/'
    files = [psf_matched + 'abell2744_bcgs_out_f275w_psf_bkg_drz.fits.gz',
             psf_matched + 'abell2744_bcgs_out_f336w_psf_bkg_drz.fits.gz',
             psf_matched + 'abell2744_bcgs_out_f435w_psf_bkg_drz.fits.gz',
             psf_matched + 'abell2744_bcgs_out_f606w_psf_bkg_drz.fits.gz',
             psf_matched + 'abell2744_bcgs_out_f814w_psf_bkg_drz.fits.gz',
             psf_matched + 'abell2744_bcgs_out_f105w_psf_bkg_drz.fits.gz',
             psf_matched + 'abell2744_bcgs_out_f125w_psf_bkg_drz.fits.gz',
             psf_matched + 'abell2744_bcgs_out_f140w_psf_bkg_drz.fits.gz',
             bcgs_out + 'abell2744_bcgs_out_f160w_bkg_drz_masked.fits.gz']
    
    models_dir = misc + 'images/bcgs_models/'
    models = [models_dir + 'abell2744_f275w_bcgs_model.fits.gz',
              models_dir + 'abell2744_f336w_bcgs_model.fits.gz',
              models_dir + 'abell2744_f435w_bcgs_model.fits.gz',
              models_dir + 'abell2744_f606w_bcgs_model.fits.gz',
              models_dir + 'abell2744_f814w_bcgs_model.fits.gz',
              models_dir + 'abell2744_f105w_bcgs_model.fits.gz',
              models_dir + 'abell2744_f125w_bcgs_model.fits.gz',
              models_dir + 'abell2744_f140w_bcgs_model.fits.gz',
              models_dir + 'abell2744_f160w_bcgs_model.fits.gz']
    
    RMS = [0.00038855623931806355, 0.0003223930286361111,
           0.0009910847767761874,  0.002203596225349592,
           0.002516569462848094,   0.03302668875071387,
           0.0078074336316639814,  0.025402163365869854,
           0.010081889]
    
    return (cluster_redshift, delta_z, delta_z, redshift_type,
            first_path, second_path, filters, segPath, files, models, RMS)

def m416_params() :
    
    cluster_redshift, sigma = 0.396, 955*u.km/u.s #sigma from Jauzac+ 2014
    delta_z = (3*sigma/speed_of_light)*(1 + cluster_redshift)
    redshift_type = 'z_spec'
    first_path, second_path = 'macs0416clu_catalogs', 'macs0416clu_v3.9'
    misc = 'misc/macs0416clu_misc/'
    
    filters = ['f225w', 'f275w', 'f336w', 'f390w', # UVIS
               'f435w', 'f475w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp', # ACS
               'f105w', 'f110w', 'f125w', 'f140w', 'f160w'] # IR
    
    segPath = misc + 'photometry/macs0416_bcgs_out_detection_seg.fits.gz'
    
    psf_matched = misc + 'images/psf_matched/'
    bcgs_out = misc + 'images/bcgs_out/'
    files = [psf_matched + 'macs0416_bcgs_out_f225w_psf_bkg_drz.fits.gz',
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
    
    models_dir = misc + 'images/bcgs_models/'
    models = [models_dir + 'macs0416_f225w_bcgs_model.fits.gz',
              models_dir + 'macs0416_f275w_bcgs_model.fits.gz',
              models_dir + 'macs0416_f336w_bcgs_model.fits.gz',
              models_dir + 'macs0416_f390w_bcgs_model.fits.gz',
              models_dir + 'macs0416_f435w_bcgs_model.fits.gz',
              models_dir + 'macs0416_f475w_bcgs_model.fits.gz',
              models_dir + 'macs0416_f606w_bcgs_model.fits.gz',
              models_dir + 'macs0416_f625w_bcgs_model.fits.gz',
              models_dir + 'macs0416_f775w_bcgs_model.fits.gz',
              models_dir + 'macs0416_f814w_bcgs_model.fits.gz',
              models_dir + 'macs0416_f850lp_bcgs_model.fits.gz',
              models_dir + 'macs0416_f105w_bcgs_model.fits.gz',
              models_dir + 'macs0416_f110w_bcgs_model.fits.gz',
              models_dir + 'macs0416_f125w_bcgs_model.fits.gz',
              models_dir + 'macs0416_f140w_bcgs_model.fits.gz',
              models_dir + 'macs0416_f160w_bcgs_model.fits.gz']
    
    RMS = [0.001022513145144539,  0.002953053080221736,
           0.01686053443037426,   0.0017191001741423658,
           0.001095424754206102,  0.0029546951779546755,
           0.0038279024768208504, 0.004080879449397976,
           0.0047054578998890475, 0.002283975917149887,
           0.004456799032220651,  0.007748914359525832,
           0.013362314701066995,  0.006921186187522245,
           0.010957062338618869,  0.0073933504]
    
    return (cluster_redshift, delta_z, delta_z, redshift_type,
            first_path, second_path, filters, segPath, files, models, RMS)

def m717_params() :
    
    cluster_redshift, sigma = 0.545, 1660*u.km/u.s # sigma from Ebeling+ 2007
    delta_z = (3*sigma/speed_of_light)*(1 + cluster_redshift)
    redshift_type = 'z_spec'
    first_path, second_path = 'macs0717clu_catalogs', 'macs0717clu_v3.9'
    misc = 'misc/macs0717clu_misc/'
    
    filters = ['f225w', 'f275w', 'f336w', 'f390w', # UVIS
               'f435w', 'f475w', 'f555w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp', # ACS
               'f105w', 'f110w', 'f125w', 'f140w', 'f160w'] # IR
    
    segPath = misc + 'photometry/macs0717_bcgs_out_detection_seg.fits.gz'
    
    psf_matched = misc + 'images/psf_matched/'
    bcgs_out = misc + 'images/bcgs_out/'
    files = [psf_matched + 'macs0717_bcgs_out_f225w_psf_bkg_drz.fits.gz',
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
    
    models_dir = misc + 'images/bcgs_models/'
    models = [models_dir + 'macs0717_f225w_bcgs_model.fits.gz',
              models_dir + 'macs0717_f275w_bcgs_model.fits.gz',
              models_dir + 'macs0717_f336w_bcgs_model.fits.gz',
              models_dir + 'macs0717_f390w_bcgs_model.fits.gz',
              models_dir + 'macs0717_f435w_bcgs_model.fits.gz',
              models_dir + 'macs0717_f475w_bcgs_model.fits.gz',
              models_dir + 'macs0717_f555w_bcgs_model.fits.gz',
              models_dir + 'macs0717_f606w_bcgs_model.fits.gz',
              models_dir + 'macs0717_f625w_bcgs_model.fits.gz',
              models_dir + 'macs0717_f775w_bcgs_model.fits.gz',
              models_dir + 'macs0717_f814w_bcgs_model.fits.gz',
              models_dir + 'macs0717_f850lp_bcgs_model.fits.gz',
              models_dir + 'macs0717_f105w_bcgs_model.fits.gz',
              models_dir + 'macs0717_f110w_bcgs_model.fits.gz',
              models_dir + 'macs0717_f125w_bcgs_model.fits.gz',
              models_dir + 'macs0717_f140w_bcgs_model.fits.gz',
              models_dir + 'macs0717_f160w_bcgs_model.fits.gz']
    
    RMS = [0.0009173335084495755,  0.0003246848223213115,
           0.00034964614334836953, 0.001443331520114308,
           0.007392274853508328,   0.011249858367578901,
           0.0015317426955391118,  0.023738641244255478,
           0.016355449526803338,   0.01626875186766334,
           0.00209748276318978,    0.014315095827882308,
           0.004786703989946258,   0.010780858656126994,
           0.005199413123335132,   0.0072590100581836495,
           0.005936407]
    
    return (cluster_redshift, delta_z, delta_z, redshift_type,
            first_path, second_path, filters, segPath, files, models, RMS)

def m1149_params() :
    
    cluster_redshift, sigma = 0.543, 1840*u.km/u.s # sigma from Ebeling+ 2007
    delta_z = (3*sigma/speed_of_light)*(1 + cluster_redshift)
    redshift_type = 'z_spec'
    first_path, second_path = 'macs1149clu_catalogs', 'macs1149clu_v3.9'
    misc = 'misc/macs1149clu_misc/'
    
    filters = ['f225w', 'f275w', 'f336w', 'f390w', # UVIS
               'f435w', 'f475w', 'f555w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp', # ACS
               'f105w', 'f110w', 'f125w', 'f140w', 'f160w'] # IR
    
    segPath = misc + 'photometry/macs1149_bcgs_out_detection_seg.fits.gz'
    
    psf_matched = misc + 'images/psf_matched/'
    bcgs_out = misc + 'images/bcgs_out/'
    files = [psf_matched + 'macs1149_bcgs_out_f225w_psf_bkg_drz.fits.gz',
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
    
    models_dir = misc + 'images/bcgs_models/'
    models = ['models/' + 'macs1149_f225w_bcgs_model.fits.gz', # user gen.
              'models/' + 'macs1149_f275w_bcgs_model.fits.gz', # user gen.
              models_dir + 'macs1149_f336w_bcgs_model.fits.gz',
              models_dir + 'macs1149_f390w_bcgs_model.fits.gz',
              models_dir + 'macs1149_f435w_bcgs_model.fits.gz',
              models_dir + 'macs1149_f475w_bcgs_model.fits.gz',
              models_dir + 'macs1149_f555w_bcgs_model.fits.gz',
              models_dir + 'macs1149_f606w_bcgs_model.fits.gz',
              models_dir + 'macs1149_f625w_bcgs_model.fits.gz',
              models_dir + 'macs1149_f775w_bcgs_model.fits.gz',
              models_dir + 'macs1149_f814w_bcgs_model.fits.gz',
              models_dir + 'macs1149_f850lp_bcgs_model.fits.gz',
              models_dir + 'macs1149_f105w_bcgs_model.fits.gz',
              models_dir + 'macs1149_f110w_bcgs_model.fits.gz',
              models_dir + 'macs1149_f125w_bcgs_model.fits.gz',
              models_dir + 'macs1149_f140w_bcgs_model.fits.gz',
              models_dir + 'macs1149_f160w_bcgs_model.fits.gz']
    
    RMS = [0.0013397469733752053, 0.00037920130212222537,
           0.0003773592403834817, 0.0016005177608797417,
           0.0005056412424501541, 0.01563980579934457,
           0.006877902429964284,  0.0016964515651718942,
           0.01774489939044516,   0.019729500832344284,
           0.00226753277782757,   0.01609277183208197,
           0.005710196804617347,  0.01754505598287007,
           0.007802336747616293,  0.008701588266849821,
           0.008023425]
    
    return (cluster_redshift, delta_z, delta_z, redshift_type,
            first_path, second_path, filters, segPath, files, models, RMS)

def a370par_params() :
    
    cluster_redshift = 0.375
    delta_z_lo = cluster_redshift - 0.2884056845219235 # min. z from clusters
    delta_z_hi = 0.5714108548187694 - cluster_redshift # max. z from clusters
    redshift_type = 'z'
    first_path, second_path = 'abell370par_catalogs', 'abell370par_v3.9'
    misc = 'misc/abell370par_misc/'
    
    filters = ['f435w', 'f606w', 'f814w', # ACS
               'f105w', 'f125w', 'f140w', 'f160w'] # IR
    
    segPath = misc + 'photometry/abell370par_bcgs_out_detection_seg.fits.gz'
    
    psf_matched = misc + 'images/psf_matched/'
    bcgs_out = misc + 'images/bcgs_out/'
    files = [psf_matched + 'abell370par_bcgs_out_f435w_psf_bkg_drz.fits.gz',
             psf_matched + 'abell370par_bcgs_out_f606w_psf_bkg_drz.fits.gz',
             psf_matched + 'abell370par_bcgs_out_f814w_psf_bkg_drz.fits.gz',
             psf_matched + 'abell370par_bcgs_out_f105w_psf_bkg_drz.fits.gz',
             psf_matched + 'abell370par_bcgs_out_f125w_psf_bkg_drz.fits.gz',
             psf_matched + 'abell370par_bcgs_out_f140w_psf_bkg_drz.fits.gz',
             bcgs_out + 'abell370par_bcgs_out_f160w_bkg_drz_masked.fits.gz']
    
    models_dir = misc + 'images/bcgs_models/'
    models = [models_dir + 'abell370par_f435w_bcgs_model.fits.gz',
              models_dir + 'abell370par_f606w_bcgs_model.fits.gz',
              models_dir + 'abell370par_f814w_bcgs_model.fits.gz',
              models_dir + 'abell370par_f105w_bcgs_model.fits.gz',
              models_dir + 'abell370par_f125w_bcgs_model.fits.gz',
              models_dir + 'abell370par_f140w_bcgs_model.fits.gz',
              models_dir + 'abell370par_f160w_bcgs_model.fits.gz']
    
    RMS = [0.0004652031764877602, 0.001064855033812716,
           0.0007438264566643753, 0.0012629274388352607,
           0.0019222345444720006, 0.0022521165138177563,
           0.0042280676]
    
    return (cluster_redshift, delta_z_lo, delta_z_hi, redshift_type,
            first_path, second_path, filters, segPath, files, models, RMS)

def a1063par_params() :
    
    cluster_redshift = 0.348
    delta_z_lo = cluster_redshift - 0.2884056845219235
    delta_z_hi = 0.5714108548187694 - cluster_redshift
    redshift_type = 'z'
    first_path, second_path = 'abell1063par_catalogs', 'abell1063par_v3.9'
    misc = 'misc/abell1063par_misc/'
    
    filters = ['f435w', 'f606w', 'f814w', # ACS
               'f105w', 'f125w', 'f140w', 'f160w'] # IR
    
    segPath = misc + 'photometry/abell1063par_bcgs_out_detection_seg.fits.gz'
    
    psf_matched = misc + 'images/psf_matched/'
    bcgs_out = misc + 'images/bcgs_out/'
    files = [psf_matched + 'abell1063par_bcgs_out_f435w_psf_bkg_drz.fits.gz',
             psf_matched + 'abell1063par_bcgs_out_f606w_psf_bkg_drz.fits.gz',
             psf_matched + 'abell1063par_bcgs_out_f814w_psf_bkg_drz.fits.gz',
             psf_matched + 'abell1063par_bcgs_out_f105w_psf_bkg_drz.fits.gz',
             psf_matched + 'abell1063par_bcgs_out_f125w_psf_bkg_drz.fits.gz',
             psf_matched + 'abell1063par_bcgs_out_f140w_psf_bkg_drz.fits.gz',
             bcgs_out + 'abell1063par_bcgs_out_f160w_bkg_drz_masked.fits.gz']
    
    models_dir = misc + 'images/bcgs_models/'
    models = [models_dir + 'abell1063par_f435w_bcgs_model.fits.gz',
              models_dir + 'abell1063par_f606w_bcgs_model.fits.gz',
              models_dir + 'abell1063par_f814w_bcgs_model.fits.gz',
              models_dir + 'abell1063par_f105w_bcgs_model.fits.gz',
              models_dir + 'abell1063par_f125w_bcgs_model.fits.gz',
              models_dir + 'abell1063par_f140w_bcgs_model.fits.gz',
              models_dir + 'abell1063par_f160w_bcgs_model.fits.gz']
    
    RMS = [0.0005126545825005789, 0.0015236742604596359,
           0.0010727261928403882, 0.0018839704626364162,
           0.003545452030027999,  0.0035276791011699285,
           0.0049573397]
    
    return (cluster_redshift, delta_z_lo, delta_z_hi, redshift_type,
            first_path, second_path, filters, segPath, files, models, RMS)

def a2744par_params() :
    
    cluster_redshift = 0.308
    delta_z_lo = cluster_redshift - 0.2884056845219235
    delta_z_hi = 0.5714108548187694 - cluster_redshift
    redshift_type = 'z'
    first_path, second_path = 'abell2744par_catalogs', 'abell2744par_v3.9'
    misc = 'misc/abell2744par_misc/'
    
    filters = ['f435w', 'f606w', 'f814w', # ACS
               'f105w', 'f125w', 'f140w', 'f160w'] # IR
    
    segPath = misc + 'photometry/abell2744par_bcgs_out_detection_seg.fits.gz'
    
    psf_matched = misc + 'images/psf_matched/'
    bcgs_out = misc + 'images/bcgs_out/'
    files = [psf_matched + 'abell2744par_bcgs_out_f435w_psf_bkg_drz.fits.gz',
             psf_matched + 'abell2744par_bcgs_out_f606w_psf_bkg_drz.fits.gz',
             psf_matched + 'abell2744par_bcgs_out_f814w_psf_bkg_drz.fits.gz',
             psf_matched + 'abell2744par_bcgs_out_f105w_psf_bkg_drz.fits.gz',
             psf_matched + 'abell2744par_bcgs_out_f125w_psf_bkg_drz.fits.gz',
             psf_matched + 'abell2744par_bcgs_out_f140w_psf_bkg_drz.fits.gz',
             bcgs_out + 'abell2744par_bcgs_out_f160w_bkg_drz_masked.fits.gz']
    
    models_dir = misc + 'images/bcgs_models/'
    models = [models_dir + 'abell2744par_f435w_bcgs_model.fits.gz',
              models_dir + 'abell2744par_f606w_bcgs_model.fits.gz',
              models_dir + 'abell2744par_f814w_bcgs_model.fits.gz',
              models_dir + 'abell2744par_f105w_bcgs_model.fits.gz',
              models_dir + 'abell2744par_f125w_bcgs_model.fits.gz',
              models_dir + 'abell2744par_f140w_bcgs_model.fits.gz',
              models_dir + 'abell2744par_f160w_bcgs_model.fits.gz']
    
    RMS = [0.000482993783866937,  0.011711077398109594,
           0.0008324629536602113, 0.0026043154700675906,
           0.0032909964401797065, 0.0043810722888363175,
           0.0036131504]
    
    return (cluster_redshift, delta_z_lo, delta_z_hi, redshift_type,
            first_path, second_path, filters, segPath, files, models, RMS)

def m416par_params() :
    
    cluster_redshift = 0.396
    delta_z_lo = cluster_redshift - 0.2884056845219235
    delta_z_hi = 0.5714108548187694 - cluster_redshift
    redshift_type = 'z'
    first_path, second_path = 'macs0416par_catalogs', 'macs0416par_v3.9'
    misc = 'misc/macs0416par_misc/'
    
    filters = ['f435w', 'f606w', 'f775w', 'f814w', 'f850lp', # ACS
               'f105w', 'f125w', 'f140w', 'f160w'] # IR
    
    segPath = misc + 'photometry/macs0416par_bcgs_out_detection_seg.fits.gz'
    
    psf_matched = misc + 'images/psf_matched/'
    bcgs_out = misc + 'images/bcgs_out/'
    files = [psf_matched + 'macs0416par_bcgs_out_f435w_psf_bkg_drz.fits.gz',
             psf_matched + 'macs0416par_bcgs_out_f606w_psf_bkg_drz.fits.gz',
             psf_matched + 'macs0416par_bcgs_out_f775w_psf_bkg_drz.fits.gz',
             psf_matched + 'macs0416par_bcgs_out_f814w_psf_bkg_drz.fits.gz',
             psf_matched + 'macs0416par_bcgs_out_f850lp_psf_bkg_drz.fits.gz',
             psf_matched + 'macs0416par_bcgs_out_f105w_psf_bkg_drz.fits.gz',
             psf_matched + 'macs0416par_bcgs_out_f125w_psf_bkg_drz.fits.gz',
             psf_matched + 'macs0416par_bcgs_out_f140w_psf_bkg_drz.fits.gz',
             bcgs_out + 'macs0416par_bcgs_out_f160w_bkg_drz_masked.fits.gz']
    
    models_dir = misc + 'images/bcgs_models/'
    models = [models_dir + 'macs0416par_f435w_bcgs_model.fits.gz',
              models_dir + 'macs0416par_f606w_bcgs_model.fits.gz',
              models_dir + 'macs0416par_f775w_bcgs_model.fits.gz',
              models_dir + 'macs0416par_f814w_bcgs_model.fits.gz',
              models_dir + 'macs0416par_f850lp_bcgs_model.fits.gz',
              models_dir + 'macs0416par_f105w_bcgs_model.fits.gz',
              models_dir + 'macs0416par_f125w_bcgs_model.fits.gz',
              models_dir + 'macs0416par_f140w_bcgs_model.fits.gz',
              models_dir + 'macs0416par_f160w_bcgs_model.fits.gz']
    
    RMS = [0.00045506415746966684, 0.0009728186213636664,
           0.03149169217228643,    0.0009347899081887372,
           0.029257469814343726,   0.0019231930413740048,
           0.002467709853135251,   0.008358303082935914,
           0.012034563]
    
    return (cluster_redshift, delta_z_lo, delta_z_hi, redshift_type,
            first_path, second_path, filters, segPath, files, models, RMS)

def m717par_params() :
    
    cluster_redshift = 0.545
    delta_z_lo = cluster_redshift - 0.2884056845219235
    delta_z_hi = 0.5714108548187694 - cluster_redshift
    redshift_type = 'z'
    first_path, second_path = 'macs0717par_catalogs', 'macs0717par_v3.9'
    misc = 'misc/macs0717par_misc/'
    
    filters = ['f435w', 'f606w', 'f814w', # ACS
               'f105w', 'f125w', 'f140w', 'f160w'] # IR
    
    segPath = misc + 'photometry/macs0717par_bcgs_out_detection_seg.fits.gz'
    
    psf_matched = misc + 'images/psf_matched/'
    bcgs_out = misc + 'images/bcgs_out/'
    files = [psf_matched + 'macs0717par_bcgs_out_f435w_psf_bkg_drz.fits.gz',
             psf_matched + 'macs0717par_bcgs_out_f606w_psf_bkg_drz.fits.gz',
             psf_matched + 'macs0717par_bcgs_out_f814w_psf_bkg_drz.fits.gz',
             psf_matched + 'macs0717par_bcgs_out_f105w_psf_bkg_drz.fits.gz',
             psf_matched + 'macs0717par_bcgs_out_f125w_psf_bkg_drz.fits.gz',
             psf_matched + 'macs0717par_bcgs_out_f140w_psf_bkg_drz.fits.gz',
             bcgs_out + 'macs0717par_bcgs_out_f160w_bkg_drz_masked.fits.gz']
    
    models_dir = misc + 'images/bcgs_models/'
    models = [models_dir + 'macs0717par_f435w_bcgs_model.fits.gz',
              models_dir + 'macs0717par_f606w_bcgs_model.fits.gz',
              models_dir + 'macs0717par_f814w_bcgs_model.fits.gz',
              models_dir + 'macs0717par_f105w_bcgs_model.fits.gz',
              models_dir + 'macs0717par_f125w_bcgs_model.fits.gz',
              models_dir + 'macs0717par_f140w_bcgs_model.fits.gz',
              models_dir + 'macs0717par_f160w_bcgs_model.fits.gz']
    
    RMS = [0.0244939716339791,    0.04543024417500927,
           0.007849304440235873,  0.0018605373200842343,
           0.0050811364834191285, 0.0029706232359012624,
           0.005247639]
    
    return (cluster_redshift, delta_z_lo, delta_z_hi, redshift_type,
            first_path, second_path, filters, segPath, files, models, RMS)

def m1149par_params() :
    
    cluster_redshift = 0.543
    delta_z_lo = cluster_redshift - 0.2884056845219235
    delta_z_hi = 0.5714108548187694 - cluster_redshift
    redshift_type = 'z'
    first_path, second_path = 'macs1149par_catalogs', 'macs1149par_v3.9'
    misc = 'misc/macs1149par_misc/'
    
    filters = ['f435w', 'f606w', 'f814w', # ACS
               'f105w', 'f125w', 'f140w', 'f160w'] # IR
    
    segPath = misc + 'photometry/macs1149par_bcgs_out_detection_seg.fits.gz'
    
    psf_matched = misc + 'images/psf_matched/'
    bcgs_out = misc + 'images/bcgs_out/'
    files = [psf_matched + 'macs1149par_bcgs_out_f435w_psf_bkg_drz.fits.gz',
             psf_matched + 'macs1149par_bcgs_out_f606w_psf_bkg_drz.fits.gz',
             psf_matched + 'macs1149par_bcgs_out_f814w_psf_bkg_drz.fits.gz',
             psf_matched + 'macs1149par_bcgs_out_f105w_psf_bkg_drz.fits.gz',
             psf_matched + 'macs1149par_bcgs_out_f125w_psf_bkg_drz.fits.gz',
             psf_matched + 'macs1149par_bcgs_out_f140w_psf_bkg_drz.fits.gz',
             bcgs_out + 'macs1149par_bcgs_out_f160w_bkg_drz_masked.fits.gz']
    
    models_dir = misc + 'images/bcgs_models/'
    models = [models_dir + 'macs1149par_f435w_bcgs_model.fits.gz',
              models_dir + 'macs1149par_f606w_bcgs_model.fits.gz',
              models_dir + 'macs1149par_f814w_bcgs_model.fits.gz',
              models_dir + 'macs1149par_f105w_bcgs_model.fits.gz',
              models_dir + 'macs1149par_f125w_bcgs_model.fits.gz',
              models_dir + 'macs1149par_f140w_bcgs_model.fits.gz',
              models_dir + 'macs1149par_f160w_bcgs_model.fits.gz']
    
    RMS = [0.0004692070016910032, 0.009543790664211372,
           0.0021566684274162314, 0.0015926369948471837,
           0.0022168557361695275, 0.0029588726664742273,
           0.0035911163]
    
    return (cluster_redshift, delta_z_lo, delta_z_hi, redshift_type,
            first_path, second_path, filters, segPath, files, models, RMS)
