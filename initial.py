
def a370_params() :
    
    cluster_redshift = 0.375
    delta_z_lo, delta_z_hi = 0.01, 0.01
    force_spec = True
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
    
    RMS = [0.006150097033513798,   0.00018895763834042654,
           0.00030218235143297724, 0.0019826714331752926,
           0.0009539727215735511,  0.002089619907261409,
           0.0018212787251443347,  0.008016906056900567,
           0.0070035497120136065,  0.0031736748573591906,
           0.006013218785180418,   0.0041546267]
    
    return (cluster_redshift, delta_z_lo, delta_z_hi, force_spec,
            first_path, second_path, filters, segPath, files, models, RMS)

def a1063_params() :
    
    cluster_redshift = 0.348
    delta_z_lo, delta_z_hi = 0.01, 0.01
    force_spec = True
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
    
    RMS = [0.001009132276649849,   0.00019468210091193593,
           0.00021583139841912794, 0.002584890399809958,
           0.0009133420729312283,  0.00489489723092036,
           0.0016408712958101048,  0.00523063870614027,
           0.0061910310910946835,  0.0016884383764566345,
           0.002285859571050604,   0.0033306179274798955,
           0.004810575815582492,   0.005576161652329551,
           0.0064130011504368615,  0.010719087]
    
    return (cluster_redshift, delta_z_lo, delta_z_hi, force_spec,
            first_path, second_path, filters, segPath, files, models, RMS)

def a2744_params() :
    
    cluster_redshift = 0.308
    delta_z_lo, delta_z_hi = 0.01, 0.01
    force_spec = True
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
    
    RMS = [0.00017821459857111745, 0.00014789605532426083,
           0.0007145838788770763,  0.001578221285547029,
           0.0018031112928780836,  0.013637986621027736,
           0.002737900021748479,   0.009309231504260962,
           0.0037151459]
    
    return (cluster_redshift, delta_z_lo, delta_z_hi, force_spec,
            first_path, second_path, filters, segPath, files, models, RMS)

def m416_params() :
    
    cluster_redshift = 0.396
    delta_z_lo, delta_z_hi = 0.01, 0.01
    force_spec = True
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
    
    RMS = [0.0005152721958643792, 0.0015232390831579111,
           0.008695907525215575,  0.0007815615517791924,
           0.0007717699865627695, 0.0020688076710073783,
           0.0025887395507373677, 0.002856857403925015,
           0.003277634974401299,  0.0016246667309102309,
           0.003120975668358854,  0.0034987608992435174,
           0.005269143193058465,  0.003089861725271567,
           0.00485683000893487,   0.0033343788]
    
    return (cluster_redshift, delta_z_lo, delta_z_hi, force_spec,
            first_path, second_path, filters, segPath, files, models, RMS)

def m717_params() :
    
    cluster_redshift = 0.545
    delta_z_lo, delta_z_hi = 0.01, 0.01
    force_spec = True
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
    
    RMS = [0.00046611518484917884, 0.00017143037839937033,
           0.00018427351205806768, 0.0006596745694923388,
           0.005280505220527691,   0.007807015679126884,
           0.0009968294917573966,  0.017375741800531513,
           0.011421332862687431,   0.01133954743386766,
           0.001586713150211549,   0.009948818178895768,
           0.002251223101431824,   0.004196418715575267,
           0.0023549107302774687,  0.0034093925114306325,
           0.0027945212]
    
    return (cluster_redshift, delta_z_lo, delta_z_hi, force_spec,
            first_path, second_path, filters, segPath, files, models, RMS)

def m1149_params() :
    
    cluster_redshift = 0.543
    delta_z_lo, delta_z_hi = 0.01, 0.01
    force_spec = True
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
    
    RMS = [0.0006290553961319055,  0.00018143297988564388,
           0.00018041860975653422, 0.0006786344022559748,
           0.00030392911762606366, 0.010212617906354983,
           0.004172468449132965,   0.0010194241574960708,
           0.011533475282206431,   0.012840495729004253,
           0.0013911904180016364,  0.010477723356620857,
           0.0025752394829416637,  0.006466929453150982,
           0.0038990153402697777,  0.003257393011206763,
           0.0040224046]
    
    return (cluster_redshift, delta_z_lo, delta_z_hi, force_spec,
            first_path, second_path, filters, segPath, files, models, RMS)

def a370par_params() :
    
    cluster_redshift = 0.375
    delta_z_lo = cluster_redshift - (0.308-0.05)
    delta_z_hi = (0.545+0.05) - cluster_redshift
    force_spec = False
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
    
    RMS = [0.0002879315725811113, 0.0006588479431452871,
           0.0004609912601621799, 0.0004973539250495655,
           0.0007551672625252938, 0.0008845389370181092,
           0.0016675465]
    
    return (cluster_redshift, delta_z_lo, delta_z_hi, force_spec,
            first_path, second_path, filters, segPath, files, models, RMS)

def a1063par_params() :
    
    cluster_redshift = 0.348
    delta_z_lo = cluster_redshift - (0.308-0.05)
    delta_z_hi = (0.545+0.05) - cluster_redshift
    force_spec = False
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
    
    RMS = [0.00032123329138275986, 0.0009547855462766231,
           0.0006668530222700777,  0.0007476069912216449,
           0.001624492568493733,   0.0013993146495282485,
           0.002274977]
    
    return (cluster_redshift, delta_z_lo, delta_z_hi, force_spec,
            first_path, second_path, filters, segPath, files, models, RMS)

def a2744par_params() :
    
    cluster_redshift = 0.308
    delta_z_lo = cluster_redshift - (0.308-0.05)
    delta_z_hi = (0.545+0.05) - cluster_redshift
    force_spec = False
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
    
    RMS = [0.00030088137886912473, 0.007504526729094637,
           0.0005139121423771429,  0.0010082562671031548,
           0.0012806185239366261,  0.0017031766639187147,
           0.0014054665]
    
    return (cluster_redshift, delta_z_lo, delta_z_hi, force_spec,
            first_path, second_path, filters, segPath, files, models, RMS)

def m416par_params() :
    
    cluster_redshift = 0.396
    delta_z_lo = cluster_redshift - (0.308-0.05)
    delta_z_hi = (0.545+0.05) - cluster_redshift
    force_spec = False
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
    
    RMS = [0.00027932040174635613, 0.0005967158652637423,
           0.02390073357236199,    0.0005743923830976937,
           0.022203974558080843,   0.0007395513438158308,
           0.0009526387529398551,  0.003227937807108707,
           0.004652325]
    
    return (cluster_redshift, delta_z_lo, delta_z_hi, force_spec,
            first_path, second_path, filters, segPath, files, models, RMS)

def m717par_params() :
    
    cluster_redshift = 0.545
    delta_z_lo = cluster_redshift - (0.308-0.05)
    delta_z_hi = (0.545+0.05) - cluster_redshift
    force_spec = False
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
    
    RMS = [0.015750940989686293,  0.029209976296150752,
           0.004962645535653485,  0.0007025983894983889,
           0.0022194734390028562, 0.0011206273251126492,
           0.0023070623]
    
    return (cluster_redshift, delta_z_lo, delta_z_hi, force_spec,
            first_path, second_path, filters, segPath, files, models, RMS)

def m1149par_params() :
    
    cluster_redshift = 0.543
    delta_z_lo = cluster_redshift - (0.308-0.05)
    delta_z_hi = (0.545+0.05) - cluster_redshift
    force_spec = False
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
    
    RMS = [0.00029077593709390443, 0.006550107220462234,
           0.0014776043499247677,  0.0006278033486495236,
           0.000873025341883191,   0.0011653710498893614,
           0.0014654497]
    
    return (cluster_redshift, delta_z_lo, delta_z_hi, force_spec,
            first_path, second_path, filters, segPath, files, models, RMS)
