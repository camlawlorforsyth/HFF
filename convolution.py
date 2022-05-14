
from astropy.convolution import convolve
from astropy.io import fits

def all_convolutions(cluster) :
    
    if cluster == 'a370' :
        models = ['misc/abell370clu_misc/images/bcgs_models/abell370_f275w_bcgs_model.fits.gz',
                  'misc/abell370clu_misc/images/bcgs_models/abell370_f336w_bcgs_model.fits.gz',
                  'misc/abell370clu_misc/images/bcgs_models/abell370_f435w_bcgs_model.fits.gz',
                  'misc/abell370clu_misc/images/bcgs_models/abell370_f475w_bcgs_model.fits.gz',
                  'misc/abell370clu_misc/images/bcgs_models/abell370_f606w_bcgs_model.fits.gz',
                  'misc/abell370clu_misc/images/bcgs_models/abell370_f625w_bcgs_model.fits.gz',
                  'misc/abell370clu_misc/images/bcgs_models/abell370_f814w_bcgs_model.fits.gz',
                  'misc/abell370clu_misc/images/bcgs_models/abell370_f105w_bcgs_model.fits.gz',
                  'misc/abell370clu_misc/images/bcgs_models/abell370_f110w_bcgs_model.fits.gz',
                  'misc/abell370clu_misc/images/bcgs_models/abell370_f125w_bcgs_model.fits.gz',
                  'misc/abell370clu_misc/images/bcgs_models/abell370_f140w_bcgs_model.fits.gz']
        
        kernels = ['misc/abell370clu_misc/kernels/abell370_bcgs_out_f275w_f160w_kernel_69p.fits.gz',
                   'misc/abell370clu_misc/kernels/abell370_bcgs_out_f336w_f160w_kernel_69p.fits.gz',
                   'misc/abell370clu_misc/kernels/abell370_bcgs_out_f435w_f160w_kernel_69p.fits.gz',
                   'misc/abell370clu_misc/kernels/abell370_bcgs_out_f475w_f160w_kernel_69p.fits.gz',
                   'misc/abell370clu_misc/kernels/abell370_bcgs_out_f606w_f160w_kernel_69p.fits.gz',
                   'misc/abell370clu_misc/kernels/abell370_bcgs_out_f625w_f160w_kernel_69p.fits.gz',
                   'misc/abell370clu_misc/kernels/abell370_bcgs_out_f814w_f160w_kernel_69p.fits.gz',
                   'misc/abell370clu_misc/kernels/abell370_bcgs_out_f105w_f160w_kernel_69p.fits.gz',
                   'misc/abell370clu_misc/kernels/abell370_bcgs_out_f110w_f160w_kernel_69p.fits.gz',
                   'misc/abell370clu_misc/kernels/abell370_bcgs_out_f125w_f160w_kernel_69p.fits.gz',
                   'misc/abell370clu_misc/kernels/abell370_bcgs_out_f140w_f160w_kernel_69p.fits.gz']
        
        filters = ['f275w', 'f336w',
                   'f435w', 'f475w', 'f606w', 'f625w', 'f814w',
                   'f105w', 'f110w', 'f125w', 'f140w']
    
    if cluster == 'a1063' :
        models = ['misc/abell1063clu_misc/images/bcgs_models/abell1063_f275w_bcgs_model.fits.gz',
                  'misc/abell1063clu_misc/images/bcgs_models/abell1063_f336w_bcgs_model.fits.gz',
                  'misc/abell1063clu_misc/images/bcgs_models/abell1063_f390w_bcgs_model.fits.gz',
                  'misc/abell1063clu_misc/images/bcgs_models/abell1063_f435w_bcgs_model.fits.gz',
                  'misc/abell1063clu_misc/images/bcgs_models/abell1063_f475w_bcgs_model.fits.gz',
                  'misc/abell1063clu_misc/images/bcgs_models/abell1063_f606w_bcgs_model.fits.gz',
                  'misc/abell1063clu_misc/images/bcgs_models/abell1063_f625w_bcgs_model.fits.gz',
                  'misc/abell1063clu_misc/images/bcgs_models/abell1063_f775w_bcgs_model.fits.gz',
                  'misc/abell1063clu_misc/images/bcgs_models/abell1063_f814w_bcgs_model.fits.gz',
                  'misc/abell1063clu_misc/images/bcgs_models/abell1063_f850lp_bcgs_model.fits.gz',
                  'misc/abell1063clu_misc/images/bcgs_models/abell1063_f105w_bcgs_model.fits.gz',
                  'misc/abell1063clu_misc/images/bcgs_models/abell1063_f110w_bcgs_model.fits.gz',
                  'misc/abell1063clu_misc/images/bcgs_models/abell1063_f125w_bcgs_model.fits.gz',
                  'misc/abell1063clu_misc/images/bcgs_models/abell1063_f140w_bcgs_model.fits.gz']
        
        kernels = ['misc/abell1063clu_misc/kernels/abell1063_bcgs_out_f275w_f160w_kernel_69p.fits.gz',
                   'misc/abell1063clu_misc/kernels/abell1063_bcgs_out_f336w_f160w_kernel_69p.fits.gz',
                   'misc/abell1063clu_misc/kernels/abell1063_bcgs_out_f390w_f160w_kernel_69p.fits.gz',
                   'misc/abell1063clu_misc/kernels/abell1063_bcgs_out_f435w_f160w_kernel_69p.fits.gz',
                   'misc/abell1063clu_misc/kernels/abell1063_bcgs_out_f475w_f160w_kernel_69p.fits.gz',
                   'misc/abell1063clu_misc/kernels/abell1063_bcgs_out_f606w_f160w_kernel_69p.fits.gz',
                   'misc/abell1063clu_misc/kernels/abell1063_bcgs_out_f625w_f160w_kernel_69p.fits.gz',
                   'misc/abell1063clu_misc/kernels/abell1063_bcgs_out_f775w_f160w_kernel_69p.fits.gz',
                   'misc/abell1063clu_misc/kernels/abell1063_bcgs_out_f814w_f160w_kernel_69p.fits.gz',
                   'misc/abell1063clu_misc/kernels/abell1063_bcgs_out_f850lp_f160w_kernel_69p.fits.gz',
                   'misc/abell1063clu_misc/kernels/abell1063_bcgs_out_f105w_f160w_kernel_69p.fits.gz',
                   'misc/abell1063clu_misc/kernels/abell1063_bcgs_out_f110w_f160w_kernel_69p.fits.gz',
                   'misc/abell1063clu_misc/kernels/abell1063_bcgs_out_f125w_f160w_kernel_69p.fits.gz',
                   'misc/abell1063clu_misc/kernels/abell1063_bcgs_out_f140w_f160w_kernel_69p.fits.gz']
        
        filters = ['f275w', 'f336w', 'f390w',
                   'f435w', 'f475w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp',
                   'f105w', 'f110w', 'f125w', 'f140w']
    
    if cluster == 'a2744' :
        models = ['misc/abell2744clu_misc/images/bcgs_models/abell2744_f275w_bcgs_model.fits.gz',
                  'misc/abell2744clu_misc/images/bcgs_models/abell2744_f336w_bcgs_model.fits.gz',
                  'misc/abell2744clu_misc/images/bcgs_models/abell2744_f435w_bcgs_model.fits.gz',
                  'misc/abell2744clu_misc/images/bcgs_models/abell2744_f606w_bcgs_model.fits.gz',
                  'misc/abell2744clu_misc/images/bcgs_models/abell2744_f814w_bcgs_model.fits.gz',
                  'misc/abell2744clu_misc/images/bcgs_models/abell2744_f105w_bcgs_model.fits.gz',
                  'misc/abell2744clu_misc/images/bcgs_models/abell2744_f125w_bcgs_model.fits.gz',
                  'misc/abell2744clu_misc/images/bcgs_models/abell2744_f140w_bcgs_model.fits.gz']
        
        kernels = ['misc/abell2744clu_misc/kernels/abell2744_bcgs_out_f275w_f160w_kernel_69p.fits.gz',
                   'misc/abell2744clu_misc/kernels/abell2744_bcgs_out_f336w_f160w_kernel_69p.fits.gz',
                   'misc/abell2744clu_misc/kernels/abell2744_bcgs_out_f435w_f160w_kernel_69p.fits.gz',
                   'misc/abell2744clu_misc/kernels/abell2744_bcgs_out_f606w_f160w_kernel_69p.fits.gz',
                   'misc/abell2744clu_misc/kernels/abell2744_bcgs_out_f814w_f160w_kernel_69p.fits.gz',
                   'misc/abell2744clu_misc/kernels/abell2744_bcgs_out_f105w_f160w_kernel_69p.fits.gz',
                   'misc/abell2744clu_misc/kernels/abell2744_bcgs_out_f125w_f160w_kernel_69p.fits.gz',
                   'misc/abell2744clu_misc/kernels/abell2744_bcgs_out_f140w_f160w_kernel_69p.fits.gz']
        
        filters = ['f275w', 'f336w',
                   'f435w', 'f606w', 'f814w',
                   'f105w', 'f125w', 'f140w']
    
    if cluster == 'm416' :
        models = ['misc/macs0416clu_misc/images/bcgs_models/macs0416_f225w_bcgs_model.fits.gz',
                  'misc/macs0416clu_misc/images/bcgs_models/macs0416_f275w_bcgs_model.fits.gz',
                  'misc/macs0416clu_misc/images/bcgs_models/macs0416_f336w_bcgs_model.fits.gz',
                  'misc/macs0416clu_misc/images/bcgs_models/macs0416_f390w_bcgs_model.fits.gz',
                  'misc/macs0416clu_misc/images/bcgs_models/macs0416_f435w_bcgs_model.fits.gz',
                  'misc/macs0416clu_misc/images/bcgs_models/macs0416_f475w_bcgs_model.fits.gz',
                  'misc/macs0416clu_misc/images/bcgs_models/macs0416_f606w_bcgs_model.fits.gz',
                  'misc/macs0416clu_misc/images/bcgs_models/macs0416_f625w_bcgs_model.fits.gz',
                  'misc/macs0416clu_misc/images/bcgs_models/macs0416_f775w_bcgs_model.fits.gz',
                  'misc/macs0416clu_misc/images/bcgs_models/macs0416_f814w_bcgs_model.fits.gz',
                  'misc/macs0416clu_misc/images/bcgs_models/macs0416_f850lp_bcgs_model.fits.gz',
                  'misc/macs0416clu_misc/images/bcgs_models/macs0416_f105w_bcgs_model.fits.gz',
                  'misc/macs0416clu_misc/images/bcgs_models/macs0416_f110w_bcgs_model.fits.gz',
                  'misc/macs0416clu_misc/images/bcgs_models/macs0416_f125w_bcgs_model.fits.gz',
                  'misc/macs0416clu_misc/images/bcgs_models/macs0416_f140w_bcgs_model.fits.gz']
        
        kernels = ['misc/macs0416clu_misc/kernels/macs0416_bcgs_out_f225w_f160w_kernel_69p.fits.gz',
                   'misc/macs0416clu_misc/kernels/macs0416_bcgs_out_f275w_f160w_kernel_69p.fits.gz',
                   'misc/macs0416clu_misc/kernels/macs0416_bcgs_out_f336w_f160w_kernel_69p.fits.gz',
                   'misc/macs0416clu_misc/kernels/macs0416_bcgs_out_f390w_f160w_kernel_69p.fits.gz',
                   'misc/macs0416clu_misc/kernels/macs0416_bcgs_out_f435w_f160w_kernel_69p.fits.gz',
                   'misc/macs0416clu_misc/kernels/macs0416_bcgs_out_f475w_f160w_kernel_69p.fits.gz',
                   'misc/macs0416clu_misc/kernels/macs0416_bcgs_out_f606w_f160w_kernel_69p.fits.gz',
                   'misc/macs0416clu_misc/kernels/macs0416_bcgs_out_f625w_f160w_kernel_69p.fits.gz',
                   'misc/macs0416clu_misc/kernels/macs0416_bcgs_out_f775w_f160w_kernel_69p.fits.gz',
                   'misc/macs0416clu_misc/kernels/macs0416_bcgs_out_f814w_f160w_kernel_69p.fits.gz',
                   'misc/macs0416clu_misc/kernels/macs0416_bcgs_out_f850lp_f160w_kernel_69p.fits.gz',
                   'misc/macs0416clu_misc/kernels/macs0416_bcgs_out_f105w_f160w_kernel_69p.fits.gz',
                   'misc/macs0416clu_misc/kernels/macs0416_bcgs_out_f110w_f160w_kernel_69p.fits.gz',
                   'misc/macs0416clu_misc/kernels/macs0416_bcgs_out_f125w_f160w_kernel_69p.fits.gz',
                   'misc/macs0416clu_misc/kernels/macs0416_bcgs_out_f140w_f160w_kernel_69p.fits.gz']
        
        filters = ['f225w', 'f275w', 'f336w', 'f390w',
                   'f435w', 'f475w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp',
                   'f105w', 'f110w', 'f125w', 'f140w']
        
    if cluster == 'm717' :
        models = ['misc/macs0717clu_misc/images/bcgs_models/macs0717_f225w_bcgs_model.fits.gz',
                  'misc/macs0717clu_misc/images/bcgs_models/macs0717_f275w_bcgs_model.fits.gz',
                  'misc/macs0717clu_misc/images/bcgs_models/macs0717_f336w_bcgs_model.fits.gz',
                  'misc/macs0717clu_misc/images/bcgs_models/macs0717_f390w_bcgs_model.fits.gz',
                  'misc/macs0717clu_misc/images/bcgs_models/macs0717_f435w_bcgs_model.fits.gz',
                  'misc/macs0717clu_misc/images/bcgs_models/macs0717_f475w_bcgs_model.fits.gz',
                  'misc/macs0717clu_misc/images/bcgs_models/macs0717_f555w_bcgs_model.fits.gz',
                  'misc/macs0717clu_misc/images/bcgs_models/macs0717_f606w_bcgs_model.fits.gz',
                  'misc/macs0717clu_misc/images/bcgs_models/macs0717_f625w_bcgs_model.fits.gz',
                  'misc/macs0717clu_misc/images/bcgs_models/macs0717_f775w_bcgs_model.fits.gz',
                  'misc/macs0717clu_misc/images/bcgs_models/macs0717_f814w_bcgs_model.fits.gz',
                  'misc/macs0717clu_misc/images/bcgs_models/macs0717_f850lp_bcgs_model.fits.gz',
                  'misc/macs0717clu_misc/images/bcgs_models/macs0717_f105w_bcgs_model.fits.gz',
                  'misc/macs0717clu_misc/images/bcgs_models/macs0717_f110w_bcgs_model.fits.gz',
                  'misc/macs0717clu_misc/images/bcgs_models/macs0717_f125w_bcgs_model.fits.gz',
                  'misc/macs0717clu_misc/images/bcgs_models/macs0717_f140w_bcgs_model.fits.gz']
        
        kernels = ['misc/macs0717clu_misc/kernels/macs0717_bcgs_out_f225w_f160w_kernel_69p.fits.gz',
                   'misc/macs0717clu_misc/kernels/macs0717_bcgs_out_f275w_f160w_kernel_69p.fits.gz',
                   'misc/macs0717clu_misc/kernels/macs0717_bcgs_out_f336w_f160w_kernel_69p.fits.gz',
                   'misc/macs0717clu_misc/kernels/macs0717_bcgs_out_f390w_f160w_kernel_69p.fits.gz',
                   'misc/macs0717clu_misc/kernels/macs0717_bcgs_out_f435w_f160w_kernel_69p.fits.gz',
                   'misc/macs0717clu_misc/kernels/macs0717_bcgs_out_f475w_f160w_kernel_69p.fits.gz',
                   'misc/macs0717clu_misc/kernels/macs0717_bcgs_out_f555w_f160w_kernel_69p.fits.gz',
                   'misc/macs0717clu_misc/kernels/macs0717_bcgs_out_f606w_f160w_kernel_69p.fits.gz',
                   'misc/macs0717clu_misc/kernels/macs0717_bcgs_out_f625w_f160w_kernel_69p.fits.gz',
                   'misc/macs0717clu_misc/kernels/macs0717_bcgs_out_f775w_f160w_kernel_69p.fits.gz',
                   'misc/macs0717clu_misc/kernels/macs0717_bcgs_out_f814w_f160w_kernel_69p.fits.gz',
                   'misc/macs0717clu_misc/kernels/macs0717_bcgs_out_f850lp_f160w_kernel_69p.fits.gz',
                   'misc/macs0717clu_misc/kernels/macs0717_bcgs_out_f105w_f160w_kernel_69p.fits.gz',
                   'misc/macs0717clu_misc/kernels/macs0717_bcgs_out_f110w_f160w_kernel_69p.fits.gz',
                   'misc/macs0717clu_misc/kernels/macs0717_bcgs_out_f125w_f160w_kernel_69p.fits.gz',
                   'misc/macs0717clu_misc/kernels/macs0717_bcgs_out_f140w_f160w_kernel_69p.fits.gz']
        
        filters = ['f225w', 'f275w', 'f336w', 'f390w',
                   'f435w', 'f475w', 'f555w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp',
                   'f105w', 'f110w', 'f125w', 'f140w']
    
    if cluster == 'm1149' :
        models = ['misc/macs1149clu_misc/images/bcgs_models/macs1149_f336w_bcgs_model.fits.gz',
                  'misc/macs1149clu_misc/images/bcgs_models/macs1149_f390w_bcgs_model.fits.gz',
                  'misc/macs1149clu_misc/images/bcgs_models/macs1149_f435w_bcgs_model.fits.gz',
                  'misc/macs1149clu_misc/images/bcgs_models/macs1149_f475w_bcgs_model.fits.gz',
                  'misc/macs1149clu_misc/images/bcgs_models/macs1149_f555w_bcgs_model.fits.gz',
                  'misc/macs1149clu_misc/images/bcgs_models/macs1149_f606w_bcgs_model.fits.gz',
                  'misc/macs1149clu_misc/images/bcgs_models/macs1149_f625w_bcgs_model.fits.gz',
                  'misc/macs1149clu_misc/images/bcgs_models/macs1149_f775w_bcgs_model.fits.gz',
                  'misc/macs1149clu_misc/images/bcgs_models/macs1149_f814w_bcgs_model.fits.gz',
                  'misc/macs1149clu_misc/images/bcgs_models/macs1149_f850lp_bcgs_model.fits.gz',
                  'misc/macs1149clu_misc/images/bcgs_models/macs1149_f105w_bcgs_model.fits.gz',
                  'misc/macs1149clu_misc/images/bcgs_models/macs1149_f110w_bcgs_model.fits.gz',
                  'misc/macs1149clu_misc/images/bcgs_models/macs1149_f125w_bcgs_model.fits.gz',
                  'misc/macs1149clu_misc/images/bcgs_models/macs1149_f140w_bcgs_model.fits.gz']
        
        kernels = ['misc/macs1149clu_misc/kernels/macs1149_bcgs_out_f336w_f160w_kernel_69p.fits.gz',
                   'misc/macs1149clu_misc/kernels/macs1149_bcgs_out_f390w_f160w_kernel_69p.fits.gz',
                   'misc/macs1149clu_misc/kernels/macs1149_bcgs_out_f435w_f160w_kernel_69p.fits.gz',
                   'misc/macs1149clu_misc/kernels/macs1149_bcgs_out_f475w_f160w_kernel_69p.fits.gz',
                   'misc/macs1149clu_misc/kernels/macs1149_bcgs_out_f555w_f160w_kernel_69p.fits.gz',
                   'misc/macs1149clu_misc/kernels/macs1149_bcgs_out_f606w_f160w_kernel_69p.fits.gz',
                   'misc/macs1149clu_misc/kernels/macs1149_bcgs_out_f625w_f160w_kernel_69p.fits.gz',
                   'misc/macs1149clu_misc/kernels/macs1149_bcgs_out_f775w_f160w_kernel_69p.fits.gz',
                   'misc/macs1149clu_misc/kernels/macs1149_bcgs_out_f814w_f160w_kernel_69p.fits.gz',
                   'misc/macs1149clu_misc/kernels/macs1149_bcgs_out_f850lp_f160w_kernel_69p.fits.gz',
                   'misc/macs1149clu_misc/kernels/macs1149_bcgs_out_f105w_f160w_kernel_69p.fits.gz',
                   'misc/macs1149clu_misc/kernels/macs1149_bcgs_out_f110w_f160w_kernel_69p.fits.gz',
                   'misc/macs1149clu_misc/kernels/macs1149_bcgs_out_f125w_f160w_kernel_69p.fits.gz',
                   'misc/macs1149clu_misc/kernels/macs1149_bcgs_out_f140w_f160w_kernel_69p.fits.gz']
        
        filters = ['f336w', 'f390w',
                   'f435w', 'f475w', 'f555w', 'f606w', 'f625w', 'f775w', 'f814w', 'f850lp',
                   'f105w', 'f110w', 'f125w', 'f140w']
    
    if cluster == 'a2744par' :
        models = ['misc/abell2744par_misc/images/bcgs_models/abell2744par_f435w_bcgs_model.fits.gz',
                  'misc/abell2744par_misc/images/bcgs_models/abell2744par_f606w_bcgs_model.fits.gz',
                  'misc/abell2744par_misc/images/bcgs_models/abell2744par_f814w_bcgs_model.fits.gz',
                  'misc/abell2744par_misc/images/bcgs_models/abell2744par_f105w_bcgs_model.fits.gz',
                  'misc/abell2744par_misc/images/bcgs_models/abell2744par_f125w_bcgs_model.fits.gz',
                  'misc/abell2744par_misc/images/bcgs_models/abell2744par_f140w_bcgs_model.fits.gz']
        
        kernels = ['misc/abell2744par_misc/kernels/abell2744par_bcgs_out_f435w_f160w_kernel_69p.fits.gz',
                   'misc/abell2744par_misc/kernels/abell2744par_bcgs_out_f606w_f160w_kernel_69p.fits.gz',
                   'misc/abell2744par_misc/kernels/abell2744par_bcgs_out_f814w_f160w_kernel_69p.fits.gz',
                   'misc/abell2744par_misc/kernels/abell2744par_bcgs_out_f105w_f160w_kernel_69p.fits.gz',
                   'misc/abell2744par_misc/kernels/abell2744par_bcgs_out_f125w_f160w_kernel_69p.fits.gz',
                   'misc/abell2744par_misc/kernels/abell2744par_bcgs_out_f140w_f160w_kernel_69p.fits.gz']
        
        filters = ['f435w', 'f606w', 'f814w',
                   'f105w', 'f125w', 'f140w']
    
    if cluster == 'm416par' :
        models = ['misc/macs0416par_misc/images/bcgs_models/macs0416par_f435w_bcgs_model.fits.gz',
                  'misc/macs0416par_misc/images/bcgs_models/macs0416par_f606w_bcgs_model.fits.gz',
                  'misc/macs0416par_misc/images/bcgs_models/macs0416par_f775w_bcgs_model.fits.gz',
                  'misc/macs0416par_misc/images/bcgs_models/macs0416par_f814w_bcgs_model.fits.gz',
                  'misc/macs0416par_misc/images/bcgs_models/macs0416par_f850lp_bcgs_model.fits.gz',
                  'misc/macs0416par_misc/images/bcgs_models/macs0416par_f105w_bcgs_model.fits.gz',
                  'misc/macs0416par_misc/images/bcgs_models/macs0416par_f125w_bcgs_model.fits.gz',
                  'misc/macs0416par_misc/images/bcgs_models/macs0416par_f140w_bcgs_model.fits.gz']
        
        kernels = ['misc/macs0416par_misc/kernels/macs0416par_bcgs_out_f435w_f160w_kernel_69p.fits.gz',
                   'misc/macs0416par_misc/kernels/macs0416par_bcgs_out_f606w_f160w_kernel_69p.fits.gz',
                   'misc/macs0416par_misc/kernels/macs0416par_bcgs_out_f775w_f160w_kernel_69p.fits.gz',
                   'misc/macs0416par_misc/kernels/macs0416par_bcgs_out_f814w_f160w_kernel_69p.fits.gz',
                   'misc/macs0416par_misc/kernels/macs0416par_bcgs_out_f850lp_f160w_kernel_69p.fits.gz',
                   'misc/macs0416par_misc/kernels/macs0416par_bcgs_out_f105w_f160w_kernel_69p.fits.gz',
                   'misc/macs0416par_misc/kernels/macs0416par_bcgs_out_f125w_f160w_kernel_69p.fits.gz',
                   'misc/macs0416par_misc/kernels/macs0416par_bcgs_out_f140w_f160w_kernel_69p.fits.gz']
        
        filters = ['f435w', 'f606w', 'f775w', 'f814w', 'f850lp',
                   'f105w', 'f125w', 'f140w']
    
    if cluster == 'm717par' :
        models = ['misc/macs0717par_misc/images/bcgs_models/macs0717par_f435w_bcgs_model.fits.gz',
                  'misc/macs0717par_misc/images/bcgs_models/macs0717par_f606w_bcgs_model.fits.gz',
                  'misc/macs0717par_misc/images/bcgs_models/macs0717par_f814w_bcgs_model.fits.gz',
                  'misc/macs0717par_misc/images/bcgs_models/macs0717par_f105w_bcgs_model.fits.gz',
                  'misc/macs0717par_misc/images/bcgs_models/macs0717par_f125w_bcgs_model.fits.gz',
                  'misc/macs0717par_misc/images/bcgs_models/macs0717par_f140w_bcgs_model.fits.gz']
        
        kernels = ['misc/macs0717par_misc/kernels/macs0717par_bcgs_out_f435w_f160w_kernel_69p.fits.gz',
                   'misc/macs0717par_misc/kernels/macs0717par_bcgs_out_f606w_f160w_kernel_69p.fits.gz',
                   'misc/macs0717par_misc/kernels/macs0717par_bcgs_out_f814w_f160w_kernel_69p.fits.gz',
                   'misc/macs0717par_misc/kernels/macs0717par_bcgs_out_f105w_f160w_kernel_69p.fits.gz',
                   'misc/macs0717par_misc/kernels/macs0717par_bcgs_out_f125w_f160w_kernel_69p.fits.gz',
                   'misc/macs0717par_misc/kernels/macs0717par_bcgs_out_f140w_f160w_kernel_69p.fits.gz']
        
        filters = ['f435w', 'f606w', 'f814w',
                   'f105w', 'f125w', 'f140w']
    
    if cluster == 'm1149par' :
        models = ['misc/macs1149par_misc/images/bcgs_models/macs1149par_f435w_bcgs_model.fits.gz',
                  'misc/macs1149par_misc/images/bcgs_models/macs1149par_f606w_bcgs_model.fits.gz',
                  'misc/macs1149par_misc/images/bcgs_models/macs1149par_f814w_bcgs_model.fits.gz',
                  'misc/macs1149par_misc/images/bcgs_models/macs1149par_f105w_bcgs_model.fits.gz',
                  'misc/macs1149par_misc/images/bcgs_models/macs1149par_f125w_bcgs_model.fits.gz',
                  'misc/macs1149par_misc/images/bcgs_models/macs1149par_f140w_bcgs_model.fits.gz']
        
        kernels = ['misc/macs1149par_misc/kernels/macs1149par_bcgs_out_f435w_f160w_kernel_69p.fits.gz',
                   'misc/macs1149par_misc/kernels/macs1149par_bcgs_out_f606w_f160w_kernel_69p.fits.gz',
                   'misc/macs1149par_misc/kernels/macs1149par_bcgs_out_f814w_f160w_kernel_69p.fits.gz',
                   'misc/macs1149par_misc/kernels/macs1149par_bcgs_out_f105w_f160w_kernel_69p.fits.gz',
                   'misc/macs1149par_misc/kernels/macs1149par_bcgs_out_f125w_f160w_kernel_69p.fits.gz',
                   'misc/macs1149par_misc/kernels/macs1149par_bcgs_out_f140w_f160w_kernel_69p.fits.gz']
        
        filters = ['f435w', 'f606w', 'f814w',
                   'f105w', 'f125w', 'f140w']
    
    for model, kernel, filt in zip(models, kernels, filters) :
        convolve_with_kernel(model, kernel,
                              'models/{}_{}_bcg_model.fits.gz'.format(
                                  cluster, filt))
    
    return

def convolve_with_kernel(image, kernel, outfile) :
    
    with fits.open(image) as hdu :
        hdr = hdu[0].header
        data = hdu[0].data
    
    with fits.open(kernel) as hduk :
        kern = hduk[0].data
    
    hdu = fits.PrimaryHDU(convolve(data, kern))
    hdu.header = hdr
    hdu.writeto(outfile)
    
    return
