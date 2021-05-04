
from astropy.io import fits

import initial

def combine(residual, model, outfile) :
    
    with fits.open(residual) as hdur :
        res_hdr = hdur[0].header
        res_data = hdur[0].data
    
    with fits.open(model) as hdum :
        model_data = hdum[0].data
    
    fits.writeto(outfile, res_data + model_data, header=res_hdr)
    
    return

def combine_cluster(residuals, models) :
    
    outfiles = []
    for model in models :
        cluster, filt, _, _ = model.split('/')[-1].split('_')
        outfile = '{}/{}_{}_res+mod.fits'.format('detec', cluster, filt)
        outfiles.append(outfile)
    
    for i in range(len(residuals)) :
        combine(residuals[i], models[i], outfiles[i])
    
    return

def combine_all() :
    
    (_, _, _, _, _, _, _, _,
     a370_files, a370_models, _) = initial.a370_params()
    (_, _, _, _, _, _, _, _,
     a1063_files, a1063_models, _) = initial.a1063_params()
    (_, _, _, _, _, _, _, _,
     a2744_files, a2744_models, _) = initial.a2744_params()
    (_, _, _, _, _, _, _, _,
     m416_files, m416_models, _) = initial.m416_params()
    (_, _, _, _, _, _, _, _,
     m717_files, m717_models, _) = initial.m717_params()
    (_, _, _, _, _, _, _, _,
     m1149_files, m1149_models, _) = initial.m1149_params()
    
    combine_cluster(a370_files, a370_models)
    combine_cluster(a1063_files, a1063_models)
    combine_cluster(a2744_files, a2744_models)
    combine_cluster(m416_files, m416_models)
    combine_cluster(m717_files, m717_models)
    combine_cluster(m1149_files, m1149_models)
    
    return
