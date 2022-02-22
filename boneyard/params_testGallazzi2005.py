
def build_all(**extras) :
    return (build_obs(**extras), build_model(**extras),
            build_sps(**extras), build_noise(**extras))

def build_model(binNum=0, fit_metallicity=1, fit_redshift=1, infile=None,
                add_neb=0, **extras) :
    
    from astropy.table import Table
    from prospect.models.priors import LogUniform, TopHat, Normal
    from prospect.models.sedmodel import SedModel
    from prospect.models.templates import TemplateLibrary
    from prospect.sources.constants import cosmo
    
    table = Table.read(infile)
    zobs = table['z'][binNum]
    max_age = cosmo.age(zobs).value
    
    model_params = TemplateLibrary['parametric_sfh'] # delay-tau model
    
    model_params['zred']['init'] = zobs
    model_params['zred']['isfree'] = bool(fit_redshift)
    model_params['zred']['prior'] = TopHat(mini=zobs-0.01, maxi=zobs+0.01)
    
    model_params['mass']['init'] = 1e7
    model_params['mass']['prior'] = LogUniform(mini=1e5, maxi=1e9)
    
    model_params['logzsol']['init'] = -0.15
    model_params['logzsol']['isfree'] = bool(fit_metallicity)
    model_params['logzsol']['prior'] = Normal(mean=-0.15, sigma=0.22)
    
    model_params['dust2']['init'] = 0.3
    model_params['dust2']['prior'] = TopHat(mini=0.0, maxi=2.0)
    
    model_params['tage']['init'] = 6.24
    model_params['tage']['prior'] = TopHat(mini=0.001, maxi=max_age)
    
    model_params['tau']['init'] = 0.29
    model_params['tau']['prior'] = LogUniform(mini=0.1, maxi=10)
    
    model_params['imf_type']['init'] = 1 # Chabrier (2003) IMF
    model_params['dust_type']['init'] = 1 # Cardelli+ (1989) MW extinction
    
    if add_neb :
        model_params.update(TemplateLibrary['nebular'])
    
    model = SedModel(model_params)
    
    return model

def build_noise(**extras) :
    return None, None

def build_obs(binNum=0, infile=None, **extras) :
    
    import numpy as np
    from astropy.table import Table
    from sedpy.observate import load_filters
    from prospect.utils.obsutils import fix_obs
    
    table = Table.read(infile)
    
    flux_columns = [col for col in table.colnames if col.endswith('_flux')]
    e_flux_columns = [col.replace('flux', 'err') for col in flux_columns]
    nPix_columns = [col.replace('flux', 'nPix') for col in flux_columns]
    use_columns = [col.replace('flux', 'use') for col in flux_columns]
    filternames = ['hff_' + col.replace('_flux', '') for col in flux_columns]
    
    fluxes = np.array(list(table[flux_columns][binNum]))/3631 # in maggies
    e_fluxes = np.array(list(table[e_flux_columns][binNum]))/3631 # in maggies
    nPixels = np.array(list(table[nPix_columns][binNum]))
    use_col = np.array(list(table[use_columns][binNum]))
    
    fluxes_per_pixel = fluxes/nPixels
    e_fluxes_per_pixel = e_fluxes/nPixels
    
    max_SNR = 20
    SNR_acceptable = (fluxes_per_pixel/e_fluxes_per_pixel < max_SNR)
    e_fluxes_per_pixel[~SNR_acceptable] = fluxes_per_pixel[~SNR_acceptable]/max_SNR
    
    obs = {}
    
    obs['filters'] = load_filters(filternames)
    obs['maggies'], obs['maggies_unc'] = fluxes_per_pixel, e_fluxes_per_pixel
    
    obs['phot_mask'] =  use_col
    obs['phot_wave'] = np.array([f.wave_effective for f in obs['filters']])
    
    obs['wavelength'], obs['spectrum'] = None, None
    obs['unc'], obs['mask'] = None, None
    
    obs['binNum'], obs['nPixels'] = binNum, nPixels
    
    obs = fix_obs(obs)
    
    return obs

def build_sps(**extras) :
    from prospect.sources import CSPSpecBasis
    sps = CSPSpecBasis(zcontinuous=1, compute_vega_mags=False)
    return sps

if __name__ == '__main__' :
    
    import sys
    from prospect import prospect_args
    from prospect.fitting import fit_model
    from prospect.io import write_results as writer
    
    parser = prospect_args.get_parser()
    
    parser.add_argument('--add_neb', type=int, default=0,
                        help='Flag to add nebular emission.')
    parser.add_argument('--binNum', type=int, default=0,
                        help='Bin/annulus of the galaxy to fit.')
    parser.add_argument('--fit_metallicity', type=int, default=1,
                        help='Flag to fit for the metallicity.')
    parser.add_argument('--fit_redshift', type=int, default=1,
                        help='Flag to fit for the redshift.')
    parser.add_argument('--infile', type=str,
                        help='File to get the photometry from.')
    
    args = parser.parse_args()
    run_params = vars(args)
    
    run_params['dynesty'] = True
    binNum = run_params['binNum']
    cluster, _, ID, _ = run_params['infile'].split('/')[-1].split('_')
    outfile = '{}/h5/{}_ID_{}_bin_{}_gallazzi'.format(cluster, cluster, ID, binNum)
    run_params['outfile'] = outfile
    
    obs, model, sps, noise = build_all(**run_params)
    
    run_params['sps_libraries'] = sps.ssp.libraries
    run_params['param_file'] = __file__
    
    if args.debug :
        sys.exit()
    
    output = fit_model(obs, model, sps, noise, **run_params)
    
    hfile = '{}.h5'.format(outfile)
    writer.write_hdf5(hfile, run_params, model, obs,
                      output['sampling'][0], output['optimization'][0],
                      tsample=output['sampling'][1],
                      toptimize=output['optimization'][1],
                      sps=sps)
    try :
        hfile.close()
    except(AttributeError) :
        pass
