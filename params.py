
def build_all(**kwargs) :
    '''
    Construct all pieces required for running a prospector fit.
    
    Parameters
    ----------
    **kwargs : kwargs
        Additional keyword arguments.
    
    Returns
    -------
    tuple
        The required pieces for running.
    
    '''
    
    return (build_obs(**kwargs), build_model(**kwargs),
            build_sps(**kwargs), build_noise(**kwargs))

def build_model(object_redshift=0.0, fixed_metallicity=None, add_duste=False,
                add_neb=False, luminosity_distance=0.0, **extras) :
    '''
    Construct a model. This method defines a number of parameter specification
    dictionaries and uses them to initialize a `models.sedmodel.SedModel`
    object.
    
    Parameters
    ----------
    object_redshift : float, optional
        If given, the model uses this redshift value. The default is 0.0.
    fixed_metallicity : bool, optional
        Boolean to set the metallicity to a fixed value. The default is None.
    add_duste : bool, optional
        Boolean to add (fixed) parameters relevant for dust emission. The
        default is False.
    add_neb : bool, optional
        Boolean to add (fixed) parameters relevant for nebular emission, and
        to turn nebular emission on. The default is False.
    luminosity_distance : float, optional
        If present, add a `lumdist` parameter to the model, and set it's value
        (in Mpc) to this. This allows one to decouple redshift from distance,
        and fit, for example, absolute magnitudes. The default is 0.0.
    **extras : kwargs
        Additional keyword arguments.
    
    Returns
    -------
    model : models.sedmodel.SedModel
        The model to use when fitting.
    
    '''
    
    from prospect.models.templates import TemplateLibrary
    from prospect.models import sedmodel
    
    # Get a basic non-parametric SFH parameter set.
    # This has 4 free parameters:
    #   'logmass', 'logzsol', 'dust2', 'logsfr_ratios'
    # and two fixed parameters:
    #   'zred'=0.1, 'mass'=1e6, 'sfh'=3, 'imf_type'=2, 'dust_type'=0,
    #   'agebins'=[[0.0, 8.0], [8.0, 9.0], [9.0, 10.0]]
    # See the python-FSPS documentation for details about most of these
    # parameters. Also look at `TemplateLibrary.describe('continuity_sfh')` to
    # view the parameters, their initial values, and the priors in detail.
    model_params = TemplateLibrary['continuity_sfh']
    
    # add lumdist parameter. If this is not added then the distance is
    # controlled by the 'zred' parameter and a WMAP9 cosmology.
    if luminosity_distance > 0 :
        model_params['lumdist'] = {'N': 1, 'isfree': False,
                                   'init': luminosity_distance, 'units':'Mpc'}
    
    # change the model parameter specifications based on some keyword arguments
    if fixed_metallicity is not None :
        # make it a fixed parameter
        model_params['logzsol']['isfree'] = False
        # and use value supplied by fixed_metallicity keyword
        model_params['logzsol']['init'] = fixed_metallicity
    
    if object_redshift != 0.0 :
        # make sure zred is fixed
        model_params['zred']['isfree'] = False
        # and set the value to the object_redshift keyword
        model_params['zred']['init'] = object_redshift
    
    if add_duste :
        # add dust emission (with fixed dust SED parameters)
        model_params.update(TemplateLibrary['dust_emission'])
    
    if add_neb :
        # add nebular emission (with fixed parameters)
        model_params.update(TemplateLibrary['nebular'])
    
    # now instantiate the model using this new dictionary of parameter specs
    model = sedmodel.SedModel(model_params)
    
    return model

def build_noise(**extras) :
    '''
    Construct a noise model.
    
    Parameters
    ----------
    **extras : kwargs
        Additional keyword arguments.
    
    Returns
    -------
    tuple
        Empty tuple containing two NoneType objects.
    
    '''
    
    return None, None

def build_obs(infile=None, filterFile=None, vorbinNum=0, **kwargs) :
    '''
    Construct the observations to use for fitting.
    
    Parameters
    ----------
    fluxes : list, optional
        The fluxes to use. The default is None.
    errs : list, optional
        The uncertainties on the above fluxes. The default is None.
    filterset : list, optional
        The filters that the observations were taken in. The default is None.
    vorbinNum : int, optional
        The vorbin number for the galaxy. The default is 0.
    **kwargs : kwargs
        Additional keyword arguments.
        
    Returns
    -------
    obs : dict
        Dictionary of observational data.
    
    '''
    
    import numpy as np
    
    from astropy.table import Table
    from sedpy.observate import load_filters
    from prospect.utils.obsutils import fix_obs
    
    with open(filterFile) as file : # open the file containing strings of the
        filters = []                # filters and create a list containing
        for line in file :          # those strings
            filt = line.strip()
            filters.append(filt)
    
    table = Table.read(infile)
    flux_cols = [string for string in table.colnames if '_flux' in string]
    err_cols = [string for string in table.colnames if '_err' in string]
    flux_table = table[flux_cols] # create a table using only the flux columns
    err_table = table[err_cols] # create a table using only the err columns
    
    fluxes = list(flux_table[vorbinNum]) # get the fluxes for the vorbin
    errs = list(err_table[vorbinNum]) # get the uncertainties for the vorbin
    
    # build output dictionary.
    obs = {}
    # this is a list of sedpy filter objects. See the
    # sedpy.observate.load_filters command for more details on its syntax.
    obs['filters'] = load_filters(filters)
    # this is a list of maggies, converted from fluxes. It should have the same
    # order as `filters` above.
    obs['maggies'] = np.array(fluxes)/3631
    # this is a list of maggie uncertainties, converted from fluxes. It should
    # have the same order as `filters` above.
    obs['maggies_unc'] = np.array(errs)/3631
    # here we mask out any NaNs or infs
    obs['phot_mask'] = np.array([True]*len(filters))
    # we have no spectrum.
    obs['wavelength'] = None
    obs['spectrum'] = None
    obs['unc'] = None
    obs['mask'] = None
    
    # add unessential bonus info. This will be stored in output
    obs['vorbinNum'] = vorbinNum
    
    # this ensures all required keys are present and adds extra useful info
    obs = fix_obs(obs)
    
    return obs

def build_sps(zcontinuous=1, compute_vega_mags=False, **extras) :
    '''
    Instantiate and return the Stellar Population Synthesis object.
    
    Parameters
    ----------
    zcontinuous : float, optional
        A value of 1 ensures that interpolation is used between SSPs to have a
        continuous metallicity parameter (`logzsol`). See python-FSPS
        documentation for details. The default is 1.
    compute_vega_mags : bool, optional
        Boolean that sets the zero points of the magnitude system. The default
        is False.
    **extras : kwargs
        Additional keyword arguments.
    
    Returns
    -------
    sps : sources.ssp_basis.FastStepBasis
        Stellar population synthesis object to use for fitting.
    
    '''
    
    from prospect.sources import FastStepBasis
    sps = FastStepBasis(zcontinuous=zcontinuous,
                        compute_vega_mags=compute_vega_mags)
    return sps

if __name__ == '__main__' :
    
    import sys
    
    from prospect import prospect_args
    from prospect.fitting import fit_model
    from prospect.io import write_results as writer
    
    # parser with default arguments
    parser = prospect_args.get_parser()
    
    # add custom arguments
    parser.add_argument('--object_redshift', type=float, default=0.0,
                        help='Redshift for the model.')
    parser.add_argument('--add_neb', action='store_true', default=False,
                        help='If set, add nebular emission in the model.')
    parser.add_argument('--add_duste', action='store_true', default=False,
                        help='If set, add dust emission to the model.')
    parser.add_argument('--luminosity_distance', type=float,
                        help='Luminosity distance in Mpc.')
    
    parser.add_argument('--infile', type=str,
                        help='File to get the photometry from.')
    parser.add_argument('--filterFile', type=str,
                        help='File containing filterset information.')
    parser.add_argument('--vorbinNum', type=int,
                        help='Vorbin number for the galaxy.')
    
    args = parser.parse_args()
    run_params = vars(args)
    
    obs, model, sps, noise = build_all(**run_params)
    
    run_params['sps_libraries'] = sps.ssp.libraries
    run_params['param_file'] = __file__
    
    if args.debug :
        sys.exit()
    
    hfile = '{}_dynesty.h5'.format(args.outfile)
    output = fit_model(obs, model, sps, noise, **run_params)
    
    writer.write_hdf5(hfile, run_params, model, obs,
                      output['sampling'][0], output['optimization'][0],
                      tsample=output['sampling'][1],
                      toptimize=output['optimization'][1],
                      sps=sps)
    
    try :
        hfile.close()
    except(AttributeError) :
        pass
    