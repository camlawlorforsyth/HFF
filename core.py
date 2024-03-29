
import os
import numpy as np

import astropy.constants as const
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.table import hstack, join, Table
import astropy.units as u
from scipy import interpolate

import plotting as plt

def build_par_file_WFC3(filt, wave, throughput) :
    
    # sedpy needs evenly spaced wavelength arrays and corresponding throughputs
    interpolated = interpolate.interp1d(wave, throughput)
    
    # evenly sample the wavelengths and use interpolated throughput values
    sampling_waves = np.arange(wave[0], wave[-1]+1, 1)
    interpolated_throughput = interpolated(sampling_waves)
    
    # write the corresponding tables
    par_table = Table([sampling_waves, interpolated_throughput],
                      names=('wave', 'thru'))
    par_table.write('hff_filters/hff_{}.par'.format(filt),
                    format='ascii.no_header', overwrite=False)
    
    return

def build_par_file_ACS(table) :
    
    cols = table.colnames
    
    # loop over all the column pairs in the table
    for i in range(0, len(cols), 2) :
        filt = cols[i].split('_')[0]
        
        # get the non-masked entries
        wave = table[cols[i]]
        wave = np.array(wave[wave > 0])
        throughput = table[cols[i+1]]
        throughput = np.array(throughput[throughput > 0])
        
        # interpolate values and write files
        build_par_file_WFC3(filt, np.float64(wave), np.float64(throughput))
    
    return

def build_filter_table(save_par=False, display=False, write=False) :
    '''
    Build the filter table which includes wavelength arrays and transmission
    arrays for all the filters included in the Hubble Frontier Fields
    DeepSpace program. This includes 17 filters in total. The initial data
    comes from various CSV files.
    
    Parameters
    ----------
    display : bool, optional
        Boolean to display the master table. The default is False.
    write : bool, optional
        Boolean to write the master table to a file. The default is False.
    
    Returns
    -------
    master_table : astropy.table.Table
        The final table which includes all wavelength and transmission arrays.
    
    '''
    
    if save_par :
        os.makedirs('hff_filters', exist_ok=True) # ensure the output directory
        # for the filter throughput files is available
    
    # HST WFC3/UVIS (Wide Field Camera 3 - UV/Visible channel)
    # F225W, F275W, F336W, F390W
    # (UV wide, UV wide, U/Stromgren u, Washington C)
    UVIS = ['f225w_UVIS_throughput.csv', 'f275w_UVIS_throughput.csv',
            'f336w_UVIS_throughput.csv', 'f390w_UVIS_throughput.csv']
    UVIS_tables = []
    for file in UVIS :
        table = Table.read('hst_filters/{}'.format(file), format='csv')
        filt = file[:5]
        wave1 = table['Chip 1 Wave (Angstroms)']
        wave2 = table['Chip 2 Wave (Angstroms)']
        wave = (wave1 + wave2)/2
        thru1 = table['Chip 1 Throughput']
        thru2 = table['Chip 2 Throughput']
        thru = (thru1 + thru2)/2
        
        if save_par :
            build_par_file_WFC3(filt, wave, thru)
        
        filt_table = Table([wave, thru], names=(filt+'_wave', filt+'_thru'))
        UVIS_tables.append(filt_table)
    UVIS_table = hstack(UVIS_tables)
    
    # HST ACS/WFC (Advanced Camera for Surveys - Wide Field Channel)
    # F435W, F475W, F555W, F606W, F625W, F775W, F814W, F850LP
    # (Johnson B, SDSS g, Johnson V, Broad V, SDSS r, SDSS i, Broad I, SDSS z)
    ACS_table = Table.read('hst_filters/ACS_WFC_throughputs.csv', format='csv')
    if save_par :
        build_par_file_ACS(ACS_table)
    
    # HST WFC3/IR (Wide Field Camera 3 - near-IR channel)
    # F105W, F110W, F125W, F140W, F160W
    # (wide Y, wide YJ, wide J, wide JH gap, H)
    IR = ['f105w_IR_throughput.csv', 'f110w_IR_throughput.csv',
          'f125w_IR_throughput.csv', 'f140w_IR_throughput.csv',
          'f160w_IR_throughput.csv']
    IR_tables = []
    for file in IR :
        table = Table.read('hst_filters/{}'.format(file), format='csv')
        filt = file[:5]
        wave = table['Wave (Angstroms)']
        thru = table['Throughput']
        
        if save_par :
            build_par_file_WFC3(filt, wave, thru)
        
        filt_table = Table([wave, thru], names=(filt+'_wave', filt+'_thru'))
        IR_tables.append(filt_table)
    IR_table = hstack(IR_tables)
    
    master_table = hstack([UVIS_table, ACS_table, IR_table])
    
    if display :
        print(master_table)
    
    if write :
        master_table.write('filter_throughputs.fits', overwrite=True)
    
    return master_table

def combine_colors(first_path, second_path, FUV_filtnum=218, U_filtnum=153,
                   V_filtnum=155, J_filtnum=161, cluster='', key='id',
                   plot_hist=False, write=False) :
    '''
    
    
    Parameters
    ----------
    first_path : string
        DESCRIPTION.
    second_path : string
        DESCRIPTION.
    FUV_filtnum : int, optional
        Number describing FUV filter color indices file. The default is 218.
    U_filtnum : int, optional
        Number describing U filter color indices file. The default is 153.
    V_filtnum : int, optional
        Number describing V filter color indices file. The default is 155.
    J_filtnum : int, optional
        Number describing J filter color indices file. The default is 161.
    key : string, optional
        The key to base the table join operation on. The default is 'id'.
    plotuvj : bool, optional
        Boolean to plot the values on a scatter plot. The default is False.
    
    Returns
    -------
    table : astropy.table.Table
        DESCRIPTION.
    
    '''
    
    direc = 'catalogs/{}/eazy/{}'.format(first_path, second_path)
    FUV_file = '{}/{}.{}.rf'.format(direc, second_path, str(FUV_filtnum))
    U_file = '{}/{}.{}.rf'.format(direc, second_path, str(U_filtnum))
    V_file = '{}/{}.{}.rf'.format(direc, second_path, str(V_filtnum))
    J_file = '{}/{}.{}.rf'.format(direc, second_path, str(J_filtnum))
    
    FUV_tab = Table.read(FUV_file, format='ascii')
    U_tab = Table.read(U_file, format='ascii')
    V_tab = Table.read(V_file, format='ascii')
    J_tab = Table.read(J_file, format='ascii')
    
    FUVU_join = join(FUV_tab, U_tab, keys=key, table_names=['FUV', 'U'],
                     uniq_col_name='{table_name}_{col_name}')
    FUVUV_join = join(FUVU_join, V_tab, keys=key, table_names=['U', 'V'],
                      uniq_col_name='{table_name}_{col_name}')
    FUVUVJ_join = join(FUVUV_join, J_tab, keys=key, table_names=['V', 'J'],
                       uniq_col_name='{table_name}_{col_name}')
    
    DM = FUVUVJ_join['FUV_DM']
    
    M_AB_FUV = -2.5*np.log10(FUVUVJ_join['L' + str(FUV_filtnum)]) + 25 - DM
    M_AB_U = -2.5*np.log10(FUVUVJ_join['L' + str(U_filtnum)]) + 25 - DM
    M_AB_V = -2.5*np.log10(FUVUVJ_join['L' + str(V_filtnum)]) + 25 - DM
    M_AB_J = -2.5*np.log10(FUVUVJ_join['L' + str(J_filtnum)]) + 25 - DM
    
    table = Table([FUVUVJ_join['id'], M_AB_FUV, M_AB_U, M_AB_V, M_AB_J],
                  names=('id', 'M_AB_FUV', 'M_AB_U', 'M_AB_V', 'M_AB_J') )
    if write :
        table.write('{}/{}_colors.fits'.format(cluster, cluster))
    
    if plot_hist :
        plt.plot_UVJ_hist(table['M_AB_V'] - table['M_AB_J'],
                          table['M_AB_U'] - table['M_AB_V'],
                          xlabel=r'$V - J$', ylabel=r'$U - V$',
                          xmin=0, xmax=1.9, ymin=0, ymax=2.4)
    
    return table

def create_blank_models() :
    '''
    Create blank "models" for the missing model images: F225W for A1063,
    F225W and F275W for MACS J1149.
    
    Returns
    -------
    None.
    
    '''
    
    a1063_dim = (5190, 4860)
    m1149_dim = (5400, 5400)
    dims = [a1063_dim, m1149_dim, m1149_dim]
    
    outfiles = ['abell1063_f225w_bcgs_model.fits.gz',
                'macs1149_f225w_bcgs_model.fits.gz',
                'macs1149_f275w_bcgs_model.fits.gz']
    
    for i in range(len(dims)) :
        data = np.zeros(dims[i])
        hdu = fits.PrimaryHDU(data)
        hdu.writeto('models/' + outfiles[i])
    
    return

def determine_initial_sample(cluster, redshift, sigma, key='id', redshift_type='z_best',
                             redshift_tol_lo=0.01, redshift_tol_hi=0.01,
                             others=False) :
    '''
    Determine the final objects to use for analysis, based on various quality
    flags included in the catalog tables. Also combine relevant data including
    FAST output with the photometric data.
    
    Parameters
    ----------
    cluster : string
        String describing the cluster which informs the path to relevant data.
    redshift : float
        The redshift of the cluster.
    key : string, optional
        The key to base the table join operation on. The default is 'id'.
    z_spec : bool, optional
        Flag to use spectroscopic redshifts. The default is True.
    plot : bool, optional
        Boolean to plot all objects on a scatter plot. The default is False.
    redshift_tol_lo : float, optional
        The low redshift tolerance to use. The default is 0.01.
    redshift_tol_hi : float, optional
        The high redshift tolerance to use. The default is 0.01.
    
    Returns
    -------
    combined[mask] : astropy.table.Table
        The intial good objects which will be used for subsequent analysis.
    
    '''
    
    # read in the data in the two files
    fastoutput = Table.read('catalogs/{}_fastoutput.fits'.format(cluster))
    photometry = Table.read('catalogs/{}_photometry.fits'.format(cluster))
    
    # join the data based on a specific key
    combined = join(fastoutput, photometry, keys=key)
    
    # use the best redshift available for the field (ie. foreground/background
    # galaxies for the cluster pointings) and cluster galaxies
    z_best = []
    for i in range(len(combined)) :
        if combined['z_spec'][i] > 0 :
            z_best.append(combined['z_spec'][i])
        else :
            z_best.append(combined['z'][i])
    combined['z_best'] = z_best
    
    # determine the environment for the galaxies
    delta_z = (3*sigma/const.c.to(u.km/u.s))*(1 + redshift)
    env = []
    for i in range(len(combined)) :
        if ((combined['z_best'][i] >= redshift - delta_z) &
            (combined['z_best'][i] <= redshift + delta_z)) :
            env.append('cluster')
        else :
            env.append('field')
    combined['env'] = env
    
    combined['cluster'] = [cluster]*len(combined)
    
    # mask the data based on a good redshift, then remove any
    # point sources and "uncertain" sources (as per Shipley+2018), then
    # remove not "OK" sources (as per Shipley+2018), then remove sources
    # that are not in the footprint of the F160W image, then remove not-massive
    # objects, and finally mask based on redshift
    mask = ((combined[redshift_type] > 0) &
            (combined['star_flag'] == 0) & (combined['use_phot'] == 1) &
            (combined['flag_F160W'] != -1) & (combined['flag_F160W'] != -99) &
            (combined['lmass'] >= 8) &
            (combined[redshift_type] >= redshift-redshift_tol_lo) &
            (combined[redshift_type] <= redshift+redshift_tol_hi))
    
    if others :
        determine_other_sample(combined, combined[mask], redshift, redshift_type,
                               redshift_tol_lo, redshift_tol_hi, cluster)
    
    return combined[mask]

def determine_final_sample(cluster, redshift, sigma, first_path, second_path,
                           key='id', FUV_filtnum=218, U_filtnum=153,
                           V_filtnum=155, J_filtnum=161,
                           redshift_tol_lo=0.01, redshift_tol_hi=0.01,
                           redshift_type='z_best', plot_all=False, plot_uvj=False,
                           write_final_objs=False, write_regions=False,
                           selection='FUVVJ', verbose=False) :
    '''
    
    
    Parameters
    ----------
    cluster : TYPE
        DESCRIPTION.
    redshift : TYPE
        DESCRIPTION.
    first_path : TYPE
        DESCRIPTION.
    second_path : TYPE
        DESCRIPTION.
    key : string, optional
        The key to base the table join operation on. The default is 'id'.
    FUV_filtnum : int, optional
        Number describing FUV filter color indices file. The default is 218.
    U_filtnum : int, optional
        Number describing U filter color indices file. The default is 153.
    V_filtnum : int, optional
        Number describing V filter color indices file. The default is 155.
    J_filtnum : int, optional
        Number describing J filter color indices file. The default is 161.
    redshift_type : str, optional
        Flag to use spectroscopic redshifts. The default is 'z_best'.
    redshift_tol_lo : float, optional
        The low redshift tolerance to use. The default is 0.01.
    redshift_tol_hi : float, optional
        The high redshift tolerance to use. The default is 0.01.
    plot_all : TYPE, optional
        DESCRIPTION. The default is False.
    plot_uvj : TYPE, optional
        DESCRIPTION. The default is False.
    write_final_objs : TYPE, optional
        DESCRIPTION. The default is False.
    write_regions : TYPE, optional
        DESCRIPTION. The default is False.
    selection : string, optional
        Determines how QGs and SFGs are determined. The default is 'FUVVJ'.
    
    Returns
    -------
    None.
    
    '''
    
    initial = determine_initial_sample(cluster, redshift, sigma, key=key,
                                       redshift_type=redshift_type,
                                       redshift_tol_lo=redshift_tol_lo,
                                       redshift_tol_hi=redshift_tol_hi)
    colors = combine_colors(first_path, second_path, FUV_filtnum, U_filtnum,
                            V_filtnum, J_filtnum, cluster, key=key)
    if len(initial) > 0 :
        final = join(initial, colors, keys=key)
        
        # determine colors and set values for the different color-color diagrams
        if selection == 'UVJ' : # using values from Muzzin et al. 2013
            xs = final['M_AB_V'] - final['M_AB_J']
            ys = final['M_AB_U'] - final['M_AB_V']
            slope, intercept, horiz, vert = 0.88, 0.69, 1.3, 1.5
        else : # using values from Leja et al. 2019b
            xs = final['M_AB_V'] - final['M_AB_J']
            ys = final['M_AB_FUV'] - final['M_AB_V']
            slope, intercept, horiz, vert = 3.24, 0.32, 3.45, 1.56
        
        first_knee, second_knee = (horiz - intercept)/slope, vert
        
        # define the quiescent region
        corner = ( ((xs <= first_knee) & (ys >= horiz)) |
                   ( ((xs >= first_knee) & (xs <= second_knee))
                     & (ys >= slope*xs + intercept)) )
        
        # create a column describing the population type
        pop_type = np.empty(len(final), dtype=str)
        pop_type[corner], pop_type[~corner] = 'Q', 'S'
        pop_type  = np.char.replace(pop_type, 'S', 'SF')
        final['pop'] = pop_type
        
        # QG_mask, SF_mask = (final['pop'] == 'Q'), (final['pop'] == 'SF')
        # QGs, SFGs = final.copy(), final.copy()
        # QGs, SFGs = QGs[QG_mask], SFGs[SF_mask]
        # print('{}: {} QGs, {} SFGs'.format(cluster, len(QGs), len(SFGs)))
        
        if plot_uvj :
            plt.plot_UVJ(final['M_AB_V'] - final['M_AB_J'],
                         final['M_AB_U'] - final['M_AB_V'],
                         xlabel=r'$V - J$', ylabel=r'$U - V$', title=cluster,
                         xmin=0, xmax=1.9, ymin=0, ymax=2.4)
        
        if write_final_objs :
            final_objs_file = '{}/{}_sample.fits'.format(cluster, cluster)
            final.write(final_objs_file)
            if verbose :
                print('Saved: Final objects to fits file.')
        
        if write_regions :
            save_regions(cluster, final)
            if verbose :
                print('Saved: Regions to file.')
    
    return

def determine_other_sample(combined, good_sample, redshift, redshift_type,
                           redshift_tol_lo, redshift_tol_hi, cluster) :
    # determine other objects that didn't satisfy the initial sample criteria
    
    # galaxies with stellar mass < 10^8 solMass, but within the redshift limits
    small = ((combined['lmass'] < 8) &
             (combined['{}'.format(redshift_type)] >= redshift-redshift_tol_lo) &
             (combined['{}'.format(redshift_type)] <= redshift+redshift_tol_hi))
    small = combined[small]
    
    # other objects beyond the redshift limits
    far = ((combined['{}'.format(redshift_type)] < redshift - redshift_tol_lo) |
           (combined['{}'.format(redshift_type)] > redshift + redshift_tol_hi))
    far = combined[far]
    
    # create the scatter plot
    xs = [far['lmass'], small['lmass'], good_sample['lmass']]
    ys = [far['{}'.format(redshift_type)], small['{}'.format(redshift_type)],
          good_sample['{}'.format(redshift_type)]]
    labels = [r'$z~\notin~z_c \pm \delta z$',
              r'$M_* < 10^{8}~M_{\odot}$',
              'final ({})'.format(len(good_sample))]
    styles = ['x', '+', 'o']
    colors = ['k', 'k', 'k']
    sizes = [75, 80, 40]
    alphas = [1, 1, 0.4]
    plt.plot_objects(xs, ys, redshift, labels, styles, colors, sizes, alphas,
                     redshift_tol_lo=redshift_tol_lo,
                     redshift_tol_hi=redshift_tol_hi,
                     xlabel=r'$\log(M_{*}/M_{\odot})$',
                     ylabel=r'$z$', title=cluster,
                     xmin=2, xmax=14, ymin=0.2, ymax=0.8)
    
    return

def exp(xx, aa, bb, cc) :
    return aa - np.exp(-(xx-bb)/cc)

def linear(xx, aa, bb) :
    return aa*xx + bb

def open_cutout(infile, exp=False, simple=False, phot=False) :
    '''
    Open a cutout image. If opening the F160W image in order to compute the
    Voroni bins, then also return the total exposure time for that image,
    to assist with the noise calculation.
    
    Parameters
    ----------
    infile : TYPE
        DESCRIPTION.
    
    Returns
    -------
    data : TYPE
        DESCRIPTION.
    dim : TYPE
        DESCRIPTION.
    photfnu : TYPE
        DESCRIPTION.
    
    '''
    with fits.open(infile) as hdu :
        hdr = hdu[0].header
        redshift = hdr['Z']
        exptime = hdr['EXPTIME']
        photfnu = hdr['PHOTFNU']
        R_e = hdr['R_e']
        sma = hdr['SMA']
        smb = hdr['SMB']
        pa = hdr['PA']
        data = hdu[0].data
        dim = data.shape
    
    if simple :
        return data
    elif exp :
        return data, exptime
    elif phot :
        return data*photfnu
    else :
        return data, dim, photfnu, R_e, redshift, sma, smb, pa

def open_image(cluster, ID, direc) :
    
    with fits.open('{}/{}_images/{}_ID_{}_{}.fits'.format(
            cluster, direc, cluster, ID, direc)) as hdu :
        image = hdu[0].data
    
    return image

def save_bin_images() :
    
    HFF = Table.read('output/tables/nbCGs.fits')
    
    for cluster, ID in zip(HFF['cluster'], HFF['ID']) :
        os.makedirs('{}/bins_images'.format(cluster), exist_ok=True)
        # ensure the output directory for the bins images is available
        
        # save images based on the annuli/binning
        bin_data = np.load('{}/bins/{}_ID_{}_annuli.npz'.format(
            cluster, cluster, ID))
        
        hdu = fits.PrimaryHDU(bin_data['image'])
        hdu.writeto('{}/bins_images/{}_ID_{}_bins.fits'.format(
            cluster, cluster, ID))
    
    return

def save_cutout(xx, yy, angular_size, data, outfile, exposure, photfnu, scale,
                rms, r_e, redshift, sma, smb, pa, vmin=None, vmax=None,
                show=False) :
    '''
    Save an individual cutout image based on RA and Dec, and given the WCS
    information. Include basic information in the header such as total
    exposure time in that filter, the inverse sensitivity (PHOTFNU), and the
    pixel scale for that filter.
    
    Parameters
    ----------
    sky_ra : TYPE
        DESCRIPTION.
    sky_dec : TYPE
        DESCRIPTION.
    angular_size : TYPE
        DESCRIPTION.
    data : TYPE
        DESCRIPTION.
    outfile : TYPE
        DESCRIPTION.
    exposure : TYPE
        DESCRIPTION.
    photfnu : TYPE
        DESCRIPTION.
    scale : TYPE
        DESCRIPTION.
    rms : TYPE
        DESCRIPTION.
    r_e : TYPE
        DESCRIPTION.
    redshift : TYPE
        DESCRIPTION.
    sma : TYPE
        DESCRIPTION.
    smb : TYPE
        DESCRIPTION.
    pa : TYPE
        DESCRIPTION.
    vmin : TYPE, optional
        DESCRIPTION. The default is None.
    vmax : TYPE, optional
        DESCRIPTION. The default is None.
    show : TYPE, optional
        DESCRIPTION. The default is False.
    
    Returns
    -------
    None.
    
    '''
    
    position = (xx, yy) # the position of the cutout's center
    size = np.ceil(angular_size.value) # size along half an axis, in pixels
    cutout = Cutout2D(data, position, 2*size) # cutout has radius of 'size'
    
    hdu = fits.PrimaryHDU(cutout.data)
    hdr = hdu.header
    hdr['Z'] = redshift
    hdr.comments['Z'] = 'object spectroscopic redshift--from table'
    hdr['EXPTIME'] = exposure
    hdr.comments['EXPTIME'] = 'exposure duration (seconds)--calculated'
    hdr['PHOTFNU'] = photfnu
    hdr.comments['PHOTFNU'] = 'inverse sensitivity, Jy*sec/electron'
    hdr['SCALE'] = scale
    hdr.comments['SCALE'] = 'Drizzle, pixel size (arcsec) of output image'
    hdr['RMS'] = rms
    hdr.comments['RMS'] = 'RMS value of science image--calculated'
    hdr['R_E'] = r_e
    hdr.comments['R_E'] = 'Radius enclosing half the total flux, pixels'
    hdr['SMA'] = sma
    hdr.comments['SMA'] = 'Semi-major axis, pixels'
    hdr['SMB'] = smb
    hdr.comments['SMB'] = 'Semi-minor axis, pixels'
    hdr['PA'] = pa
    hdr.comments['PA'] = 'Posn. angle of SMA, count.-clock. from East'
    hdu.writeto(outfile)
    
    if show :
        plt.display_image_simple(cutout.data, vmin=vmin, vmax=vmax)
    
    return

def save_flux_fraction_images() :
    
    HFF = Table.read('output/tables/nbCGs.fits')
    
    for cluster, ID in zip(HFF['cluster'], HFF['ID']) :
        # ensure the output directory for the flux fraction images is available
        os.makedirs('{}/fluxfrac_images'.format(cluster), exist_ok=True)
        
        # save flux fraction images based on F160W
        flux = open_cutout('{}/cutouts/{}_ID_{}_f160w.fits'.format(
            cluster, cluster, ID), simple=True)
        
        segmap = open_cutout('{}/cutouts/{}_ID_{}_segmap.fits'.format(
            cluster, cluster, ID), simple=True)
        
        bins = open_image(cluster, ID, 'bins')
        
        # mask pixels associated with other galaxies, that aren't in an annulus
        flux[(segmap > 0) & (segmap != ID)] = 0
        flux[np.isnan(bins)] = 0
        
        hdu = fits.PrimaryHDU(flux/np.sum(flux))
        hdu.writeto('{}/fluxfrac_images/{}_ID_{}_fluxfrac.fits'.format(
            cluster, cluster, ID))
    
    return

def save_regions(cluster, table) :
    '''
    Save region files containing the locations of the different populations.
    
    Parameters
    ----------
    cluster : string
        String describing the cluster which informs the path to relevant data.
    table : astropy.table.Table
        Table of final objects which will be used for subsequent analysis.
    
    Returns
    -------
    None.
    
    '''
    
    os.makedirs('{}/regions'.format(cluster), exist_ok=True) # ensure the
        # output directory is available
    
    q_region_file = '{}/regions/{}_final_QGs_R_e.reg'.format(cluster, cluster)
    sf_region_file = '{}/regions/{}_final_SFGs_R_e.reg'.format(cluster, cluster)
    bCG_region_file = '{}/regions/{}_final_bCGs_R_e.reg'.format(cluster, cluster)
    bCG_sqr_file = '{}/regions/{}_bCGs_sqr.reg'.format(cluster, cluster)
    bCG_loc_file = '{}/regions/{}_bCGs_loc.reg'.format(cluster, cluster)
    
    first = '# Region file format: DS9 version 4.1\n'
    q_second = ('global color=red width=2 select=1 ' +
                'edit=1 move=1 delete=1 include=1\n')
    sf_second = ('global color=blue width=2 select=1 ' +
                 'edit=1 move=1 delete=1 include=1\n')
    bCG_second = ('global color=green width=2 select=1 ' +
                  'edit=1 move=1 delete=1 include=1\n')
    third = 'fk5\n'
    
    with open(q_region_file, 'a+') as file :
        file.write(first)
        file.write(q_second)
        file.write(third)
        for i in range(len(table)) :
            if (table['pop'][i] == 'Q') & (table['id'][i] < 20000) :
                string = 'circle({},{},{}") # {}\n'.format(
                    str(table['ra'][i]), str(table['dec'][i]),
                    str(table['flux_radius'][i]*0.06), str(table['id'][i]))
                file.write(string)
    
    with open(sf_region_file, 'a+') as file :
        file.write(first)
        file.write(sf_second)
        file.write(third)
        for i in range(len(table)) :
            if (table['pop'][i] == 'SF') & (table['id'][i] < 20000) :
                string = 'circle({},{},{}") # {}\n'.format(
                    str(table['ra'][i]), str(table['dec'][i]),
                    str(table['flux_radius'][i]*0.06), str(table['id'][i]))
                file.write(string)
    
    with open(bCG_region_file, 'a+') as file :
        file.write(first)
        file.write(bCG_second)
        file.write(third)
        for i in range(len(table)) :
            if (table['id'][i] > 20000) :
                string = 'circle({},{},{}") # {}\n'.format(
                    str(table['ra'][i]), str(table['dec'][i]),
                    str(table['flux_radius'][i]*0.06), str(table['id'][i]))
                file.write(string)
    
    with open(bCG_sqr_file, 'a+') as file :
        file.write(first)
        file.write(q_second)
        file.write('physical\n')
        for i in range(len(table)) :
            if (table['id'][i] > 20000) :
                string = 'box({},{},{},{},0) # {}\n'.format(
                    str(table['x'][i]), str(table['y'][i]),
                    str(table['flux_radius'][i]*10),
                    str(table['flux_radius'][i]*10), str(table['id'][i]))
                file.write(string)
    
    with open(bCG_loc_file, 'a+') as file :
        file.write(first)
        file.write(bCG_second)
        file.write(third)
        for i in range(len(table)) :
            if (table['id'][i] > 20000) :
                string = 'text({},{}) # text={{{}}}\n'.format(
                    str(table['ra'][i]), str(table['dec'][i]), str(table['id'][i]))
                file.write(string)
    
    return

def tanh(xx, aa, bb, cc, dd) :
    return aa + bb*np.tanh((xx - cc)/dd)
