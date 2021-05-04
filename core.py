
import os
import numpy as np

from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.table import hstack, join, Table
from scipy import interpolate, stats

import plotting as plt

def build_par_file_WFC3(filt, wave, throughput) :
    
    # sedpy needs evenly spaced wavelength arrays and corresponding throughputs
    interpolated = interpolate.interp1d(wave, throughput)
    start, end = wave[0], wave[-1]
    
    # evenly sample the wavelengths and use interpolated throughput values
    wave_new = np.arange(start, end+1, 1)
    throughput_new = interpolated(wave_new)
    
    # write the corresponding tables
    par_table = Table([wave_new, throughput_new], names=('wave', 'thru'))
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
        
        # pad the wavelength arrays and throughputs with leading/trailing zeros
        wave_front = np.arange(wave[0]-100, wave[0], 1)
        wave_end = np.arange(wave[-1], wave[-1]+100, 1)
        zeros_front = np.zeros(wave_front.shape)
        zeros_end = np.zeros(wave_end.shape)
        
        # create final arrays
        wave = np.concatenate([wave_front, wave, wave_end])
        throughput = np.concatenate([zeros_front, throughput, zeros_end])
        
        # interpolate values and write files
        build_par_file_WFC3(filt, np.float64(wave), throughput)
    
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
    
    filters = os.listdir('hst_filters/')
    
    # HST WFC3/UVIS (Wide Field Camera 3 - UV/Visible channel)
    # F225W, F275W, F336W, F390W
    # (UV wide, UV wide, U/Stromgren u, Washington C)
    UVIS = filters[6:]
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
    ACS = filters[0]
    ACS_table = Table.read('hst_filters/{}'.format(ACS), format='csv')
    if save_par :
        build_par_file_ACS(ACS_table)
    
    # HST WFC3/IR (Wide Field Camera 3 - near-IR channel)
    # F105W, F110W, F125W, F140W, F160W
    # (wide Y, wide YJ, wide J, wide JH gap, H)
    IR = filters[1:6]
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
                   V_filtnum=155, J_filtnum=161, key='id', plot_hist=False) :
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

def determine_final_flags(table) :
    '''
    Determine the final flag for each object in `table`.
    
    Parameters
    ----------
    table : astropy.table.Table
        The table to determine the final flags for.
    
    Returns
    -------
    final_flag : list
        The final flags for the objects in the table.
    
    '''
    
    # determine the band flags that are in the catalog
    band_flags = [string for string in table.colnames if 'flag_F' in string]
    
    # use only those columns to create a new subtable
    flag_table = Table([table[band_flag] for band_flag in band_flags],
                       names=tuple(band_flags))
    
    # determine the appropriate final flag for each galaxy
    final_flag = []
    for row in flag_table.iterrows() :
        '''
        row = np.array(row)
        # mask out any zeros
        non_zero = row[row != 0]
        mode, count = stats.mode(non_zero)
        if len(mode) == 0 :
            final_flag.append(0)
        else :
            # print(mode[0])
            final_flag.append(mode[0])
        '''
        final_flag.append(0)
    
    # print(final_flag)
    # table['final_flag'] = final_flag
    # set_of_sums = set(combined['final_flag'])
    # print(np.sort(list(set_of_sums)))
    
    return final_flag

def determine_final_objects(cluster, redshift, key='id',
                            z_spec=True, plot=False,
                            redshift_tol_lo=0.01, redshift_tol_hi=0.01) :
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
    final_objs : astropy.table.Table
        The final objects which will be used for subsequent analysis.
    
    '''
    
    # read in the data in the two files
    fastoutput = Table.read('catalogs/{}_fastoutput.fits'.format(cluster))
    photometry = Table.read('catalogs/{}_photometry.fits'.format(cluster))
    
    # join the data based on a specific key
    combined = join(fastoutput, photometry, keys=key)
    
    # mask the data based on a good spectroscopic redshift
    if z_spec :
        mask = combined['z_spec'] > 0
        combined = combined[mask]
    else :
        mask = combined['z'] > 0
        combined = combined[mask]
    
    # remove any point sources and "uncertain" sources, as per Shipley+ (2018)
    mask = combined['star_flag'] == 0
    combined = combined[mask]
    
    # remove not "OK" sources, as per Shipley+ (2018)
    mask = combined['use_phot'] == 1
    combined = combined[mask]
    
    # remove sources that are not in the footprint of the F160W image
    mask = (combined['flag_F160W'] != -1) & (combined['flag_F160W'] != -99)
    combined = combined[mask]
    
    # mask the modelled-and-subtracted bCGs
    # mask = combined['bandtotal'] != 'bcg'
    # combined = combined[mask]
    
    # mask specific galaxies from certain clusters
    if cluster == 'm717' :
        mask = ((combined['id'] != 777) & (combined['id'] != 5551) &
                (combined['id'] != 5890) & (combined['id'] != 5956) &
                (combined['id'] != 6308) & (combined['id'] != 6335))
        # the above IDs all have no synethsized FUV, U, V, J color indices
    if cluster == 'a1063par' :
        mask = (combined['id'] != 789) & (combined['id'] != 4690)
        # ID 789 is too faint and resulted in no annuli
        # ID 4690 is an errant detection on the edge of the field
        combined = combined[mask]
    
    # determine the final flags for all galaxies
    final_flags = determine_final_flags(combined)
    combined['final_flag'] = final_flags
    
    # plot the final objects that will be used
    if z_spec :
        final_objs_mask = ((combined['final_flag']==0) &
                           (combined['lmass'] >= 8) &
                           (combined['z_spec'] >= redshift - redshift_tol_lo) &
                           (combined['z_spec'] <= redshift + redshift_tol_hi))
    else :
        final_objs_mask = ((combined['final_flag']==0) &
                           (combined['lmass'] >= 8) &
                           (combined['z'] >= redshift - redshift_tol_lo) &
                           (combined['z'] <= redshift + redshift_tol_hi))
    final_objs = combined[final_objs_mask]
    final_objs_lmass = final_objs['lmass']
    if z_spec :
        final_objs_z = final_objs['z_spec']
    else :
        final_objs_z = final_objs['z']
    length_of_final = len(final_objs)
    
    # plot other objects with no flags, but with a stellar mass < 10^8 solMass
    if z_spec :
        small_mass_mask = ((combined['final_flag']==0) &
                           (combined['lmass'] < 8) &
                           (combined['z_spec'] >= redshift - redshift_tol_lo) &
                           (combined['z_spec'] <= redshift + redshift_tol_hi))
    else :
        small_mass_mask = ((combined['final_flag']==0) &
                           (combined['lmass'] < 8) &
                           (combined['z'] >= redshift - redshift_tol_lo) &
                           (combined['z'] <= redshift + redshift_tol_hi))
    small_mass_objs = combined[small_mass_mask]
    small_mass_objs_lmass = small_mass_objs['lmass']
    if z_spec :
        small_mass_objs_z = small_mass_objs['z_spec']
    else :
        small_mass_objs_z = small_mass_objs['z']
    
    # plot other objects with no flags, but beyond the redshift limits
    if z_spec :
        no_flags_mask = ((combined['final_flag']==0) &
                         ((combined['z_spec'] < redshift - redshift_tol_lo) |
                          (combined['z_spec'] > redshift + redshift_tol_hi)))
    else :
        no_flags_mask = ((combined['final_flag']==0) &
                         ((combined['z'] < redshift - redshift_tol_lo) |
                          (combined['z'] > redshift + redshift_tol_hi)))
    no_flags_objs = combined[no_flags_mask]
    no_flags_objs_lmass = no_flags_objs['lmass']
    if z_spec :
        no_flags_objs_z = no_flags_objs['z_spec']
    else :
        no_flags_objs_z = no_flags_objs['z']
    
    # create the scatter plot
    if plot :
        xs = [no_flags_objs_lmass, small_mass_objs_lmass, final_objs_lmass]
        ys = [no_flags_objs_z, small_mass_objs_z, final_objs_z]
        labels = [r'$z~\notin~z_c \pm \delta z$',
                  r'$M_* < 10^{8}~M_{\odot}$',
                  'final ({})'.format(length_of_final)]
        styles = ['x', '+', 'o']
        colors = ['k', 'k', 'k']
        sizes = [75, 80, 40]
        alphas = [1, 1, 0.4]
        plt.plot_objects(xs, ys, redshift, labels, styles, colors,sizes,alphas,
                         redshift_tol_lo=redshift_tol_lo,
                         redshift_tol_hi=redshift_tol_hi,
                         xlabel=r'$\log(M_{*}/M_{\odot})$',
                         ylabel=r'$z$', title=cluster,
                         xmin=2, xmax=14, ymin=0.2, ymax=0.8)
    
    return final_objs

def determine_finalObjs_w_color(cluster, redshift, first_path, second_path,
                                key='id', FUV_filtnum=218, U_filtnum=153,
                                V_filtnum=155, J_filtnum=161,
                                redshift_tol_lo=0.01, redshift_tol_hi=0.01,
                                z_spec=True, plot_all=False, plot_uvj=False,
                                write_final_objs=False, write_regions=False,
                                selection='FUVVJ') :
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
    z_spec : bool, optional
        Flag to use spectroscopic redshifts. The default is True.
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
    
    final_objs = determine_final_objects(cluster, redshift, key=key,
                                         z_spec=z_spec, plot=plot_all,
                                         redshift_tol_lo=redshift_tol_lo,
                                         redshift_tol_hi=redshift_tol_hi)
    colors = combine_colors(first_path, second_path, FUV_filtnum, U_filtnum,
                            V_filtnum, J_filtnum, key=key)
    final = join(final_objs, colors, keys=key)
    
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
    
    if plot_uvj :
        plt.plot_UVJ(final['M_AB_V'] - final['M_AB_J'],
                     final['M_AB_U'] - final['M_AB_V'],
                     xlabel=r'$V - J$', ylabel=r'$U - V$', title=cluster,
                     xmin=0, xmax=1.9, ymin=0, ymax=2.4)
    
    if write_final_objs :
        final_objs_file = '{}/{}_final_objects.fits'.format(cluster, cluster)
        final.write(final_objs_file)
    
    if write_regions :
        save_regions(cluster, final)
    
    return

def open_cutout(infile) :
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
        photfnu = hdr['PHOTFNU']
        R_e = hdr['R_e']
        sma = hdr['SMA']
        smb = hdr['SMB']
        pa = hdr['PA']
        data = hdu[0].data
        dim = data.shape
    
    return data, dim, photfnu, R_e, redshift, sma, smb, pa

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
    
    q_region_file = '{}/{}_final_QGs_R_e.reg'.format(cluster, cluster)
    sf_region_file = '{}/{}_final_SFGs_R_e.reg'.format(cluster, cluster)
    bCG_region_file = '{}/{}_final_bCGs_R_e.reg'.format(cluster, cluster)
    
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
            if (table['pop'][i] == 'Q') & (table['bandtotal'][i] != 'bcg') :
                string = 'circle({},{},{}") # {}\n'.format(
                    str(table['ra'][i]), str(table['dec'][i]),
                    str(table['flux_radius'][i]*0.06), str(table['id'][i]))
                file.write(string)
    
    with open(sf_region_file, 'a+') as file :
        file.write(first)
        file.write(sf_second)
        file.write(third)
        for i in range(len(table)) :
            if (table['pop'][i] == 'SF') & (table['bandtotal'][i] != 'bcg') :
                string = 'circle({},{},{}") # {}\n'.format(
                    str(table['ra'][i]), str(table['dec'][i]),
                    str(table['flux_radius'][i]*0.06), str(table['id'][i]))
                file.write(string)
    
    with open(bCG_region_file, 'a+') as file :
        file.write(first)
        file.write(bCG_second)
        file.write(third)
        for i in range(len(table)) :
            if (table['pop'][i] == 'Q') & (table['bandtotal'][i] == 'bcg') :
                string = 'circle({},{},{}") # {}\n'.format(
                    str(table['ra'][i]), str(table['dec'][i]),
                    str(table['flux_radius'][i]*0.06), str(table['id'][i]))
                file.write(string)
    
    return
