
import os
import numpy as np

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.table import hstack, join, Table
import astropy.units as u
from scipy import stats

import plotting as plt

def build_par_file_WFC3(filt, wave, throughput) :
    
    kfilter = len(wave)*['KFILTER']
    par_table = Table([kfilter, wave, throughput],
                      names=('kfilt', 'wave', 'thru'))
    par_table.write('hff_filters/hff_' + filt + '.par',
                    format='ascii.no_header', overwrite=False)
    
    return

def build_par_file_ACS(table) :
    
    cols = table.colnames
    for i in range(0, len(cols), 2) :
        filt = cols[i].split('_')[0]
        kfilter = len(table[cols[i]])*['KFILTER']
        par_table = Table([kfilter, table[cols[i]], table[cols[i+1]]],
                          names=('kfilt', 'wave', 'thru'))
        par_table.write('hff_filters/hff_' + filt + '.par',
                        format='ascii.no_header', overwrite=False)
    
    return

def build_filter_table(display=False, write=False) :
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
        table = Table.read('hst_filters/' + file, format='csv')
        filt = file[:5]
        wave1 = table['Chip 1 Wave (Angstroms)']
        wave2 = table['Chip 2 Wave (Angstroms)']
        wave = (wave1 + wave2)/2
        thru1 = table['Chip 1 Throughput']
        thru2 = table['Chip 2 Throughput']
        thru = (thru1 + thru2)/2
        # build_par_file_WFC3(filt, wave, thru)
        filt_table = Table([wave, thru], names=(filt+'_wave', filt+'_thru'))
        UVIS_tables.append(filt_table)
    UVIS_table = hstack(UVIS_tables)
    
    # HST ACS/WFC (Advanced Camera for Surveys - Wide Field Channel)
    # F435W, F475W, F555W, F606W, F625W, F775W, F814W, F850LP
    # (Johnson B, SDSS g, Johnson V, Broad V, SDSS r, SDSS i, Broad I, SDSS z)
    ACS = filters[0]
    ACS_table = Table.read('hst_filters/' + ACS, format='csv')
    # build_par_file_ACS(ACS_table)
    
    # HST WFC3/IR (Wide Field Camera 3 - near-IR channel)
    # F105W, F110W, F125W, F140W, F160W
    # (wide Y, wide YJ, wide J, wide JH gap, H)
    IR = filters[1:6]
    IR_tables = []
    for file in IR :
        table = Table.read('hst_filters/' + file, format='csv')
        filt = file[:5]
        wave = table['Wave (Angstroms)']
        thru = table['Throughput']
        # build_par_file_WFC3(filt, wave, thru)
        filt_table = Table([wave, thru], names=(filt+'_wave', filt+'_thru'))
        IR_tables.append(filt_table)
    IR_table = hstack(IR_tables)
    
    master_table = hstack([UVIS_table, ACS_table, IR_table])
    
    if display :
        print(master_table)
    
    if write :
        master_table.write('filter_throughputs.fits', overwrite=True)
    
    return master_table

def combine_UVJ(first_path, second_path, U_filtnum, V_filtnum, J_filtnum,
                plotuvj_hist=False) :
    '''
    
    
    Parameters
    ----------
    first_path : string
        DESCRIPTION.
    second_path : string
        DESCRIPTION.
    U_filtnum : int
        DESCRIPTION.
    V_filtnum : int
        DESCRIPTION.
    J_filtnum : int
        DESCRIPTION.
    plotuvj : bool, optional
        Boolean to plot the values on a scatter plot. The default is False.
    
    Returns
    -------
    table : astropy.table.Table
        DESCRIPTION.
    
    '''
    
    direc = 'catalogs/' + first_path + '/eazy/' + second_path + '/'
    U_file = direc + second_path + ('.%s' % U_filtnum) + '.rf'
    V_file = direc + second_path + ('.%s' % V_filtnum) + '.rf'
    J_file = direc + second_path + ('.%s' % J_filtnum) + '.rf'
    
    U_tab = Table.read(U_file, format='ascii')
    V_tab = Table.read(V_file, format='ascii')
    J_tab = Table.read(J_file, format='ascii')
    
    UV_join = join(U_tab, V_tab, keys='id')
    UVJ_join = join(UV_join, J_tab, keys='id')
    
    DM = UVJ_join['DM']
    
    M_AB_U = -2.5*np.log10(UVJ_join['L' + str(U_filtnum)]) + 25 - DM
    M_AB_V = -2.5*np.log10(UVJ_join['L' + str(V_filtnum)]) + 25 - DM
    M_AB_J = -2.5*np.log10(UVJ_join['L' + str(J_filtnum)]) + 25 - DM
    
    table = Table([UVJ_join['id'], M_AB_U, M_AB_V, M_AB_J],
                  names=('id', 'M_AB_U', 'M_AB_V', 'M_AB_J') )
    
    if plotuvj_hist :
        plt.plot_UVJ_hist(table['M_AB_V'] - table['M_AB_J'],
                          table['M_AB_U'] - table['M_AB_V'],
                          xlabel=r'$V - J$', ylabel=r'$U - V$',
                          xmin=0, xmax=1.9, ymin=0, ymax=2.4)
    
    return table

def determine_final_objects(cluster, key, redshift, plot=False) :
    '''
    Determine the final objects to use for analysis, based on various quality
    flags included in the catalog tables. Also combine relevant data including
    FAST output with the photometric data.
    
    Parameters
    ----------
    cluster : string
        String describing the cluster which informs the path to relevant data.
    key : string
        String describing the key which is used to join the data.
    redshift : float
        The redshift of the cluster.
    plot : bool, optional
        Boolean to plot all objects on a scatter plot. The default is False.
    
    Returns
    -------
    final_objs : astropy.table.Table
        The final objects which will be used for subsequent analysis.
    
    '''
    
    # read in the data in the two files
    fastoutput = Table.read('catalogs/' + cluster + '_cluster_fastoutput.fits')
    photometry = Table.read('catalogs/' + cluster + '_cluster_photometry.fits')
    
    # join the data based on a specific key
    combined = join(fastoutput, photometry, keys=key)
    
    # mask the data based on a good spectroscopic redshift
    mask = combined['z_spec'] > 0
    combined = combined[mask]
    
    # determine the band flags that are in the catalog
    band_flags = [string for string in combined.colnames if 'flag_F' in string]
    
    # use only those colmns to create a new subtable
    flag_table = Table([combined[band_flag] for band_flag in band_flags],
                       names=tuple(band_flags))
    
    # determine the appropriate final flag for each galaxy
    final_flag = []
    for row in flag_table.iterrows() :
        row = np.array(row)
        # mask out any zeros
        non_zero = row[row != 0]
        mode, count = stats.mode(non_zero)
        if len(mode) == 0 :
            final_flag.append(0)
        else :
            # print(mode[0])
            final_flag.append(mode[0])
    # print(final_flag)
    combined['final_flag'] = final_flag
    # set_of_sums = set(combined['final_flag'])
    # print(np.sort(list(set_of_sums)))
    
    # plot the final objects that will be used
    final_objs_mask = ((combined['final_flag']==0) &
                       (combined['lmass'] >= 8) &
                       (combined['z_spec'] >= redshift-0.05) &
                       (combined['z_spec'] <= redshift+0.05))
    final_objs = combined[final_objs_mask]
    final_objs_lmass = final_objs['lmass']
    final_objs_z = final_objs['z_spec']
    length_of_final = len(final_objs)
    
    # plot other objects with no flags, but with a stellar mass < 10^8 solMass
    small_mass_mask = ((combined['final_flag']==0) &
                       (combined['lmass'] < 8) &
                       (combined['z_spec'] >= redshift-0.05) &
                       (combined['z_spec'] <= redshift+0.05))
    small_mass_objs = combined[small_mass_mask]
    small_mass_objs_lmass = small_mass_objs['lmass']
    small_mass_objs_z = small_mass_objs['z_spec']
    
    # plot other objects with no flags, but beyond the redshift limits
    no_flags_mask = ((combined['final_flag']==0) &
                     ((combined['z_spec'] < redshift-0.05) |
                      (combined['z_spec'] > redshift+0.05)))
    no_flags_objs = combined[no_flags_mask]
    no_flags_objs_lmass = no_flags_objs['lmass']
    no_flags_objs_z = no_flags_objs['z_spec']
    
    # plot objects which are likely BCGs
    likely_BCGs_mask = combined['final_flag']==4
    likely_BCGs = combined[likely_BCGs_mask]
    likely_BCGs_lmass = likely_BCGs['lmass']
    likely_BCGs_z = likely_BCGs['z_spec']
    
    # plot all other types of objects: final_flag in [-99, -1, 2, 3]
    others_mask = ((combined['final_flag']==-99) |
                   (combined['final_flag']==-1) |
                   (combined['final_flag']==1) |
                   (combined['final_flag']==2) |
                   (combined['final_flag']==3))
    others = combined[others_mask]
    others_lmass = others['lmass']
    others_z = others['z_spec']
    
    if plot :
        # create the scatter plot
        xs = [others_lmass, likely_BCGs_lmass, no_flags_objs_lmass,
              small_mass_objs_lmass, final_objs_lmass]
        ys = [others_z, likely_BCGs_z, no_flags_objs_z,
              small_mass_objs_z, final_objs_z]
        labels = ['flags', 'BCGs', r'$z~\notin~z_c \pm 0.05$',
                  r'$M_* < 10^{8}~M_{\odot}$', 'final']
        styles = ['v', '^', 's', 's', 'o']
        colors = ['orange', 'm', 'brown', 'green', 'cyan']
        plt.plot_objects(xs, ys, redshift, length_of_final, labels, styles,
                         colors, xlabel=r'$\log(M_{*}/M_{\odot})$',
                         ylabel=r'$z$', title=cluster,
                         xmin=2, xmax=14, ymin=0, ymax=0.9)
    
    return final_objs

def determine_finalObjs_w_UVJ(cluster, key, redshift, first_path, second_path,
                              U_filtnum, V_filtnum, J_filtnum, plotuvj=False,
                              write_final_objs=False, write_regions=False) :
    '''
    
    
    Parameters
    ----------
    cluster : TYPE
        DESCRIPTION.
    key : TYPE
        DESCRIPTION.
    redshift : TYPE
        DESCRIPTION.
    first_path : TYPE
        DESCRIPTION.
    second_path : TYPE
        DESCRIPTION.
    U_filtnum : TYPE
        DESCRIPTION.
    V_filtnum : TYPE
        DESCRIPTION.
    J_filtnum : TYPE
        DESCRIPTION.
    plotuvj : TYPE, optional
        DESCRIPTION. The default is False.
    write_final_objs : TYPE, optional
        DESCRIPTION. The default is False.
    write_regions : TYPE, optional
        DESCRIPTION. The default is False.
    
    Returns
    -------
    None.
    
    '''
    
    final_objs = determine_final_objects(cluster, key, redshift)
    UVJ = combine_UVJ(first_path, second_path, U_filtnum, V_filtnum, J_filtnum,
                      plotuvj_hist=False)
    final = join(final_objs, UVJ, keys='id')
    
    if write_final_objs :
        final.write(cluster + '/' + cluster + '_final_objects.fits')
    
    if write_regions :
        region_file = (cluster + '/' + cluster + 
                       '_final_objects_flux_radius.reg')
        first = '# Region file format: DS9 version 4.1\n'
        second = ('global color=red width=2 select=1 ' +
                  'edit=1 move=1 delete=1 include=1\n')
        third = 'fk5\n'
        with open(region_file, 'a+') as file :
            file.write(first)
            file.write(second)
            file.write(third)
            for i in range(len(final)) :
                string = ('circle(' + str(final['ra'][i]) + ',' +
                          str(final['dec'][i]) + ',' +
                          str(final['flux_radius'][i]*0.06) + '") # ' +
                          str(final['id'][i]) + '\n' )
                file.write(string)
    
    if plotuvj :
        plt.plot_UVJ(final['M_AB_V'] - final['M_AB_J'],
                     final['M_AB_U'] - final['M_AB_V'],
                     xlabel=r'$V - J$', ylabel=r'$U - V$', title=cluster,
                     xmin=0, xmax=1.9, ymin=0, ymax=2.4)
    
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
        photfnu = hdr['PHOTFNU']
        R_e = hdr['R_e']
        data = hdu[0].data
        dim = data.shape
    
    return data, dim, photfnu, R_e

def save_cutout(sky_ra, sky_dec, angular_size, data, wcs, outfile, exposure,
                photfnu, scale, rms, r_e, vmin=None, vmax=None, show=False) :
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
    wcs : TYPE
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
    
    position = SkyCoord(sky_ra, sky_dec, unit='deg', frame='icrs')
    size = u.Quantity(angular_size.value, u.arcsec) # size along each axis
    
    cutout = Cutout2D(data, position, 2*size, wcs=wcs) # cutout will have
        # radius of 'size'
    
    hdu = fits.PrimaryHDU(cutout.data)
    hdr = hdu.header
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
    hdu.writeto(outfile)
    
    if show :
        plt.display_image_simple(cutout.data, vmin=vmin, vmax=vmax)
    
    return

# still necessary?
# determine_final_objects('abell_370', 'id', 0.375)
# determine_final_objects('abell_S1063', 'id', 0.348)
# determine_final_objects('abell_2744', 'id', 0.308, plot=True)
# determine_final_objects('macs_0416', 'id', 0.392)
# determine_final_objects('macs_0717', 'id', 0.5458)
# determine_final_objects('macs_1149', 'id', 0.543)

# still necessary?
# combine_UVJ('abell2744clu_catalogs', 'abell2744clu_v3.9', 153, 155, 161)
# plt.hst_transmission_curves(build_filter_table())
# build_filter_table()
