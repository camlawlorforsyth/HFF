
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import cm
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

currentFig = 1

def display_annuli(data, xy, rins, eta, theta, bad='black', cbar_label='',
                   cmap=cm.gray, norm=LogNorm(), vmin=None, vmax=None) :
    '''
    Display the given image with various ellipses.
    
    Parameters
    ----------
    data : numpy.ndarray
        The data that will displayed.
    xy : tuple
        Coordinates of the ellipse centers.
    rins : list
        List of semi-major axis values for the inner edge of each ellipse.
    eta : float
        Ellipticity of the ellipses.
    theta : float
        Position angle of the ellipses in degrees, measured anti-clockwise.
    bad : string, optional
        The 'bad' value to use when displaying the image. Default is 'black'.
    cbar_label : string, optional
        The colorbar label. Default is Empty.
    cmap : matplotlib.colors.LinearSegmentedColormap, optional
        The colormap to use. The default is cm.gray.
    norm : matplotlib.colors.LogNorm, optional
        The normalization to use. The default is LogNorm().
    vmin : float, optional
        The minimum value to use for colormap scaling. The default is None.
    vmax : float, optional
        The maximum value to use for colormap scaling. The default is None.
    
    Returns
    -------
    None.
    
    '''
    import matplotlib.patches as mpatches
    
    global currentFig
    fig = plt.figure(currentFig)
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    cmap = cmap
    cmap.set_bad(bad, 1)
    
    frame = ax.imshow(data, origin='lower', cmap=cmap, norm=norm,
                      vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(frame)
    cbar.set_label(cbar_label)
    
    for rin in rins :
        ellipse = mpatches.Ellipse(xy, 2*rin, 2*(1-eta)*rin, theta,
                                   edgecolor='red', facecolor='none')
        ax.add_patch(ellipse)
    
    plt.tight_layout()
    plt.show()
    
    return

def display_image_with_wcs(data, wcs, vmin=None, vmax=None) :
    '''
    Display the given image data with the corresponding WCS.
    
    Parameters
    ----------
    data : numpy.ndarray
        The data that will displayed.
    wcs : astropy.wcs.WCS
        The world coordinate system to use to display the data.
    vmin : float, optional
        The minimum value to use for colormap scaling. The default is None.
    vmax : float, optional
        The maximum value to use for colormap scaling. The default is None.
    
    Returns
    -------
    None.
    
    '''
    
    global currentFig
    fig = plt.figure(currentFig)
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111, projection=wcs)
    
    ax.coords[0].set_axislabel('Right Ascension')
    ax.coords[1].set_axislabel('Declination')
    
    cmap = cm.gray
    cmap.set_bad('black', 1)
    
    frame = ax.imshow(data, origin='lower', norm=LogNorm(), cmap=cmap,
                      vmin=vmin, vmax=vmax)
    
    plt.tight_layout()
    plt.show()
    
    return

def display_image_simple(data, bad='black', cbar_label='', cmap=cm.gray, 
                         norm=LogNorm(), vmin=None, vmax=None) :
    '''
    Display the given image data without the WCS.
    
    Parameters
    ----------
    data : numpy.ndarray
        The data that will displayed.
    bad : string, optional
        The 'bad' value to use when displaying the image. Default is 'black'.
    cbar_label : string, optional
        The colorbar label. Default is Empty.
    cmap : matplotlib.colors.LinearSegmentedColormap, optional
        The colormap to use. The default is cm.gray.
    vmin : float, optional
        The minimum value to use for colormap scaling. The default is None.
    vmax : float, optional
        The maximum value to use for colormap scaling. The default is None.
    
    Returns
    -------
    None.
    
    '''
    
    global currentFig
    fig = plt.figure(currentFig)
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    cmap = cmap
    cmap.set_bad(bad, 1)
    
    frame = ax.imshow(data, origin='lower', cmap=cmap, norm=norm,
                      vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(frame)
    cbar.set_label(cbar_label)
    
    # for i in range(len(xs)) :
    #     plt.arrow(dim[1]/2, dim[0]/2, xs[i], ys[i], color='r')
    # plt.plot(xs+dim[1]/2, ys+dim[0]/2, 'ro')
    
    plt.tight_layout()
    plt.show()
    
    # imgs = ax.get_images()
    # if len(imgs) > 0 :
    #    vmin, vmax = imgs[0].get_clim()
    #    print(vmin, vmax)
    
    return

def histogram(data, label, title=None, bins=None, log=False) :
    
    global currentFig
    fig = plt.figure(currentFig)
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    if bins and not log :
        ax.hist(data, bins=bins)
    elif bins and log :
        ax.hist(data, bins=bins, log=log)
    elif log and not bins :
        ax.hist(data, log=log)
    else :
        ax.hist(data)
    
    ax.set_xlabel('%s' % label, fontsize = 15)
    if title :
        if title[0] == 'a' :
            title = 'Abell ' + title[1:]
        if title[0] == 'm' :
            title = 'MACS J' + title[1:]
    ax.set_title(title, fontsize=18)
    
    plt.tight_layout()
    plt.show()
    
    return

def hst_transmission_curves(table, plot=True) :
    '''
    Prepare the data to plot the transmission curves for all the filters
    included in the Hubble Frontier Fields Deepspace program. This includes
    17 filters in total.
    
    Parameters
    ----------
    table : astropy.table.Table
        The table which includes all wavelength and transmission arrays.
    plot : bool, optional
        Boolean to plot the resulting filter curves. The default is True.
    
    Returns
    -------
    None.
    
    '''
    
    # wavelength arrays
    f225w_wave = table['f225w_wave']
    f275w_wave = table['f275w_wave']
    f336w_wave = table['f336w_wave']
    f390w_wave = table['f390w_wave']
    f435w_wave = table['f435w_wave']
    f475w_wave = table['f475w_wave']
    f555w_wave = table['f555w_wave']
    f606w_wave = table['f606w_wave']
    f625w_wave = table['f625w_wave']
    f775w_wave = table['f775w_wave']
    f814w_wave = table['f814w_wave']
    f850lp_wave = table['f850lp_wave']
    f105w_wave = table['f105w_wave']
    f110w_wave = table['f110w_wave']
    f125w_wave = table['f125w_wave']
    f140w_wave = table['f140w_wave']
    f160w_wave = table['f160w_wave']
    
    # transmission arrays
    f225w_thru = table['f225w_thru']
    f275w_thru = table['f275w_thru']
    f336w_thru = table['f336w_thru']
    f390w_thru = table['f390w_thru']
    f435w_thru = table['f435w_thru']
    f475w_thru = table['f475w_thru']
    f555w_thru = table['f555w_thru']
    f606w_thru = table['f606w_thru']
    f625w_thru = table['f625w_thru']
    f775w_thru = table['f775w_thru']
    f814w_thru = table['f814w_thru']
    f850lp_thru = table['f850lp_thru']
    f105w_thru = table['f105w_thru']
    f110w_thru = table['f110w_thru']
    f125w_thru = table['f125w_thru']
    f140w_thru = table['f140w_thru']
    f160w_thru = table['f160w_thru']
    
    wavelengths = [f225w_wave, f336w_wave, f475w_wave, f625w_wave, f775w_wave,
                   f850lp_wave, f125w_wave,
                   f275w_wave, f390w_wave, f555w_wave, f105w_wave, f140w_wave,
                   f435w_wave, f606w_wave, f814w_wave, f110w_wave, f160w_wave]
    filters = [f225w_thru, f336w_thru, f475w_thru, f625w_thru, f775w_thru,
               f850lp_thru, f125w_thru,
               f275w_thru, f390w_thru, f555w_thru, f105w_thru, f140w_thru,
               f435w_thru, f606w_thru, f814w_thru, f110w_thru, f160w_thru]
    labels = ['F225W', 'F336W', 'F475W', 'F625W', 'F775W', 'F850LP', 'F125W',
              'F275W', 'F390W', 'F555W', 'F105W', 'F140W',
              'F435W', 'F606W', 'F814W', 'F110W', 'F160W']
    colors = ['hotpink', 'mediumorchid', 'royalblue', 'cyan', 'springgreen',
              'greenyellow', 'darkorange',
              'violet', 'darkviolet', 'dodgerblue', 'yellow', 'orangered',
              'darkslateblue', 'deepskyblue', 'lime', 'gold', 'red']
    text_x = [2100, 3125, 4500, 6050, 7500, 8700, 12250,
              2500, 3700, 5150, 10300, 13850,
              4100, 5650, 7900, 11300, 15200] # horizontal label positions
    text_y = [1.22, 1.28, 1.35, 1.4, 1.4, 1.28, 1.45,
              0.62, 0.7, 0.75, 0.8, 0.85,
              0.15, 0.2, 0.15, 0.225, 0.25] # vertical label positions
    
    if plot :
        plot_transmission_curves(wavelengths, filters, labels, colors,
                                 text_x, text_y, r'Wavelength ($\rm \AA$)',
                                 'Integrated System Throughput',
                                 xmin=1750, xmax=17500, ymin=0, ymax=1.8,
                                 ytick_labels=[0.0, 0.2, 0.4, 0.6,
                                               0.2, 0.4, 0.6, 0.2, 0.4, 0.6])
    
    return

def plot_colorcolor_multi(xs, ys, labels, lengths, colors, markers, sizes,
                          alphas, contour_xs, contour_ys, contour_zs, version,
                          xlabel=None, ylabel=None, plot_divider=True,
                          xmin=None, xmax=None, ymin=None, ymax=None, loc=0) :
    '''
    Plot a standard color-color diagram. This includes the UVJ diagram, as well
    as the FUVVJ diagram.
    
    Parameters
    ----------
    xs : list
        DESCRIPTION.
    ys : list
        DESCRIPTION.
    xlabel : TYPE, optional
        DESCRIPTION. The default is None.
    ylabel : TYPE, optional
        DESCRIPTION. The default is None.
    plot_divider : bool, optional
        Flag to plot the quiescent region dividing line. The default is True.
    xmin : TYPE, optional
        DESCRIPTION. The default is None.
    xmax : TYPE, optional
        DESCRIPTION. The default is None.
    ymin : TYPE, optional
        DESCRIPTION. The default is None.
    ymax : TYPE, optional
        DESCRIPTION. The default is None.
    loc : float, optional
        Location of the legend. The default is 0.
    
    Returns
    -------
    None.
    
    '''
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(9.5, 7))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    # plot the diagonal region that contains quiescent galaxies
    if version == 'FUVVJ' :
        slope, intercept, horiz, vert = 3.24, 0.32, 3.45, 1.56
    else : # UVJ is the default
        slope, intercept, horiz, vert = 0.88, 0.69, 1.3, 1.5
    
    first_knee, second_knee = (horiz - intercept)/slope, vert
    divide_x = np.linspace(first_knee, second_knee, 1000)
    divide_y = slope*divide_x + intercept
    if plot_divider :
        ax.hlines(horiz, xmin, first_knee, ls='-', color='k', zorder=5)
        ax.plot(divide_x, divide_y, 'k-', zorder=5)
        ax.vlines(vert, slope*second_knee + intercept, ymax,
                  ls='-', color='k', zorder=5)
    
    for i in range(len(xs)) :
        ax.scatter(xs[i], ys[i], s=sizes[i], color=colors[i],
                   marker=markers[i], alpha=alphas[i],
                   label='{} ({})'.format(labels[i], lengths[i]))
    
    for i in range(len(contour_xs)) :
        ax.contour(contour_xs[i], contour_ys[i], contour_zs[i], colors='k',
                   alpha=0.1, zorder=2)
    
    # these are two star-forming (according to the FUVVJ) bCGs from A370 and M416
    # ax.scatter([-20.313091957851302 - -21.698459182796444,
    #             -22.290556427098757 - -23.441648553116536],
    #            [-16.005818274869203 - -20.313091957851302,
    #             -19.444035537051526 - -22.290556427098757], s=150, color='k')
    
    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.legend(loc=loc, facecolor='whitesmoke', framealpha=1,
              fontsize=15)
    
    plt.tight_layout()
    plt.show()
    
    return

def plot_objects(array_of_xs, array_of_ys, redshift, labels, markers, colors,
                 sizes, alphas, redshift_tol_lo=0.05, redshift_tol_hi=0.05,
                 xlabel=None, ylabel=None, title=None,
                 xmin=None, xmax=None, ymin=None, ymax=None) :
    '''
    Plot all the objects included in the catalog. This includes final objects
    and other objects that have certain quality flags. Also marks cluster
    boundaries based on redshift.
    
    Parameters
    ----------
    array_of_xs : list
        Arrays of object stellar masses.
    array_of_ys : list
        Arrays of object redshifts.
    redshift : float
        The redshift of the cluster.
    labels : list
        Labels for each object type.
    markers : list
        Marker types to use when plotting.
    colors : list
        Colors to use when plotting.
    redshift_tol_lo : float, optional
        The low redshift tolerance to use. The default is 0.05.
    redshift_tol_hi : float, optional
        The high redshift tolerance to use. The default is 0.05.
    xlabel : string, optional
        X-axis label. The default is None.
    ylabel : string, optional
        Y-axis label. The default is None.
    title : string, optional
        Title of the plot. The default is None.
    xmin : float, optional
        X-axis lower limit. The default is None.
    xmax : float, optional
        X-axis upper limit. The default is None.
    ymin : float, optional
        Y-axis lower limit. The default is None.
    ymax : float, optional
        Y-axis upper limit. The default is None.
    
    Returns
    -------
    None.
    
    '''
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(11, 9))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    # global currentFig
    # fig, ax = plt.subplots()
    # currentFig += 1
    # plt.clf()
    
    for i in range(len(array_of_xs)) :
        ax.scatter(array_of_xs[i], array_of_ys[i], s=sizes[i], alpha=alphas[i],
                   c=colors[i], marker=markers[i], label=labels[i], zorder=i+3)
    
    ax.axhline(redshift, color='k', ls=':', label=r'$z$={}'.format(redshift))
    ax.axhline(redshift - redshift_tol_lo, color='b', ls=':', zorder=0,
               label=r'$z \pm \delta z$')
    ax.axhline(redshift + redshift_tol_hi, color='b', ls=':', zorder=0)
    
    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)
    if title :
        if title[0] == 'a' :
            title = 'Abell ' + title[1:]
        if title[0] == 'm' :
            title = 'MACS J' + title[1:]
    ax.set_title(title, fontsize=18)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.legend(loc='upper right', facecolor='whitesmoke', framealpha=1,
              fontsize=15) # bbox_to_anchor=(1.01, 1)
    
    # inset axes
    axins = inset_axes(ax, width='40%', height='30%', loc='upper left',
                       borderpad=4)
    axins.axhline(redshift, color='k', ls=':')
    axins.axhline(redshift - redshift_tol_lo, color='b', ls=':', zorder=0)
    axins.axhline(redshift + redshift_tol_hi, color='b', ls=':', zorder=0)
    for i in range(len(array_of_xs)) :
        axins.scatter(array_of_xs[i], array_of_ys[i], s=sizes[i],
                      alpha=alphas[i], c=colors[i], marker=markers[i],
                      zorder=i+3)
        axins.set_xlim(7, 12.5)
        axins.set_ylim(redshift-2*redshift_tol_lo, redshift+2*redshift_tol_hi)
    
    
    plt.tight_layout()
    plt.show()
    
    return

def plot_transmission_curves(wavelengths, throughputs, labels, colors,
                             text_x, text_y, xlabel, ylabel, ytick_labels=None,
                             xmin=None, xmax=None, ymin=None, ymax=None,
                             figsizewidth=15, figsizeheight=6.75) :
    '''
    Plot the transmission curves for all included filters, with corresponding
    labels and colors.
    
    Parameters
    ----------
    wavelengths : list
        Arrays of wavelengths for each filter.
    throughputs : list
        Arrays of throughputs for each filter.
    labels : list
        Labels for each filter.
    colors : list
        Colors for each filter.
    text_x : list
        Horizontal label positions.
    text_y : list
        Vertical label positions.
    xlabel : string
        X-axis label. The default is None.
    ylabel : string
        Y-axis label. The default is None.
    ytick_labels : list, optional
        List of values to use as labels for the Y-axis. The default is None.
    xmin : float, optional
        X-axis lower limit. The default is None.
    xmax : float, optional
        X-axis upper limit. The default is None.
    ymin : float, optional
        Y-axis lower limit. The default is None.
    ymax : float, optional
        Y-axis upper limit. The default is None.
    figsizewidth : TYPE, optional
        DESCRIPTION. The default is 15.
    figsizeheight : TYPE, optional
        DESCRIPTION. The default is 6.75.
    
    Returns
    -------
    None.
    
    '''
    
    global currentFig        
    fig = plt.figure(currentFig, figsize=(figsizewidth, figsizeheight))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    for i in range(7) :
        # ax.plot(wavelengths[i], 1.2+throughputs[i], '-', color=colors[i],
        #         label=labels[i])
        ax.fill_between(wavelengths[i], 1.2, 1.2+throughputs[i],
                        color=colors[i], alpha=0.4)
        ax.text(text_x[i], text_y[i], s=labels[i])
    for i in range(7, 12) :
        # ax.plot(wavelengths[i], 0.6+throughputs[i], '-', color=colors[i],
        #         label=labels[i])
        ax.fill_between(wavelengths[i], 0.6, 0.6+throughputs[i],  
                        color=colors[i], alpha=0.4)
        ax.text(text_x[i], text_y[i], s=labels[i])
    for i in range(12, len(wavelengths)) :
        # ax.plot(wavelengths[i], throughputs[i], '-', color=colors[i],
        #         label=labels[i])
        ax.fill_between(wavelengths[i], 0, throughputs[i],
                        color=colors[i], alpha=0.4)
        ax.text(text_x[i], text_y[i], s=labels[i])
    
    ax.axhline(0.0, ls='-', color='k', lw=1.5)
    ax.axhline(0.6, ls='-', color='k', lw=1.5)
    ax.axhline(1.2, ls='-', color='k', lw=1.5)
    
    if ytick_labels is not None :
        ax.set_yticks(np.arange(0, 2, 0.2))
        ax.set_yticklabels(np.asarray(ytick_labels, dtype=np.str))
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)
    
    plt.tight_layout()
    plt.show()
    
    return

def plot_UVJ(xs, ys, xlabel=None, ylabel=None, title=None,
             xmin=None, xmax=None, ymin=None, ymax=None) :
    '''
    Plot the standard U-V versus V-J color-color diagram.
    
    Parameters
    ----------
    xs : TYPE
        DESCRIPTION.
    ys : TYPE
        DESCRIPTION.
    xlabel : TYPE, optional
        DESCRIPTION. The default is None.
    ylabel : TYPE, optional
        DESCRIPTION. The default is None.
    title : TYPE, optional
        DESCRIPTION. The default is None.
    xmin : TYPE, optional
        DESCRIPTION. The default is None.
    xmax : TYPE, optional
        DESCRIPTION. The default is None.
    ymin : TYPE, optional
        DESCRIPTION. The default is None.
    ymax : TYPE, optional
        DESCRIPTION. The default is None.
    
    Returns
    -------
    None.
    
    '''
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(7.5, 5.5))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    slope, intercept, horiz, vert = 0.88, 0.69, 1.3, 1.5
    first_knee, second_knee = (horiz - intercept)/slope, vert
    
    ax.hlines(horiz, xmin, first_knee, ls='-', color='k', zorder=5)
    divide_x = np.linspace(first_knee, second_knee, 1000)
    divide_y = slope*divide_x + intercept
    ax.plot(divide_x, divide_y, 'k-')
    ax.vlines(vert, slope*second_knee + intercept, ymax,
              ls='-', color='k', zorder=5)
    
    # the corner region where the quiescent galaxies reside
    corner = ( ((xs <= first_knee) & (ys >= horiz)) |
               ( ((xs >= first_knee) & (xs <= second_knee))
                 & (ys >= slope*xs + intercept)) )
    
    # divide the objects into their regions
    q_x, q_y = xs[corner], ys[corner]
    sf_x, sf_y = xs[~corner], ys[~corner]
    q_length, sf_length = np.sum(corner), np.sum(~corner)
    
    ax.scatter(q_x, q_y, color='r', marker='o',
               label='Quiescent ({})'.format(q_length))
    ax.scatter(sf_x, sf_y, color='b', marker='o',
               label='Star Forming ({})'.format(sf_length))
    
    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)
    if title :
        if title[0] == 'a' :
            title = 'Abell ' + title[1:]
        if title[0] == 'm' :
            title = 'MACS J' + title[1:]
    ax.set_title(title, fontsize=18)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.legend(loc='upper left', facecolor='whitesmoke', framealpha=1,
              fontsize=15)
    
    plt.tight_layout()
    plt.show()
    
    return

def plot_UVJ_hist(xs, ys, xlabel=None, ylabel=None, title=None,
                  xmin=None, xmax=None, ymin=None, ymax=None) :
    '''
    Plot the standard U-V versus V-J color-color diagram as a 2D histogram.
    
    Parameters
    ----------
    xs : TYPE
        DESCRIPTION.
    ys : TYPE
        DESCRIPTION.
    xlabel : TYPE, optional
        DESCRIPTION. The default is None.
    ylabel : TYPE, optional
        DESCRIPTION. The default is None.
    title : TYPE, optional
        DESCRIPTION. The default is None.
    xmin : TYPE, optional
        DESCRIPTION. The default is None.
    xmax : TYPE, optional
        DESCRIPTION. The default is None.
    ymin : TYPE, optional
        DESCRIPTION. The default is None.
    ymax : TYPE, optional
        DESCRIPTION. The default is None.
    
    Returns
    -------
    None.
    
    '''
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(8, 6))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    ax.scatter(xs, ys, color='b', marker='o')
    
    ax.hlines(1.3, xmin, (1.3-0.59)/0.88, ls='-', color='k', zorder=5)
    divide_x = np.linspace((1.3-0.59)/0.88, 1.5, 100)
    divide_y = 0.88*divide_x + 0.59
    ax.plot(divide_x, divide_y, 'k-')
    ax.vlines(1.5, 0.88*1.5+0.59, ymax, ls='-', color='k', zorder=5)
    
    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)
    if title :
        if title[0] == 'a' :
            title = 'Abell ' + title[1:]
        if title[0] == 'm' :
            title = 'MACS J' + title[1:]
    ax.set_title(title, fontsize=18)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    plt.tight_layout()
    plt.show()
    
    return
