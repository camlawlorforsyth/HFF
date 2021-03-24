
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import cm
from matplotlib.colors import LogNorm
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', category=UserWarning)

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

def histogram(data, label, bins=None) :
    
    global currentFig
    fig = plt.figure(currentFig)
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    if bins :
        ax.hist(data, bins=bins, log=True)
    else :
        ax.hist(data, log=True)
    
    ax.set_xlabel('%s' % label, fontsize = 15)
    
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
    f225w_wave = table['F225W_wave']
    f275w_wave = table['F275W_wave']
    f336w_wave = table['F336W_wave']
    f390w_wave = table['F390W_wave']
    f435w_wave = table['F435W_wave']
    f475w_wave = table['F475W_wave']
    f555w_wave = table['F555W_wave']
    f606w_wave = table['F606W_wave']
    f625w_wave = table['F625W_wave']
    f775w_wave = table['F775W_wave']
    f814w_wave = table['F814W_wave']
    f850lp_wave = table['F850LP_wave']
    f105w_wave = table['F105W_wave']
    f110w_wave = table['F110W_wave']
    f125w_wave = table['F125W_wave']
    f140w_wave = table['F140W_wave']
    f160w_wave = table['F160W_wave']
    
    # transmission arrays
    f225w_thru = table['F225W_thru']
    f275w_thru = table['F275W_thru']
    f336w_thru = table['F336W_thru']
    f390w_thru = table['F390W_thru']
    f435w_thru = table['F435W_thru']
    f475w_thru = table['F475W_thru']
    f555w_thru = table['F555W_thru']
    f606w_thru = table['F606W_thru']
    f625w_thru = table['F625W_thru']
    f775w_thru = table['F775W_thru']
    f814w_thru = table['F814W_thru']
    f850lp_thru = table['F850LP_thru']
    f105w_thru = table['F105W_thru']
    f110w_thru = table['F110W_thru']
    f125w_thru = table['F125W_thru']
    f140w_thru = table['F140W_thru']
    f160w_thru = table['F160W_thru']
    
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
                                 xmin=1750, xmax=17500, ymin=0, ymax=1.8)
    
    return

def plot_objects(array_of_xs, array_of_ys, redshift, length, labels, markers,
                 colors, xlabel=None, ylabel=None, title=None,
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
    length : float
        The number of final objects to be analyzed.
    labels : list
        Labels for each object type.
    markers : list
        Marker types to use when plotting.
    colors : list
        Colors to use when plotting.
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
    fig = plt.figure(currentFig, figsize=(9.5, 5.5))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    for i in range(len(array_of_xs)) :
        if colors[i] == 'cyan':
            edgecolor = 'k'
            new_label = (labels[i] + ' (%g)') % length
            size = 100
        else :
            edgecolor = None
            new_label = labels[i]
            size = 60
        ax.scatter(array_of_xs[i], array_of_ys[i], s=size,
                   c=colors[i], marker=markers[i], edgecolors=edgecolor,
                   label = new_label, zorder=i+3)
    
    ax.axhline(redshift, color='k', ls='--', label=r'$z$=%g' %redshift)
    ax.axhline(redshift-0.05, color='b', ls='--', zorder=0,
               label=r'$z \pm 0.05$')
    ax.axhline(redshift+0.05, color='b', ls='--', zorder=0)
    
    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)
    if title :
        strings = title.split('_')
        if strings[0] == 'abell' :
            plot_title = 'Abell ' + strings[1]
        if strings[0] == 'macs' :
            plot_title = 'MACS J' + strings[1]
    ax.set_title(plot_title, fontsize=18)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.legend(bbox_to_anchor=(1.01, 1),
              facecolor='whitesmoke', framealpha=1, fontsize=15,
              loc='upper left')
    
    plt.tight_layout()
    plt.show()
    
    return

def plot_transmission_curves(wavelengths, throughputs, labels, colors,
                             text_x, text_y, xlabel, ylabel,
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
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)
    
    plt.tight_layout()
    plt.show()
    # plt.legend()
    
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
        strings = title.split('_')
        if strings[0] == 'abell' :
            plot_title = 'Abell ' + strings[1]
        if strings[0] == 'macs' :
            plot_title = 'MACS J' + strings[1]
    ax.set_title(plot_title, fontsize=18)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
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
        strings = title.split('_')
        if strings[0] == 'abell' :
            plot_title = 'Abell ' + strings[1]
        if strings[0] == 'macs' :
            plot_title = 'MACS J' + strings[1]
    ax.set_title(plot_title, fontsize=18)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    plt.tight_layout()
    plt.show()
    
    return
