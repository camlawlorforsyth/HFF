
import copy
import numpy as np
import matplotlib.pyplot as plt

import corner
from matplotlib import cm
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import prospect.io.read_results as reader

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
    
    cmap = copy.copy(cmap)
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

def display_cutouts(cutout_data, nrows, ncols, filters, flags, outfile,
                    save=False) :
    '''
    Display all the cutouts for a given galaxy, saving the resulting file if
    desired.
    
    Parameters
    ----------
    cutout_data : list
        List of cutout data, with the same length as filters.
    nrows : int
        Number of subplot rows.
    ncols : int
        Number of subplot columns.
    filters : list
        List of filters used for the observations.
    flags : list
        List of flag values, one per filter.
    outfile : string
        The name of the output file.
    save : bool, optional
        Flag to save the complete image. The default is False.
    
    Returns
    -------
    None.
    
    '''
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(2.5*ncols, 2.5*nrows))
    currentFig += 1
    plt.clf()
    
    position = range(1, len(cutout_data) + 1)
    
    # old version
    # cmap = copy.copy(cm.inferno)
    # cmap.set_bad('white', 1)
    # cmap.set_under('white', 1)
    
    for i in range(len(cutout_data)) :
        if flags[i] == 0 :
            fw = 'normal'
        else :
            fw = 'bold'
        
        ax = fig.add_subplot(nrows, ncols, position[i])
        data = cutout_data[i]
        
        # new version
        pos = np.ma.masked_array(data, data <= 0)
        neg = np.ma.masked_array(data, data > 0)
        ax.imshow(pos, origin='lower', norm=LogNorm(), cmap='inferno')
        ax.imshow(-neg, origin='lower', norm=LogNorm(), cmap='gray')
        
        # old version
        # ax.imshow(data, origin='lower', norm=LogNorm(), cmap=cmap)
        
        ax.axis('off')
        ax.set_title('{} : {}'.format(filters[i], flags[i]), fontsize=13,
                     fontweight=fw)
    
    if save :
        plt.tight_layout()
        plt.savefig(outfile, bbox_inches='tight')
        plt.close()
    else :
        plt.tight_layout()
        plt.show()
    
    return

def display_hists(cutout_data, nrows, ncols, filters, medians, num_bins,
                  outfile, save=False) :
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(2.4*ncols, 2.4*nrows))
    currentFig += 1
    plt.clf()
    
    position = range(1, len(cutout_data) + 1)
    
    for i in range(len(cutout_data)) :
        ax = fig.add_subplot(nrows, ncols, position[i])
        data = cutout_data[i]
        
        ax.hist(data, bins=num_bins[i], color='k', histtype='step')
        
        ax.set_title('{}:{:.2g}'.format(filters[i], medians[i]), fontsize=13)
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=8)
    
    plt.tight_layout()
    
    if save :
        plt.savefig(outfile, bbox_inches='tight')
        plt.close()
    else :
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
    
    cmap = copy.copy(cm.gray)
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
    
    cmap = copy.copy(cmap)
    cmap.set_bad(bad, 1)
    
    frame = ax.imshow(data, origin='lower', cmap=cmap, norm=norm,
                      vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(frame)
    cbar.set_label(cbar_label, fontsize=15)
    
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

def histogram(data, label, title=None, bins=None, log=False, histtype='bar',
              vlines=[], colors=[], labels=[]) :
    
    global currentFig
    fig = plt.figure(currentFig)
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    if bins and not log :
        ax.hist(data, bins=bins, color='k', histtype=histtype)
    elif bins and log :
        ax.hist(data, bins=bins, log=log, color='k', histtype=histtype)
    elif log and not bins :
        ax.hist(data, log=log, color='k', histtype=histtype)
    else :
        ax.hist(data, histtype=histtype)
    
    if len(vlines) > 0 :
        for i in range(len(vlines)) :
            ax.axvline(vlines[i], ls='--', color=colors[i], lw=1, alpha=0.5,
                        label=labels[i])
    
    ax.set_xlabel('{}'.format(label), fontsize = 15)
    if title :
        if title[0] == 'a' :
            title = 'Abell ' + title[1:]
        if title[0] == 'm' :
            title = 'MACS J' + title[1:]
    ax.set_title(title, fontsize=18)
    
    if len(vlines) > 0 :
        ax.legend(loc='upper left', facecolor='whitesmoke', framealpha=1,
                  fontsize=15)
    
    plt.tight_layout()
    plt.show()
    
    return

def histogram_multi(data, label, bins=[], log=False, histtype='step',
                    colors=[], labels=[], styles=[],
                    xmin=None, xmax=None, ymin=None, ymax=None) :
    
    global currentFig
    fig = plt.figure(currentFig)
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    if bins and not log :
        for i in range(len(data)) :
            ax.hist(data[i], bins=bins[i], color=colors[i], histtype=histtype,
                    label=labels[i], linestyle=styles[i])
    elif bins and log :
        for i in range(len(data)) :
            ax.hist(data[i], bins=bins[i], color=colors[i], histtype=histtype,
                    label=labels[i], log=log, linestyle=styles[i])
    elif log and not bins :
        for i in range(len(data)) :
            ax.hist(data[i], color=colors[i], histtype=histtype, log=log,
                    label=labels[i], linestyle=styles[i])
    else :
        for i in range(len(data)) :
            ax.hist(data[i], color=colors[i], histtype=histtype,
                    label=labels[i], linestyle=styles[i])
    
    ax.set_xlabel('{}'.format(label), fontsize = 15)
    # ax.set_xscale('log')
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.legend(loc='upper right', facecolor='whitesmoke', framealpha=1,
              fontsize=15)
    
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
    text_x = [2000, 2975, 4350, 5900, 7350, 8650, 12100,
              2375, 3525, 4975, 10150, 13650,
              3950, 5575, 7700, 11100, 15050] # horizontal label positions
    text_y = [1.32, 1.42, 1.35, 1.4, 1.4, 1.26, 1.45,
              0.75, 0.68, 0.75, 0.8, 0.85,
              0.13, 0.2, 0.15, 0.225, 0.25] # vertical label positions
    
    if plot :
        plot_transmission_curves(wavelengths, filters, labels, colors,
                                 text_x, text_y, r'Wavelength ($\mathrm{\AA}$)',
                                 'Integrated System Throughput',
                                 xmin=1750, xmax=17500, ymin=0, ymax=1.8,
                                 ytick_labels=[0.0, 0.2, 0.4, 0.6,
                                               0.2, 0.4, 0.6, 0.2, 0.4, 0.6])
    
    return

def plot_chains(result, results_type='dynesty') :
    
    if results_type == 'emcee' :
        chosen = np.random.choice(result['run_params']['nwalkers'], size=10,
                                  replace=False)
        tracefig = reader.traceplot(result, figsize=(11, 6), chains=chosen)
    else :
        tracefig = reader.traceplot(result, figsize=(11, 6))
    
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

def plot_corner(samples, labels, npar, result, ranges=None,
                figsizewidth=9, figsizeheight=9, outfile=None, save=False) :
    
    '''
    cornerfig = reader.subcorner(result, start=0, thin=1,
                                  fig=plt.subplots(npar, npar,
                                    figsize=(1.5*npar, 1.5*npar))[0],
                                      range=ranges)
    '''
    
    cornerfig = corner.corner(samples, labels=labels, range=ranges,
                              fig=plt.subplots(npar, npar,
                                               figsize=(2*npar, 2*npar))[0],
                              fill_contours=True, plot_datapoints=False,
                              plot_density=False, quantiles=[0.16, 0.5, 0.84],
                              show_titles=True)
    
    if save :
        plt.savefig(outfile)
        plt.close()
    else :
        plt.show()
    
    return

def plot_degeneracies(list_of_xs, list_of_xerrs_lo, list_of_xerrs_hi,
                      list_of_ys, list_of_yerrs_lo, list_of_yerrs_hi,
                      list_of_best_ys, list_of_slopes, list_of_intercepts,
                      q_xi, q_yi, q_z, sf_xi, sf_yi, sf_z, V_J, FUV_V,
                      labels='', xlabel=None, list_of_ylabels=None, title=None,
                      xmin=None, xmax=None, outfile=None, loc='best',
                      figsizewidth=16, figsizeheight=9, save=False) :
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(figsizewidth, figsizeheight))
    currentFig += 1
    plt.clf()
    
    spec = gridspec.GridSpec(ncols=2, nrows=3, figure=fig, width_ratios=[2, 1])
    
    xs = np.linspace(xmin, xmax, 1000)
    
    ax1 = fig.add_subplot(spec[0, 0])
    ax2 = fig.add_subplot(spec[1, 0], sharex=ax1)
    ax3 = fig.add_subplot(spec[2, 0], sharex=ax1)
    ax4 = fig.add_subplot(spec[0:2, 1])
    
    ax1.plot(list_of_xs[0], list_of_best_ys[0], '--', color='darkred',
             label=labels[0])
    ax1.plot(xs, list_of_slopes[0]*xs + list_of_intercepts[0], ':',
             color='maroon', label=labels[1])
    ax1.errorbar(list_of_xs[0], list_of_ys[0],
                 xerr=(list_of_xerrs_lo[0], list_of_xerrs_hi[0]),
                 yerr=(list_of_yerrs_lo[0], list_of_yerrs_hi[0]),
                 color='red', ecolor='k', elinewidth=0.7, linestyle='-',
                 marker='.', markeredgecolor='k', markerfacecolor='k',
                 markersize=11, label=labels[2], zorder=3)
    ax1.set_ylabel(list_of_ylabels[0], fontsize=15)
    ax1.tick_params(axis='x', which='major', labelbottom=False)
    
    ax2.plot(list_of_xs[1], list_of_best_ys[1], '--', color='darkblue')
    ax2.plot(xs, list_of_slopes[1]*xs + list_of_intercepts[1], ':',
             color='navy')
    ax2.errorbar(list_of_xs[1], list_of_ys[1],
                 xerr=(list_of_xerrs_lo[1], list_of_xerrs_hi[1]),
                 yerr=(list_of_yerrs_lo[1], list_of_yerrs_hi[1]),
                 color='blue', ecolor='k', elinewidth=0.7, linestyle='-',
                 marker='.', markeredgecolor='k', markerfacecolor='k',
                 markersize=11, zorder=3)
    ax2.set_ylabel(list_of_ylabels[1], fontsize=15)
    ax2.tick_params(axis='x', which='major', labelbottom=False)
    
    ax3.plot(list_of_xs[2], list_of_best_ys[2], '--', color='darkgrey')
    ax3.plot(xs, list_of_slopes[2]*xs + list_of_intercepts[2], ':', color='k')
    ax3.errorbar(list_of_xs[2], list_of_ys[2],
                 xerr=(list_of_xerrs_lo[2], list_of_xerrs_hi[2]),
                 yerr=(list_of_yerrs_lo[2], list_of_yerrs_hi[2]),
                 color='grey', ecolor='k', elinewidth=0.7, linestyle='-',
                 marker='.', markeredgecolor='k', markerfacecolor='k',
                 markersize=11, zorder=3)
    ax3.set_ylabel(list_of_ylabels[2], fontsize=15)
    ax3.set_xlabel(xlabel, fontsize=15)
    
    # location on the FUVVJ diagram
    ax4.contour(q_xi, q_yi, q_z, colors='darkred', alpha=0.2, zorder=2)
    ax4.contour(sf_xi, sf_yi, sf_z, colors='darkblue', alpha=0.2, zorder=2)
    ax4.set_xlabel(r'V$-$J')
    ax4.set_ylabel(r'FUV$-$V')
    ax4.set_xlim(0, 2.1)
    ax4.set_ylim(0, 8.4)
    ax4.plot(V_J, FUV_V, marker='o', markerfacecolor='red',
             markeredgecolor='k', markersize=9)
    slope, intercept, horiz, vert = 3.24, 0.32, 3.45, 1.56
    first_knee, second_knee = (horiz - intercept)/slope, vert
    divide_x = np.linspace(first_knee, second_knee, 1000)
    divide_y = slope*divide_x + intercept
    ax4.hlines(horiz, xmin, first_knee, ls='-', color='k', zorder=5)
    ax4.plot(divide_x, divide_y, 'k-', zorder=5)
    ax4.vlines(vert, slope*second_knee + intercept, 8.4,
               ls='-', color='k', zorder=5)
    
    ax1.set_title(title, fontsize=18)
    ax1.set_xlim(xmin, xmax)
    ax1.legend(loc=loc, facecolor='whitesmoke', framealpha=1, fontsize=11)
    
    plt.tight_layout()
    
    if save :
        plt.savefig(outfile, bbox_inches='tight')
        plt.close()
    else :
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
        if title == 'a1063' :
            title = 'Abell S' + title[1:]
        elif title == 'm1149' :
            title = 'MACS J' + title[1:]
        elif title[0] == 'a' :
            title = 'Abell ' + title[1:]
        elif title[0] == 'm' :
            title = 'MACS J0' + title[1:]
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

def plot_sed(xs, list_of_ys, labels, colors, outfile, xlabel=None, ylabel=None,
             xmin=None, xmax=None, ymin=None, ymax=None, save=False,
             figsizewidth=12, figsizeheight=9, style='--', xlog=False) :
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(figsizewidth, figsizeheight))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    for i in range(len(list_of_ys)) :
        ax.plot(xs, list_of_ys[i], style, color=colors[i], label=labels[i])
    
    ax.set_yscale('log')
    if xlog :
        ax.set_xscale('log')
    
    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.legend(facecolor='whitesmoke', framealpha=1, fontsize=8)
    
    plt.tight_layout()
    
    if save :
        plt.savefig(outfile, bbox_inches='tight')
        plt.close()
    else :
        plt.show()
    
    return

def plot_sed_alone(obs) :
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(16, 9))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    waves, mask = obs['phot_wave'], obs['phot_mask']
    fluxes, e_fluxes = obs['maggies'], obs['maggies_unc']
    
    xmin, xmax = np.min(waves)*0.8, np.max(waves)/0.8
    ymin, ymax = np.max([np.min(fluxes), 1e-15])*0.8, np.max(fluxes)/0.4
    
    ax.plot(waves, fluxes, label='all observed photometry',
            marker='o', markersize=12, alpha=0.8, ls='', lw=3,
            color='slateblue')
    
    ax.errorbar(waves[mask], fluxes[mask], yerr=e_fluxes[mask],
                label='photometry to fit', marker='o', markersize=8, alpha=0.8,
                ls='', lw=3, ecolor='tomato', mfc='none', mec='tomato', mew=3)
    
    for filt in obs['filters'] :
        wave, trans = filt.wavelength.copy(), filt.transmission.copy()
        trans = trans / trans.max()
        trans = 10**(0.2*(np.log10(ymax/ymin)))*trans*ymin
        ax.loglog(wave, trans, lw=3, color='gray', alpha=0.7)
    
    ax.set(xscale='log', yscale='log', xlim=(xmin, xmax), ylim=(ymin, ymax))
    ax.set_xlabel(r'Wavelength ($\rm \AA$)', fontsize=15)
    ax.set_ylabel('Flux Density (maggies)', fontsize=15)
    ax.legend(loc='upper left', facecolor='whitesmoke', framealpha=1,
              fontsize=15)
    
    plt.tight_layout()
    plt.show()
    
    return

def plot_sed_and_model(obs, model, init_phot, init_spec, sps) :
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(16, 9))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    waves, mask = obs['phot_wave'], obs['phot_mask']
    fluxes, e_fluxes = obs['maggies'], obs['maggies_unc']
    
    xmin, xmax = np.min(waves)*0.8, np.max(waves)/0.8
    ymin, ymax = np.max([np.min(fluxes), 1e-15])*0.8, np.max(fluxes)/0.4
    
    ax.plot(waves, fluxes, label='all observed photometry',
            marker='o', markersize=12, alpha=0.8, ls='', lw=3,
            color='slateblue')
    
    ax.errorbar(waves[mask], fluxes[mask], yerr=e_fluxes[mask],
                label='photometry to fit', marker='o', markersize=8, alpha=0.8,
                ls='', lw=3, ecolor='tomato', mfc='none', mec='tomato', mew=3)
    
    ax.plot(waves, init_phot, label='initial model photometry',
            marker='s', markersize=10, alpha=0.6, ls='', lw=3, mfc='none',
            mec='gray', mew=3)
    
    model_waves = sps.wavelengths*(1.0 + model.params.get('zred', 0.0))
    ax.plot(model_waves, init_spec, label='model spectrum', lw=0.7,
            color='navy', alpha=0.7)
    
    ax.set(xscale='log', yscale='log')
    ax.set_xlabel(r'Wavelength ($\rm \AA$)', fontsize=15)
    ax.set_ylabel('Flux Density (maggies)', fontsize=15)
    ax.legend(loc='upper left', facecolor='whitesmoke', framealpha=1,
              fontsize=15)
    
    plt.tight_layout()
    plt.show()
    
    return

def plot_sed_and_models(obs, model, init_phot, init_spec, opt_phot, opt_spec,
                        sps) :
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(16, 9))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    waves, mask = obs['phot_wave'], obs['phot_mask']
    fluxes, e_fluxes = obs['maggies'], obs['maggies_unc']
    
    xmin, xmax = np.min(waves)*0.8, np.max(waves)/0.8
    ymin, ymax = np.max([np.min(fluxes), 1e-15])*0.8, np.max(fluxes)/0.4
    
    ax.plot(waves, fluxes, label='all observed photometry',
            marker='o', markersize=12, alpha=0.8, ls='', lw=3,
            color='slateblue')
    
    ax.errorbar(waves[mask], fluxes[mask], yerr=e_fluxes[mask],
                label='photometry to fit', marker='o', markersize=8, alpha=0.8,
                ls='', lw=3, ecolor='tomato', mfc='none', mec='tomato', mew=3)
    
    ax.plot(waves, init_phot, label='initial model photometry',
            marker='s', markersize=10, alpha=0.6, ls='', lw=3, mfc='none',
            mec='gray', mew=3)
    
    ax.plot(waves, opt_phot, label='optimized model photometry',
            marker='s', markersize=10, alpha=1, ls='', lw=3, mfc='none',
            mec='c', mew=3)
    
    model_waves = sps.wavelengths*(1.0 + model.params.get('zred', 0.0))
    ax.plot(model_waves, init_spec, label='initial model spectrum', lw=0.7,
            color='navy', alpha=0.7)
    
    ax.plot(model_waves, opt_spec, label='optimized model spectrum', lw=0.7,
            color='c', alpha=1)
    
    ax.set(xscale='log', yscale='log')
    ax.set_xlabel(r'Wavelength ($\rm \AA$)', fontsize=15)
    ax.set_ylabel('Flux Density (maggies)', fontsize=15)
    ax.legend(loc='upper left', facecolor='whitesmoke', framealpha=1,
              fontsize=15)
    
    plt.tight_layout()
    plt.show()
    
    return

def plot_sed_from_fit(waves, fluxes, e_fluxes, mask, map_spec, map_phot,
                      model_waves, chisq='', lo=None, hi=None, title=None,
                      outfile=None, save=False) :
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(16, 9))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    xmin, xmax = 0.8*np.min(waves), 1.2*np.max(waves)
    if np.min(fluxes) < 0 :
        ymin = 0.8*np.min(map_phot)
    else :
        ymin = 0.8*np.min([np.min(fluxes), np.min(map_phot)])
    ymax = 1.2*np.max(fluxes)
    
    ax.plot(waves, fluxes, label='all observed photometry',
            marker='o', markersize=12, alpha=0.8, ls='', lw=3,
            color='slateblue')
    
    ax.errorbar(waves[mask], fluxes[mask], yerr=e_fluxes[mask],
                label='photometry to fit', marker='o', markersize=8, alpha=0.8,
                ls='', lw=3, ecolor='tomato', mfc='none', mec='tomato', mew=3)
    
    ax.plot(waves, map_phot, label='MAP model photometry',
            marker='s', markersize=10, alpha=1, ls='', lw=3, mfc='none',
            mec='r', mew=3)
    
    ax.plot(model_waves, map_spec, label='MAP model spectrum', lw=0.7,
            color='navy', alpha=0.7)
    
    # ax.fill_between(model_waves, lo, hi, color='lightgrey', alpha=0.2)
    
    ax.annotate(r'$\chi^2_{\rm reduced}$ = ' + '{:.2f}'.format(chisq),
                (0.75, 0.25), fontsize=15, xycoords='axes fraction')
    
    ax.set(xscale='log', yscale='log', xlim=(xmin, xmax), ylim=(ymin, ymax))
    
    ax.set_title(title, fontsize=18)
    ax.set_xlabel(r'Wavelength ($\rm \AA$)', fontsize=15)
    ax.set_ylabel('Flux Density (maggies)', fontsize=15)
    ax.legend(loc='upper left', facecolor='whitesmoke', framealpha=1,
              fontsize=15)
    
    plt.tight_layout()
    
    if save :
        plt.savefig(outfile, bbox_inches='tight')
        plt.close()
    else :
        plt.show()
    
    return

def plot_simple(xs, ys, yerr, xerr=None, label='',
                xlabel=None, ylabel=None, title=None,
                xmin=None, xmax=None, ymin=None, ymax=None,
                figsizewidth=9, figsizeheight=6) :
    
    global currentFig
    fig = plt.figure(currentFig, figsize=(figsizewidth, figsizeheight))
    currentFig += 1
    plt.clf()
    ax = fig.add_subplot(111)
    
    ax.plot(xs, ys, 'r-', zorder=1)
    ax.errorbar(xs, ys, xerr=xerr, yerr=yerr, fmt='k.', label=label, zorder=2,
                elinewidth=1)
    # ax.plot(xs, ys, 'ko', label=label, zorder=2)
    
    ax.set_yscale('log')
    # ax.set_xscale('log')
    
    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    if label :
        ax.legend(facecolor='whitesmoke', framealpha=1, fontsize=15)
    
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
        ax.text(text_x[i], text_y[i], s=labels[i], fontsize=15)
    for i in range(7, 12) :
        # ax.plot(wavelengths[i], 0.6+throughputs[i], '-', color=colors[i],
        #         label=labels[i])
        ax.fill_between(wavelengths[i], 0.6, 0.6+throughputs[i],  
                        color=colors[i], alpha=0.4)
        ax.text(text_x[i], text_y[i], s=labels[i], fontsize=15)
    for i in range(12, len(wavelengths)) :
        # ax.plot(wavelengths[i], throughputs[i], '-', color=colors[i],
        #         label=labels[i])
        ax.fill_between(wavelengths[i], 0, throughputs[i],
                        color=colors[i], alpha=0.4)
        ax.text(text_x[i], text_y[i], s=labels[i], fontsize=15)
    
    ax.axhline(0.0, ls='-', color='k', lw=1.5)
    ax.axhline(0.6, ls='-', color='k', lw=1.5)
    ax.axhline(1.2, ls='-', color='k', lw=1.5)
    
    if ytick_labels is not None :
        ax.set_yticks(np.arange(0, 2, 0.2))
        ax.set_yticklabels(np.asarray(ytick_labels, dtype=np.str))
    
    ax.tick_params(axis='both', which='major', labelsize=15)
    
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    
    ax.set_xlabel(xlabel, fontsize=22)
    ax.set_ylabel(ylabel, fontsize=22)
    
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
