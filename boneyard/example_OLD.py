
import numpy as np

from astropy.table import Table
from astropy.visualization import make_lupton_rgb
from matplotlib import cm

from core import open_cutout
import plotting as plt

import warnings
warnings.filterwarnings('ignore')

cluster = 'm416'
ID = 2876
filt = 'f160w'
val = 10

binDir = '{}/bins'.format(cluster)
cutoutDir = '{}/cutouts'.format(cluster)
outfile = '{}/photometry/{}_ID_{}_photometry.fits'.format(cluster, cluster, ID)

'''
sci_file = '{}/{}_ID_{}_{}.fits'.format(cutoutDir, cluster, ID, filt)
noise_file = '{}/{}_ID_{}_{}_noise.fits'.format(cutoutDir, cluster, ID, filt)
segmap_file = '{}/{}_ID_{}_segmap.fits'.format(cutoutDir, cluster, ID)

# science and segmentation map
sci, dim, photfnu, r_e, redshift, sma, smb, pa = open_cutout(sci_file)
noise, _, _, _, _, _, _, _ = open_cutout(noise_file)
segMap, _, _, _, _, _, _, _ = open_cutout(segmap_file)
# plt.display_image_simple(sci, cmap=cm.spring)
# plt.display_image_simple(segMap, norm=None)

# science image masked based on the segmentation map
new_sci = sci.copy()
new_sci[(segMap > 0) & (segMap != ID)] = 0

newnew_sci = sci.copy()
newnew_sci[(segMap > 0) & (segMap != ID)] = np.nan
# plt.display_image_simple(newnew_sci, norm=None, cmap=cm.spring)

# the annulus mask
bin_path = '{}/{}_ID_{}_annuli.npz'.format(binDir, cluster, ID)
bin_data = np.load(bin_path)
bins_image = bin_data['image']
# plt.display_image_simple(bins_image, norm=None, cmap=cm.gist_rainbow,
#                          cbar_label='Bin Number')

# annulus mask, only showing the last annulus
temp_bins = bins_image.copy()
temp_bins[temp_bins != val] = 0
# plt.display_image_simple(temp_bins, norm=None, vmax=11)

temp_sci = new_sci.copy()
temp_sci[bins_image != val] = np.nan
# plt.display_image_simple(temp_sci, norm=None, cmap=cm.spring)

pix_sci = temp_sci.copy()
# pix_sci[(pix_sci > 0) | (pix_sci < 0)] = np.nan
pix_sci[pix_sci != 0] = np.nan
pix_sci[pix_sci == 0] = 1
# plt.display_image_simple(pix_sci, norm=None, cmap=cm.spring, vmax=1.08)
# print(np.nansum(pix_sci))
'''

'''
ys = 2.5*np.log10(phot['f625w_flux']/phot['f475w_flux'])
yerr = 2.5/np.log(10)*np.sqrt(np.square(phot['f625w_err']/phot['f625w_flux']) +
                              np.square(phot['f475w_err']/phot['f475w_flux']))
# plt.plot_simple(phot['sma']/R_e, ys, yerr, 'k', '',
#                 xlabel=r'Radius ($R_{\rm e}$)',
#                 ylabel=r'$(m_{\mathrm{F435W}} - m_{\mathrm{F625W}}) \approx (U-V)$')
'''

'''
f225, _, _, _, _, _, _, _ = open_cutout('{}/{}_ID_{}_{}.fits'.format(cutoutDir, cluster, ID, 'f225w'))
f275, _, _, _, _, _, _, _ = open_cutout('{}/{}_ID_{}_{}.fits'.format(cutoutDir, cluster, ID, 'f275w'))
f336, _, _, _, _, _, _, _ = open_cutout('{}/{}_ID_{}_{}.fits'.format(cutoutDir, cluster, ID, 'f336w'))
f390, _, _, _, _, _, _, _ = open_cutout('{}/{}_ID_{}_{}.fits'.format(cutoutDir, cluster, ID, 'f390w'))

f435, _, _, _, _, _, _, _ = open_cutout('{}/{}_ID_{}_{}.fits'.format(cutoutDir, cluster, ID, 'f435w'))
f475, _, pfnu_475, _, _, _, _, _ = open_cutout('{}/{}_ID_{}_{}.fits'.format(cutoutDir, cluster, ID, 'f475w'))
f606, _, _, _, _, _, _, _ = open_cutout('{}/{}_ID_{}_{}.fits'.format(cutoutDir, cluster, ID, 'f606w'))
f625, _, pfnu_625, _, _, _, _, _ = open_cutout('{}/{}_ID_{}_{}.fits'.format(cutoutDir, cluster, ID, 'f625w'))
f775, _, _, _, _, _, _, _ = open_cutout('{}/{}_ID_{}_{}.fits'.format(cutoutDir, cluster, ID, 'f775w'))
f814, _, _, _, _, _, _, _ = open_cutout('{}/{}_ID_{}_{}.fits'.format(cutoutDir, cluster, ID, 'f814w'))
f850, _, _, _, _, _, _, _ = open_cutout('{}/{}_ID_{}_{}.fits'.format(cutoutDir, cluster, ID, 'f850lp'))

f105, _, _, _, _, _, _, _ = open_cutout('{}/{}_ID_{}_{}.fits'.format(cutoutDir, cluster, ID, 'f105w'))
f110, _, _, _, _, _, _, _ = open_cutout('{}/{}_ID_{}_{}.fits'.format(cutoutDir, cluster, ID, 'f110w'))
f125, _, _, _, _, _, _, _ = open_cutout('{}/{}_ID_{}_{}.fits'.format(cutoutDir, cluster, ID, 'f125w'))
f140, _, _, _, _, _, _, _ = open_cutout('{}/{}_ID_{}_{}.fits'.format(cutoutDir, cluster, ID, 'f140w'))
f160, _, _, _, _, _, _, _ = open_cutout('{}/{}_ID_{}_{}.fits'.format(cutoutDir, cluster, ID, 'f160w'))

cbar_label = r'$(m_{\mathrm{F435W}} - m_{\mathrm{F625W}}) \approx (U-V)$'
UminusV = 2.5*np.log10((pfnu_625*f625)/(pfnu_475*f475))
UminusV[(segMap > 0) & (segMap != ID)] = np.nan
# plt.display_image_simple(UminusV, norm=None, cbar_label=cbar_label,
#                          cmap=cm.RdBu_r, vmin=0, vmax=2.4)

temp = UminusV.flatten()
temp = temp[np.isfinite(temp)]
bins = int(np.round(np.sqrt(len(temp))))
# plt.histogram(temp, cbar_label, bins=bins, histtype='step')

# rgb = make_lupton_rgb(f105+f110+f125+f140+f160,
#                       f435+f475+f606+f625+f775+f814+f850,
#                       f225+f275+f336+f390, Q=15, stretch=0.3)
# plt.display_image_simple(rgb, norm=None)
'''
