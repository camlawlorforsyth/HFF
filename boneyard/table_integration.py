
import numpy as np
from astropy.table import Table
import astropy.units as u

infile = 'm717/m717_ID_2156_photometry.fits'

table = Table.read(infile)
binNum = len(table) - 1

sma = table['sma'][binNum]
smb = table['smb'][binNum]
flux = np.sum(table['flux'])
err = np.sqrt(np.sum(np.square(table['err'])))
SN = flux/err
nPixels = np.sum(table['nPixels'])
width = table['width'][binNum]
PA = table['PA'][binNum]
R_e = (table['R_e'][binNum])*(table['R_e'].unit)
z = table['z'][binNum]
lumDist = (table['lumDist'][binNum])*(table['lumDist'].unit)

f225w_flux = np.sum(table['f225w_flux'])
f225w_err = np.sqrt(np.sum(np.square(table['f225w_err'])))
f225w_nPix = np.sum(table['f225w_nPix'])

f275w_flux = np.sum(table['f275w_flux'])
f275w_err = np.sqrt(np.sum(np.square(table['f275w_err'])))
f275w_nPix = np.sum(table['f275w_nPix'])

f336w_flux = np.sum(table['f336w_flux'])
f336w_err = np.sqrt(np.sum(np.square(table['f336w_err'])))
f336w_nPix = np.sum(table['f336w_nPix'])

f390w_flux = np.sum(table['f390w_flux'])
f390w_err = np.sqrt(np.sum(np.square(table['f390w_err'])))
f390w_nPix = np.sum(table['f390w_nPix'])

f435w_flux = np.sum(table['f435w_flux'])
f435w_err = np.sqrt(np.sum(np.square(table['f435w_err'])))
f435w_nPix = np.sum(table['f435w_nPix'])

f475w_flux = np.sum(table['f475w_flux'])
f475w_err = np.sqrt(np.sum(np.square(table['f475w_err'])))
f475w_nPix = np.sum(table['f475w_nPix'])

f555w_flux = np.sum(table['f555w_flux'])
f555w_err = np.sqrt(np.sum(np.square(table['f555w_err'])))
f555w_nPix = np.sum(table['f555w_nPix'])

f606w_flux = np.sum(table['f606w_flux'])
f606w_err = np.sqrt(np.sum(np.square(table['f606w_err'])))
f606w_nPix = np.sum(table['f606w_nPix'])

f625w_flux = np.sum(table['f625w_flux'])
f625w_err = np.sqrt(np.sum(np.square(table['f625w_err'])))
f625w_nPix = np.sum(table['f625w_nPix'])

f775w_flux = np.sum(table['f775w_flux'])
f775w_err = np.sqrt(np.sum(np.square(table['f775w_err'])))
f775w_nPix = np.sum(table['f775w_nPix'])

f814w_flux = np.sum(table['f814w_flux'])
f814w_err = np.sqrt(np.sum(np.square(table['f814w_err'])))
f814w_nPix = np.sum(table['f814w_nPix'])

f850lp_flux = np.sum(table['f850lp_flux'])
f850lp_err = np.sqrt(np.sum(np.square(table['f850lp_err'])))
f850lp_nPix = np.sum(table['f850lp_nPix'])

f105w_flux = np.sum(table['f105w_flux'])
f105w_err = np.sqrt(np.sum(np.square(table['f105w_err'])))
f105w_nPix = np.sum(table['f105w_nPix'])

f110w_flux = np.sum(table['f110w_flux'])
f110w_err = np.sqrt(np.sum(np.square(table['f110w_err'])))
f110w_nPix = np.sum(table['f110w_nPix'])

f125w_flux = np.sum(table['f125w_flux'])
f125w_err = np.sqrt(np.sum(np.square(table['f125w_err'])))
f125w_nPix = np.sum(table['f125w_nPix'])

f140w_flux = np.sum(table['f140w_flux'])
f140w_err = np.sqrt(np.sum(np.square(table['f140w_err'])))
f140w_nPix = np.sum(table['f140w_nPix'])

f160w_flux = np.sum(table['f160w_flux'])
f160w_err = np.sqrt(np.sum(np.square(table['f160w_err'])))
f160w_nPix = np.sum(table['f160w_nPix'])

values = [[0], [sma], [smb], [flux], [err], [SN], [nPixels], [width],
          [PA], [R_e], [z], [lumDist],
          [f225w_flux], [f225w_err], [f225w_nPix],
          [f275w_flux], [f275w_err], [f275w_nPix],
          [f336w_flux], [f336w_err], [f336w_nPix],
          [f390w_flux], [f390w_err], [f390w_nPix],
          [f435w_flux], [f435w_err], [f435w_nPix],
          [f475w_flux], [f475w_err], [f475w_nPix],
          [f555w_flux], [f555w_err], [f555w_nPix],
          [f606w_flux], [f606w_err], [f606w_nPix],
          [f625w_flux], [f625w_err], [f625w_nPix],
          [f775w_flux], [f775w_err], [f775w_nPix],
          [f814w_flux], [f814w_err], [f814w_nPix],
          [f850lp_flux], [f850lp_err], [f850lp_nPix],
          [f105w_flux], [f105w_err], [f105w_nPix],
          [f110w_flux], [f110w_err], [f110w_nPix],
          [f125w_flux], [f125w_err], [f125w_nPix],
          [f140w_flux], [f140w_err], [f140w_nPix],
          [f160w_flux], [f160w_err], [f160w_nPix]]

integrated_table = Table(values, names = tuple(table.columns))
integrated_table.write('m717/m717_ID_2156_integrated.fits')
