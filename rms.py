
import numpy as np

import core
import plotting

data = core.open_cutout('a2744/cutouts/a2744_ID_3859_f275w.fits', simple=True)

exp = 22720.0

plotting.display_image_simple(data*exp, norm=None)

vals = exp*data.flatten()

plotting.histogram(vals, 'vals')

rms = np.sqrt(np.nanmean(np.square(vals)))
print(rms)
