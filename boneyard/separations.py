
import numpy as np

from astropy.coordinates import SkyCoord

a370 = SkyCoord('02h39m52.9s', '-01d34m36.5s')
a370_par = SkyCoord('02h40m13.4s', '-01d37m32.8s')

a1063 = SkyCoord('22h48m44.4s', '-44d31m48.5s')
a1063_par = SkyCoord('22h49m17.7s', '-44d32m43.8s')

a2744 = SkyCoord('00h14m21.2s', '-30d23m50.1s')
a2744_par = SkyCoord('00h13m53.6s', '-30d22m54.3s')

m416 = SkyCoord('04h16m08.9s', '-24d04m28.7s')
m416_par = SkyCoord('04h16m33.1s', '-24d06m48.7s')

m717 = SkyCoord('07h17m34.0s', '+37d44m49.0s')
m717_par = SkyCoord('07h17m17.0s', '+37d49m47.3s')

m1149 = SkyCoord('11h49m36.3s', '+22d23m58.1s')
m1149_par = SkyCoord('11h49m40.5s', '+22d18m02.3s')

clus = np.array([a370, a1063, a2744, m416, m717, m1149])
pars = np.array([a370_par, a1063_par, a2744_par, m416_par, m717_par, m1149_par])

for i in range(len(clus)) :
    sep = clus[i].separation(pars[i]).arcmin
    print('{:.2f}'.format(sep))
