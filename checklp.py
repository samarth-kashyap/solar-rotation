import numpy as np
from astropy.io import fits
from scipy.integrate import simps
from math import sqrt, pi

lp0 = fits.open('lp_n00.fits')[0].data
lp1 = fits.open('lp_n01.fits')[0].data
lp2 = fits.open('lp_n02.fits')[0].data
lp3 = fits.open('lp_n03.fits')[0].data
lp4 = fits.open('lp_n04.fits')[0].data
lp5 = fits.open('lp_n05.fits')[0].data
lp6 = fits.open('lp_n06.fits')[0].data

th = np.linspace(2, 178, 89)*pi/180.
costh = np.cos(th)

l00 = simps(lp0*lp0, costh)
l11 = simps(lp1*lp1, costh)
l22 = simps(lp2*lp2, costh)
l33 = simps(lp3*lp3, costh)
l44 = simps(lp4*lp4, costh)
l55 = simps(lp5*lp5, costh)

l15 = simps(lp1*lp5, costh)
l23 = simps(lp2*lp3, costh)
print(l00, 2.)
print(l11, 2./3)
print(l22, 2./5)
print(l33, 2./7)
print(l44, 2./9)
print(l55, 2./11)

print(l15, 0.)
print(l23, 0.)


