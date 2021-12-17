from __future__ import print_function
import scipy.interpolate as interp
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import os

def writefitsfile(a, fname):
    print("Writing "+fname)
    from astropy.io import fits
    hdu = fits.PrimaryHDU()
    hdu.data = a
    hdu.writeto(fname, overwrite=True)
    return None

r = np.squeeze(fits.open('radius.fits')[0].data)
w1mesh = np.squeeze(fits.open('w1rot-hmi.fits')[0].data)
w3mesh = np.squeeze(fits.open('w3rot-hmi.fits')[0].data)
w5mesh = np.squeeze(fits.open('w5rot-hmi.fits')[0].data)
rmesh = np.squeeze(fits.open('rad-hmi.fits')[0].data)

rmin_idx = np.where(r > rmesh[0])[0][0]
rmax_idx = np.where(r > rmesh[-1])[0][0]

wint1 = interp.interp1d(rmesh, w1mesh)
wint3 = interp.interp1d(rmesh, w3mesh)
wint5 = interp.interp1d(rmesh, w5mesh)
w1 = np.zeros(r.shape[0])
w3 = np.zeros(r.shape[0])
w5 = np.zeros(r.shape[0])

w1[rmin_idx:rmax_idx] = wint1(r[rmin_idx:rmax_idx])
w1[rmax_idx:] = w1[rmax_idx-1]
w1[:rmin_idx] = w1[rmin_idx]

w3[rmin_idx:rmax_idx] = wint3(r[rmin_idx:rmax_idx])
w3[rmax_idx:] = w3[rmax_idx-1]
w3[:rmin_idx] = w3[rmin_idx]

w5[rmin_idx:rmax_idx] = wint5(r[rmin_idx:rmax_idx])
w5[rmax_idx:] = w5[rmax_idx-1]
w5[:rmin_idx] = w5[rmin_idx]

writefitsfile(w1, '/home/g.samarth/rotation/w1_interp-hmi.fits')
writefitsfile(w3, '/home/g.samarth/rotation/w3_interp-hmi.fits')
writefitsfile(w5, '/home/g.samarth/rotation/w5_interp-hmi.fits')
