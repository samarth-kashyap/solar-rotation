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

# r = np.squeeze(fits.open('radius.fits')[0].data)
r = np.loadtxt('input_files/r_jesper.dat')
w1mesh = np.squeeze(fits.open('output_files/w1rot-hmi.fits')[0].data)
w3mesh = np.squeeze(fits.open('output_files/w3rot-hmi.fits')[0].data)
w5mesh = np.squeeze(fits.open('output_files/w5rot-hmi.fits')[0].data)
rmesh = np.squeeze(fits.open('output_files/rad-hmi.fits')[0].data)


try:
    rmin_idx = np.where(r >= rmesh[0])[0][0]
except IndexError:
    rmin_idx = 0

try:
    rmax_idx = np.where(r >= rmesh[-1])[0][0]
except IndexError:
    rmax_idx = len(r)

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

writefitsfile(w1, 'output_files/w1_interp-hmi.fits')
writefitsfile(w3, 'output_files/w3_interp-hmi.fits')
writefitsfile(w5, 'output_files/w5_interp-hmi.fits')
