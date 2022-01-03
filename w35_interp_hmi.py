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


def fix_extrapolation(wk, interp_func):
    wk[rmin_idx:rmax_idx] = interp_func(r[rmin_idx:rmax_idx])
    wk[rmax_idx:] = wk[rmax_idx-1]
    wk[:rmin_idx] = wk[rmin_idx]
    return wk

# r = np.squeeze(fits.open('radius.fits')[0].data)
r = np.loadtxt('input_files/r_jesper.dat')
w1mesh = np.squeeze(fits.open('output_files/w1rot-hmi.fits')[0].data)
w3mesh = np.squeeze(fits.open('output_files/w3rot-hmi.fits')[0].data)
w5mesh = np.squeeze(fits.open('output_files/w5rot-hmi.fits')[0].data)
e1mesh = np.squeeze(fits.open('output_files/w1err-hmi.fits')[0].data)
e3mesh = np.squeeze(fits.open('output_files/w3err-hmi.fits')[0].data)
e5mesh = np.squeeze(fits.open('output_files/w5err-hmi.fits')[0].data)
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
eint1 = interp.interp1d(rmesh, e1mesh)
eint3 = interp.interp1d(rmesh, e3mesh)
eint5 = interp.interp1d(rmesh, e5mesh)

w1 = np.zeros(r.shape[0])
w3 = np.zeros(r.shape[0])
w5 = np.zeros(r.shape[0])
e1 = np.zeros(r.shape[0])
e3 = np.zeros(r.shape[0])
e5 = np.zeros(r.shape[0])

w1 = fix_extrapolation(w1, wint1)
w3 = fix_extrapolation(w3, wint3)
w5 = fix_extrapolation(w5, wint5)
e1 = fix_extrapolation(e1, eint1)
e3 = fix_extrapolation(e3, eint3)
e5 = fix_extrapolation(e5, eint5)

writefitsfile(w1, 'output_files/w1_interp-hmi.fits')
writefitsfile(w3, 'output_files/w3_interp-hmi.fits')
writefitsfile(w5, 'output_files/w5_interp-hmi.fits')

writefitsfile(e1, 'output_files/e1_interp-hmi.fits')
writefitsfile(e3, 'output_files/e3_interp-hmi.fits')
writefitsfile(e5, 'output_files/e5_interp-hmi.fits')
