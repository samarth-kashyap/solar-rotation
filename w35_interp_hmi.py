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
#---------------------- reading the fits files ------------------------------
r = np.loadtxt('input_files/r_jesper.dat')
rmesh = np.squeeze(fits.open('output_files/rad-hmi.fits')[0].data)

w1mesh = np.squeeze(fits.open('output_files/w1rot-hmi.fits')[0].data)
w3mesh = np.squeeze(fits.open('output_files/w3rot-hmi.fits')[0].data)
w5mesh = np.squeeze(fits.open('output_files/w5rot-hmi.fits')[0].data)

w1up_mesh = np.squeeze(fits.open('output_files/w1up-hmi.fits')[0].data)
w3up_mesh = np.squeeze(fits.open('output_files/w3up-hmi.fits')[0].data)
w5up_mesh = np.squeeze(fits.open('output_files/w5up-hmi.fits')[0].data)

w1lo_mesh = np.squeeze(fits.open('output_files/w1lo-hmi.fits')[0].data)
w3lo_mesh = np.squeeze(fits.open('output_files/w3lo-hmi.fits')[0].data)
w5lo_mesh = np.squeeze(fits.open('output_files/w5lo-hmi.fits')[0].data)

e1mesh = np.squeeze(fits.open('output_files/w1err-hmi.fits')[0].data)
e3mesh = np.squeeze(fits.open('output_files/w3err-hmi.fits')[0].data)
e5mesh = np.squeeze(fits.open('output_files/w5err-hmi.fits')[0].data)
#-----------------------------------------------------------------------------

try:
    rmin_idx = np.where(r >= rmesh[0])[0][0]
except IndexError:
    rmin_idx = 0

try:
    rmax_idx = np.where(r >= rmesh[-1])[0][0]
except IndexError:
    rmax_idx = len(r)

#----------------------- creating interpolants --------------------------------
wint1 = interp.interp1d(rmesh, w1mesh)
wint3 = interp.interp1d(rmesh, w3mesh)
wint5 = interp.interp1d(rmesh, w5mesh)

wint1up = interp.interp1d(rmesh, w1up_mesh)
wint3up = interp.interp1d(rmesh, w3up_mesh)
wint5up = interp.interp1d(rmesh, w5up_mesh)

wint1lo = interp.interp1d(rmesh, w1lo_mesh)
wint3lo = interp.interp1d(rmesh, w3lo_mesh)
wint5lo = interp.interp1d(rmesh, w5lo_mesh)

eint1 = interp.interp1d(rmesh, e1mesh)
eint3 = interp.interp1d(rmesh, e3mesh)
eint5 = interp.interp1d(rmesh, e5mesh)
#------------------------------------------------------------------------------

w1 = np.zeros(r.shape[0])
w3 = np.zeros(r.shape[0])
w5 = np.zeros(r.shape[0])

w1up = np.zeros(r.shape[0])
w3up = np.zeros(r.shape[0])
w5up = np.zeros(r.shape[0])

w1lo = np.zeros(r.shape[0])
w3lo = np.zeros(r.shape[0])
w5lo = np.zeros(r.shape[0])

e1 = np.zeros(r.shape[0])
e3 = np.zeros(r.shape[0])
e5 = np.zeros(r.shape[0])

#------------------------------------------------------------------------------

w1 = fix_extrapolation(w1, wint1)
w3 = fix_extrapolation(w3, wint3)
w5 = fix_extrapolation(w5, wint5)

w1up = fix_extrapolation(w1up, wint1up)
w3up = fix_extrapolation(w3up, wint3up)
w5up = fix_extrapolation(w5up, wint5up)

w1lo = fix_extrapolation(w1lo, wint1lo)
w3lo = fix_extrapolation(w3lo, wint3lo)
w5lo = fix_extrapolation(w5lo, wint5lo)

e1 = fix_extrapolation(e1, eint1)
e3 = fix_extrapolation(e3, eint3)
e5 = fix_extrapolation(e5, eint5)
#------------------------------------------------------------------------------

writefitsfile(w1, 'output_files/w1_interp-hmi.fits')
writefitsfile(w3, 'output_files/w3_interp-hmi.fits')
writefitsfile(w5, 'output_files/w5_interp-hmi.fits')

writefitsfile(w1up, 'output_files/w1up_interp-hmi.fits')
writefitsfile(w3up, 'output_files/w3up_interp-hmi.fits')
writefitsfile(w5up, 'output_files/w5up_interp-hmi.fits')

writefitsfile(w1lo, 'output_files/w1lo_interp-hmi.fits')
writefitsfile(w3lo, 'output_files/w3lo_interp-hmi.fits')
writefitsfile(w5lo, 'output_files/w5lo_interp-hmi.fits')

writefitsfile(e1, 'output_files/e1_interp-hmi.fits')
writefitsfile(e3, 'output_files/e3_interp-hmi.fits')
writefitsfile(e5, 'output_files/e5_interp-hmi.fits')
