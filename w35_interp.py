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
	try:
		os.system("rm "+fname) # removes already existing fits file
	except:
		pass
	hdu.data = a
	hdu.writeto(fname)
	return 0;

r = np.squeeze(fits.open('radius.fits')[0].data)
w1w = np.squeeze(fits.open('w1rot.fits')[0].data)
w3w = np.squeeze(fits.open('w3rot.fits')[0].data)
w5w = np.squeeze(fits.open('w5rot.fits')[0].data)
w3r = np.squeeze(fits.open('rad.fits')[0].data)

print("w3Rmin = %5.5e, w3Rmax = %5.5e" %(w3r.min(), w3r.max()))
print("Rmin = %5.5e, Rmax = %5.5e" %(r.min(), r.max()))
wint1 = interp.interp1d(w3r, w1w)
wint3 = interp.interp1d(w3r, w3w)
wint5 = interp.interp1d(w3r, w5w)
w1 = np.zeros(r.shape[0])
w3 = np.zeros(r.shape[0])
w5 = np.zeros(r.shape[0])

w1[130:7050] = wint1(r[130:7050])
w1[7050:] = w1[7049]
w1[:130] = w1[130]

w3[130:7050] = wint3(r[130:7050])
w3[7050:] = w3[7049]
w3[:130] = w3[130]

w5[130:7050] = wint5(r[130:7050])
w5[7050:] = w5[7049]
w5[:130] = w5[130]

plt.plot(r, w1)
plt.plot(r, w3)
plt.plot(r, w5)
plt.show()

writefitsfile(w1, '/home/samarth/rotation/w1_interp.fits')
writefitsfile(w3, '/home/samarth/rotation/w3_interp.fits')
writefitsfile(w5, '/home/samarth/rotation/w5_interp.fits')
