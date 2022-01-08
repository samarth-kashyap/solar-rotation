from __future__ import print_function
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from math import sqrt, pi

OM = 2.096367060263e-05
unit_conv = 1e-9/OM

#---------------------------------- reading files ------------------------------------
rhmi = fits.open('output_files/rad-hmi.fits')[0].data

w1hmi = fits.open('output_files/w1rot-hmi.fits')[0].data*unit_conv
w3hmi = fits.open('output_files/w3rot-hmi.fits')[0].data*unit_conv
w5hmi = fits.open('output_files/w5rot-hmi.fits')[0].data*unit_conv

w1interp = fits.open('output_files/w1_interp-hmi.fits')[0].data.flatten()*unit_conv
w3interp = fits.open('output_files/w3_interp-hmi.fits')[0].data.flatten()*unit_conv
w5interp = fits.open('output_files/w5_interp-hmi.fits')[0].data.flatten()*unit_conv

w1up_interp = fits.open('output_files/w1up_interp-hmi.fits')[0].data.flatten()*unit_conv
w3up_interp = fits.open('output_files/w3up_interp-hmi.fits')[0].data.flatten()*unit_conv
w5up_interp = fits.open('output_files/w5up_interp-hmi.fits')[0].data.flatten()*unit_conv

w1lo_interp = fits.open('output_files/w1lo_interp-hmi.fits')[0].data.flatten()*unit_conv
w3lo_interp = fits.open('output_files/w3lo_interp-hmi.fits')[0].data.flatten()*unit_conv
w5lo_interp = fits.open('output_files/w5lo_interp-hmi.fits')[0].data.flatten()*unit_conv

e1interp = fits.open('output_files/e1_interp-hmi.fits')[0].data.flatten()*unit_conv
e3interp = fits.open('output_files/e3_interp-hmi.fits')[0].data.flatten()*unit_conv
e5interp = fits.open('output_files/e5_interp-hmi.fits')[0].data.flatten()*unit_conv
rinterp = np.loadtxt('input_files/r_jesper.dat').flatten()
#----------------------------------------------------------------------------------------
wnew = np.zeros((3, len(rinterp)))
wup_new = np.zeros((3, len(rinterp)))
wlo_new = np.zeros((3, len(rinterp)))
enew = np.zeros((3, len(rinterp)))

wnew[0, :] = w1interp
wnew[1, :] = w3interp
wnew[2, :] = w5interp

wup_new[0, :] = w1up_interp
wup_new[1, :] = w3up_interp
wup_new[2, :] = w5up_interp

wlo_new[0, :] = w1lo_interp
wlo_new[1, :] = w3lo_interp
wlo_new[2, :] = w5lo_interp

enew[0, :] = e1interp
enew[1, :] = e3interp
enew[2, :] = e5interp

fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(5, 10))
axs = axs.flatten()
axs[0].plot(rhmi, w1hmi, 'b', label='HMI')
axs[0].plot(rinterp, w1interp, '--r', label='Interpolated')
axs[0].plot(rinterp, w1up_interp, '--g', label='Interpolated - up')
axs[0].plot(rinterp, w1lo_interp, '--b', label='Interpolated - lo')
axs[0].set_title('w1')
axs[0].set_xlabel('r/R')
axs[0].legend()

axs[1].plot(rhmi, w3hmi, 'b', label='HMI')
axs[1].plot(rinterp, w3interp, '--r', label='Interpolated')
axs[1].plot(rinterp, w3up_interp, '--g', label='Interpolated - up')
axs[1].plot(rinterp, w3lo_interp, '--b', label='Interpolated - lo')
axs[1].set_title('w3')
axs[1].set_xlabel('r/R')
axs[1].legend()

axs[2].plot(rhmi, w5hmi, 'b', label='HMI')
axs[2].plot(rinterp, w5interp, '--r', label='Interpolated')
axs[2].plot(rinterp, w5up_interp, '--g', label='Interpolated - up')
axs[2].plot(rinterp, w5lo_interp, '--b', label='Interpolated - lo')
axs[2].set_title('w5')
axs[2].set_xlabel('r/R')
axs[2].legend()
fig.tight_layout()
fig.show()

np.savetxt("output_files/w_hmi.dat", wnew)
np.savetxt("output_files/wup_hmi.dat", wup_new)
np.savetxt("output_files/wlo_hmi.dat", wlo_new)
np.savetxt("output_files/err_hmi.dat", enew)
