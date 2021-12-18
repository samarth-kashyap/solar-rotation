from __future__ import print_function
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from math import sqrt, pi

OM = 2.096367060263e-05
unit_conv = 1e-9/OM

w1hmi = fits.open('output_files/w1rot-hmi.fits')[0].data*unit_conv
w3hmi = fits.open('output_files/w3rot-hmi.fits')[0].data*unit_conv
w5hmi = fits.open('output_files/w5rot-hmi.fits')[0].data*unit_conv
rhmi = fits.open('output_files/rad-hmi.fits')[0].data

w1interp = fits.open('output_files/w1_interp-hmi.fits')[0].data.flatten()*unit_conv
w3interp = fits.open('output_files/w3_interp-hmi.fits')[0].data.flatten()*unit_conv
w5interp = fits.open('output_files/w5_interp-hmi.fits')[0].data.flatten()*unit_conv
# rinterp = fits.open('input_files/radius.fits')[0].data.flatten()
rinterp = np.loadtxt('input_files/r_jesper.dat').flatten()

wnew = np.zeros((3, len(rinterp)))
wnew[0, :] = w1interp
wnew[1, :] = w3interp
wnew[2, :] = w5interp

fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(5, 10))
axs = axs.flatten()
axs[0].plot(rhmi, w1hmi, 'b', label='HMI')
axs[0].plot(rinterp, w1interp, '--r', label='Interpolated')
axs[0].set_title('w1')
axs[0].set_xlabel('r/R')
axs[0].legend()

axs[1].plot(rhmi, w3hmi, 'b', label='HMI')
axs[1].plot(rinterp, w3interp, '--r', label='Interpolated')
axs[1].set_title('w3')
axs[1].set_xlabel('r/R')
axs[1].legend()

axs[2].plot(rhmi, w5hmi, 'b', label='HMI')
axs[2].plot(rinterp, w5interp, '--r', label='Interpolated')
axs[2].set_title('w5')
axs[2].set_xlabel('r/R')
axs[2].legend()
fig.tight_layout()
fig.show()

np.savetxt("output_files/w_hmi.dat", wnew)
