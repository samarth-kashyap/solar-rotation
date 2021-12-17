from __future__ import print_function
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from math import sqrt, pi

w1hmi = fits.open('w1rot-hmi.fits')[0].data
w3hmi = fits.open('w3rot-hmi.fits')[0].data
w5hmi = fits.open('w5rot-hmi.fits')[0].data
rhmi = fits.open('rad-hmi.fits')[0].data

w1interp = fits.open('w1_interp-hmi.fits')[0].data.flatten()
w3interp = fits.open('w3_interp-hmi.fits')[0].data.flatten()
w5interp = fits.open('w5_interp-hmi.fits')[0].data.flatten()
rinterp = fits.open('radius.fits')[0].data.flatten()

fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(5, 15))
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
fig.show()
