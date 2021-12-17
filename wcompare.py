from __future__ import print_function
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from math import sqrt, pi

w1sgk = fits.open('w1rot.fits')[0].data
w3sgk = fits.open('w3rot.fits')[0].data
w5sgk = fits.open('w5rot.fits')[0].data
rsgk = fits.open('rad.fits')[0].data

w1antia = np.loadtxt('rotgongs414.1d23')[:,1]
w3antia = np.loadtxt('rotgongs414.3d23n')[:,1]
w5antia = np.loadtxt('rotgongs414.5d23n')[:,1]
rantia = np.loadtxt('rotgongs414.5d23n')[:,0]

maxratio = abs(w3antia).max()/(abs(w3sgk)/rsgk).max()
print("mag ratio, s3 = %10.5e" %(maxratio))
plt.plot(rsgk, abs(w3sgk)/rsgk, 'r')
plt.plot(rsgk, abs(w3sgk)/rsgk*maxratio, 'g')
plt.plot(rantia, abs(w3antia), 'b')
plt.show()
plt.close()

maxratio = abs(w5antia).max()/(abs(w5sgk)/rsgk).max()
print("mag ratio, s5 %10.5e" %(maxratio))
plt.plot(rsgk, abs(w5sgk)/rsgk, 'r')
plt.plot(rsgk, abs(w5sgk)/rsgk*maxratio, 'g')
plt.plot(rantia, abs(w5antia), 'b')
plt.show()
plt.close()

maxratio = abs(w1antia).max()/(abs(w1sgk)/rsgk).max()
print("mag ratio, s1 = %10.5e" %(maxratio))
plt.plot(rsgk, abs(w1sgk)/rsgk, 'r')
plt.plot(rsgk, abs(w1sgk)/rsgk*maxratio, 'g')
plt.plot(rantia, abs(w1antia), 'b')
plt.show()
