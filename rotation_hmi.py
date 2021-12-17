import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from math import sqrt, pi
from scipy.integrate import simps

def writefitsfile(a, fname, overwrite=False):
    from astropy.io import fits
    hdu = fits.PrimaryHDU()
    hdu.data = a
    hdu.writeto(fname, overwrite=overwrite)

def pLegSave(l, x, filePrefix, overwrite=False):
    LP0 = np.ones(x.shape);
    LP1 = x;
    LP2 = LP0;
    filename = filePrefix+str(0).zfill(2)+'.fits'
    writefitsfile(LP0, filename, overwrite=overwrite)
    filename = filePrefix+str(1).zfill(2)+'.fits'
    writefitsfile(LP1, filename, overwrite=overwrite)
    for i in range(2,l):
        LP2 = ((2*i - 1)*x*LP1 - (i-1)*LP0)/float(i);
        filename = filePrefix+str(i).zfill(2)+'.fits'
        writefitsfile(LP2, filename, overwrite=overwrite)
        LP0 = LP1;
        LP1 = LP2;

if __name__=="__main__":
    rmesh = np.loadtxt('rmesh.hmi')[::4]
    lenr = len(rmesh)

    rot2d_upper = np.loadtxt('rot2d.hmi')
    rot2d_lower = rot2d_upper[:, ::-1][:, :-1]
    rot2d = np.zeros((lenr, rot2d_upper.shape[1]*2-1))

    lent = rot2d_upper.shape[1]
    rot2d[:, :lent-1] = rot2d_lower
    rot2d[:, lent-1:] = rot2d_upper

    tmesh = np.arange(rot2d_upper.shape[1])
    theta = 90 - tmesh*15./8
    theta = np.append(theta[::-1][:-1], 180. - theta)
    theta = theta*pi/180.

    costh = np.cos(theta)
    sinth = np.sin(theta)

    writefitsfile(rot2d, 'rot2dfull-hmi.fits', overwrite=True)
    pLegSave(15, costh, 'lp_n', overwrite=True)

    omega0 = np.zeros(lenr)
    omega2 = np.zeros(lenr)
    omega4 = np.zeros(lenr)
    omega6 = np.zeros(lenr)
    omega8 = np.zeros(lenr)
    omega10 = np.zeros(lenr)

    w1 = np.zeros(lenr)
    w3 = np.zeros(lenr)
    w5 = np.zeros(lenr)

    lp0 = fits.open('lp_n00.fits')[0].data
    lp2 = fits.open('lp_n02.fits')[0].data
    lp4 = fits.open('lp_n04.fits')[0].data
    lp6 = fits.open('lp_n06.fits')[0].data
    lp8 = fits.open('lp_n08.fits')[0].data
    lp10 = fits.open('lp_n10.fits')[0].data

    for i in range(lenr):
        omega0[i] = simps(rot2d[i, :]*lp0, costh)*(2*0+1.)/2.
        omega2[i] = simps(rot2d[i, :]*lp2, costh)*(2*2+1.)/2.
        omega4[i] = simps(rot2d[i, :]*lp4, costh)*(2*4+1.)/2.
        omega6[i] = simps(rot2d[i, :]*lp6, costh)*(2*6+1.)/2.
        omega8[i] = simps(rot2d[i, :]*lp8, costh)*(2*8+1.)/2.
        omega10[i] = simps(rot2d[i, :]*lp10, costh)*(2*10+1.)/2.
    
    w1 = 2*sqrt(pi/3.)*rmesh*(omega0 - omega2/5.)
    w3 = 2*sqrt(pi/7.)*rmesh*(omega2/5. - omega4/9.)
    w5 = 2*sqrt(pi/11)*rmesh*(omega4/9.)

    writefitsfile(w1, 'w1rot-hmi.fits', overwrite=True)
    writefitsfile(w3, 'w3rot-hmi.fits', overwrite=True)
    writefitsfile(w5, 'w5rot-hmi.fits', overwrite=True)
    writefitsfile(rmesh, 'rad-hmi.fits', overwrite=True)

    plt.plot(rmesh, w1)
    plt.show()
