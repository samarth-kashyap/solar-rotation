import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from math import sqrt, pi
from scipy.integrate import simps

def writefitsfile(a, fname):
    from astropy.io import fits
    hdu = fits.PrimaryHDU()
    hdu.data = a
    hdu.writeto(fname)

def pLegSave(l, x, filePrefix):
    LP0 = np.ones(x.shape);
    LP1 = x;
    LP2 = LP0;
    filename = filePrefix+str(0).zfill(2)+'.fits'
    writefitsfile(LP0, filename)
    filename = filePrefix+str(1).zfill(2)+'.fits'
    writefitsfile(LP1, filename)
    for i in range(2,l):
        LP2 = ((2*i - 1)*x*LP1 - (i-1)*LP0)/float(i);
        filename = filePrefix+str(i).zfill(2)+'.fits'
        writefitsfile(LP2, filename)
        LP0 = LP1;
        LP1 = LP2;

if __name__=="__main__":
    temp = np.loadtxt('rot2d.gong8')
    r = temp[:,0]
    rotLen = temp.shape[1]-1
    rotLen = 2*rotLen -1

    theta = np.linspace(2, 178, rotLen)
    theta = theta*pi/180.
    costh = np.cos(theta)
    sinth = np.sin(theta)

    rot = np.zeros((temp.shape[0], rotLen))
    for i in xrange(temp.shape[1]-1):
        rot[:, rotLen/2+i] = temp[:, i+1]
        rot[:, rotLen/2-i] = temp[:, i+1]

    writefitsfile(rot, 'rot2dfull.fits')
    pLegSave(15, costh, 'lp_n')

    omega0 = np.zeros(temp.shape[0])
    omega2 = np.zeros(temp.shape[0])
    omega0new = np.zeros(temp.shape[0])
    omega2new = np.zeros(temp.shape[0])
    omega4 = np.zeros(temp.shape[0])
    omega6 = np.zeros(temp.shape[0])
    omega8 = np.zeros(temp.shape[0])
    omega10 = np.zeros(temp.shape[0])
    w1 = np.zeros(temp.shape[0])
    w3 = np.zeros(temp.shape[0])
    w5 = np.zeros(temp.shape[0])

    lp0 = fits.open('lp_n00.fits')[0].data
    lp2 = fits.open('lp_n02.fits')[0].data
    lp4 = fits.open('lp_n04.fits')[0].data
    lp6 = fits.open('lp_n06.fits')[0].data
    lp8 = fits.open('lp_n08.fits')[0].data
    lp10 = fits.open('lp_n10.fits')[0].data

    for i in xrange(temp.shape[0]):
        omega0[i] = simps(rot[i,:]*lp0, costh)*(2*0+1.)/2.
        omega2[i] = simps(rot[i,:]*lp2, costh)*(2*2+1.)/2.
        omega4[i] = simps(rot[i,:]*lp4, costh)*(2*4+1.)/2.
        omega6[i] = simps(rot[i,:]*lp6, costh)*(2*6+1.)/2.
        omega8[i] = simps(rot[i,:]*lp8, costh)*(2*8+1.)/2.
        omega10[i] = simps(rot[i,:]*lp10, costh)*(2*10+1.)/2.
    
    plt.plot(omega0)
    plt.plot(omega0new)
    plt.show()
    w1 = 2*sqrt(pi/3.)*r*(omega0 - omega2/5.)
    w3 = 2*sqrt(pi/7.)*r*(omega2/5. - omega4/9.)
    w5 = 2*sqrt(pi/11)*r*(omega4/9.)

    writefitsfile(w1, 'w1rot.fits')
    writefitsfile(w3, 'w3rot.fits')
    writefitsfile(w5, 'w5rot.fits')
    writefitsfile(r, 'rad.fits')
