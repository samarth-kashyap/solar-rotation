import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from math import sqrt, pi
from scipy.integrate import simps
from pyshtools import legendre
from scipy.special import legendre as scipy_legendre

NAX = np.newaxis

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

def get_both_hemispheres(w2d_upper):
    w2d_lower = w2d_upper[:, ::-1][:, 1:]
    w2d = np.zeros((lenr, w2d_upper.shape[1]*2-1))

    lent = w2d_upper.shape[1]
    w2d[:, :lent] = w2d_upper
    w2d[:, lent:] = w2d_lower
    return w2d


def wsr_to_omegas(wsr, r, s):
    return wsr/r * np.sqrt((2*s+1)/4./np.pi)


def wsr_to_omega(wsr):
    omega1 = wsr_to_omegas(wsr[0], rmesh, 1)
    omega3 = wsr_to_omegas(wsr[1], rmesh, 3)
    omega5 = wsr_to_omegas(wsr[2], rmesh, 5)
    omega = omega1[:, NAX] * dlegp[:, 0][NAX, :]
    omega += omega3[:, NAX] * dlegp[:, 1][NAX, :] 
    omega += omega5[:, NAX] * dlegp[:, 2][NAX, :]
    return omega

def wsr_to_omega_sbd(wsr):
    omega1 = wsr_to_omegas(wsr[0], rmesh, 1)
    omega3 = wsr_to_omegas(wsr[1], rmesh, 3)
    omega5 = wsr_to_omegas(wsr[2], rmesh, 5)
    Omegasr = [omega1, omega3, omega5]
    # carrying out the sum over angular degree s
    s = np.array([1, 3, 5])
    lens = len(s)
    Omega_r_theta = 0
    for s_ind in range(lens):
        s = 2 * s_ind + 1    # we want odd s only. s = 1, 3, ...
        
        # creating the legendre polynomial. P_s(x)
        leg_poly = scipy_legendre(s)
        # derivative with respect to its argument dP_s(x)/dx
        leg_poly_deriv = leg_poly.deriv()
        Omega_r_theta += Omegasr[s_ind][:, NAX] * leg_poly_deriv(costh)[NAX, :]
    return Omega_r_theta


def esr_to_err_sbd(wsr):
    omega1 = wsr_to_omegas(wsr[0], rmesh, 1)
    omega3 = wsr_to_omegas(wsr[1], rmesh, 3)
    omega5 = wsr_to_omegas(wsr[2], rmesh, 5)
    Omegasr = [omega1, omega3, omega5]
    # carrying out the sum over angular degree s
    s = np.array([1, 3, 5])
    lens = len(s)
    Omega_r_theta = 0
    for s_ind in range(lens):
        s = 2 * s_ind + 1    # we want odd s only. s = 1, 3, ...
        
        # creating the legendre polynomial. P_s(x)
        leg_poly = scipy_legendre(s)
        # derivative with respect to its argument dP_s(x)/dx
        leg_poly_deriv = leg_poly.deriv()
        Omega_r_theta += (Omegasr[s_ind][:, NAX] * leg_poly_deriv(costh)[NAX, :])**2
    return np.sqrt(Omega_r_theta)



if __name__=="__main__":
    rmesh = np.loadtxt('input_files/rmesh.hmi')[::4]
    rot2d = np.loadtxt('input_files/rot2d.hmi')
    err2d = np.loadtxt('input_files/err2d.hmi')
    lenr = len(rmesh)
    
    tmesh = np.arange(rot2d.shape[1])
    theta = 90 - tmesh*15./8
    theta = np.append(-theta, theta[::-1][1:]) + 90.
    theta = theta*pi/180.

    rot2d = get_both_hemispheres(rot2d)
    err2d = get_both_hemispheres(err2d)
    err2d = err2d**2 #the linear operator acts on variance

    rot2dup = rot2d + err2d
    rot2dlo = rot2d - err2d

    costh = np.cos(theta)
    sinth = np.sin(theta)

    legp, dlegp = legendre.PlBar_d1(5, costh[0])/np.sqrt(2)
    for i in range(1, len(costh)):
        lp1, dp1 = legendre.PlBar_d1(5, costh[i])/np.sqrt(2)
        dlegp = np.vstack((dlegp, dp1))
        legp = np.vstack((legp, lp1))
    dlegp = dlegp[:, 1::2]

    writefitsfile(rot2d, 'output_files/rot2dfull-hmi.fits', overwrite=True)
    pLegSave(15, costh, 'data_files/lp_n', overwrite=True)

    #--------------------------- initializing the new profiles ----------------------
    omega0 = np.zeros(lenr)
    omega2 = np.zeros(lenr)
    omega4 = np.zeros(lenr)

    omega0up = np.zeros(lenr)
    omega2up = np.zeros(lenr)
    omega4up = np.zeros(lenr)

    omega0lo = np.zeros(lenr)
    omega2lo = np.zeros(lenr)
    omega4lo = np.zeros(lenr)

    err0 = np.zeros(lenr)
    err2 = np.zeros(lenr)
    err4 = np.zeros(lenr)

    w1 = np.zeros(lenr)
    w1up = np.zeros(lenr)
    w1lo = np.zeros(lenr)

    w3 = np.zeros(lenr)
    w3up = np.zeros(lenr)
    w3lo = np.zeros(lenr)

    w5 = np.zeros(lenr)
    w5up = np.zeros(lenr)
    w5lo = np.zeros(lenr)
    #--------------------------------------------------------------------------------

    ord_list = np.arange(0, 12, 2)
    lp_list = []
    for ord in ord_list:
        _lp = fits.open(f'data_files/lp_n{ord:02d}.fits')[0].data
        lp_list.append(_lp)

    lp0 = fits.open('data_files/lp_n00.fits')[0].data
    lp2 = fits.open('data_files/lp_n02.fits')[0].data
    lp4 = fits.open('data_files/lp_n04.fits')[0].data

    lp0 = scipy_legendre(0)(costh)
    lp2 = scipy_legendre(2)(costh)
    lp4 = scipy_legendre(4)(costh)

    for i in range(lenr):
        omega0[i] = simps(rot2d[i, :]*lp0, costh)*(2*0+1.)/2.
        omega2[i] = simps(rot2d[i, :]*lp2, costh)*(2*2+1.)/2.
        omega4[i] = simps(rot2d[i, :]*lp4, costh)*(2*4+1.)/2.

        omega0up[i] = simps(rot2dup[i, :]*lp0, costh)*(2*0+1.)/2.
        omega2up[i] = simps(rot2dup[i, :]*lp2, costh)*(2*2+1.)/2.
        omega4up[i] = simps(rot2dup[i, :]*lp4, costh)*(2*4+1.)/2.

        omega0lo[i] = simps(rot2dlo[i, :]*lp0, costh)*(2*0+1.)/2.
        omega2lo[i] = simps(rot2dlo[i, :]*lp2, costh)*(2*2+1.)/2.
        omega4lo[i] = simps(rot2dlo[i, :]*lp4, costh)*(2*4+1.)/2.
        
        err0[i] = simps(err2d[i, :]*lp0**2 * sinth**2, theta)* ((2*0+1.)/2.)**2
        err2[i] = simps(err2d[i, :]*lp2**2 * sinth**2, theta)* ((2*2+1.)/2.)**2
        err4[i] = simps(err2d[i, :]*lp4**2 * sinth**2, theta)* ((2*4+1.)/2.)**2
    
    w1 = 2*sqrt(pi/3.)*rmesh*(omega0/1. - omega2/5.)
    w3 = 2*sqrt(pi/7.)*rmesh*(omega2/5. - omega4/9.)
    w5 = 2*sqrt(pi/11)*rmesh*(omega4/9.)

    w1up = 2*sqrt(pi/3.)*rmesh*(omega0up/1. - omega2up/5.)
    w3up = 2*sqrt(pi/7.)*rmesh*(omega2up/5. - omega4up/9.)
    w5up = 2*sqrt(pi/11)*rmesh*(omega4up/9.)

    w1lo = 2*sqrt(pi/3.)*rmesh*(omega0lo/1. - omega2lo/5.)
    w3lo = 2*sqrt(pi/7.)*rmesh*(omega2lo/5. - omega4lo/9.)
    w5lo = 2*sqrt(pi/11)*rmesh*(omega4lo/9.)

    e1 = (2*sqrt(pi/3.))**2 *rmesh**2 *(err0/(1.**2) - err2/(5.**2))
    e3 = (2*sqrt(pi/7.))**2 *rmesh**2 *(err2/(5.**2) - err4/(9.**2))
    e5 = (2*sqrt(pi/11))**2 *rmesh**2 *(err4/(9.**2))

    e1 = np.sqrt(e1)
    e3 = np.sqrt(e3)
    e5 = np.sqrt(e5)

    w1sig = w1up - w1
    w3sig = w3up - w3
    w5sig = w5up - w5
    wsr_sig = [abs(w1sig), abs(w3sig), abs(w5sig)]

    #--------------------------- storing the profiles ----------------------
    writefitsfile(rmesh, 'output_files/rad-hmi.fits', overwrite=True)

    writefitsfile(w1, 'output_files/w1rot-hmi.fits', overwrite=True)
    writefitsfile(w3, 'output_files/w3rot-hmi.fits', overwrite=True)
    writefitsfile(w5, 'output_files/w5rot-hmi.fits', overwrite=True)

    writefitsfile(w1up, 'output_files/w1up-hmi.fits', overwrite=True)
    writefitsfile(w3up, 'output_files/w3up-hmi.fits', overwrite=True)
    writefitsfile(w5up, 'output_files/w5up-hmi.fits', overwrite=True)

    writefitsfile(w1lo, 'output_files/w1lo-hmi.fits', overwrite=True)
    writefitsfile(w3lo, 'output_files/w3lo-hmi.fits', overwrite=True)
    writefitsfile(w5lo, 'output_files/w5lo-hmi.fits', overwrite=True)

    writefitsfile(w1sig, 'output_files/w1err-hmi.fits', overwrite=True)
    writefitsfile(w1sig, 'output_files/w3err-hmi.fits', overwrite=True)
    writefitsfile(w1sig, 'output_files/w5err-hmi.fits', overwrite=True)
    #------------------------------------------------------------------------

    plt.plot(rmesh, w1)
    plt.fill_between(rmesh, w1-e1, w1+e1, alpha=0.5)
    plt.show()

