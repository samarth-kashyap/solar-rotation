import os
import numpy as np
from astropy.io import fits
from math import sqrt, pi
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.special import legendre as scipy_legendre

current_dir = os.path.dirname(os.path.realpath(__file__))
package_dir = os.path.dirname(current_dir)
input_dir = f"{package_dir}/input_files"
output_dir = f"{package_dir}/output_files"
data_dir = f"{package_dir}/input_files"
NAX = np.newaxis


def writefitsfile(a, fname, overwrite=False):
    from astropy.io import fits
    hdu = fits.PrimaryHDU()
    hdu.data = a
    hdu.writeto(fname, overwrite=overwrite)


def get_both_hemispheres(w2d_upper):
    """Mirrors the rotation data from top hemisphere to the bottom hemisphere

    Inputs:
    -------
    w2d_upper - np.ndarray(ndim=2, dtype=float)
        Rotation profile in the top hemisphere

    Returns:
    --------
    w2d - np.ndarray(ndim=2, dtype=float)
        Rotation profile on the full sphere
    """
    lenr = w2d_upper.shape[0]
    lent = w2d_upper.shape[1]

    w2d_lower = w2d_upper[:, ::-1][:, 1:]
    w2d = np.zeros((lenr, w2d_upper.shape[1]*2-1))

    w2d[:, :lent] = w2d_upper
    w2d[:, lent:] = w2d_lower
    return w2d

#---------------------------((( TEST FUNCTIONS --------------------
"""
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
"""
#------------------------------ TEST FUNCTIONS ))) --------------------


def get_omega024(rot_profile):
    omega0 = np.zeros(lenr)
    omega2 = np.zeros(lenr)
    omega4 = np.zeros(lenr)

    for i in range(lenr):
        omega0[i] = simps(rot_profile[i, :]*lp0, costh)*(2*0+1.)/2.
        omega2[i] = simps(rot_profile[i, :]*lp2, costh)*(2*2+1.)/2.
        omega4[i] = simps(rot_profile[i, :]*lp4, costh)*(2*4+1.)/2.

    return omega0, omega2, omega4


def get_w135(rot_profile, fprefix='rot'):
    omega0, omega2, omega4 = get_omega024(rot_profile)
    w1 = 2*sqrt(pi/3.)*rmesh*(omega0/1. - omega2/5.)
    w3 = 2*sqrt(pi/7.)*rmesh*(omega2/5. - omega4/9.)
    w5 = 2*sqrt(pi/11)*rmesh*(omega4/9.)
    writefitsfile(w1, f'{output_dir}/w1{fprefix}-hmi.fits', overwrite=True)
    writefitsfile(w3, f'{output_dir}/w3{fprefix}-hmi.fits', overwrite=True)
    writefitsfile(w5, f'{output_dir}/w5{fprefix}-hmi.fits', overwrite=True)
    return w1, w3, w5


def load_data():
    """Reads hemispherical rotation data and returns full rotation profiles"""

    # Reading radial-mesh, rotation profile and error
    rmesh = np.loadtxt(f'{input_dir}/rmesh.hmi')[::4]
    rot2d = np.loadtxt(f'{input_dir}/rot2d.hmi')
    err2d = np.loadtxt(f'{input_dir}/err2d.hmi')
    lenr = len(rmesh)

    # converting hemispherical theta-mesh to full spherical mesh
    tmesh = np.arange(rot2d.shape[1])
    theta = 90 - tmesh*15./8
    theta = np.append(-theta, theta[::-1][1:]) + 90.
    theta = theta*pi/180.

    rot2d = get_both_hemispheres(rot2d)
    err2d = get_both_hemispheres(err2d)
    err2d = err2d**2 #the linear operator acts on variance
    writefitsfile(rot2d, f'{output_dir}/rot2dfull-hmi.fits', overwrite=True)
    return (rmesh, theta), (rot2d, err2d)


def get_legpoly(theta):
    costh = np.cos(theta)
    lp0 = scipy_legendre(0)(costh)
    lp2 = scipy_legendre(2)(costh)
    lp4 = scipy_legendre(4)(costh)
    return lp0, lp2, lp4



if __name__=="__main__":
    (rmesh, theta), (rot2d, err2d) = load_data()
    lp0, lp2, lp4 = get_legpoly(theta)

    rot2dup = rot2d + err2d
    rot2dlo = rot2d - err2d
    #--------------------------- initializing the new profiles ----------------------
    err0 = np.zeros(lenr)
    err2 = np.zeros(lenr)
    err4 = np.zeros(lenr)
   #--------------------------------------------------------------------------------
    w1, w3, w5 = get_w135(rot2d, fprefix='rot')
    w1up, w3up, w5up = get_w135(rot2dup, fprefix='up')
    w1lo, w3lo, w5lo = get_w135(rot2dlo, fprefix='lo')

    for i in range(lenr):
        err0[i] = simps(err2d[i, :]*lp0**2 * sinth**2, theta)* ((2*0+1.)/2.)**2
        err2[i] = simps(err2d[i, :]*lp2**2 * sinth**2, theta)* ((2*2+1.)/2.)**2
        err4[i] = simps(err2d[i, :]*lp4**2 * sinth**2, theta)* ((2*4+1.)/2.)**2
    
    e1 = np.sqrt(2*sqrt(pi/3.))**2 *rmesh**2 *(err0/(1.**2) - err2/(5.**2))
    e3 = np.sqrt(2*sqrt(pi/7.))**2 *rmesh**2 *(err2/(5.**2) - err4/(9.**2))
    e5 = np.sqrt(2*sqrt(pi/11))**2 *rmesh**2 *(err4/(9.**2))

    w1sig = w1up - w1
    w3sig = w3up - w3
    w5sig = w5up - w5
    wsr_sig = [abs(w1sig), abs(w3sig), abs(w5sig)]

    #--------------------------- storing the profiles ----------------------
    writefitsfile(rmesh, f'{output_dir}/rad-hmi.fits', overwrite=True)
    writefitsfile(w1sig, 'output_files/w1err-hmi.fits', overwrite=True)
    writefitsfile(w1sig, 'output_files/w3err-hmi.fits', overwrite=True)
    writefitsfile(w1sig, 'output_files/w5err-hmi.fits', overwrite=True)
    #------------------------------------------------------------------------

    plt.plot(rmesh, w1)
    plt.fill_between(rmesh, w1-e1, w1+e1, alpha=0.5)
    plt.show()
