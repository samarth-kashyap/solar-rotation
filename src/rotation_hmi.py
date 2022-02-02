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
    print(f"Writing {fname}")


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


def get_omegak(rot_profile, smax):
    costh = np.cos(theta)
    klist = np.arange(smax+2)[::2]
    lplist = []
    omegak = []
    for k in klist:
        lplist.append(scipy_legendre(k)(costh))
        omegak.append(np.zeros(lenr))

    for ik, k in enumerate(klist):
        for ir in range(lenr):
            omegak[ik][ir] = simps(rot_profile[ir, :]*lplist[ik], costh)*(2*k+1.)/2.
    return omegak


def get_ws(rot_profile, smax, fprefix='rot'):
    omegak = get_omegak(rot_profile, smax)
    klist = np.arange(smax+2)[::2]
    slist = np.arange(smax)[::2] + 1
    ws = []
    for iess, s in enumerate(slist):
        prefac = -2*sqrt(pi/(2*s+1))
        omfac1 = 2*klist[iess] + 1
        omfac2 = omfac1 + 4
        ws.append(prefac*rmesh*(omegak[iess]/omfac1 - omegak[iess+1]/omfac2))
        writefitsfile(ws[iess], f'{output_dir}/w{s}{fprefix}-hmi.fits', overwrite=True)
    return ws


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

    tharr = np.arange(0, 90, 15)
    err_list = []
    for th in tharr:
        absdiff = abs(th*np.pi/180 - theta)
        thidx = np.argmin(abs(th*np.pi/180 - theta))
        if absdiff[thidx] == 0:
            err_list.append(err2d[:, thidx])
        else:
            err_plus = err2d[:, thidx+1]
            dth_plus = abs(theta[thidx+1] - th*np.pi/180)
            err_minus = err2d[:, thidx-1]
            dth_minus = abs(theta[thidx-1] - th*np.pi/180)
            dth = abs(theta[thidx+1] - theta[thidx-1])
            err_list.append((err_plus*dth_plus +
                             err_minus*dth_minus)/dth)

    rot2d = get_both_hemispheres(rot2d)
    err2d = get_both_hemispheres(err2d)
    err2d = err2d**2 #the linear operator acts on variance
    writefitsfile(rmesh, f'{output_dir}/rad-hmi.fits', overwrite=True)
    writefitsfile(rot2d, f'{output_dir}/rot2dfull-hmi.fits', overwrite=True)
    np.save(f"{output_dir}/err1d-hmi.npy", np.array(err_list))
    return (rmesh, theta), (rot2d, err2d)


if __name__=="__main__":
    smax = 5
    (rmesh, theta), (rot2d, err2d) = load_data()
    lenr = len(rmesh)

    rot2dup = rot2d + err2d
    rot2dlo = rot2d - err2d

    ws = get_ws(rot2d, smax, fprefix='rot')
    ws_up = get_ws(rot2dup, smax, fprefix='up')
    ws_lo = get_ws(rot2dlo, smax, fprefix='lo')
    wsig = []

    for i in range(len(ws)):
        sval = int(2*i+1)
        wsig.append(abs(ws_up[i] - ws[i]))
        writefitsfile(wsig[i], f'{output_dir}/w{sval}err-hmi.fits', overwrite=True)

    plt.plot(rmesh, ws[0])
    plt.fill_between(rmesh, ws[0]-wsig[0], ws[0]+wsig[0], alpha=0.5)
    plt.show()
