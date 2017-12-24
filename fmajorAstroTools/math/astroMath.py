import numpy as np
import ctypes
from ctypes import c_int, c_double, c_short
import warnings
import os
from fmajorAstroUtils.utils import messager
from scipy.interpolate import interp1d
import glob
import time
import unittest
import astropy.io.fits as pyfits
from scipy.interpolate import LSQUnivariateSpline
from scipy.interpolate import UnivariateSpline
from numpy.polynomial.legendre import Legendre


c_int_p = ctypes.POINTER(c_int)
c_short_p = ctypes.POINTER(c_short)
c_double_p = ctypes.POINTER(c_double)


class MathWarning(UserWarning):
    pass
#warnings.simplefilter("always", MathWarning)
warnings.simplefilter("ignore", MathWarning)

logger = messager.Messager(warningClass=MathWarning)

#region: how to use ctypes
#   so = CDLL(os.path.join(os.path.dirname(__file__), "src" ,"extract_row.so"))
#   extract_row_c = so.extract_row
#   np.ctypeslib.as_ctypes
#endregion: how to use ctypes

_lib = ctypes.CDLL(os.path.join(os.path.dirname(__file__), "lib", "astroMath.so"))
iterstat1d_c = _lib.iterstat1d
iterstat3d_c = _lib.iterstat3d
median_c = _lib.median
median_c.restype = c_double
median_c.argtypes = [c_double_p, c_int]
linearInterpWithError_c = _lib.linearInterpWithError
fillHollowness_c = _lib.fillHollowness
fillHollowness_c.restype = c_int

def fillHollowness(mask):
    """input bad mask, output filled bad mask"""
    mask = mask.copy()
    mask = mask.astype(np.int64)
    ny, nx = mask.shape
    mask_c = np.ctypeslib.as_ctypes(mask)
    fillHollowness_c(mask_c, ny, nx)
    return mask

# As comments for C version
# all mask is "badmask"
def iterstat1d(data, invvar=None, sigma=None, mask=None, sigrej=None, maxiter=10,
               useMedian=False, maxRej=None, stdRej=False, maxDev=None):
    """
        mask is bad mask
        the priority is invvar > sigma > mask
        if use invvar or sigma, will do weighted mean using invvar as the weight
            sigma is converted to invvar
        if have (input) sigma, use input sigma to do rejection
            else, calculate a sigma to do rejection (stdRej = 1)
        you can also force stdRej using stdRej=1
        badmask = newBadMask | mask
    """
    #region: setups
    if data.size==1:
        if invvar is not None:# use invvar to get badmask
            return {"mean": data, "std":np.nan, "median": data,'niter':0,
                    "badmask":invvar==0,
                    "newBadMask": np.zeros(data.shape, dtype=np.bool)}
        elif sigma is not None:# use sigma to get badmask
            return {"mean": data, "std":np.nan, "median": data,'niter':0,
                    "badmask": np.isfinite(sigma),
                    "newBadMask": np.zeros(data.shape, dtype=np.bool)}
        elif mask is not None:# use mask as badmask
            return {"mean": data, "std":np.nan, "median": data,'niter':0,
                    "badmask":mask, "newBadMask": np.zeros(data.shape, dtype=np.bool)}
        else:
            return {"mean": data, "std":np.nan, "median": data,'niter':0,
                    "badmask": np.zeros(data.shape, dtype=np.bool),
                    "newBadMask": np.zeros(data.shape, dtype=np.bool)}

    # data size larger than 1
    if sigrej is None:
        if sigma is None and invvar is None:# calculate std to do reject
            sigrej = 3
        else:# use self sigma for each point to do reject
            sigrej = 5
    if invvar is not None: # invvar is the 1st
        assert data.shape==invvar.shape
        assert np.all(np.isfinite(invvar))
        badmask = (invvar==0)
    elif sigma is not None: # then follow sigma
        invvar = 1./(sigma**2)
        assert data.shape==invvar.shape
        assert np.all(np.isfinite(invvar))
        badmask = (invvar==0)
    elif mask is not None: # or only have mask (and use stdRej)
        assert mask.shape == data.shape
        assert mask.dtype==np.bool, "mask dtype must be bool"
        badmask = mask.copy()
        stdRej = True
    else:# all points are good (and use stdRej)
        badmask = np.zeros(data.shape, dtype=np.bool)
        stdRej = True
    #endregion: setups

    # record initBadmask and initNanIndex
    initBadMask = badmask.copy()
    data = data.copy()
    initNanIndex = np.isnan(data) # record init nan position
    data[initNanIndex] = 0 # to remove NaN

    #!! first cycle
    #!! here data number must > 1
    allBad = False
    goodmask = np.logical_not(badmask)
    goodN = np.sum(goodmask)
    if goodN==0:
        mean = np.nan
        std = np.nan
        allBad = True
    else: # goodN > 0
        if invvar is not None:# weight mean
            mean = np.sum(data*invvar)/np.sum(invvar)
        else:# normal mean
            mean = np.sum(data*goodmask)/goodN
        if goodN>1:
            std = np.sqrt(np.sum(goodmask * (data-mean)**2) / (goodN - 1))
        else:
            std = np.nan
    if np.isfinite(std): # goodN>1
        if useMedian:
            median = np.median(data[goodmask])
            diff = np.abs(data - median)
        else:
            diff = np.abs(data - mean)
        if stdRej:# use std to reject
            badness = diff/std - sigrej
            badness = badness * (badness>0)
        else:# use self sigma to reject
            badness = diff*np.sqrt(invvar) - sigrej
            badness = badness * (badness>0)
        if maxDev is not None:# use maxDev to reject
            diff = (diff - maxDev)/maxDev
            diff = diff * (diff>0)
            badness = badness + diff
        newBadMask = badness>0

        newBadMask[badmask] = False# only leave new badmask
        newBadMaskNum = newBadMask.sum()
    else:
        newBadMaskNum = 0

    #print("first\n", mean,std, diff/std)
    niter = 0
    #print("\t python: iter: {} mean: {:.10f} std: {:.10f} badnumber: {} newBadNumber: {}".format(niter, mean, std, badmask.sum(), newBadMaskNum));
    #!! main loop
    while niter<maxiter and newBadMaskNum:
        niter += 1
        if maxRej is not None:# reject less than maxRej in each round
            newBadIndex = np.where(newBadMask)[0]
            newBadValue = badness[newBadMask]
            if maxRej==1:
                newBadValueIndexTodo = np.argmax(newBadValue)
            else:
                todo = newBadMaskNum if maxRej > newBadMaskNum else maxRej
                auxIndex = np.argsort(newBadValue)
                newBadValueIndexTodo = auxIndex[-todo:]

            newBadIndex = newBadIndex[newBadValueIndexTodo]
            #print("\t\t{}:{}".format(newBadIndex, badness[newBadIndex]))
            badmask[newBadIndex] = True
        else: # reject all
            #print("\t\t{}:{}".format(newBadMask, badness[newBadIndex]))
            badmask[newBadMask] = True
        #!! after select bad mask
        #!! calculate mean and std
        goodmask = np.logical_not(badmask)
        goodN = np.sum(goodmask)
        if goodN==0:
            mean = np.nan
            std = np.nan
            allBad = True
            break
        else: # goodN > 0
            if invvar is not None:
                mean = np.sum(data*invvar)/np.sum(invvar)
            else:
                mean = np.sum(data*goodmask)/goodN
            if goodN>1:
                std = np.sqrt(np.sum(goodmask * (data-mean)**2) / (goodN - 1))
            else:
                std = np.nan
                break
        if useMedian:
            median = np.median(data[goodmask])
            diff = np.abs(data - median)
        else:
            diff = np.abs(data - mean)
        if stdRej:
            badness = diff/std - sigrej
            badness = badness * (badness>0)
        else:
            badness = diff*np.sqrt(invvar) - sigrej
            badness = badness * (badness>0)
        if maxDev is not None:
            diff = (diff - maxDev)/maxDev
            diff = diff * (diff>0)
            badness = badness + diff
        newBadMask = badness>0
        #!! calculate new badMask
        newBadMask[badmask] = False
        newBadMaskNum = newBadMask.sum()

        #print("\t python: iter: {} mean: {:.10f} std: {:.10f} badnumber: {} newBadNumber: {}".format(niter, mean, std, badmask.sum(), newBadMaskNum));

    newBadMask = (badmask.astype(np.int) - initBadMask.astype(np.int)).astype(np.bool)
    result = {"mean": mean, "std":std, "badmask": badmask,
              "niter":niter, "newBadMask": newBadMask}
    if niter==maxiter:
        logger.warn("iterstat reach max iteration number")
        goodmask = np.logical_not(badmask)

    if allBad:
        result["median"] = np.nan
    else:
        result["median"] = np.median(data[goodmask])

    return result
# C version of last function
def iterstat1dC(data, invvar=None, sigma=None, mask=None, sigrej=None, hsigrej=None, lsigrej=None,
                maxiter=10, useMedian=False, maxRej=None, stdRej=False, maxDev=None):
    """
        mask is bad mask
        the priority is invvar > sigma > mask
        if use invvar or sigma, will do weighted mean using invvar as the weight
            sigma is converted to invvar
        if have (input) sigma, use input sigma to do rejection
            else, calculate a sigma to do rejection (stdRej = 1)
        you can also force stdRej using stdRej=1
        badmask = newBadMask | mask
    """
    #region: setups
    if data.size==1:
        if invvar is not None:# use invvar to get badmask
            return {"mean": data, "std":np.nan, "median": data,'niter':0,
                    "badmask":invvar==0,
                    "newBadMask": np.zeros(data.shape, dtype=np.bool)}
        elif sigma is not None:# use sigma to get badmask
            return {"mean": data, "std":np.nan, "median": data,'niter':0,
                    "badmask": np.isfinite(sigma),
                    "newBadMask": np.zeros(data.shape, dtype=np.bool)}
        elif mask is not None:# use mask as badmask
            return {"mean": data, "std":np.nan, "median": data,'niter':0,
                    "badmask":mask, "newBadMask": np.zeros(data.shape, dtype=np.bool)}
        else:
            return {"mean": data, "std":np.nan, "median": data,'niter':0,
                    "badmask": np.zeros(data.shape, dtype=np.bool),
                    "newBadMask": np.zeros(data.shape, dtype=np.bool)}

    # data size larger than 1
    if sigrej is None:
        if sigma is None and invvar is None:
            sigrej = 3
        else:
            sigrej = 5
    if hsigrej is None:
        hsigrej = sigrej
    if lsigrej is None:
        lsigrej = sigrej

    if invvar is not None:
        assert data.shape==invvar.shape
        assert np.all(np.isfinite(invvar))
        badmask = (invvar==0)
    elif sigma is not None:
        invvar = 1./(sigma**2)
        assert data.shape==invvar.shape
        assert np.all(np.isfinite(invvar))
        badmask = (invvar==0)
    elif mask is not None:
        assert mask.shape == data.shape
        assert mask.dtype==np.bool, "mask dtype must be bool"
        invvar = np.ones(data.shape, dtype=np.float64) * -1
        badmask = mask.copy()
        stdRej = True
    else:
        invvar = np.ones(data.shape, dtype=np.float64) * -1
        badmask = np.zeros(data.shape, dtype=np.bool)
        stdRej = True

    if maxRej is None:
        maxRej = -1
    if maxDev is None:
        maxDev = -1

    initBadMask = badmask.copy()
    data = data.copy()
    initNanIndex = np.isnan(data) # record init nan position
    data[initNanIndex] = 0 # to remove NaN
    #endregion: setups
    #region: data 2 ctypes
    data = data.astype(np.float64, copy=True, order="C")
    invvar = invvar.astype(np.float64, order="C")
    initBadMask = initBadMask.astype(np.short, order="C")
    badmask = badmask.astype(np.short, order="C")
    newBadMask = np.zeros(data.shape, dtype=np.short, order="C")
    result = np.zeros(3, dtype=np.float64, order="C")

    data_c = data.ctypes.data_as(c_double_p)
    invvar_c = invvar.ctypes.data_as(c_double_p)
    initBadMask_c = initBadMask.ctypes.data_as(c_short_p)
    badmask_c = badmask.ctypes.data_as(c_short_p)
    newBadMask_c = newBadMask.ctypes.data_as(c_short_p)
    result_c = result.ctypes.data_as(c_double_p)

    #endregion: data 2 ctypes
    #int iterstat1d(const double data[], int M,
    #           const double invvar[], double lsigrej, double hsigrej,
    #           int maxiter, short const initmask[], short badmask[], short newBadMask[],
    #           short calSigma, short useMedian, double maxdev,
    #           double result[]) {
    iterstat1d_c.argtypes = [c_double_p, c_int,
                             c_double_p, c_double, c_double,
                             c_int, c_short_p, c_short_p, c_short_p,
                             c_short, c_short, c_int, c_double,
                             c_double_p]
    # main c function
    niter = iterstat1d_c(data_c, data.size,
                         invvar_c, lsigrej, hsigrej,
                         maxiter, initBadMask_c, badmask_c, newBadMask_c,
                         stdRej, useMedian, maxRej, maxDev,
                         result_c)
    mean = result[0]
    std  = result[1]
    theMedian = result[2]

    result = {"mean": mean, "std":std, "badmask": badmask.astype(np.bool),
              "niter":niter, "newBadMask": newBadMask.astype(np.bool)}
    if niter==maxiter:
        logger.warn("iterstat reach max iteration number")
    result["median"] = theMedian

    return result

def iterstat3d(data, invvar=None, sigma=None, mask=None, sigrej=None, maxiter=10,
               useMedian=False, maxRej=None, stdRej=False, maxDev=None):
    """the priority is invvar > sigma > mask
    """
    #region: setups
    assert data.ndim==3

    # data size larger than 1
    if sigrej is None:
        if sigma is None and invvar is None:
            sigrej = 3
        else:
            sigrej = 5

    if invvar is not None:
        assert data.shape==invvar.shape
        assert np.all(np.isfinite(invvar))
    elif sigma is not None:
        invvar = 1./(sigma**2)
        assert data.shape==invvar.shape
        assert np.all(np.isfinite(invvar))
    elif mask is not None:
        assert mask.shape == data.shape
        assert mask.dtype==np.bool, "mask dtype must be bool"
        badmask = mask.copy()
        stdRej = True
    else:
        badmask = np.zeros(data.shape, dtype=np.bool)
        stdRej = True
    #endregion: setups
    #region: data 2 ctypes
    data = data.copy()
    initNanIndex = np.isnan(data) # record init nan position
    data[initNanIndex] = 0 # to remove NaN
    a, b, c = data.shape # N, ny, nx

    #!! ./test/testArrayReshape.py
    #!! show some numpy reshape exam that can help to understand thing below
    """
        a = np.arange(120).reshape((4,5,6))
        b = a.T

        z, m, n = a.shape
        for eachM in range(m):
            for eachN in range(n):
                print(a[:, eachM, eachN])
                assert np.all(a[:, eachM, eachN]==b[eachN, eachM, :])
    """
    data = data.T
    if invvar is not None:
        invvar = invvar.T
    else:
        badmask = badmask.T

    result_mean = np.empty((c,b))
    result_std  = np.empty((c,b))
    result_badmask  = np.empty((c,b,a), dtype=np.bool)
    result_newBadMask  = np.empty((c,b,a), dtype=np.bool)

    result_median = np.empty((c,b))
    #endregion: data 2 ctypes

    for eachC in range(c):
        for eachB in range(b):
            #print("{}/{} {}/{}".format(eachC, c, eachB, b))
            #print(data[eachC, eachB])
            if invvar is not None:
                aux = iterstat1dC(data[eachC, eachB], invvar=invvar[eachC, eachB], sigrej=sigrej, maxiter=maxiter, useMedian=useMedian, maxRej=maxRej, stdRej=stdRej, maxDev=maxDev)
            else:
                aux = iterstat1dC(data[eachC, eachB], mask=badmask[eachC, eachB], sigrej=sigrej, maxiter=maxiter, useMedian=useMedian, maxRej=maxRej, stdRej=stdRej, maxDev=maxDev)
            #!! record result
            result_mean[eachC, eachB] = aux["mean"]
            result_std[eachC, eachB] = aux["std"]
            result_badmask[eachC, eachB, :] = aux["badmask"][:]
            result_newBadMask[eachC, eachB, :] = aux["newBadMask"][:]
            result_median[eachC, eachB] = aux["median"]

    result = {"mean": result_mean.T, "std":result_std.T,
              "badmask": result_badmask.T,
              "newBadMask": result_newBadMask.T}
    result["median"] = result_median.T

    return result

def iterstat3dC(data, invvar=None, sigma=None, mask=None, sigrej=None, hsigrej=None, lsigrej=None,
                maxiter=10, useMedian=False, maxRej=None, stdRej=False, maxDev = None):
    """the priority is invvar > sigma > mask
    """
    #region: setups
    assert data.ndim==3

    # data size larger than 1
    if sigrej is None:
        if sigma is None and invvar is None:
            sigrej = 3
        else:
            sigrej = 5
    if hsigrej is None:
        hsigrej = sigrej
    if lsigrej is None:
        lsigrej = sigrej

    if invvar is not None:
        assert data.shape==invvar.shape
        assert np.all(np.isfinite(invvar))
    elif sigma is not None:
        invvar = 1./(sigma**2)
        assert data.shape==invvar.shape
        assert np.all(np.isfinite(invvar))
    elif mask is not None:
        assert mask.shape == data.shape
        assert mask.dtype==np.bool, "mask dtype must be bool"
        badmask = mask.copy()
        stdRej = True
    else:
        badmask = np.zeros(data.shape, dtype=np.bool)
        stdRej = True
    if maxRej is None:
        maxRej = -1
    if maxDev is None:
        maxDev = -1

    #endregion: setups

    #region: data 2 ctypes
    data = data.copy()
    initNanIndex = np.isnan(data) # record init nan position
    data[initNanIndex] = 0 # to remove NaN
    a, b, c = data.shape

    #!! ./test/testArrayReshape.py
    #!! show some numpy reshape exam that can help to understand thing below
    data = data.T
    data = data.astype(np.float64, order="C")
    if invvar is not None:
        invvar = invvar.T
        invvar = invvar.astype(np.float64, order="C")
        badmask = (invvar==0)
        badmask = badmask.astype(np.short, order="C")
        initBadMask = (invvar==0)
        initBadMask = initBadMask.astype(np.short, order="C")
    else:
        badmask = badmask.T
        badmask = badmask.astype(np.short, order="C")
        invvar = np.zeros(badmask.shape, dtype=np.float64, order="C")
        invvar[:,:,0] = -1
        initBadMask = badmask.copy()

    newBadMask = np.zeros(data.shape, dtype=np.short, order="C")
    cresult = np.zeros((c,b,3), dtype=np.float64, order="C")

    # c,b,a
    data_c = np.ctypeslib.as_ctypes(data)
    invvar_c = np.ctypeslib.as_ctypes(invvar)
    initBadMask_c = np.ctypeslib.as_ctypes(initBadMask)
    badmask_c = np.ctypeslib.as_ctypes(badmask)
    newBadMask_c = np.ctypeslib.as_ctypes(newBadMask)
    # c,b,3
    result_c = np.ctypeslib.as_ctypes(cresult)

    doubleClass = data_c.__class__
    resultClass = result_c.__class__
    shortClass = badmask_c.__class__

    #endregion: data 2 ctypes

    #void iterstat3d(const double data[], int c, int b, int M,
    #           const double invvar[], double lsigrej, double hsigrej,
    #           int maxiter, short const initmask[], short badmask[], short newBadMask[],
    #           short stdRej, short useMedian, int maxRej,
    #           double result[]) {
    iterstat3d_c.argtypes = [doubleClass, c_int, c_int, c_int,
                             doubleClass, c_double, c_double,
                             c_int, shortClass, shortClass, shortClass,
                             c_short, c_short, c_int, c_double,
                             resultClass]
    #print("data:", data.shape)
    #print("initBadMask:", initBadMask.shape)
    #print("badmask:", badmask.shape)
    #print("newBadMask:", newBadMask.shape)
    iterstat3d_c(data_c, c, b, a,
                 invvar_c, lsigrej, hsigrej,
                 maxiter, initBadMask_c, badmask_c, newBadMask_c,
                 stdRej, useMedian, maxRej, maxDev,
                 result_c)

    cresult = cresult.T
    badmask = badmask.T.astype(np.bool)
    newBadMask = newBadMask.T.astype(np.bool)

    result = {"mean": cresult[0], "std":cresult[1],
              "badmask": badmask, "newBadMask": newBadMask}
    result["median"] = cresult[2]

    return result

def linearInterpWithError(x, y, sigma, xx):
    """sigmasigma[i] = sqrt(hi * sigma[ii] * sigma[ii] + (1-hi) * sigma[ii+1] * sigma[ii+1])"""
    assert np.all(np.diff(x)>0)
    assert np.all(np.diff(xx)>0)
    assert x.ndim==1
    M = xx.size
    x = x.copy()
    xx = xx.copy()
    yy = np.zeros(M)
    y = y.copy()
    sigma = sigma.copy()
    sigmasigma = np.zeros(M)

    goodPoint = np.isfinite(sigma)
    y = y[goodPoint]
    x = x[goodPoint]
    sigma = sigma[goodPoint]
    N = x.size

    x = x.astype(np.float64, order="C")
    xx = xx.astype(np.float64, order="C")
    y = y.astype(np.float64, order="C")
    yy = yy.astype(np.float64, order="C")
    sigma = sigma.astype(np.float64, order="C")
    sigmasigma = sigmasigma.astype(np.float64, order="C")

    x_c = np.ctypeslib.as_ctypes(x)
    y_c = np.ctypeslib.as_ctypes(y)
    sigma_c = np.ctypeslib.as_ctypes(sigma)
    xx_c = np.ctypeslib.as_ctypes(xx)
    yy_c = np.ctypeslib.as_ctypes(yy)
    sigmasigma_c = np.ctypeslib.as_ctypes(sigmasigma)

    #void linearInterpWithError(double x[], double y[], double sigma[], int N,
    #                     double xx[], double yy[], double sigmasigma[], int M)
    linearInterpWithError_c.argtypes = [c_double_p, c_double_p, c_double_p, c_int,
                                  c_double_p, c_double_p, c_double_p, c_int]
    linearInterpWithError_c(x_c, y_c, sigma_c, N,
                      xx_c, yy_c, sigmasigma_c, M)
    return {"x": xx, "y": yy, "sigma": sigmasigma}

def gaussian(x, A, mu, sigma):
    return A*np.exp(-(x-mu)**2/(2*sigma**2))/(sigma*np.sqrt(2*np.pi))

def Lorentzian(x, A, mu, w):
    return A/(1+((x-mu)/(w/2))**2)

def medianC(A):
    N = A.size
    data = A.astype(np.float64)
    data_c = np.ctypeslib.as_ctypes(data)
    return median_c(data_c, N)

def optimalFixExtract1D(x, y, var, profile):
    assert x.shape==y.shape==var.shape==profile.shape
    #!! boxsum
    boxsum = y.sum()
    boxsumvar = var.sum()
    #!! first calculation
    mask = np.isfinite(var)
    goodN = mask.sum()
    tFlux = np.sum(profile)
    denom = np.sum((profile[mask]**2/var[mask]))
    C = np.sum(y[mask]*profile[mask]/var[mask])/denom
    varC = 1/denom
    chis = (C*profile[mask]-y[mask])**2/var[mask]
    chiss = np.sqrt(chis)
    chi2 = np.sum(chis)
    reduce_chi2 = chi2/(goodN-1)
    flux = C*tFlux
    fluxvar = varC*tFlux**2
    #!! detect bad pixel
    result = iterstat1dC(chiss, mask=np.logical_not(mask), sigrej=3, useMedian=1)
    chissig = result['std']
    newmask = np.logical_not(result['badmask'])
    newbadmask = result['newBadMask']
    badpixelFlag = False
    print(("".join(["{:9.3e},"]*chiss.size)).format(*(chiss.tolist())))
    print(("".join(["{:^10d}"]*newmask.size)).format(*(np.logical_not(newmask).tolist())))
    if not np.all(mask==newmask):
        mask = newmask
        goodN = mask.sum()
        denom = np.sum((profile[mask]**2/var[mask]))
        C = np.sum(y[mask]*profile[mask]/var[mask])/denom
        varC = 1/denom
        chis = (C*profile[mask]-y[mask])**2/var[mask]
        chi2 = np.sum(chis)
        reduce_chi2 = chi2/(goodN-1)
        flux = C*tFlux
        fluxvar = varC*tFlux**2
        badpixelFlag = True


    return {"C":C, "varC":varC, "outy":profile*C,
            "flux":flux, "fluxvar":fluxvar,
            "boxsum":boxsum, "boxsumvar":boxsumvar,
            "chiss":chiss, 'chissig':chissig,
            "chi2":chi2, "reduce_chi2":reduce_chi2,
            "newmask":newmask, 'newBadMask':newbadmask,
            "badpixelFlag":badpixelFlag
            }
def optimalFixExtract2D(xx, yy, varvar, profile):
    ny, nx = xx.shape
    assert xx.shape==yy.shape==varvar.shape
    assert profile.size == nx
    result_C = np.zeros(ny)
    result_varC = np.ones(ny) * np.inf
    result_outy = np.zeros((ny, nx))
    result_flux = np.zeros(ny) * np.nan
    result_fluxvar = np.zeros(ny)
    result_boxsum = np.zeros(ny) * np.nan
    result_boxsumvar = np.zeros(ny) * np.inf
    result_chiss = np.zeros((ny, nx))
    result_chissig = np.zeros(ny)
    result_chi2 = np.zeros(ny)
    result_reduce_chi2 = np.zeros(ny)
    result_newmask = np.zeros((ny, nx))
    result_newBadMask = np.zeros((ny, nx))
    result_badpixelFlag = np.zeros((ny, nx))

    for yindex in np.arange(ny):
        x = xx[yindex]
        y = yy[yindex]
        var = varvar[yindex]
        #!! boxsum
        boxsum = y.sum()
        boxsumvar = var.sum()
        #!! first calculation
        mask = np.isfinite(var)
        goodN = mask.sum()
        tFlux = np.sum(profile)
        denom = np.sum((profile[mask]**2/var[mask]))
        C = np.sum(y[mask]*profile[mask]/var[mask])/denom
        varC = 1/denom
        chis = (C*profile-y)**2/var
        chiss = np.sqrt(chis)
        chi2 = np.sum(chis)
        reduce_chi2 = chi2/(goodN-1)
        flux = C*tFlux
        fluxvar = varC*tFlux**2
        #!! detect bad pixel
        result = iterstat1dC(chiss, mask=np.logical_not(mask), sigrej=3, useMedian=1)
        chissig = result['std']
        newmask = np.logical_not(result['badmask'])
        newbadmask = result['newBadMask']
        badpixelFlag = False
        #print(("".join(["{:9.3e},"]*chiss.size)).format(*(chiss.tolist())))
        #print(("".join(["{:^10d}"]*newmask.size)).format(*(np.logical_not(newmask).tolist())))
        if not np.all(mask==newmask):
            mask = newmask
            goodN = mask.sum()
            denom = np.sum((profile[mask]**2/var[mask]))
            C = np.sum(y[mask]*profile[mask]/var[mask])/denom
            varC = 1/denom
            chis = (C*profile[mask]-y[mask])**2/var[mask]
            chi2 = np.sum(chis)
            reduce_chi2 = chi2/(goodN-1)
            flux = C*tFlux
            fluxvar = varC*tFlux**2
            badpixelFlag = True
        result_C[yindex]            = C
        result_varC[yindex]         = varC
        result_outy[yindex,:]       = (profile*C)[:]
        result_flux[yindex]         = flux
        result_fluxvar[yindex]      = fluxvar
        result_boxsum[yindex]       = boxsum
        result_boxsumvar[yindex]    = boxsumvar
        result_chiss[yindex,:]      = chiss[:]
        result_chissig[yindex]      = chissig
        result_chi2[yindex]         = chi2
        result_reduce_chi2[yindex]  = reduce_chi2
        result_newmask[yindex,:]    = newmask[:]
        result_newBadMask[yindex,:] = newbadmask[:]
        result_badpixelFlag[yindex] = badpixelFlag


    return {"C":result_C, "varC":result_varC, "outy":result_outy,
            "flux":result_flux, "fluxvar":result_fluxvar,
            "boxsum":result_boxsum, "boxsumvar":result_boxsumvar,
            "chiss":result_chiss, 'chissig':result_chissig,
            "chi2":result_chi2, "reduce_chi2":result_reduce_chi2,
            "newmask":result_newmask, 'newBadMask':result_newBadMask,
            "badpixelFlag":result_badpixelFlag
            }

class TestMath(unittest.TestCase):
    #@unittest.skip("skip")
    def setUp(self):
        pass
    def tearDown(self):
        pass
    @unittest.skip("skip")
    def test_median(self):
        N = 999999
        print()
        for i in range(N):
            M = np.random.randint(1,9999)
            data = np.random.uniform(0,1,M)
            #print("{}/{} size:{}".format(i,N,M))
            self.assertEqual(np.median(data), medianC(data))
    def floatEqual(self, a, b):
        if np.isnan(a) and np.isnan(b):
            return True
        tol = 1e-10
        return np.abs(a-b)<tol
    def iterstat1d_test(self, doPlot=False, doPrint=False):
        #region: set data and figures
        self.count+=1
        secondwave=1
        sys.stdout.write("{}\r".format(self.count))
        sys.stdout.flush()
        if doPlot:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            c = cm.get_cmap("rainbow")
            c = c(np.linspace(0,1,9))
        N = np.random.randint(1,200)
        # test for iterstat1d
        data = np.random.randn(N)+10
        sigma = (np.random.random(N)-0.5)/2+1
        badmask = np.zeros(data.shape, dtype=np.bool)
        badPointMask = np.zeros(data.shape, dtype=np.bool)

        if data.size>1:
            nbad = np.random.randint(0, np.min([data.size, 10]))
            # add bad points
            badIndex = np.random.randint(0, high=data.size-1, size=nbad)
            for eachIndex in badIndex:
                sigma[eachIndex] = np.inf
                badmask[eachIndex] = 1
                badPointMask[eachIndex] = 1

            # add singular points
            nbad = np.random.randint(0, np.min([data.size, 20]))
            badIndex = np.random.randint(0, high=data.size-1, size=nbad)
            for eachIndex in badIndex:
                data[eachIndex] += (np.random.random()-0.5)*20

            # add nan points
            nbad = np.random.randint(0, np.min([data.size, 5]))
            badIndex = np.random.randint(0, high=data.size-1, size=nbad)
            for eachIndex in badIndex:
                data[eachIndex] = np.nan
                sigma[eachIndex] = np.inf
                badmask[eachIndex] = 1

        if doPlot:
            plt.scatter(range(data.size), data)
            plt.errorbar(range(data.size), data, ls="None", yerr=sigma)
            plt.scatter(np.where(badPointMask)[0], data[badPointMask], color="r")
        #endregion: set data and figures
        #region: tests
        if doPrint:
            print()
        r0 = iterstat1d(data,  maxRej=4)
        if doPrint:
            print("r0:   {} +- {} @ {} {}".format(r0["mean"], r0["std"], r0["niter"], r0["median"]))
        rt0 = iterstat1dC(data,  maxRej=4)
        if doPrint:
            print("rt0:  {} +- {} @ {} {}".format(rt0["mean"], rt0["std"], rt0["niter"], rt0["median"]))
        self.assertTrue(self.floatEqual(r0["mean"],   rt0["mean"]))
        self.assertTrue(self.floatEqual(r0["std"],    rt0["std"]))
        self.assertTrue(self.floatEqual(r0["niter"],  rt0["niter"]))
        self.assertTrue(self.floatEqual(r0["median"], rt0["median"]))
        self.assertTrue(np.all(r0["badmask"]==rt0["badmask"]))
        self.assertTrue(np.all(r0["newBadMask"]==rt0["newBadMask"]))

        r1 = iterstat1d(data, mask=badmask,  maxRej=3)
        if doPrint:
            print("r1:   {} +- {} @ {} {}".format(r1["mean"], r1["std"], r1["niter"], r1["median"]))
        rt1 = iterstat1dC(data, mask=badmask,  maxRej=3)
        if doPrint:
            print("rt1:  {} +- {} @ {} {}".format(rt1["mean"], rt1["std"], rt1["niter"], rt1["median"]))
        self.assertTrue(self.floatEqual(r1["mean"],   rt1["mean"]))
        self.assertTrue(self.floatEqual(r1["std"],    rt1["std"]))
        self.assertTrue(self.floatEqual(r1["niter"],  rt1["niter"]))
        self.assertTrue(self.floatEqual(r1["median"], rt1["median"]))
        self.assertTrue(np.all(r1["badmask"]==rt1["badmask"]))
        self.assertTrue(np.all(r1["newBadMask"]==rt1["newBadMask"]))

        r2 = iterstat1d(data, sigma=sigma,  maxRej=2)
        if doPrint:
            print("r2:   {} +- {} @ {} {}".format(r2["mean"], r2["std"], r2["niter"], r2["median"]))
        rt2 = iterstat1dC(data, sigma=sigma,  maxRej=2)
        if doPrint:
            print("rt2:  {} +- {} @ {} {}".format(rt2["mean"], rt2["std"], rt2["niter"], rt2["median"]))
        self.assertTrue(np.all(r2["badmask"]==rt2["badmask"]))
        self.assertTrue(np.all(r2["newBadMask"]==rt2["newBadMask"]))
        self.assertTrue(self.floatEqual(r2["mean"],   rt2["mean"]))
        self.assertTrue(self.floatEqual(r2["std"],    rt2["std"]))
        self.assertTrue(self.floatEqual(r2["niter"],  rt2["niter"]))
        self.assertTrue(self.floatEqual(r2["median"], rt2["median"]))

        r3 = iterstat1d(data, invvar=1/sigma**2,  maxRej=3, maxDev=2)
        if doPrint:
            print("r3:   {} +- {} @ {} {}".format(r3["mean"], r3["std"], r3["niter"], r3["median"]))
        rt3 = iterstat1dC(data, invvar=1/sigma**2,  maxRej=3, maxDev=2)
        if doPrint:
            print("rt3:  {} +- {} @ {} {}".format(rt3["mean"], rt3["std"], rt3["niter"], rt3["median"]))
        self.assertTrue(np.all(r3["badmask"]==rt3["badmask"]))
        self.assertTrue(np.all(r3["newBadMask"]==rt3["newBadMask"]))
        self.assertTrue(self.floatEqual(r3["mean"],   rt3["mean"]))
        self.assertTrue(self.floatEqual(r3["std"],    rt3["std"]))
        self.assertTrue(self.floatEqual(r3["niter"],  rt3["niter"]))
        self.assertTrue(self.floatEqual(r3["median"], rt3["median"]))

        r4 = iterstat1d(data, sigma=sigma, invvar=1/sigma**2,  maxRej=4)
        if doPrint:
            print("r4:   {} +- {} @ {} {}".format(r4["mean"], r4["std"], r4["niter"], r4["median"]))
        rt4 = iterstat1dC(data, sigma=sigma, invvar=1/sigma**2,  maxRej=4)
        if doPrint:
            print("rt4:  {} +- {} @ {} {}".format(rt4["mean"], rt4["std"], rt4["niter"], rt4["median"]))
        self.assertTrue(np.all(r4["badmask"]==rt4["badmask"]))
        self.assertTrue(np.all(r4["newBadMask"]==rt4["newBadMask"]))
        self.assertTrue(self.floatEqual(r4["mean"],   rt4["mean"]))
        self.assertTrue(self.floatEqual(r4["std"],    rt4["std"]))
        self.assertTrue(self.floatEqual(r4["niter"],  rt4["niter"]))
        self.assertTrue(self.floatEqual(r4["median"], rt4["median"]))

        r5 = iterstat1d(data, invvar=1/sigma**2, maxiter=2,  maxRej=5)
        if doPrint:
            print("r5:   {} +- {} @ {} {}".format(r5["mean"], r5["std"], r5["niter"], r5["median"]))
        rt5 = iterstat1dC(data, invvar=1/sigma**2, maxiter=2,  maxRej=5)
        if doPrint:
            print("rt5:  {} +- {} @ {} {}".format(rt5["mean"], rt5["std"], rt5["niter"], rt5["median"]))
        self.assertTrue(np.all(r5["badmask"]==rt5["badmask"]))
        self.assertTrue(np.all(r5["newBadMask"]==rt5["newBadMask"]))
        self.assertTrue(self.floatEqual(r5["mean"],   rt5["mean"]))
        self.assertTrue(self.floatEqual(r5["std"],    rt5["std"]))
        self.assertTrue(self.floatEqual(r5["niter"],  rt5["niter"]))
        self.assertTrue(self.floatEqual(r5["median"], rt5["median"]))

        r6 = iterstat1d(data, invvar=1/sigma**2, maxiter=2,  stdRej=True, maxRej=1)
        if doPrint:
            print("r6:   {} +- {} @ {} {}".format(r6["mean"], r6["std"], r6["niter"], r6["median"]))
        rt6 = iterstat1dC(data, invvar=1/sigma**2, maxiter=2,  stdRej=True, maxRej=1)
        if doPrint:
            print("rt6:  {} +- {} @ {} {}".format(rt6["mean"], rt6["std"], rt6["niter"], rt6["median"]))
        self.assertTrue(np.all(r6["badmask"]==rt6["badmask"]))
        self.assertTrue(np.all(r6["newBadMask"]==rt6["newBadMask"]))
        self.assertTrue(self.floatEqual(r6["mean"],   rt6["mean"]))
        self.assertTrue(self.floatEqual(r6["std"],    rt6["std"]))
        self.assertTrue(self.floatEqual(r6["niter"],  rt6["niter"]))
        self.assertTrue(self.floatEqual(r6["median"], rt6["median"]))

        r8 = iterstat1d(data,  maxRej=1, maxiter=20)
        if doPrint:
            print("r8:   {} +- {} @ {} {}".format(r8["mean"], r8["std"], r8["niter"], r8["median"]))
        rt8 = iterstat1dC(data,  maxRej=1, maxiter=20)
        if doPrint:
            print("rt8:  {} +- {} @ {} {}".format(rt8["mean"], rt8["std"], rt8["niter"], rt8["median"]))
        self.assertTrue(np.all(r8["badmask"]==rt8["badmask"]))
        self.assertTrue(np.all(r8["newBadMask"]==rt8["newBadMask"]))
        self.assertTrue(self.floatEqual(r8["mean"],   rt8["mean"]))
        self.assertTrue(self.floatEqual(r8["std"],    rt8["std"]))
        self.assertTrue(self.floatEqual(r8["niter"],  rt8["niter"]))
        self.assertTrue(self.floatEqual(r8["median"], rt8["median"]))

        r7 = iterstat1d(data, invvar=1/sigma**2, maxRej=1, maxiter=20)
        if doPrint:
            print("r7:   {} +- {} @ {} {}".format(r7["mean"], r7["std"], r7["niter"], r7["median"]))
        rt7 = iterstat1dC(data, invvar=1/sigma**2, maxRej=1, maxiter=20)
        if doPrint:
            print("rt7:  {} +- {} @ {} {}".format(rt7["mean"], rt7["std"], rt7["niter"], rt7["median"]))
        self.assertTrue(np.all(r7["badmask"]==rt7["badmask"]))
        self.assertTrue(np.all(r7["newBadMask"]==rt7["newBadMask"]))
        self.assertTrue(self.floatEqual(r7["mean"],   rt7["mean"]))
        self.assertTrue(self.floatEqual(r7["std"],    rt7["std"]))
        self.assertTrue(self.floatEqual(r7["niter"],  rt7["niter"]))
        self.assertTrue(self.floatEqual(r7["median"], rt7["median"]))

        if doPlot:
            plt.plot([0,data.size-1], [r1["mean"], r1["mean"]], label="r1", c=c[1])
            plt.plot([0,data.size-1], [r1["mean"]+r1["std"], r1["mean"]+r1["std"]], "--", label="r1", c=c[1])
            plt.plot([0,data.size-1], [r1["mean"]-r1["std"], r1["mean"]-r1["std"]], "--", label="r1", c=c[1])
            plt.plot([0,data.size-1], [r2["mean"], r2["mean"]], label="r2", c=c[2])
            plt.plot([0,data.size-1], [r2["mean"]+r2["std"], r2["mean"]+r2["std"]], "--", label="r2", c=c[2])
            plt.plot([0,data.size-1], [r2["mean"]-r2["std"], r2["mean"]-r2["std"]], "--", label="r2", c=c[2])
            plt.plot([0,data.size-1], [r3["mean"], r3["mean"]], label="r3", c=c[3])
            plt.plot([0,data.size-1], [r3["mean"]+r3["std"], r3["mean"]+r3["std"]], "--", label="r3", c=c[3])
            plt.plot([0,data.size-1], [r3["mean"]-r3["std"], r3["mean"]-r3["std"]], "--", label="r3", c=c[3])
            plt.plot([0,data.size-1], [r4["mean"], r4["mean"]], label="r4", c=c[4])
            plt.plot([0,data.size-1], [r4["mean"]+r4["std"], r4["mean"]+r4["std"]], "--", label="r4", c=c[4])
            plt.plot([0,data.size-1], [r4["mean"]-r4["std"], r4["mean"]-r4["std"]], "--", label="r4", c=c[4])
            plt.plot([0,data.size-1], [r5["mean"], r5["mean"]], label="r5", c=c[5])
            plt.plot([0,data.size-1], [r5["mean"]+r5["std"], r5["mean"]+r5["std"]], "--", label="r5", c=c[5])
            plt.plot([0,data.size-1], [r5["mean"]-r5["std"], r5["mean"]-r5["std"]], "--", label="r5", c=c[5])
            plt.plot([0,data.size-1], [r6["mean"], r6["mean"]], label="r6", c=c[6])
            plt.plot([0,data.size-1], [r6["mean"]+r6["std"], r6["mean"]+r6["std"]], "--", label="r6", c=c[6])
            plt.plot([0,data.size-1], [r6["mean"]-r6["std"], r6["mean"]-r6["std"]], "--", label="r6", c=c[6])
            plt.plot([0,data.size-1], [r7["mean"], r7["mean"]], label="r7", c=c[6])
            plt.plot([0,data.size-1], [r7["mean"]+r7["std"], r7["mean"]+r7["std"]], "--", label="r7", c=c[7])
            plt.plot([0,data.size-1], [r7["mean"]-r7["std"], r7["mean"]-r7["std"]], "--", label="r7", c=c[7])
            plt.plot([0,data.size-1], [r8["mean"], r8["mean"]], label="r8", c=c[6])
            plt.plot([0,data.size-1], [r8["mean"]+r8["std"], r8["mean"]+r8["std"]], "--", label="r8", c=c[8])
            plt.plot([0,data.size-1], [r8["mean"]-r8["std"], r8["mean"]-r8["std"]], "--", label="r8", c=c[8])
        #!! second wave
        if secondwave:
            r0 = iterstat1d(data, useMedian=True, maxRej=4)
            if doPrint:
                print("second wave")
                print("r0:   {} +- {} @ {} {}".format(r0["mean"], r0["std"], r0["niter"], r0["median"]))
            rt0 = iterstat1dC(data, useMedian=True, maxRej=4)
            if doPrint:
                print("rt0:  {} +- {} @ {} {}".format(rt0["mean"], rt0["std"], rt0["niter"], rt0["median"]))
            self.assertTrue(self.floatEqual(r0["mean"],   rt0["mean"]))
            self.assertTrue(self.floatEqual(r0["std"],    rt0["std"]))
            self.assertTrue(self.floatEqual(r0["niter"],  rt0["niter"]))
            self.assertTrue(self.floatEqual(r0["median"], rt0["median"]))
            self.assertTrue(np.all(r0["badmask"]==rt0["badmask"]))
            self.assertTrue(np.all(r0["newBadMask"]==rt0["newBadMask"]))

            r1 = iterstat1d(data, mask=badmask, useMedian=True, maxRej=3)
            if doPrint:
                print("r1:   {} +- {} @ {} {}".format(r1["mean"], r1["std"], r1["niter"], r1["median"]))
            rt1 = iterstat1dC(data, mask=badmask, useMedian=True, maxRej=3)
            if doPrint:
                print("rt1:  {} +- {} @ {} {}".format(rt1["mean"], rt1["std"], rt1["niter"], rt1["median"]))
            self.assertTrue(self.floatEqual(r1["mean"],   rt1["mean"]))
            self.assertTrue(self.floatEqual(r1["std"],    rt1["std"]))
            self.assertTrue(self.floatEqual(r1["niter"],  rt1["niter"]))
            self.assertTrue(self.floatEqual(r1["median"], rt1["median"]))
            self.assertTrue(np.all(r1["badmask"]==rt1["badmask"]))
            self.assertTrue(np.all(r1["newBadMask"]==rt1["newBadMask"]))

            r2 = iterstat1d(data, sigma=sigma, useMedian=True, maxRej=2)
            if doPrint:
                print("r2:   {} +- {} @ {} {}".format(r2["mean"], r2["std"], r2["niter"], r2["median"]))
            rt2 = iterstat1dC(data, sigma=sigma, useMedian=True, maxRej=2)
            if doPrint:
                print("rt2:  {} +- {} @ {} {}".format(rt2["mean"], rt2["std"], rt2["niter"], rt2["median"]))
            self.assertTrue(np.all(r2["badmask"]==rt2["badmask"]))
            self.assertTrue(np.all(r2["newBadMask"]==rt2["newBadMask"]))
            self.assertTrue(self.floatEqual(r2["mean"],   rt2["mean"]))
            self.assertTrue(self.floatEqual(r2["std"],    rt2["std"]))
            self.assertTrue(self.floatEqual(r2["niter"],  rt2["niter"]))
            self.assertTrue(self.floatEqual(r2["median"], rt2["median"]))

            r3 = iterstat1d(data, invvar=1/sigma**2, useMedian=True, maxRej=3, maxDev=2)
            if doPrint:
                print("r3:   {} +- {} @ {} {}".format(r3["mean"], r3["std"], r3["niter"], r3["median"]))
            rt3 = iterstat1dC(data, invvar=1/sigma**2, useMedian=True, maxRej=3, maxDev=2)
            if doPrint:
                print("rt3:  {} +- {} @ {} {}".format(rt3["mean"], rt3["std"], rt3["niter"], rt3["median"]))
            self.assertTrue(np.all(r3["badmask"]==rt3["badmask"]))
            self.assertTrue(np.all(r3["newBadMask"]==rt3["newBadMask"]))
            self.assertTrue(self.floatEqual(r3["mean"],   rt3["mean"]))
            self.assertTrue(self.floatEqual(r3["std"],    rt3["std"]))
            self.assertTrue(self.floatEqual(r3["niter"],  rt3["niter"]))
            self.assertTrue(self.floatEqual(r3["median"], rt3["median"]))

            r4 = iterstat1d(data, sigma=sigma, invvar=1/sigma**2, useMedian=True, maxRej=4)
            if doPrint:
                print("r4:   {} +- {} @ {} {}".format(r4["mean"], r4["std"], r4["niter"], r4["median"]))
            rt4 = iterstat1dC(data, sigma=sigma, invvar=1/sigma**2, useMedian=True, maxRej=4)
            if doPrint:
                print("rt4:  {} +- {} @ {} {}".format(rt4["mean"], rt4["std"], rt4["niter"], rt4["median"]))
            self.assertTrue(np.all(r4["badmask"]==rt4["badmask"]))
            self.assertTrue(np.all(r4["newBadMask"]==rt4["newBadMask"]))
            self.assertTrue(self.floatEqual(r4["mean"],   rt4["mean"]))
            self.assertTrue(self.floatEqual(r4["std"],    rt4["std"]))
            self.assertTrue(self.floatEqual(r4["niter"],  rt4["niter"]))
            self.assertTrue(self.floatEqual(r4["median"], rt4["median"]))

            r5 = iterstat1d(data, invvar=1/sigma**2, maxiter=2, useMedian=True, maxRej=5)
            if doPrint:
                print("r5:   {} +- {} @ {} {}".format(r5["mean"], r5["std"], r5["niter"], r5["median"]))
            rt5 = iterstat1dC(data, invvar=1/sigma**2, maxiter=2, useMedian=True, maxRej=5)
            if doPrint:
                print("rt5:  {} +- {} @ {} {}".format(rt5["mean"], rt5["std"], rt5["niter"], rt5["median"]))
            self.assertTrue(np.all(r5["badmask"]==rt5["badmask"]))
            self.assertTrue(np.all(r5["newBadMask"]==rt5["newBadMask"]))
            self.assertTrue(self.floatEqual(r5["mean"],   rt5["mean"]))
            self.assertTrue(self.floatEqual(r5["std"],    rt5["std"]))
            self.assertTrue(self.floatEqual(r5["niter"],  rt5["niter"]))
            self.assertTrue(self.floatEqual(r5["median"], rt5["median"]))

            r6 = iterstat1d(data, invvar=1/sigma**2, maxiter=2, useMedian=True, stdRej=True, maxRej=1)
            if doPrint:
                print("r6:   {} +- {} @ {} {}".format(r6["mean"], r6["std"], r6["niter"], r6["median"]))
            rt6 = iterstat1dC(data, invvar=1/sigma**2, maxiter=2, useMedian=True, stdRej=True, maxRej=1)
            if doPrint:
                print("rt6:  {} +- {} @ {} {}".format(rt6["mean"], rt6["std"], rt6["niter"], rt6["median"]))
            self.assertTrue(np.all(r6["badmask"]==rt6["badmask"]))
            self.assertTrue(np.all(r6["newBadMask"]==rt6["newBadMask"]))
            self.assertTrue(self.floatEqual(r6["mean"],   rt6["mean"]))
            self.assertTrue(self.floatEqual(r6["std"],    rt6["std"]))
            self.assertTrue(self.floatEqual(r6["niter"],  rt6["niter"]))
            self.assertTrue(self.floatEqual(r6["median"], rt6["median"]))

            r8 = iterstat1d(data, useMedian=True, maxRej=1, maxiter=20)
            if doPrint:
                print("r8:   {} +- {} @ {} {}".format(r8["mean"], r8["std"], r8["niter"], r8["median"]))
            rt8 = iterstat1dC(data, useMedian=True, maxRej=1, maxiter=20)
            if doPrint:
                print("rt8:  {} +- {} @ {} {}".format(rt8["mean"], rt8["std"], rt8["niter"], rt8["median"]))
            self.assertTrue(np.all(r8["badmask"]==rt8["badmask"]))
            self.assertTrue(np.all(r8["newBadMask"]==rt8["newBadMask"]))
            self.assertTrue(self.floatEqual(r8["mean"],   rt8["mean"]))
            self.assertTrue(self.floatEqual(r8["std"],    rt8["std"]))
            self.assertTrue(self.floatEqual(r8["niter"],  rt8["niter"]))
            self.assertTrue(self.floatEqual(r8["median"], rt8["median"]))

            r7 = iterstat1d(data, invvar=1/sigma**2, useMedian=True, maxRej=1, maxiter=20)
            if doPrint:
                print("r7:   {} +- {} @ {} {}".format(r7["mean"], r7["std"], r7["niter"], r7["median"]))
            rt7 = iterstat1dC(data, invvar=1/sigma**2, useMedian=True, maxRej=1, maxiter=20)
            if doPrint:
                print("rt7:  {} +- {} @ {} {}".format(rt7["mean"], rt7["std"], rt7["niter"], rt7["median"]))
            self.assertTrue(np.all(r7["badmask"]==rt7["badmask"]))
            self.assertTrue(np.all(r7["newBadMask"]==rt7["newBadMask"]))
            self.assertTrue(self.floatEqual(r7["mean"],   rt7["mean"]))
            self.assertTrue(self.floatEqual(r7["std"],    rt7["std"]))
            self.assertTrue(self.floatEqual(r7["niter"],  rt7["niter"]))
            self.assertTrue(self.floatEqual(r7["median"], rt7["median"]))


        #endregion: tests
        if doPlot:
            plt.legend()
            plt.scatter(np.where(r3["newBadMask"])[0], data[r3["newBadMask"]], c="red", s=300, marker="x")
            plt.show(0)
    @unittest.skip("skip")
    def test_iterstat1d(self):
        N=99999
        self.count = 0
        print()
        for i in range(N):
            self.iterstat1d_test(doPrint=0)
    @unittest.skip("skip")
    def test_iterstat1d_2(self):
        data = [0.008567, 0.09059, 0.2519, 0.5506, 0.2691, 2.149, 0.3507, 0.7372, 0.08919, 1.125, 1.912, 0.5372, 0.9009, 0.9384, 0.4078, 0.3988, 35.2, 39.56, 0.2773, 0.06454, 0.03878]
        data = np.array(data)
        mask = data<0
        result = iterstat1d(data, mask=mask, maxiter=1)
        newBadMask =  result['newBadMask']
        print()
        print(("".join(["{:9.3e},"]*data.size)).format(*(data.tolist())))
        print(("".join(["{:^10d}"]*newBadMask.size)).format(*(newBadMask.tolist())))
        print(result)
    @unittest.skip("skip")
    def test_linearInterpWithError(self):
        N=100+1
        x = np.arange(N) + (np.random.random(N)-0.5)
        xx = np.linspace(-100,N, 30000)
        f = lambda x:10*np.sin(0.02*x+np.random.random(1)) + 6*np.cos(0.03*x+np.random.random(1)) + 0.001*(x-500)*2
        y = f(x)
        sigma = np.sqrt(np.abs(y) + 10)
        result = linearInterpWithError(x, y, sigma, xx)
        yy = result["y"]
        sigmasigma = result["sigma"]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(x, y, color='red')
        ax.scatter(x, sigma, color='green')
        ax.plot(xx, yy, "r-")
        ax.plot(xx, sigmasigma, "g-")
        plt.show(0)
    def optimalFixExtract1D_test(self, doPlot=False, doPrint=False, cosmicSize=0):
        mpl.rcParams['figure.subplot.top'] = 0.8
        mpl.rcParams['figure.subplot.bottom']= 0.05
        mpl.rcParams['figure.subplot.left']  = 0.05
        mpl.rcParams['figure.subplot.right'] = 0.95
        N = 41
        xstart = -(N-1)/2
        xend = (N-1)/2
        xmid = int((N-1)/2)
        halfwidth = 10
        xxstart = xmid-halfwidth
        xxend = xmid+halfwidth+1

        x = np.linspace(xstart,xend,N)
        xdense = np.linspace(xstart, xend, 999)
        sigma = np.random.uniform(1,3)
        A = np.random.uniform(50,5000)
        y = gaussian(x, A, 0, sigma)
        RDNOISE = 2
        in_var = np.abs(y)+RDNOISE
        noise = np.random.randn(N) * np.sqrt(in_var)
        xx = x[xxstart:xxend]
        yy = y[xxstart:xxend]
        y = y  + noise
        ydense = gaussian(xdense, A, 0, sigma)
        iny = y[xxstart:xxend]
        inin_var = in_var[xxstart:xxend]
        if cosmicSize>0:
            maxpeak = yy.max()
            pos = np.random.randint(xxstart, xxend)
            cosmic = np.random.uniform(0.8*maxpeak, 3*maxpeak, cosmicSize)
            y[pos:pos+cosmicSize] = cosmic[:]
        result = optimalFixExtract1D(xx,iny,inin_var, yy)
        C = result["C"]
        varC = result["varC"]
        outprofile = result["outy"]
        flux = result['flux']
        fluxvar = result['fluxvar']
        boxsum = result['boxsum']
        boxsumvar = result['boxsumvar']
        chi2 = result['chi2']
        chiss = result['chiss']
        chissig = result['chissig']
        reduce_chi2 = result['reduce_chi2']
        newBadMask = result['newBadMask']
        newmask = result['newmask']

        badpixelFlag = result['badpixelFlag']

        if doPlot:
            title = []
            title += ["A:     {:15.7f} sigma: {:15.7f}".format(A, sigma)]
            title += ["flux:  {:15.7f} fluxvar:  {:15.7f} {:15.7}~{:<15.7} C:{:15.7f} varC:{:15.7f}".
                       format(flux, fluxvar, flux-np.sqrt(fluxvar), flux+np.sqrt(fluxvar), C, varC)]
            title += ["boxsum:{:15.7f} boxsumvar:{:15.7f} {:15.7}~{:<15.7}".
                      format(boxsum, boxsumvar, boxsum-np.sqrt(boxsumvar), boxsum+np.sqrt(boxsumvar))]
            title += ["chi2:  {:15.7f} reduce_chi2: {:15.7f}".format(chi2, reduce_chi2)]
            title += ["chissig: {}".format(chissig)]
            title += ["chiss: "+ "".join((["{:5.3e} "]*chiss.size)).format(*chiss)]
            title += ["newmask: "+ "".join((["{:2d}"]*newmask.size)).format(*newmask)]
            title += ["newbadmask: "+ "".join((["{:2d}"]*newBadMask.size)).format(*newBadMask)]
            title = "\n".join(title)
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_title(title, loc='left')
            ax.scatter(x,y, edgecolor="None", s=20)
            ax.plot(xdense, ydense, color="green")
            ax.plot(xx, outprofile, color="red")
            if badpixelFlag:
                ax.scatter(xx[newBadMask], iny[newBadMask], marker='x', s=40)
            plt.show(0)
    #@unittest.skip("skip")
    def test_optimalFixExtract1D(self, doPlot=True):
        self.optimalFixExtract1D_test(doPlot=1, doPrint=1, cosmicSize=3)
        self.optimalFixExtract1D_test(doPlot=1, doPrint=1, cosmicSize=3)
        self.optimalFixExtract1D_test(doPlot=1, doPrint=1, cosmicSize=3)
        self.optimalFixExtract1D_test(doPlot=1, doPrint=1, cosmicSize=3)
        self.optimalFixExtract1D_test(doPlot=1, doPrint=1, cosmicSize=3)
    @unittest.skip("skip")
    def iterstat3d_test(self):
        obj = pyfits.open("./test/test3dstack.fits")
        data = obj[0].data
        var = obj[1].data
        mask = np.ones(data.shape)<0
        N, ny, nx = data.shape
        result = iterstat3d(data, invvar=1/var, sigrej=10)
        p = [pyfits.PrimaryHDU(result["mean"])]
        d = []
        for i in range(N):
            originData = data[i]
            originMaskedData = originData.copy()
            originMaskedData[result["badmask"][i]]=np.nan
            d.append(pyfits.ImageHDU(originData))
            d.append(pyfits.ImageHDU(originMaskedData))
        last = [pyfits.ImageHDU(np.logical_not(result["badmask"]).sum(axis=0))]
        fits  = pyfits.HDUList(p+d+last)
        fits.writeto("test/test_iterstat3d.fits", clobber=True)
        return result
    @unittest.skip("skip")
    def iterstat3dC_test(self):
        obj = pyfits.open("./test/test3dstack.fits")
        data = obj[0].data
        var = obj[1].data
        mask = np.ones(data.shape)<0
        N, ny, nx = data.shape
        result = iterstat3dC(data, invvar=1/var, sigrej=10)
        p = [pyfits.PrimaryHDU(result["mean"])]
        d = []
        for i in range(N):
            originData = data[i]
            originMaskedData = originData.copy()
            originMaskedData[result["badmask"][i]]=np.nan
            d.append(pyfits.ImageHDU(originData))
            d.append(pyfits.ImageHDU(originMaskedData))
        last = [pyfits.ImageHDU(np.logical_not(result["badmask"]).sum(axis=0))]
        fits  = pyfits.HDUList(p+d+last)
        fits.writeto("test/test_iterstat3dC.fits", clobber=True)
        return result
    @unittest.skip("skip")
    def test_iterstat3d_C(self):
        result0 = self.iterstat3d_test()
        result1 = self.iterstat3dC_test()
        self.assertTrue(np.all(result0['badmask']==result1['badmask']))
        self.assertTrue(np.all(result0['newBadMask']==result1['newBadMask']))
        t0=result0['mean']
        t1=result1['mean']
        self.assertTrue(np.all((t0==t1)|((~np.isfinite(t0))&(~np.isfinite(t1)))))
        t0=result0['std']
        t1=result1['std']
        self.assertTrue(np.all((t0==t1)|((~np.isfinite(t0))&(~np.isfinite(t1)))))
    @unittest.skip("skip")
    def test_fillHollowness(self):
        data = np.zeros((100,100))
        data_result = np.zeros((100,100))
        pos = np.random.randint(1,98, (10,2))
        N = 10
        for eachy, eachx in pos.tolist():
            data_result[eachy, eachx] = 1
            direc = np.random.randint(0,3)
            if direc==0:
                data[eachy-1, eachx] = 1
                data[eachy, eachx+1] = 1
                data[eachy, eachx-1] = 1
                data_result[eachy-1, eachx] = 1
                data_result[eachy, eachx+1] = 1
                data_result[eachy, eachx-1] = 1
            elif direc==1:
                data[eachy+1, eachx] = 1
                data[eachy, eachx+1] = 1
                data[eachy, eachx-1] = 1
                data_result[eachy+1, eachx] = 1
                data_result[eachy, eachx+1] = 1
                data_result[eachy, eachx-1] = 1
            elif direc==2:
                data[eachy+1, eachx] = 1
                data[eachy-1, eachx] = 1
                data[eachy, eachx-1] = 1
                data_result[eachy+1, eachx] = 1
                data_result[eachy-1, eachx] = 1
                data_result[eachy, eachx-1] = 1
            elif direc==3:
                data[eachy+1, eachx] = 1
                data[eachy-1, eachx] = 1
                data[eachy, eachx+1] = 1
                data_result[eachy+1, eachx] = 1
                data_result[eachy-1, eachx] = 1
                data_result[eachy, eachx+1] = 1
        newdata = fillHollowness(data)
        fig = plt.figure()
        plt.imshow(newdata)
        fig = plt.figure()
        plt.imshow(data_result)
        plt.show(0)
        #self.assertTrue(np.all(newdata==data_result))

def fixpixel(data, mask, kind="cubic", bad=None):
    """mask is good mask"""
    assert mask.dtype==np.bool
    ny, nx = mask.shape
    bady = np.unique(np.where(np.logical_not(mask))[0])
    result = data.copy()
    if kind in ['slinear','quadratic','cubic']:
        fill_value = np.nan
    else:
        fill_value = 'extrapolate'
    for eachy in bady:
        eachxx = np.arange(nx)
        eachyy = data[eachy]
        eachmask = mask[eachy]
        try:
            eachobj = interp1d(eachxx[eachmask], eachyy[eachmask], kind=kind, bounds_error=0, fill_value=fill_value)
            result[eachy,:] = eachobj(eachxx)[:]
        except Exception as e:
            try:
                eachobj = interp1d(eachxx[eachmask], eachyy[eachmask], kind='linear', bounds_error=0, fill_value=fill_value)
                result[eachy,:] = eachobj(eachxx)[:]
            except Exception as ee:
                print(e)
                print(ee)
                if bad is not None:
                    result[eachy,:] = bad
                #import ipdb
                #ipdb.set_trace()
                print()
    return result

def splineWithReject(wave, data, std,
                  tpoints=1500,
                  mindev = 2.5,
                  sigrej = None,
                  sigrejMin = 7,
                  maxEachReject=None,
                  maxiter=500):
    goodmask = np.isfinite(std)&np.isfinite(data)
    finalIndex = np.where(goodmask)[0]
    finalmask = np.zeros(wave.shape, dtype=np.bool)
    if maxEachReject is None:
        maxEachReject = int(wave.size/9999)+1
    if sigrej is None:
        sigrej = np.linspace(20,sigrejMin,maxiter)
    wave = wave[goodmask]
    data = data[goodmask]
    std  =  std[goodmask]
    assert np.all(std>0)
    rejections = []
    niter = 0
    while niter<maxiter:
        thissigrej = sigrej[niter]
        niter += 1
        #!! get ts
        snr = data/std
        snr[snr<0] = -snr[snr<0]
        thisgoodmask = snr>0
        snr = np.power(snr, 2./3)
        if not np.all(np.isfinite(snr)):
            import ipdb
            ipdb.set_trace()
        snrsum = snr.sum()
        cumsum = np.cumsum(snr)
        tsums = np.linspace(0,snrsum,tpoints+2)[1:-1]
        tindexs = []
        for eachcum in tsums:
            tindexs.append(cumsum.searchsorted(eachcum))
        tindex = np.array(tindexs)
        ts = wave[tindex]
        if mindev is not None:
            goodmask = np.concatenate((np.diff(ts)>mindev, np.array([True])))
            ts = ts[goodmask]
        #!! do fitting and reject
        try:
            #thisfit = LSQUnivariateSpline(wave, data, ts, 1/std**2)
            thisfit = LSQUnivariateSpline(wave, data, ts)
        except Exception as e:
            print(e)
            import ipdb
            ipdb.set_trace()
        #thisfit = UnivariateSpline(wave, data, w=1/std**2, s=s)
        #ts = thisfit.get_knots()
        thisy = thisfit(wave)
        thisdiffstd = np.abs(thisy-data)/std
        thisgoodmask = thisdiffstd<thissigrej
        if (~thisgoodmask).sum()>maxEachReject:
            sortIndex = np.argsort(thisdiffstd)[:-maxEachReject]
            thisgoodmask[sortIndex]=1
        rejections.append((~thisgoodmask).sum())
        wave = wave[thisgoodmask]
        data = data[thisgoodmask]
        std  =  std[thisgoodmask]
        finalIndex = finalIndex[thisgoodmask]
    finalmask[finalIndex] = 1
    return {"mask":finalmask,
            'fit':thisfit,
            "niter":niter,
            "ts":ts,
            'rejections':rejections}





if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib import cm
    import sys
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMath)
    unittest.TextTestRunner(verbosity=2).run(suite)


#def iterstat1dC(data, sigma=None, invvar=None, sigrej=None, maxiter=10, mask=None, median=False):
    #"""mask is badmask"""
    #if data.size==1:
        #if invvar is not None:
            #if mask is not None:
                #return {"mean": data, "std":None, "median": data, "mask":(invvar==0)*mask}
            #else:
                #return {"mean": data, "std":None, "median": data, "mask":mask}
        #if sigma is not None:
            #if mask is not None:
                #return {"mean": data, "std":None, "median": data, "mask":np.logical_not(np.isfinite(sigma))*mask}
            #else:
                #return {"mean": data, "std":None, "median": data, "mask":mask}
    #if sigrej is None:
        #if sigma is None and invvar is None:
            #sigrej = 3
        #else:
            #sigrej = 5

    #if sigma is not None and invvar is not None:
        #assert data.shape==invvar.shape
        #assert np.all(np.isfinite(invvar))
        #logger.warn("both sigma and invvar set, use invvar and ignore sigma")
        #badmask = (invvar==0)
        #calSigma = False
    #elif invvar is None and sigma is not None:
        #invvar = 1./(sigma**2)
        #assert data.shape==invvar.shape
        #assert np.all(np.isfinite(invvar))
        #badmask = (invvar==0)
        #calSigma = False
    #elif invvar is not None and sigma is None:
        #assert data.shape==invvar.shape
        #badmask = (invvar==0)
        #calSigma = False
    #else:# both sigma and invvar not set, leave invvar None, and set all data good
        #invvar = np.ones(data.shape, dtype=np.float64) * -1
        #badmask = np.zeros(data.shape, dtype=np.short)
        #calSigma = True
    #if mask is not None:
        #assert mask.dtype==np.bool, "mask dtype must be bool"
        #badmask[mask] = 1
    #else:
        #if calSigma:# no mask and sigma and invvar, init mask is all good
            #mask = np.zeros(data.shape, dtype=np.bool)
        #else:
            #mask = (invvar==0) # init mask with invvar

    #initBadMask = badmask.copy()
    #data = data.copy()
    #initNanIndex = np.isnan(data) # record init nan position
    #data[initNanIndex] = 0 # to remove NaN

    #data = data.astype(np.float64, copy=True, order="C")
    #invvar = invvar.astype(np.float64, order="C")
    #initBadMask = initBadMask.astype(np.short, order="C")
    #badmask = badmask.astype(np.short, order="C")
    #newBadMask = np.zeros(data.shape, dtype=np.short, order="C")
    #result = np.zeros(3, dtype=np.float64, order="C")

    #data_c = data.ctypes.data_as(c_double_p)
    #invvar_c = invvar.ctypes.data_as(c_double_p)
    #initBadMask_c = initBadMask.ctypes.data_as(c_short_p)
    #badmask_c = badmask.ctypes.data_as(c_short_p)
    #newBadMask_c = newBadMask.ctypes.data_as(c_short_p)
    #result_c = result.ctypes.data_as(c_double_p)

    ##int iterstat1d(const double data[], int M,
    ##           const double invvar[], double sigrej,
    ##           int maxiter, short const initmask[], short badmask[], short newBadMask[],
    ##           short calSigma, short calMedian,
    ##           double result[]) {
    #iterstat1d_c.argtypes = [c_double_p, c_int,
                             #c_double_p, c_double,
                             #c_int, c_short_p, c_short_p, c_short_p,
                             #c_short, c_short,
                             #c_double_p]
    ## main c function
    #niter = iterstat1d_c(data_c, data.size,
                         #invvar_c, sigrej,
                         #maxiter, initBadMask_c, badmask_c, newBadMask_c,
                         #calSigma, median,
                         #result_c)
    #mean = result[0]
    #std  = result[1]
    #theMedian = result[2]

    #result = {"mean": mean, "std":std, "badmask": badmask.astype(np.bool), "niter":niter, "newBadMask": newBadMask.astype(np.bool)}
    #if niter==maxiter:
        #logger.warn("iterstat reach max iteration number")
    #if median:
        #result.update({"median":theMedian})

    #return result
## need debug...
#def iterstat2d(data, sigma=None, invvar=None, sigrej=None, maxiter=10, mask=None, median=False):
    #"""mask is badmask"""

    #if sigma is not None and invvar is not None:
        #assert data.shape==invvar.shape
        #assert np.all(np.isfinite(invvar))
        #logger.warn("both sigma and invvar set, use invvar and ignore sigma")
        #badmask = (invvar==0)
    #elif invvar is None and sigma is not None:
        #invvar = 1./(sigma**2)
        #assert data.shape==invvar.shape
        #assert np.all(np.isfinite(invvar))
        #badmask = (invvar==0)
    #elif invvar is not None and sigma is None:
        #assert data.shape==invvar.shape
        #badmask = (invvar==0)
    #else:# both sigma and invvar not set, leave invvar None, and set all data good
        #invvar = np.array([None] * data.shape[1])
        #badmask = np.zeros(data.shape, dtype=np.bool)

    #if mask is not None:
        #assert mask.shape == data.shape
        #assert mask.dtype==np.bool, "mask dtype must be bool"
    #else:
        #mask = np.array([None] * data.shape[1])

    #initNanIndex = np.isnan(data) # record init nan position

    #data = data.T
    #data[initNanIndex] = 0 # to remove NaN
    #invvar = invvar.T
    #mask = mask.T

    #allResults = []
    #for i in range(data.shape[0]):
        #eachData, eachInvvar, eachMask = data[i], invvar[i], mask[i]
        #allResults.append(iterstat1d(eachData, invvar=eachInvvar,
                                     #sigrej=sigrej, maxiter=maxiter, mask=eachMask, median=median))

    #mean = np.array([each["mean"] for each in allResults])
    #std = np.array([each["std"] for each in allResults])
    #badmask = np.array([each["badmask"] for each in allResults])
    #niter = np.array([each["niter"] for each in allResults])
    #newBadMask = np.array([each["newBadMask"] for each in allResults])

    #result = {"mean": mean.T, "std":std.T, "badmask": badmask.T, "niter":niter, "newBadMask": newBadMask.T}

    #if median:
        #allMedian = np.array([each["median"] for each in allResults])
        #result.update({"median": np.array(allMedian).T})

    #data = data.T
    #data[initNanIndex] = np.nan
    #return result
#def iterstat2d_test():
    ##region: set data and figures
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #c = cm.get_cmap("rainbow")
    #c = c(np.linspace(0,1,6))
    ## test for iterstat2d
    #m = 300
    #n = 50

    #x = np.arange(m)
    #y = 2*np.sin(0.141*x) + 3.1*np.cos(0.312*x)
    #data2d = []
    #badmask2d = []
    #sigma2d = []
    #badPointMask2d = []
    #for i in range(n):
        #data = y + np.random.randn(x.size)
        #sigma = (np.random.random(x.size)-0.5)/2+1
        #badmask = np.zeros(data.shape, dtype=np.bool)
        #badPointMask = np.zeros(data.shape, dtype=np.bool)

        ## add bad points
        #badIndex = np.random.randint(0, high=data.size-1, size=30)
        #for eachIndex in badIndex:
            #sigma[eachIndex] = np.inf
            #badmask[eachIndex] = 1
            #badPointMask[eachIndex] = 1

        ## add singular points
        #badIndex = np.random.randint(0, high=data.size-1, size=60)
        #for eachIndex in badIndex:
            #data[eachIndex] += (np.random.random()-0.5)*20

        ## add nan points
        #badIndex = np.random.randint(0, high=data.size-1, size=10)
        #for eachIndex in badIndex:
            #data[eachIndex] = np.nan
            #sigma[eachIndex] = np.inf
            #badmask[eachIndex] = 1

        #data2d.append(data)
        #sigma2d.append(sigma)
        #badmask2d.append(badmask)
        #badPointMask2d.append(badPointMask)

    #data2d = np.array(data2d)
    #sigma2d = np.array(sigma2d)
    #badmask2d = np.array(badmask2d)
    #badPointMask2d = np.array(badPointMask2d)

    #for data,sigma,badPointMask in zip(data2d, sigma2d,badPointMask2d):
        #plt.scatter(range(data.size), data)
        #plt.errorbar(range(data.size), data, ls="None", yerr=sigma)
        #plt.scatter(np.where(badPointMask)[0], data[badPointMask], color="r")
    ##endregion: set data and figures
    ##region: tests
    #r0 = iterstat2d(data2d, calMedian=True)
    ##print("r0:   {} +- {} @ {} {}".format(r0["mean"], r0["std"], r0["niter"], r0["median"]))
    ##rt0 = iterstat1dC(data2d, calMedian=True)
    ##print("rt0:  {} +- {} @ {} {}".format(rt0["mean"], rt0["std"], rt0["niter"], rt0["median"]))
    ##assert np.all(r0["badmask"]==rt0["badmask"])
    ##assert np.all(r0["newBadMask"]==rt0["newBadMask"])
    #r1 = iterstat2d(data2d, mask=badmask2d, calMedian=True)
    ##print("r1:   {} +- {} @ {} {}".format(r1["mean"], r1["std"], r1["niter"], r1["median"]))
    ##rt1 = iterstat1dC(data2d, mask=badmask, calMedian=True)
    ##print("rt1:  {} +- {} @ {} {}".format(rt1["mean"], rt1["std"], rt1["niter"], rt1["median"]))
    ##assert np.all(r1["badmask"]==rt1["badmask"])
    ##assert np.all(r1["newBadMask"]==rt1["newBadMask"])
    #r2 = iterstat2d(data2d, sigma=sigma2d, calMedian=True)
    ##print("r2:   {} +- {} @ {} {}".format(r2["mean"], r2["std"], r2["niter"], r2["median"]))
    ##rt2 = iterstat1dC(data2d, sigma=sigma, calMedian=True)
    ##print("rt2:  {} +- {} @ {} {}".format(rt2["mean"], rt2["std"], rt2["niter"], rt2["median"]))
    ##assert np.all(r2["badmask"]==rt2["badmask"])
    ##assert np.all(r2["newBadMask"]==rt2["newBadMask"])
    #r3 = iterstat2d(data2d, invvar=1/sigma2d**2, calMedian=True)
    ##print("r3:   {} +- {} @ {} {}".format(r3["mean"], r3["std"], r3["niter"], r3["median"]))
    ##rt3 = iterstat1dC(data2d, invvar=1/sigma**2, calMedian=True)
    ##print("rt3:  {} +- {} @ {} {}".format(rt3["mean"], rt3["std"], rt3["niter"], rt3["median"]))
    ##assert np.all(r3["badmask"]==rt3["badmask"])
    ##assert np.all(r3["newBadMask"]==rt3["newBadMask"])
    #r4 = iterstat2d(data2d, sigma=sigma2d, invvar=1/sigma2d**2, calMedian=True)
    ##print("r4:   {} +- {} @ {} {}".format(r4["mean"], r4["std"], r4["niter"], r4["median"]))
    ##rt4 = iterstat1dC(data2d, sigma=sigma, invvar=1/sigma**2, calMedian=True)
    ##print("rt4:  {} +- {} @ {} {}".format(rt4["mean"], rt4["std"], rt4["niter"], rt4["median"]))
    ##assert np.all(r4["badmask"]==rt4["badmask"])
    ##assert np.all(r4["newBadMask"]==rt4["newBadMask"])
    #r5 = iterstat2d(data2d, invvar=1/sigma2d**2, maxiter=2, calMedian=True)
    ##print("r5:   {} +- {} @ {} {}".format(r5["mean"], r5["std"], r5["niter"], r5["median"]))
    ##rt5 = iterstat1dC(data2d, invvar=1/sigma**2, maxiter=2, calMedian=True)
    ##print("rt5:  {} +- {} @ {} {}".format(rt5["mean"], rt5["std"], rt5["niter"], rt5["median"]))
    ##assert np.all(r5["badmask"]==rt5["badmask"])
    ##assert np.all(r5["newBadMask"]==rt5["newBadMask"])
    #plt.plot(x, r1["mean"], label="r1", c=c[1])
    #plt.plot(x, r1["mean"]+r1["std"], "--", label="r1", c=c[1])
    #plt.plot(x, r1["mean"]-r1["std"], "--", label="r1", c=c[1])
    #plt.plot(x, r2["mean"], label="r2", c=c[2])
    #plt.plot(x, r2["mean"]+r2["std"], "--", label="r2", c=c[2])
    #plt.plot(x, r2["mean"]-r2["std"], "--", label="r2", c=c[2])
    #plt.plot(x, r3["mean"], label="r3", c=c[3])
    #plt.plot(x, r3["mean"]+r3["std"], "--", label="r3", c=c[3])
    #plt.plot(x, r3["mean"]-r3["std"], "--", label="r3", c=c[3])
    #plt.plot(x, r4["mean"], label="r4", c=c[4])
    #plt.plot(x, r4["mean"]+r4["std"], "--", label="r4", c=c[4])
    #plt.plot(x, r4["mean"]-r4["std"], "--", label="r4", c=c[4])
    #plt.plot(x, r5["mean"], label="r5", c=c[5])
    #plt.plot(x, r5["mean"]+r5["std"], "--", label="r5", c=c[5])
    #plt.plot(x, r5["mean"]-r5["std"], "--", label="r5", c=c[5])
    ##endregion: tests
    #plt.legend()
    #for i in range(n):
        #plt.scatter(np.where(r3["newBadMask"][i])[0], data2d[i][r3["newBadMask"][i]], c="red", s=300, marker="x")
    #fig.tight_layout()
    #plt.show(0)
## reject with error
#def reject(data, model, badmask=None,
           #sigma=None, invvar=None,
           #lower=None, upper=None, maxdev=None,
           #maxrej=None
           #):
    #"""reject 1d data"""

    #if sigma is None:
        #if invvar is None:
            #sigma = std(data)

    #if data.ndim==1:
        #pass
# old method backup
#def optimalFixExtract1D(x, y, var, profile, maxBad = 3, doPrint=False):
    #assert x.shape==y.shape==var.shape==profile.shape
    #allGood = np.isfinite(var)
    #N = x.size
    #niter = N*(maxBad) + int(maxBad*(maxBad+1)/2)

    #records = []
    #for nbad in range(0,maxBad+1):
        #for startIndex in range(0, N-nbad+1):
            #mask = allGood.copy()
            #mask[startIndex:startIndex+nbad] = 0
            #goodN = mask.sum()
            #tFlux = np.sum(profile[mask])
            #denom = np.sum((profile[mask]**2/var[mask]))
            #C = np.sum(y[mask]*profile[mask]/var[mask])/denom
            #varC = 1/denom
            #chi2 = np.sum((C*profile[mask]-y[mask])**2/var[mask])
            #reduce_chi2 = chi2/(goodN-1)
            #flux = C*tFlux
            #fluxvar = varC*tFlux**2
            #if doPrint:
                #print("flux:  {:15.7f} fluxvar:  {:15.7f} {:15.7}~{:<15.7} C:{:15.7f} varC:{:15.7f}".
                        #format(flux, fluxvar, flux-np.sqrt(fluxvar), flux+np.sqrt(fluxvar), C, varC))
                #print("chi2:  {:15.7f} reduce_chi2: {:15.7f}".format(chi2, reduce_chi2))
                #print("mask:  {}".format(mask.astype(np.int)))
            #records.append((goodN, flux, fluxvar, chi2, reduce_chi2))
            #if nbad==0:
                #break
    #boxsum = y.sum()
    #boxsumvar = var.sum()
    #if doPrint:
        #print("boxsum:{:15.7f} boxsumvar:{:15.7f} {:15.7}~{:<15.7}".
                  #format(boxsum, boxsumvar, boxsum-np.sqrt(boxsumvar), boxsum+np.sqrt(boxsumvar)))
        #lgoodN, lflux, lfluxvar, lchi2, lreduce_chi2 = list(zip(*records))
        #print(lflux)
        #print(lreduce_chi2)
    #return {"C":C, "varC":varC, "outy":profile*C,
            #"flux":flux, "fluxvar":fluxvar,
            #"boxsum":boxsum, "boxsumvar":boxsumvar,
            #"chi2":chi2, "reduce_chi2":reduce_chi2}
