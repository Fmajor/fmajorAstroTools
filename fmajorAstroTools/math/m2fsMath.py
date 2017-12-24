from __future__ import (absolute_import, division, print_function, unicode_literals)
import numpy as np
import ctypes
from ctypes import c_int, c_double, c_short
import warnings
import os
#from ..utils import messager
from fmajorAstroUtils.utils import messager
import glob
import time
import ipdb
import astropy.io.fits as pyfits


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

so = ctypes.CDLL(os.path.join(os.path.dirname(__file__), "lib", "astroMath.so"))
#so = ctypes.CDLL("./astroMath.so")
iterstat1d_c = so.iterstat1d
iterstat3d_c = so.iterstat3d
linearInterpWithError_c = so.linearInterpWithError

# As comments for C version
# all mask is "badmask"
def iterstat1d(data, invvar=None, sigma=None, mask=None, sigrej=None, maxiter=10,
               calMedian=False, maxRej=None, stdRej=False, maxDev=None):
    """the priority is invvar > sigma > mask
    """
    #region: setups
    if data.size==1:
        if invvar is not None:# use invvar
            return {"mean": data, "std":np.nan, "median": data,
                    "badmask":invvar==0,
                    "newBadMask": np.zero(data.shape, dtype=np.bool)}
        elif sigma is not None:
            return {"mean": data, "std":np.nan, "median": data,
                    "badmask": np.isfinite(sigma),
                    "newBadMask": np.zeros(data.shape, dtype=np.bool)}
        elif mask is not None:
            return {"mean": data, "std":None, "median": data,
                    "mask":mask, "newBadMask": np.zeros(data.shape, dtype=np.bool)}
        else:
            return {"mean": data, "std":None, "median": data,
                    "mask": np.zeros(data.shape, dtype=np.bool),
                    "newBadMask": np.zeros(data.shape, dtype=np.bool)}

    # data size larger than 1
    if sigrej is None:
        if sigma is None and invvar is None:
            sigrej = 3
        else:
            sigrej = 5

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
        badmask = mask.copy()
        stdRej = True
    else:
        badmask = np.zeros(data.shape, dtype=np.bool)
        stdRej = True

    #endregion: setups

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
        if invvar is not None:
            mean = np.sum(data*invvar)/np.sum(invvar)
        else:
            mean = np.sum(data*goodmask)/goodN
        if goodN>1:
            std = np.sqrt(np.sum(goodmask * (data-mean)**2) / (goodN - 1))
        else:
            std = np.nan
    if np.isfinite(std):
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

        newBadMask[badmask] = False
        newBadMaskNum = newBadMask.sum()
    else:
        newBadMaskNum = 0

    niter = 0
    #print("\t python: iter: {} mean: {:.10f} std: {:.10f} badnumber: {} newBadNumber: {}".format(niter, mean, std, badmask.sum(), newBadMaskNum));
    #!! main loop
    while niter<maxiter and newBadMaskNum:
        niter += 1
        if maxRej is not None:
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

        newBadMask[badmask] = False
        newBadMaskNum = newBadMask.sum()

        #print("\t python: iter: {} mean: {:.10f} std: {:.10f} badnumber: {} newBadNumber: {}".format(niter, mean, std, badmask.sum(), newBadMaskNum));

    newBadMask = (badmask.astype(np.int) - initBadMask.astype(np.int)).astype(np.bool)
    result = {"mean": mean, "std":std, "badmask": badmask, "niter":niter, "newBadMask": newBadMask}
    if niter==maxiter:
        logger.warn("iterstat reach max iteration number")
    if calMedian:
        goodmask = np.logical_not(badmask)
        #print(np.sort(data[goodmask]))
        #print("\t\t\tpython: ", data[goodmask].sum(), goodmask.sum())
        if allBad:
            result["median"] = np.nan
        else:
            result["median"] = np.median(data[goodmask])

    return result

# C version of last function
def iterstat1dC(data, invvar=None, sigma=None, mask=None, sigrej=None, hsigrej=None, lsigrej=None,
                maxiter=10, calMedian=False, maxRej=None, stdRej=False, maxDev=None):
    """the priority is invvar > sigma > mask
    """
    #region: setups
    if data.size==1:
        if invvar is not None:# use invvar
            return {"mean": data, "std":np.nan, "median": data,
                    "badmask":invvar==0,
                    "newBadMask": np.zero(data.shape, dtype=np.bool)}
        elif sigma is not None:
            return {"mean": data, "std":np.nan, "median": data,
                    "badmask": np.isfinite(sigma),
                    "newBadMask": np.zeros(data.shape, dtype=np.bool)}
        elif mask is not None:
            return {"mean": data, "std":None, "median": data,
                    "mask":mask, "newBadMask": np.zeros(data.shape, dtype=np.bool)}
        else:
            return {"mean": data, "std":None, "median": data,
                    "mask": np.zeros(data.shape, dtype=np.bool),
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
    #           short calSigma, short calMedian, double maxdev,
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
                         stdRej, calMedian, maxRej, maxDev,
                         result_c)
    mean = result[0]
    std  = result[1]
    theMedian = result[2]

    result = {"mean": mean, "std":std, "badmask": badmask, "niter":niter, "newBadMask": newBadMask}
    if niter==maxiter:
        logger.warn("iterstat reach max iteration number")
    if calMedian:
        result["median"] = theMedian

    return result

def iterstat3d(data, invvar=None, sigma=None, mask=None, sigrej=None, maxiter=10,
               calMedian=False, maxRej=None, stdRej=False, maxDev=None):
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
    a, b, c = data.shape

    #!! ./test/testArrayReshape.py
    #!! show some numpy reshape exam that can help to understand thing below
    data = data.T
    if invvar is not None:
        invvar = invvar.T
    else:
        badmask = badmask.T

    result_mean = np.empty((c,b))
    result_std  = np.empty((c,b))
    result_badmask  = np.empty((c,b,a), dtype=np.bool)
    result_newBadMask  = np.empty((c,b,a), dtype=np.bool)

    if calMedian:
        result_median = np.empty((c,b))
    #endregion: data 2 ctypes

    for eachC in range(c):
        for eachB in range(b):
            print("{}/{} {}/{}".format(eachC, c, eachB, b))
            if invvar is not None:
                aux = iterstat1dC(data[eachC, eachB], invvar=invvar[eachC, eachB], sigrej=sigrej, maxiter=maxiter, calMedian=calMedian, maxRej=maxRej, stdRej=stdRej, maxDev=maxDev)
            else:
                aux = iterstat1dC(data[eachC, eachB], mask=badmask[eachC, eachB], sigrej=sigrej, maxiter=maxiter, calMedian=calMedian, maxRej=maxRej, stdRej=stdRej, maxDev=maxDev)
            #!! record result
            result_mean[eachC, eachB] = aux["mean"]
            result_std[eachC, eachB] = aux["std"]
            result_badmask[eachC, eachB, :] = aux["badmask"][:]
            result_newBadMask[eachC, eachB, :] = aux["newBadMask"][:]
            if calMedian:
                result_median[eachC, eachB] = aux["median"]

    result = {"mean": result_mean.T, "std":result_std.T,
              "badmask": result_badmask.T,
              "newBadMask": result_newBadMask.T}
    if calMedian:
        result["median"] = result_median.T

    return result

def iterstat3dC(data, invvar=None, sigma=None, mask=None, sigrej=None, hsigrej=None, lsigrej=None,
                maxiter=10, calMedian=False, maxRej=None, stdRej=False, maxDev = None):
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
    #           short stdRej, short calMedian, int maxRej,
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
                 stdRej, calMedian, maxRej, maxDev,
                 result_c)

    cresult = cresult.T
    badmask = badmask.T.astype(np.bool)
    newBadMask = newBadMask.T.astype(np.bool)

    result = {"mean": cresult[0], "std":cresult[1],
              "badmask": badmask, "newBadMask": newBadMask}
    if calMedian:
        result["median"] = cresult[2]

    return result

def linearInterpWithError(x, y, sigma, xx):
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

# reject with error
def reject(data, model, badmask=None,
           sigma=None, invvar=None,
           lower=None, upper=None, maxdev=None,
           maxrej=None
           ):
    """reject 1d data"""

    if sigma is None:
        if invvar is None:
            sigma = std(data)

    if data.ndim==1:
        pass


import matplotlib.pyplot as plt
from matplotlib import cm
def iterstat1d_test():
    #region: set data and figures
    fig = plt.figure()
    ax = fig.add_subplot(111)
    c = cm.get_cmap("rainbow")
    c = c(np.linspace(0,1,7))
    # test for iterstat1d
    data = np.random.randn(80)+10
    sigma = (np.random.random(80)-0.5)/2+1
    badmask = np.zeros(data.shape, dtype=np.bool)
    badPointMask = np.zeros(data.shape, dtype=np.bool)

    # add bad points
    badIndex = np.random.randint(0, high=data.size-1, size=10)
    for eachIndex in badIndex:
        sigma[eachIndex] = np.inf
        badmask[eachIndex] = 1
        badPointMask[eachIndex] = 1

    # add singular points
    badIndex = np.random.randint(0, high=data.size-1, size=20)
    for eachIndex in badIndex:
        data[eachIndex] += (np.random.random()-0.5)*20

    # add nan points
    badIndex = np.random.randint(0, high=data.size-1, size=5)
    for eachIndex in badIndex:
        data[eachIndex] = np.nan
        sigma[eachIndex] = np.inf
        badmask[eachIndex] = 1

    plt.scatter(range(data.size), data)
    plt.errorbar(range(data.size), data, ls="None", yerr=sigma)
    plt.scatter(np.where(badPointMask)[0], data[badPointMask], color="r")
    #endregion: set data and figures
    #region: tests
    r0 = iterstat1d(data, calMedian=True, maxRej=4)
    print("r0:   {} +- {} @ {} {}".format(r0["mean"], r0["std"], r0["niter"], r0["median"]))
    rt0 = iterstat1dC(data, calMedian=True, maxRej=4)
    print("rt0:  {} +- {} @ {} {}".format(rt0["mean"], rt0["std"], rt0["niter"], rt0["median"]))
    assert np.all(r0["badmask"]==rt0["badmask"])
    assert np.all(r0["newBadMask"]==rt0["newBadMask"])

    r1 = iterstat1d(data, mask=badmask, calMedian=True, maxRej=3)
    print("r1:   {} +- {} @ {} {}".format(r1["mean"], r1["std"], r1["niter"], r1["median"]))
    rt1 = iterstat1dC(data, mask=badmask, calMedian=True, maxRej=3)
    print("rt1:  {} +- {} @ {} {}".format(rt1["mean"], rt1["std"], rt1["niter"], rt1["median"]))
    assert np.all(r1["badmask"]==rt1["badmask"])
    assert np.all(r1["newBadMask"]==rt1["newBadMask"])

    r2 = iterstat1d(data, sigma=sigma, calMedian=True, maxRej=2)
    print("r2:   {} +- {} @ {} {}".format(r2["mean"], r2["std"], r2["niter"], r2["median"]))
    rt2 = iterstat1dC(data, sigma=sigma, calMedian=True, maxRej=2)
    print("rt2:  {} +- {} @ {} {}".format(rt2["mean"], rt2["std"], rt2["niter"], rt2["median"]))
    assert np.all(r2["badmask"]==rt2["badmask"])
    assert np.all(r2["newBadMask"]==rt2["newBadMask"])

    r3 = iterstat1d(data, invvar=1/sigma**2, calMedian=True, maxRej=3, maxDev=2)
    print("r3:   {} +- {} @ {} {}".format(r3["mean"], r3["std"], r3["niter"], r3["median"]))
    rt3 = iterstat1dC(data, invvar=1/sigma**2, calMedian=True, maxRej=3, maxDev=2)
    print("rt3:  {} +- {} @ {} {}".format(rt3["mean"], rt3["std"], rt3["niter"], rt3["median"]))
    assert np.all(r3["badmask"]==rt3["badmask"])
    assert np.all(r3["newBadMask"]==rt3["newBadMask"])

    r4 = iterstat1d(data, sigma=sigma, invvar=1/sigma**2, calMedian=True, maxRej=4)
    print("r4:   {} +- {} @ {} {}".format(r4["mean"], r4["std"], r4["niter"], r4["median"]))
    rt4 = iterstat1dC(data, sigma=sigma, invvar=1/sigma**2, calMedian=True, maxRej=4)
    print("rt4:  {} +- {} @ {} {}".format(rt4["mean"], rt4["std"], rt4["niter"], rt4["median"]))
    assert np.all(r4["badmask"]==rt4["badmask"])
    assert np.all(r4["newBadMask"]==rt4["newBadMask"])

    r5 = iterstat1d(data, invvar=1/sigma**2, maxiter=2, calMedian=True, maxRej=5)
    print("r5:   {} +- {} @ {} {}".format(r5["mean"], r5["std"], r5["niter"], r5["median"]))
    rt5 = iterstat1dC(data, invvar=1/sigma**2, maxiter=2, calMedian=True, maxRej=5)
    print("rt5:  {} +- {} @ {} {}".format(rt5["mean"], rt5["std"], rt5["niter"], rt5["median"]))
    assert np.all(r5["badmask"]==rt5["badmask"])
    assert np.all(r5["newBadMask"]==rt5["newBadMask"])

    r6 = iterstat1d(data, invvar=1/sigma**2, maxiter=2, calMedian=True, stdRej=True, maxRej=1)
    print("r6:   {} +- {} @ {} {}".format(r6["mean"], r6["std"], r6["niter"], r6["median"]))
    rt6 = iterstat1dC(data, invvar=1/sigma**2, maxiter=2, calMedian=True, stdRej=True, maxRej=1)
    print("rt6:  {} +- {} @ {} {}".format(rt6["mean"], rt6["std"], rt6["niter"], rt6["median"]))
    assert np.all(r6["badmask"]==rt6["badmask"])
    assert np.all(r6["newBadMask"]==rt6["newBadMask"])

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
    #endregion: tests
    plt.legend()
    plt.scatter(np.where(r3["newBadMask"])[0], data[r3["newBadMask"]], c="red", s=300, marker="x")
    plt.show(0)
def iterstat2d_test():
    #region: set data and figures
    fig = plt.figure()
    ax = fig.add_subplot(111)
    c = cm.get_cmap("rainbow")
    c = c(np.linspace(0,1,6))
    # test for iterstat2d
    m = 300
    n = 50

    x = np.arange(m)
    y = 2*np.sin(0.141*x) + 3.1*np.cos(0.312*x)
    data2d = []
    badmask2d = []
    sigma2d = []
    badPointMask2d = []
    for i in range(n):
        data = y + np.random.randn(x.size)
        sigma = (np.random.random(x.size)-0.5)/2+1
        badmask = np.zeros(data.shape, dtype=np.bool)
        badPointMask = np.zeros(data.shape, dtype=np.bool)

        # add bad points
        badIndex = np.random.randint(0, high=data.size-1, size=30)
        for eachIndex in badIndex:
            sigma[eachIndex] = np.inf
            badmask[eachIndex] = 1
            badPointMask[eachIndex] = 1

        # add singular points
        badIndex = np.random.randint(0, high=data.size-1, size=60)
        for eachIndex in badIndex:
            data[eachIndex] += (np.random.random()-0.5)*20

        # add nan points
        badIndex = np.random.randint(0, high=data.size-1, size=10)
        for eachIndex in badIndex:
            data[eachIndex] = np.nan
            sigma[eachIndex] = np.inf
            badmask[eachIndex] = 1

        data2d.append(data)
        sigma2d.append(sigma)
        badmask2d.append(badmask)
        badPointMask2d.append(badPointMask)

    data2d = np.array(data2d)
    sigma2d = np.array(sigma2d)
    badmask2d = np.array(badmask2d)
    badPointMask2d = np.array(badPointMask2d)

    for data,sigma,badPointMask in zip(data2d, sigma2d,badPointMask2d):
        plt.scatter(range(data.size), data)
        plt.errorbar(range(data.size), data, ls="None", yerr=sigma)
        plt.scatter(np.where(badPointMask)[0], data[badPointMask], color="r")
    #endregion: set data and figures
    #region: tests
    r0 = iterstat2d(data2d, calMedian=True)
    #print("r0:   {} +- {} @ {} {}".format(r0["mean"], r0["std"], r0["niter"], r0["median"]))
    #rt0 = iterstat1dC(data2d, calMedian=True)
    #print("rt0:  {} +- {} @ {} {}".format(rt0["mean"], rt0["std"], rt0["niter"], rt0["median"]))
    #assert np.all(r0["badmask"]==rt0["badmask"])
    #assert np.all(r0["newBadMask"]==rt0["newBadMask"])
    r1 = iterstat2d(data2d, mask=badmask2d, calMedian=True)
    #print("r1:   {} +- {} @ {} {}".format(r1["mean"], r1["std"], r1["niter"], r1["median"]))
    #rt1 = iterstat1dC(data2d, mask=badmask, calMedian=True)
    #print("rt1:  {} +- {} @ {} {}".format(rt1["mean"], rt1["std"], rt1["niter"], rt1["median"]))
    #assert np.all(r1["badmask"]==rt1["badmask"])
    #assert np.all(r1["newBadMask"]==rt1["newBadMask"])
    r2 = iterstat2d(data2d, sigma=sigma2d, calMedian=True)
    #print("r2:   {} +- {} @ {} {}".format(r2["mean"], r2["std"], r2["niter"], r2["median"]))
    #rt2 = iterstat1dC(data2d, sigma=sigma, calMedian=True)
    #print("rt2:  {} +- {} @ {} {}".format(rt2["mean"], rt2["std"], rt2["niter"], rt2["median"]))
    #assert np.all(r2["badmask"]==rt2["badmask"])
    #assert np.all(r2["newBadMask"]==rt2["newBadMask"])
    r3 = iterstat2d(data2d, invvar=1/sigma2d**2, calMedian=True)
    #print("r3:   {} +- {} @ {} {}".format(r3["mean"], r3["std"], r3["niter"], r3["median"]))
    #rt3 = iterstat1dC(data2d, invvar=1/sigma**2, calMedian=True)
    #print("rt3:  {} +- {} @ {} {}".format(rt3["mean"], rt3["std"], rt3["niter"], rt3["median"]))
    #assert np.all(r3["badmask"]==rt3["badmask"])
    #assert np.all(r3["newBadMask"]==rt3["newBadMask"])
    r4 = iterstat2d(data2d, sigma=sigma2d, invvar=1/sigma2d**2, calMedian=True)
    #print("r4:   {} +- {} @ {} {}".format(r4["mean"], r4["std"], r4["niter"], r4["median"]))
    #rt4 = iterstat1dC(data2d, sigma=sigma, invvar=1/sigma**2, calMedian=True)
    #print("rt4:  {} +- {} @ {} {}".format(rt4["mean"], rt4["std"], rt4["niter"], rt4["median"]))
    #assert np.all(r4["badmask"]==rt4["badmask"])
    #assert np.all(r4["newBadMask"]==rt4["newBadMask"])
    r5 = iterstat2d(data2d, invvar=1/sigma2d**2, maxiter=2, calMedian=True)
    #print("r5:   {} +- {} @ {} {}".format(r5["mean"], r5["std"], r5["niter"], r5["median"]))
    #rt5 = iterstat1dC(data2d, invvar=1/sigma**2, maxiter=2, calMedian=True)
    #print("rt5:  {} +- {} @ {} {}".format(rt5["mean"], rt5["std"], rt5["niter"], rt5["median"]))
    #assert np.all(r5["badmask"]==rt5["badmask"])
    #assert np.all(r5["newBadMask"]==rt5["newBadMask"])
    plt.plot(x, r1["mean"], label="r1", c=c[1])
    plt.plot(x, r1["mean"]+r1["std"], "--", label="r1", c=c[1])
    plt.plot(x, r1["mean"]-r1["std"], "--", label="r1", c=c[1])
    plt.plot(x, r2["mean"], label="r2", c=c[2])
    plt.plot(x, r2["mean"]+r2["std"], "--", label="r2", c=c[2])
    plt.plot(x, r2["mean"]-r2["std"], "--", label="r2", c=c[2])
    plt.plot(x, r3["mean"], label="r3", c=c[3])
    plt.plot(x, r3["mean"]+r3["std"], "--", label="r3", c=c[3])
    plt.plot(x, r3["mean"]-r3["std"], "--", label="r3", c=c[3])
    plt.plot(x, r4["mean"], label="r4", c=c[4])
    plt.plot(x, r4["mean"]+r4["std"], "--", label="r4", c=c[4])
    plt.plot(x, r4["mean"]-r4["std"], "--", label="r4", c=c[4])
    plt.plot(x, r5["mean"], label="r5", c=c[5])
    plt.plot(x, r5["mean"]+r5["std"], "--", label="r5", c=c[5])
    plt.plot(x, r5["mean"]-r5["std"], "--", label="r5", c=c[5])
    #endregion: tests
    plt.legend()
    for i in range(n):
        plt.scatter(np.where(r3["newBadMask"][i])[0], data2d[i][r3["newBadMask"][i]], c="red", s=300, marker="x")
    fig.tight_layout()
    plt.show(0)
def iterstat3d_test():
    todoFits = glob.glob("test/sci*")
    data = []
    sigma = []
    for eachFile in todoFits:
        obj = pyfits.open(eachFile)
        data.append(obj["sci"].data)
        sigma.append(obj["std"].data)
    data = np.array(data)
    N = len(data)
    sigma = np.array(sigma)
    result = iterstat3d(data, sigma=sigma, calMedian=True, maxRej=1, stdRej=True)
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
    fits.writeto("test/test_iterstat3d.fits", clobber='True')
def iterstat3dC_test():
    todoFits = glob.glob("test/sci*")
    data = []
    sigma = []
    for eachFile in todoFits:
        obj = pyfits.open(eachFile)
        data.append(obj["sci"].data)
        sigma.append(obj["std"].data)
    data = np.array(data)
    N = len(data)
    sigma = np.array(sigma)
    result = iterstat3dC(data, sigma=sigma, calMedian=True, maxRej=1, stdRej=1, maxDev=700)
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
    fits.writeto("test/iterstat3dC.fits", clobber='True')
def linearInterpWithError_test():
    N=100+1
    x = np.arange(N) + (np.random.random(N)-0.5)
    xx = np.linspace(-100,N, 30000)
    #sigma = np.ones(N) * 2
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


if __name__=="__main__":
    #for i in range(1000):
        #iterstat1d_test()

    #iterstat2d_test()
    #iterstat3d_test()
    #iterstat3dC_test()

    linearInterpWithError_test()






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
