import numpy as np
import time
import unittest
import os
import ctypes
from ctypes import c_int, c_double, c_short, c_long
import ipdb
c_int_p = ctypes.POINTER(c_int)
c_short_p = ctypes.POINTER(c_short)
c_long_p = ctypes.POINTER(c_long)
c_double_p = ctypes.POINTER(c_double)
so = ctypes.CDLL(os.path.join(os.path.dirname(__file__), "lib", "astroMath.so"))
yAx_c = so.yAx

#np.polynomial.Legendre
#np.polynomial.Chebyshev

eps = np.finfo(np.float).eps

def Lnx(order, x):
    """Legendre(order, x), x must be a 1D array"""
    N = len(x)
    assert order>=1
    result  = np.zeros((order, N))
    resultP = np.zeros((order, N))
    result[0,:] = 1
    resultP[0,:] = 0
    if order>=2:
        result[1,:] = x[:]
        resultP[1,:] = 1
    for n in range(2,order):
        result[n,:] = ((2*n-1)*x*result[n-1,:] - (n-1)*result[n-2,:])/n
        resultP[n,:] = ((2*n-1)*(result[n-1,:]+x*resultP[n-1,:]) - (n-1)*resultP[n-2,:])/n
    return result, resultP

def Cnx(order, x):
    """Chebyshev(order, x), x must be a 1D array"""
    N = len(x)
    assert order>=1
    result = np.zeros((order, N))
    resultP = np.zeros((order, N))
    result[0,:] = 1
    resultP[0,:] = 0
    if order>=2:
        result[1,:]  = x[:]
        resultP[1,:] = 1
    for n in range(2,order):
        result[n,:]  = 2*x*result[n-1,:] - result[n-2,:]
        resultP[n,:] = 2*(result[n-1,:]+x*resultP[n-1,:]) - resultP[n-2,:]
    return result, resultP

class Polynomial2D(object):
    def __init__(self, xtype, ytype, coeffs, xlims=None, ylims=None):
        self.xtype = xtype
        self.ytype = ytype
        self.yOrder, self.xOrder = coeffs.shape
        self.coeffs = np.matrix(coeffs)
        self.xlims = xlims
        self.ylims = ylims
        self.xlen = xlims[1]-xlims[0]
        self.ylen = ylims[1]-ylims[0]
        if self.xtype=="L":
            self.Px = lambda _:Lnx(self.xOrder, _)
        elif self.xtype=="C":
            self.Px = lambda _:Cnx(self.xOrder, _)
        else:
            raise ValueError("xtype should be L or C")
        if self.ytype=="L":
            self.Py = lambda _:Lnx(self.yOrder, _)
        elif self.ytype=="C":
            self.Py = lambda _:Cnx(self.yOrder, _)
        else:
            raise ValueError("ytype should be L or C")

    def call_py(self, x, y):
        """x,y should all be 1D array"""
        assert len(x)==len(y)
        x,y = self.mapCoord(x,y)
        self.px = np.matrix(self.Px(x)[0])
        self.py = np.matrix(self.Py(y)[0]).T
        result = self.py*self.coeffs*self.px
        #print()
        #print(self.py.shape, "*", self.coeffs.shape, "*",
        #      self.px.shape, "=", result.shape)
        return np.diag(result)

    def __call__(self, x, y):
        """x,y should all be 1D array, calculate in c"""
        assert len(x)==len(y)
        x,y = self.mapCoord(x,y)
        self.px = (self.Px(x)[0].T).copy()
        self.py = (self.Py(y)[0].T).copy()
        return self.yAx(self.py, self.coeffs.A, self.px)

    def partialx(self, x, y):
        assert len(x)==len(y)
        x,y = self.mapCoord(x,y)
        self.px = (self.Px(x)[1].T).copy()/self.xlen*2
        self.py = (self.Py(y)[0].T).copy()
        return self.yAx(self.py, self.coeffs.A, self.px)

    def partialy(self, x, y):
        assert len(x)==len(y)
        x,y = self.mapCoord(x,y)
        self.px = (self.Px(x)[0].T).copy()
        self.py = (self.Py(y)[1].T).copy()/self.ylen*2
        return self.yAx(self.py, self.coeffs.A, self.px)

    def yAx(self, y, A, x):
        N = y.shape[0]
        y = y.astype(np.float64)
        x = x.astype(np.float64)
        A = A.astype(np.float64)
        y_c = np.ctypeslib.as_ctypes(y)
        x_c = np.ctypeslib.as_ctypes(x)
        A_c = np.ctypeslib.as_ctypes(A)
        result = np.zeros(N, dtype=np.float64)
        result_c = np.ctypeslib.as_ctypes(result)
        yAx_c.argtypes = [
                y_c.__class__, A_c.__class__, x_c.__class__,
                c_int, c_int, c_int,
                result_c.__class__
            ]
        yAx_c(y_c, A_c, x_c,
              N, self.xOrder, self.yOrder,
              result_c)
        #for i in range(N):
        #    for j in range(self.yOrder):
        #        assert y[i,j]==y_c[i][j], (i,j)
        #for i in range(N):
        #    for j in range(self.xOrder):
        #        assert x[i,j]==x_c[i][j], (i,j)
        return result

    def mapCoord(self, x, y):
        """map x and y data into [-1, 1]
           if self.xlims and self.ylims exists
        """
        if self.xlims is None:
            if not np.all((-1<x)&(x<1)):
                raise ValueError("no xlims and x not in [-1,1]")
        else:
            x = (x-self.xlims[0])/(self.xlims[1]-self.xlims[0])*2-1
        if self.ylims is None:
            if not np.all((-1<y)&(y<1)):
                raise ValueError("no ylims and y not in [-1,1]")
        else:
            y = (y-self.ylims[0])/(self.ylims[1]-self.ylims[0])*2-1
        return x,y

class Legendre2D(Polynomial2D):
    def __init__(self, coeffs, xlims=None, ylims=None):
        super().__init__("L", "L", coeffs, xlims, ylims)
class Chebyshev2D(Polynomial2D):
    def __init__(self, coeffs, xlims=None, ylims=None):
        super().__init__("C", "C", coeffs, xlims, ylims)


if __name__ == '__main__':
    testN = 9999
    tol=5000
    class TestPolynomial(unittest.TestCase):
        def setUp(self):
            self.starttime = None
            self.midtime = None
            self.endtime = None
        def tearDown(self):
            if self.starttime is not None:
                print("{}s {}s".format(self.midtime-self.starttime, self.endtime-self.midtime))
        def test_Lnx(self):
            x = np.random.uniform(-1,1,testN*10)
            N = 10
            result = np.zeros((N,len(x)))
            resultP = np.zeros((N,len(x)))
            testResult,testResultP = Lnx(N,x)
            for i in range(N):
                result[i,:] = np.polynomial.legendre.Legendre.basis(i)(x)
                resultP[i,:] = np.polynomial.Legendre([0]*i+[1]).deriv()(x)

            diff = np.abs(result-testResult)
            max = diff.max()
            maxarg = np.argmax(max)
            maxDiffValue = diff.flatten()[maxarg]
            maxValue = result.flatten()[maxarg]
            if not maxValue:
                maxValue = testResult.flatten()[maxarg]
            self.assertTrue(diff.max()<eps*tol, msg="{}>{}, relative:{}/{} count:{}/{}".
                    format(diff.max(), eps*tol, maxDiffValue,maxValue,(diff>eps*tol).sum(), testN))
            diff = np.abs(resultP-testResultP)
            max = diff.max()
            maxarg = np.argmax(max)
            maxDiffValue = diff.flatten()[maxarg]
            maxValue = resultP.flatten()[maxarg]
            if not maxValue:
                maxValue = testResultP.flatten()[maxarg]
            self.assertTrue(diff.max()<eps*tol, msg="{}>{}, relative:{}/{} count:{}/{}".
                    format(diff.max(), eps*tol, maxDiffValue,maxValue,(diff>eps*tol).sum(), testN))
            #diff = np.abs(result-testResult)
            #self.assertTrue(diff.max()<eps*tol, msg="{}>{}".format(diff.max(), eps*tol))
            #diff = np.abs(resultP-testResultP)
            #self.assertTrue(diff.max()<eps*tol, msg="{}>{}".format(diff.max(), eps*tol))
        def test_Cnx(self):
            x = np.random.uniform(-1,1,testN*10)
            N = 10
            result = np.zeros((N,len(x)))
            resultP = np.zeros((N,len(x)))
            testResult,testResultP = Cnx(N,x)
            for i in range(N):
                result[i,:] = np.polynomial.chebyshev.Chebyshev.basis(i)(x)
                resultP[i,:] = np.polynomial.Chebyshev([0]*i+[1]).deriv()(x)
            diff = np.abs(result-testResult)
            max = diff.max()
            maxarg = np.argmax(max)
            maxDiffValue = diff.flatten()[maxarg]
            maxValue = result.flatten()[maxarg]
            if not maxValue:
                maxValue = testResult.flatten()[maxarg]
            self.assertTrue(diff.max()<eps*tol, msg="{}>{}, relative:{}/{} count:{}/{}".
                    format(diff.max(), eps*tol, maxDiffValue,maxValue,(diff>eps*tol).sum(), testN))

            diff = np.abs(resultP-testResultP)
            max = diff.max()
            maxarg = np.argmax(max)
            maxDiffValue = diff.flatten()[maxarg]
            maxValue = resultP.flatten()[maxarg]
            if not maxValue:
                maxValue = testResultP.flatten()[maxarg]
            self.assertTrue(diff.max()<eps*tol, msg="{}>{}, relative:{}/{} count:{}/{}".
                    format(diff.max(), eps*tol, maxDiffValue, maxValue,(diff>eps*tol).sum(), testN))
            #diff = np.abs(result-testResult)
            #self.assertTrue(diff.max()<eps*tol, msg="{}>{}".format(diff.max(), eps*tol))
            #diff = np.abs(resultP-testResultP)
            #self.assertTrue(diff.max()<eps*tol, msg="{}>{}".format(diff.max(), eps*tol))
        def test_Legendre2D(self):
            self.starttime=time.time()
            coeffs = np.random.randn(6,8)*100
            x = np.random.uniform(-1,1,testN)
            y = np.random.uniform(-1,1,testN)
            l = Legendre2D(coeffs)
            result_py = l.call_py(x,y)
            self.midtime=time.time()
            result = l(x,y)
            self.endtime=time.time()
            diff = np.abs(result-result_py)
            max = diff.max()
            maxarg = np.argmax(max)
            maxDiffValue = diff.flatten()[maxarg]
            maxValue = result.flatten()[maxarg]
            self.assertTrue(diff.max()<eps*tol, msg="{}>{}, relative:{} count:{}/{}".
                    format(diff.max(), eps*tol, maxDiffValue/maxValue,(diff>eps*tol).sum(), testN))
        def test_Legendre2D_2(self):
            coeffs = np.random.randn(6,8)*100
            x = np.random.uniform(-1000,1000,testN)
            y = np.random.uniform(-2000,2000,testN)
            l = Legendre2D(coeffs, xlims=[-1000,1000], ylims=[-2000,2000])
            result_py = l.call_py(x,y)
            result = l(x,y)
            diff = np.abs(result-result_py)
            max = diff.max()
            maxarg = np.argmax(max)
            maxDiffValue = diff.flatten()[maxarg]
            maxValue = result.flatten()[maxarg]
            self.assertTrue(diff.max()<eps*tol, msg="{}>{}, relative:{} count:{}/{}".
                    format(diff.max(), eps*tol, maxDiffValue/maxValue,(diff>eps*tol).sum(), testN))
        def test_Chebyshev2D(self):
            self.starttime=time.time()
            coeffs = np.random.randn(6,8)*100
            x = np.random.uniform(-1,1,testN)
            y = np.random.uniform(-1,1,testN)
            l = Chebyshev2D(coeffs)
            result_py = l.call_py(x,y)
            self.midtime=time.time()
            result = l(x,y)
            self.endtime=time.time()
            #diff = np.abs(result-result_py)
            #self.assertTrue(diff.max()<eps*tol, msg="{}>{}, count:{}/{}".
            #        format(diff.max(), eps*tol, (diff>eps*tol).sum(), testN))
            diff = np.abs(result-result_py)
            max = diff.max()
            maxarg = np.argmax(max)
            maxDiffValue = diff.flatten()[maxarg]
            maxValue = result.flatten()[maxarg]
            self.assertTrue(diff.max()<eps*tol, msg="{}>{}, relative:{} count:{}/{}".
                    format(diff.max(), eps*tol, maxDiffValue/maxValue,(diff>eps*tol).sum(), testN))
        def test_Chebyshev2D_2(self):
            coeffs = np.random.randn(6,8)*100
            x = np.random.uniform(-1000,1000,testN)
            y = np.random.uniform(-2000,2000,testN)
            l = Chebyshev2D(coeffs, xlims=[-1000,1000], ylims=[-2000,2000])
            result_py = l.call_py(x,y)
            result = l(x,y)
            diff = np.abs(result-result_py)
            max = diff.max()
            maxarg = np.argmax(max)
            maxDiffValue = diff.flatten()[maxarg]
            maxValue = result.flatten()[maxarg]
            self.assertTrue(diff.max()<eps*tol, msg="{}>{}, relative:{} count:{}/{}".
                    format(diff.max(), eps*tol, maxDiffValue/maxValue,(diff>eps*tol).sum(), testN))
        def test_ChebyshevPartial(self):
            coeffs = np.array([[0,2],
                               [4,0]])
            x = np.random.uniform(-1,1,testN)
            y = np.random.uniform(-1,1,testN)
            l = Chebyshev2D(coeffs)
            resultx = l.partialx(x,y)
            resulty = l.partialy(x,y)
            self.assertTrue(np.all(resultx==2), msg="")
            self.assertTrue(np.all(resulty==4), msg="")

    suite = unittest.TestLoader().loadTestsFromTestCase(TestPolynomial)
    unittest.TextTestRunner(verbosity=2).run(suite)
