from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline,InterpolatedUnivariateSpline,LSQUnivariateSpline
from scipy import interpolate
import ipdb

N=100+1
xx = np.arange(N) + (np.random.random(N)-0.5)
xxx = np.linspace(-100,N, 30000)
w = np.ones(N) * 2
f = lambda x:10*np.sin(0.02*x+np.random.random(1)) + 6*np.cos(0.03*x+np.random.random(1)) + 0.001*(x-500)*2
yy = f(xx)

yk1 = InterpolatedUnivariateSpline(xx, yy, w=w, k=1)
yk2 = InterpolatedUnivariateSpline(xx, yy, w=w, k=2)
yk31 = InterpolatedUnivariateSpline(xx, yy, w=w, k=3)
yk3 = LSQUnivariateSpline(xx, yy, xx[2:-2], w=w, k=3)
yk311 = LSQUnivariateSpline(xx, yy, xx[1:-1], w=w, k=3)

def pp(xxx, obj, color="red", s=40):
    print("=======================")
    knots = obj.get_knots()
    xxx = np.sort(np.concatenate((xxx,knots)))
    plt.scatter(knots, obj(knots), color=color, s=s )
    plt.plot(xxx, obj(xxx), color=color)

def p123():
    pp(xxx, yk3, color="blue", s=150)
    pp(xxx, yk2, color="green", s=80)
    pp(xxx, yk1, color="red", s=20)
    plt.plot(xx,yy, "k.", markerSize=3)
    plt.show(0)

def testCoeff(obj):
    fig = plt.figure()
    print("diff:")
    print(xx.shape)
    #print(np.diff(yy)/np.diff(xx))
    pp(xxx, obj, color="red", s=20)
    print("nodes:")
    print(obj.get_knots().shape)
    print("coeff:")
    print(obj.get_coeffs().shape)
    print(obj.get_residual())
    plt.plot(xx,yy, "k.", markerSize=3)
    plt.show(0)

# see https://en.wikipedia.org/wiki/Spline_interpolation
def testWiki():
    yk3 = UnivariateSpline(xx, yy, k=2, s=0)
    knots = yk3.get_knots()
    coeffs = yk3.get_coeffs()
    #ipdb.set_trace()
    r1 = yk3(xxx)
    r2 = gety(xxx, xx, yy, coeffs)
    plt.plot(xx,yy, "k.", markerSize=10)
    plt.plot(xxx, r1, "g.-", markerSize=1)
    plt.plot(xxx, r2, "r.-", markerSize=1)
    plt.show(0)

# len(x) = N = n + 1
# point number: 0~n + 1
def gety(X, x, y, k):
    #!! np.searchsorted([5,6,7,8,9], 7, side="right") = 3
    n = x.size-1
    i = np.searchsorted(x, X, side='right')
    rightOut = (X > x[n])
    leftOut  = (i==0)
    i[rightOut] = n
    t = (X-x[i-1])/(x[i]-x[i-1])
    a = k[i-1]*(x[i]-x[i-1]) - (y[i]-y[i-1])
    b = -k[i]*(x[i]-x[i-1]) + (y[i]-y[i-1])
    q = (1-t)*y[i-1] + t*y[i] + t*(1-t)*(a*(1-t)+b*t)
    #leftoutResult = k[0]*(x[1]-x[0])*(X-x[0])
    #rightoutResult= k[n]*(x[n]-x[n-1])*(X-x[n])
    #q[leftOut] = leftoutResult[leftOut]
    #q[rightOut] = rightoutResult[rightOut]
    return q

def toStr(l):
    return [str(each) for each in l]

def hackinUnivariateSpline(obj):
    print("data len: {}".format(xx.size))
    print("x")
    print(" ".join(toStr(xx.tolist())))
    print("y")
    print(" ".join(toStr(yy.tolist())))
    print("len knots: {}".format(obj.get_knots().size))
    print(" ".join(toStr(obj.get_knots().tolist())))
    print("len coeffs: {}".format(obj.get_coeffs().size))
    print(" ".join(toStr(obj.get_coeffs().tolist())))
    data = obj._data
    xb = data[3]
    xe = data[4]
    x = data[0]
    y = data[1]
    k = data[5]
    s = data[6]
    n = data[7]
    t = data[8]
    c = data[9]
    print("len data: ", len(data))
    print("x[0], x[1], x[-2], x[-1]:\n\t{} {} {} {}".format(x[0],x[1],x[-2],x[-1]))
    print("xb:{}, xe:{}".format(xb, xe))
    print("k:{}, s:{}".format(k, s))
    print("n(?)")
    print(n)
    print("t(?)")
    print(t)
    print("c(?)")
    print(c)
    tck = interpolate.splrep(x,y,s=0)
    ipdb.set_trace()

hackinUnivariateSpline(yk3)

#testCoeff(yk1)
#testCoeff(yk2)
#testCoeff(yk3)
#p123()
#testWiki()
