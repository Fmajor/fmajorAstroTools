import numpy as np
import tempfile
import time
import unittest
import copy
import os
import ctypes
import astropy.io.fits as fits
from ctypes import c_int, c_double, c_short, c_long, POINTER, Structure
import io
from io import StringIO
import sys
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.collections import PatchCollection
import itertools
c_int_p    = POINTER(c_int)
c_short_p  = POINTER(c_short)
c_long_p   = POINTER(c_long)
c_double_p = POINTER(c_double)
#p = ctypes.create_string_buffer("Hello", 10)
null_ptr = POINTER(c_int)()
_lib = ctypes.CDLL(os.path.join(os.path.dirname(__file__), "lib", "astroMath.so"))

eps = np.finfo(np.float).eps
#.argtypes =
#.restype =
#y = y.astype(np.float64)
#y_c = np.ctypeslib.as_ctypes(y)

#region: structure and functions
class struct_Vertex(Structure):
    pass
struct_Vertex.__slots__ = [
    'x',
    'y',
    'prev',
    'next',
    'flag',
]
struct_Vertex._fields_ = [
    ('x', c_double),
    ('y', c_double),
    ('prev', POINTER(struct_Vertex)),
    ('next', POINTER(struct_Vertex)),
    ('flag', c_int),
]
class struct_Polygon(Structure):
    pass
struct_Polygon.__slots__ = [
    'P',
    'num',
]
struct_Polygon._fields_ = [
    ('P', POINTER(struct_Vertex)),
    ('num', c_int),
]
class struct_Vector2(Structure):
    pass
struct_Vector2.__slots__ = [
    'x',
    'y',
]
struct_Vector2._fields_ = [
    ('x', c_double),
    ('y', c_double),
]
class struct_Vector3(Structure):
    pass
struct_Vector3.__slots__ = [
    'x',
    'y',
    'z',
]
struct_Vector3._fields_ = [
    ('x', c_double),
    ('y', c_double),
    ('z', c_double),
]
struct_Vertex_p = POINTER(struct_Vertex)
struct_Polygon_p = POINTER(struct_Polygon)
struct_Vector2_p = POINTER(struct_Vector2)
struct_Vector3_p = POINTER(struct_Vector3)


TriangleTrianglePolygon_c = _lib.TriangleTrianglePolygon
TriangleTrianglePolygon_c.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double]
TriangleTrianglePolygon_c.restype = POINTER(struct_Polygon)

TriangleSquarePolygon_c = _lib.TriangleSquarePolygon
TriangleSquarePolygon_c.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double]
TriangleSquarePolygon_c.restype = POINTER(struct_Polygon)

areaTriangleSquare_c = _lib.areaTriangleSquare
areaTriangleSquare_c.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, POINTER(c_double), POINTER(c_double), POINTER(c_double)]
areaTriangleSquare_c.restype = None
areaPolygon_c = _lib.areaPolygon
areaPolygon_c.argtypes = [POINTER(struct_Polygon)]
areaPolygon_c.restype = c_double

newPolygon_c = _lib.newPolygon
newPolygon_c.argtypes = []
newPolygon_c.restype = POINTER(struct_Polygon)

#!! struct Vertex *addVertex(struct Polygon *obj, double x, double y, int flag, int index)
addVertex_c = _lib.addVertex
addVertex_c.argtypes = [POINTER(struct_Polygon), c_double, c_double, c_int, c_int]
addVertex_c.restype = POINTER(struct_Vertex)

#!! struct Vertex *delVertex(struct Polygon *obj,int index)
delVertex_c = _lib.delVertex
delVertex_c.argtypes = [POINTER(struct_Polygon), c_int]
delVertex_c.restype = POINTER(struct_Vertex)

#!! int showPolygon(struct Polygon *obj)
showPolygon_c = _lib.showPolygon
showPolygon_c.argtypes = [POINTER(struct_Polygon)]
showPolygon_c.restype = c_int

#!! struct Vertex* newVertex(double x, double y, int flag)
newVertex_c = _lib.newVertex
newVertex_c.argtypes = [c_double, c_double, c_int]
newVertex_c.restype = POINTER(struct_Vertex)

#!! struct Vertex* intersect(struct Vertex *A, struct Vertex *B
#!!                          struct Vertex *C, struct Vertex *D)
intersect_c = _lib.intersect
intersect_c.argtypes = [POINTER(struct_Vertex), POINTER(struct_Vertex),
                        POINTER(struct_Vertex), POINTER(struct_Vertex)]
intersect_c.restype = POINTER(struct_Vertex)


cross3_c = _lib.cross3
cross3_c.argtypes = [POINTER(struct_Vector3), POINTER(struct_Vector3)]
cross3_c.restype = POINTER(struct_Vector3)
norm3_c = _lib.norm3
norm3_c.argtypes = [POINTER(struct_Vector3)]
norm3_c.restype = c_double
dot3_c = _lib.dot3
dot3_c.argtypes = [POINTER(struct_Vector3), POINTER(struct_Vector3)]
dot3_c.restype = c_double
newVector3_c = _lib.newVector3
newVector3_c.argtypes = [c_double, c_double, c_double]
newVector3_c.restype = POINTER(struct_Vector3)
neg3_c = _lib.neg3
neg3_c.argtypes = [POINTER(struct_Vector3)]
neg3_c.restype = POINTER(struct_Vector3)
minus3_c = _lib.minus3
minus3_c.argtypes = [POINTER(struct_Vector3), POINTER(struct_Vector3)]
minus3_c.restype = POINTER(struct_Vector3)
plus3_c = _lib.plus3
plus3_c.argtypes = [POINTER(struct_Vector3), POINTER(struct_Vector3)]
plus3_c.restype = POINTER(struct_Vector3)

norm2_c = _lib.norm2
norm2_c.argtypes = [POINTER(struct_Vector2)]
norm2_c.restype = c_double
dot2_c = _lib.dot2
dot2_c.argtypes = [POINTER(struct_Vector2), POINTER(struct_Vector2)]
dot2_c.restype = c_double
newVector2_c = _lib.newVector2
newVector2_c.argtypes = [c_double, c_double]
newVector2_c.restype = POINTER(struct_Vector2)
neg2_c = _lib.neg2
neg2_c.argtypes = [POINTER(struct_Vector2)]
neg2_c.restype = POINTER(struct_Vector2)
minus2_c = _lib.minus2
minus2_c.argtypes = [POINTER(struct_Vector2), POINTER(struct_Vector2)]
minus2_c.restype = POINTER(struct_Vector2)
plus2_c = _lib.plus2
plus2_c.argtypes = [POINTER(struct_Vector2), POINTER(struct_Vector2)]
plus2_c.restype = POINTER(struct_Vector2)

get4map_c = _lib.get4map
get4map_c.restype = c_int

searchsorted_c = _lib.searchsorted
searchsorted_c.argtypes = [c_int, c_double_p, c_double]
searchsorted_c.restype = c_int

interpolate2D_c = _lib.interpolate2D
interpolate1D_c = _lib.interpolate1D
#endregion: structure and functions

#region: py functions
def searchsorted(array, value):
    array = array.astype(np.float64)
    n = array.size
    array_c = np.ctypeslib.as_ctypes(array)
    return searchsorted_c(n, array_c, value)

def newPolygon():
    return newPolygon_c()
def addVertex(obj, x, y, flag, index):
    return addVertex_c(obj, x, y, flag, index)
def delVertex(obj, index):
    return delVertex_c(obj, index)
def showPolygon(obj):
    return showPolygon_c(obj)
def newVertex(x, y, flag):
    return newVertex_c(x, y, flag)
def intersect(a, b, c, d):
    return intersect_c(a, b, c, d)

def cross3(a,b):
    return cross3_c(a,b)
def norm3(v):
    return norm3_c(v)
def dot3(a,b):
    return dot3_c(a,b)
def newVector3(x,y,z):
    return newVector3_c(x,y,z)
def plus3(a,b):
    return plus3_c(a,b)
def minus3(a,b):
    return minus3_c(a,b)
def neg3(v):
    return neg3_c(v)

def norm2(v):
    return norm2_c(v)
def dot2(a,b):
    return dot2_c(a,b)
def newVector2(x,y):
    return newVector2_c(x,y)
def plus2(a,b):
    return plus2_c(a,b)
def minus2(a,b):
    return minus2_c(a,b)
def neg2(v):
    return neg2_c(v)

def TriangleTrianglePolygon(Ax, Ay, Bx, By, Cx, Cy, Dx, Dy, Ex, Ey, Fx, Fy):
    return TriangleTrianglePolygon_c(Ax, Ay, Bx, By, Cx, Cy, Dx, Dy, Ex, Ey, Fx, Fy)

def TriangleSquarePolygon(x1, y1, x2, y2, x3, y3, x4, y4, Dx, Dy, Ex, Ey, Fx, Fy):
    return TriangleSquarePolygon_c(x1, y1, x2, y2, x3, y3, x4, y4, Dx, Dy, Ex, Ey, Fx, Fy)

def areaTriangleSquare(x1, y1, x2, y2, x3, y3, x4, y4, Dx, Dy, Ex, Ey, Fx, Fy):
    fullArea    = c_double()
    overlapArea = c_double()
    squareArea = c_double()
    areaTriangleSquare_c(x1, y1, x2, y2, x3, y3, x4, y4, Dx, Dy, Ex, Ey, Fx, Fy,
                       ctypes.byref(fullArea), ctypes.byref(overlapArea),
                       ctypes.byref(squareArea))
    return fullArea.value, overlapArea.value

def get4map(xx, yy, mask):
    xx = xx.astype(np.float64)
    yy = yy.astype(np.float64)
    mask = (mask>0).astype(np.int64)
    assert xx.shape==yy.shape==mask.shape
    #print("mask:=================")
    #print(mask)
    #print("mask:=================")
    shape = xx.shape
    ny, nx = shape
    xx_c = np.ctypeslib.as_ctypes(xx)
    yy_c = np.ctypeslib.as_ctypes(yy)
    mask_c = np.ctypeslib.as_ctypes(mask)
    ULx = np.zeros(shape, dtype=np.float64) * np.nan
    ULy = np.zeros(shape, dtype=np.float64) * np.nan
    URx = np.zeros(shape, dtype=np.float64) * np.nan
    URy = np.zeros(shape, dtype=np.float64) * np.nan
    LLx = np.zeros(shape, dtype=np.float64) * np.nan
    LLy = np.zeros(shape, dtype=np.float64) * np.nan
    LRx = np.zeros(shape, dtype=np.float64) * np.nan
    LRy = np.zeros(shape, dtype=np.float64) * np.nan
    ULx_c = np.ctypeslib.as_ctypes(ULx)
    ULy_c = np.ctypeslib.as_ctypes(ULy)
    URx_c = np.ctypeslib.as_ctypes(URx)
    URy_c = np.ctypeslib.as_ctypes(URy)
    LLx_c = np.ctypeslib.as_ctypes(LLx)
    LLy_c = np.ctypeslib.as_ctypes(LLy)
    LRx_c = np.ctypeslib.as_ctypes(LRx)
    LRy_c = np.ctypeslib.as_ctypes(LRy)
    #get4map_c.argtypes = [c_int, c_int, POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double)]
    get4map_c(nx, ny,
              mask_c,
              xx_c, yy_c,
              ULx_c, ULy_c,
              URx_c, URy_c,
              LLx_c, LLy_c,
              LRx_c, LRy_c,
        )
    return ULx, ULy, URx, URy, LLx, LLy, LRx, LRy

def get4map_plot(xx, yy, mask, doprint=False):
    if doprint:
        print(xx[0:10,0:10]); print(yy[0:10,0:10]);print(mask[0:10,0:10])
    print(xx.shape, yy.shape, mask.shape)
    ULx, ULy, URx, URy, LLx, LLy, LRx, LRy = get4map(xx,yy,mask)
    if doprint:
        print(ULx[0:10,0:10]);print(ULy[0:10,0:10]);
        print(URx[0:10,0:10]);print(URy[0:10,0:10]);
        print(LLx[0:10,0:10]);print(LLy[0:10,0:10]);
        print(LRx[0:10,0:10]);print(LRy[0:10,0:10]);
    meandx = np.nanmedian(np.diff(ULx))
    meandy = np.nanmedian(np.diff(ULy, axis=0))
    meand = np.sqrt(meandx**2+meandy**2)
    if doprint:
        print(meandx, meandy, meand)

    import matplotlib as mpl
    mpl.rcParams['figure.subplot.bottom']= 0.05
    mpl.rcParams['figure.subplot.hspace']= 0.05
    mpl.rcParams['figure.subplot.left']  = 0.04
    mpl.rcParams['figure.subplot.right'] = 1
    mpl.rcParams['figure.subplot.top']   = 1
    fig = plt.figure()
    ax = fig.add_subplot(111)
    xypairs = [(ULx, ULy), (URx, URy), (LLx, LLy), (LRx, LRy)]
    dirs  = [[-1,1], [1,1], [-1,-1], [1,-1]]
    texts = ['1', '2', '3', '4']
    ax.scatter(xx[mask], yy[mask], s=100, color='black')
    ny, nx = xx.shape
    for j in range(ny):
        for i in range(nx):
            if mask[j,i]:
                if j<ny-1 and mask[j+1, i]:
                    ax.plot([xx[j,i], xx[j+1,i]],
                            [yy[j,i], yy[j+1,i]],
                            color="black")
                if i<nx-1 and mask[j, i+1]:
                    ax.plot([xx[j,i], xx[j,i+1]],
                            [yy[j,i], yy[j,i+1]],
                            color="black")
    ax.scatter(xx[~mask], yy[~mask], s=50, color='black', marker='x')
    ax.scatter(ULx, ULy, s=30, edgecolor="None", color="red")
    ax.scatter(LLx, LLy, s=30, edgecolor="None", color="green")
    ax.scatter(LRx, LRy, s=30, edgecolor="None", color="blue")
    ax.scatter(URx, URy, s=30, edgecolor="None", color="pink")
    for i in range(4):
        thisxx, thisyy = xypairs[i]
        thisdx, thisdy = dirs[i]
        eachtext = texts[i]
        thisdx *= meandx/10
        thisdy *= meandy/10
        for eachx, eachy in zip(thisxx.flatten(), thisyy.flatten()):
            ax.arrow(eachx, eachy, thisdx, thisdy,
                      head_width = 0.03*meand,
                      head_length = 0.1*meand,
                      )
            #ax.text(eachx, eachy, eachtext)
    plt.show(0)

def get4map_plot_large(xx, yy, mask):
    ULx, ULy, URx, URy, LLx, LLy, LRx, LRy = get4map(xx,yy,mask)
    ULdx = ULx - xx
    ULdy = ULy - yy
    URdx = URx - xx
    URdy = URy - yy
    LLdx = LLx - xx
    LLdy = LLy - yy
    LRdx = LRx - xx
    LRdy = LRy - yy
    from fmajorAstroUtils.DS10 import ds9
    d = ds9.ds10('test')
    d.newFrame
    d.set_np2arr(xx)
    d.newFrame
    d.set_np2arr(yy)
    d.newFrame
    d.set_np2arr(ULdx)
    d.newFrame
    d.set_np2arr(ULdy)
    d.newFrame
    d.set_np2arr(URdx)
    d.newFrame
    d.set_np2arr(URdy)
    d.newFrame
    d.set_np2arr(LLdx)
    d.newFrame
    d.set_np2arr(LLdy)
    d.newFrame
    d.set_np2arr(LRdx)
    d.newFrame
    d.set_np2arr(LRdy)

def interpolate2D(inxx, inyy, inmask, influx, invar, outx, outy):
    assert inxx.shape==inyy.shape==inmask.shape==influx.shape==invar.shape
    inny, innx   = inxx.shape

    inxx = inxx.astype(np.float64)
    inyy = inyy.astype(np.float64)
    influx = influx.astype(np.float64)
    invar  = invar.astype(np.float64)
    inmask = (inmask>0).astype(np.int64)
    outx = outx.astype(np.float64)
    outy = outy.astype(np.float64)
    outxx, outyy = np.meshgrid(outx, outy)
    outxx = outxx.astype(np.float64)
    outyy = outyy.astype(np.float64)
    outny, outnx = outxx.shape

    inshape = inxx.shape
    inoutmask = np.zeros(inshape, dtype=np.float64)
    inoutmask_c = np.ctypeslib.as_ctypes(inoutmask)
    inULx = np.zeros(inshape, dtype=np.float64) * np.nan
    inULy = np.zeros(inshape, dtype=np.float64) * np.nan
    inURx = np.zeros(inshape, dtype=np.float64) * np.nan
    inURy = np.zeros(inshape, dtype=np.float64) * np.nan
    inLLx = np.zeros(inshape, dtype=np.float64) * np.nan
    inLLy = np.zeros(inshape, dtype=np.float64) * np.nan
    inLRx = np.zeros(inshape, dtype=np.float64) * np.nan
    inLRy = np.zeros(inshape, dtype=np.float64) * np.nan
    inULx_c = np.ctypeslib.as_ctypes(inULx)
    inULy_c = np.ctypeslib.as_ctypes(inULy)
    inURx_c = np.ctypeslib.as_ctypes(inURx)
    inURy_c = np.ctypeslib.as_ctypes(inURy)
    inLLx_c = np.ctypeslib.as_ctypes(inLLx)
    inLLy_c = np.ctypeslib.as_ctypes(inLLy)
    inLRx_c = np.ctypeslib.as_ctypes(inLRx)
    inLRy_c = np.ctypeslib.as_ctypes(inLRy)

    inxx_c = np.ctypeslib.as_ctypes(inxx)
    inyy_c = np.ctypeslib.as_ctypes(inyy)

    influx_c = np.ctypeslib.as_ctypes(influx)
    invar_c = np.ctypeslib.as_ctypes(invar)
    inmask_c = np.ctypeslib.as_ctypes(inmask)

    outshape = outxx.shape
    outULx = np.zeros(outshape, dtype=np.float64) * np.nan
    outULy = np.zeros(outshape, dtype=np.float64) * np.nan
    outURx = np.zeros(outshape, dtype=np.float64) * np.nan
    outURy = np.zeros(outshape, dtype=np.float64) * np.nan
    outLLx = np.zeros(outshape, dtype=np.float64) * np.nan
    outLLy = np.zeros(outshape, dtype=np.float64) * np.nan
    outLRx = np.zeros(outshape, dtype=np.float64) * np.nan
    outLRy = np.zeros(outshape, dtype=np.float64) * np.nan
    outULx_c = np.ctypeslib.as_ctypes(outULx)
    outULy_c = np.ctypeslib.as_ctypes(outULy)
    outURx_c = np.ctypeslib.as_ctypes(outURx)
    outURy_c = np.ctypeslib.as_ctypes(outURy)
    outLLx_c = np.ctypeslib.as_ctypes(outLLx)
    outLLy_c = np.ctypeslib.as_ctypes(outLLy)
    outLRx_c = np.ctypeslib.as_ctypes(outLRx)
    outLRy_c = np.ctypeslib.as_ctypes(outLRy)

    outxx_c = np.ctypeslib.as_ctypes(outxx)
    outyy_c = np.ctypeslib.as_ctypes(outyy)
    outflux = np.zeros(outshape, dtype=np.float64)
    outvar  = np.zeros(outshape, dtype=np.float64)
    outmask = np.zeros(outshape, dtype=np.float64)
    outflux_c = np.ctypeslib.as_ctypes(outflux)
    outvar_c  = np.ctypeslib.as_ctypes(outvar)
    outmask_c  = np.ctypeslib.as_ctypes(outmask)
    outx_c  = np.ctypeslib.as_ctypes(outx)
    outy_c  = np.ctypeslib.as_ctypes(outy)

    success = interpolate2D_c(
                    inshape[1], inshape[0],
                    inxx_c,   inyy_c,
                    influx_c, invar_c,
                    inmask_c, inoutmask_c,
                    inULx_c,  inULy_c,
                    inURx_c,  inURy_c,
                    inLLx_c,  inLLy_c,
                    inLRx_c,  inLRy_c,
                    outshape[1], outshape[0],
                    outxx_c, outyy_c,
                    outflux_c, outvar_c,
                    outULx_c,  outULy_c,
                    outURx_c,  outURy_c,
                    outLLx_c,  outLLy_c,
                    outLRx_c,  outLRy_c,
                    outmask_c)
    if success<0:
        raise RuntimeError("have error in interpolate2D")
    debug2D=0
    if debug2D:
        from fmajorAstroUtils.DS10 import ds9
        d = ds9.ds10('test')
        d.newFrame
        d.set_np2arr(inULx)
        d.newFrame
        d.set_np2arr(inULy)
        d.newFrame
        d.set_np2arr(inURx)
        d.newFrame
        d.set_np2arr(inURy)
        d.newFrame
        d.set_np2arr(inLRx)
        d.newFrame
        d.set_np2arr(inLRy)
        d.newFrame
        d.set_np2arr(inLLx)
        d.newFrame
        d.set_np2arr(inLLy)
        d.newFrame

        d.set_np2arr(outULx)
        d.newFrame
        d.set_np2arr(outULy)
        d.newFrame
        d.set_np2arr(outURx)
        d.newFrame
        d.set_np2arr(outURy)
        d.newFrame
        d.set_np2arr(outLRx)
        d.newFrame
        d.set_np2arr(outLRy)
        d.newFrame
        d.set_np2arr(outLLx)
        d.newFrame
        d.set_np2arr(outLLy)
    return {"outflux":outflux,
            "outvar":outvar,
            "outmask":outmask,
            "inoutmask":inoutmask}

def interpolate1D(inx, inmask, influx, invar, outx):
    assert inx.shape==inmask.shape==influx.shape==invar.shape
    innx,  = inx.shape

    inx = inx.astype(np.float64)
    influx = influx.astype(np.float64)
    invar  = invar.astype(np.float64)
    inmask = (inmask>0).astype(np.int64)
    outx = outx.astype(np.float64)
    outnx, = outx.shape

    inshape = inx.shape
    inoutmask = np.zeros(inshape, dtype=np.float64)
    inoutmask_c = np.ctypeslib.as_ctypes(inoutmask)
    inLx = np.zeros(inshape, dtype=np.float64) * np.nan
    inRx = np.zeros(inshape, dtype=np.float64) * np.nan
    inLx_c = np.ctypeslib.as_ctypes(inLx)
    inRx_c = np.ctypeslib.as_ctypes(inRx)

    inx_c = np.ctypeslib.as_ctypes(inx)

    influx_c = np.ctypeslib.as_ctypes(influx)
    invar_c = np.ctypeslib.as_ctypes(invar)
    inmask_c = np.ctypeslib.as_ctypes(inmask)

    outshape = outx.shape
    outLx = np.zeros(outshape, dtype=np.float64) * np.nan
    outRx = np.zeros(outshape, dtype=np.float64) * np.nan
    outLx_c = np.ctypeslib.as_ctypes(outLx)
    outRx_c = np.ctypeslib.as_ctypes(outRx)

    outx_c = np.ctypeslib.as_ctypes(outx)
    outflux = np.zeros(outshape, dtype=np.float64)
    outvar  = np.zeros(outshape, dtype=np.float64)
    outmask = np.zeros(outshape, dtype=np.float64)
    outflux_c = np.ctypeslib.as_ctypes(outflux)
    outvar_c  = np.ctypeslib.as_ctypes(outvar)
    outmask_c  = np.ctypeslib.as_ctypes(outmask)
    outx_c  = np.ctypeslib.as_ctypes(outx)
    success = interpolate1D_c(
                    inshape[0],
                    inx_c,
                    influx_c, invar_c,
                    inmask_c, inoutmask_c,
                    inLx_c,  inRx_c,
                    outshape[0],
                    outx_c,
                    outflux_c, outvar_c,
                    outLx_c,  outRx_c,
                    outmask_c)
    if success<0:
        raise RuntimeError("have error in interpolate2D")
    return {"outflux":outflux,
            "outvar":outvar,
            "outmask":outmask,
            "inoutmask":inoutmask}

def interpolate1DDensity(inx, inmask, influx, invar, outx):
    assert inx.shape==inmask.shape==influx.shape==invar.shape
    innx,  = inx.shape

    inx = inx.astype(np.float64)
    influx = influx.astype(np.float64)
    invar  = invar.astype(np.float64)
    inmask = (inmask>0).astype(np.int64)
    outx = outx.astype(np.float64)
    outnx, = outx.shape

    inshape = inx.shape
    inoutmask = np.zeros(inshape, dtype=np.float64)
    inoutmask_c = np.ctypeslib.as_ctypes(inoutmask)
    inLx = np.zeros(inshape, dtype=np.float64) * np.nan
    inRx = np.zeros(inshape, dtype=np.float64) * np.nan
    inLx_c = np.ctypeslib.as_ctypes(inLx)
    inRx_c = np.ctypeslib.as_ctypes(inRx)

    inx_c = np.ctypeslib.as_ctypes(inx)

    influx_c = np.ctypeslib.as_ctypes(influx)
    invar_c = np.ctypeslib.as_ctypes(invar)
    inmask_c = np.ctypeslib.as_ctypes(inmask)

    outshape = outx.shape
    outLx = np.zeros(outshape, dtype=np.float64) * np.nan
    outRx = np.zeros(outshape, dtype=np.float64) * np.nan
    outLx_c = np.ctypeslib.as_ctypes(outLx)
    outRx_c = np.ctypeslib.as_ctypes(outRx)

    outx_c = np.ctypeslib.as_ctypes(outx)
    outflux = np.zeros(outshape, dtype=np.float64)
    outvar  = np.zeros(outshape, dtype=np.float64)
    outmask = np.zeros(outshape, dtype=np.float64)
    outflux_c = np.ctypeslib.as_ctypes(outflux)
    outvar_c  = np.ctypeslib.as_ctypes(outvar)
    outmask_c  = np.ctypeslib.as_ctypes(outmask)
    outx_c  = np.ctypeslib.as_ctypes(outx)
    success = interpolate1DDensity_c(
                    inshape[0],
                    inx_c,
                    influx_c, invar_c,
                    inmask_c, inoutmask_c,
                    inLx_c,  inRx_c,
                    outshape[0],
                    outx_c,
                    outflux_c, outvar_c,
                    outLx_c,  outRx_c,
                    outmask_c)
    if success<0:
        raise RuntimeError("have error in interpolate2D")
    return {"outflux":outflux,
            "outvar":outvar,
            "outmask":outmask,
            "inoutmask":inoutmask}

#endregion: py func

tol = 10
PRINT_TEST_DEBUG = 1
if __name__ == '__main__':
    class MainTest(unittest.TestCase):
        def setUp(self):
            self.starttime = None
            self.midtime = None
            self.endtime = None
            self.stdout = StringIO()
            sys.stdout = self.stdout
        def tearDown(self):
            sys.stdout = sys.__stdout__
            if PRINT_TEST_DEBUG:
                print(self.stdout.getvalue())
            if self.starttime is not None:
                print("{}s {}s".format(self.midtime-self.starttime, self.endtime-self.midtime))
        def getPolygon(self, polygon):
            number = polygon.contents.num
            if not number:
                return None, None, None
            thisV = polygon.contents.P
            x = []
            y = []
            flag = []
            for i in range(number):
                x.append(thisV.contents.x)
                y.append(thisV.contents.y)
                flag.append(thisV.contents.flag)
                thisV = thisV.contents.next
            x.append(x[0])
            y.append(y[0])
            flag.append(flag[0])
            return x,y,flag
        #@unittest.skipIf(not PRINT_TEST_DEBUG, 'not print test debug')
        @unittest.skip('debug by eye')
        def test_0_Polygon(self):
            polygon = newPolygon()
            showPolygon(polygon)
            addVertex(polygon, 0,0,0,0)
            addVertex(polygon, 1,1,1,0)
            addVertex(polygon, 2,2,2,0)
            addVertex(polygon, 3,3,3,0)
            addVertex(polygon, 4,4,4,0)
            showPolygon(polygon)
            addVertex(polygon, 5,5,5,4)
            howPolygon(polygon)
            addVertex(polygon, 6,6,6,6)
            showPolygon(polygon)
            addVertex(polygon, 6,6,6,8)
            delVertex(polygon, 0)
            showPolygon(polygon)
            delVertex(polygon, 5)
            showPolygon(polygon)
            delVertex(polygon, 4)
            showPolygon(polygon)
            delVertex(polygon, 3)
            showPolygon(polygon)
            delVertex(polygon, 1)
            showPolygon(polygon)
            delVertex(polygon, 0)
            showPolygon(polygon)
            delVertex(polygon, 0)
            showPolygon(polygon)
            delVertex(polygon, 0)
        @unittest.skip('skip')
        def test_0_searchsorted(self):
            N = 9999
            for i in range(N):
                xx = np.arange(100)
                deltax = np.random.uniform(-0.2,0.2,100)
                xx = xx + deltax
                value = np.random.uniform(10,90)
                result = searchsorted(xx, value)
                result_np = xx.searchsorted(value)
                self.assertEqual(result, result_np)
        @unittest.skip('debug by eye')
        def test_1_intersec(self):
            plot = 1
            def test4points(ax,ay,bx,by,cx,cy,dx,dy):
                va = newVertex(ax,ay,0)
                vb = newVertex(bx,by,0)
                vc = newVertex(cx,cy,0)
                vd = newVertex(dx,dy,0)
                vv = intersect(va,vb,vc,vd)
                vx = vv.contents.x
                vy = vv.contents.y
                if plot:
                    fig = plt.figure()
                    axis = fig.add_subplot(111)
                    t = np.linspace(-10,10,9999)
                    l1x = ax + (bx-ax)*t
                    l1y = ay + (by-ay)*t
                    l2x = cx + (dx-cx)*t
                    l2y = cy + (dy-cy)*t
                    axis.scatter(ax, ay, s=10, color='green')
                    axis.scatter(bx, by, s=10, color='green')
                    axis.scatter(cx, cy, s=10, color='red')
                    axis.scatter(dx, dy, s=10, color='red')
                    axis.scatter(vx, vy, s=15, color='black')
                    axis.plot(l1x, l1y, color="green")
                    axis.plot(l2x, l2y, color='red')
                    xcen = (ax+bx+cx+dx+vx)/5
                    ycen = (ay+by+cy+dy+vy)/5
                    xsize = np.max([ax,bx,cx,dx,vx])-np.min([ax,bx,cx,dx,vx])
                    ysize = np.max([ay,by,cy,dy,vy])-np.min([ay,by,cy,dy,vy])
                    xlim = [xcen-xsize, xcen+xsize]
                    ylim = [ycen-ysize, ycen+ysize]
                    axis.set_xlim(xlim)
                    axis.set_ylim(ylim)
                    plt.show(0)
            def randomTest():
                min = -10
                max =  10
                test4points(
                    np.random.uniform(min,max),
                    np.random.uniform(min,max),
                    np.random.uniform(min,max),
                    np.random.uniform(min,max),
                    np.random.uniform(min,max),
                    np.random.uniform(min,max),
                    np.random.uniform(min,max),
                    np.random.uniform(min,max),
                )
            test4points(0,0,1,1,1,0,0,1)
            for i in range(20):
                randomTest()
        @unittest.skip('skip')
        def test_2_vectors(self):
            def equal3(A,V):
                if np.abs(A[0]-V.contents.x)>tol*eps:
                    return False
                if np.abs(A[1]-V.contents.y)>tol*eps:
                    return False
                if np.abs(A[2]-V.contents.z)>tol*eps:
                    return False
                return True
            def equal2(A,V):
                if np.abs(A[0]-V.contents.x)>tol*eps:
                    return False
                if np.abs(A[1]-V.contents.y)>tol*eps:
                    return False
                return True
            for i in range(100):
                ax = np.random.randn()
                ay = np.random.randn()
                az = np.random.randn()
                bx = np.random.randn()
                by = np.random.randn()
                bz = np.random.randn()
                Aa = np.array([ax, ay, az])
                Ab = np.array([bx, by, bz])
                Va = newVector3(ax,ay,az)
                Vb = newVector3(bx,by,bz)
                self.assertTrue(equal3(-Aa, neg3(Va)))
                self.assertTrue(equal3(-Ab, neg3(Vb)))
                self.assertTrue(np.abs(np.sqrt((Aa*Aa).sum())-norm3(Va))<tol*eps )
                self.assertTrue(np.abs(np.sqrt((Ab*Ab).sum())-norm3(Vb))<tol*eps )
                self.assertTrue(equal3(Aa+Ab, plus3(Va,Vb)))
                self.assertTrue(equal3(Aa-Ab, minus3(Va,Vb)))
                self.assertTrue(np.abs((Aa*Ab).sum()-dot3(Va, Vb))<tol*eps)
                self.assertTrue(equal3(np.cross(Aa, Ab), cross3(Va, Vb)))
            for i in range(100):
                ax = np.random.randn()
                ay = np.random.randn()
                bx = np.random.randn()
                by = np.random.randn()
                Aa = np.array([ax, ay])
                Ab = np.array([bx, by])
                Va = newVector2(ax,ay)
                Vb = newVector2(bx,by)
                self.assertTrue(equal2(-Aa, neg2(Va)))
                self.assertTrue(equal2(-Ab, neg2(Vb)))
                self.assertTrue(np.abs(np.sqrt((Aa*Aa).sum())-norm2(Va))<tol*eps )
                self.assertTrue(np.abs(np.sqrt((Ab*Ab).sum())-norm2(Vb))<tol*eps )
                self.assertTrue(equal2(Aa+Ab, plus2(Va,Vb)))
                self.assertTrue(equal2(Aa-Ab, minus2(Va,Vb)))
                self.assertTrue(np.abs((Aa*Ab).sum()-dot2(Va, Vb))<tol*eps)
        @unittest.skip('skip')
        def test_3_triangle_triangle_intersec(self):
            debug_dir = "test_triangle_triangle_intersec"
            os.makedirs(debug_dir, exist_ok=True)
            for i in range(1000):
                print("process ", i)
                min = -10
                max =  10
                fig = plt.figure()
                ax = fig.add_subplot(111)
                inPara = [np.random.uniform(-10,10) for i in range(12)]
                inParaA = np.array(inPara)
                Axs = inParaA[[0,2,4]]
                Ays = inParaA[[1,3,5]]
                Bxs = inParaA[[6,8,10]]
                Bys = inParaA[[7,9,11]]
                Axsp = inParaA[[0,2,4,0]]
                Aysp = inParaA[[1,3,5,1]]
                Bxsp = inParaA[[6,8,10,6]]
                Bysp = inParaA[[7,9,11,7]]
                ax.scatter(Axs, Ays, s=10, color='green')
                ax.scatter(Bxs, Bys, s=10, color='red')
                ax.plot(Axsp, Aysp, color="green")
                ax.plot(Bxsp, Bysp, color="red")
                ax.text(Axs[0], Ays[0], "A")# ({:7.3f}, {:7.3f})".format(Axs[0], Ays[0]))
                ax.text(Axs[1], Ays[1], "B")# ({:7.3f}, {:7.3f})".format(Axs[1], Ays[1]))
                ax.text(Axs[2], Ays[2], "C")# ({:7.3f}, {:7.3f})".format(Axs[2], Ays[2]))
                ax.text(Bxs[0], Bys[0], "D")# ({:7.3f}, {:7.3f})".format(Bxs[0], Bys[0]))
                ax.text(Bxs[1], Bys[1], "E")# ({:7.3f}, {:7.3f})".format(Bxs[1], Bys[1]))
                ax.text(Bxs[2], Bys[2], "F")# ({:7.3f}, {:7.3f})".format(Bxs[2], Bys[2]))
                poly = TriangleTrianglePolygon(*inPara)
                Cxs, Cys, flags = self.getPolygon(poly)
                if Cxs is not None:
                    ps = PatchCollection([patches.Polygon(np.array(list(zip(Cxs, Cys))), True)], alpha=0.4)
                    ax.add_collection(ps)
                fig.savefig(os.path.join(debug_dir, "{}.png".format(i)))
                plt.close(fig)
        @unittest.skip('skip')
        def test_4_triangle_square(self):
            doPlot = 0
            N=999999
            if doPlot:
                debug_dir = "test_triangle_square_intersec"
                os.makedirs(debug_dir, exist_ok=True)
            count = 0
            while count<N:
                if doPlot:
                    print("process ", count)
                min = -10
                max =  10
                Axcen, Aycen = np.random.uniform(min*0.8, max*0.8),np.random.uniform(min*0.8, max*0.8)
                AHalfSize = np.random.uniform(0.1, 8)
                Axs  = [Axcen-AHalfSize, Axcen-AHalfSize, Axcen+AHalfSize, Axcen+AHalfSize]
                Axsp = [Axcen-AHalfSize, Axcen-AHalfSize, Axcen+AHalfSize, Axcen+AHalfSize, Axcen-AHalfSize]
                Ays  = [Aycen+AHalfSize, Aycen-AHalfSize, Aycen-AHalfSize, Aycen+AHalfSize]
                Aysp = [Aycen+AHalfSize, Aycen-AHalfSize, Aycen-AHalfSize, Aycen+AHalfSize, Aycen+AHalfSize]
                Bxs = [np.random.uniform(min,max) for i in range(3)]
                Bys = [np.random.uniform(min,max) for i in range(3)]
                Bxsp = copy.copy(Bxs)
                Bxsp.append(Bxs[0])
                Bysp = copy.copy(Bys)
                Bysp.append(Bys[0])
                inParas = [each for eachtuple in itertools.chain(zip(Axs, Ays), zip(Bxs, Bys))
                                    for each in eachtuple]
                poly = TriangleSquarePolygon(*inParas)
                Cxs, Cys, flags = self.getPolygon(poly)
                fullArea, overlapArea = areaTriangleSquare(*inParas)
                if overlapArea==0.:
                    if np.random.uniform(0,1)<0.9:# only have 1/10 of all the empty area
                        continue
                count += 1
                if doPlot:
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    ax.scatter(Axs, Ays, s=10, color='green')
                    ax.scatter(Bxs, Bys, s=10, color='red')
                    ax.plot(Axsp, Aysp, color="green")
                    ax.plot(Bxsp, Bysp, color="red")
                    ax.text(Axs[0], Ays[0], "1")# ({:7.3f}, {:7.3f})".format(Axs[0], Ays[0]))
                    ax.text(Axs[1], Ays[1], "2")# ({:7.3f}, {:7.3f})".format(Axs[1], Ays[1]))
                    ax.text(Axs[2], Ays[2], "3")# ({:7.3f}, {:7.3f})".format(Axs[2], Ays[2]))
                    ax.text(Axs[3], Ays[3], "4")# ({:7.3f}, {:7.3f})".format(Axs[3], Ays[3]))

                    ax.text(Bxs[0], Bys[0], "D")# ({:7.3f}, {:7.3f})".format(Bxs[0], Bys[0]))
                    ax.text(Bxs[1], Bys[1], "E")# ({:7.3f}, {:7.3f})".format(Bxs[1], Bys[1]))
                    ax.text(Bxs[2], Bys[2], "F")# ({:7.3f}, {:7.3f})".format(Bxs[2], Bys[2]))
                    if Cxs is not None:
                        ps = PatchCollection([patches.Polygon(np.array(list(zip(Cxs, Cys))), True)], alpha=0.4)
                        ax.add_collection(ps)
                    ax.set_title("fullArea: {:10.5f}, overlapArea: {:10.5f}, percent: {:7.3f}".
                            format(fullArea, overlapArea, overlapArea/fullArea*100.0))
                    fig.savefig(os.path.join(debug_dir, "{}.png".format(count)))
                    plt.close(fig)
        @unittest.skip('skip')
        def test_5_get4map_simple(self):
            m = 5
            n = 6
            x = np.arange(3,3+n)
            y = np.arange(6,6+m)
            deltad = 0.2
            deltax = np.linspace(0,0.8,n)
            deltadeltax = np.random.uniform(-deltad,deltad, (m,n))
            deltay = np.linspace(0,0.7,m).reshape(m,1)
            deltadeltay = np.random.uniform(-deltad,deltad, (m,n))
            xx,yy = np.meshgrid(x,y)
            #print(xx.shape, yy.shape, deltax.shape,
            #      deltay.shape, deltadeltax.shape, deltadeltay.shape)

            xx = xx+deltay
            yy = yy+deltax
            xx = xx+deltadeltax
            yy = yy+deltadeltay
            mask = np.ones(xx.shape, dtype=np.bool)
            get4map_plot(xx, yy, mask)
        @unittest.skip('skip')
        def test_6_get4map_test1(self):
            obj  = fits.open('./test/map_temp.fits')

            data = obj[1].data
            var  = obj[2].data
            mask = obj[3].data.astype(np.bool)
            pinholemap = obj[4].data
            pixelmap = obj[5].data
            wavemap  = obj[6].data
            sky = obj[7].data

            data[~mask] = np.nan
            var[~mask] = np.nan
            pinholemap[~mask] = np.nan
            pixelmap[~mask] = np.nan
            wavemap[~mask] = np.nan
            sky[~mask] = np.nan

            ycut = slice(260, 276)
            xcut = slice(25, 86)
            ny, nx = data[ycut, xcut].shape
            wrapper = 5
            ny += 2*wrapper
            nx += 2*wrapper
            xcutcut = slice(wrapper, -wrapper)
            ycutcut = slice(wrapper, -wrapper)

            data_       = data[ycut, xcut]
            var_        = var[ycut, xcut]
            mask_       = mask[ycut, xcut]
            pinholemap_ = pinholemap[ycut, xcut]
            pixelmap_   = pixelmap[ycut, xcut]
            wavemap_    = wavemap[ycut, xcut]
            sky_        = sky[ycut, xcut]

            data       = np.zeros((ny, nx))*np.nan
            var        = np.zeros((ny, nx))*np.nan
            mask       = np.zeros((ny, nx), dtype=np.bool)
            pinholemap = np.zeros((ny, nx))*np.nan
            pixelmap   = np.zeros((ny, nx))*np.nan
            wavemap    = np.zeros((ny, nx))*np.nan
            sky        = np.zeros((ny, nx))*np.nan

            data[ycutcut, xcutcut]       = data_[:,:]
            var[ycutcut, xcutcut]        = var_[:,:]
            mask[ycutcut, xcutcut]       = mask_[:,:]
            pinholemap[ycutcut, xcutcut] = pinholemap_[:,:]
            pixelmap[ycutcut, xcutcut]   = pixelmap_[:,:]
            wavemap[ycutcut, xcutcut]    = wavemap_[:,:]
            sky[ycutcut, xcutcut]        = sky_[:,:]

            #from fmajorAstroUtils.DS10 import ds9
            #d = ds9.ds10('test')
            #d.newFrame
            #d.set_np2arr(data)
            #d.newFrame
            #d.set_np2arr(var)
            #d.newFrame
            #d.set_np2arr(pinholemap)
            #d.newFrame
            #d.set_np2arr(pixelmap)
            #d.newFrame
            #d.set_np2arr(wavemap)
            #d.newFrame
            #d.set_np2arr(sky)

            ny, nx = data.shape
            y = np.arange(ny)
            x = np.arange(nx)
            xx,yy = np.meshgrid(x,y)
            get4map_plot(pinholemap, pixelmap, mask)
        @unittest.skip('skip')
        def test_6_get4map_test2(self):
            obj  = fits.open('./test/map_temp.fits')

            data = obj[1].data
            var  = obj[2].data
            mask = obj[3].data.astype(np.bool)
            pinholemap = obj[4].data
            pixelmap = obj[5].data
            wavemap  = obj[6].data
            sky = obj[7].data

            data[~mask] = np.nan
            var[~mask] = np.nan
            pinholemap[~mask] = np.nan
            pixelmap[~mask] = np.nan
            wavemap[~mask] = np.nan
            sky[~mask] = np.nan

            ycut = slice(260, 276)
            xcut = slice(25, 86)
            ny, nx = data[ycut, xcut].shape
            wrapper = 5
            ny += 2*wrapper
            nx += 2*wrapper
            xcutcut = slice(wrapper, -wrapper)
            ycutcut = slice(wrapper, -wrapper)

            data       = data[ycut, xcut]
            var        = var[ycut, xcut]
            mask       = mask[ycut, xcut]
            pinholemap = pinholemap[ycut, xcut]
            pixelmap   = pixelmap[ycut, xcut]
            wavemap    = wavemap[ycut, xcut]
            sky        = sky[ycut, xcut]

            ny, nx = data.shape
            y = np.arange(ny)
            x = np.arange(nx)
            xx,yy = np.meshgrid(x,y)
            get4map_plot(pinholemap, pixelmap, mask)
        @unittest.skip('skip')
        def test_6_get4map_test3(self):
            '''2d diff test'''
            obj  = fits.open('./test/map_temp.fits')

            data = obj[1].data
            var  = obj[2].data
            mask = obj[3].data.astype(np.bool)
            pinholemap = obj[4].data
            pixelmap = obj[5].data
            wavemap  = obj[6].data
            sky = obj[7].data

            data[~mask] = np.nan
            var[~mask] = np.nan
            pinholemap[~mask] = np.nan
            pixelmap[~mask] = np.nan
            wavemap[~mask] = np.nan
            sky[~mask] = np.nan

            #from fmajorAstroUtils.DS10 import ds9
            #d = ds9.ds10('test')
            #d.newFrame
            #d.set_np2arr(data)
            #d.newFrame
            #d.set_np2arr(var)
            #d.newFrame
            #d.set_np2arr(pinholemap)
            #d.newFrame
            #d.set_np2arr(pixelmap)
            #d.newFrame
            #d.set_np2arr(wavemap)
            #d.newFrame
            #d.set_np2arr(sky)

            ny, nx = data.shape
            y = np.arange(ny)
            x = np.arange(nx)
            xx,yy = np.meshgrid(x,y)
            get4map_plot_large(pinholemap, pixelmap, mask)

        #@unittest.skip('skip new method')
        def test_6_get4map_test4_new_method(self):
            obj = fits.open("./test/map_temp2.fits")

            pinholemap = obj['pinholemap'].data
            pixelmap = obj['pixelmap'].data

            xcut = slice(None, 80)
            ycut = slice(None, 80)

            xcut = slice(43,  107)
            ycut = slice(279, 317)

            xx = pinholemap[ycut, xcut]
            yy =   pixelmap[ycut, xcut]
            mask = np.isfinite(xx)
            get4map_plot(xx, yy, mask)

        def get4map_range(self, xcut, ycut):
            obj  = fits.open('./test/map_temp.fits')

            data = obj[1].data
            var  = obj[2].data
            mask = obj[3].data.astype(np.bool)
            pinholemap = obj[4].data
            pixelmap = obj[5].data
            wavemap  = obj[6].data
            sky = obj[7].data

            data[~mask] = np.nan
            var[~mask] = np.nan
            pinholemap[~mask] = np.nan
            pixelmap[~mask] = np.nan
            wavemap[~mask] = np.nan
            sky[~mask] = np.nan

            #ycut = slice(896)
            #xcut = slice(131)
            ny, nx = data[ycut, xcut].shape
            wrapper = 5
            ny += 2*wrapper
            nx += 2*wrapper
            xcutcut = slice(wrapper, -wrapper)
            ycutcut = slice(wrapper, -wrapper)

            data       = data[ycut, xcut]
            var        = var[ycut, xcut]
            mask       = mask[ycut, xcut]
            pinholemap = pinholemap[ycut, xcut]
            pixelmap   = pixelmap[ycut, xcut]
            wavemap    = wavemap[ycut, xcut]
            sky        = sky[ycut, xcut]

            ny, nx = data.shape
            y = np.arange(ny)
            x = np.arange(nx)
            xx,yy = np.meshgrid(x,y)
            get4map_plot(pinholemap, pixelmap, mask)
        @unittest.skip('skip')
        def test_7_get4map_range(self):
            ycut = slice(None, 25)
            xcut = slice(None, 61)
            self.get4map_range(xcut, ycut)
            ycut = slice(900,950)
            xcut = slice(163)
            self.get4map_range(xcut, ycut)

        def get4map_flux(self, xcut=None, ycut=None, wrapper=5, full=False, name=""):
            obj  = fits.open('./test/map_temp.fits')

            data = obj[1].data
            var  = obj[2].data
            mask = obj[3].data.astype(np.bool)
            pinholemap = obj[4].data
            pixelmap = obj[5].data
            wavemap  = obj[6].data
            sky = obj[7].data

            outx = np.arange(-5,data.shape[1]+5+1)
            outy = np.arange(-5,data.shape[0]+5+1)

            data[~mask] = np.nan
            var[~mask] = np.nan
            pinholemap[~mask] = np.nan
            pixelmap[~mask] = np.nan
            wavemap[~mask] = np.nan
            sky[~mask] = np.nan

            if not full:
                ny, nx = data[ycut, xcut].shape
                ny += 2*wrapper
                nx += 2*wrapper
                xcutcut = slice(wrapper, -wrapper)
                ycutcut = slice(wrapper, -wrapper)

                data_       = data[ycut, xcut]
                var_        = var[ycut, xcut]
                mask_       = mask[ycut, xcut]
                pinholemap_ = pinholemap[ycut, xcut]
                pixelmap_   = pixelmap[ycut, xcut]
                wavemap_    = wavemap[ycut, xcut]
                sky_        = sky[ycut, xcut]

                data       = np.zeros((ny, nx))*np.nan
                var        = np.zeros((ny, nx))*np.nan
                mask       = np.zeros((ny, nx), dtype=np.bool)
                pinholemap = np.zeros((ny, nx))*np.nan
                pixelmap   = np.zeros((ny, nx))*np.nan
                wavemap    = np.zeros((ny, nx))*np.nan
                sky        = np.zeros((ny, nx))*np.nan

                data[ycutcut, xcutcut]       = data_[:,:]
                var[ycutcut, xcutcut]        = var_[:,:]
                mask[ycutcut, xcutcut]       = mask_[:,:]
                pinholemap[ycutcut, xcutcut] = pinholemap_[:,:]
                pixelmap[ycutcut, xcutcut]   = pixelmap_[:,:]
                wavemap[ycutcut, xcutcut]    = wavemap_[:,:]
                sky[ycutcut, xcutcut]        = sky_[:,:]

            result = interpolate2D(pinholemap, pixelmap, mask, data, var, outx, outy);
            outflux = result['outflux']
            outvar = result['outvar']
            outmask = result['outmask']
            inoutmask = result['inoutmask']
            result = interpolate2D(pinholemap, pixelmap, mask, sky, var, outx, outy);
            outsky = result['outflux']
            from fmajorAstroUtils.DS10 import ds9
            d = ds9.ds10(name)
            d.newFrame
            d.set_np2arr(data)
            d.newFrame
            d.set_np2arr(var)
            d.newFrame
            d.set_np2arr(pinholemap)
            d.newFrame
            d.set_np2arr(pixelmap)
            d.newFrame
            d.set_np2arr(wavemap)
            d.newFrame
            d.set_np2arr(sky)
            d.newFrame
            d.set_np2arr(outflux)
            d.newFrame
            d.set_np2arr(outsky)
            d.newFrame
            d.set_np2arr(outvar)
            d.newFrame
            d.set_np2arr(outmask)
            d.newFrame
            d.set_np2arr(inoutmask)
            initFlux = np.nansum(data)
            newFlux = np.nansum(outflux)
            initvar = np.nansum(var)
            newvar = np.nansum(outvar)
            initsky = np.nansum(sky)
            newsky = np.nansum(outsky)
            print("{:15.7f} ?= {:15.7f}, resdual: {}".
                    format(initFlux, newFlux, (newFlux-initFlux)/initFlux))
            print("{:15.7f} ?= {:15.7f}, resdual: {}".
                    format(initvar, newvar, (newvar-initvar)/initvar))
            print("{:15.7f} ?= {:15.7f}, resdual: {}".
                    format(initsky, newsky, (newsky-initsky)/initsky))
        #@unittest.skip('skip')
        def test_8_get4map_flux(self):
            ycut = slice(None, 500)
            xcut = slice(None, 100)
            self.get4map_flux(xcut, ycut, name="testpart")
            self.get4map_flux(full=True, name="testfull")

        def gaussion(self, x, A, mu, sigma):
            return A*np.exp(-0.5*(x-mu)**2/(2*sigma**2))/(np.sqrt(2*np.pi*sigma))
        #@unittest.skip('skip')
        def test_9_interpolate1D(self):
            N=1
            doPlot = 1
            for i in range(N):
                xstart = 0
                xend = 100
                inN = 101
                step = (xend+1-xstart)/inN
                wrapper = 5
                offset = 0.5
                inx = np.linspace(xstart, xend, inN)
                outx = np.arange(xstart-wrapper, xend+wrapper, step)+offset
                dense = np.linspace(xstart-wrapper, xend+wrapper, 9999)
                influx = np.random.randn(inN)
                Ng = np.random.randint(1,4)
                SNR = 3
                print("insize:{} outsize:{}".format(inN, outx.size))
                if doPlot:
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                for j in range(Ng):
                    sigma = np.random.uniform(0.4,2)
                    mu = np.random.uniform(xstart+10*sigma, xend-10*sigma)
                    A = np.random.uniform(10*SNR,100*SNR)
                    yy = self.gaussion(inx, A, mu, sigma)
                    influx += yy
                    ax.plot([mu, mu], [0,A], '--', lw=1)
                if doPlot:
                    ax.scatter(inx, influx, edgecolor='None', s=50, marker='x')
                    #ax.plot(inx, influx, color='black')
                invar = np.abs(influx)
                from scipy.interpolate import interp1d
                interpLinear = interp1d(inx, influx, kind='linear',
                                        bounds_error=False, fill_value=0)
                outLinear = interpLinear(outx)
                interpSpline = interp1d(inx, influx, kind='cubic',
                                        bounds_error=False, fill_value=0)
                outSpline = interpSpline(outx)
                mask = np.ones(inN, dtype=np.bool)
                outFlux = interpolate1D(inx, mask, influx, invar, outx)['outflux']
                if doPlot:
                    ax.scatter(outx, outLinear, color='b', edgecolor="None", s=15)
                    ax.plot(dense, interpLinear(dense), 'b--', lw=0.5)
                    ax.scatter(outx, outSpline, color='g', edgecolor="None", s=15)
                    ax.plot(dense, interpSpline(dense), 'g--', lw=0.5)
                    ax.scatter(outx, outFlux, edgecolor="None", marker="+", s=60)
            if doPlot:
                plt.show(0)

        @unittest.skip('skip')
        def test_testdata(selfname=""):
            obj  = fits.open('./test/map_temp.fits')

            data = obj[1].data
            var  = obj[2].data
            mask = obj[3].data.astype(np.bool)
            pinholemap = obj[4].data
            pixelmap = obj[5].data
            wavemap  = obj[6].data
            sky = obj[7].data

            outx = np.arange(-5,data.shape[1]+5+1)
            outy = np.arange(-5,data.shape[0]+5+1)

            data[~mask] = np.nan
            var[~mask] = np.nan
            pinholemap[~mask] = np.nan
            pixelmap[~mask] = np.nan
            wavemap[~mask] = np.nan
            sky[~mask] = np.nan

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.scatter(pixelmap.flatten(), wavemap.flatten(),
                       s=1, edgecolor="None")
            plt.show()


    suite = unittest.TestLoader().loadTestsFromTestCase(MainTest)
    unittest.TextTestRunner(verbosity=2).run(suite)

