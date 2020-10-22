from __future__ import (print_function, absolute_import, division)
try:
    input = raw_input
except NameError:
    pass
from   pyds9 import DS9
import astropy.io.fits as pyfits
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
import re
import copy
import sys
import os
import getopt
import time
from .goodFigure import GoodFigure
import pdb
import pickle
from scipy.optimize import curve_fit
from astropy import wcs as WCS
import warnings
warnings.simplefilter('ignore', WCS.FITSFixedWarning)

def gaussian(x, A, mu, sigma, C):
    return A * np.exp(- (x-mu)**2/(2 * sigma**2)) + C

#!! todo 
'''
    change region size
'''
if sys.version_info.major==3:
    unicode = str


#<== mode setter
modeList = ["none", "region", "crosshair",
            "colorbar", "pan", "zoom", "rotate"]
class ModeSetter(object):
    def __init__(self, father):
        self._father = father
    @property
    def pan(self):
        self._father.set("mode pan")
    @property
    def rotate(self):
        self._father.set("mode rotate")
    @property
    def none(self):
        self._father.set("mode none")
    @property
    def region(self):
        self._father.set("mode region")
    @property
    def zoom(self):
        self._father.set("mode zoom")
    @property
    def crosshair(self):
        self._father.set("mode crosshair")
    @property
    def crop(self):
        self._father.set("mode crop")

#==> mode setter

#<== match setter
matchLevel2 = ["wcs", "image", "physical", "amplifier", "detector"]
submatchStr="""
@property
def {match}(self):
    if self._name=="slice" and '{match}' in ["physical", "amplifier", "detector"]:
        print("unchanged")
    else:
        self._father.set("match {{}} {match}".format(self._name))
"""
class SubMatchSetter(object):
    def __init__(self, father, name):
        self._father = father
        self._name = name

    for eachSubMatch in matchLevel2:
        exec(submatchStr.format(match=eachSubMatch), globals(), locals())

matchSetterStr="""
@property
def {matchname}(self):
    self._father.set("match {match}")
"""
matchLevel1 = ["bin", "scale", "colorbar", "smooth", "scalelimits"]
matchLevel1Name = ["bin", "scale", "colorbar", "smooth", "limits"]
matchLevel12 = ["frame", "crosshair", "crop", "slice"]

class MatchSetter(object):
    def __init__(self, father):
        self._father = father
        for eachSubMatch in matchLevel12:
            exec("self.{f} = SubMatchSetter(self._father,'{f}')".format(f=eachSubMatch), globals(), locals())

    for eachMatch,eachMatchName in zip(matchLevel1, matchLevel1Name):
        todo = matchSetterStr.format(match = eachMatch, matchname=eachMatchName)
        exec(todo, globals(), locals())

    @property
    def limits(self):
        thisLimit = self._father.get("scale limists")
        self._father.set("match scale limits {}".format(thisLimit))
    @limits.setter
    def limits(self, value):
        if len(value)==2:
            try:
                float(value[0])
                float(value[1])
                self.set("match scale limits {} {}".\
                        format(value[0], value[1]))
                return
            except:
                raise Exception("invaild value or data type for {}: {}, type {}".\
                    format(sys._getframe().f_code.co_name, value, type(value)))

#==> match setter

#<== lock setter
lockLevel2 = ["wcs", "none", "image", "physical", "amplifier", "detector"]
sublockStr="""
@property
def {lock}(self):
    if self._name=="slice" and '{lock}' in ["physical", "amplifier", "detector"]:
        print("unchanged")
    else:
        self._father.set("lock {{}} {lock}".format(self._name))
"""
class SubLockSetter(object):
    def __init__(self, father, name):
        self._father = father
        self._name = name

    for eachSubLock in lockLevel2:
        exec(sublockStr.format(lock=eachSubLock), globals(), locals())
    @property
    def get(self):
        return self._father.get("lock {}".format(self._name))

lockSetterStr="""
@property
def {lockname}(self):
    result = self._father.get("lock {lock}")
    if result[0]=="y":
        self._father.set("lock {lock} off")
        return 0
    elif result[0]=="n":
        self._father.set("lock {lock} on")
        return 1
    else:
        print(result)
        raise Exception("unknown result from {lock}")

"""

lockSetterStrOn="""
@property
def {lockname}_on(self):
    result = self._father.set("lock {lock} on")
"""
lockSetterStrOff="""
@property
def {lockname}_off(self):
    result = self._father.set("lock {lock} off")
"""

lockLevel1 = ["bin", "scale", "colorbar", "smooth","scalelimits"]
lockLevel1Name = ["bin", "scale", "colorbar", "smooth", "limits"]
lockLevel12 = ["frame", "crosshair", "crop", "slice"]

class LockSetter(object):
    def __init__(self, father):
        self._father = father
        for eachSubLock in lockLevel12:
            exec("self.{f} = SubLockSetter(self._father,'{f}')".format(f=eachSubLock), globals(), locals())


    for eachLock, eachLockName in zip(lockLevel1, lockLevel1Name):
        todo = lockSetterStr.format(lock = eachLock, lockname=eachLockName)
        exec(todo, globals(), locals())
    for eachLock, eachLockName in zip(lockLevel1, lockLevel1Name):
        todo = lockSetterStrOn.format(lock = eachLock, lockname=eachLockName)
        exec(todo, globals(), locals())
    for eachLock, eachLockName in zip(lockLevel1, lockLevel1Name):
        todo = lockSetterStrOff.format(lock = eachLock, lockname=eachLockName)
        exec(todo, globals(), locals())



    @property
    def get(self):
        allResult = []
        for eachLock in lockLevel1:
            result = self._father.get("lock {}".format(eachLock))
            if result[0]=="y":
                allResult.append(eachLock)
            elif result[0]!="n":
                print(result)
                raise Exception("unknown result from {}".format(eachLock))
        return allResult

#==> lock setter

#<== regions
regionTypes = ["circle", "ellipse", "box", "polygon", "point",
               "line", "vector", "text", "ruler", "compass",
               "projection", "annulus", "ellipse", "box", "panda",
               "epanda", "bpanda", "composite"]

str_in_parentheses = r"\((.*?)\)"
re_in_parentheses = re.compile(str_in_parentheses)
str_before_parentheses = r"(\w*)\("
re_before_parentheses = re.compile(str_before_parentheses)

str_tag = r"tag=\{(.*?)\}"
re_tag = re.compile(str_tag)

str_inside = r"\{(.*?)\}"
re_inside = re.compile(str_inside)

class Region(object):
    def __init__(self):
        pass
        #self.defaultFont = {"fontname":"times", "size":"12", "bold":"normal", "italic":"roman"}
        #self.defaultFontKeys = ["fontname", "size", "bold", "italic"]
    def dict2eq(self, d):# in d, text and tag and fond already inside ''
        tagList = d.get("tag",[])
        result = ",".join(["{}={}".format(eachKey, d[eachKey]) for eachKey in d if eachKey != "tag"])
        resultTag = ",".join(["tag={}".format(eachTag) for eachTag in tagList])
        if resultTag:
            if result:
                return ",".join([result, resultTag])
            else:
                return resultTag
        else:
            return result
    def processRegionLine(self,s):
        if "file format" in s:
            return s
        if s in ["image", "physical"]:
            return s
        if s == "wcs;":
            return "wcs"
        if "global" in s:
            command="global"
            params=""
            afterCommand=s
        else:
            aux = s.split("#")
            if "(" in aux[0] and ")" in aux[0]: # command in first part
                command = re_before_parentheses.findall(aux[0])[0]
                params = re_in_parentheses.findall(aux[0])[0].split(",")
                if len(aux)==2:
                    afterCommand = aux[1].strip()
                else:
                    afterCommand = ""
            elif "(" in aux[1] and ")" in aux[1]:
                command = re_before_parentheses.findall(aux[1])[0]
                params = re_in_parentheses.findall(aux[1])[0].split(",")
                afterCommand = aux[1].split(")")[-1].strip()
            else:
                raise Exception("string {} not understood".format(s))
        configAux = afterCommand.split("=")
        totalN = len(configAux)-1
        if totalN == 0:
            return command, params, {}
        configPairs = [""] * totalN
        haveTag = False
        for i,each in enumerate(configAux):
            thisPair = each.split(" ")
            #print(i, thisPair)
            N = len(thisPair)
            if i>0:
                if i<totalN:
                    thisKey = thisPair[-1].strip()
                    if thisKey=="tag":
                        haveTag=True
                    configPairs[i] = [thisKey, "?"]
                    configPairs[i-1][1] = " ".join(thisPair[0:-1])
                else:
                    configPairs[i-1][1] = each
            else:
                thisKey = thisPair[0].strip()
                if thisKey=="tag":
                    haveTag=True
                configPairs[i] = [thisKey, "?"]
        #print(configPairs)
        if haveTag:
            tagList = [eachPair[1][1:-1] for eachPair in configPairs if eachPair[0]=="tag"]
        configs = dict(configPairs)
        if haveTag:
            configs["tag"] = tagList
        if "font" in configs.keys():
            values = configs["font"][1:-1].split()
            assert len(values)==4
            configs["font"] = {self.defaultFontKeys[i]:eachValue for i,eachValue in enumerate(values)}
        if "text" in configs.keys():
            configs["text"] = configs["text"][1:-1]
        if command=="global":
            return command, configs

        return command, params, configs
    def genRegionCommand(self, command, params, configs={}, **kwargs):
        myParams = [str(each) for each in params]
        front = "{command}({params})".format(command=command, params = ",".join(myParams))
        c = {}
        c.update({"tag":"pyds9"})
        c.update(configs)
        if kwargs is not None:
            c.update(kwargs)
        fontKeys = list(set(self.defaultFontKeys).intersection(set(c.keys())))
        if len(fontKeys)>0 or "font" in c.keys():
            if "font" in c.keys():
                if type(c["font"])==str or type(c["font"])==unicode:
                    thisFontConfig = dict([(eachKey, eachValue)
                        for eachKey, eachValue in
                             zip(self.defaultFontKeys, c["font"].split())])
                elif type(c["font"])==dict:
                    thisFontConfig = copy.copy(c["font"])
            else:
                thisFontConfig = copy.copy(self.defaultFont)
            for eachKey in fontKeys:
                thisFontConfig[eachKey] = c[eachKey]
            #print(c)
            #print(thisFontConfig)
            #print(self.defaultFontKeys)
            thisFontStr = " ".join([thisFontConfig[eachkey] for eachkey in self.defaultFontKeys])
            c["font"] = "'{}'".format(thisFontStr)
            for eachKey in fontKeys:
                c.pop(eachKey)
        if "tag" in c.keys():
            if type(c["tag"])==str or type(c["tag"])==unicode:
                aux = ["'{}'".format(each) for each in c["tag"].split(",")]
                c["tag"] = aux
            elif type(c["tag"])==list or type(c["tag"])==tuple:
                aux = ["'{}'".format(each) for each in c["tag"]]
                c["tag"] = aux
            else:
                raise Exception("known type of tag: {}".format(c["tag"]))
        if "text" in c.keys():
            c["text"] = "'{}'".format(c["text"])
        end = self.dict2eq(c)
        return "regions command \"{sys};{front} # {end}\"".format(front=front, end=end, sys=self.rsys)
    def addRegion(self, command, params, configs={}, **kwargs):
        c = copy.copy(configs)
        c.update(kwargs)
        command = self.genRegionCommand(command, params, c)
        print(command)
        try:
            self.set(command)
        except:
            print("error command:", command)
    def _resampleLine(self, params):
        x1,y1,x2,y2,width = params
        x1 = float(x1); y1 = float(y1);
        x2 = float(x2); y2 = float(y2);
        width = float(width)
        if self.rsys == "wcs":
            w = WCS.WCS(self.header)
            ((x1,y1), (x2,y2))= w.wcs_world2pix([[x1, y1], [x2, y2]], 1)
            print('use wcs: {} ==> {}'.format(params, [x1,y1, x2,y2, width]))
        p1 = np.array([x1,y1], dtype=np.float64)
        p2 = np.array([x2,y2], dtype=np.float64)
        dp = p2-p1
        length = np.sqrt(np.sum((p2-p1)**2))
        dpN = dp/length
        if width<2:
            N = 1
        else:
            N = int(width)
        idpN = np.array([-dpN[1], dpN[0]])
        iLength = int(length)+1
        #step = np.linspace(0, 1, iLength)
        step = np.arange(iLength)
        result = np.zeros((iLength,N), dtype=np.dtype([("x", np.float64), ("y", np.float64)]))
        iresult = np.zeros((iLength,N), dtype=np.dtype([("x", np.int64), ("y", np.int64)]))
        for i in range(N):
            xx = x1 + step * dpN[0]
            yy = y1 + step * dpN[1]
            result["x"][:,i] = xx
            result["y"][:,i] = yy
            iresult["x"][:,i] = np.round(xx)-1
            iresult["y"][:,i] = np.round(yy)-1
            x1 += idpN[0]
            y1 += idpN[1]
        return iresult, result# result[:,0]

#==> regions

def getDimStr(header):
    extName = header.get("extname", "")
    dim = header.get("NAXIS", 0)
    dims = ", ".join([str(header.get("NAXIS"+str(each),"?")) for each in range(1, dim+1)])
    obj = header.get("object", "")
    if not dim:
        return "{} {} {}".format(extName, dim, obj)
    else:
        return "{} {} [{}] {}".format(extName, dim, dims, obj)

#<== help strings
h="""Use:
    d.hb    some basic commands (about edit, view, zoom .. etc)
    d.hr    commands about regions and plots
    d.hc    documents about the control panel of 1D array (you must install pyqt4)
    d.hpro  commands to gen projections
    d.hf    commands to do operations for frames
"""
hb=""" in the following command list, "A==>B==>C" means this command is equivalent to the menu "A", sub menu "B" and subsub menu "C" in the ds9 window.

# xpa
    d.set(command), d.get(command): exec ds9 xpa command directly
# switch commands: this commands will have no output, only do switching job
    d.blink_on, d.blink_off: turn on and turn off blinking
    d.top: bring the ds9 window to top and focus on it
    d.single: 'frame==>single'
    d.tile:   'frame==>tile'
    d.minmax:  'scale==>minmax'
    d.zscale:  'scale==>zscale'
    d.match.?: 'frame==>match==>?'
    d.lock.?:  'frame==>lock==>?'
    d.zoomfit: set zoom to fit

# getter commands: type the command in ipython console and get results
    d.newFrame: open a new frame
    d.frames: all existing frames
    d.header: header string for current frame
    d.data:   data for current frame (could be 2d array, 1d array or numpy table)
    d.info:   list information about pan, rotate, zoom, window, scale bin crop, srosshair, mode and orientation.
    d.blink:  the blinking status
    d.shape:  d.data.shape
    d.pyfits: return the astropy.io.fits.hdu.image.PrimaryHDU of current frame
    d.xy:     the mouse pointer because a bigdot, click on image to get [x, y]
    d.kxy:    the mouse pointer because a bigdot, press a key with mouse on image to get [key, x, y]
    d.file    get filename of current frame

# getter and setter commands: you can use "command" to get value and "command=valeu" to set value, the 'getter' and 'setter' will tell you what you get and change using this command
    d.mode:   getter&setter: the mode in "Edit" menu, value could be none, region, crosshair, colorbar, pan, zoom, rotate, crop, catlog and exam
    d.frame:  getter&setter: the number of current frame
    d.blink:  getter&setter: the blinking interval in seconds
    d.rotate: getter&setter: the rotation degree of current frame
    d.zoom:   getter&setter: the zoom factor of current frame
    d.bin:    getter&setter: the x and y bin size of current frame
    d.pan:    getter: the (x, y) center position of current frame
              setter: the (delta_x, delta_y) you want to pan
    d.panto:  getter: the (x, y) center position of current frame
              setter: the (x, y) position you want to pan to
    d.nan:    getter&setter: the color of the nan pixel
    d.limits: getter&setter: the scale limits value
    d.scale:  getter&setter: the scale method, could be linear, log, power, sqrt, squared, asinh, sinh, histogram
    d.zc:     getter&setter: the zscale contrast
    d.crop:   getter&setter: (cropCenter_x, cropCenter_y, cropSize_x, cropSize_y)
    d.cropnon:getter: set no crop
    d.cropbox:getter&setter: (bottomLeft_x, bottomLeft, upperRight_x, upperRight_y)
    d.wcs_format: getter&setter: wcs format, can be degrees or sexagesimal
    d.width:  getter&setter: ds9 window width
    d.height: getter&setter: ds9 window height
    d.window: getter&setter: (width, height)
    d.fo    : getter: list all frames
              setter: close all other frames except for this frames
    d.c     : getter: get all comments in *.fits.comment
              setter: add comments for *.fits[extname] to *.fits.comment
    d.co    : getter: get all comments in *.fits.comment
              setter: set unique comment for *.fits[extname] to *.fits.comment
# functions
    d.ndarrayList(l): give a list of numpy 2d array to show them in new frames
"""

hr=""" # commands about regions: with getter and setter
    d.rsys:     getter&setter: current region system, could be image, physical and wcs
    d.rcolor:   getter&setter: current region color, could be Black, White, Red, Green, Blue, Cyan, Magenta and Yellow
# commands about regions: only with getter
    d.region:   region string for current frame
    d.regions:  region info as a dict for current frame
    d.sregion:  selected region string for current frame
    d.sregions: info of selected region as a dict for current frame
# commands about regions: do operations on region
    d.ms:          give a new tag for all selected regions in this frame
    d.si:          show info about selected regions (as string in a list)
    d.sa:          select all regions in this frame
    d.sat:         delete all tagged regions in this frame
    d.ct:          setter: copy selected regions to frame(s)
    d.ca:          getter: get names for all unique tagged regions in this frame
                   setter: copy all tagged regions to other frame(s)
    d.us:          getter: try to find regions that have same tag name of the selected regions in other frame
                   setter: update position of regions that have same tag name of the selected regions in other frames
    d.usa:         d.usa=True is equivalent to d.us=d.frames
    d.saut:        same as d.ca, but also select them
    d.sat:         select all tagged regions is this frame (some regions may have the same tag name)
    d.ds:          d.ds=True will delete all selected regions in this frame
    d.cleanRegion: delete all regions in this frame
# commands about plot:
    d.ps:   plot or update selected projections
    d.usp:  getter: same as d.us (in help of region)
            setter: same as d.us, and also plot these projections
    d.usap: setter: d.usap=True is equivalent to d.usp=d.frames
"""
hpro="""
# commands to generate projections in ds9:
    d.p:  plot a large projections with tag in the middle of this frame
    d.ph: plot a large projections horizontally
    d.pv: plot a large projections vertically
    d.phi:select two points in ds9 window and plot a horizontal projection within the xlimits and use the y average
    d.pvi:select two points in ds9 window and plot a vertical projection within the ylimits and use the x average
    d.pkci:
        hit two 'h' in the image to gen a horizontal projection (like d.phi)
        hit two 'v' in the image to gen a horizontal projection (like d.phi)
        hit 'q' to return
    d.lx(values): draw horizontal lines at values, values cat be a value or a list
    d.ly(values): draw vertical lines at values, values cat be a value or a list
"""
hc="""
# control panel (you must install pyqt4)
    x limits:
        manually set with "start" and "end"
        use U button to update, use T button to use tight x limits
    y limits:
        there are four method:
            auto: default  ylimits set by metplotlib
            updown: manually set up and down
            updownRej: manually set up and down, or add additional margin to the auto results
            percentile: set ylimits based on percentile
        use the U button to update ylimits
    plots status:
        three select field (this status are used later)
            V means visual
            A means active
            C means focus
        use X button to delete this plot
    keyboard shortcut in the plot window:
        e: open legend
        d: open free cursor
        c: open a cursor attached to the "focus" plot
    change data or plot scale for all "active" plot
        y: how to change the value. Examples:
            "y/y.max()" will normalize all plots
            "x" will make all plots become "y=x"
        style: change style of the plots. Examples:
            "color=blue" will change the color to blue
            "color=blue,linestyle=dashed" will alse change the linestyle to dashed
        scatter: change style for makers. Examples:
            "s=10"
            "marker=x"
            "edgecolor=none"
"""
hf=""" # calculation between frames:
    d.minus:    getter: list all frames
                setter:
                    d.minus=1,2 will calculate {frame1} - {frame2} and put the results in a new frame
                    d.minus=1,2,5 will calculate {frame1} - {frame2} and put the results into frame 5
    d.math:     getter: get help string
                setter:
                    d.math='{1}-{2}*1.1' will calculate {frame1} - {frame2}*1.1 and put the result in a new frame
                    d.math='{1}-{2}*1.1',3 will calculate {frame1} - {frame2}*1.1 and put the result in frame3
    d.vcon:    get all other images that have same x shape
    d.hcon:    get all other images that have same y shape
    d.vconDo:  d.vconDo(frames, s=2, newFrame=None) will vertically concatenate frames, with a margin size of 2
    d.vconDo:  d.hconDo(frames, s=2, newFrame=None) will vertically concatenate frames, with a margin size of 2
"""
#==> help strings

class ds10(DS9, Region):#<==
    """status get and set
    --------------------
    """#==>
    def __init__(self, para=None, figShow=False, ds9Show=True, qt=True, configs={}, **kwargs):
        if ds9Show:
            if para is not None and para:
                super(ds10, self).__init__(para)
            else:
                super(ds10, self).__init__()
        self._data = None
        #!! mode setter
        self.m = ModeSetter(self)
        #!! lock setter
        self.lock = LockSetter(self)
        #!! match setter
        self.match = MatchSetter(self)

        self.configs = configs
        self.configs.update(kwargs)
        # Default value
        self.rconfig={"format":"ds9", "system":"image", 'color':'Cyan'}
        self.defaultFont = {"fontname":"times", "size":"12", "bold":"normal", "italic":"roman"}
        self.defaultFontKeys = ["fontname", "size", "bold", "italic"]
        self._plotType = 'plot'

        #!! my startup configs
        width = kwargs.get("width")
        height = kwargs.get("height")
        if width is not None:
            self.width = width
        if height is not None:
            self.height = height

        if ds9Show:
            #self.set("cmap Heat")
            self.set("lock scale")
            self.set("lock colorbar")
            self.set("region shape projection")
            #self.tile
            self.zscale
            self.lock.frame.image
            self.nan="red"
            #self.m.none
            self.rsys = 'image'
            self.rcolor = 'magenta'

        # plot
        #self._initMarkIndex()
        self._frameMarkIndex = {}
        self.goodfigure = GoodFigure(name='1', parent=self, show=figShow, qt=qt)
        self.plots = self.goodfigure.plots
    #<== only get
    @property
    def frames(self):
        return list(map(int, self.get("frame all").split()))
    @property
    def info(self):
        print("cd: {}".format(self.get("cd")))
        print("file: {}".format(self.get("file")))
        print("frame: {}".format(self.get("frame")))
        print("pan: {} | ".format(self.get("pan")),end="")
        print("rotate: {} | ".format(self.get("rotate")),end="")
        print("zoom: {} {} | ".format(self.get("zoom x"), self.get("zoom y")),end="")
        print("windows: {} {}".format(self.get("width"), self.get("height")))
        print("scale: {} | ".format(self.get("scale")),end="")
        print("bin: {} | ".format(self.get("bin factor")),end="")
        print("crop: {} | ".format(self.get("crop")),end="")
        print("crosshair: {}".format(self.get("crosshair")),)
        print("mode: {}: | ".format(self.get("mode")),end="")
        print("nan: {} | ".format(self.get("nan")),end="")
        print("orient: {} | ".format(self.get("orient")),end="")
    @property
    def file(self):
      return self.get('file')
    @property
    def dfile(self):
      f = self.get('file')
      if '[' in f:
        filename = ''.join(f.split('[')[:-1])
        extname = f.split('[')[-1]
        if not extname.endswith(']'):
          raise Exception('bad filename: {}'.format(f))
        extname = extname[:-1]
        return {'filename': filename, 'extname': extname}
      else:
        return {'filename': filename, 'extname': None}
    #==>
    #<== set and get
    @property
    def c(self):
      filename = self.dfile['filename']
      comment_name = filename + '.comment'
      if os.path.exists(comment_name):
        with open(comment_name) as f:
          result = f.readlines()
        return result
      else:
        return None
    @c.setter
    def c(self, value):
      dfile = self.dfile
      filename = dfile['filename']
      extname = dfile['extname']
      comment_name = filename + '.comment'
      with open(comment_name, 'a') as f:
        exists = f.readlines()
        if not len(exists):
          f.write('{} ==> {}'.format(extname, value))
        else:
          f.write('\n{} ==> {}'.format(extname, value))
    @property
    def co(self):
      return self.c
    @co.setter
    def co(self, value):
      dfile = self.dfile
      filename = dfile['filename']
      extname = dfile['extname']
      comment_name = filename + '.comment'
      if os.path.exists(comment_name):
        with open(comment_name, 'r') as f:
          exists = f.readlines()
        with open(comment_name, 'w') as f:
          others = list(filter(lambda _:not _.startswith('{} ==> '.format(extname)), exists))
          others.append('{} ==> {}'.format(extname, value))
          f.write('\n'.join(others))
      else:
        with open(comment_name, 'w') as f:
          f.write('{} ==> {}'.format(extname, value))

    @property
    def mode(self):
        return self.get("mode")
    @mode.setter
    def mode(self, value):
        self.set("mode {}".format(value))
    @property
    def frame(self):
        return int(self.get("frame"))
    @frame.setter
    def frame(self, value):
        try:
            self.set("frame {}".format(int(value)))
        except Exception as e:
            print(e)
            raise Exception("invaild value or data type for {}: {}, type {}".\
                    format(sys._getframe().f_code.co_name, value, type(value)))
    @property
    def rotate(self):
        return float(self.get("rotate"))
    @rotate.setter
    def rotate(self, value):
        try:
            self.set("rotate {}".format(float(value)-self.rotate))
        except Exception as e:
            print(e)
            raise Exception("invaild value or data type for {}: {}, type {}".\
                    format(sys._getframe().f_code.co_name, value, type(value)))
    @property
    def blink(self):
        result = self.get("blink")
        if result=="yes":
            return 1
        elif result=="no":
            return 0
    @blink.setter
    def blink(self, value):
        assert value, isinstance(value, int) or isinstance(value, float)
        result = self.set("blink interval {}".format(value))
    @property
    def blink_on(self):
        self.set("blink yes")
    @property
    def blink_off(self):
        self.set("blink no")

    @property
    def zoom(self):
        aux = self.get("zoom")
        if aux:
            aux = aux.split()
            if len(aux)==1:
                return float(aux[0])
            else:
                return float(aux[0]), float(aux[1])
        else:
            return None
    @zoom.setter
    def zoom(self, value):
        try:
            if isinstance(value, int) or isinstance(value, float):
                self.set("zoom to {}".format(value))
            else:
                self.set("zoom to {} {}".format(value[0], value[1]))
        except Exception as e:
            print(e)
            raise Exception("invaild value or data type for {}: {}, type {}".\
                    format(sys._getframe().f_code.co_name, value, type(value)))
    @property
    def zoomfit(self):
        self.set("zoom to fit")
    @property
    def bin(self):
        aux = self.get("bin factor")
        if aux:
            aux = aux.split()
            if len(aux)==1:
                return float(aux[0])
            else:
                return float(aux[0]), float(aux[1])
        else:
            return None
    @bin.setter
    def bin(self, value):
        try:
            if isinstance(value, int) or isinstance(value, float):
                self.set("bin factor {}".format(value))
            else:
                self.set("bin factor {} {}".format(value[0], value[1]))
        except Exception as e:
            print(e)
            raise Exception("invaild value or data type for {}: {}, type {}".\
                    format(sys._getframe().f_code.co_name, value, type(value)))
    @property
    def pan(self):
        temp = self.get("pan").split()
        return np.array([float(each) for each in temp])
    @pan.setter
    def pan(self, value):
        current = self.pan
        if len(value)==2:
            self.set("pan {} {}".\
                    format(float(value[0]), float(value[1])))
        else:
            raise Exception("invaild value or data type for {}: {}, type {}".\
                format(sys._getframe().f_code.co_name, value, type(value)))
    @property
    def panto(self):
        temp = self.get("pan").split()
        return np.array([float(each) for each in temp])
    @pan.setter
    def panto(self, value):
        current = self.pan
        if len(value)==2:
            try:
                self.set("pan {} {}".\
                        format(float(value[0])-current[0], float(value[1])-current[1]))
            except:
                raise Exception("invaild value or data type for {}: {}, type {}".\
                    format(sys._getframe().f_code.co_name, value, type(value)))
        elif len(value)==3:
                self.set("pan to {} {} {}".format(float(value[0]),
                                               float(value[1]), value[2]))
        elif len(value)==4:
                self.set("pan to {} {} {} {}".format(float(value[0]),
                                                  float(value[1]), value[2], value[3]))
    @property
    def nan(self):
        return self.get("nan")
    @nan.setter
    def nan(self, value):
        self.set("nan {}".format(value))
    @property
    def limits(self):
        temp = self.get("scale limits").split()
        return np.array([float(each) for each in temp])
    @limits.setter
    def limits(self, value):
        current = [int(each) for each in self.limits]
        if len(value)==2:
            try:
                float(value[0])
                float(value[1])
                self.set("scale limits {} {}".\
                        format(value[0], value[1]))
                return
            except:
                raise Exception("invaild value or data type for {}: {}, type {}".\
                    format(sys._getframe().f_code.co_name, value, type(value)))
    @property
    def scale(self):
        return self.get("scale")
    @scale.setter
    def scale(self, value):
        valueList = ["linear","log","pow","sqrt","squared","asinh","sinh","histequ"]
        if value in valueList:
            self.set("scale {}".format(value))
            return
        raise Exception("invaild value or data type for {}: {}, type {}\n\tvalue must be a string in {}".\
                format(sys._getframe().f_code.co_name, value, type(value), valueList))
    @property
    def zc(self):
        return float(self.get('zscale contrast'))
    @zc.setter
    def zc(self, value):
        self.set('zscale contrast {}'.format(value))
    @property
    def cropnone(self):
        x,y = self.shape
        self.cropbox=1,1,x,y
    @property
    def crop(self):
        result = self.get("crop").split()
        return [float(each) for each in result]
    @crop.setter
    def crop(self, value):
        if value is None:
            self.set("crop reset")
            return
        elif len(value) == 4:
            try:
                int(value[0]);float(value[0])
                int(value[1]);float(value[1])
                int(value[2]);float(value[2])
                int(value[3]);float(value[3])
                self.set("crop {} {} {} {}".format(value[0], value[1], value[2], value[3]))
                return
            except Exception as e:
                print(e)
        raise Exception("invaild value or data type for {}: {}, type {}".\
                format(sys._getframe().f_code.co_name, value, type(value)))
    @property
    def cropbox(self):
        result = self.get("crop").split()
        result = list(map(float, result))
        x, y, dx, dy = result
        return [x-(dx-1)/2, y-(dy-1)/2, x+(dx-1)/2, y+(dy-1)/2]
    @cropbox.setter
    def cropbox(self, value):
        if value is None:
            self.set("crop reset")
            return
        elif len(value) == 4:
            try:
                int(value[0]);float(value[0])
                int(value[1]);float(value[1])
                int(value[2]);float(value[2])
                int(value[3]);float(value[3])
                #self.set("crop {} {} {} {}".format(value[0], value[1], value[2], value[3]))
                # x1, y1, x2, y2
                # y is the first index in numpy
                self.set("crop {} {} {} {}".
                            format( (value[2]+value[0])/2, (value[3]+value[1])/2,
                                    (value[2]-value[0]+1), (value[3]-value[1]+1) ))
                return
            except Exception as e:
                print(e)
        raise Exception("invaild value or data type for {}: {}, type {}".\
                format(sys._getframe().f_code.co_name, value, type(value)))
    @property
    def width(self):
        return float(self.get("width"))
    @width.setter
    def width(self, value):
        assert value, isinstance(value, int) or isinstance(value, float)
        self.set("width {}".format(value))
    @property
    def height(self):
        return float(self.get("height"))
    @height.setter
    def height(self, value):
        assert value, isinstance(value, int) or isinstance(value, float)
        self.set("height {}".format(value))
    @property
    def window(self):
        height = self.height
        width = self.width
        return width, height
    @window.setter
    def window(self, value):
        assert type(value) in [list, tuple]
        if len(value)==2:
            try:
                self.width = int(value[0])
                self.height = int(value[1])
                return
            except Exception as e:
                print(e)
        raise Exception("invaild value or data type for {}: {}, type {}".\
                    format(sys._getframe().f_code.co_name, value, type(value)))
    @property
    def wcs_format(self):
        return self.get('wcs format')
    @wcs_format.setter
    def wcs_format(self, value):
        if value not in ['degrees', 'sexagesimal']:
            raise ValueError('wcs format can only be degrees or sexagesimal')
        else:
            self.set('wcs format {}'.format(value))

    #==> set and get
    #<== only set
    @property
    def top(self):
        self.set("raise")
    @property
    def single(self):
        self.set("single")
    @property
    def tile(self):
        self.set("tile")
    @property
    def zscale(self):
        self.set("zscale")
    @property
    def minmax(self):
        self.set("minmax")
    #==> only set
    #<== about data
    @property
    def shape(self):
        "return shape of (x, y), note that the readin pyfits data is of shape (y,x)"
        data = self.get_arr2np()
        return data.T.shape
    @property
    def data(self):
        self._thisData = self.pyfits.data
        return self._thisData
    @property
    def pyfits(self):
        self._thisPyfits = self.get_fits()[0]
        return self._thisPyfits
    @property
    def header(self):
        self._thisHeader=self.get_fits()[0].header
        return self._thisHeader
    def fits(self, name, frame=None):
        "open a file in new frame"
        try:#!! test if there has no data frame yet
            data = self.data
        except:
            self.set("frame delete")
        frames = self.frames
        aux = [int(each) for each in frames]
        if len(aux)==0:
            startNumber = 1
        else:
            startNumber = np.array(aux).max()+1
        print("add new frame: {}".format(startNumber))
        self.frame = startNumber
        if frame is not None:
            thisName = "{}[{}]".format(name, frame)
            print("\t", thisName, " ", getDimStr(pyfits.open(name)[frame].header))
            self.set("fits \"{}\"".format(os.path.abspath(thisName)))
        else:
            self._data = pyfits.open(name)
            N = len(self._data)
            for i in range(N):
                thisName = "{}[{}]".format(name, i)
                print("\t", thisName, "", getDimStr(pyfits.open(name)[i].header))
                self.set("fits \"{}\"".format(thisName))
        return startNumber
    @property
    def newFrame(self):
        try:
            data = self.data
        except:
            self.set("frame delete")
        frames = self.frames
        aux = [int(each) for each in frames]
        if len(aux)==0:
            startNumber = 1
        else:
            startNumber = np.array(aux).max()+1
        print("add new frame: {}".format(startNumber))
        self.frame = startNumber
        return startNumber
    @property
    def xy(self):
        result = self.get("iexam coordinate image").split()
        return [float(result[0]), float(result[1])]
    @property
    def kxy(self):
        result = self.get("iexam key").split()
        return [result[0], float(result[1]), float(result[2])]
    def ndarrayList(self, l):
        if isinstance(l, list):
            for eachArray in l:
                try:#!! test if there has no data frame yet
                    data = self.data
                except:
                    self.set("frame delete")
                frames = self.frames
                aux = [int(each) for each in frames]
                if len(aux)==0:
                    startNumber = 1
                else:
                    startNumber = np.array(aux).max()+1
                print("add new frame: {}".format(startNumber))
                self.frame = startNumber
                self.set_np2arr(eachArray)
        else:
            raise Exception("please input a list of ndarray")

    #==> about data
    #<== about region
    @property
    def rsys(self):
        return self.rconfig["system"]
    @rsys.setter
    def rsys(self, value):
        if value not in ['wcs', 'image', 'physical']:
            raise ValueError('region system should be wcs, image or physical')
        else:
            self.rconfig["system"] = value
            self.set('regions system {}'.format(value))
    @property
    def rcolor(self):
        return self.rconfig["color"]
    @rcolor.setter
    def rcolor(self, value):
        clist = ['black', 'white', 'red', 'green', 'blue', 'cyan', 'magenta', 'yellow']
        if value not in clist:
            raise ValueError('color should be in {}'.format(clist))
        self.rconfig["color"] = value

    @property
    def rformat(self):
        return self.get('regions format')
    @property
    def region(self):
        'information about all regions'
        self.set("regions format " + self.rconfig["format"])
        self.set("regions system " + self.rconfig["system"])
        self._thisRegion = self.get("regions list").split("\n")
        temp = [each for each in self._thisRegion if each]
        self._thisRegion = temp
        self._thisTags = set()
        self._thisTagCount = {}
        self._thisTagObjs = {}
        for each in self._thisRegion:
            if "tag" in each:
                temp = re_tag.findall(each)
                for eachResult in temp:
                    if eachResult in self._thisTags:
                        self._thisTagCount[eachResult] += 1
                        self._thisTagObjs[eachResult].append(each)
                    else:
                        self._thisTags.add(eachResult)
                        self._thisTagCount[eachResult] = 1
                        self._thisTagObjs[eachResult] = [each]
        self._thisUniqueTags = [each for each in self._thisTags
                                    if self._thisTagCount[each]==1]
        return self._thisRegion
    @property
    def regions(self):
        'region info in dict'
        self.set("regions format " + self.rconfig["format"])
        self.set("regions system " + self.rconfig["system"])
        self._thisRegions = {}
        for eachType in regionTypes:
            self._thisRegions[eachType] = []
        self._thisRegions["other"] = []
        temp = [each for each in self.get("regions list").split("\n") if each]
        for eachTerm in temp:
            flag = False
            for eachType in regionTypes:
                if (eachType+"(") in eachTerm.lower():
                    self._thisRegions[eachType].append(eachTerm)
                    flag = True
                    break
            if not flag:
                self._thisRegions["other"].append(eachTerm)
        temp=[]
        for eachTerm in self._thisRegions["panda"]:
            if "epanda" in eachTerm:
                self._thisRegions["epanda"].append(eachTerm)
            else:
                temp.append(eachTerm)
        self._thisRegions["panda"] = temp
        for eachKey in list(self._thisRegions.keys()):
            if not self._thisRegions[eachKey]:
                self._thisRegions.pop(eachKey)
            else:
                self._thisRegions[eachKey] = \
                        [self.processRegionLine(eachLine)
                                for eachLine in self._thisRegions[eachKey]]

        return self._thisRegions
    @property
    def sregion(self):
        'information about selected regions'
        self.set("regions format " + self.rconfig["format"])
        self.set("regions system " + self.rconfig["system"])
        self._thisRegion = self.get("regions selected").split("\n")
        temp = [each for each in self._thisRegion if each]
        self._thisRegion = temp
        return self._thisRegion
    @property
    def sregions(self):
        self.set("regions format " + self.rconfig["format"])
        self.set("regions system " + self.rconfig["system"])
        self._thisRegions = {}
        for eachType in regionTypes:
            self._thisRegions[eachType] = []
        self._thisRegions["other"] = []
        temp = [each for each in self.get("regions selected").split("\n") if each]
        for eachTerm in temp:
            flag = False
            for eachType in regionTypes:
                if (eachType+"(") in eachTerm.lower():
                    self._thisRegions[eachType].append(eachTerm)
                    flag = True
                    break
            if not flag:
                self._thisRegions["other"].append(eachTerm)
        temp=[]
        for eachTerm in self._thisRegions["panda"]:
            if "epanda" in eachTerm:
                self._thisRegions["epanda"].append(eachTerm)
            else:
                temp.append(eachTerm)
        self._thisRegions["panda"] = temp
        for eachKey in list(self._thisRegions.keys()):
            if not self._thisRegions[eachKey]:
                self._thisRegions.pop(eachKey)
            else:
                self._thisRegions[eachKey] = \
                        [self.processRegionLine(eachLine)
                                for eachLine in self._thisRegions[eachKey]]
        return self._thisRegions
    @property
    def cleanRegion(self):
        self.set("regions deleteall")
    def addProjection(self, params, configs={}, **kwargs):
        "addProjection((y1,x1, y2,x2, width), tag=[], text='')"
        assert len(params)==5
        c = {}; c.update(configs)
        c.update(kwargs)
        command = self.genRegionCommand("projection", params, configs=c)
        try:
            self.set(command)
        except:
            print("error command:", command)
    def addCircle(self, params, configs={}, **kwargs):
        "addCircle((y,x, r), tag=[], text='')"
        assert len(params)==3
        c = {}; c.update(configs)
        c.update(kwargs)
        command = self.genRegionCommand("circle", params, configs=c)
        try:
            self.set(command)
        except:
            print("error command:", command)
    def addText(self, xyText, configs={}, **kwargs):
        "addText((x,y,'text'), tag=[])"
        assert len(xyText)==3
        c = {}; c.update(configs)
        c.update(kwargs)
        c.update({"text":xyText[2]})
        command = self.genRegionCommand("text", xyText[:2], configs=c)
        try:
            self.set(command)
        except:
            print("error command:", command)
    def addPoint(self, params, configs={"point":"X"}, **kwargs):
        "addPoint((x,y), point='X')"
        assert len(params)==2
        c = {}; c.update(configs)
        c.update(kwargs)
        command = self.genRegionCommand("point", params, configs=c)
        try:
            self.set(command)
        except:
            print("error command:", command)
    #==> about region
    #<== some utils
    def markAllPixel(self,size=20):
        "make a rect region full of points"
        data = self.data
        m,n = self.pan
        allDataPair = [(i,j)
                for i in range(int(m-size), int(m+size))
                for j in range(int(n-size), int(n+size))]
        for eachPair in allDataPair:
            self.addPoint(eachPair)
    def getMarkIndex(self, prefix=None):
        frame = self.frame
        if prefix is None:
            if not self._frameMarkIndex.get(frame, {}).get("default", 0):
                if not self._frameMarkIndex.get(frame, {}):
                    self._frameMarkIndex[frame]={}
                self._frameMarkIndex[frame]["default"]=1
                return 1
            else:
                self._frameMarkIndex[frame]["default"] += 1
                return self._frameMarkIndex[frame]["default"]
        else:
            if not self._frameMarkIndex.get(frame, {}).get(prefix, 0):
                if not self._frameMarkIndex.get(frame, {}):
                    self._frameMarkIndex[frame]={}
                self._frameMarkIndex[frame][prefix]=1
                return "{}-{}".format(prefix, 1)
            else:
                self._frameMarkIndex[frame][prefix] += 1
                return "{}-{}".format(prefix, self._frameMarkIndex[frame][prefix])
    def _initMarkIndex(self): # todo
        self.region
        uniqueTags = self._thisUniqueTags
        print(uniqueTags)
    #==> some utils
    #<== gen projection in ds9
    @property
    def p(self):
        'generate a marked projection in the middle of the frame'
        p = self.pan
        frame = self.frame
        zoom = self.zoom
        p1 = p - 70./zoom
        p2 = p + 70./zoom
        self.region
        allTag = self._thisTags
        tagName = "{}:{}".format(frame, self.getMarkIndex())
        while tagName in allTag:
            tagName = "{}:{}".format(frame, self.getMarkIndex())
        self.addProjection((p1[0],p1[1],p2[0],p2[1],0), tag=tagName, text=tagName)
        self.m.region
        self.set("regions group {} select".format(tagName))
    @property
    def ph(self):
        'generate a marked projection in the middle of the frame'
        zoom = self.zoom
        p = self.pan
        frame = self.frame
        p1 = p - 70./zoom
        p1[1] += 70./zoom
        p2 = p + 70./zoom
        p2[1] -= 70./zoom
        self.region
        allTag = self._thisTags
        tagName = "{}:{}".format(frame, self.getMarkIndex())
        while tagName in allTag:
            tagName = "{}:{}".format(frame, self.getMarkIndex())
        self.addProjection((p1[0],p1[1],p2[0],p2[1],0), tag=tagName, text=tagName)
        self.m.region
        self.set("regions group {} select".format(tagName))
    @property
    def pv(self):
        'generate a marked projection in the middle of the frame'
        zoom = self.zoom
        p = self.pan
        frame = self.frame
        p1 = p - 70./zoom
        p1[0] += 70./zoom
        p2 = p + 70./zoom
        p2[0] -= 70./zoom
        self.region
        allTag = self._thisTags
        tagName = "{}:{}".format(frame, self.getMarkIndex())
        while tagName in allTag:
            tagName = "{}:{}".format(frame, self.getMarkIndex())
        self.addProjection((p1[0],p1[1],p2[0],p2[1],0), tag=tagName, text=tagName)
        self.m.region
        self.set("regions group {} select".format(tagName))
    @property
    def phi(self):
        'generate a marked projection in the middle of the frame'
        self.region
        frame = self.frame
        allTag = self._thisTags
        tagName = "{}:{}".format(frame, self.getMarkIndex())
        while tagName in allTag:
            tagName = "{}:{}".format(frame, self.getMarkIndex())
        self.top
        p1 = self.xy
        print("point: {}".format(p1))
        p2 = self.xy
        print("point: {}".format(p2))
        middle = (p1[1]+p2[1])/2
        self.addProjection((p1[0],middle,p2[0],middle,0), tag=tagName, text=tagName)
        self.m.region
        self.set("regions group {} select".format(tagName))
    @property
    def pvi(self):
        'generate a marked projection in the middle of the frame'
        self.region
        frame = self.frame
        allTag = self._thisTags
        tagName = "{}:{}".format(frame, self.getMarkIndex())
        while tagName in allTag:
            tagName = "{}:{}".format(frame, self.getMarkIndex())
        self.top
        p1 = self.xy
        print("point: {}".format(p1))
        p2 = self.xy
        print("point: {}".format(p2))
        middle = (p1[0]+p2[0])/2
        self.addProjection((middle,p1[1],middle,p2[1],0), tag=tagName, text=tagName)
        self.m.region
        self.set("regions group {} select".format(tagName))
    @property
    def pki(self):
        'generate a marked projection in the middle of the frame'
        self.region
        frame = self.frame
        allTag = self._thisTags
        tagName = "{}:{}".format(frame, self.getMarkIndex())
        while tagName in allTag:
            tagName = "{}:{}".format(frame, self.getMarkIndex())
        kDict = {"v":"vertical", "h":"horizontal", "q":"quit"}
        self.top
        p1 = self.kxy
        print("type: {} point: {}".format(kDict.get(p1[0], "unknown"), p1[1:]))
        if p1[0]=="q":
            return 0
        p2 = self.kxy
        print("type: {} point: {}".format(kDict.get(p2[0], "unknown"), p2[1:]))
        if p1[0]==p2[0]=="v":
            middle = (p1[1]+p2[1])/2
            self.addProjection((middle, p1[2], middle, p2[2], 0), tag=tagName, text=tagName)
            self.m.region
            self.set("regions group {} select".format(tagName))
        elif p1[0]==p2[0]=="h":
            middle = (p1[2]+p2[2])/2
            self.addProjection((p1[1], middle, p2[1], middle, 0), tag=tagName, text=tagName)
            self.m.region
            self.set("regions group {} select".format(tagName))
        elif p2[0]=="q":
            return 0
        else:
            print("two key must be the same!!")
    @property
    def pkci(self):
        while True:
            result = self.pki
            if result is not None:
                return

    def lx(self, values, **kwargs):
        "draw horizontal line(s) : lx(xind) or lx([xind1, xind2, ..])"
        thisFrame = self.frame
        c = copy.copy(kwargs)
        xMax, yMax = self.shape
        if isinstance(values, int):
            todo = [values]
        elif isinstance(values, list) or isinstance(values, tuple):
            todo = values
        else:
            raise Exception("input not understood: {}".format(values))
        for eachValue in todo:
            eachTag = "{}:row_{}".format(thisFrame, eachValue)
            c.update({"tag": eachTag})
            self.addProjection((0,eachValue,yMax,eachValue,0), **c)
            self.set("regions group {} select".format(eachTag))
    def ly(self, values, **kwargs):
        "draw vertical lines(s) : ly(yind) or ly([yind1, yind2, ..])"
        thisFrame = self.frame
        c = copy.copy(kwargs)
        xMax, yMax = self.shape
        if isinstance(values, int):
            todo = [values]
        elif isinstance(values, list) or isinstance(values, tuple):
            todo = values
        else:
            raise Exception("input not understood: {}".format(values))
        for eachValue in todo:
            eachTag = "{}:column_{}".format(thisFrame, eachValue)
            c.update({"tag": eachTag})
            c.update({"tag": eachTag})
            self.addProjection((eachValue,0,eachValue,xMax,0), **c)
            self.set("regions group {} select".format(eachTag))
    @property
    def pro(self):
        "show all tags for projections in current frame"
        aux = self.regions.get("projection",[])
        result = [each[2].get("tag",[]) for each in aux]
        return result
    #==> print projection
    #<== select, modify tagged region
    @property
    def ms(self):
        "mark selected: give a tag to the selected projection if it do not have one"
        allRegions = self.region
        allTag = self._thisTags
        allGoodRegions=[]
        for eachRegion in self.sregion:
            aux = self.processRegionLine(eachRegion)
            if isinstance(aux, tuple) or isinstance(aux, list):
                if aux[0]!="global":
                    allGoodRegions.append(aux)
        if not allGoodRegions:
            raise Exception("no select regions!")
        self.set("regions delete select")
        toSelect  = []
        for eachGoodRegion in allGoodRegions:
            regionType, params, configs = eachGoodRegion
            frame = self.frame
            aux = configs.get("tag", [])
            thisUniqueTag = []
            for eachTag in self._thisUniqueTags:
                if eachTag in aux:
                    thisUniqueTag.append(eachTag)
            N=len(thisUniqueTag) # have no tag
            if N==0:
                tagName = "{}:{}".format(frame, self.getMarkIndex())
                while tagName in allTag:
                    tagName = "{}:{}".format(frame, self.getMarkIndex())
                print("new marked region:", regionType, params, configs, tagName)
            else:
                tagName = thisUniqueTag[0]
            toSelect.append(tagName)
            self.addRegion(regionType, params, configs=configs, tag=tagName, text=tagName)
        for eachTag in toSelect:
            self.set("regions group {} select".format(eachTag))
        self.set('regions color {}'.format(self.rcolor))
    @ms.setter
    def ms(self, prefix):
        "mark selected: give a tag to the selected projection if it do not have one"
        allRegions = self.region
        allTag = self._thisTags
        allGoodRegions=[]
        for eachRegion in self.sregion:
            aux = self.processRegionLine(eachRegion)
            if isinstance(aux, tuple) or isinstance(aux, list):
                if aux[0]!="global":
                    allGoodRegions.append(aux)
        if not allGoodRegions:
            raise Exception("no select regions!")
        self.set("regions delete select")
        toSelect  = []
        for eachGoodRegion in allGoodRegions:
            regionType, params, configs = eachGoodRegion
            frame = self.frame
            aux = configs.get("tag", [])
            tagName = "{}:{}".format(frame, self.getMarkIndex(prefix=prefix))
            toSelect.append(tagName)
            print(regionType, params, configs, tagName)
            self.addRegion(regionType, params, configs=configs, tag=tagName, text=tagName)
        for eachTag in toSelect:
            self.set("regions group {} select".format(eachTag))
    @property
    def si(self):
        'show info about the selected projection'
        r = self.sregion
        todo = []
        for each in r:
            if "file format" in each:
                continue
            if each in ["image", "physical", "wcs;"]:
                continue
            if "global" in each:
                continue
            todo.append(each)
        if len(todo)<1:
            raise Exception("no select regions!")
        else:
            return todo

    @property
    def ct(self):
        return None
    @ct.setter
    def ct(self, values):
        "copy selected regions to frame(s)"
        if type(values)==int:
            toFrames = [values]
        else:
            toFrames = values
        thisFrame = self.frame
        todoList = []
        for eachLine in self.sregion:
            aux = self.processRegionLine(eachLine)
            if len(aux) == 3 and isinstance(aux, (tuple, list)):
                todoList.append(aux)
        for toFrame in toFrames:
            self.frame = toFrame
            for eachTodo in todoList:
                regionType, params, configs = eachTodo
                tagName = configs.get("tag",["{}:?".format(thisFrame)])
                tagName = tagName[0]
                self.addRegion(regionType, params, configs=configs, tag=tagName, text=tagName)
        self.frame = thisFrame

    @property
    def us(self):
        'show info about the selected projection'
        allRegions = self.si
        allTags = []
        thisFrame = self.frame
        allFrameTags = {"_current:{}".format(thisFrame):allTags}
        for eachRegion in allRegions:
            aux = self.processRegionLine(eachRegion)
            if len(aux)==3:
                tagName = aux[2]["tag"][0]
                allTags.append(tagName)
        allFrames = self.frames
        allFrames.pop(allFrames.index(thisFrame))
        for eachFrame in allFrames:
            thisTags=[]
            self.frame = eachFrame
            allRegions = self.region
            for eachTag in allTags:
                for eachRegion in allRegions:
                    if "tag={{{}}}".format(eachTag) in eachRegion:
                        thisTags.append(eachTag)
            if len(thisTags)>0:
                allFrameTags.update({eachFrame:thisTags})

        self.frame = thisFrame
        return allFrameTags
    @us.setter
    def us(self, values):
        '''update selectd regions to frame(s)'''
        if type(values)==int:
            toFrames = [values]
        else:
            toFrames = values
        thisFrame = self.frame
        todoList = []
        for eachLine in self.sregion:
            aux = self.processRegionLine(eachLine)
            if len(aux) == 3:
                todoList.append(aux)
        unknownTag = "{}:?".format(thisFrame)
        for toFrame in toFrames:
            self.frame = toFrame
            self.region
            allUniqueTag2 = self._thisTags
            for eachTodo in todoList:
                regionType, params, configs = eachTodo
                tagName = configs.get("tag",[unknownTag])
                tagName = tagName[0]
                if tagName in allUniqueTag2 and tagName!=unknownTag:
                    self.set("regions group {} delete".format(tagName))
                self.addRegion(regionType, params, configs=configs, tag=tagName, text=tagName)
        self.frame = thisFrame
    @property
    def usa(self,confirm):
      return None
    @usa.setter
    def usa(self,confirm):
        '''update selectd regions to all frame(s)'''
        # t = input("do you want to update all selected regions IN ALL FRAMES? [0=no, 1=yes]")
        if str(confirm)==True:
            self.us = self.frames
    @property
    def ca(self):
        'show what to "copy all"'
        self.region
        allTag = self._thisTags
        return allTag
    @ca.setter
    def ca(self, values):
        '''copy all projections with unique tag to frames
        '''
        if type(values)==int:
            toFrames = [values]
        else:
            toFrames = values
        thisFrame = self.frame
        self.region
        allUniqueTag = self._thisTags
        allTodo = []
        for each in allUniqueTag:
            result = self.processRegionLine(self._thisTagObjs[each][0])
            thisTag = result[2]["tag"][0]
            allTodo.append((result[0], result[1], result[2], thisTag))
        print(allTodo)
        for toFrame in toFrames:
            self.frame = toFrame
            self.region
            allUniqueTag2 = self._thisTags
            for eachParas in allTodo:
                regionType, params, configs, tagName = eachParas
                print(tagName)
                if tagName in allUniqueTag2:
                    self.set("regions group {} delete".format(tagName))
                self.addRegion(regionType, params, configs=configs, tag=tagName, text=tagName)
        self.frame = thisFrame
    @property
    def saut(self):
        'select all projections with unique tag in this frame'
        self.region
        allUniqueTag = self._thisUniqueTags
        self.set("regions select none")
        for eachTag in allUniqueTag:
            self.set("regions select group {}".format(eachTag))
        return allUniqueTag
    @property
    def sa(self):
        self.set("regions select all")
    @property
    def sat(self):
        'select all tagged region in this frame'
        self.region
        allUniqueTag = self._thisTags
        self.set("regions select none")
        for eachTag in allUniqueTag:
            self.set("regions select group {}".format(eachTag))
        return allUniqueTag
    @property
    def ds(self):
      return None
    @ds.setter
    def ds(self,confirm):
        'delete selected regions in this frame'
        # t = input("do you want to delete all selected regions? [0=no, 1=yes]")
        if confirm==True:
            self.set("regions delete select")
    #==> select, modify tagged region
    #<== about figure plot
    @property
    def ps(self):
        "plot or update current selected projection"
        self.ms
        self.region
        allTag = self._thisTags
        allInfo = []
        for eachProjection in self.sregions["projection"]:
            try:
                aux, params, configs = eachProjection
            except:
                raise Exception("no select projection!")
            frame = self.frame
            aux = configs.get("tag",[])
            thisUniqueTag = []
            for eachTag in self._thisUniqueTags:
                if eachTag in aux:
                    thisUniqueTag.append(eachTag)
            N=len(thisUniqueTag)
            if N==0:
                tagName = "{}:{}".format(frame, self.getMarkIndex())
                self.set("regions delete select")
                self.addProjection(params, configs=configs, tag=tagName, text=tagName)
                self.set("regions group {} select".format(tagName))
            else:
                tagName = thisUniqueTag[0]

            # calculate
            data = self.data.T
            print('input:', params)
            iIndex, index = self._resampleLine(params)
            try:
                self._idataOri = idataOri = data[iIndex["x"], iIndex["y"]]
            except IndexError:
                raise Exception('please put your projections inside the image')
            self._idata = idata = idataOri.mean(axis=1)
            N = idata.size
            xx = np.arange(N)
            thisName = "{}@{}".format(tagName, frame)
            info = self.goodfigure.addPlot(self._plotType, thisName, xx, idata)
            allInfo.append(info)

            #self.m.region
            #self.set("regions group {} select".format(tagName))
        #return allInfo
    @ps.setter
    def ps(self, value):
        "plot or update current selected projection"
        self.ms = value
        self.region
        allTag = self._thisTags
        allInfo = []
        for eachProjection in self.sregions["projection"]:
            try:
                aux, params, configs = eachProjection
            except:
                raise Exception("no select projection!")
            frame = self.frame
            aux = configs.get("tag",[])
            thisUniqueTag = []
            for eachTag in self._thisUniqueTags:
                if eachTag in aux:
                    thisUniqueTag.append(eachTag)
            N=len(thisUniqueTag)
            if N==0:
                tagName = "{}:{}".format(frame, self.getMarkIndex())
                self.set("regions delete select")
                self.addProjection(params, configs=configs, tag=tagName, text=tagName)
                self.set("regions group {} select".format(tagName))
            else:
                tagName = thisUniqueTag[0]

            # calculate
            data = self.data.T
            iIndex, index = self._resampleLine(params)
            try:
                idataOri = data[iIndex["x"], iIndex["y"]]
            except IndexError:
                raise Exception('please put your projections inside the image')
            idata = idataOri.mean(axis=1)
            N = idata.size
            xx = np.arange(N)
            thisName = "{}@{}".format(tagName, frame)
            info = self.goodfigure.addPlot(self._plotType, thisName, xx, idata)
            allInfo.append(info)

            #self.m.region
            #self.set("regions group {} select".format(tagName))
        #return allInfo
    @property
    def usp(self):
        'get info of projections in all frames'
        self.ms
        allRegions = self.si
        allTags = []
        thisFrame = self.frame
        allFrameTags = {"_current:{}".format(thisFrame):allTags}
        for eachRegion in allRegions:
            aux = self.processRegionLine(eachRegion)
            if len(aux)==3:
                tagName = aux[2]["tag"][0]
                allTags.append(tagName)
        allFrames = self.frames
        allFrames.pop(allFrames.index(thisFrame))
        for eachFrame in allFrames:
            thisTags=[]
            self.frame = eachFrame
            allRegions = self.region
            for eachTag in allTags:
                for eachRegion in allRegions:
                    if "tag={{{}}}".format(eachTag) in eachRegion:
                        thisTags.append(eachTag)
            if len(thisTags)>0:
                allFrameTags.update({eachFrame:thisTags})

        self.frame = thisFrame
        return allFrameTags
    @usp.setter
    def usp(self, values):
        'update selected projection (in frames?) in figure'
        if type(values)==int:
            toFrames = [values]
        else:
            toFrames = values
        thisFrame = self.frame
        self.ps
        todoList = []
        for eachLine in self.sregion:
            aux = self.processRegionLine(eachLine)
            if len(aux) == 3:
                todoList.append(aux)
        unknownTag = "{}:?".format(thisFrame)
        for toFrame in toFrames:
            self.frame = toFrame
            self.region
            allUniqueTag2 = self._thisTags
            for eachTodo in todoList:
                regionType, params, configs = eachTodo
                tagName = configs.get("tag",[unknownTag])
                tagName = tagName[0]
                if tagName in allUniqueTag2 and tagName!=unknownTag:
                    self.set("regions group {} delete".format(tagName))
                self.addRegion(regionType, params, configs=configs, tag=tagName, text=tagName)
                self.set("regions group {} select".format(tagName))
                self.ps
        self.frame = thisFrame
    @property
    def usap(self):
      return None
    @usap.setter
    def usap(self,confirm):
        '''update selectd regions to all frame(s) and update their figure'''
        # t = input("do you want to update all selected regions and their FIGURE IN ALL FRAMES? [0=no, 1=yes]")
        if confirm == True:
            self.usp = self.frames
    @property
    def ps_step(self):
        try:
            self._plotType='step'
            self.ps
        finally:
            self._plotType='plot'
    @ps_step.setter
    def ps_step(self, value):
        try:
            self._plotType='step'
            self.ps = value
        finally:
            self._plotType='plot'
    @property
    def usp_step(self):
        try:
            self._plotType='step'
            self.usp
        finally:
            self._plotType='plot'
    @usp_step.setter
    def usp_step(self, values):
        try:
            self._plotType='step'
            self.usp = values
        finally:
            self._plotType='plot'
    @property
    def usap_step(self):
        try:
            self._plotType='step'
            self.usap=True
        finally:
            self._plotType='plot'

    def save(self, filename):
        with open(filename, "wb") as f:
            pickle.dump(self.goodfigure.fig, f)

    #==> about figure plot
    #<== about fitting
    @property
    def gaussianFit(self):
        """fit a gaussian in figure plot"""
        #datas = self.goodfigure.fig.ginput(2)
        focus = self.goodfigure.plots[self.goodfigure._focus]["buttons"]["focus"]
        self.goodfigure.dynamic.openGetData(2)
        temp = input("please click twice to get x boundary and type Enter\n")
        datas = self.goodfigure.dynamic.xyData
        if len(datas)<2:
            raise Exception("you need to click twice to get the fitting boundary")
        xLow = datas[0][0]
        xHigh = datas[1][0]
        allx = self.goodfigure.dynamic.x
        ally = self.goodfigure.dynamic.y
        index = np.logical_and(xLow<allx, allx<xHigh)
        x = allx[index]
        y = ally[index]
        popt, pcov = curve_fit(gaussian, x, y, p0=(1, x.mean(), 1, 0))
        xx = np.linspace(allx.min(), allx.max(), len(allx)*10)
        yy = gaussian(xx, *popt)
        print("A: {}, mu: {}, sigma: {}, C: {}\n".format(*popt))
        self.goodfigure.addPlot("plot", "fitting" , xx, yy, color = "red")
        focus.toggle()
    #==> about fitting
    #<== todo
    @property
    def slidep(self): # equal to reset
        "todo..."
        self.region
        allTag = self._thisTags
        allInfo = []
        for eachProjection in self.sregions["projection"]:
            try:
                aux, params, configs = eachProjection
            except:
                raise Exception("no select projection!")
            frame = self.frame
            aux = configs.get("tag",[])
            thisUniqueTag = []
            for eachTag in self._thisUniqueTags:
                if eachTag in aux:
                    thisUniqueTag.append(eachTag)
            N=len(thisUniqueTag)
            if N==0:
                tagName = "{}:{}".format(frame, self.getMarkIndex())
                self.set("regions delete select")
                configs.update({"tag":tagName, "text": tagName})
                self.addProjection(params, configs=configs)
                self.set("regions group {} select".format(tagName))
            else:
                tagName = thisUniqueTag[0]

            xy = [float(each) for each in params[0:4]]
            p = np.array([(xy[0]+xy[2])/2., (xy[1]+xy[3])/2.])
            zoom = self.zoom
            pStart = p - 70./zoom
            pEnd = p + 70./zoom

            startTag = "{}:start".format(tagName)
            self.set("region select none")
            self.set("region group {} select".format(startTag))
            if len(self.sregions["circle"])==0:
                self.addPoint(pStart[0], pStart[1], tag=[startTag, tagName], text=startTag)
            else:
                aux, points, aux1 = self.sregions["circle"][-1]
                pStart = np.array([float(points[0]), float(points[1])])

            endTag = "{}:end".format(tagName)
            self.set("region select none")
            self.set("region group {} select".format(endTag))
            if len(self.sregions["circle"])==0:
                self.addPoint(pEnd[0], pEnd[1], tag=[endTag, tagName], text=endTag)
            else:
                aux, points, aux1 = self.sregions["circle"][-1]
                pEnd = np.array([float(points[0]), float(points[1])])

            # calculate
            data = self.data.T
            iIndex, index = self._resampleLine(params)
            try:
                idataOri = data[iIndex["x"], iIndex["y"]]
            except IndexError:
                raise Exception('please put your projections inside the image')
            idata = idataOri.mean(axis=1)
            N = idata.size
            xx = np.arange(N)
            thisName = "{}@{}".format(tagName, frame)
            info = self.goodfigure.addSlidePlot("plot", thisName, pStart, pEnd, [1,1], xx, idata)
            allInfo.append(info)

            #self.m.region
            self.set("regions group {} select".format(tagName))
        return allInfo

    def printProjection(self, params):
        data = self.data
        iIndex, index = self._resampleLine(params)
        datas = data[iIndex["x"], iIndex["y"]].mean(axis=1)
        N = data.size
        xx = np.arange(N)
    #==> todo
    #<== math operation between frames
    @property
    def minus(self):
        return self.frames
    @minus.setter
    def minus(self, values):
        if len(values)==2:
            A=values[0]
            B=values[1]
            C=None
        elif len(values)==3:
            A=values[0]
            B=values[1]
            C=values[2]
        else:
            raise Exception("input shoud be tuple of 2 or 3 frame")
        allFrames = self.frames
        for each in [A, B]:
            assert (each in allFrames), each
        self.frame = A
        dataA = self.data
        self.frame = B
        dataB = self.data
        dataC = dataA - dataB
        print("shape:", dataC.shape)
        if C is None:
            self.set("frame new")
        else:
            self.frame = C
        newFrame = self.frame
        thisObj = pyfits.PrimaryHDU()
        thisObj.data = dataC
        thisObj.header["EXTNAME"] = "[{}]-[{}]".format(A,B)
        thisObj.header["INHERIT"] = False
        thisObj.header["OBJECT"] = "ds10.minus result"
        self.set_pyfits(pyfits.HDUList([thisObj]))
        print("result in ", newFrame)
    @property
    def fo(self):
        return self.frames
    @fo.setter
    def fo(self, values):
        "frame only"
        allFrames = self.frames
        if isinstance(values, int):
            values = [values]
        for each in values:
            assert (each in allFrames), each
        toClose = set(allFrames).difference(set([each for each in values]))
        for eachToClose in toClose:
            self.frame = eachToClose
            self.set("frame delete")
    @property
    def math(self):
        return "do math! example: d.math='{56}-{55}', 59"
    @math.setter
    def math(self, values):
        "math between frames: example: d.math='{56}-{55}', 59"
        if type(values) in [str, unicode]:
            exp = values
            newFrame = None
        elif type(values) in [tuple, list] and len(values)==2:
            exp = values[0]
            newFrame = values[1]
        else:
            raise Exception("error input: {}".format(values))
        originExp = exp
        allFrames = self.frames
        allThisFrames = set(re_inside.findall(exp))
        d=locals()
        for each in allThisFrames:
            assert int(each) in allFrames, "no frame {}".format(each)
        for each in allThisFrames:
            self.frame = each
            toExec = "data{}=self.data".format(each)
            #print(toExec)
            exec(toExec, globals(), d)
        for each in allThisFrames:
            exp = exp.replace("{{{}}}".format(each), "data{}".format(each))

        exp = "result = " + exp
        print(exp)
        exec(exp, globals(), d)
        result = d["result"]
        if newFrame is None:
            self.set("frame new")
            newFrame = self.frame
        else:
            self.frame = newFrame
        thisObj = pyfits.PrimaryHDU()
        thisObj.data = result
        thisObj.header["EXTNAME"] = originExp
        thisObj.header["INHERIT"] = False
        thisObj.header["OBJECT"] = "ds10.math result"
        self.set_pyfits(pyfits.HDUList([thisObj]))
        print("result in ", newFrame)

    @property
    def vcon(self):
        thisFrame = self.frame
        frames = self.frames
        ny, nx = self.data.shape
        todo = []
        for eachFrame in frames:
            self.frame = eachFrame
            tny, tnx = self.data.shape
            if tnx==nx and eachFrame!=thisFrame:
                todo.append(eachFrame)
        self.frame=thisFrame
        print(thisFrame)
        print(todo)

    def vconDo(self, frames, s=2, newFrame=None, nan="red", exps=None):
        assert isinstance(frames, list)
        if exps is not None:
            assert isinstance(exps, list)
            assert len(frames)==len(exps)
        else:
            exps=[""] * len(frames)
        if len(frames)>1:
            self.frame = frames[0]
            data = self.data
            if exps[0]:
                todo = "data=" + exps[0].format("data")
                print("in frame {}: {}".format(frames[0], todo))
                exec(todo, globals(), locals())
            conData = data
            for eachFrame, eachExp in zip(frames[1:], exps[1:]):
                self.frame = eachFrame
                data = self.data
                if eachExp:
                    todo = "data=" + eachExp.format("data")
                    print("in frame {}: {}".format(eachFrame, todo))
                    exec(todo, globals(), locals())
                aux = np.ones((s, conData.shape[1])) * np.nan
                conData = np.concatenate((conData, aux, data), axis=0)
            if newFrame is None:
                self.newFrame
            else:
                self.frame = newFrame
                print("go to frame: {}".format(newFrame))

            strFrames = [str(each) for each in frames]
            thisObj = pyfits.PrimaryHDU()
            thisObj.data = conData
            thisObj.header["EXTNAME"] = ",".join(strFrames)
            thisObj.header["OBJECT"] = "ds9 combine image"
            self.set_pyfits(pyfits.HDUList([thisObj]))
            self.nan=nan
            print("\tconcatenate {}".format(", ".join(strFrames)) )

    @property
    def hcon(self):
        thisFrame = self.frame
        frames = self.frames
        ny, nx = self.data.shape
        todo = []
        for eachFrame in frames:
            self.frame = eachFrame
            tny, tnx = self.data.shape
            if tny==ny and eachFrame!=thisFrame:
                todo.append(eachFrame)
        self.frame=thisFrame
        print(thisFrame)
        print(todo)

    def hconDo(self, frames, s=2, newFrame=None, nan="red", exps=None):
        assert isinstance(frames, list)
        if exps is not None:
            assert isinstance(exps, list)
            assert len(frames)==len(exps)
        else:
            exps=[""] * len(frames)
        if len(frames)>1:
            self.frame = frames[0]
            data = self.data
            if exps[0]:
                todo = "data=" + exps[0].format("data")
                print("in frame {}: {}".format(frames[0], todo))
                exec(todo, globals(), locals())
            conData = data
            for eachFrame, eachExp in zip(frames[1:], exps[1:]):
                self.frame = eachFrame
                data = self.data
                if eachExp:
                    todo = "data=" + eachExp.format("data")
                    print("in frame {}: {}".format(eachFrame, todo))
                    exec(todo, globals(), locals())
                aux = np.ones((conData.shape[0], s)) * np.nan
                conData = np.concatenate((conData, aux, data), axis=1)
            if newFrame is None:
                self.newFrame
            else:
                self.frame = newFrame
                print("go to frame: {}".format(newFrame))

            strFrames = [str(each) for each in frames]
            thisObj = pyfits.PrimaryHDU()
            thisObj.data = conData
            thisObj.header["EXTNAME"] = ",".join(strFrames)
            thisObj.header["OBJECT"] = "ds9 combine image"
            self.set_pyfits(pyfits.HDUList([thisObj]))
            self.nan=nan
            print("\tconcatenate {}".format(", ".join(strFrames)) )

    #==> math operation between frames
    #<== help strings
    @property
    def h(self):
        print(h)
    @property
    def ha(self):
        print(h)
        print(hb)
        print(hr)
        print(hpro)
        print(hc)
        print(hf)
    @property
    def hb(self):
        print(hb)
    @property
    def hr(self):
        print(hr)
    @property
    def hpro(self):
        print(hpro)
    @property
    def hc(self):
        print(hc)
    @property
    def hf(self):
        print(hf)
    #==> help strings
