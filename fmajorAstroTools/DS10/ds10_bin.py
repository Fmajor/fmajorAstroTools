# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function, unicode_literals)
import os
import sys
import getopt
import copy
import re
import numpy as np
import ipdb
import pickle
#import pyfits
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from fmajorAstroTools.DS10 import ds9
from fmajorAstroTools.utils import imheader
p = plt

usageStr="""ds10 filenames_group0 [options] [filenames_group1 [options]] ...
    -n: norm the plot
    -i: initial plot(not use wcs)
    -e: do expression before plot. Example: -e"{}/1000"
    -a: use all frames in fits files
    --name: create a new ds9 window or connect to a named ds9 window
    --alter: show fits alternatively.
        Example: ds10 1.fits 2.fits -f1,2
            show 1.fits[1] 2.fits[1] 1.fits[2] 2.fits[2]
                instead of 1.fits[1,2] 2.fits[1,2]
    -f: set frames to use. Example: -f1,2,3 -f1~3,5~7

    -b: filter files to open using header information
    --boolFrames: specify which frame to be use for bool filter (default is all frames)
    -B: filter extensions to open using header information

    --width: width for the ds9 window
    --height: width for the ds9 window
    --regions: open region file
    --selfRegions: open region file with the same name
    --exec: make exec after plot
    --exit: not stay in python interactive mode
    --notOpen: not open 2D image and 1D array, just print filter result

    --: no options for args before

    -p: set plot style. Example: -p'color=red,ls=dashed', for 1D plot
    --xlim: '(start,end)|tight' for 1D plot
    --ylim: 'updown(down,up)' for 1D plot
"""

# some functions are in development
usageStrFull="""ds10 filenames_group0 [options] [filenames_group1 [options]]
    -i: initial plot(not use wcs)
    -e: do expression before plot. Example: -e"{}/1000"
    -a: use all frames in fits files
    --name: create a new ds9 window or connect to a named ds9 window
    --alter: show fits alternatively.
        Example: ds10 1.fits 2.fits -f1,2
            show 1.fits[1] 2.fits[1] 1.fits[2] 2.fits[2]
                instead of 1.fits[1,2] 2.fits[1,2]
    -f: set frames to use. Example: -f1,2,3 -f1~3,5~7

    -b: filter files to open using header information
    -B: filter extensions to open using header information

    --width: width for the ds9 window
    --height: width for the ds9 window
    --regions: open region file
    --selfRegions: open region file with the same name
    --exec: make exec after plot
    --exit: not stay in python interactive mode
    --notOpen: not open 2D image and 1D array, just print filter result
    --: no options for args before

    -n: norm the plot, for 1D array
    -p: set plot style. Example: -p'color=red,ls=dashed', for 1D plot
    --xlim: '(start,end)|tight' for 1D plot
    --ylim: 'updown(down,up)' for 1D plot
    --force1d: force to plot the 2d array
    --vc: concatenate vertically
    --hc: concatenate horizontally
"""


# TODO: falter

def usage():
    print(usageStr)

def wcs2xy(data, header, ydata, rotateWCS="False"):
    assert len(data.shape)==1
    wcs = WCS(header)
    print(wcs)
    pixel_x = np.arange(data.size)
    pixel_y = np.zeros((data.size,))
    #ipdb.set_trace()
    if wcs.naxis==2:
        axisToUse = 1
        mask = (wcs.wcs.crval==1.0)|(np.abs(np.diag(wcs.wcs.cd))==1.0)
        print(mask)
        if mask.sum()==1:
            if mask[1]:
                axisToUse = 1-axisToUse

        if rotateWCS: axisToUse = 1-axisToUse
        if axisToUse: # for vertical spectrum
            pixel = np.array([pixel_y, pixel_x], dtype=np.int).T
            xdata = wcs.all_pix2world(pixel, 1)[:,1].flatten()
        else:
            pixel = np.array([pixel_x, pixel_y], dtype=np.int).T
            xdata = wcs.all_pix2world(pixel, 1)[:,0].flatten()
    else:
        xdata = wcs.all_pix2world(pixel_x.reshape((data.size,1)),1).flatten()
    inds = np.argsort(xdata)
    xdata = xdata[inds]
    ydata = ydata[inds]
    print("xdata range:", xdata[0], xdata[-1])
    return xdata, ydata

def orderedArgvParser(args, options="", long_options=[]):#<==
    """demo:
    todo = "-abc --test3 arg1 arg2 -ab -c -e123 -f'123' -g 123 --test= --test1='asd' -h '123' arg3 -bc -e 'lala' 'arg4' -d --test3 "
    result = orderedArgvParser(todo.split(), "abcde:f:g:h:", ["test=","test1=","test2","test3"])
    print(todo)
    for each in result:
        print(each)
    print("\\n\\n")
    opts, args = getopt.gnu_getopt(todo.split(),
            "abcde:f:g:h:", ["test=","test1=","test2","test3"])
    print(opts)
    print(args)
    """#==>
    #<== parse options
    l_op = set()
    s_op = set()
    l_lop = set()
    s_lop = set()
    options+="!"
    aux = options.split(":")
    for each in aux: # for short options
        if len(each)==1:
            l_op.add(each)
        else:
            l_op.add(each[-1])
            s_op.update({eachStr for eachStr in each[:-1]})
    l_op.remove("!")

    for each in long_options:
        if each[-1]=="=":
            l_lop.add(each[:-1])
        else:
            s_lop.add(each)
    #print(l_op, s_op)
    #print(l_lop, s_lop)
    #==> parse options
    #<== parse args
    result = []
    thisOpts = []
    thisArgs = []
    isLastArg = True
    nextMustBeArg=False
    ofArg = ''
    for op in args:
        #print("procesing", op)
        if op[:2]=="--":
            isLastArg = False
            if ofArg:
                thisOpts.append(("-{}".format(ofArg), ""))
                ofArg=""
            op = op[2:]
            if not op:
                isLastArg = False
                continue
            if "=" in op: #!! long option
                ind = op.index("=")
                name = op[:ind]
                value = op[ind+1:]
                if name not in l_lop:
                    raise Exception("error long_options: --{}".format(op))
                thisOpts.append(("--{}".format(name), value))
            else:
                if op not in s_lop:
                    raise Exception("error long_options: --{}".format(op))
                else:
                    thisOpts.append(("--{}".format(op), ""))
        elif op[:1]=="-": #!! short option
            isLastArg = False
            if ofArg:
                thisOpts.append(("-{}".format(ofArg), ""))
                ofArg=""
            op = op[1:]
            if op[0] in l_op:
                if len(op)>1:#!! short option with arg like this: -k123
                    name = op[0]
                    value = op[1:]
                    thisOpts.append(("-{}".format(name), value))
                else:#!! may be like: -f 123
                    ofArg = op
            elif op[0] in s_op:
                if len(op)>1:#!! short option array like -abckdef
                    for each in op:
                        if each in s_op:
                            thisOpts.append(("-{}".format(each), ""))
                        else:
                            raise Exception("unknown short config in -{}".format(op))
                else:
                    thisOpts.append(("-{}".format(op), ""))
            else:
                raise Exception("unknown short config in -{}".format(op))

        else:#!! process args
            if ofArg:
                value = op
                thisOpts.append(("-{}".format(ofArg), value))
                ofArg=""
                isLastArg = False
            else:
                if not isLastArg: #!! last is not Arg but options, so push all this* in to all*
                    result.append((thisOpts, thisArgs))
                    thisOpts = []
                    thisArgs = []
                thisArgs.append(op)
                isLastArg = True
    if len(thisArgs) or len(thisOpts):
        result.append((thisOpts, thisArgs))

    return result
    #==> parse args

def getValue(string):
    """get value from string, if it's like
        if it's like 'blabla' or "blabla", it's string
            or we will try to return int and float
    """
    if len(string)>2:
        if (string[0]=="'" and string[-1]=="'") or (string[0]=='"' and string[-1]=='"'):
            return string[1:-1]
    try:
        int(string)
    except:
        try:
            result = float(string)
            return result
        except:
            return string
    return int(string)
def getDict(input):
    allCommands = input.split(",")
    result={}
    for each in allCommands:
        aux = each.split("=")
        assert len(aux)==2, "plotDict format error!"
        if aux[1]:
            result[aux[0]]=getValue(aux[1])
    return result

def main(argv):
    #<== setups
    argv = argv[1:]
    reInBracket = re.compile(r"\((.*),(.*)\)")
    xlim = ["tight", {}]
    ylim = None
    #==> setups

    print(argv)
    fileArgOpts = orderedArgvParser(argv, "af:re:niB:b:p:",
            ["boolFrames=", "xlim=", "ylim=", "force1d", "exit", "notOpen",
             "vc", "hc", "name=", "exec=", "savefig=", "notShow",
             "alter", "width=", "height=", "regions=", "selfRegions"])

    ds9Name="ds10"
    width = None
    height = None
    file_ext_confs = []
    fileExt_names = []
    execStr = []
    savefig = ""
    notShow = False
    for eachOpts, eachFiles in fileArgOpts:
        thisAlter = False
        imheaderConfig = {
            "frames":     "0",
            "allFrames":  False,
            "boolFrames": "",
            "boolExp":    "",
            "extBoolExp":    "",
            "mode":       "feq"
        }
        eachConfig = {
            "rotateWCS":False,
            "normPlot":False,
            "useWcs":True,
            "expStr":"",
            "plotConfig":{},
            "force1d":False,
            "concatenate":"",
            "regions":"",
            "selfRegions":False,
        }
        notOpen = False
        for op, value in eachOpts:
            if op=="-r":
                eachConfig["rotateWCS"] = True
            elif op == "-b":
                imheaderConfig["boolExp"] = value
            elif op == "-B":
                imheaderConfig["extBoolExp"] = value
            elif op == "--boolFrames":
                imheaderConfig["boolFrames"] = value
            elif op == "-f":
                imheaderConfig["frames"] = value
            elif op == "-a":
                imheaderConfig["allFrames"] = True
            elif op=="-e":
                if value[0]=="'" and value[-1]=="'": value = value[1:-1]
                if value[0]=='"' and value[-1]=='"': value = value[1:-1]
                eachConfig["expStr"] = value
            elif op=="-n":
                eachConfig["normPlot"] = True
            elif op=="-i":
                eachConfig["useWcs"] = False
            elif op=="-p":
                eachConfig["plotConfig"] = getDict(value)
            elif op=="--xlim":
                if value=="tight":
                    xlim = ["xlimitTight", None]
                continue
                result = reInBracket.findall(value)
                if len(result)!=1:
                    raise Exception("unknown config: {}:{}".format(op, value))
                else:
                    xlim = ["xlimit", {"start":float(result[0][0]), "end":float(result[0][1])}]
            elif op=="--ylim":
                result = reInBracket.findall(value)
                if len(result)!=1:
                    raise Exception("unknown config: {}:{}".format(op, value))
                else:
                    ylim = ["updown", {"down":float(result[0][0]), "up":float(result[0][1])}]
            elif op=="--force1d":
                eachConfig["force1d"] = True
            elif op=="--savefig":
                savefig = value
            elif op=="--notShow":
                notShow = True
            elif op=="--vc":
                eachConfig["concatenate"] = "v"
            elif op=="--hc":
                eachConfig["concatenate"] = "h"
            elif op=="--name":
                ds9Name = value
            elif op=="--exec":
                execStr = value
            elif op=="--regions":
                eachConfig["regions"] = value
            elif op=="--selfRegions":
                eachConfig["selfRegions"] = True
            elif op=="--alter":
                thisAlter = True
            elif op=="--width":
                width = int(value)
            elif op=="--height":
                height = int(value)
            elif op=="--notOpen":
                notOpen = True

        fileInfo = imheader.main(eachFiles, imheaderConfig)
        todoFiles = fileInfo["todoFiles"]
        todoFrames = fileInfo["todoFrames"]
        if thisAlter:
            fileN = len(todoFiles)
            assert fileN==len(todoFrames)
            allFramesLen = [len(each) for each in todoFrames]
            assert np.all(np.diff(allFramesLen)==0), "todoFrame for each Fits file do not have the same length"
            allFramesLen = allFramesLen[0]
            aux_todoFiles = todoFiles * allFramesLen
            aux_todoFrames = ["?"] * fileN * allFramesLen
            index=0
            for frameIndex in range(allFramesLen):
                for fileIndex in range(fileN):
                    aux_todoFrames[index] = [todoFrames[fileIndex][frameIndex]]
                    index+=1
            todoFiles = aux_todoFiles
            todoFrames = aux_todoFrames

        file_ext_confs.append((todoFiles, todoFrames, eachConfig))
        fileExt_names.append(fileInfo["names"])

    if len(argv)==0:
        usage()
        d = ds9.ds10(ds9Name, figShow=False, width=width, height=height)
        return d

    d = None
    plot1d = False
    if savefig: notShow = True
    if notShow:
        plt.switch_backend("Agg")
    for i, (eachGroupFile, eachGroupFrames, eachConfig) in enumerate(file_ext_confs):
        eachExp = eachConfig.get("expStr")
        eachPlotConfig = eachConfig.get("plotConfig")
        eachConcatenate = eachConfig.get("concatenate")
        print("============================")
        print("Group ", i)
        fileExt_name_count = -1
        for each in fileExt_names[i]:
            print(each)
        if eachExp:
            print("\tdo ", eachExp)
        print("----------------------------")
        if notOpen:
            continue

        conNum = 0
        conFrame=[]
        for j, (eachFile, eachFrames) in enumerate(zip(eachGroupFile, eachGroupFrames)):
            obj = pyfits.open(eachFile)
            if eachConfig.get("normPlot", 0): #!! get max value for norm
                eachDatas = [obj[eachFrame].data for eachFrame in eachFrames if obj[eachFrame].data is not None]
                eachMax = [np.nanmax(each) for each in eachDatas]
                if len(eachMax)>0:
                    maxInFits = np.nanmax(eachMax)
                else:
                    maxInFits = 1
            else:
                maxInFits = 1
            force1d = eachConfig.get("force1d",False)
            regionFile = [eachConfig.get("regions", "")]
            if eachConfig["selfRegions"]:
                eachFilePath, eachFileName = os.path.split(eachFile)
                regionFile.append(os.path.join(eachFilePath, eachFileName.split(".fits")[0]+".reg"))
            #print(eachConfig)
            for eachFrame in eachFrames:
                fileExt_name_count += 1
                thisObj = obj[eachFrame]
                if not thisObj.header.get("XTENSION")=='BINTABLE' and thisObj.header.get("NAXIS",0)>1 and not force1d:
                    if d is None:
                        d = ds9.ds10(ds9Name, figShow=False, width=width, height=height)
                        d.t=[]
                        d.tn=[]
                    if eachExp:
                        data = thisObj.data
                        header = thisObj.header
                        newFrame = d.newFrame
                        d.set_pyfits(pyfits.HDUList([obj[0], thisObj]))
                        todo = "data="+eachExp.format("data")
                        print("exec({})".format(todo))
                        l=locals()
                        exec(todo, globals(), l)
                        thisObj.data = l['data']
                        thisObj.header = l['header']
                        thisName = "{}[{}]".format(os.path.basename(eachFile), eachFrame)
                        message = "{} after {}".\
                                format(thisName, eachExp)
                        thisObj.header["EXTNAME"] = message
                        d.set_pyfits(pyfits.HDUList([obj[0], thisObj]))
                        print("\t", thisName, "", ds9.getDimStr(thisObj.header))
                    else:
                        newFrame = d.fits(eachFile, eachFrame)
                    print(regionFile)
                    for eachReg in regionFile:
                        if os.path.exists(eachReg):
                            print("\tload region file {}".format(eachReg))
                            d.set("regions {}".format(eachReg))
                        elif eachReg:
                            print("\tregion file {} not exists".format(eachReg))
                    #!! process concatenate
                    if eachConcatenate and conNum==0:
                        conData = thisObj.data
                        conNum += 1
                        conFrame.append(str(newFrame))
                    elif eachConcatenate=="v":
                        aux = np.ones((2, conData.shape[1])) * np.nan
                        conData = np.concatenate((conData, thisObj.data), axis=0)
                        conNum += 1
                        conFrame.append(str(newFrame))
                    elif eachConcatenate=="h":
                        aux = np.ones((conData.shape[0],2)) * np.nan
                        conData = np.concatenate((conData, aux, thisObj.data), axis=1)
                        conNum += 1
                        conFrame.append(str(newFrame))
                    else:
                        pass
                elif not thisObj.header.get("XTENSION")=='BINTABLE' and force1d:
                    if d is None:
                        d = ds9.ds10(ds9Name, figShow=False, width=width, height=height)
                        d.t=[]
                        d.tn=[]
                    data = thisObj.data
                    shape = data.shape
                    print(shape)
                    if shape[0]>shape[1]:
                        data = data.T
                    for k,eachLine in enumerate(data):
                        x = np.arange(eachLine.size)
                        y = eachLine/maxInFits
                        if eachExp:
                            todo = "y="+eachExp.format("y")
                            print("exec({})".format(todo))
                            exec(todo, globals())
                        d.goodfigure.addPlot("simplePlot", "{}[{}][{}]".format(eachFile, eachFrame, k), x, y, **eachPlotConfig)
                else:# plot 1D
                    if d is None:
                        if notShow:
                            d = ds9.ds10(ds9Name, figShow=False, width=width, height=height, ds9Show=False, qt=False)
                        else:
                            d = ds9.ds10(ds9Name, figShow=False, width=width, height=height)
                        d.t=[]
                        d.tn=[]
                    if thisObj.header.get("NAXIS",0)==1:
                        plot1d = True
                        if eachConfig.get("useWcs",0):
                            x,y = wcs2xy(thisObj.data, thisObj.header,
                                    thisObj.data/maxInFits,
                                    rotateWCS=eachConfig.get("rotateWCS"))
                        else:
                            x = np.arange(thisObj.data.size)
                            y = thisObj.data/maxInFits
                        if eachExp:
                            todo = "y="+eachExp.format("y")
                            print("exec({})".format(todo))
                            l=locals()
                            exec(todo, globals(), l)
                            y=l['y']
                        d.goodfigure.addPlot("simplePlot", "{}[{}]".format(eachFile, eachFrame), x, y, **eachPlotConfig)

                    if thisObj.header.get("XTENSION")=='BINTABLE':
                        d.t.append(thisObj)
                        d.tn.append(fileExt_names[i][fileExt_name_count])
        if conNum>0:
            newFrame = d.newFrame
            thisObj = pyfits.PrimaryHDU()
            thisObj.data = conData
            thisObj.header["EXTNAME"] = ",".join(conFrame)
            thisObj.header["OBJECT"] = "ds10 combine image"
            if not notOpen:
                d.set_pyfits(pyfits.HDUList([thisObj]))
            d.nan="red"
            print("\tconcatenate {}".format(", ".join(conFrame)) )

    if execStr:
        print("exec({})".format(execStr))
        exec(execStr, globals(), locals())

    if plot1d:# set xlim and ylim
        if xlim is not None:
            d.goodfigure._setLimitType(xlim[0], xlim[1])
        if ylim is not None:
            print(ylim)
            d.goodfigure._setLimitType(ylim[0], ylim[1])
        if savefig:
            with open(savefig, "wb") as f:
                pickle.dump(d.goodfigure.fig, f)

    return d

if __name__=="__main__":
    d=main(sys.argv)
