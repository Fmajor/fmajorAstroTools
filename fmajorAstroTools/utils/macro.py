# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function, unicode_literals)
# Copyright 2016 Wu Jin <wujin_astro@pku.edu.cn>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
''' Notes
    vim:set foldmethod=indent
    use vimrc config 'customFolding.vim' from https://gist.github.com/Fmajor to
        make reading the code much more easier!
'''
''' Update history
    2016-08-12 add the license
'''
''' Bugs or Cautions:
    None
'''
import warnings
import numpy as np
import astropy.io.fits as pyfits
from . import messager

def genResults(namePairs, l):
    return {eachname:l[eachvar] for eachvar, eachname in namePairs}

def cleanComments(header):
    keys = list(header)
    if len(keys)>0:
        while not header[-1] and keys[-1]=="COMMENT":
            header.pop()
            keys.pop()
            if not len(keys):
                break

class macroWarning(UserWarning):
    pass
warnings.simplefilter("always", macroWarning)

class macroError(Exception):
    pass


#!! gen list of hdu list 
genHDUListListMacro='''
try:
    unicode
except NameError:
    unicode = str
HDUListType = [astropy.io.fits.hdu.hdulist.HDUList]
output = configs.get("output", None)
rerun = configs.get("rerun", False)
{addLogger}
if type(images) not in [list, tuple]:
    raise logger.Exception("images should be a list", errorClass=macro.macroError)
hduListList = []
for eachImage in images:
    if type(eachImage) in HDUListType:
        hduListList.append(eachImage)
    elif type(eachImage) in [str, unicode]:
        hduListList.append(pyfits.open(eachImage))
    else:
        raise logger.Exception("images should be a list which contains HDULists",
                                    errorClass=macro.macroError)

headerInfoFrame = configs.get("headerInfoFrame", "first")
if headerInfoFrame not in ["first", "every"]:
    raise logger.Exception("unknown headerInfoFrame: {{}}".format(headerInfoFrame))
'''.format(addLogger=messager.initLoggerMacro)

#!! macro for beginImageProcessMacro
# gen hduList
# parameters:
#   makeZeroExt : {True, False}
#   headerInfoFrame : {"every", "first"}
#!! use   exec(macro.imageHeaderPair) # gen hdulist
genHDUListMacro= '''
try:
    unicode
except NameError:
    unicode = str
{addLogger}
HDUListType = [astropy.io.fits.hdu.hdulist.HDUList]
HDUType = [astropy.io.fits.hdu.image.PrimaryHDU, astropy.io.fits.hdu.image.ImageHDU]
output = configs.get("output", None)
rerun = configs.get("rerun", False)
imageType = type(image)
if imageType in [str, unicode]:
    imageType = "file"
elif imageType in [list, tuple]:# only accept list of numpy array
    imageType = "list"
    for each in image:
        if not isinstance(each, np.ndarray):
            raise logger.Exception("only accept np.ndarray in the list", errorClass=macro.macroError)
elif imageType in HDUListType:
    imageType = "hdulist"
elif imageType in HDUType:
    imageType = "hdu"
else:
    raise logger.Exception("error imageType: {{}}".format(imageType), errorClass=macro.macroError)

makeZeroExt = configs.get("makeZeroExt", False)

if imageType=="list":
    imageHeaderPair = [(each, pyfits.Header()) for each in image]
    primaryHDU = [(None, pyfits.Header())]
    hduList = pyfits.HDUList(primaryHDU + imageHeaderPair)
elif imageType=="file":
    obj = pyfits.open(image)
    if makeZeroExt:
        if len(obj)==1:
            primaryHDU = [pyfits.PrimaryHDU(None, obj[0].header)]
            secondHDU = [pyfits.ImageHDU(obj[0].data, pyfits.Header())]
        else:
            primaryHDU = [pyfits.PrimaryHDU(None, pyfits.Header())]
            secondHDU = [pyfits.ImageHDU(eachObj.data, eachObj.header) for eachObj in obj]
        hduList = pyfits.HDUList(primaryHDU + secondHDU)
    else:
        hduList = obj
elif imageType == "hdu":
    imageHeaderPair = [(image.data, image.header)]
    primaryHDU = [(None, pyfits.Header())]
    hduList = pyfits.HDUList(primaryHDU + imageHeaderPair)
elif imageType =="hdulist":
    hduList = image
else:
    raise Exception("should not be here... needs debug")

headerInfoFrame = configs.get("headerInfoFrame", "first")
if headerInfoFrame not in ["first", "every"]:
    raise logger.Exception("unknown headerInfoFrame: {{}}".format(headerInfoFrame))
'''.format(addLogger=messager.initLoggerMacro)


#!! before this, we should have
#!!    resultImageHeaderPair, allHeaders

addHistyroMacro="""
for eachHistory in historys:
    if type(eachHistory) in [list, tuple]:
        for eachFrame, eachHistoryStr in zip(hduList, eachHistory):
            macro.cleanComments(eachFrame.header):
            eachFrame.header.add_history(eachHistoryStr)
    elif type(eachHistory) in [str, unicode]:
        for eachFrame in hduList:
            macro.cleanComments(eachFrame.header):
            eachFrame.header.add_history(eachHistory)
"""

#!! gen outputs
genOutputMacro='''
if output is not None:
    if not os.path.exists(output) or rerun:
        hduList.writeto(output)
'''


"""notes for debug system
    to debug:
        if debug<[1|2|3|4]:

    maxDebugDepth = ?
    debugLevel
        -1      silence
         0      no debug info
         1      normal debug info (default)
         2      more debug, but the program still run without stop
         3      much more debug, but the program still run without stop
         4      the program will be block by things like matplotlib
"""

