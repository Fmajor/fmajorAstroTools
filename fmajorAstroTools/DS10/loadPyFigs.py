#!/usr/local/bin/ipython3
# -*- coding: utf-8 -*-

import sys
import os
import pickle
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from fmajorAstroTools.DS10.debugFigure import DebugFigure
import ipdb
import getopt

def loadFile(filename):
    with open(filename, 'rb') as f:
        obj = pickle.load(f)
        if type(obj)==list:
            return obj
        elif type(obj)==matplotlib.figure.Figure:
            return [obj]
        elif isinstance(obj, dict):
            thisClass = obj.pop("class")
            return [thisClass(**obj)]
        else:
            raise Exception("unknown data type of {}", type(obj))

opts, args = getopt.gnu_getopt(sys.argv[1:], "",
        ["debug", "gui", "class"])
filenames = args

debug = False
gui = False
for op, value in opts:
    if op == "--debug":
        debug = True
    if op == "--gui":
        gui = True

if len(filenames)<1:
    print("usage: loadPyFig filename1 filename2 ...")
    sys.exit(1)

for each in filenames:
    if not os.path.exists(each):
        print("invalid file name:{}".format(each))

allObj = []
for each in filenames:
    allObj += loadFile(each)

fig = plt.gcf()
ax = plt.gca()
if gui:
    debugFigure = DebugFigure(ax)
elif debug:
    plt.show(0)
else:
    plt.show()
