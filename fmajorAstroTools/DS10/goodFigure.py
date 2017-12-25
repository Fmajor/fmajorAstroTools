# -*- coding: utf-8 -*-
from __future__ import (print_function, absolute_import, division)
import sys
import time
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import copy
import numpy as np
from ..utils.messager import cPrint
from scipy.optimize import curve_fit
import pdb

def gaussian(x, A, mu, sigma, C):
    return A * np.exp(- (x-mu)**2/(2 * sigma**2)) + C

try:
    import PyQt5.QtWidgets as QtWidgets
    from PyQt5.QtWidgets import QWidget
    from PyQt5.QtWidgets import QCheckBox
    from PyQt5.QtWidgets import QApplication
    from PyQt5.QtWidgets import QRadioButton
    from PyQt5.QtWidgets import QHBoxLayout
    from PyQt5.QtWidgets import QVBoxLayout
    from PyQt5.QtWidgets import QTextEdit
    from PyQt5.QtWidgets import QPushButton
    from PyQt5.QtWidgets import QLabel
    from PyQt5.QtWidgets import QGridLayout
    from PyQt5.QtWidgets import QDialog
    from PyQt5.QtWidgets import QDialogButtonBox
    from PyQt5.QtWidgets import QInputDialog
    from PyQt5.QtWidgets import QSizePolicy
    from PyQt5 import QtGui
    from PyQt5 import QtCore
    from PyQt5.QtCore import Qt
    haveQt = 1
except ImportError as e:
    print(e)
    print("have no PyQt")
    haveQt = 0

temp = plt.rcParams["keymap.back"].index("c")
if temp>-1:
    plt.rcParams["keymap.back"].pop(temp)


class GoodFigure(object):
    def __init__(self, name="good plot", parent=None, show=True,
                 ion=True, figsize=None, logfile=None, fig=None, qt=True):
        self.name = name
        self.qt = qt
        self.app = None
        if fig is None:
            self.fig = plt.figure(self.name, figsize=figsize)
            self.ax = self.fig.add_subplot(111)
        else:
            self.fig = fig
            self.ax = self.fig.get_axes()[0]
            if haveQt and self.qt:
                self.app = QApplication([])
                self._controlWindow = PlotControl(self, app=self.app)
                #self._mainWindow = PlotControlMainWindow(self, app=self.app)
                #self._controlWindow = self._mainWindow.main_widget

        #!! disable number key mapping
        manager, canvas = self.fig.canvas.manager, self.fig.canvas
        canvas.mpl_disconnect(manager.key_press_handler_id)

        self.ion = ion
        self.dynamic = PlotDynamic(self.ax, show=show, ion=ion, logfile=logfile)
        self.plots = {}
        self.parent=parent
        self.logfile = logfile
        # active data pair
        self._activeNumber = 0
        self._activeNames = []

        self.scatterConfig = {"marker":"o", "s":5, "edgecolor": "none"}
        self._colorMapLength = 20
        temp = cm.get_cmap("rainbow")
        self.colorMap = temp(np.linspace(0,1,self._colorMapLength))
        # used in fitting
        self._cidOnkey = self.fig.canvas.mpl_connect("key_press_event", self.on_key)

    @property
    def colorMapLength(self):
        return self._colorMapLength
    @colorMapLength.setter
    def colorMapLength(self, value):
        self._colorMapLength = value
        temp = cm.get_cmap("rainbow")
        self.colorMap = temp(np.linspace(0,1,self._colorMapLength))

    def _getColorMap(self, N):
        if N<0: N=0
        if N>self.colorMapLength - 1: N = self.colorMapLength-1
        return self.colorMap[N]

    def _setLimitType(self, name, configs):
        print("limit type: {}".format(name))
        if name == "auto":
            self.ax.set_ylim(auto=True)
            self.ax.relim()
            self.ax.autoscale_view(True,True,True)
            self.fig.canvas.draw()
        elif name == "updown":
            thisTop = configs["up"]
            thisBottom = configs["down"]
            if np.isfinite(thisTop):
                self.ax.set_ylim(top=thisTop)
            if np.isfinite(thisBottom):
                self.ax.set_ylim(bottom=thisBottom)
            print("top: {}, bottom: {}".format(thisTop, thisBottom))
            self.fig.canvas.draw()
        elif name == "updownRej":
            self._originData = {}
            thisUp = configs["up"]
            thisDown = configs["down"]
            mup = configs["mup"]
            mdown = configs["mdown"]
            for eachKey in self.plots:
                thisPlot = self.plots[eachKey]["obj"]
                thisx = self.plots[eachKey]["x"]
                thisy = self.plots[eachKey]["y"]
                self._originData[eachKey] = thisy
                tempy = thisy.copy()
                if np.isfinite(thisUp):
                    tempy[tempy>thisUp] = np.nan
                if np.isfinite(thisDown):
                    tempy[tempy<thisDown] = np.nan
                thisPlot.set_ydata(tempy)
                thisPlot.set_xdata(thisx)
            self.ax.set_ylim(auto=True)
            self.ax.set_xlim(auto=True)
            self.ax.relim()
            self.ax.autoscale_view(True,True,True)
            self.fig.canvas.draw()
            ydown, yup = self.ax.get_ylim()
            ydown -= mdown
            yup += mup
            self.ax.set_ylim((ydown, yup))
            for eachKey in self.plots:
                thisPlot = self.plots[eachKey]["obj"]
                thisPlot.set_ydata(self._originData[eachKey])
            print("top: {}, bottom: {}".format(yup, ydown))
            self.fig.canvas.draw()
        elif name == "percentile":
            up   = -np.inf
            down = np.inf
            up1 = configs["up1"]
            up2 = configs["up2"]
            up3 = configs["up3"]
            up4 = configs["up4"]
            setup = np.all(np.isfinite([up1, up2, up3, up4]))
            down1 = configs["down1"]
            down2 = configs["down2"]
            down3 = configs["down3"]
            down4 = configs["down4"]
            setdown = np.all(np.isfinite([down1, down2, down3, down4]))
            for eachKey in self.plots:
                thisy = self.plots[eachKey]["y"]
                if setup:
                    thisUp = (np.percentile(thisy, up1) + up2) * up3 + up4
                else:
                    thisUp = -np.inf
                if setdown:
                    thisDown = (np.percentile(thisy, down1) + down2) * down3 + down4
                else:
                    thisDown = np.inf
                if thisUp > up:
                    up = thisUp
                if thisDown < down:
                    down = thisDown
            if setup:
                self.ax.set_ylim(top=up)
            if setdown:
                self.ax.set_ylim(bottom=down)

            print("top: {}, bottom: {}".format(up, down))
            self.fig.canvas.draw()
        elif name == "xlimit":
            thisStart = configs["start"]
            thisEnd = configs["end"]
            if np.isfinite(thisStart):
                self.ax.set_xlim(left=thisStart)
            if np.isfinite(thisEnd):
                self.ax.set_xlim(right=thisEnd)
            print("start: {}, end: {}".format(thisStart, thisEnd))
            self.fig.canvas.draw()
        elif name == "xlimitTight":
            xmin = []
            xmax = []
            for eachKey in self.plots:
                thisx = self.plots[eachKey].get("x", None)
                if thisx is not None:
                    xmin.append(np.nanmin(thisx))
                    xmax.append(np.nanmax(thisx))
            if len(xmin):
                thisStart = np.min(xmin)
                self.ax.set_xlim(left=thisStart)
            else:
                thisStart = np.nan
            if len(xmax):
                thisEnd = np.max(xmax)
                self.ax.set_xlim(right=thisEnd)
            else:
                thisEnd = np.nan
            print("start: {}, end: {}".format(thisStart, thisEnd))
            self.fig.canvas.draw()

    def _setPlotStyle(self, config):
        if not self._activeNumber:
            return
        print("all actives:", self._activeNames)
        print("set style to", config)
        for eachActive in copy.copy(self._activeNames):
            eachPlot = self.plots[eachActive]
            eachType = eachPlot["type"]
            eachName = eachPlot["name"]
            cPrint("===> {} {}".format(eachActive, eachName), "y b")
            eachX    = eachPlot["x"]
            eachY    = eachPlot["y"]
            lastConfig = eachPlot["config"]
            scatterConfig = eachPlot["scatterConfig"]
            if scatterConfig is not None:
                config["scatterConfig"] = scatterConfig
            self._controlWindow.removeLine(eachName)()
            self.app.processEvents()
            print("update style: ", eachName)
            self.addPlot(eachType, eachName, eachX, eachY, **config)

    def _setScatterStyle(self, config):
        if not self._activeNumber:
            return
        print("all actives:", self._activeNames)
        print("set scatter style to", config)
        for eachActive in copy.copy(self._activeNames):
            eachPlot = self.plots[eachActive]
            eachType = eachPlot["type"]
            eachName = eachPlot["name"]
            #print("===>", eachActive, eachName)
            cPrint("===> {} {}".format(eachActive, eachName), "y b")
            eachX    = eachPlot["x"]
            eachY    = eachPlot["y"]
            lastConfig = eachPlot["config"]
            lastScatterConfig = eachPlot["scatterConfig"]
            print(lastScatterConfig)
            if lastScatterConfig is None:
                lastConfig["scatterConfig"] = config
            else:
                lastScatterConfig.update(config)
                lastConfig["scatterConfig"] = lastScatterConfig
            self._controlWindow.removeLine(eachName)()
            self.app.processEvents()
            print("update style: ", eachName)
            print(lastConfig)
            self.addPlot(eachType, eachName, eachX, eachY, **lastConfig)

    def _setPlotValue(self, text):
        if not self._activeNumber:
            return
        print("all actives:", self._activeNames)
        print("set value to", text)
        for eachActive in copy.copy(self._activeNames):
            eachPlot = self.plots[eachActive]
            eachType = eachPlot["type"]
            eachName = eachPlot["name"]
            config= eachPlot["config"]
            #print("===>", eachActive, eachName)
            cPrint("===> {} {}".format(eachActive, eachName), "y b")
            x = eachPlot["x"]
            y = eachPlot["y"]
            lastConfig = eachPlot["config"]
            scatterConfig = eachPlot["scatterConfig"]
            if eachPlot.get("originY") is None: # should be the origin version
                originY = y
                originX = x
            else:
                originY = eachPlot["originY"]
                originX = eachPlot["originX"]
            if scatterConfig is not None:
                config["scatterConfig"] = scatterConfig
            if text.lower()=="reset":
                if eachPlot.get("originY") is None: # should be the origin version
                    print("it is the origin version, noting to update")
                    continue
                y = eachPlot["originY"]
            else:
                y = eval(text)
            self._controlWindow.removeLine(eachName)()
            self.app.processEvents()
            print("update value: ", eachName)
            print(eachType, eachName, config)
            self.addPlot(eachType, eachName, x, y, **config)
            eachPlot = self.plots[eachActive]
            eachPlot["originX"] = originX
            eachPlot["originY"] = originY

    def addPlot(self, type, name, x, y, **kwargs):
        if not self.dynamic._show and self.ion:
            plt.show()
        if self.app is None:
            if haveQt and self.qt:
                self.app = QApplication([])
                self._controlWindow = PlotControl(self, app=self.app)
                print('open control window')
                #self._mainWindow = PlotControlMainWindow(self, app=self.app)
                #self._controlWindow = self._mainWindow.main_widget
                self.dynamic._show = True
        if name in self.plots:
            print("update: {}".format(name))
            thisPlot = self.plots[name]["obj"]
            thisPlot.set_xdata(x)
            thisPlot.set_ydata(y)
            color = thisPlot.get_color()
            scatterConfig = self.plots[name].get("scatterConfig")
            if scatterConfig is None:
                c = copy.copy(self.scatterConfig)
                scatterConfig = copy.deepcopy(c)
            else:
                c = copy.copy(scatterConfig)
            c.update({"color":color})
            if self.plots[name]["objScatter"] is not None:
                self.plots[name]["objScatter"].remove()
                thisScatter = self.plots[name]["objScatter"]  = self.ax.scatter(x, y, **c)
            self.ax.relim()
            self.ax.autoscale_view(True,True,True)
            self.fig.canvas.draw()
            thisDict = thisPlot
            if haveQt and self.qt:
                buttons = self.plots[name]['buttons']
                buttons["focus"].toggle()
            #thisPlot = self.plots[name]["buttons"]["focus"].toggle()
        else:
            print("add new figure: {}".format(name))
            kwargs.update({"label":name})
            scatterConfig = kwargs.get("scatterConfig", None)
            if scatterConfig is not None:
                kwargs.pop("scatterConfig")

            #!! if not set color, use rainbow color
            customColor = kwargs.get("color", None)
            if customColor is None:
                customColor = kwargs.get("c", None)
            if customColor is None:
                try:
                    index = int(name.split("@")[-1])
                    customColor = self._getColorMap(index)
                except Exception as e:
                    customColor = None

            if type=="plot":
                c=copy.copy(kwargs)
                if customColor is not None:
                    c.update({"c":customColor})
                c.update(kwargs)
                thisPlot, = self.ax.plot(x, y, **c)
                if scatterConfig is None:
                    c = copy.copy(self.scatterConfig)
                    scatterConfig = copy.deepcopy(c)
                else:
                    c = copy.copy(scatterConfig)
                color = thisPlot.get_color()
                c.update({"color":color})
                c.update(kwargs);c.pop("label")
                thisScatter = self.ax.scatter(x, y, **c)
            elif type=="step":
                c={"where":"post"}
                if customColor is not None:
                    c.update({"c":customColor})
                c.update(kwargs)
                thisPlot, = self.ax.step(x, y, **c)
                if scatterConfig is None:
                    c = copy.copy(self.scatterConfig)
                    scatterConfig = copy.deepcopy(c)
                else:
                    c = copy.copy(scatterConfig)
                color = thisPlot.get_color()
                c.update({"color":color})
                c.update(kwargs);c.pop("label")
                thisScatter = self.ax.scatter(x, y, **c)
            elif type=="simplePlot":
                c=copy.copy(kwargs)
                if customColor is not None:
                    c.update({"c":customColor})
                c.update(kwargs)
                thisPlot, = self.ax.plot(x, y, **c)
                if scatterConfig is None:
                    c = copy.copy(self.scatterConfig)
                    scatterConfig = copy.deepcopy(c)
                else:
                    c = copy.copy(scatterConfig)
                color = thisPlot.get_color()
                c.update({"color":color})
                c.update(kwargs);c.pop("label")
                thisScatter = None
            self.plots[name] = {}
            thisDict = self.plots[name]
            thisDict["type"] = type
            thisDict["obj"] = thisPlot
            thisDict["objScatter"] = thisScatter
            thisDict["name"] = name
            thisDict["x"] = x
            thisDict["y"] = y
            thisDict["config"] = kwargs
            thisDict["scatterConfig"] = scatterConfig
            if haveQt and self.qt:
                buttons = self._controlWindow.addLine(name)
                self.app.processEvents()
                thisDict["buttons"] = buttons
        # update active values: this has been done in app (the checkbox is auto selected)
        #self._activeNames.append(name)
        #self._activeNumber += 1
        self.dynamic.update(x,y, "x={{:.2f}}, {}={{:.2f}}".format(name))
        return thisDict

    #!! ??
    def addSlidePlot(self, type, name, pStart, pEnd, pDatas, reset, x, y, **kwargs):
        if name in self.plots:
            print("update: {}".format(name))
            thisObj = self.plots[name]
            thisPlot = self.plots[name]["obj"]
            if reset[0]:
                thisObj["startPoint"] = pStart
                print("update start point: {}".format(pStart))
            if reset[1]:
                thisObj["endPoint"] = pStart
                print("update end point: {}".format(pEnd))
            if np.any(reset):
                thisObj["pDatas"] = pDatas

            thisPlot.set_xdata(x)
            thisPlot.set_ydata(y)
            color = thisPlot.get_color()
            c = copy.copy(self.scatterConfig)
            c.update({"color":color})
            if self.plots[name]["objScatter"] is not None:
                self.plots[name]["objScatter"].remove()
                thisScatter = self.plots[name]["objScatter"]  = self.ax.scatter(x, y, **c)
            self.ax.relim()
            self.ax.autoscale_view(True,True,True)
            self.fig.canvas.draw()
            thisDict = thisPlot
        else:
            print("add new figure: {}".format(name))
            kwargs.update({"label":name})
            if type=="plot":
                c=copy.copy(kwargs)
                c.update(kwargs)
                thisPlot, = self.ax.plot(x, y, **c)
                c = copy.copy(self.scatterConfig)
                color = thisPlot.get_color()
                c.update({"color":color})
                c.update(kwargs);c.pop("label")
                thisScatter = self.ax.scatter(x, y, **c)
            elif type=="step":
                c={"where":"post"}
                c.update(kwargs)
                thisPlot, = self.ax.step(x, y, **c)
                c = copy.copy(self.scatterConfig)
                color = thisPlot.get_color()
                c.update({"color":color})
                c.update(kwargs);c.pop("label")
                thisScatter = self.ax.scatter(x, y, **c)
            self._controlWindow.addLine(name)
            self.app.processEvents()
            self.plots[name] = {}
            thisDict = self.plots[name]
            thisDict["type"] = type
            thisDict["obj"] = thisPlot
            thisDict["objScatter"] = thisScatter
            thisDict["name"] = name
            thisDict["x"] = x
            thisDict["y"] = y
        self.dynamic.update(x,y, "x={{:.2f}}, {}={{:.2f}}".format(name))
        return thisDict

    def _removePlot(self, name): # be called when click X
        if self.dynamic.cursor: # close cursor when delete plot
            self.dynamic.closeCursor()
        if self.dynamic.legend:
            self.dynamic.closeLegend()
        if name in self.plots: # or it may came from a invaild call
            print("remove plot: {}".format(name))
            if name in self._activeNames:
                self._activeNames.pop(self._activeNames.index(name))
                self._activeNumber -= 1
            thisPlot = self.plots[name]["obj"]
            thisScatter = self.plots[name]["objScatter"]
            self.plots.pop(name)
            thisPlot.remove()
            if thisScatter is not None:
                thisScatter.remove()
            self.ax.relim()
            self.ax.autoscale_view(True,True,True)
            self.fig.canvas.draw()
            self.app.processEvents()
        else:
            raise Exception("invaild call of this functions???")
    def _setFocus(self, name): # be called when click focus
        if name in self.plots:
            print("change focus to {}".format(name))
            self._focus = name
            self.dynamic.update(self.plots[name]["x"], self.plots[name]["y"],
                    "x={{:.2f}}, {}={{:.2f}}".format(name))
        else:
            print("change focus to None")
            self._focus = None
            self.dynamic.update(None, None, "x={{:.2f}}, None={{:.2f}}".format(name))
    def _setActive(self, name, status):
        if name not in self._activeNames and status:
            print("set {} active on".format(name))
            self._activeNames.append(name)
            self._activeNumber += 1
        elif name in self._activeNames and not status:
            print("set {} active off".format(name))
            self._activeNames.pop(self._activeNames.index(name))
            self._activeNumber -= 1
        else:
            raise Exception("information in plot and control window not match!")
    def _setVisible(self, name, status):
        assert name in self.plots, "name {} not in plots".format(name)
        self.plots[name]["obj"].set_visible(status)
        if self.plots[name]["objScatter"] is not None:
            self.plots[name]["objScatter"].set_visible(status)
        self.fig.canvas.draw()

    def gaussianFitWithBoundary(self, event):
        """fit a gaussian in figure plot
            F is band to this function
                after F is pressed, click is also band to this function
        """
        if self.dynamic.xyDataLength==-1:
            self.dynamic.openGetData(2)
            self._cidFitting = self.fig.canvas.mpl_connect('button_press_event', self.gaussianFitWithBoundary)
            print("click twice to get x boundary")
            return
        elif len(self.dynamic.xyData)<self.dynamic.xyDataLength:
            return
        elif len(self.dynamic.xyData)==self.dynamic.xyDataLength:
            focus = self.plots[self._focus]["buttons"]["focus"]
            datas = self.dynamic.xyData
            xLow = datas[0][0]
            xHigh = datas[1][0]
            allx = self.dynamic.x
            ally = self.dynamic.y
            index = np.logical_and(xLow<allx, allx<xHigh)
            x = allx[index]
            y = ally[index]
            if len(x)<5:
                print("ind:{}, x:{}, y:{}, bad fit".format(np.where(index)[0], x, y))
                return
            popt, pcov = curve_fit(gaussian, x, y, p0=(1, x.mean(), 1, 0))
            xx = np.linspace(allx.min(), allx.max(), len(allx)*10)
            yy = gaussian(xx, *popt)
            toPrintStr = "xLow: {:5}, xHigh: {:5}, A: {:10}, mu: {:.4}, sigma: {:.4}, C: {:.4}".format(xLow, xHigh, *popt)
            print(toPrintStr)
            if self.logfile is not None:
                toPrintStr = "{:15.13} {:10.8} {:15.10} {:15.10}".format(popt[1], popt[2], popt[0], popt[3])
                with open(self.logfile, "a") as f:
                    f.write(toPrintStr+"\n")
            self.addPlot("plot", "fitting" , xx, yy, color = "red")
            focus.toggle()
            self.fig.canvas.mpl_disconnect(self._cidFitting)
            self.dynamic.xyDataLength = -1
        else:
            raise Exception("why here?? need debug")

    def gaussianFit(self, event):
        """fit a gaussian in figure plot
            F is band to this function
                after F is pressed, click is also band to this function
        """
        focus = self.plots[self._focus]["buttons"]["focus"]
        allx = self.dynamic.x
        ally = self.dynamic.y
        x, y = event.xdata, event.ydata
        if x is None:
            return
        else:
            print(x,y)
        try:
            up = np.concatenate(([False],(ally[1:]-ally[:-1])>0))
            down = np.concatenate(((ally[1:]-ally[:-1])<0, [False]))
            peaks = np.where((up.astype(np.int)+down.astype(np.int))==2)[0]
            index = np.argmin(np.abs(allx[peaks] - x))
            peakIndex = peaks[index]
            xLow  = peakIndex
            xHigh = peakIndex
            while up[xLow]:
                xLow -= 1
            if xLow>0:
                xLow -= 1
            while down[xHigh]:
                xHigh += 1
            if xHigh<allx.size-1:
                xHigh +=1
            xx = np.arange(allx.size)
            index = np.logical_and(xLow<xx, xx<xHigh)
            x = allx[index]
            y = ally[index]
            if len(x)<5:
                print("ind:{}, x:{}, y:{}, bad fit".format(np.where(index)[0], x, y))
                return
            popt, pcov = curve_fit(gaussian, x, y, p0=(1, x.mean(), 1, 0))
            xx = np.linspace(allx.min(), allx.max(), len(allx)*10)
            yy = gaussian(xx, *popt)
        except Exception as e:
            print(e)
            import ipdb
            ipdb.set_trace()
        toPrintStr = "xLow: {:5}, xHigh: {:5}, A: {:10}, mu: {:.4}, sigma: {:.4}, C: {:.4}".format(xLow, xHigh, *popt)
        print(toPrintStr)
        if self.logfile is not None:
            toPrintStr = "{:15.13} {:10.8} {:15.10} {:15.10}".format(popt[1], popt[2], popt[0], popt[3])
            with open(self.logfile, "a") as f:
                f.write(toPrintStr+"\n")
        self.addPlot("plot", "fitting" , xx, yy, color = "red")
        focus.toggle()

    def toggleHideVisible(self):
        for eachName in self._activeNames:
            visible = self.plots[eachName]["buttons"]["visible"]
            visible.toggle()

    def getPoints(self, event):
        """fit a gaussian in figure plot
            F is band to this function
                after F is pressed, click is also band to this function
        """
        if self.dynamic.xyDataLength==-1:
            self.dynamic.openGetData(2)
            self._cidFitting = self.fig.canvas.mpl_connect('button_press_event', self.getPoints)
            self._datatype = event.key
            print("click twice to get x data pair")
            return
        elif len(self.dynamic.xyData)<self.dynamic.xyDataLength:
            return
        elif len(self.dynamic.xyData)==self.dynamic.xyDataLength:
            datas = self.dynamic.xyData
            toPrintStr = "{:20.15f} {:20.15f} {:4}".format(datas[0][0], datas[1][0], self._datatype)
            print(toPrintStr)
            if self.logfile is not None:
                with open(self.logfile, "a") as f:
                    f.write(toPrintStr+"\n")
            self.fig.canvas.mpl_disconnect(self._cidFitting)
            self._datatype = None
            self.dynamic.xyDataLength = -1
        else:
            raise Exception("why here?? need debug")

    def on_key(self, event):
        if event.key=="F":
            self.gaussianFitWithBoundary(None)
        if event.key=="f":
            self.gaussianFit(event)
        if event.key=="z":
            self.toggleHideVisible()
        if event.key in ["1", "2", "3", "4", "5", "6", "7", "8", "9"]:
            self.getPoints(event)

#!! use an axis to init a PlotDynamic
class PlotDynamic(object):
    def __init__(self, ax, show=True,
                 cursor_x_format=".2f", cursor_y_format=".2f",
                 ion=True, logfile=None):
        self.ax = ax
        self.fig = ax.get_figure()
        self.fig.subplots_adjust(left=0.08, right=0.95, top=0.95, bottom=0.05)
        self._cidOnkey = self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        self._show=show
        self.logfile = logfile
        if ion:
            plt.ion()
        if show:
            plt.show()
        # active data pair
        #self.x = np.array([]) # for the mostly used single cursor
        #self.y = np.array([]) # for the mostly used single cursor
        self.x = None
        self.y = None
        self.focusButton = None
        self._cursor_x_format=cursor_x_format
        self._cursor_y_format=cursor_y_format
        self._cursorStr = "x={{:{x}}}, y={{:{y}}}".\
                format(x=self._cursor_x_format, y=self._cursor_y_format)
        # init cursor
        self._cursor = False
        self.ly = None
        # init legend
        self._legend = -1
        self._legendObj = None
        self.xyData = []
        self.xyDataLength = -1
    def _mouse_move_cursor(self, event):
        if not event.inaxes: return
        xdata = event.xdata
        indx = np.searchsorted(self.x, [xdata])[0]
        if indx<self.x.size:
            x = self.x[indx]
            y = self.y[indx]
            self.ly.set_xdata(x)
            self.lx.set_ydata(y)
            self.ax.set_title( self._cursorStr.format(x,y)+ " index: {}".format(indx) )
        else:
            self.ax.set_title( self._cursorStr.format(self.x[-1],self.y[-1]) )
        self.fig.canvas.draw()
    def _mouse_move_simple_cursor(self, event):
        if not event.inaxes: return
        x, y = event.xdata, event.ydata
        self.lx.set_ydata(y )
        self.ly.set_xdata(x )
        self.ax.set_title( 'x=%1.2f, y=%1.2f'%(x,y) )
        self.fig.canvas.draw()

    #<=== get data and fit
    def _mouse_get_data(self, event):
        if not event.inaxes: return
        x, y = event.xdata, event.ydata
        self.xyData.append((x,y))
        print("x:{} y:{}".format(x, y))
        if len(self.xyData) == self.xyDataLength:
            self.closeGetData()

    def openGetData(self, N=1):
        self.xyData = []
        self.xyDataLength = N
        self._cidGetData = self.fig.canvas.mpl_connect(\
                'button_press_event', self._mouse_get_data)
        if not self._cursor:
            self._cursor = 1 - self._cursor
            success = self.openCursor()
            if not success:
                self._cursor = 1 - self._cursor

    def closeGetData(self):
        self.fig.canvas.mpl_disconnect(self._cidGetData)
        if self._cursor:
            self._cursor = 1 - self._cursor
            self.closeCursor()
        return self.xyData
    #===> todos

    def recordPosition(self, event):
        x, y = event.xdata, event.ydata
        if x is None:
            return
        toPrintStr = "{:15.13} {:10.8}".format(x, y)
        print(toPrintStr)
        if self.logfile is not None:
            with open(self.logfile, "a") as f:
                f.write(toPrintStr+"\n")
    def recordError(self, event):
        x, y = event.xdata, event.ydata
        if x is None:
            return
        toPrintStr = "!!XXXXX!! {:15.13} {:10.8}".format(x, y)
        print(toPrintStr)
        if self.logfile is not None:
            with open(self.logfile, "a") as f:
                f.write(toPrintStr+"\n")
    def on_key(self, event):
        if event.key=="c":  # normal cursor
            if self.ly is None:
                x, y = event.xdata, event.ydata
                if x is None:
                    x = np.mean(self.ax.get_xlim())
                    y = np.mean(self.ax.get_ylim())
                self.ly = self.ax.axvline(x, color='k')  # the vert line
                self.ly.set_alpha(0)
                self.lx = self.ax.axhline(y, color='k')  # the vert line
                self.lx.set_alpha(0)
            self._cursor = 1 - self._cursor
            if self._cursor:
                success = self.openCursor()
                if not success:
                    self._cursor = 1 - self._cursor
            else:
                self.closeCursor()
        elif event.key=="e":
            if self._legend==-1:
                self._legend=1
                self._legendObj = self.ax.legend(loc=0)
                self._legendObj.set_visible(True)
                self.fig.canvas.draw()
            else:
                self._legend = 1 - self._legend
                if self._legend == 1:
                    self._legendObj = self.ax.legend(loc=0)
                if self._legendObj is not None:
                    self._legendObj.set_visible(self._legend)
                self.fig.canvas.draw()
        elif event.key=="d": # simple cursor
            if self.ly is None:
                x, y = event.xdata, event.ydata
                if x is None:
                    x = np.mean(self.ax.get_xlim())
                    y = np.mean(self.ax.get_ylim())
                self.ly = self.ax.axvline(x, color='k')  # the vert line
                self.ly.set_alpha(0)
                self.lx = self.ax.axhline(y, color='k')  # the vert line
                self.lx.set_alpha(0)
            self._cursor = 1 - self._cursor
            if self._cursor:
                self._cidcursor = self.fig.canvas.mpl_connect('motion_notify_event', self._mouse_move_simple_cursor)
                self.lx.set_alpha(1)
                self.ly.set_alpha(1)
                self.fig.canvas.draw()
            else:
                self.fig.canvas.mpl_disconnect(self._cidcursor)
                self.lx.set_alpha(0)
                self.ly.set_alpha(0)
                self.ax.set_title("")
                self.fig.canvas.draw()
        elif event.key=="D":
            if self.ly is None:
                x, y = event.xdata, event.ydata
                if x is None:
                    x = np.mean(self.ax.get_xlim())
                    y = np.mean(self.ax.get_ylim())
                self.ly = self.ax.axvline(x, color='k')  # the vert line
                self.ly.set_alpha(0)
                self.lx = self.ax.axhline(y, color='k')  # the vert line
                self.lx.set_alpha(0)
            self._cursor = 1 - self._cursor
            if self._cursor:
                self._cidcursor = self.fig.canvas.mpl_connect('motion_notify_event', self._mouse_move_cursor)
                self.lx.set_alpha(1)
                self.ly.set_alpha(1)
                self.fig.canvas.draw()
            else:
                self.fig.canvas.mpl_disconnect(self._cidcursor)
        if event.key=="w":
            self.recordPosition(event)
        if event.key=="x":
            self.recordError(event)

    def openCursor(self):
        if self.x is not None:
            self._cursor = 1
            self._cidcursor = self.fig.canvas.mpl_connect(\
                    'motion_notify_event', self._mouse_move_cursor)
            self.ly.set_alpha(1)
            self.lx.set_alpha(1)
            self.fig.canvas.draw()
            return True
        else:
            return False
    def closeCursor(self):
        self._cursor = 0
        self.fig.canvas.mpl_disconnect(self._cidcursor)
        self.ly.set_alpha(0)
        self.lx.set_alpha(0)
        self.ax.set_title("")
        self.fig.canvas.draw()
    def closeLegend(self):
        self._legend = 0
        if self._legendObj is not None:
            self._legendObj.set_visible(self._legend)

    def update(self, x, y, title=None):
        if title is None:
            self._cursorStr = "x={{:{x}}}, y={{:{y}}}".\
                    format(x=self._cursor_x_format, y=self._cursor_y_format)
        else:
            self._cursorStr = title
        self.x = x
        self.y = y
        self.fig.canvas.draw()
    @property
    def cursor(self):
        return self._cursor
    @property
    def legend(self):
        return self._legend

#<== some functions
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
    if not input:
        return {}
    allCommands = input.split(",")
    result={}
    for each in allCommands:
        aux = each.split("=")
        assert len(aux)==2, "plotDict format error!"
        if aux[1]:
            result[aux[0]]=getValue(aux[1])
    return result
#==> some functions

if haveQt:
    class PlotControlMainWindow(QtWidgets.QMainWindow):
         def __init__(self, parent = None, app=None):
             #QtWidgets.QMainWindow.__init__(self)
             super(PlotControlMainWindow, self).__init__()
             self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
             self.main_widget = PlotControl(parent=parent, app=app)

    #!! plot control gui
    class PlotControl(QWidget):
        def __init__(self, parent = None, app=None):
            super(PlotControl, self).__init__()
            self.setWindowTitle('control')
            self.setGeometry(20, 30, 180, 300)
            self.labelFormat = "{:>10}"
            self._focusName = "None"
            self.bodys={}
            self.app = app
            self.limitConfigs={}
            self._focusGroup = QtWidgets.QButtonGroup()
            self._limitGroup = QtWidgets.QButtonGroup()
            self.parent = parent

            self.genLayout()
            self.show()

        #<== gen layout
        def genLayout(self):
            self.rootVbox = QVBoxLayout()
            self.rootVbox.setContentsMargins(5,5,5,5)
            self.rootVbox.setAlignment(Qt.AlignTop)
            self.rootVbox.setSpacing(0)

            self.headerVbox = QVBoxLayout()
            self.setLayout(self.rootVbox)

            self.headerVbox.setContentsMargins(0,0,0,0)
            #self.headerVbox.setMargin(0)
            self.headerVbox.setSpacing(0)
            self.headerWidget = QWidget()
            self.headerWidget.setLayout(self.headerVbox)
            self.headerWidget.setContentsMargins(0,0,0,0)
            self.rootVbox.addWidget(self.headerWidget)

            #!! gen layout
            self.limitControl()
            self.header()

            self.lineVbox = QVBoxLayout()
            self.lineVbox.setContentsMargins(0,0,0,0)
            #self.lineVbox.setMargin(0)
            self.lineVbox.setSpacing(0)
            self.lineVbox.setAlignment(Qt.AlignTop)
            self.lineWidget = QWidget()
            self.lineWidget.setLayout(self.lineVbox)
            self.lineWidget.setContentsMargins(0,0,0,0)

            self.scroll = QtWidgets.QScrollArea()
            #self.scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
            self.scroll.setWidgetResizable(True)
            self.scroll.setWidget(self.lineWidget)
            self.scroll.setContentsMargins(0,0,0,0)
            self.rootVbox.addWidget(self.scroll)

            if self.parent is None:
                for i in range(20):
                    self.addLine(str(i))

            self.footVbox = QVBoxLayout()
            self.footVbox.setAlignment(Qt.AlignTop)
            self.footVbox.setSpacing(0)
            self.footVbox.setContentsMargins(0,0,0,0)
            self.footWidget = QWidget()
            self.footWidget.setLayout(self.footVbox)
            self.rootVbox.addWidget(self.footWidget)
            self.footer()

        def limitControl(self):
            vbox = QVBoxLayout()
            vbox.setSpacing(0)
            box = QHBoxLayout()
            box1 = QHBoxLayout()
            vbox.addLayout(box); vbox.addLayout(box1);
            self.headerVbox.addLayout(vbox)
            # label and text field
            lbl = QtWidgets.QLabel("start", self)
            self._x_limit_start = QtWidgets.QLineEdit(self)
            aux = self.limitConfigs.get("_x_limit_start","")
            if aux:
                self._x_limit_start.setText(aux)
            else:
                self._x_limit_start.setText("-")
            box.addWidget(lbl)
            box.addWidget(self._x_limit_start)
            box.addStretch(1)
            # label and text field
            lbl = QtWidgets.QLabel("end", self)
            self._x_limit_end = QtWidgets.QLineEdit(self)
            aux = self.limitConfigs.get("_x_limit_end","")
            if aux:
                self._x_limit_end.setText(aux)
            else:
                self._x_limit_end.setText("-")
            box.addWidget(lbl)
            box.addWidget(self._x_limit_end)
            box.addStretch(1)
            # U button
            self._x_limit_button_update = QtWidgets.QPushButton("U", self)
            self._x_limit_button_update.clicked.connect(self._x_limit_button_f)
            box1.addWidget(self._x_limit_button_update)
            # T button
            self._x_limit_button_tight = QtWidgets.QPushButton("T", self)
            self._x_limit_button_tight.clicked.connect(self._x_limit_button_ft)
            box1.addWidget(self._x_limit_button_tight)

            # ylimit control
            vbox = QVBoxLayout()
            vbox.setSpacing(5)
            hbox = QHBoxLayout()
            vbox.addLayout(hbox)
            self.headerVbox.addLayout(vbox)
            hbox_auto = QVBoxLayout()
            hbox_updown = QVBoxLayout()
            hbox_updownRej = QVBoxLayout()
            hbox_percentile = QVBoxLayout()
            vvbox1 = QVBoxLayout()
            vvbox2 = QVBoxLayout()
            hbox.addLayout(vvbox1)
            hbox.addLayout(vvbox2)
            vvbox1.addLayout(hbox_auto)
            vvbox2.addLayout(hbox_updown)
            vvbox1.addLayout(hbox_updownRej)
            vvbox2.addLayout(hbox_percentile)
            hhbox = QHBoxLayout()
            hhbox.setContentsMargins(0,0,0,0)
            #hhbox.setMargin(0)
            self._limitLaylout = hhbox
            vbox.addLayout(hhbox)
            rb = QRadioButton('Auto', self)
            rb.toggle()
            self._limitType = "auto"
            self._limitGroup.addButton(rb)
            rb.toggled.connect(self.limits("auto"))
            hbox_auto.addWidget(rb)

            rb = QRadioButton('updown', self)
            rb.toggled.connect(self.limits("updown"))
            self._limitGroup.addButton(rb)
            hbox_updown.addWidget(rb)

            rb = QRadioButton('updownRej', self)
            rb.toggled.connect(self.limits("updownRej"))
            self._limitGroup.addButton(rb)
            hbox_updownRej.addWidget(rb)

            rb = QRadioButton('percentile', self)
            rb.toggled.connect(self.limits("percentile"))
            self._limitGroup.addButton(rb)
            hbox_percentile.addWidget(rb)
        def header(self):
            hbox = QHBoxLayout()
            hbox.setContentsMargins(0,5,0,0)
            hbox.setSpacing(10)
            self.headerVbox.addLayout(hbox)
            hbox1 = QHBoxLayout()
            hbox2 = QHBoxLayout()
            hbox1.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
            hbox2.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
            hbox.addLayout(hbox1)
            hbox.addLayout(hbox2)
            # visible
            t = QLabel(' V', self)
            hbox1.addWidget(t)
            # active
            t = QLabel(' A', self)
            hbox1.addWidget(t)
            # cursor
            t = QLabel(' C', self)
            hbox1.addWidget(t)
            # name
            t = QLabel(self.labelFormat.format('name'), self)
            hbox2.addWidget(t)
            # close
            t = QLabel('      X      ', self)
            hbox2.addWidget(t)
            ## second line
            hbox = QHBoxLayout()
            hbox.setSpacing(10)
            hbox.setContentsMargins(0,0,0,0)
            #hbox.setMargin(0)
            hbox.allObjs = {}
            self.headerVbox.addLayout(hbox)
            hbox1 = QHBoxLayout()
            hbox2 = QHBoxLayout()
            hbox1.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
            hbox2.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
            hbox.addLayout(hbox1)
            hbox.addLayout(hbox2)
            # visible all
            cb = QCheckBox('', self)
            cb.stateChanged.connect(self.changeAllVisible)
            hbox1.addWidget(cb)
            hbox.allObjs["visible"] = cb
            # active all
            cb = QCheckBox('', self)
            cb.stateChanged.connect(self.changeAllActive)
            hbox1.addWidget(cb)
            hbox.allObjs["active"] = cb
            # cursor none
            self._focus_None = QRadioButton('', self)
            self._focus_None.toggle()
            self._focus_None.toggled.connect(self.focus("None"))
            hbox1.addWidget(self._focus_None)
            self._focusGroup.addButton(self._focus_None)
            hbox.allObjs["focus"] = self._focus_None
            # name
            t = QLabel(self.labelFormat.format('all'), self)
            hbox2.addWidget(t)
            hbox.allObjs["name"] = t
            # close
            pb = QPushButton('X', self)
            pb.clicked.connect(self.removeAllLine)
            hbox2.addWidget(pb)
            hbox.allObjs["close"] = t
        #!!  self.lineVbox.addLayout(hbox), self.bodys[name] = hbox
        def addLine(self, name):
            hbox = QHBoxLayout()
            hbox.setContentsMargins(0,0,0,0)
            #hbox.setMargin(0)
            hbox.setSpacing(10)
            self.lineVbox.addLayout(hbox)
            self.bodys[name] = hbox
            buttons = {}

            hbox.plotname = name
            hbox.allObjs = {}

            hbox1 = QHBoxLayout(); hbox2 = QHBoxLayout()
            hbox1.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
            hbox2.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
            hbox.addLayout(hbox1); hbox.addLayout(hbox2)
            # visible
            cb = QCheckBox('', self)
            cb.stateChanged.connect(self.visible(name))
            cb.toggle()
            hbox1.addWidget(cb)
            hbox.allObjs["visible"] = cb
            buttons["visible"] = cb
            # active
            cb = QCheckBox('', self)
            cb.stateChanged.connect(self.active(name))
            cb.toggle()
            hbox1.addWidget(cb)
            hbox.allObjs["active"] = cb
            buttons["active"] = cb
            # for cursor
            rb = QRadioButton('', self)
            rb.toggled.connect(self.focus(name))
            self._focusGroup.addButton(rb)
            rb.toggle()
            hbox1.addWidget(rb)
            hbox.allObjs["focus"] = rb
            buttons["focus"] = rb
            # name
            nl = QLabel(self.labelFormat.format(name), self)
            nl.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
            hbox2.addWidget(nl)
            hbox.allObjs["name"] = nl
            # delete
            pb = QPushButton('X', self)
            pb.clicked.connect(self.removeLine(name))
            hbox2.addWidget(pb)
            hbox.allObjs["delete"] = pb
            return buttons
        def footer(self):
            # get data
            # self.text = QTextEdit()
            # self.rootVbox.addWidget(self.text)
            # label and text field
            #!! plot value update
            hbox = QHBoxLayout()
            hbox.setContentsMargins(0,0,0,0)
            #hbox.setMargin(0)
            hbox.setSpacing(10)
            self.footVbox.addLayout(hbox)
            lbl = QtWidgets.QLabel("y", self)
            self._plotValue = QtWidgets.QLineEdit(self)
            self._plotValue.setText("")
            self._plotValue.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                          QtWidgets.QSizePolicy.Fixed)
            hbox.addWidget(lbl)
            hbox.addWidget(self._plotValue)
            # U button
            self._plotValue_update = QtWidgets.QPushButton("U", self)
            self._plotValue_update.clicked.connect(self._plotValue_button_f)
            hbox.addWidget(self._plotValue_update)
            #!! plot style update
            hbox = QHBoxLayout()
            hbox.setContentsMargins(0,0,0,0)
            #hbox.setMargin(0)
            hbox.setSpacing(10)
            self.footVbox.addLayout(hbox)
            lbl = QtWidgets.QLabel("style", self)
            self._plotStyle = QtWidgets.QLineEdit(self)
            self._plotStyle.setText("")
            self._plotStyle.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                          QtWidgets.QSizePolicy.Fixed)
            hbox.addWidget(lbl)
            hbox.addWidget(self._plotStyle)
            # U button
            self._plotStyle_update = QtWidgets.QPushButton("U", self)
            self._plotStyle_update.clicked.connect(self._plotStyle_button_f)
            hbox.addWidget(self._plotStyle_update)
            #!! scatter style update
            hbox = QHBoxLayout()
            hbox.setContentsMargins(0,0,0,0)
            #hbox.setMargin(0)
            hbox.setSpacing(10)
            self.footVbox.addLayout(hbox)
            lbl = QtWidgets.QLabel("scatter", self)
            self._plotScatter = QtWidgets.QLineEdit(self)
            self._plotScatter.setText("")
            self._plotScatter.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                          QtWidgets.QSizePolicy.Fixed)
            hbox.addWidget(lbl)
            hbox.addWidget(self._plotScatter)
            # U button
            self._plotScatter_update = QtWidgets.QPushButton("U", self)
            self._plotScatter_update.clicked.connect(self._plotScatter_button_f)
            hbox.addWidget(self._plotScatter_update)

        #==> gen layout

        #<== todo
        def addSlideLine(self, name, maxValue=100):#!! just a demo
            Vbox = QVBoxLayout()
            hbox = QHBoxLayout()
            hhbox = QHBoxLayout()
            hhhbox = QHBoxLayout()
            Vbox.addLayout(hbox)
            Vbox.addLayout(hhbox)
            Vbox.addLayout(hhhbox)
            Vbox.plotname = name
            self.bodys[name] = Vbox
            self.lineVbox.addLayout(Vbox)
            Vbox.allObjs = {}
            hbox1 = QHBoxLayout()
            hbox2 = QHBoxLayout()
            hbox1.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
            hbox2.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
            hbox.addLayout(hbox1)
            hbox.addLayout(hbox2)
            # visible
            cb = QCheckBox('', self)
            cb.toggle()
            cb.stateChanged.connect(self.visible(name))
            hbox1.addWidget(cb)
            Vbox.allObjs["visible"] = cb
            # active
            cb = QCheckBox('', self)
            cb.stateChanged.connect(self.active(name))
            cb.toggle()
            hbox1.addWidget(cb)
            Vbox.allObjs["active"] = cb
            # for cursor
            rb = QRadioButton('', self)
            rb.toggled.connect(self.focus(name))
            self._focusGroup.addButton(rb)
            rb.toggle()
            hbox1.addWidget(rb)
            Vbox.allObjs["focus"] = rb
            # name
            nl = QLabel(self.labelFormat.format(name), self)
            nl.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
            hbox2.addWidget(nl)
            Vbox.allObjs["name"] = nl
            # delete
            pb = QPushButton('X', self)
            pb.clicked.connect(self.removeLine(name))
            hbox2.addWidget(pb)
            Vbox.allObjs["delete"] = pb

            # use inner data
            cb = QCheckBox('', self)
            cb.toggle()
            #cb.stateChanged.connect(self.visible(name))
            hhbox.addWidget(cb)
            Vbox.allObjs["innerData"] = cb
            # slide
            sld = QtWidgets.QSlider(QtCore.Qt.Horizontal, self)
            sld.valueChanged.connect(self.slide(name))
            sld.setMinimum(0)
            sld.setMaximum(maxValue)
            sld.setTickInterval(1)
            sld.setSingleStep(1)
            hhbox.addWidget(sld)
            Vbox.allObjs["slide"] = sld
            # reset
            pb = QPushButton('R', self)
            #pb.clicked.connect(self.removeLine(name))
            hhbox.addWidget(pb)
            Vbox.allObjs["reset"] = pb
            # update
            pb = QPushButton('U', self)
            #pb.clicked.connect(self.removeLine(name))
            hhbox.addWidget(pb)
            Vbox.allObjs["update"] = pb

            lbx1 = QtWidgets.QLabel("x", self)
            x1 = QtWidgets.QLineEdit(self)
            lby1 = QtWidgets.QLabel("y", self)
            y1 = QtWidgets.QLineEdit(self)
            pb1 = QPushButton('U', self)
            lbx2 = QtWidgets.QLabel("x", self)
            x2 = QtWidgets.QLineEdit(self)
            lby2 = QtWidgets.QLabel("y", self)
            y2 = QtWidgets.QLineEdit(self)
            pb2 = QPushButton('U', self)
            Vbox.allObjs["updatePoint1"] = pb
            Vbox.allObjs["updatePoint2"] = pb
            Vbox.allObjs["x1"] = x1
            Vbox.allObjs["x2"] = x2
            Vbox.allObjs["y1"] = y1
            Vbox.allObjs["y2"] = y2
            hhhbox.addWidget(lbx1)
            hhhbox.addWidget(x1)
            hhhbox.addWidget(lby1)
            hhhbox.addWidget(y1)
            hhhbox.addWidget(pb1)
            hhhbox.addWidget(lbx2)
            hhhbox.addWidget(x2)
            hhhbox.addWidget(lby2)
            hhhbox.addWidget(y2)
            hhhbox.addWidget(pb2)
        def slide(self, name):
            def setter(value):
                print("name: {}, value:{}".format(name, value))
            return setter

        #==> todo

        #<== buttun function
        def limits(self, name):
            def setter(state):
                if state:
                    if self._limitType!= name:
                        self._limitType=name
                        self.clearLayout(self._limitLaylout)
                        box = self._limitLaylout
                        if self._limitType == "auto":
                            if self.parent is not None:
                                self.parent._setLimitType("auto", None)
                        elif self._limitType == "updown":
                            lbl = QtWidgets.QLabel("down", self)
                            self._limit_down = QtWidgets.QLineEdit(self)
                            aux = self.limitConfigs.get("_limit_down","")
                            if aux:
                                self._limit_down.setText(aux)
                            else:
                                self._limit_down.setText("-")
                            box.addWidget(lbl)
                            box.addWidget(self._limit_down)
                            box.addStretch(1)
                            lbl = QtWidgets.QLabel("up", self)
                            self._limit_up = QtWidgets.QLineEdit(self)
                            aux = self.limitConfigs.get("_limit_up","")
                            if aux:
                                self._limit_up.setText(aux)
                            else:
                                self._limit_up.setText("-")
                            box.addWidget(lbl)
                            box.addWidget(self._limit_up)
                            box.addStretch(1)
                            self._limit_button_update = QtWidgets.QPushButton("U", self)
                            self._limit_button_update.clicked.connect(self._limit_button_f)
                            box.addWidget(self._limit_button_update)
                        elif self._limitType == "updownRej":
                            vbox = QVBoxLayout()
                            box.addLayout(vbox)
                            hbox1  = QHBoxLayout()
                            hbox2  = QHBoxLayout()
                            hbox3  = QHBoxLayout()
                            vbox.addLayout(hbox1)
                            vbox.addLayout(hbox2)
                            vbox.addLayout(hbox3)
                            lbl = QtWidgets.QLabel("down", self)
                            self._limit_reject_down = QtWidgets.QLineEdit(self)
                            aux = self.limitConfigs.get("_limit_reject_down","")
                            if aux:
                                self._limit_reject_down.setText(aux)
                            else:
                                self._limit_reject_down.setText("-")
                            hbox1.addWidget(lbl)
                            hbox1.addWidget(self._limit_reject_down)
                            hbox1.addStretch(1)
                            lbl = QtWidgets.QLabel("up", self)
                            self._limit_reject_up = QtWidgets.QLineEdit(self)
                            aux = self.limitConfigs.get("_limit_reject_up","")
                            if aux:
                                self._limit_reject_up.setText(aux)
                            else:
                                self._limit_reject_up.setText("-")
                            hbox1.addWidget(lbl)
                            hbox1.addWidget(self._limit_reject_up)
                            hbox1.addStretch(1)
                            self._reject_limit_button_update = QtWidgets.QPushButton("U", self)
                            self._reject_limit_button_update.clicked.connect(self._reject_limit_button_f)
                            #qle.textChanged[str].connect(self.onChanged)
                            hbox3.addWidget(self._reject_limit_button_update)

                            lbl1 = QtWidgets.QLabel("margins", self)
                            lbl2 = QtWidgets.QLabel("d", self)
                            self._limit_reject_down_margin = QtWidgets.QLineEdit(self)
                            aux = self.limitConfigs.get("_limit_reject_down_margin","")
                            if aux:
                                self._limit_reject_down_margin.setText(aux)
                            else:
                                self._limit_reject_down_margin.setText("0")
                            hbox2.addWidget(lbl1)
                            hbox2.addWidget(lbl2)
                            hbox2.addWidget(self._limit_reject_down_margin)
                            hbox2.addStretch(1)
                            lbl = QtWidgets.QLabel("u", self)
                            self._limit_reject_up_margin = QtWidgets.QLineEdit(self)
                            aux = self.limitConfigs.get("_limit_reject_up_margin","")
                            if aux:
                                self._limit_reject_up_margin.setText(aux)
                            else:
                                self._limit_reject_up_margin.setText("0")
                            hbox2.addWidget(lbl)
                            hbox2.addWidget(self._limit_reject_up_margin)
                        elif self._limitType == "percentile":
                            vbox = QVBoxLayout()
                            vbox.setSpacing(0)
                            box.addLayout(vbox)
                            hbox1  = QHBoxLayout()
                            hbox2  = QHBoxLayout()
                            hbox3  = QHBoxLayout()
                            vbox.addLayout(hbox1)
                            vbox.addLayout(hbox2)
                            vbox.addLayout(hbox3)
                            hbox1.addWidget(QtWidgets.QLabel("( P(",self))
                            self._up_1 = QtWidgets.QLineEdit(self)
                            aux = self.limitConfigs.get("_up_1","")
                            if aux:
                                self._up_1.setText(aux)
                            else:
                                self._up_1.setText("100")
                                self.limitConfigs["_up_1"] = "100"
                            hbox1.addWidget(self._up_1)
                            hbox1.addWidget(QtWidgets.QLabel(")+",self))
                            self._up_2 = QtWidgets.QLineEdit(self)
                            aux = self.limitConfigs.get("_up_2","")
                            if aux:
                                self._up_2.setText(aux)
                            else:
                                self._up_2.setText("0")
                                self.limitConfigs["_up_2"] = "0"
                            hbox1.addWidget(self._up_2)
                            hbox1.addWidget(QtWidgets.QLabel(")x",self))
                            self._up_3 = QtWidgets.QLineEdit(self)
                            aux = self.limitConfigs.get("_up_3","")
                            if aux:
                                self._up_3.setText(aux)
                            else:
                                self._up_3.setText("1")
                                self.limitConfigs["_up_3"] = "0"
                            hbox1.addWidget(self._up_3)
                            hbox1.addWidget(QtWidgets.QLabel("+",self))
                            self._up_4 = QtWidgets.QLineEdit(self)
                            aux = self.limitConfigs.get("_up_4","")
                            if aux:
                                self._up_4.setText(aux)
                            else:
                                self._up_4.setText("0")
                                self.limitConfigs["_up_4"] = "0"
                            hbox1.addWidget(self._up_4)

                            hbox2.addWidget(QtWidgets.QLabel("( P(",self))
                            self._down_1 = QtWidgets.QLineEdit(self)
                            aux = self.limitConfigs.get("_down_1","")
                            if aux:
                                self._down_1.setText(aux)
                            else:
                                self._down_1.setText("0")
                                self.limitConfigs["_down_1"] = "0"
                            hbox2.addWidget(self._down_1)
                            hbox2.addWidget(QtWidgets.QLabel(")+",self))
                            self._down_2 = QtWidgets.QLineEdit(self)
                            aux = self.limitConfigs.get("_down_2","")
                            if aux:
                                self._down_2.setText(aux)
                            else:
                                self._down_2.setText("0")
                                self.limitConfigs["_down_2"] = "0"
                            hbox2.addWidget(self._down_2)
                            hbox2.addWidget(QtWidgets.QLabel(")x",self))
                            self._down_3 = QtWidgets.QLineEdit(self)
                            aux = self.limitConfigs.get("_down_3","")
                            if aux:
                                self._down_3.setText(aux)
                            else:
                                self._down_3.setText("1")
                                self.limitConfigs["_down_3"] = "0"
                            hbox2.addWidget(self._down_3)
                            hbox2.addWidget(QtWidgets.QLabel("+",self))
                            self._down_4 = QtWidgets.QLineEdit(self)
                            aux = self.limitConfigs.get("_down_4","")
                            if aux:
                                self._down_4.setText(aux)
                            else:
                                self._down_4.setText("0")
                                self.limitConfigs["_down_4"] = "0"
                            hbox2.addWidget(self._down_4)
                            self._percentile_button_update = QtWidgets.QPushButton("U", self)
                            self._percentile_button_update.clicked.connect(self._percentile_button_f)
                            hbox3.addWidget(self._percentile_button_update)
                    self.app.processEvents()
            return setter
        def visible(self, name):
            setterName = "setVisible"
            def setter(state):
                if self.parent is not None:
                    try:
                        if state == QtCore.Qt.Checked:
                            print("{} visible on".format(name))
                            self.parent._setVisible(name, True)
                        else:
                            print("{} visible off".format(name))
                            self.parent._setVisible(name, False)
                    except:
                        raise Exception("failed to {}".format(setterName))
                else:
                    if state == QtCore.Qt.Checked:
                        print("{}: {} {}".format(setterName, name, True))
                    else:
                        print("{}: {} {}".format(setterName, name, False))
            return setter
        def active(self, name):
            setterName = "setActive"
            def setter(state):
                if self.parent is not None:
                    if state == QtCore.Qt.Checked:
                       self.parent._setActive(name, True)
                    else:
                       self.parent._setActive(name, False)
                else:
                    if state == QtCore.Qt.Checked:
                        print("{}: {} {}".format(setterName, name, True))
                    else:
                        print("{}: {} {}".format(setterName, name, False))
            return setter
        def focus(self, name):
            setterName = "setFocus"
            def setter(state):
                if self.parent is not None:
                    try:
                        if state:
                            self._focusName = name
                            self.parent._setFocus(name)
                    except:
                        raise Exception("failed to {}".format(setterName))
                else:
                    if state:
                        self._focusName = name
                        print("{}: {} {}".format(setterName, name, True))
                    else:
                        print("{}: {} {}".format(setterName, name, False))
            return setter
        def changeAllVisible(self, state):
            boxType = "visible"
            thisKeys = list(self.bodys.keys())
            if state == QtCore.Qt.Checked:
                for eachName in thisKeys:
                    widget = self.getCheckBoxObj(eachName, boxType)
                    if widget is not None:
                        widget.setChecked(True)
            else:
                for eachName in thisKeys:
                    widget = self.getCheckBoxObj(eachName, boxType)
                    if widget is not None:
                        widget.setChecked(False)
        def changeAllActive(self, state):
            boxType = "active"
            thisKeys = list(self.bodys.keys())
            if state == QtCore.Qt.Checked:
                for eachName in thisKeys:
                    widget = self.getCheckBoxObj(eachName, boxType)
                    if widget is not None:
                        widget.setChecked(True)
            else:
                for eachName in thisKeys:
                    widget = self.getCheckBoxObj(eachName, boxType)
                    if widget is not None:
                        widget.setChecked(False)
        def removeAllLine(self):
            ok = SimpleDialog.ok("Do you want to delete all the plots?")
            #text, ok = QtWidgets.QInputDialog.getText(self, 'Input Dialog', 'Enter your name:')
            if ok:
                for eachKey in list(self.bodys.keys()):
                    print("remove:", eachKey)
                    self.bodys[eachKey].allObjs['delete'].click()
        def _x_limit_button_f(self):
            start = self._x_limit_start.text()
            end = self._x_limit_end.text()
            try:
                start = float(start)
                self.limitConfigs["_x_limit_start"] = str(start)
            except:
                self.limitConfigs["_x_limit_start"] = "-"
                start = np.nan
            try:
                end = float(end)
                self.limitConfigs["_x_limit_end"] = str(end)
            except:
                self.limitConfigs["_x_limit_end"] = "-"
                end = np.nan
            configs = {"start":start, "end":end}
            print(configs)
            if self.parent is not None:
                self.parent._setLimitType("xlimit", configs)
        def _x_limit_button_ft(self):
            if self.parent is not None:
                self.parent._setLimitType("xlimitTight", None)
        def _limit_button_f(self):
            up =   self._limit_up.text()
            down = self._limit_down.text()
            try:
                up = float(up)
                self.limitConfigs["_limit_up"] = str(up)
            except:
                self.limitConfigs["_limit_up"] = "-"
                up = np.nan
            try:
                down = float(down)
                self.limitConfigs["_limit_down"] = str(down)
            except:
                self.limitConfigs["_limit_down"] = "-"
                down = np.nan
            configs = {"up":up, "down":down}
            print(configs)
            if self.parent is not None:
                self.parent._setLimitType("updown", configs)
        def _reject_limit_button_f(self):
            up =   self._limit_reject_up.text()
            down = self._limit_reject_down.text()
            mup =   self._limit_reject_up_margin.text()
            mdown = self._limit_reject_down_margin.text()
            try:
                up = float(up)
                self.limitConfigs["_limit_reject_up"] = str(up)
            except:
                self.limitConfigs["_limit_reject_up"] = "-"
                up = np.nan
            try:
                down = float(down)
                self.limitConfigs["_limit_reject_down"] = str(down)
            except:
                self.limitConfigs["_limit_reject_down"] = "-"
                down = np.nan
            try:
                mdown = float(mdown)
                self.limitConfigs["_limit_reject_down_margin"] = str(mdown)
            except:
                self.limitConfigs["_limit_reject_down_margin"] = "-"
                mdown = 0
            try:
                mup = float(mup)
                self.limitConfigs["_limit_reject_up_margin"] = str(mup)
            except:
                self.limitConfigs["_limit_reject_up_margin"] = "-"
                mup = 0
            configs = {"up":up, "down":down, "mup":mup, "mdown":mdown}
            print(configs)
            if self.parent is not None:
                self.parent._setLimitType("updownRej", configs)
        def _percentile_button_f(self):
            up1 = self._up_1.text()
            up2 = self._up_2.text()
            up3 = self._up_3.text()
            up4 = self._up_4.text()
            down1 = self._down_1.text()
            down2 = self._down_2.text()
            down3 = self._down_3.text()
            down4 = self._down_4.text()
            try:
                up1 = float(up1)
                assert 0<=up1<=100
                self.limitConfigs["_up_1"] = str(up1)
            except:
                up1 = np.nan
                self.limitConfigs["_up_1"] = "-"
            try:
                up2 = float(up2)
                self.limitConfigs["_up_2"] = str(up2)
            except:
                self.limitConfigs["_up_2"] = "-"
                up2 = np.nan
            try:
                up3 = float(up3)
                self.limitConfigs["_up_3"] = str(up3)
            except:
                self.limitConfigs["_up_3"] = "-"
                up3 = np.nan
            try:
                up4 = float(up4)
                self.limitConfigs["_up_4"] = str(up4)
            except:
                self.limitConfigs["_up_4"] = "-"
                up4 = np.nan
            try:
                down1 = float(down1)
                assert 0<=down1<=100
                self.limitConfigs["_down_1"] = str(down1)
            except:
                self.limitConfigs["_down_1"] = "-"
                down1 = np.nan
            try:
                down2 = float(down2)
                self.limitConfigs["_down_2"] = str(down2)
            except:
                down2 = np.nan
                self.limitConfigs["_down_2"] = '-'
            try:
                down3 = float(down3)
                self.limitConfigs["_down_3"] = str(down3)
            except:
                down3 = np.nan
                self.limitConfigs["_down_3"] = "-"
            try:
                down4 = float(down4)
                self.limitConfigs["_down_4"] = str(down4)
            except:
                down4 = np.nan
                self.limitConfigs["_down_4"] = "-"
            configs = {"up1":up1, "up2":up2, "up3":up3, "up4":up4,
                       "down1":down1, "down2":down2, "down3":down3, "down4":down4}
            print(configs)
            if self.parent is not None:
                self.parent._setLimitType("percentile", configs)
        def _plotStyle_button_f(self):
            text = str(self._plotStyle.text())
            config = getDict(text)
            if self.parent is not None:
                self.parent._setPlotStyle(config)
        def _plotScatter_button_f(self):
            text = str(self._plotScatter.text())
            config = getDict(text)
            if self.parent is not None:
                self.parent._setScatterStyle(config)

        def _plotValue_button_f(self):
            text = str(self._plotValue.text())
            if self.parent is not None and text:
                self.parent._setPlotValue(text)

        #==> buttun function

        #<== layout tools
        def clearLayout(self, layout):
            for i in reversed(range(layout.count())):
                thisWidget = layout.itemAt(i).widget()
                thisLayout = layout.itemAt(i).layout()
                thisSpace  = layout.itemAt(i).spacerItem()
                if thisWidget is not None:
                    thisWidget.setParent(None)
                elif thisLayout is not None:
                    self.clearLayout(thisLayout)
                elif thisSpace is not None:
                    pass
                else:
                    raise Exception("unknown object: type({}) = {}".format(\
                            layout.itemAt(i), type(layout.itemAt(i))))
        def removeLine(self, name):
            def remove():
                self.clearLayout(self.bodys[name])
                self.lineVbox.removeItem(self.bodys[name])
                self.bodys.pop(name)
                if self._focusName==name:
                    self._focus_None.toggle()
                    self.focusName = "None"
                if self.parent is not None:
                    self.parent._removePlot(name)
            return remove
        def getCheckBoxObj(self, name, boxName):
            "ex.getCheckBoxStatus('1111','visible')"
            if boxName not in ["visible", "active", "focus"]:
                raise("unknown box name: {}".format(boxName))
            if name not in list(self.bodys.keys()):
                raise("unknown name: {}".format(name))
            if boxName=="visible":
                result =  self.bodys[name].itemAt(0).itemAt(0).widget()
                return result
            if boxName=="active":
                result =  self.bodys[name].itemAt(0).itemAt(1).widget()
                return result
            if boxName=="focus":
                result =  self.bodys[name].itemAt(0).itemAt(2).widget()
                return result

        #==> layout tools

    class SimpleDialog(QDialog):
        def __init__(self, text="", parent = None):
            super(SimpleDialog, self).__init__(parent)

            layout = QVBoxLayout(self)
            t = QLabel(text, self)
            layout.addWidget(t)

            # OK and Cancel buttons
            self.buttons = QDialogButtonBox(
                QDialogButtonBox.Ok | QDialogButtonBox.Cancel,
                Qt.Horizontal, self)
            layout.addWidget(self.buttons)
            self.buttons.accepted.connect(self.accept)
            self.buttons.rejected.connect(self.reject)

        # static method to create the dialog and return (date, time, accepted)
        @staticmethod
        def ok(text="", parent = None):
            dialog = SimpleDialog(text, parent)
            result = dialog.exec_()
            return QDialog.Accepted == result

if __name__ == '__main__':
    if haveQt:
        app = QApplication(sys.argv)
        ex = PlotControl(app=app)
