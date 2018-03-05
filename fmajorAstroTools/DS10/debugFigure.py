import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.rcParams['keymap.save'] = ''

class DebugFigure(object):
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
        self._raw_title = self.ax.get_title()
        self._title = self._raw_title

        self._doSelect = False
        self._selectLine = None
    def _mouse_move_cursor(self, event):
        if not event.inaxes: return
        xdata = event.xdata
        indx = np.searchsorted(self.x, [xdata])[0]
        if indx<self.x.size:
            x = self.x[indx]
            y = self.y[indx]
            self.ly.set_xdata(x)
            self.lx.set_ydata(y)
            self.ax.set_title(self.title + '\n' + self._cursorStr.format(x,y)+ " index: {}".format(indx) )
        else:
            self.ax.set_title( self.title + '\n' +self._cursorStr.format(self.x[-1],self.y[-1]) )
        self.fig.canvas.draw()
    def _mouse_move_simple_cursor(self, event):
        if not event.inaxes: return
        x, y = event.xdata, event.ydata
        self.lx.set_ydata(y)
        self.ly.set_xdata(x)
        self.ax.set_title('{}\nx={:1.2f}, y={:1.2f}'.format(self.title, x, y))
        self.fig.canvas.draw()

    def _mouse_move_select(self, event):
        if not self._doSelect: return
        if not event.inaxes: return
        x, y = event.xdata, event.ydata
        self.selectLine(x, y)

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
                self.ax.set_title(self.title)
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
        elif event.key=="s":
            self._doSelect = not self._doSelect
            if self._doSelect:
                self._cidselect = self.fig.canvas.mpl_connect('motion_notify_event', self._mouse_move_select)
            else:
                self.fig.canvas.mpl_disconnect(self._cidselect)
        elif event.key=="S":
            if self._selectLine:
                self._selectLine.set_lw(self._lastLW)
            self.fig.canvas.draw()
        elif event.key=="w":
            self.recordPosition(event)
        elif event.key=="x":
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
        self.ax.set_title(self.title)
        self.fig.canvas.draw()
    def closeLegend(self):
        self._legend = 0
        if self._legendObj is not None:
            self._legendObj.set_visible(self._legend)

    def selectLine(self, x, y):
        minDis = np.inf
        if self._selectLine:
            self._selectLine.set_lw(self._lastLW)
        for eachline in self.ax.lines:
            thisxy = eachline.get_xydata()
            mindis = ((thisxy - np.array([[x, y]]))**2).sum(axis=1).min()
            if mindis < minDis:
                minDis = mindis
                self._selectLine = eachline
                self.x = eachline.get_xdata()
                self.y = eachline.get_ydata()
        self._lastLW = self._selectLine.get_lw()
        self._selectLine.set_lw(3)
        self.title = '{}\nlabel: {}'.format(self._raw_title, self._selectLine.get_label())
        self.ax.set_title(self.title)
        self.fig.canvas.draw()

if __name__ == "__main__":
    fig = plt.figure()
    ax = fig.add_subplot(111)

    for i in range(10):
        ax.plot(np.arange(20), np.random.rand(20), label="{}".format(i))

    gf = DebugFigure(ax)
