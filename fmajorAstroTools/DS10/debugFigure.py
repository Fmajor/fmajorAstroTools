import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.rcParams['keymap.save'] = ''

tt = None
class DebugFigure(object):
    def __init__(self, fig, show=True,
                 cursor_x_format=".2f", cursor_y_format=".2f",
                 ion=True, logfile=None):
        self.fig = fig
        #self.fig.subplots_adjust(left=0.08, right=0.95, top=0.95, bottom=0.05)
        self._cidOnkey = self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        self._show=show
        self.logfile = logfile
        if ion:
            plt.ion()
        if show:
            plt.show()
        self.data = {}
        self._cursor = False
        self._doSelect = False
        for eachax in self.fig.axes:
            thisd = self.data[eachax] = {}

            thisd['ax'] = eachax
            thisd['x'] = None
            thisd['y'] = None
            thisd['focusButton'] = None
            thisd['_cursor_x_format'] = cursor_x_format
            thisd['_cursor_y_format'] = cursor_y_format
            thisd['_cursorStr'] = "x={{:{x}}}, y={{:{y}}}".\
                    format(x=thisd['_cursor_x_format'], y=thisd['_cursor_y_format'])
            # init cursor
            thisd['ly'] = None
            # init legend
            thisd['_legend'] = -1
            thisd['_legendObj'] = None
            thisd['xyData'] = []
            thisd['xyDataLength'] = -1
            thisd['_raw_title'] = eachax.get_title()
            thisd['title'] = thisd['_raw_title']

            thisd['_selectLine'] = None
    def _mouse_move_cursor(self, event):
        if not event.inaxes: return
        d = self.data[event.inaxes]
        if d['x'] is None:
            return
        if d['ly'] is None:
            x, y = event.xdata, event.ydata
            if x is None:
                x = np.mean(d['ax'].get_xlim())
                y = np.mean(d['ax'].get_ylim())
            d['ly'] = d['ax'].axvline(x, color='k')  # the vert line
            d['lx'] = d['ax'].axhline(y, color='k')  # the vert line
        xdata = event.xdata
        d['lx'].set_alpha(1)
        d['ly'].set_alpha(1)
        sx = d['x']
        sy = d['y']
        indx = np.searchsorted(sx, [xdata])[0]
        if indx<sx.size:
            x = sx[indx]
            y = sy[indx]
            d['ly'].set_xdata(x)
            d['lx'].set_ydata(y)
            d['ax'].set_title(d['title'] + '\n' + d['_cursorStr'].format(x,y)+ " index: {}".format(indx) )
        else:
            self.ax.set_title( d['title'] + '\n' + d['_cursorStr'].format(sx[-1], sy[-1]) )
        self.fig.canvas.draw()
    def _mouse_move_simple_cursor(self, event):
        if not event.inaxes: return
        d = self.data[event.inaxes]
        if d['ly'] is None:
            x, y = event.xdata, event.ydata
            if x is None:
                x = np.mean(d['ax'].get_xlim())
                y = np.mean(d['ax'].get_ylim())
            d['ly'] = d['ax'].axvline(x, color='k')  # the vert line
            d['lx'] = d['ax'].axhline(y, color='k')  # the vert line
        x, y = event.xdata, event.ydata
        d['lx'].set_alpha(1)
        d['ly'].set_alpha(1)
        d['lx'].set_ydata(y)
        d['ly'].set_xdata(x)
        d['ax'].set_title('{}\nx={:1.2f}, y={:1.2f}'.format(d['title'], x, y))
        self.fig.canvas.draw()
    def _mouse_move_select(self, event):
        if not self._doSelect: return
        if not event.inaxes: return
        self.selectLine(event)

    def recordPosition(self, event):
        if not event.inaxes: return
        x, y = event.xdata, event.ydata
        toPrintStr = "{:15.13} {:10.8}".format(x, y)
        print(toPrintStr)
        if self.logfile is not None:
            with open(self.logfile, "a") as f:
                f.write(toPrintStr+"\n")
    def recordError(self, event):
        if not event.inaxes: return
        x, y = event.xdata, event.ydata
        toPrintStr = "!!XXXXX!! {:15.13} {:10.8}".format(x, y)
        print(toPrintStr)
        if self.logfile is not None:
            with open(self.logfile, "a") as f:
                f.write(toPrintStr+"\n")
    def on_key(self, event):
        if event.key=="c":  # normal cursor
            self._cursor = not self._cursor
            if self._cursor:
                success = self.openCursor(event)
                if not success:
                    self._cursor = not self._cursor
            else:
                self.closeCursor()
        elif event.key=="e":
            for eachax in self.data:
                d = self.data[eachax]
                if d['_legend'] == -1:
                    d['_legend'] = 1
                    d['_legendObj'] = eachax.legend(loc=0)
                    d['_legendObj'].set_visible(True)
                else:
                    d['_legend'] = not d['_legend']
                    if d['_legend'] == 1:
                        d['_legendObj'] = eachax.legend(loc=0)
                    if d['_legendObj'] is not None:
                        d['_legendObj'].set_visible(d['_legend'])
            self.fig.canvas.draw()
        elif event.key=="d": # simple cursor
            self._cursor = not self._cursor
            if self._cursor:
                self._cidcursor = self.fig.canvas.mpl_connect('motion_notify_event', self._mouse_move_simple_cursor)
                if event.inaxes:
                    d = self.data[event.inaxes]
                    if d['ly'] is None:
                        x, y = event.xdata, event.ydata
                        if x is None:
                            x = np.mean(d['ax'].get_xlim())
                            y = np.mean(d['ax'].get_ylim())
                        d['ly'] = d['ax'].axvline(x, color='k')  # the vert line
                        d['lx'] = d['ax'].axhline(y, color='k')  # the vert line
                        self.fig.canvas.draw()
            else:
                self.closeCursor()
        elif event.key=="s":
            self._doSelect = not self._doSelect
            if self._doSelect:
                self._cidselect = self.fig.canvas.mpl_connect('motion_notify_event', self._mouse_move_select)
            else:
                self.fig.canvas.mpl_disconnect(self._cidselect)
        elif event.key=="S":
            if d['_selectLine']:
                d['_selectLine'].set_lw(d['_lastLW'])
            d['ax'].set_title(d['_raw_title'])
            self.fig.canvas.draw()
        elif event.key=="w":
            self.recordPosition(event)
        elif event.key=="x":
            self.recordError(event)

    def openCursor(self, event):
        success = False
        for eachax in self.data:
            d = self.data[eachax]
            if d['x'] is not None:
                success = True
                if event.inaxes is eachax and d['ly'] is None:
                    x, y = event.xdata, event.ydata
                    if x is None:
                        x = np.mean(d['ax'].get_xlim())
                        y = np.mean(d['ax'].get_ylim())
                    d['ly'] = d['ax'].axvline(x, color='k')  # the vert line
                    d['lx'] = d['ax'].axhline(y, color='k')  # the vert line
                    self.fig.canvas.draw()
        if success:
            self._cursor = 1
            self._cidcursor = self.fig.canvas.mpl_connect(\
                    'motion_notify_event', self._mouse_move_cursor)
            self.fig.canvas.draw()
        return success

    def closeCursor(self):
        self.fig.canvas.mpl_disconnect(self._cidcursor)
        for eachax in self.data:
            d = self.data[eachax]
            if d['ly']:
                d['ly'].set_alpha(0)
                d['lx'].set_alpha(0)
                d['ax'].set_title(d['title'])
        self.fig.canvas.draw()

    def selectLine(self, event):
        x, y = event.xdata, event.ydata
        d = self.data[event.inaxes]
        minDis = np.inf
        thisSelectLine = None
        if d['_selectLine']:
            d['_selectLine'].set_lw(d['_lastLW'])
        for eachline in d['ax'].lines:
            if not eachline.__dict__.get('_selectable'):
                continue
            thisxy = eachline.get_xydata()
            mindis = ((thisxy - np.array([[x, y]]))**2).sum(axis=1).min()
            if mindis < minDis:
                minDis = mindis
                thisSelectLine = eachline
                d['x'] = eachline.get_xdata()
                d['y'] = eachline.get_ydata()
        if thisSelectLine:
            d['_selectLine'] = thisSelectLine
            d['_lastLW'] = d['_selectLine'].get_lw()
            d['_selectLine'].set_lw(3)
            d['title'] = '{}\nlabel: {}'.format(d['_raw_title'], d['_selectLine'].get_label())
            d['ax'].set_title(d['title'])
            self.fig.canvas.draw()

if __name__ == "__main__":
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    for i in range(10):
        ax1.plot(np.arange(20), np.random.rand(20), label="{}".format(i))[0].__dict__['_selectable'] = True
        ax2.plot(np.arange(20), np.random.rand(20), label="{}".format(i))[0].__dict__['_selectable'] = True

    gf = DebugFigure(fig)
