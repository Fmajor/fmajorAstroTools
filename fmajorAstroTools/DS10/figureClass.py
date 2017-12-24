import numpy as np
import matplotlib.pyplot as plt

class DebugFigure():
    def __init__(self, figure=None, x=None):
        self.fig = figure
        self.ax = self.fig.get_axes()[0]
        self.x = x
        self._cidOnkey = self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        self._cursor_x_format=".2f"
        self._cursor_y_format=".2f"
        self._cursorStr = "x={{:{x}}}, y={{:{y}}}".\
                format(x=self._cursor_x_format, y=self._cursor_y_format)
        self._cursor = 0
    def on_key(self, event):
        if event.key=="c":  # normal cursor
            if self._cursor==0:
                print("open index")
                self._cursor = 1
                self._cidcursor = self.fig.canvas.mpl_connect(\
                        'motion_notify_event', self._mouse_move_cursor)
            else:
                print("close index")
                self._cursor = 0
                self.fig.canvas.mpl_disconnect(self._cidcursor)
                self.ax.set_xlabel("")
                self.fig.canvas.draw()

    def _mouse_move_cursor(self, event):
        if not event.inaxes: return
        xdata = event.xdata
        y = event.ydata
        indx = np.searchsorted(self.x, [xdata])[0]
        if indx<self.x.size:
            x = self.x[indx]
            self.ax.set_xlabel( self._cursorStr.format(x,y)+ " index: {}".format(indx) )
        else:
            self.ax.set_xlabel( self._cursorStr.format(self.x[-1],self.y[-1]) )
        self.fig.canvas.draw()
