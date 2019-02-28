#!/usr/bin/env python

from matplotlib.lines import Line2D

class Line:
    def __init__(self, axes, x0 = None, b0 = None, b1 = None):
        self.axes = axes
        if x0 is None:
            x0 = self.axes.get_xlim()[0]
        self.x0 = x0
        if b0 is None:
            b0 = sum(axes.get_ylim()) / 2.0
        self.b0 = b0
        if b1 is None:
            b1 = 0.0
        self.b1 = b1

        xlim = (self.x0, self.axes.get_xlim()[1])
        self.artist = Line2D(xlim, (0.0, 0.0), color = 'b')
        axes.add_artist(self.artist)

        self.xr = self.x0

        self.update()

        self.press = None

        self.press_handler = self.axes.figure.canvas.mpl_connect(
            'button_press_event', self.onPress)
        self.release_heandler = self.axes.figure.canvas.mpl_connect(
            'button_release_event', self.onRelease)
        self.motion_handler = self.axes.figure.canvas.mpl_connect(
            'motion_notify_event', self.onMotion)

    def update(self):
        xlim = (self.x0, self.axes.get_xlim()[1])
        ydata = [self.b0 + self.b1 * x for x in xlim]
        self.artist.set_data(xlim, ydata)
        self.artist.figure.canvas.draw()

    def transformToData(self, mouseevent):
        transform = self.axes.get_xaxis_transform().inverted()
        xt = transform.transform_point((mouseevent.x, 0.0))[0]
        transform = self.axes.get_yaxis_transform().inverted()
        yt = transform.transform_point((0.0, mouseevent.y))[1]
        return xt, yt

    def onPress(self, e):
        if e.inaxes == self.axes:
            lineContains, point = self.artist.contains(e)
            if lineContains:
                self.press = e.button,

    def onMotion(self, e):
        if self.press is not None:
            xt, yt = self.transformToData(e)

            if self.press[0] == 3:
                yr = self.b0 + self.b1 * self.xr
                self.b1 = (yt - yr) / (xt - self.xr)
            self.b0 = yt - self.b1 * xt

            if self.press[0] == 1:
                self.xr = xt

            self.update()

    def onRelease(self, e):
        self.press = None
        self.update()
