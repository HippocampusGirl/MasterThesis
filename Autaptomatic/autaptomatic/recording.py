#!/usr/bin/env python

import numpy as np

import matplotlib
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar

from matplotlib.figure import Figure
from matplotlib.lines import Line2D

import matplotlib.cm as cm
cmap = matplotlib.cm.get_cmap('prism')

import wx

class Recording:
    def __init__(self, name = None, x = None, y = None, io = None, enableSelection = True):
        self.name = name
        self.x = x
        self.y = y
        self.baseline = np.nan
        self.io = io
        self.enableSelection = enableSelection
        self.isSelected = [True for i in range(self.y.shape[0])]

    def makePlot(self, frame):
        self.figure = Figure()
        self.canvas = FigureCanvasWxAgg(frame.scrolled_panel, -1, self.figure)
        frame.sizer.Add(self.canvas, 1.0, wx.EXPAND)

        self.navigation_toolbar = NavigationToolbar(self.canvas)

        self.axes = self.figure.add_subplot(111)
        self.axes.set_title(self.name)
        self.axes.set_xlabel(self.x.dimensionality)
        self.axes.set_ylabel(self.y.dimensionality)

        self.lines = []
        for i in range(self.y.shape[0]):
            line = Line2D(self.x, self.y[i, :])
            if self.isSelected[i]:
                line.set_color(cmap(i * 3))
            else:
                line.set_color((0.2, 0.2, 0.2, 0.5))
            self.lines.append(line)
            self.axes.add_line(line)

        self.axes.autoscale_view()
        frame.update()

        if self.enableSelection:
            self.press_handler = self.figure.canvas.mpl_connect(
                'button_press_event', self.onPress)

    def onPress(self, e):
        if e.inaxes == self.axes:
            for i, line in enumerate(self.lines):
                lineContains, _ = line.contains(e)
                if lineContains:
                    self.isSelected[i] = not self.isSelected[i]
                    if self.isSelected[i]:
                        line.set_color(cmap(i * 3))
                    else:
                        line.set_color((0.2, 0.2, 0.2, 0.5))
                    return
