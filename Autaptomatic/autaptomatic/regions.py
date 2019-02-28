#!/usr/bin/env python

import numpy as np
import quantities as pq

import matplotlib
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar

from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

import matplotlib.cm as cm
cmap = matplotlib.cm.get_cmap('prism')

class RegionRecording:
    def __init__(self, recording, initialRegions = None):
        self.recording = recording
        self.initialRegions = initialRegions

    def makePlot(self, frame, enableAddAndDelete = True):
        self.recording.makePlot(frame)
        self.regions = Regions(self.recording.axes, self.initialRegions,
            recording = self.recording, enableAddAndDelete = enableAddAndDelete)

    def getRegions(self):
        regions = self.regions.getRegions()
        return [(x0 * self.recording.x.units, x1 * self.recording.x.units) for x0, x1 in regions]

    def getRecordings(self):
        isSelected = np.asarray(self.recording.isSelected)
        if np.all(isSelected):
            return None
        else:
            return np.where(isSelected)

class Regions:
    def __init__(self, axes, regions = None, pickradius = 5.0, recording = None, enableAddAndDelete = True):
        self.axes = axes
        self.rectangles = []
        self.pickradius = pickradius
        self.recording = recording

        self.enableAddAndDelete = enableAddAndDelete

        xmin = None
        xmax = None
        ymin = None
        ymax = None
        if regions is not None:
            for x0, x1 in regions:
                self.addRegion(x0, x1)

                mask = np.logical_and(self.recording.x >= x0, self.recording.x <= x1)
                if xmax is None:
                    xmax = x1
                else:
                    xmax = max(xmax, x1)
                if xmin is None:
                    xmin = x0
                else:
                    xmin = min(xmin, x0)

                ymin_ = np.min(self.recording.y[:, mask])
                if ymin is None:
                    ymin = ymin_
                else:
                    ymin = min(ymin, ymin_)

                ymax_ = np.max(self.recording.y[:, mask])
                if ymax is None:
                    ymax = ymax_
                else:
                    ymax = max(ymax, ymax_)

            if xmax is None:
                xmax = self.axes.get_xlim()[1]
            if xmin is None:
                xmin = self.axes.get_xlim()[0]
            if ymax is None:
                ymax = self.axes.get_ylim()[1]
            if ymin is None:
                ymin = self.axes.get_ylim()[0]

            delta = xmax - xmin
            self.axes.set_xlim(xmin - delta / 10., xmax + delta / 10.)
            delta = ymax - ymin
            self.axes.set_ylim(ymin - delta / 5., ymax + delta / 5.)

        self.press = None

        self.press_handler = self.axes.figure.canvas.mpl_connect(
            'button_press_event', self.onPress)
        self.release_heandler = self.axes.figure.canvas.mpl_connect(
            'button_release_event', self.onRelease)
        self.motion_handler = self.axes.figure.canvas.mpl_connect(
            'motion_notify_event', self.onMotion)

    def addRegion(self, x0, x1):
        ylim = self.axes.get_ylim()
        color = cmap(len(self.rectangles) * 3)
        color = (color[0], color[1], color[2], 0.4)
        # import pdb; pdb.set_recording()
        rectangle = Rectangle((x0, ylim[0]), x1 - x0, ylim[1] - ylim[0],
            color = color)
        self.rectangles.append(rectangle)
        self.axes.add_patch(rectangle)
        return rectangle

    def getRegions(self):
        return [(rectangle.get_x(), rectangle.get_x() + rectangle.get_width()) for rectangle in self.rectangles]

    def update(self):
        self.axes.figure.canvas.draw()

    def transformToData(self, mouseevent):
        transform = self.axes.get_xaxis_transform().inverted()
        xt = transform.transform_point((mouseevent.x, 0.0))[0]
        transform = self.axes.get_yaxis_transform().inverted()
        yt = transform.transform_point((0.0, mouseevent.y))[1]
        return xt, yt

    def onPress(self, e):
        if e.inaxes == self.axes:
            xt, yt = self.transformToData(e)
            transform = self.axes.get_xaxis_transform().inverted()
            pickradius = transform.transform_point((self.pickradius, 0.0))[0] \
                - transform.transform_point((0.0, 0.0))[0]
            for rectangle in self.rectangles:
                rectContains, _ = rectangle.contains(e)
                if rectContains:
                    if e.button == 3:
                        if self.enableAddAndDelete:
                            rectangle.remove()
                            self.rectangles.remove(rectangle)
                            self.update()
                    else:
                        x0 = rectangle.get_x()
                        x1 = x0 + rectangle.get_width()

                        if abs(xt - x0) < pickradius:
                            self.press = rectangle,1,
                        elif abs(xt - x1) < pickradius:
                            self.press = rectangle,2,
                        else:
                            self.press = rectangle,3,xt - x0,xt - x1,
                    return
            if self.enableAddAndDelete:
                rectangle = self.addRegion(xt, xt + 2 * pickradius)
                self.press = rectangle,2,

    def onMotion(self, e):
        if e.inaxes == self.axes:
            if self.press is not None:
                rectangle = self.press[0]
                xt, yt = self.transformToData(e)
                x0 = rectangle.get_x()
                x1 = x0 + rectangle.get_width()
                if self.press[1] == 1:
                    x0 = xt
                elif self.press[1] == 2:
                    x1 = xt
                elif self.press[1] == 3:
                    x0 = xt + self.press[2]
                    x1 = xt + self.press[3]
                rectangle.set_x(min(x0, x1))
                rectangle.set_width(abs(x1 - x0))
                self.update()

    def onRelease(self, e):
        self.press = None
        self.update()
