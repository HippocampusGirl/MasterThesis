#!/usr/bin/env python

import numpy as np
import quantities as pq
import scipy.optimize as optimize
import scipy.ndimage as ndimage

from random import random

from autaptomatic.regions import Regions
from autaptomatic.line import Line

def logLikelihood(y, x, s):
    return -len(y) / 2.0 * np.log(s) - np.sum(np.square(y - x)) / (2.0 * s)

def curve(x, a, q, p):
    return -1.0 * a * np.power(x, q) * np.exp(-x / p)

def fitCurve(y, m = 0.0, x = None, niter = 2):
    def cost(solution):
        a = solution[0]
        q = solution[1]
        p = solution[2]
        b = solution[3]
        n = solution[4]
        s = solution[5]

        yhat = curve(x, a, q, p) + b * x + n + m
        return -logLikelihood(y, yhat, s)

    solution = np.asarray([10.0, 1.0, 50.0, 0.0, 0.0, 1.0])
    bnds = [(1.0, None), (1.0, None), (10.0, None), (None, None), (None, None), (1.0, None)]
    res = optimize.basinhopping(cost, solution, niter = niter,
        minimizer_kwargs = {'bounds': bnds, 'options': {'disp': False}})

    return res.x

def acceptanceProbability(diff, temperature):
    return np.exp(-diff / temperature)

def anneal(cost, neighbor, solution):
    old_cost = cost(solution)
    T = 1e4
    T_min = 0.00001
    alpha = 0.8
    while T > T_min:
        new_solution = neighbor(solution)
        new_cost = cost(new_solution)
        diff = new_cost - old_cost
        if diff < 0 or acceptanceProbability(diff, T) > random():
            solution = new_solution
            old_cost = new_cost
        T *= alpha
    return solution, old_cost

class SucroseRecording:
    def __init__(self, recording, x0 = None, x1 = None, b0 = None, b1 = None):
        self.recording = recording
        self.x0 = x0
        self.x1 = x1
        self.b0 = b0
        self.b1 = b1
        self.line = None
        self.regions = None

    def makePlot(self, frame):
        self.recording.makePlot(frame)
        self.regions = Regions(self.recording.axes, [(float(self.x0), float(self.x1))],
            enableAddAndDelete = False, recording = self.recording)
        self.line = Line(self.recording.axes, b0 = float(self.b0), b1 = float(self.b1))

    def getParams(self):
        if self.line is not None:
            self.b0 = self.line.b0
            self.b1 = self.line.b1
        if self.regions is not None:
            self.x0 = self.regions.getRegions()[0][0]
            self.x1 = self.regions.getRegions()[0][1]
        return self.x0, self.b0, self.b1, self.x1

    def fit(self, subsampling = 100):
        y = np.ravel(self.recording.y[0, :])
        y = ndimage.filters.gaussian_filter1d(y, 1000.0 / (2.0 * np.sqrt(2.0 * np.log(2.0))))
        y = np.mean(np.reshape(y, (-1, subsampling)), axis = 1)

        m = np.mean(y[:20])

        def neighbor(solution):
            solution = solution + 10 * np.random.normal(size = solution.size)

            if solution[0] >= len(y) - 1:
                solution[0] = len(y) - 2

            if solution[0] <= 1:
                solution[0] = 1

            if solution[1] <= 1e-6:
                solution[1] = 1e-6

            return solution

        def cost(solution):
            b = int(solution[0])
            s = solution[1]

            yy = y[b:]
            x = np.arange(len(yy))
            curveParameters = fitCurve(yy, m = m, x = x)

            a = curveParameters[0]
            q = curveParameters[1]
            p = curveParameters[2]
            bb = curveParameters[3]
            n = curveParameters[4]
            yhat = curve(x, a, q, p) + bb * x + n + m

            return -logLikelihood(y[:b], m, s) + -logLikelihood(yy, yhat, s)

        solution = np.asarray([300.0, 1.0])
        solution, _ = anneal(cost, neighbor, solution)

        yy = y[int(solution[0]):]
        x = np.arange(len(yy))
        curveParameters = fitCurve(yy, m = m, x = x, niter = 10)

        a = curveParameters[0]
        q = curveParameters[1]
        p = curveParameters[2]
        bb = curveParameters[3]
        n = curveParameters[4]

        xx = np.asarray(self.recording.x[int(solution[0] * subsampling)::subsampling])

        self.b1 = (bb * x[0] - bb * x[len(x)-1]) / (xx[0] - xx[len(xx)-1])
        self.b0 = (bb * x[0] + m + n) - self.b1 * xx[0]

        self.x0 = self.recording.x[int(solution[0] * subsampling)]
        self.x1 = np.max(self.recording.x)
