#!/usr/bin/env python

import quantities as pq
import numpy as np

import os
import sys
    
import argparse

import wx
from .frame import Frame
from .regions import RegionRecording
from .recording import Recording
from .events import MeasureEvents
from .measure import get_measure

def main():
    parser = argparse.ArgumentParser(description = 'AUTAPTOMATIC')

    parser.add_argument('-average', action = 'store_true',
                        help = 'average selected recordings')
    parser.add_argument('-baseline', metavar = ('START', 'END'), nargs = 2,
                        type = float, help = 'specify time window [s]')
    parser.add_argument('-cut', metavar = ('START', 'END'), nargs = 2,
                        type = float, help = 'specify time window [s]')

    analysis = parser.add_argument_group(title = 'analysis type')
    groupgroup = analysis.add_mutually_exclusive_group(required = True)
    groupgroup.add_argument('-regions', action = 'store_true')
    groupgroup.add_argument('-events', action = 'store_true')
    groupgroup.add_argument('-sucrose', action = 'store_true')

    action = parser.add_argument_group(title = 'action')
    groupgroup = action.add_mutually_exclusive_group(required = True)
    groupgroup.add_argument('-specify', action = 'store_true')
    groupgroup.add_argument('-check', action = 'store_true')
    groupgroup.add_argument('-measure', nargs = '+')

    parser.add_argument('input', metavar = 'FILE', nargs = '+',
                        help = 'input files')

    args = parser.parse_args()

    def writeProgress(i, n):
        progress = float(i + 1) / float(n) * 100.
        sys.stderr.write('\r%f%%' % progress)
        sys.stderr.flush()

    recordings = []

    if args.events:
        rows = []
        colNames = []
        import pandas as pd
        for ii, f in enumerate(args.input):
            name = os.path.splitext(os.path.basename(f))[0]
            amp, charge, freq = MeasureEvents(f)
            if amp is not None:
                colNames = ['NAME', 'AMPLITUDE_%s' % amp.units, 'CHARGE_%s' % charge.units, 'FREQUENCY_%s' % freq.units]
                print(name, float(amp), float(charge), float(freq))
                rows.append([name, float(amp), float(charge), float(freq)])
        dataFrame = pd.DataFrame({key:value for key, value in zip(colNames, zip(*rows))}, columns = colNames)
        dataFrame.to_excel('measurements.xlsx', sheet_name = 'measurements')
    else:
        import neo
        
        times = None
        for ii, f in enumerate(args.input):
            name = os.path.splitext(os.path.basename(f))[0]
            try:
                io = neo.get_io(f)
            except ValueError:
                print(f)
            block = io.read_block(lazy = False)
            signals = sum([s.analogsignals for s in block.segments], [])
            for s in signals:
                s.t_start = 0.0 * pq.second
            isOk = True
            if args.cut:
                try:
                    signals = [s.time_slice(args.cut[0] * pq.second, args.cut[1] * pq.second)
                        for s in signals]
                except:
                    pass
            for s in signals:
                times_ = s.times.rescale(pq.second)
                if times is None:
                    times = times_
                else:
                    if not np.array_equal(times, times_):
                        recordings.append(None)
                        sys.stderr.write('\rEXCLUDE:%s\n' % name)
                        sys.stderr.flush()
                        isOk = False
                        break
            if isOk:
                signals = np.asarray([np.ravel(s.as_quantity()) for s in signals]) * signals[0].units
                recordings.append(Recording(name = f, x = times, y = signals, io = io))
            writeProgress(ii, len(args.input))

        args.input = [inputfile for i, inputfile in enumerate(args.input) if recordings[i] is not None]
        recordings = [recording for recording in recordings if recording is not None]

        if args.regions:
            if args.specify:
                avg = np.mean(np.concatenate([recording.y for recording in recordings], axis = 0), axis = 0) \
                    * recordings[0].y.units
                avgRecording = Recording(name = "Average", x = times, y = avg[None, :], enableSelection = False)

                app = wx.App(False)
                aframe = Frame()

                psc = RegionRecording(avgRecording)
                psc.makePlot(aframe)

                aframe.Show()
                app.SetTopWindow(aframe)
                app.MainLoop()

                initialRegions = psc.getRegions()

                for f in args.input:
                    name = os.path.splitext(os.path.basename(f))[0]
                    oo = os.path.join(os.path.dirname(f), '%s.regions' % name)
                    np.savetxt(oo, np.asarray(initialRegions))

            elif args.check:
                app = wx.App(False)
                aframe = Frame()

                pscs = []
                for f, recording in zip(args.input, recordings):
                    initialRegions = None

                    name = os.path.splitext(os.path.basename(f))[0] #+ \
                        #os.path.basename(recording.io._axon_info["sProtocolPath"].decode())
                    oo = os.path.join(os.path.dirname(f), '%s.regions' % name)

                    if os.path.exists(oo):
                        initialRegions = np.loadtxt(oo).tolist()
                        if type(initialRegions[0]) is not list:
                            initialRegions = [initialRegions]

                    oo = os.path.join(os.path.dirname(f), '%s.selection' % name)

                    if os.path.exists(oo):
                        selection = np.loadtxt(oo)
                        recording.isSelected = [i in selection for i in range(recording.y.shape[0])]

                    psc = RegionRecording(recording, initialRegions = initialRegions)
                    psc.makePlot(aframe, enableAddAndDelete = False)
                    pscs.append(psc)

                aframe.Show()
                app.SetTopWindow(aframe)
                app.MainLoop()

                for f, psc in zip(args.input, pscs):
                    name = os.path.splitext(os.path.basename(f))[0]

                    oo = os.path.join(os.path.dirname(f), '%s.regions' % name)
                    checkedRegions = psc.getRegions()
                    np.savetxt(oo, np.asarray(checkedRegions))

                    selected = psc.getRecordings()
                    if selected is not None:
                        oo = os.path.join(os.path.dirname(f), '%s.selection' % name)
                        np.savetxt(oo, selected, '%02d')
        elif args.sucrose:
            
            from autaptomatic.sucrose import SucroseRecording

            if args.specify:
                for i, (f, recording) in enumerate(zip(args.input, recordings)):
                    suc = SucroseRecording(recording)
                    suc.fit()

                    name = os.path.splitext(os.path.basename(f))[0]
                    oo = os.path.join(os.path.dirname(f), '%s.sucrose' % name)
                    np.savetxt(oo, np.asarray(suc.getParams()))

                    writeProgress(i, len(args.input))

            elif args.check:
                app = wx.App(False)
                aframe = Frame()

                sucroses = []
                for f, recording in zip(args.input, recordings):
                    name = os.path.splitext(os.path.basename(f))[0]
                    oo = os.path.join(os.path.dirname(f), '%s.sucrose' % name)
                    params = np.loadtxt(oo)

                    recording.enableSelection = False
                    sucrose = SucroseRecording(recording, x0 = params[0], x1 = params[3],
                        b0 = params[1], b1 = params[2])
                    sucrose.makePlot(aframe)
                    sucroses.append(sucrose)

                aframe.Show()
                app.SetTopWindow(aframe)
                app.MainLoop()

                for f, sucrose in zip(args.input, sucroses):
                    name = os.path.splitext(os.path.basename(f))[0]
                    oo = os.path.join(os.path.dirname(f), '%s.sucrose' % name)
                    np.savetxt(oo, np.asarray(sucrose.getParams()))
        if args.measure:
            for f, recording in zip(args.input, recordings):
                name = os.path.splitext(os.path.basename(f))[0]
                oo = os.path.join(os.path.dirname(f), '%s.selection' % name)

                if os.path.exists(oo):
                    selection = np.loadtxt(oo).astype(int)
                    # recording.isSelected = [i in selection for i in range(recording.y.shape[0])]
                    recording.y = recording.y[selection, :]

            for recording in recordings:
                if len(recording.y.shape) == 1:
                    recording.y = recording.y[None, :]

            if args.average:
                for recording in recordings:
                    if recording.y.shape[0] > 1:
                        recording.y = np.mean(recording.y, axis = 0)[None, :]

            import pandas as pd
            colNames = None
            rows = []

            if args.regions:

                if args.baseline is not None:
                    x0, x1 = args.baseline
                    for recording in recordings:
                        recording.baseline = np.mean(recording.y[:, np.logical_and(recording.x >= x0, recording.x <= x1)], axis = 1)[:, None]
                        recording.y -= recording.baseline

                for ii, (f, recording) in enumerate(zip(args.input, recordings)):
                    name = os.path.splitext(os.path.basename(f))[0]
                    oo = os.path.join(os.path.dirname(f), '%s.regions' % name)
                    
                    if not os.path.isfile(oo):
                        continue
                        
                    checkedRegions = np.loadtxt(oo).tolist()
                    if type(checkedRegions[0]) is not list:
                        checkedRegions = [checkedRegions]
                    
                    colNames_ = ['NAME', 'PROTOCOL', 'TIMESTAMP', 'BASELINE %s' % np.mean(recording.baseline).units]
                    row = [f, 
                        os.path.basename(recording.io._axon_info["sProtocolPath"].decode()),
                        recording.io._axon_info["rec_datetime"].strftime("%H:%M:%S.%f"), float(np.mean(recording.baseline))]
                    for measure in args.measure:
                        for i, (x0, x1) in enumerate(checkedRegions):
                            mask = np.logical_and(recording.x >= x0, recording.x <= x1)
                            yy = recording.y[:, mask]

                            for j in range(yy.shape[0]):
                                value = get_measure(measure, recording, mask, yy[j, :])
                                colNames_.append('%s_REGION_%02d_SERIES_%02d_%s' % (measure, i, j, value.units))
                                row.append(float(value))
                    if colNames is None:
                        colNames = colNames_

                    # if np.all(np.asarray([a == b for a, b in zip(colNames, colNames_)])):
                    while len(row) < len(colNames):
                        row.append(np.nan)
                    rows.append(row)

                    writeProgress(ii, len(args.input))

            elif args.sucrose:
                for ii, (f, recording) in enumerate(zip(args.input, recordings)):
                    name = os.path.splitext(os.path.basename(f))[0]
                    oo = os.path.join(os.path.dirname(f), '%s.sucrose' % name)
                    params = np.loadtxt(oo)

                    colNames_ = ['NAME', 'PROTOCOL', 'TIMESTAMP']
                    row = [f, 
                        os.path.basename(recording.io._axon_info["sProtocolPath"].decode()),
                        recording.io._axon_info["rec_datetime"].strftime("%H:%M:%S.%f")]
                    for measure in args.measure:
                        x0 = params[0]
                        x1 = params[3] 
                        mask = np.logical_and(recording.x >= x0, recording.x <= x1)
                        yy = recording.y[0, mask][None, :]
                        yy -= (params[1] + np.asarray(recording.x[mask] * params[2])) * yy.units

                        if args.baseline is not None:
                            xx0, xx1 = args.baseline
                            xx0 += x0
                            xx1 += x0
                            yy -= np.mean(recording.y[:, np.logical_and(recording.x >= xx0, recording.x <= xx1)], axis = 1)[:, None]

                        for j in range(yy.shape[0]):                            
                            value = get_measure(measure, recording, mask, yy[j, :])
                            
                            colNames_.append('%s_SERIES_%02d_%s' % (measure,  j, value.units))
                            row.append(float(value))
                    if colNames is None:
                        colNames = colNames_

                    assert np.all(np.asarray([a == b for a, b in zip(colNames, colNames_)]))

                    rows.append(row)
                    writeProgress(ii, len(args.input))

            dataFrame = pd.DataFrame({key:value for key, value in zip(colNames, zip(*rows))}, columns = colNames)
            dataFrame.to_excel('measurements.xlsx', sheet_name = 'measurements')
