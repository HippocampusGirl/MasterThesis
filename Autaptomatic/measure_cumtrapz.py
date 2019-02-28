#!/usr/bin/env python3

import neo
import quantities as pq
import numpy as np
import pandas as pd
from scipy.integrate import cumtrapz
import scipy.ndimage.filters as filters

from multiprocessing import Pool
from neo.rawio.axonrawio import StructFile

import matplotlib.pyplot as plt

import os
from os import path as op
import sys

from argparse import ArgumentParser

ap = ArgumentParser(description = "AUTAPTOMATIC")

ap.add_argument("-sucrose", action = "store_true")
ap.add_argument("-baseline", action = "store_true")
ap.add_argument("-average", action = "store_true")
ap.add_argument("input", metavar = "FILE", nargs = "+",
                    help = "input files")

args = ap.parse_args()

def get_extra_header(filename):
    header_description = [
        ("lEpochPulsePeriod", 2136, "20i"),
        ("lEpochPulseWidth", 2216, "20i"),
        ("nWaveformEnable", 2296, "2h"),
        ("nWaveformSource", 2300, "2h"),
        ("nInterEpisodeLevel", 2304, "2h"),
        ("nEpochType", 2308, "20h"),
        ("fEpochInitLevel", 2348, "20f"),
        ("fEpochLevelInc", 2428, "20f"),
        ("lEpochInitDuration", 2508, "20i"),
        ("lEpochDurationInc", 2588, "20i")
    ]
    
    header = {}
    
    with open(filename, "rb") as fid:
        f = StructFile(fid)
        for key, offset, fmt in header_description:
            val = f.read_f(fmt, offset=offset)
            if len(val) == 1:
                header[key] = val[0]
            else:
                header[key] = np.array(val)
    
    return header

def get_protocol(fname):
    io = neo.get_io(fname)
    _axon_info2 = get_extra_header(fname)
    block = io.read_block(lazy = False)
    signals = sum([s.analogsignals for s in block.segments], [])
    
    if args.sucrose:
        signals = [signals[0]]
    
    times = signals[0].times.rescale(pq.second)
    
    for s in signals:
        s.t_start = 0.0 * pq.second
    sampling_rate = signals[0].sampling_rate
        
    signals = np.asarray([np.ravel(s.as_quantity()) for s in signals]) * signals[0].units
    
    n_samples = io._axon_info["lNumSamplesPerEpisode"]
    
    digital = np.ones(n_samples) * io._axon_info["nDigitalHolding"]
    analog = np.ones(n_samples) * -70
    
    nEpochs = sum(io._axon_info["nEpochType"] != 0)
    i_last = int(n_samples * 15625 / 10 ** 6)
    for i in range(nEpochs):
        if io._axon_info["nEpochType"][i] == 1:
            i_begin = i_last
            i_end = i_last + io._axon_info["lEpochInitDuration"][i] 
                # + io._axon_info["lEpochDurationInc"][i] * epiNum # ignore for now
            dif = i_end - i_begin
            analog[i_begin:i_end] += np.ones((dif)) * \
                io._axon_info["fEpochInitLevel"][i]
            digital[i_begin:i_end] += np.ones((dif)) * \
                io._axon_info["nDigitalValue"][i]
                # + epoch["fEpochLevelInc"] * epiNum) # ignore for now
            i_last = i_end
        elif io._axon_info["nEpochType"][i] == 3:
            i_duration = io._axon_info["lEpochInitDuration"][i] 
            i_period = _axon_info2["lEpochPulsePeriod"][i]
            dif = _axon_info2["lEpochPulseWidth"][i]
            for i_offset in range(0, i_duration, i_period):
                i_begin = i_last + i_offset
                i_end = i_begin + dif
                analog[i_begin:i_end] += np.ones((dif)) * \
                    io._axon_info["fEpochInitLevel"][i]
            i_last += i_duration
    
    i_start = np.where(np.diff(analog) > 0)
    t_start = (i_start / sampling_rate).simplified
    
    t_end = t_start * 1.0
    if t_start.size > 1:
        t_end += np.min(np.diff(t_start))
    else:
        dur = (times[-1] - t_start)
        t_end += dur
    
    return t_start, t_end

##
##

times_ = None

colNames = ["NAME", "FILENAME", "PROTOCOL", "TIMESTAMP", "SERIES", "DATA"]
nMeta = len(colNames) - 1
dtypes = ["U64", "U64", "U64", "U64", "i4", "f4"]
columns = None

for ii, fname in enumerate(args.input):
    try:
        name = os.path.splitext(os.path.basename(fname))[0]
        oo = os.path.join(os.path.dirname(fname), "%s.regions" % name)

        print("%d/%d %s" % (ii, len(args.input), name))
        
        if not os.path.isfile(oo):
            continue

        t_start, t_end = get_protocol(fname)

        checkedRegions = np.loadtxt(oo).tolist()
        if type(checkedRegions[0]) is not list:
            checkedRegions = [checkedRegions]

        io = neo.get_io(fname)
        block = io.read_block(lazy = False)
        signals = sum([s.analogsignals for s in block.segments], [])
        
        if args.sucrose:
            signals = [signals[0]]

        for s in signals:
            s.t_start = 0.0 * pq.second

        sampling_rate = signals[0].sampling_rate

        times = signals[0].times.rescale(pq.second)
        
        if times_ is None:
            times_ = times
        else:
            if not np.array_equal(times, times_):
                continue
                
        # # import pdb; pdb.set_trace()
        # sigma = 2e5 * 0.53002 / 4 / 1e6 * float(sampling_rate)
        # smoothed_signals = [filters.gaussian_filter(s, sigma) * s.units for s in signals]
        # 
        # # high pass filter
        # filtered_signals = [s - t for s, t in zip(signals, smoothed_signals)]

        signals = np.asarray([np.ravel(s.as_quantity()) for s in signals]) * signals[0].units
        
        if args.average:
            signals = np.mean(signals, axis = 0)[None, :]
            
        indices = np.asarray(np.zeros_like(times, dtype = np.bool))
        
        for region in checkedRegions:
            indices = np.logical_or(indices, np.logical_and(times >= region[0], times <= region[1]))
            
        if args.baseline is not None:
            # i_baseline = int(((np.ravel(t_start)[0] - 5 * pq.millisecond) * sampling_rate).simplified)
            # baseline = np.mean(signals[:, :i_baseline], axis = 1)[:, None] 
            i_baseline = int(((checkedRegions[-1][-1] * pq.second + 1 * pq.second) * sampling_rate).simplified)
            baseline = np.mean(signals[:, i_baseline:(i_baseline+int(sampling_rate))], axis = 1)[:, None] 
            signals -= baseline
        # import pdb; pdb.set_trace()

        indices = np.logical_or(indices, times >= checkedRegions[-1][1])

        regionsignals = signals[:, indices]
        
        # plt.plot(regionsignals.T)
        # plt.show()
        
        ct = (cumtrapz(regionsignals, x = times[:indices.sum()]) * regionsignals.units * times.units).rescale("pA*s")
        ct = np.asarray(ct)[:, ::10]
        
        # # import pdb; pdb.set_trace()
        # plt.plot(ct.T)
        # plt.show()
        
        if columns is None:
            columns = [[] for i in range(ct.shape[1])]
        
        for i in range(ct.shape[0]):
            columns[0].append(name)
            columns[1].append(fname)
            columns[2].append(os.path.basename(io._axon_info["sProtocolPath"].decode()))
            columns[3].append(io._axon_info["rec_datetime"].strftime("%H:%M:%S.%f"))
            columns[4].append(i)
            for j in range(nMeta, len(columns)):
                if j < ct.shape[1]:
                    columns[j].append(ct[i, j - nMeta])
                else:
                    columns[j].append(np.nan)
    except ValueError:
        pass

nData = len(columns) - nMeta

header = colNames[:-1] + ["%s_%d" % (colNames[-1], i) for i in range(nData)]

dtypes = dtypes[:-1] + [dtypes[-1] for i in range(nData)]

# import pdb; pdb.set_trace()
data = np.array(list(zip(*columns)), dtype = list(zip(header, dtypes)))

np.savetxt("cumtrapz.csv", data, 
    fmt = ["%s" for i in range(nMeta)] + ["%.18e" for i in range(nData)],
    header = " ".join(header), comments = "")

# dataFrame = pd.DataFrame({key:value for key, value in zip(colNames, zip(*rows))}, columns = colNames)
# dataFrame.to_csv("cumtrapz.csv")





