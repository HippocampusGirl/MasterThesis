#!/usr/bin/env python3

import neo
import quantities as pq
import numpy as np
import pandas as pd

from multiprocessing import Pool

from neo.rawio.axonrawio import StructFile

import os
from os import path as op
import sys

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


from argparse import ArgumentParser

ap = ArgumentParser(description = "AUTAPTOMATIC")

ap.add_argument("-sucrose", action = "store_true")
ap.add_argument("input", metavar = "FILE", nargs = "+",
                    help = "input files")

args = ap.parse_args()

def process(fname):
    print(fname)
    
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
    
    if args.sucrose:
        n_samples /= 2
        n_samples = int(n_samples)
    
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
        if dur > 25 * pq.millisecond:
            dur = 25 * pq.millisecond
        t_end += dur
    
    t_start += 3 * pq.millisecond
    
    # import pdb; pdb.set_trace()
    fprefix, _ = op.splitext(fname)
    np.savetxt(fprefix + ".regions", np.vstack((t_start, t_end)).T)
    
# from matplotlib import pyplot as plt
# plt.plot(np.arange(n_samples), analog)
# # plt.plot(np.arange(n_samples), np.ravel(signals))
# plt.plot(np.arange(n_samples), signals.mean(0))
# # plt.plot(np.arange(n_samples), digital)
# plt.show()
    
    # print(io)

for i in args.input:
    try:
        process(i)
    except ValueError:
        pass

# rows = Pool().map(process, args.input)
# rows = sum(rows, [])