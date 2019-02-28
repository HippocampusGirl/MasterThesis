#!/usr/bin/env python3

import numpy as np
import quantities as pq
import neo
from scipy import signal

from neo.rawio.axonrawio import StructFile

import matplotlib.pyplot as plt

import os

global sucrose
sucrose = False

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
        ("lEpochDurationInc", 2588, "20i"),
        ("nDigitalTrainValue", 2682, "20h")
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
    global sucrose

    io = neo.get_io(fname)
    _axon_info2 = get_extra_header(fname)
    block = io.read_block(lazy = False)
    signals = sum([s.analogsignals for s in block.segments], [])
    
    if sucrose:
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
                _axon_info2["nDigitalTrainValue"][i]
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
    
    return analog, digital

start = 30000
end = 100000

# sucrose = True
# 
# f = "/Users/lw/Dropbox/Univ/masterThesis/data/Ephys/20180726/0012_0007.abf"
# io = neo.get_io(f)
# block = io.read_block()
# scr = np.ravel(sum([s.analogsignals for s in block.segments], []))
# scr = scr[start:end] - np.mean(scr[start:start+5000])
# 
# a, d = get_protocol(f)
# 
# # f = "/Users/lw/Dropbox/Univ/masterThesis/data/Ephys/20180803B/0081_0004.abf"
# # io = neo.get_io(f)
# # block = io.read_block()
# # kd = np.ravel(sum([s.analogsignals for s in block.segments], []))
# # kd = kd[start:end] - np.mean(kd[start:start+5000])
# 
# d[np.abs(d - d[0]) < 0.5] = np.nan
# 
# f, (ax1, ax2) = plt.subplots(2, sharex = True, sharey = True)
# ax1.plot(scr[::10], linestyle = "-", linewidth = 0.33)
# ax1.plot(d[start:end][::10], linestyle = "-", linewidth = 0.33)
# # ax2.plot(kd, linestyle = "-", linewidth = 0.33)
# ax2.plot([0, 1000], [0, 0], 'k-', lw=2)
# ax2.plot([0, 0], [0, -250], 'k-', lw=2)
# 
# f.savefig("figure3a.pdf")
# 
def loadavg1(f):
    name = os.path.splitext(os.path.basename(f))[0]
    io = neo.get_io(f)
    block = io.read_block()
    signals = sum([s.analogsignals for s in block.segments], [])
    for s in signals:
        s.t_start = 0.0 * pq.second
    oo = os.path.join(os.path.dirname(f), '%s.regions' % name)
    regions = np.loadtxt(oo).tolist()
    region = regions[0]
    
    a, d = get_protocol(f)
    a[np.abs(a - a[0]) < 0.5] = np.nan

    oo = os.path.join(os.path.dirname(f), '%s.selection' % name)
    if os.path.exists(oo):
        selection = np.loadtxt(oo)
        signals = [signals[int(i)] for i in selection]

    times = None
    for s in signals:
        times_ = s.times.rescale(pq.second)
        if times is None:
            times = times_
        else:
            assert(np.array_equal(times, times_))

    signals = np.asarray([np.ravel(s.as_quantity()) for s in signals]) * signals[0].units
    signals = signals - np.mean(signals[:, :1500], axis = 1)[:, None]
    signals = np.mean(signals, axis = 0)
    
    indices = np.logical_and(times > region[0], times < region[1])
    indices[:np.min(np.where(signals[indices] < 0)[0])] = False
    
    first = np.where(indices)[0][0]
    indices[first-50:first] = True
        
    # import pdb; pdb.set_trace()

    signals = signals[indices]
    signals[:50] = 0 * signals.units
    # signals = np.concatenate((np.full(50, 0), signals))
    
    signals[signals > 0] = 0 * signals.units
    
    a = a[indices]

    return signals, a

# f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex = True, sharey = True)
# with open("figure6a.txt", "r") as list:
#     listitems = list.readlines()
#     def plotitem(ax, item):
#         scr = loadavg1(os.path.join("/Users/lw/Dropbox/Univ/masterThesis/data/Ephys/", item.strip()))
#         ax.plot(scr, linestyle = "-", linewidth = 0.33)
#     plotitem(ax1, listitems[0])
#     plotitem(ax1, listitems[1])
#     plotitem(ax2, listitems[2])
#     plotitem(ax2, listitems[3])
#     plotitem(ax3, listitems[4])
#     plotitem(ax3, listitems[5])
# ax4.plot([0, 100], [0, 0], 'k-', lw=1)
# ax4.plot([0, 0], [0, -1000], 'k-', lw=1)
# f.savefig("figure6a.pdf")

# f, (ax1, ax2, ax4) = plt.subplots(3, sharex = True, sharey = True)
# with open("figure1a.txt", "r") as list:
#     listitems = list.readlines()
#     def plotitem(ax, item):
#         f = os.path.join("/Users/lw/Dropbox/Univ/masterThesis/data/Ephys/", item.strip())
#         scr, a = loadavg1(f)
#         ax.plot(scr, linestyle = "-", linewidth = 0.33)
#         ax.plot(a, linestyle = "-", linewidth = 0.33)
#         # import pdb; pdb.set_trace()
#     plotitem(ax1, listitems[0])
#     plotitem(ax1, listitems[1])
#     # plotitem(ax2, listitems[2])
#     # plotitem(ax2, listitems[3])
# ax4.plot([0, 100], [0, 0], 'k-', lw=1)
# ax4.plot([0, 0], [0, -1000], 'k-', lw=1)
# f.savefig("figure1a.pdf")

# # 
# 
# f = "/Users/lw/Dropbox/Univ/masterThesis/data/Ephys/20180802A/0045_0000.abf"
# scr = loadavg1(f)
# f = "/Users/lw/Dropbox/Univ/masterThesis/data/Ephys/20180802A/0045_0001.abf"
# nbqx = loadavg1(f)
# 
# f, (ax1, ax2) = plt.subplots(2, sharex = True, sharey = True)
# ax1.plot(scr, linestyle = "-", linewidth = 0.33)
# ax1.plot(nbqx, linestyle = "-", linewidth = 0.33)
# ax2.plot([0, 100], [0, 0], 'k-', lw=1)
# ax2.plot([0, 0], [0, -1000], 'k-', lw=1)
# 
# f.savefig("figure1a.pdf")
# 
def loadavg2(f):
    name = os.path.splitext(os.path.basename(f))[0]
    io = neo.get_io(f)
    block = io.read_block()
    signals = sum([s.analogsignals for s in block.segments], [])
    for s in signals:
        s.t_start = 0.0 * pq.second
    oo = os.path.join(os.path.dirname(f), '%s.regions' % name)
    regions = np.loadtxt(oo).tolist()
    
    a, d = get_protocol(f)
    a[np.abs(a - a[0]) < 0.5] = np.nan

    oo = os.path.join(os.path.dirname(f), '%s.selection' % name)
    if os.path.exists(oo):
        selection = np.loadtxt(oo)
        signals = [signals[int(i)] for i in selection]

    times = None
    for s in signals:
        times_ = s.times.rescale(pq.second)
        if times is None:
            times = times_
        else:
            assert(np.array_equal(times, times_))

    signals = np.asarray([np.ravel(s.as_quantity()) for s in signals]) * signals[0].units
    signals = signals - np.mean(signals[:, :1500], axis = 1)[:, None]
    signals = np.mean(signals, axis = 0)

    last = 0
    
    signalsout = None
    aout = None
    
    for region in regions:
        indices = np.zeros(signals.shape, dtype = bool)
        indices[np.logical_and(times > region[0], times < region[1])] = True
        # indices[:np.min(np.where(signals[indices] < 0)[0])] = False
        
        extraindices = np.zeros(signals.shape, dtype = bool)
        first = np.where(indices)[0][0]
        extraindices[first-50:first] = True
        if last > 0:
            extraindices[last:first] = True
        extraindices = np.logical_and(extraindices, np.logical_not(indices))
        
        last = np.where(indices)[0][-1]
        
        allindices = np.logical_or(extraindices, indices)
        
        signals_ = signals[allindices]
        signals_[extraindices[allindices]] = 0 * signals.units
        signals_ -= signals[-1]
        signals_[extraindices[allindices]] = 0 * signals.units
        signals_[signals_ > 0] = 0 * signals.units
        
        a_ = a[allindices]
        
        if signalsout is None:
            signalsout = signals_
            aout = a_
        else:
            signalsout = np.hstack((signalsout, signals_))
            aout = np.hstack((aout, a_))
        
    return signalsout, aout

# f, (ax1, ax2, ax3) = plt.subplots(3, sharex = True, sharey = True)
# with open("figure4a.txt", "r") as list:
#     listitems = list.readlines()
#     def plotitem(ax, item):
#         scr, a = loadavg2(os.path.join("/Users/lw/Dropbox/Univ/masterThesis/data/Ephys/", item.strip()))
#         ax.plot(scr, linestyle = "-", linewidth = 0.33)
#         ax.plot(a, linestyle = "-", linewidth = 0.33)
#     plotitem(ax1, listitems[0])
#     # plotitem(ax2, listitems[1])
# ax3.plot([0, 100], [0, 0], 'k-', lw=1)
# ax3.plot([0, 0], [0, -1000], 'k-', lw=1)
# 
# f.savefig("figure4a.pdf")




# f, (ax1, ax2) = plt.subplots(2, sharex = True, sharey = True)
# with open("figure2a.txt", "r") as list:
#     listitems = list.readlines()
#     for item in listitems:
#         scr = loadavg1(os.path.join("/Users/lw/Dropbox/Univ/masterThesis/data/Ephys/", item.strip()))
#         ax1.plot(scr, linestyle = "-", linewidth = 0.33)
# ax2.plot([0, 100], [0, 0], 'k-', lw=1)
# ax2.plot([0, 0], [0, -1000], 'k-', lw=1)
# 
# f.savefig("figure2a.pdf")

# f, (ax1, ax2) = plt.subplots(2, sharex = True, sharey = True)
# with open("figure2a.txt", "r") as list:
#     listitems = list.readlines()
#     for item in listitems:
#         scr = loadavg1(os.path.join("/Users/lw/Dropbox/Univ/masterThesis/data/Ephys/", item.strip()))
#         ax1.plot(-scr / scr.min(), linestyle = "-", linewidth = 0.33)
# ax2.plot([0, 100], [0, 0], 'k-', lw=1)
# ax2.plot([0, 0], [0, -1], 'k-', lw=1)
# 
# f.savefig("figure2a.pdf")
# f, (ax1, ax2) = plt.subplots(2, sharex = True, sharey = True)
# with open("figure2a_.txt", "r") as list:
#     listitems = list.readlines()
#     for item in listitems:
#         scr = loadavg1(os.path.join("/Users/lw/Dropbox/Univ/masterThesis/data/Ephys/", item.strip()))
#         ax1.plot(-scr / scr.min(), linestyle = "-", linewidth = 0.33)
# ax2.plot([0, 100], [0, 0], 'k-', lw=1)
# ax2.plot([0, 0], [0, -1], 'k-', lw=1)
# 
# f.savefig("figure2a_.pdf")


# 
# def loadavg3(f):
#     name = os.path.splitext(os.path.basename(f))[0]
#     io = neo.get_io(f)
#     block = io.read_block(lazy = False, cascade = True)
#     signals = sum([s.analogsignals for s in block.segments], [])
#     for s in signals:
#         s.t_start = 0.0 * pq.second
#     oo = os.path.join(os.path.dirname(f), '%s.regions' % name)
#     regions = np.loadtxt(oo).tolist()
# 
#     oo = os.path.join(os.path.dirname(f), '%s.selection' % name)
#     if os.path.exists(oo):
#         selection = np.loadtxt(oo)
#         signals = [signals[int(i)] for i in selection]
# 
#     times = None
#     for s in signals:
#         times_ = s.times.rescale(pq.second)
#         if times is None:
#             times = times_
#         else:
#             assert(np.array_equal(times, times_))
# 
#     signals = np.asarray([np.ravel(s.as_quantity()) for s in signals]) * signals[0].units
#     signals = signals - np.mean(signals[:, :1500], axis = 1)[:, None]
#     signals = np.mean(signals, axis = 0)
#     signals0 = signals[np.logical_and(times > regions[0][0], times < regions[0][1])]
#     signals1 = signals[np.logical_and(times > regions[1][0], times < regions[1][1])]
#     signals2 = signals[np.logical_and(times > regions[2][0], times < regions[2][1])]
#     signals3 = signals[np.logical_and(times > regions[3][0], times < regions[3][1])]
#     signals4 = signals[np.logical_and(times > regions[4][0], times < regions[4][1])]
# 
#     signals0 -= np.mean(signals0[len(signals0)-50:])
#     signals1 -= np.mean(signals1[len(signals1)-50:])
#     signals2 -= np.mean(signals2[len(signals2)-50:])
#     signals3 -= np.mean(signals3[len(signals3)-50:])
#     signals4 -= np.mean(signals4[len(signals4)-50:])
# 
#     # import pdb; pdb.set_trace()
#     # signals0 = signals0[np.min(np.where(signals0 < 0)[0]):]
#     # signals1 = signals1[np.min(np.where(signals1 < 0)[0]):]
#     # signals2 = signals2[np.min(np.where(signals2 < 0)[0]):]
#     # signals3 = signals3[np.min(np.where(signals3 < 0)[0]):]
#     # signals4 = signals4[np.min(np.where(signals4 < 0)[0]):]
#     signals = np.concatenate((np.full(50, 0), signals0, signals1, signals2, signals3, signals4))
# 
#     return signals
# 
# f = "50Hz/170823_Scr_05_0.abf"
# scr = loadavg3(f)
# f = "50Hz/170814_Kd100_01_1.abf"
# kd100 = align(scr, loadavg3(f))
# f = "50Hz/170816_Kd200_02_0.abf"
# kd200 = align(scr, loadavg3(f))
# f = "50Hz/170807_Kd300_03_0.abf"
# kd300 = align(scr, loadavg3(f))
# f = "50Hz/170821_Kd400_04_0.abf"
# kd400 = align(scr, loadavg3(f))
# 
# f, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, sharex = True, sharey = True)
# ax1.plot(scr, linestyle = "-", linewidth = 0.33)
# ax2.plot(kd100, linestyle = "-", linewidth = 0.33)
# ax3.plot(kd200, linestyle = "-", linewidth = 0.33)
# ax4.plot(kd300, linestyle = "-", linewidth = 0.33)
# ax5.plot(kd400, linestyle = "-", linewidth = 0.33)
# ax6.plot([0, 250], [0, 0], 'k-', lw=2)
# ax6.plot([0, 0], [0, -500], 'k-', lw=2)
# 
# f.savefig("figures/50Hz_example.pdf")

def loadsignals(f):
    name = os.path.splitext(os.path.basename(f))[0]
    io = neo.get_io(f)
    block = io.read_block(lazy = False)
    signals = sum([s.analogsignals for s in block.segments], [])
    for s in signals:
        s.t_start = 0.0 * pq.second
    oo = os.path.join(os.path.dirname(f), '%s.regions' % name)
    regions = np.loadtxt(oo).tolist()

    oo = os.path.join(os.path.dirname(f), '%s.selection' % name)
    if os.path.exists(oo):
        selection = np.loadtxt(oo)
        signals = [signals[int(i)] for i in selection]
        
    a, d = get_protocol(f)
    a[np.abs(a - a[0]) < 0.5] = np.nan

    times = None
    for s in signals:
        times_ = s.times.rescale(pq.second)
        if times is None:
            times = times_
        else:
            assert(np.array_equal(times, times_))
            
    signals = np.asarray([np.ravel(s.as_quantity()) for s in signals]) * signals[0].units
    signals = signals - np.mean(signals[:, :50], axis = 1)[:, None]
    signals0 = signals
    signals0[:, np.logical_or(times < regions[0], times > regions[1])] = 0 * signals0.units
    # signals0 = signals[:, np.logical_and(times > regions[0], times < regions[1]+5)]
    signals0 = signals0 - np.mean(signals0[:, -100:], axis = 1)[:, None]
    
    n = 10000 // signals.shape[1]
    
    signals0 = signals0[:n, :]
    # import pdb; pdb.set_trace()
    signals0[signals0 > 0] = 0 * signals0.units

    signals = np.ravel(signals0)
    a = np.tile(a, n)
    
    # import pdb; pdb.set_trace()

    return signals, a

f = "20180726/0005_0009.abf"
scr, a = loadsignals(os.path.join("/Users/lw/Dropbox/Univ/masterThesis/data/Ephys/", f))

f = "20180726/0005_0011.abf"
kd100, a1 = loadsignals(os.path.join("/Users/lw/Dropbox/Univ/masterThesis/data/Ephys/", f))

f, (ax1, ax2, ax6) = plt.subplots(3, sharex = True, sharey = True)
ax1.plot(a + 100, linestyle = "-", linewidth = 0.33)
ax1.plot(scr, linestyle = "-", linewidth = 0.33)
ax2.plot(a1 + 100, linestyle = "-", linewidth = 0.33)
ax2.plot(kd100, linestyle = "-", linewidth = 0.33)
ax6.plot([0, 1000], [0, 0], 'k-', lw=2)
ax6.plot([0, 0], [0, -1000], 'k-', lw=2)

# # plt.show()
f.savefig("figure4a_.pdf")
