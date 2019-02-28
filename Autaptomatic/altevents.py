#!/usr/bin/env python3

import neo
import quantities as pq
import numpy as np
import scipy.signal as signal
import scipy.ndimage.filters as filters
from multiprocessing import Pool
import pandas as pd

import matplotlib.pyplot as plt

import os
import sys

from argparse import ArgumentParser

ap = ArgumentParser(description = "AUTAPTOMATIC")

ap.add_argument("input", metavar = "FILE", nargs = "+",
                    help = "input files")

args = ap.parse_args()

colNames = ["FILENAME", "PROTOCOL", "TIME", "CHARGE"]
rows = []

def process(f):
    try:
        # f = "/Users/lw/Dropbox/Univ/masterThesis/data/Ephys/20180803A/0059_0000.abf"
        print(f)
        io = neo.get_io(f)
        block = io.read_block(lazy = False)
        signals = sum([s.analogsignals for s in block.segments], [])

        for s in signals:
            s.t_start = 0.0 * pq.second

        sampling_rate = signals[0].sampling_rate

        # low pass filter
        sigma = 1000 * 0.53002 / 4 / 1e6 * float(sampling_rate)
        filtered_signals = [np.ravel(filters.gaussian_filter(s, sigma)) for s in signals]

        # skip one second
        skip = int(1 * sampling_rate)
        filtered_signals = [s[skip:] for s in filtered_signals]

        # smooth
        sigma = 1e5 * 0.53002 / 4 / 1e6 * float(sampling_rate)
        smoothed_signals = [np.ravel(filters.gaussian_filter(s, sigma)) for s in filtered_signals]

        # high pass filter
        filtered_signals = [s - t
            for s, t in zip(filtered_signals, smoothed_signals)]

        # threshold
        thres = 1.96 * np.mean([np.std(s) for s in filtered_signals])

        # delay embedding
        filtered_signals = [s[:, None ]for s in filtered_signals]
        y = np.concatenate(filtered_signals, axis = 1)

        y_dim = y.shape[0]
        emb_dim = int((10*pq.millisecond*sampling_rate).simplified)
        emb_offset = int((2*pq.millisecond*sampling_rate).simplified)
        m = int((y_dim - (emb_dim - 1)) / emb_offset)

        if m < 0:
            return []

        indices = np.repeat([np.arange(emb_dim)], m, axis = 0)
        indices += (np.arange(m).astype(np.int64) * emb_offset)[:, None]

        onset = indices[:, 0]
        onset = np.repeat(onset[:, None], y.shape[1], axis = 1)
        onset += np.arange(y.shape[1])[None, :] * int(signals[0].duration * sampling_rate)
        onset += (np.arange(y.shape[1])[None, :] + 1) * skip
        onset = np.ravel(onset.T)

        y_emb = np.asarray(y[indices])
        y_emb = np.transpose(y_emb, (2, 0, 1)).reshape((-1, emb_dim))

        filt = np.logical_and(np.abs(y_emb[:,  0]) < thres, 
                              np.abs(y_emb[:, -1]) < thres)

        y_emb = y_emb[filt, :]
        onset = onset[filt]

        x = (np.arange(emb_dim) / sampling_rate).simplified

        charge = -pq.trapz(y_emb * signals[0].units, x = x).rescale("pA*s") # picocoulomb
        
        filt = charge > -10000
        
        charge = charge[filt]
        onset = onset[filt]

        np.savetxt("charge.txt", charge)

        # plt.hist(np.asarray(charge), bins = 1000)
        # plt.show()
        
        onset = (onset / sampling_rate).rescale("s")
        
        rows = []
        for i in range(len(charge)):
            row = [f, os.path.basename(io._axon_info["sProtocolPath"].decode()),
                float(onset[i]), float(charge[i])]
            rows.append(row)
        return rows
    except ValueError:
        return []

rows = Pool().map(process, args.input)
rows = sum(rows, [])

dataFrame = pd.DataFrame({key:value for key, value in zip(colNames, zip(*rows))}, columns = colNames)
dataFrame.to_csv("events.csv.gz", compression = "gzip", index = False)
