#!/usr/bin/env python3

import neo
import quantities as pq
import numpy as np
import scipy.signal as signal
import scipy.ndimage.filters as filters
from multiprocessing import Pool
import pandas as pd

import os
import sys

from argparse import ArgumentParser

ap = ArgumentParser(description = "AUTAPTOMATIC")

ap.add_argument("input", metavar = "FILE", nargs = "+",
                    help = "input files")

args = ap.parse_args()

colNames = ["FILENAME", "PROTOCOL", "TIME", "AMPLITUDE", "CHARGE", "RISE_TIME"]
rows = []

def process(f):
    try:
        print(f)
        io = neo.get_io(f)
        block = io.read_block(lazy = False)
        signals = sum([s.analogsignals for s in block.segments], [])

        for s in signals:
            s.t_start = 0.0 * pq.second
        sampling_rate = signals[0].sampling_rate

        # filter
        sigma = 1000 * 0.53002 / 4 / 1e6 * float(sampling_rate)
        filtered_signals = [np.ravel(filters.gaussian_filter(s, sigma)) for s in signals]

        # skip one second
        filtered_signals = [s[int(3 * sampling_rate):] for s in filtered_signals]

        # define template
        baseline = 5.0 * pq.millisecond
        length = 10.0 * pq.millisecond

        amplitude = 2.0 * pq.picoampere / 2
        rise = 0.5 * pq.millisecond
        decay = 3.0 * pq.millisecond

        n_baseline = (baseline * sampling_rate).simplified
        n_samples = (length * sampling_rate).simplified
        t = np.linspace(0, float(length), float(n_samples)) * pq.millisecond

        template = (1 - np.exp(-t/rise)) * np.exp(-t/decay)
        template = np.concatenate([np.zeros((int(n_baseline),)), template])
        template = template * amplitude / np.max(template)
        template = np.asarray(template)

        # delay embedding
        filtered_signals = [s[:, None ]for s in filtered_signals]
        y = np.concatenate(filtered_signals, axis = 1)

        y_dim = y.shape[0]
        emb_dim = template.shape[0]
        m = y_dim - (emb_dim - 1)
        
        if m < 0:
            return []
        
        indices = np.repeat([np.arange(emb_dim)], m, axis = 0)
        indices += np.arange(m)[:, None]

        y_emb = np.asarray(y[indices])

        # event detection criterion
        s_y = np.sum(y_emb, axis = 1)
        s_y2 = np.sum(np.square(y_emb), axis = 1)
        s_ey = np.sum(y_emb * template[None, :, None], axis = 1)
        s_e = np.sum(template)
        s_e2 = np.sum(np.square(template))

        n = len(template)
        s = (s_ey - s_e * s_y / n) / (s_e2 - s_e * s_e / n)
        c = (s_y - s * s_e) / n

        # sse = (y_emb - s[:, None, :] * template[None, :, None])
        sse = s_y2 + np.square(s) * s_e2 + np.square(c) * n \
            - 2 * (s * s_ey + c * s_y - s * c * s_e)

        se = np.sqrt(sse / (n - 1))
        dc = s / se

        # event separation
        thres = -3
        events = [np.where(dc[:, i] < thres)[0] for i in range(dc.shape[1])]

        separation = 10.0 * pq.millisecond
        n_samples_separation = int((separation * sampling_rate).simplified)

        filtered_y_emb_onset = []
        filtered_y_emb = []
        for i in range(dc.shape[1]):
            ii = np.insert(np.diff(events[i]) > n_samples_separation, 0, True)
            if len(events[i]) > 0:
                ii = events[i][ii]
                indices = np.repeat([np.arange(emb_dim)], len(ii), axis = 0)
                indices += np.asarray(ii)[:, None]
                filtered_y_emb_onset.append(y.shape[0] * i + int(sampling_rate) * (i+1) + ii)
                filtered_y_emb.append(y[indices, i])

        filtered_y_emb_onset = np.concatenate(filtered_y_emb_onset)        
        filtered_y_emb = np.vstack(filtered_y_emb)        

        s_y = np.sum(filtered_y_emb, axis = 1)
        s_y2 = np.sum(np.square(filtered_y_emb), axis = 1)
        s_ey = np.sum(filtered_y_emb * template[None, :], axis = 1)

        s = (s_ey - s_e * s_y / n) / (s_e2 - s_e * s_e / n)
        c = (s_y - s * s_e) / n

        filtered_y_emb -= c[:, None]

        # event filtering
        base = np.mean(filtered_y_emb[:, :int(n_baseline/2)], axis = 1)

        filt = np.ones((3,)) / 3
        mov_avg = np.apply_along_axis(lambda m: np.convolve(m, filt, mode="same"), axis = 1, arr = filtered_y_emb)

        peak = np.min(mov_avg, axis = 1)
        peak_ind = np.argmin(mov_avg, axis = 1)

        lo = 0.2 * (peak - base)
        hi = 0.8 * (peak - base)

        m = filtered_y_emb.shape[0]

        rise_time = np.zeros((m,)) * pq.second
        for i in range(m):
            v = np.abs(filtered_y_emb[i, :peak_ind[i]] - base[i])
            indices = np.arange(peak_ind[i])
            indices[v >= np.abs(lo[i])] = 0
            inner_lo_ind = np.max(indices)
            indices = np.arange(peak_ind[i])
            indices[v <= np.abs(hi[i])] = 1e6
            inner_hi_ind = np.min(indices)
            rise_time[i] = ((inner_hi_ind-inner_lo_ind) / sampling_rate).simplified

        mid = 0.5 * (peak - base)
        half_width = np.zeros((m,)) * pq.second
        for i in range(m):
            v = np.abs(filtered_y_emb[i, :peak_ind[i]] - base[i])
            indices = np.arange(peak_ind[i])
            indices[v >= np.abs(mid[i])] = 0
            left_ind = np.max(indices)
            v = np.abs(filtered_y_emb[i, peak_ind[i]:] - base[i])
            indices = np.arange(n-peak_ind[i])
            indices[v >= np.abs(mid[i])] = 1e6
            right_ind = np.min(indices) + peak_ind[i]
            half_width[i] = ((right_ind-left_ind) / sampling_rate).simplified

        filt = np.logical_and(\
            np.logical_and(rise_time < 1.5 * pq.millisecond, \
                rise_time > 0.15 * pq.millisecond), \
            np.logical_and(half_width < 5 * pq.millisecond, \
                half_width > 0.5 * pq.millisecond))

        print("Detected %i events" % np.sum(filt))

        filtered_y_emb = filtered_y_emb[filt, :]
        filtered_y_emb_onset = filtered_y_emb_onset[filt]

        base = np.mean(filtered_y_emb[:, :int(n_baseline/2)], axis = 1)
        filtered_y_emb -= base[:, None]
        filtered_y_emb = filtered_y_emb * pq.picoampere

        t = np.linspace(0, float(length+baseline), num = int(n_samples+n_baseline)) * pq.millisecond

        onset = (filtered_y_emb_onset / sampling_rate).rescale("s")

        # keep only event for measurement, zero everything else
        trunc_y_emb = filtered_y_emb
        trunc_y_emb[:, :int(n_baseline)] = 0 * pq.picoampere
        for i in range(len(onset)):
            trunc_y_emb[i, int(n_baseline+2*half_width[i]*sampling_rate):] = 0 * pq.picoampere

        amp = -np.min(trunc_y_emb, axis = 1).rescale("pA")
        charge = -pq.trapz(y = trunc_y_emb, x = t, axis = 1).rescale("pA*s")
        rise_time = rise_time[filt].rescale("ms")

        rows = []
        for i in range(len(onset)):
            row = [f, os.path.basename(io._axon_info["sProtocolPath"].decode()),
                float(onset[i]), float(amp[i]), float(charge[i]), float(rise_time[i])]
            rows.append(row)
        return rows
    except ValueError:
        return []

rows = Pool().map(process, args.input)
rows = sum(rows, [])

dataFrame = pd.DataFrame({key:value for key, value in zip(colNames, zip(*rows))}, columns = colNames)
dataFrame.to_csv("events.csv", index = False)
