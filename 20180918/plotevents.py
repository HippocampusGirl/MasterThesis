#!/usr/bin/env python3

import neo
import quantities as pq
import numpy as np
import scipy.signal as signal
import scipy.ndimage.filters as filters

import matplotlib.pyplot as plt

import os
import sys

def measure(f):
    io = neo.get_io(os.path.join("/Users/lw/Dropbox/Univ/masterThesis/data/Ephys/", f))
    block = io.read_block(lazy = False)
    signals = sum([s.analogsignals for s in block.segments], [])

    for s in signals:
        s.t_start = 0.0 * pq.second
    sampling_rate = signals[0].sampling_rate
    
    signals = [s[int(sampling_rate*3):] for s in signals]

    # filter
    sigma = 1000 * 0.53002 / 4 / 1e6 * float(sampling_rate)
    filtered_signals = [np.ravel(filters.gaussian_filter(s, sigma)) for s in signals]

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
    filtered_signals = [s[:, None ] for s in filtered_signals]
    y = np.concatenate(filtered_signals, axis = 1)

    y_dim = y.shape[0]
    emb_dim = template.shape[0]
    m = y_dim - (emb_dim - 1)
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
    
    amplitude = (y_emb - y_emb[:, :int(n_baseline), :].mean(axis = 1)[:, None, :]).min(axis = 1)
    
    # import pdb; pdb.set_trace()
    
    # dc *= -1

    # event separation
    thres = -3
    minamp = -10
    events = [np.where( \
        np.logical_and(dc[:, i] < thres, amplitude[:, i] < minamp))[0] for i in range(dc.shape[1])]

    separation = 3.0 * pq.millisecond
    n_samples_separation = int((separation * sampling_rate).simplified)
    
    filtered_y_emb = None
    for i in range(dc.shape[1]):
        if len(events[i]) > 0:
            ii = np.insert(np.diff(events[i]) > n_samples_separation, 0, True)
            events[i] = events[i][ii]

    # signals = [s[int(sampling_rate):] for s in signals]

    return filtered_signals, events

slen = 40000
wlen = 4000

def getwindow(e, ts, target):
    jmin = None, None
    vmin = 5000
    kmin = np.inf
    nntarget = 0
    for i in range(0, len(e), 1):
        for j in range(0, slen-wlen, 10):
            nn = np.sum(np.logical_and(e[i] > j, e[i] < j + wlen))
            k = np.square(float(nn) - target)
            if k < kmin and nn > 0:
                nntarget = nn
                kmin = k
    
    for i in range(0, len(e), 1):
        for j in range(0, slen-wlen, 10):
            scr = ts[i][j:j+wlen]
            scr = scr - np.mean(scr)
            nn = np.sum(np.logical_and(e[i] > j, e[i] < j + wlen))
            # if nn == nntarget:
            if nn == nntarget and np.max(scr) < vmin:
                vmin = np.max(scr)
                print(nn, vmin)
                jmin = i, j
    return jmin

# f = "20180726/0006_0000.abf"
f = "20180727/0028_0000.abf"
scr, scre = measure(f)
# import pdb; pdb.set_trace()
i, j = getwindow(scre, scr, 4)
scr = scr[i][j:j+wlen]
scr = scr - np.mean(scr)

print()
# f = "20180726/0006_0001.abf"
f = "20180910/0136_0001.abf"
kd100, kd100e = measure(f)
i, j = getwindow(kd100e, kd100, 4)
kd100 = kd100[i][j:j+wlen+250]
kd100 = kd100 - np.mean(kd100)

print()
f = "20180924/0159_0003.abf"
kd200, kd200e = measure(f)
i, j = getwindow(kd200e, kd200, 4)
kd200 = kd200[i][j:j+wlen+250]
kd200 = kd200 - np.mean(kd200)

print()
f = "20180910/0136_0005.abf"
kd300, kd300e = measure(f)
i, j = getwindow(kd300e, kd300, 4)
kd300 = kd300[i][j:j+wlen+250]
kd300 = kd300 - np.mean(kd300)

print()
f = "20180910/0136_0010.abf"
kd400, kd400e = measure(f)
i, j = getwindow(kd400e, kd400, 4)
kd400 = kd400[i][(j-250):j+wlen]
kd400 = kd400 - np.mean(kd400)

print()
f = "20180830/0116_0007.abf"
kd500, kd500e = measure(f)
i, j = getwindow(kd500e, kd500, 1)
kd500 = kd500[i][(j-250):j+wlen]
kd500 = kd500 - np.mean(kd500)

# plt.plot(scr)
# plt.plot(kd100)
# plt.plot(kd200)
# plt.plot(kd300)
# plt.plot(kd400)
# plt.show()

f, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(7, sharex = True, sharey = True)
ax1.plot(scr[::1], linestyle = "-", linewidth = 0.33)
ax2.plot(kd100[::1], linestyle = "-", linewidth = 0.33)
ax3.plot(kd200[::1], linestyle = "-", linewidth = 0.33)
ax4.plot(kd300[::1], linestyle = "-", linewidth = 0.33)
ax5.plot(kd400[::1], linestyle = "-", linewidth = 0.33)
ax6.plot(kd500[::1], linestyle = "-", linewidth = 0.33)
ax7.plot([0, 250./1.], [0, 0], 'k-', lw=2)
ax7.plot([0, 0], [0, -50], 'k-', lw=2)

f.savefig("figure5a_.pdf")
