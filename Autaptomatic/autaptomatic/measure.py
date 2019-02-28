

import quantities as pq
import numpy as np

from scipy.signal import find_peaks

def get_measure(measure, recording, mask, y):
    value = np.nan
    if measure == 'AMPLITUDE':
        value = -np.min(y)
    elif measure == 'CHARGE':
        value = -pq.trapz(y, x = recording.x[mask])
    elif measure == 'RISE_TIME':     
        filt = np.ones((3,)) / 3
        mov_avg = np.convolve(y, filt, mode="same")
        # 
        base = 0
        peak = np.min(mov_avg)
        peak_ind = np.argmin(mov_avg)
        
        value = np.nan * pq.millisecond
        try:
            # 
            lo = 0.2 * (peak - base)
            hi = 0.8 * (peak - base)
            # 
            v = y[:peak_ind]
            indices = np.arange(peak_ind)
            indices[v <= lo] = 0
            inner_lo_ind = np.max(indices)
            indices = np.arange(peak_ind)
            indices[v >= hi] = 1e6
            inner_hi_ind = np.min(indices)
            
            value = ((inner_hi_ind-inner_lo_ind) * (recording.x[1] - recording.x[0])).rescale("ms")
            # import pdb; pdb.set_trace()
        except ValueError:
            pass
    elif measure == 'N_PEAKS':       
        peaks, _ = find_peaks(-y, width = 10)
        value = peaks.size
        if value < 1:
            value = 1
        value = value * pq.dimensionless
        
        
    return value