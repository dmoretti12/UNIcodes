# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 16:50:29 2022

@author: Davide
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


plt.rcParams['figure.figsize'] = [10,6]
plt.rcParams.update({'font.size': 18})
plt.style.use('seaborn')
stZ = pd.read_csv('C1XArapuca_CRT_250MHz_00000.txt', sep=" ", header=5, index_col = None, encoding = 'unicode_escape')

#stZ = read(filenameZ)
streams = [stZ.copy()]
traces = []
for st in streams:
    tr = st[0].trim(otime, otime+120)
    traces.append(tr)
    
delta = stZ[0].stats.delta
starttime = np.datetime64(stZ[0].stats.starttime)
endtime = np.datetime64(stZ[0].stats.endtime)
signalZ = traces[0].data/10**6
minsignal, maxsignal = signalZ.min(), signalZ.max()

t = traces[0].times("utcdatetime") 

## Compute Fourier Transform
n = len(t)
fhat = np.fft.fft(signalZ, n) #computes the fft
psd = fhat * np.conj(fhat)/n
freq = (1/(delta*n)) * np.arange(n) #frequency array
idxs_half = np.arange(1, np.floor(n/2), dtype=np.int32) #first half index
psd_real = np.abs(psd[idxs_half]) #amplitude for first half


## Filter out noise
sort_psd = np.sort(psd_real)[::-1]
# print(len(sort_psd))
threshold = sort_psd[300]
psd_idxs = psd > threshold #array of 0 and 1
psd_clean = psd * psd_idxs #zero out all the unnecessary powers
fhat_clean = psd_idxs * fhat #used to retrieve the signal

signal_filtered = np.fft.ifft(fhat_clean) #inverse fourier transform


## Visualization
fig, ax = plt.subplots(4,1)
ax[0].plot(t, signalZ, color='b', lw=0.5, label='Noisy Signal')
ax[0].set_xlabel('t axis')
ax[0].set_ylabel('Accn in Gal')
ax[0].legend()

ax[1].plot(freq[idxs_half], np.abs(psd[idxs_half]), color='b', lw=0.5, label='PSD noisy')
ax[1].set_xlabel('Frequencies in Hz')
ax[1].set_ylabel('Amplitude')
ax[1].legend()

ax[2].plot(freq[idxs_half], np.abs(psd_clean[idxs_half]), color='r', lw=1, label='PSD clean')
ax[2].set_xlabel('Frequencies in Hz')
ax[2].set_ylabel('Amplitude')
ax[2].legend()

ax[3].plot(t, signal_filtered, color='r', lw=1, label='Clean Signal Retrieved')
ax[3].set_ylim([minsignal, maxsignal])
ax[3].set_xlabel('t axis')
ax[3].set_ylabel('Accn in Gal')
ax[3].legend()

plt.subplots_adjust(hspace=0.6)
plt.savefig('real-signal-analysis.png', bbox_inches='tight', dpi=300)