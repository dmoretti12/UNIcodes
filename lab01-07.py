# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 17:05:21 2022

@author: Davide
"""

from functions import read_long_txt_time, denoising, norm_fit, signal, isto_carica, tot, ampl_pe, mov_av, multi_fit, gauss, multi_gauss, time_res
import numpy as np
import matplotlib.pyplot as plt
import math, os
from scipy.optimize import curve_fit
from scipy.stats import norm

#%%
plt.close("all") #chiude tutte le finestre

def expo(x, a, b):
    return np.exp(a*x + b)

def retta(x, m, q):
    return m*x + q

#  DATA lab SiPM 01-07-22

# os.chdir(r"H:\Il mio Drive\WORK\20220701_darkcount_check")
file = 'C1TraceC00002.txt'   #apri file txt e guarda linee
righe_iniz = 5
puntixwvf = 125001

a, t, a62, t62 = read_long_txt_time(file, righe_iniz, puntixwvf)

#%% manipulating array

q=a.flatten()
q=np.delete(q, len(q)-1)
q = np.reshape(q, (100, 1250))
tt=t.flatten()
tt=np.delete(tt, len(t)-1)
tt = np.reshape(tt, (100, 1250))

mu = np.empty(len(q))
std = np.empty(len(q))
tre = np.empty(len(q))
for i in range(0, len(q)):

    mu[i], std[i] = norm.fit(q[i])
    q[i] = q[i] - mu[i]
    # tre[i] = 0.3 #mu[i] + 3*std[i]


#%%denoising
from functions import signal
plt.close('all')
qff = denoising(q, q, 40)
qf = mov_av(q, qff, 6)
tre.fill(0.0115)


carica, ss, peak, peakpos, sop, sot, e = signal(q, q, qff, tt, tre, 15, 80)

# for i in range(13, 24):
#     plt.figure()
#     plt.plot(tt[i], q[i], '-')
#     plt.plot(tt[i], qff[i], '-')
#     plt.show()
    
    
#%%
# ssa = np.array(0)
# peaka = np.array(0)
# ca = np.array(0)

# for i in range(0, len(q)):
#     ssa = np.append(ssa, ss[i])
#     peaka = np.append(peaka, peak[i])
#     ca = np.append(ca, carica[i])
# ssa = np.delete(ssa, 0)
# peaka = np.delete(peaka, 0)
# ca = np.delete(ca, 0)

plt.figure()
plt.title('Charge histogram 250 MHz')
plt.ylabel('Counts')
plt.xlabel('Charge [ns.V]')
plt.hist(carica.flatten(), np.linspace(-0.25e-8, 0.75e-8, 100), density=False)
#plt.savefig('ch_histo250', dpi=300)
plt.show()

