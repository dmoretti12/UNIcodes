# -*- coding: utf-8 -*-
"""
Created on Sat May 21 10:00:56 2022

@author: Davide
"""

from functions import read_binary, norm_fit, ampl_pe, tot
import numpy as np
import matplotlib.pyplot as plt
#import math

plt.close('all')

file="wave0.dat"
header=12   #in 32 bit, header=6 in 16 bit
vpp_adc = 2   #V_peak-to-peak of adc
adcbit=14
adc_res = vpp_adc/2**(adcbit)

a, t, a62, t62 = read_binary(file, header, adc_res)

# plt.figure()
# for i in range(0, 100):
#     plt.hist2d(t[i], a[i], bins=(500, 500), range=None, density=False, weights=None, cmin=None, cmax=None)
# plt.show()

#%%
mu, std = norm_fit(a, a, t, inizio_sig=3.5)

a0 = np.empty((len(a), len(a[0])))

for i in range (0, len(a)):
    
    a0[i] = a[i] - mu[i]
    
asa = np.zeros(0)   # array che contiene indice wvf con media molto diversa dalle altre
for i in range(0, len(a)):
    if mu[i] > 0.214 or mu[i] < 0.207:
        asa = np.append(asa, i)
        # print(i, mu[i])
    
a0b = np.delete(a, asa, axis=0)
tb = np.delete(t, asa, axis=0)
#%% 
t_ns = tb * 1000

carica, picco = ampl_pe(a0b, a0b, t_ns, t_min=int(4*1000/4), t_max=int(4.5*1000/4))
#c_tot, ss_tot, p_tot, _ = tot(a0, a0, t, t_min=int(4*1000/4), t_max=int(4.6*1000/4), threshold=0.02)
# deltaT = np.zeros(len(a))
# for i in range (0, len(a)):
#     deltaT[i] = (ss_tot[i,1] - ss_tot[i,0])*4e-03
    
plt.figure()
plt.title('Amplitude vs Charge 250 MHz')
plt.ylabel('Amplitude [V]')
plt.xlabel('Charge [V.ns]')
plt.plot(carica, picco, '.')
plt.grid(True)
#plt.savefig('amp_ch62_C1X', dpi=300)
plt.show()

# rifai la baseline tagliando valori più bassi  e più alti di tot in ampiezza

#%%
#HISTOGRAM 2D (color for the z) OF ALL THE WAVEFORMS (PERCISTENCY HISTOGRAM) (2500 bin per x, 1000 per y)
nbin=1000
xbi = np.arange(0, 10, 10/(2.5*nbin))
ybi = np.arange(0.195, 0.33, (0.33-0.195)/nbin)

plt.figure()
# for i in range(0, 100):
h, xb, yb = np.histogram2d(t.flatten(), a.flatten(), bins=(xbi, ybi), range=[[0, 10], [0.195, 0.33]], normed=None, weights=None, density=None)
h = np.log(h.T)
#X, Y = np.meshgrid(xb, yb)
plt.pcolormesh(xb, yb, h, cmap='jet')
plt.xlabel('Time [us]')
plt.ylabel('Amplitude [V]')
plt.colorbar()
plt.show()



