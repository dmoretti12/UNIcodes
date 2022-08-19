# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 17:45:28 2022

@author: Davide
"""
from functions import read_binary, norm_fit, ampl_pe, isto_carica, denoising, mov_av, signal
import numpy as np
import matplotlib.pyplot as plt
#import math

plt.close('all')

file="0_wave0_34V40_20ADC_external_4led1.dat"
header=12   #in 32 bit, header=6 in 16 bit
vpp_adc = 2   #V_peak-to-peak of adc
adcbit=12
adc_res = vpp_adc/2**(adcbit)

a, t, a62, t62 = read_binary(file, header, adc_res)

#%%  BASELINE

mu, std = norm_fit(a, a, t, inizio_sig=3.5)

a0 = np.empty((len(a), len(a[0])))

for i in range (0, len(a)):
    
    a0[i] = a[i] - mu[i]
    
# asa = np.zeros(0)   # array che contiene indice wvf con media molto diversa dalle altre
# for i in range(0, len(a)):
#     if mu[i] > 0.214 or mu[i] < 0.207:
#         asa = np.append(asa, i)
#         # print(i, mu[i])
    
# a0b = np.delete(a, asa, axis=0)
# tb = np.delete(t, asa, axis=0)

#%%

t_ns = t * 1000

carica, picco = ampl_pe(a0, a0, t_ns, t_min=int(7.97*1000/4), t_max=int(8.3*1000/4))

plt.figure()
plt.title('Amplitude vs Charge 250 MHz')
plt.ylabel('Amplitude [V]')
plt.xlabel('Charge [V.ns]')
plt.plot(carica, picco, '.')
plt.grid(True)
#plt.savefig('run0', dpi=300)
plt.show()


#%%
#HISTOGRAM 2D (color for the z) OF ALL THE WAVEFORMS (PERCISTENCY HISTOGRAM) 
nbin=1000
xbi = np.arange(0, 10, 10/(2.5*nbin))
ybi = np.arange(0.195, 0.33, (0.33-0.195)/(nbin/4))

plt.figure()
h, xb, yb = np.histogram2d(t.flatten(), a.flatten(), bins=(xbi, ybi), range=[[0, 10], [0.195, 0.33]], normed=None, weights=None, density=None)
h = np.log(h.T)
#X, Y = np.meshgrid(xb, yb)
plt.pcolormesh(xb, yb, h, cmap='jet')
plt.xlabel('Time [us]')
plt.ylabel('Amplitude [V]')
plt.colorbar()
plt.show()

#%% CHARGE HISTO 16 dec miniarapuca 250 MHz
from functions import multi_fit

t_min = int((7.96*1000)/4)   #4 ns per point
t_max = int((8.6*1000)/4)
mm = np.array(0)
# for j in range (0, 200):
    

_bin=150

n, b, charge = isto_carica(a, a0, t_ns, t_min, t_max, _bin)

plt.figure()
plt.title('Charge histogram (run0) 250 MHz')
plt.ylabel('Counts')
plt.xlabel('Charge [ns.V]')
plt.hist(charge, b, density=False)
#plt.savefig('ch_histo250', dpi=300)
plt.show()

x = np.linspace(-0.35, 1.9, 1000)
# p0 = [2, -0.025, 0.08, 2, 0.14, 0.10, 1.2, 0.3]
p0 = [500, -0.025, 0.08, 500, 0.14, 0.10, 250, 0.3] # density=False in both functions

#p0=[1.8, -0.05, 0.08, 1.6, 0.13, 0.07, 1, 0.3]   #parametri molto giusti

numg = 3
for i in range(0, numg):
    p0 = np.append(p0, p0[6]/(3**(i+1)))
    
p0 = np.array(p0)

pop, cov, SNR, chi2 = multi_fit(numg, p0, n, b, charge, x)

print(pop[7]-pop[4], SNR, "chi2 = ", chi2)#, j)



#%% controllo per lab01-07

qff = denoising(a, a)
qf = mov_av(a, qff)

tre=np.zeros(len(a))
tre.fill(0.3)
carica, ss, peak, peakpos = signal(a, a, qf, t, tre, 20, 40)

ssa = np.array(0)
peaka = np.array(0)
ca = np.array(0)

for i in range(0, len(a)):
    ssa = np.append(ssa, ss[i])
    peaka = np.append(peaka, peak[i])
    ca = np.append(ca, carica[i])
ssa = np.delete(ssa, 0)
peaka = np.delete(peaka, 0)
ca = np.delete(ca, 0)


plt.figure()
plt.title('Charge histogram 250 MHz')
plt.ylabel('Counts')
plt.xlabel('Charge [ns.V]')
plt.hist(ca, 100, density=False)
#plt.savefig('ch_histo250', dpi=300)
plt.show()






