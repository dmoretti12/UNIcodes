# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 16:00:45 2022

@author: Davide
"""

from functions import read_binary, norm_fit, ampl_pe, isto_carica, denoising, mov_av, signal
import numpy as np
import matplotlib.pyplot as plt
#import math
import os

plt.close('all')
os.chdir(r"C:\Users\Davide\Desktop\20220726_board_at_warm_test\run0_34V4_external_pd2")
file="0_wave0_34V4_external_pd2.dat"
header=12   #in 32 bit, header=6 in 16 bit
vpp_adc = 2   #V_peak-to-peak of adc
adcbit=12
adc_res = vpp_adc/2**(adcbit)

a, t, a62, t62 = read_binary(file, header, adc_res)

os.chdir(r"H:\Il mio Drive\WORK")
#%%  BASELINE

mu, std = norm_fit(a, a, t, inizio_sig=3.5)

a0 = np.empty((len(a), len(a[0])))

for i in range (0, len(a)):
    
    a0[i] = a[i] - mu[i]


t_ns = t * 1000

qff = denoising(a, a0, 50)
qf = mov_av(a, a0, 3)


#%% HISTOGRAM 2D (color for the z) OF ALL THE WAVEFORMS (PERCISTENCY HISTOGRAM) (2500 bin per x, 1000 per y)
nbin=225
xbi = np.arange(0, 10, 10/(11.1*nbin))
ybi = np.arange(0.43, 0.54, (0.54-0.43)/nbin)

plt.figure()
# for i in range(0, 100):
h, xb, yb = np.histogram2d(t.flatten(), qff.flatten(), bins=(xbi, ybi), range=[[0, 10], [0.43, 0.54]], normed=None, weights=None, density=None)
h = np.log(h.T)
#X, Y = np.meshgrid(xb, yb)
plt.pcolormesh(xb, yb, h, cmap='jet')
plt.xlabel('Time [us]')
plt.ylabel('Amplitude [V]')
plt.colorbar()
plt.show()

#%%

from functions import multi_fit
    
t_min = int((5.028*1000)/4)   #4 ns per point
t_max = int((5.324*1000)/4)
    

_bin=133 #115 #80+j #143  

n, b, charge = isto_carica(a0, qff, t_ns, t_min, t_max, _bin)

plt.figure()
plt.title('Charge histogram 26-07 data 250 MHz')
plt.ylabel('Counts')
plt.xlabel('Charge [ns.V]')
plt.hist(charge, b, density=False)
#plt.savefig('ch_histo250', dpi=300)
plt.show()
#%%
x = np.linspace(-2.5, 12.5, 1000)

p0 = [400, -0.025, 0.5, 400, 1.7, 0.5, 250, 3.4] # density=False in both functions

numg = 3
for i in range(0, numg):
    p0 = np.append(p0, p0[6]/(3**(i+1)))
    
p0 = np.array(p0)

pop, cov, SNR, chi2 = multi_fit(numg, p0, n, b, charge, x)

print('mu1-mu0 =', round(pop[7]-pop[4], 3), ',  SNR = ', round(SNR, 2), ",  chi2 =", round(chi2, 2))#, j)


#%% optimizing integration time
from functions import multi_fit  
#   TOGLI PLOT DALLA FUNZIONE MULTI_FIT

start1 = 5.100    #divisibili per 4???
start2 = 5.252
fin1 = 4.852
fin2 = 5.652

    
_bin=133
    
x = np.linspace(-2.5, 12.5, 1000)

p0 = [400, -0.025, 0.5, 400, 1.7, 0.5, 250, 3.4] # density=False in both functions

numg = 3
for i in range(0, numg):
    p0 = np.append(p0, p0[6]/(3**(i+1)))
    
p0 = np.array(p0)


sn = np.array(0)
mu1 = np.array(0)
mu1dif = np.array(0)
time = np.array(0)

while True:

    t_min = int((start1*1000)/4)   #4 ns per point
    t_max = int((start2*1000)/4)
        
    n, b, charge = isto_carica(a0, qff, t_ns, t_min, t_max, _bin)

    pop, cov, SNR, chi2 = multi_fit(numg, p0, n, b, charge, x)
    
    time = np.append(time, (t_max-t_min)*4)   #indici --> ns
    mu1 = np.append(mu1, pop[4])
    mu1dif = np.append(mu1dif, pop[7]-pop[4])
    sn = np.append(sn, SNR)
    # if SNR>3.026: print(start1, start2)
    
    if start1-0.004 > fin1:
        start1 = start1 - 0.004
    if start2+0.004 < fin2:
        start2 = start2 + 0.004
    else: break 

time = np.delete(time, 0)  #indici --> ns
mu1 = np.delete(mu1, 0)
mu1dif = np.delete(mu1dif, 0)
sn = np.delete(sn, 0)

plt.figure()
plt.plot(time, sn, '-')
plt.title('SNR')
plt.ylabel('SNR')
plt.xlabel('Integration time [ns]')
plt.show()

plt.figure()
plt.plot(time, mu1, '-')
plt.title('mu1')
plt.ylabel('mu1 [ns.V]')
plt.xlabel('Integration time [ns]')
plt.show()

plt.figure()
plt.plot(time, mu1dif, '-')
plt.title('mu2-mu1')
plt.ylabel('mu2-mu1 [ns.V]')
plt.xlabel('Integration time [ns]')
plt.show()

#%% optimizing binning
#   TOGLI PLOT DALLA FUNZIONE MULTI_FIT

from functions import multi_fit
    
t_min = int((5.028*1000)/4)   #4 ns per point
t_max = int((5.324*1000)/4)


x = np.linspace(-2.5, 12.5, 1000)

p0 = [400, -0.025, 0.5, 400, 1.7, 0.5, 250, 3.4] # density=False in both functions

numg = 3
for i in range(0, numg):
    p0 = np.append(p0, p0[6]/(3**(i+1)))
    
p0 = np.array(p0)

sn = np.array(0)
mu1 = np.array(0)
mu1dif = np.array(0)
bn = np.array(0)

for i in range(70, 300):
    
    n, b, charge = isto_carica(a0, a0, t_ns, t_min, t_max, i)

    pop, cov, SNR, chi2 = multi_fit(numg, p0, n, b, charge, x)
    
    bn = np.append(bn, i)   #indici --> ns
    mu1 = np.append(mu1, pop[4])
    mu1dif = np.append(mu1dif, pop[7]-pop[4])
    sn = np.append(sn, SNR)


bn = np.delete(bn, 0)  #indici --> ns
mu1 = np.delete(mu1, 0)
mu1dif = np.delete(mu1dif, 0)
sn = np.delete(sn, 0)

plt.figure()
plt.plot(bn, sn, '-')
plt.title('SNR')
plt.ylabel('SNR')
plt.xlabel('#bin')
plt.show()

plt.figure()
plt.plot(bn, mu1, '-')
plt.title('mu1')
plt.ylabel('mu1 [ns.V]')
plt.xlabel('#bin')
plt.show()

plt.figure()
plt.plot(bn, mu1dif, '-')
plt.title('mu2-mu1')
plt.ylabel('mu2-mu1 [ns.V]')
plt.xlabel('#bin')
plt.show()





