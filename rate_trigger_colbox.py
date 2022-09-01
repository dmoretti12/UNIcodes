# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 17:02:54 2022

@author: Davide
"""

from functions import zm, read_binary, read_long_txt, norm_fit, ampl_pe, isto_carica, denoising, mov_av, signal, center_bins, baseline, plot
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})
#import math
import os
from scipy.stats import norm
from scipy.optimize import curve_fit



plt.close('all')
#%% DATI AGOSTO 22 COLDBOX


os.chdir("H:\Il mio Drive\WORK\AUG22COLDBOX")
file="0_wave1_light_leakage_test_led_flange_covered.dat"   # 1139.5 Hz/mm^2 ev 410.2 k  - tre=37 15,100, filt = 7 17)
# file="0_wave1_light_leakage_test_led_fibers_connected.dat"   #  1099.0 Hz/mm^2 ev 395.6 k - tre=37 15,100, filt = 7 17)
# file="0_wave1_light_leakage_test_led_fibers_connected_dcem_digital_off.dat"  #  910.6 Hz/mm^2 ev 327.8 k - tre 37 15,100, filt= 7 17)


header=12
a, t, _, _, _ = read_binary(file, header, adc_res=1, ns=2) #ns tra i punti, 250Mhz=4ns, 500Mhz=2ns
os.chdir("H:\Il mio Drive\WORK")

#%% FILTER

q = mov_av(a, a, 7)


#%%  BASELINE
from functions import baseline

qf, mu = baseline(q, threshold=35, bins=100)

for i in range(0, len(a)):
    a[i] = a[i] - mu[i]
    
qff = mov_av(a, a, 17)

#%%  SIGNALS
from functions import signal
tre = 37
prima = 15
dopo = 100
time_after_pulse = int(345/2)
lim_amp = 450
debug=True
num_plot=2
num_ev_expected=600

carica, e = signal(a, qf, qff, t, tre, prima, dopo, time_after_pulse, lim_amp, debug, num_plot, num_ev_expected, 0, 0)

#  %time o  %%time per tempo esecuzione riga o cella direttamente sulla console (iPyhton)

#%% CHARGE HIST

charge = np.zeros(0)
for i in range(0, len(a)):
    c=np.zeros(len(np.where(carica[i]!=0)))
    c = carica[i, np.where(carica[i]!=0)]
    charge = np.append(charge, c)
        
    
plt.figure(figsize=[10, 6])
ax=plt.gca()
ax.tick_params(bottom=True, top=True, left=True, right=True)
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
plt.ylabel('Counts')
plt.xlabel('Charge [ns.ADC]')
plt.hist(charge, np.arange(-3000, 45000, 250), density=False)
plt.savefig('coldbox_aug_rand_trig_darkcounsearch.png', dpi=600)
# plt.show()

# TRIGGER RATE
rate = sum(e)/( len(a)*len(a[0])*2e-9*(20*36))  #Hz
print(round(rate,1), 'Hz/mm^2', 'ev', round(sum(e)*1e-3, 1),'k')





















#%% DATI AUG22 SPE

os.chdir("H:\Il mio Drive\WORK\AUG22COLDBOX")
header=12
file1="0_wave1_led_7V_50ns_digital_on.dat"  # 633.2 Hz/mm^2  ev 45.6 k     carica = 5101.268 ,  SNR =  4.17 ,  chi2 = 1.58 
file2='1_wave1_led_7V_50ns_digital_on.dat' #  647.7 Hz/mm^2  ev 46.6 k    carica = 5184.381 ,  SNR =  4.13 ,  chi2 = 1.3
file3='2_wave1_led_7V_50ns_digital_on.dat' #  647.0 Hz/mm^2   ev 46.6 k     carica = 5042.97 ,  SNR =  4.13 ,  chi2 = 1.64
a1, _, _, _, _ = read_binary(file1, header, adc_res=1, ns=4) #ns tra i punti, 250Mhz=4ns, 500Mhz=2ns
a2, _, _, _, _ = read_binary(file2, header, adc_res=1, ns=4) #ns tra i punti, 250Mhz=4ns, 500Mhz=2ns
a3, _, _, _, _ = read_binary(file3, header, adc_res=1, ns=4) #ns tra i punti, 250Mhz=4ns, 500Mhz=2ns
os.chdir("H:\Il mio Drive\WORK")

a = np.zeros((len(a1)*3, len(a1[0])))
np.concatenate((a1, a2, a3), axis=0, out=a)  #carica = 5109.719 ,  SNR =  4.14 , chi2=3.89
#carica = 5177.131 ,  SNR =  4.57 ,  chi2 = 0.98  con double FIT
# 642.6 Hz/mm^2
# ev 138.8 k
# t = np.arange(0, 4*len(a1[0]), 4).reshape(len(a),len(a[0]))
t= np.zeros((len(a),len(a[0])))
for i in range(0,len(t)):
    t[i] = np.arange(0, 4*len(a1[0]), 4)
del a1,a2,a3



#%% FILTER
from functions import mov_av

q = mov_av(a, a, 7)


#%% BASELINE                       
from functions import baseline

qf, mu = baseline(q, threshold=30, bins=100)

for i in range(0, len(a)):
    a[i] = a[i] - mu[i]

qff = mov_av(a, a, 17)

#%%
#HISTOGRAM 2D (color for the z) OF ALL THE WAVEFORMS (PERCISTENCY HISTOGRAM) (2500 bin per x, 1000 per y)

xbi = np.arange(4000, 7000, 4)
ybi = np.arange(-100, 450, 3)

plt.figure()
h, xb, yb = np.histogram2d(t.flatten(), qf.flatten(), bins=(xbi, ybi), normed=None, weights=None, density=None)
h = np.log(h.T)
#X, Y = np.meshgrid(xb, yb)
plt.pcolormesh(xb, yb, h, cmap='jet')
plt.xlabel('Time [ns]')
plt.ylabel('Amplitude [ADC]')
plt.colorbar()
plt.show()


#%%
from functions import multi_fit, isto_carica
    
t_min = int(5320/4)   #4 ns per point
t_max = int(5550/4)
    

_bin = 200 
x = -6000
y = 65000

n, b, charge = isto_carica(a, qf, t, t_min, t_max, _bin, x, y)

plt.figure()
plt.title('Charge histogram 26-07 data 250 MHz')
plt.ylabel('Counts')
plt.xlabel('Charge [ns.ADC]')
plt.hist(charge, b, density=False)
#plt.savefig('ch_histo250', dpi=300)
plt.show()

#%%
from functions import multi_fit

x = np.linspace(-10000, 70000, 1000000)

p0 = [170, -2000, 2000, 240, 3800, 2000, 200, 10000] # density=False in both functions

numg = 8
for i in range(0, numg):
    p0 = np.append(p0, p0[6]/(3**(i+1)))   
p0 = np.array(p0)

pop, cov, SNR, chi2 = multi_fit(numg, p0, n, b, charge, x, fitdue=True, optimiz=False)

print('carica =', round(pop[7]-pop[4], 3), ',  SNR = ', round(SNR, 2), ",  chi2 =", round(chi2, 2))#, j)

    
#%% SIGNALS
from functions import signal

tre = 37
prima = 3
dopo = 54     
time_after_pulse = 100   #400 ns
lim_amp = 450
debug=True
w=2
n=50

carica, e = signal(qf, qf, qff, t, tre, prima, dopo, time_after_pulse, lim_amp, debug, w, n, t_min, t_max)

#  SISTEMA THRESHOLD, PRIMA E DOPO PER OTTENERE BUON GRAFICO O EXCLUSION WINDOW O TIME_AFTER_PULSE!!!!
#  O META' BASELINE bc there's always signal at some point

#%% CHARGE HIST

charge = np.zeros(0)
for i in range(0, len(qf)):
    c=np.zeros(len(np.where(carica[i]!=0)))
    c = carica[i, np.where(carica[i]!=0)]
    charge = np.append(charge, c)
        
    
plt.figure(figsize=[8.5, 6])
ax=plt.gca()
ax.tick_params(bottom=True, top=True, left=True, right=True)
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
plt.ylabel('Counts')
plt.xlabel('Charge [ns.ADC]')
plt.hist(charge, np.arange(-3500, 45000, 300), density=False)
plt.savefig('colbox_aug_spe_dark_counts_rate.png', dpi=600)
# plt.show()

# TRIGGER RATE
rate = sum(e)/( len(a)*len(qf[0])*4e-9*(20*36))  #Hz
print(round(rate,1), 'Hz/mm^2')
print('ev', round(sum(e)*1e-3,1),'k')



























#%%   DATILED 17th DEC

os.chdir(r"H:\Il mio Drive\WORK")
filez="C3MiniArapuca_PlateA_250MHz_LED_width20ns_ampl29V_Vbias46V_00000.txt" 
# 853.7 Hz/mm^2
# ev 5.2 k

a, t = read_long_txt(filez, 3, 938, 4)  # 4ns = 250 Mhz

a=a*1e3  #mV
#%% FILTER
from functions import mov_av

q = mov_av(a, a, 7)


#%% BASELINE                       
from functions import baseline

qf, mu = baseline(q, threshold=30, bins=100)

puliti=np.ones(len(a))
k = np.zeros((len(a), len(a[0])))
for i in range(0, len(a)):
    if min(qf[i,0:600])<-3 or min(qf[i,750:938])<-3:
        puliti[i]=0
    a[i] = a[i] - mu[i]

    
qff = mov_av(a, a, 17)

#%%  VALID EVENTS

k = np.zeros((int(sum(puliti)), len(a[0])))
kf = np.zeros((int(sum(puliti)), len(a[0])))

j=0
for i in range(0, len(a)):
    if puliti[i]==1:
        k[j] = qf[i]
        kf[j] = qff[i]
        j+=1
t0 = t[0:int(sum(puliti))]

#%% HISTOGRAM 2D (color for the z) OF ALL THE WAVEFORMS (PERCISTENCY HISTOGRAM)
xbi = np.linspace(0, 3748, 938)
ybi = np.linspace(-15, 35, 500)

plt.figure(figsize=[8,6])
h, xb, yb = np.histogram2d(t.flatten(), qf.flatten(), bins=(xbi, ybi), normed=None, weights=None, density=None)
h = np.log(h.T)
#X, Y = np.meshgrid(xb, yb)
ax=plt.gca()
ax.tick_params(bottom=True, top=True, left=True, right=True)
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
plt.pcolormesh(xb, yb, h, cmap='jet')
plt.xlabel('Time [ns]')
plt.ylabel('Amplitude [mV]')
plt.colorbar()
plt.show()


#%% CHARGE HIST
from functions import multi_fit, isto_carica
    
t_min = int(2760/4)   #4 ns per point
t_max = int(2840/4)
    

_bin =55
x = -250
y = 1000
n, b, charge = isto_carica(k, k, t0, t_min, t_max, _bin, x, y)


plt.figure(figsize=[8,6])
ax=plt.gca()
ax.tick_params(bottom=True, top=True, left=True, right=True)
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
plt.ylabel('Counts')
plt.xlabel('Charge [ns.V]')
plt.hist(charge, b, density=False)
#plt.savefig('ch_histo250', dpi=300)
plt.show()

#%% MULTIFIT
from functions import multi_fit

x = np.linspace(-250, 1000, 1000)

p0 = [160, 0, 50, 120, 110, 30, 80, 50] # density=False in both functions

numg = 0
for i in range(0, numg):
    p0 = np.append(p0, p0[6]/(3**(i+1)))   
p0 = np.array(p0)

pop, cov, SNR, chi2 = multi_fit(numg, p0, n, b, charge, x, fitdue=False, optimiz=False)

print('carica =', round(pop[7]-pop[4], 3), ',  SNR = ', round(SNR, 2), ",  chi2 =", round(chi2, 2))#, j)
#%% OPTIMIZATION

start1 = 2740    #divisibili per 4???
start2 = 2820
fin1 = 2700
fin2 = 3100  


sn = np.array(0)
mu1 = np.array(0)
mu1dif = np.array(0)
time = np.array(0)

while True:

    t_min = int((start1)/4)   #4 ns per point
    t_max = int((start2)/4)
        
    n, b, charge = isto_carica(k, k, t0, t_min, t_max, _bin, -250, 100)

    pop, cov, SNR, chi2 = multi_fit(numg, p0, n, b, charge, x, fitdue=False, optimiz=True)
    
    time = np.append(time, (t_max-t_min)*4)   #indici --> ns
    mu1 = np.append(mu1, pop[4])
    mu1dif = np.append(mu1dif, pop[7]-pop[4])
    sn = np.append(sn, SNR)
    # if SNR>3.026: print(start1, start2)
    
    if start1-4 > fin1:
        start1 = start1 - 4
    if start2+4 < fin2:
        start2 = start2 + 4
    else: break 

time = np.delete(time, 0)  #indici --> ns
# bn = np.array(0)

# for i in range(70, 300):
    
#     n, b, charge = isto_carica(k, k, t0, t_min, t_max, i, -250, 100)

#     pop, cov, SNR, chi2 = multi_fit(numg, p0, n, b, charge, x, fitdue=False, optimiz=True)
    
#     bn = np.append(bn, i)   #indici --> ns
#     mu1 = np.append(mu1, pop[4])
#     mu1dif = np.append(mu1dif, pop[7]-pop[4])
#     sn = np.append(sn, SNR)


# bn = np.delete(bn, 0)  #indici --> ns
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









#%% SIGNALS
from functions import signal

tre = 1.3
prima = 7
dopo = 18
time_after_pulse = 50   #400 ns
lim_amp = 450
debug=False
w=6
n=50
t_min = int(2760/4)   #4 ns per point
t_max = int(2840/4)

carica, e = signal(k, k, kf, t0, tre, prima, dopo, time_after_pulse, lim_amp, debug, w, n, t_min, t_max)

#  SISTEMA THRESHOLD, PRIMA E DOPO PER OTTENERE BUON GRAFICO O EXCLUSION WINDOW O TIME_AFTER_PULSE!!!!
#  O META' BASELINE bc there's always signal at some point

#%% CHARGE HIST

charge = np.zeros(0)
for i in range(0, len(k)):
    c=np.zeros(len(np.where(carica[i]!=0)))
    c = carica[i, np.where(carica[i]!=0)]
    charge = np.append(charge, c)
        
    
plt.figure(figsize=[8.5, 6])
ax=plt.gca()
ax.tick_params(bottom=True, top=True, left=True, right=True)
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
plt.ylabel('Counts')
plt.xlabel('Charge [mV.ns]')
plt.hist(charge, np.arange(-10, 650, 11), density=False)
#plt.savefig('ch_histo250', dpi=300)
plt.show()

# TRIGGER RATE
rate = sum(e)/( len(k)*len(k[0])*4e-9*(20*36))  #Hz
print(round(rate,1), 'Hz/mm^2')
print('ev', round(sum(e)*1e-3,1),'k')












#%%

file="wave2.dat"   # carica = 1345.185 ,  SNR =  2.51 ,  chi2 = 2.67
header=12     # 1194.1 Hz/mm^2   ev 193.5 k
a, t, _, _, _ = read_binary(file, header, adc_res=1, ns=4) #ns tra i punti, 250Mhz=4ns, 500Mhz=2ns


#%% FILTER

q = mov_av(a, a, 7)


#%%  BASELINE
from functions import baseline

qf, mu = baseline(q, threshold=35, bins=100)

for i in range(0, len(a)):
    a[i] = a[i] - mu[i]
    
qff = mov_av(a, a, 17)


#%%
#HISTOGRAM 2D (color for the z) OF ALL THE WAVEFORMS (PERCISTENCY HISTOGRAM) (2500 bin per x, 1000 per y)

xbi = np.arange(3900, 5000, 4)
ybi = np.arange(-100, 200, 0.4)

plt.figure()
h, xb, yb = np.histogram2d(t.flatten(), qf.flatten(), bins=(xbi, ybi), normed=None, weights=None, density=None)
h = np.log(h.T)
#X, Y = np.meshgrid(xb, yb)
plt.pcolormesh(xb, yb, h, cmap='jet')
plt.xlabel('Time [ns]')
plt.ylabel('Amplitude [ADC]')
plt.colorbar()
plt.show()


#%%
from functions import multi_fit, isto_carica
    
t_min = int(4080/4)   # 4074
t_max = int(4190/4)   # 4240
    

_bin = 80
x = -2550
y = 12000

n, b, charge = isto_carica(a, qf, t, t_min, t_max, _bin, x, y)

plt.figure()
plt.title('Charge histogram 26-07 data 250 MHz')
plt.ylabel('Counts')
plt.xlabel('Charge [ns.ADC]')
plt.hist(charge, b, density=False)
#plt.savefig('ch_histo250', dpi=300)
plt.show()
#%%  MULTIFIT

from functions import multi_fit

x = np.linspace(-2550, 12000, 1000)

p0 = [1100, -400, 800, 950, 1400, 600, 500, 3000] # density=False in both functions

numg = 2
for i in range(0, numg):
    p0 = np.append(p0, p0[6]/(3**(i+1)))   
p0 = np.array(p0)

pop, cov, SNR, chi2 = multi_fit(numg, p0, n, b, charge, x, fitdue=False, optimiz=False)

print('carica =', round(pop[7]-pop[4], 3), ',  SNR = ', round(SNR, 2), ",  chi2 =", round(chi2, 2))#, j)
#%%OPTIMIZATION
start1 = 4120    #divisibili per 4???
start2 = 4150
fin1 = 4000
fin2 = 4350  


sn = np.array(0)
mu1 = np.array(0)
mu1dif = np.array(0)
time = np.array(0)

while True:

    t_min = int((start1)/4)   #4 ns per point
    t_max = int((start2)/4)
        
    n, b, charge = isto_carica(qf, qf, t, t_min, t_max, _bin, -5550, 25000)

    pop, cov, SNR, chi2 = multi_fit(numg, p0, n, b, charge, x, fitdue=False, optimiz=True)
    
    time = np.append(time, (t_max-t_min)*4)   #indici --> ns
    mu1 = np.append(mu1, pop[4])
    mu1dif = np.append(mu1dif, pop[7]-pop[4])
    sn = np.append(sn, SNR)
    # if SNR>2.45: print(start1, start2)
    
    if start1-4 > fin1:
        start1 = start1 - 4
    if start2+4 < fin2:
        start2 = start2 + 4
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









#%%  SIGNALS
from functions import signal
tre = 12.5
prima = 3
dopo = 25
time_after_pulse = int(345/4)
lim_amp = 450
debug=False
num_plot=4
num_ev_expected=60
t1=int(4050/4)
t2=int(4250/4)

carica, e = signal(a, qf, qff, t, tre, prima, dopo, time_after_pulse, lim_amp, debug, num_plot, num_ev_expected, t1, t2)

#  %time o  %%time per tempo esecuzione riga o cella direttamente sulla console (iPyhton)

#%% CHARGE HIST

charge = np.zeros(0)
for i in range(0, len(a)):
    c=np.zeros(len(np.where(carica[i]!=0)))
    c = carica[i, np.where(carica[i]!=0)]
    charge = np.append(charge, c)
        
    
plt.figure(figsize=[8.5, 6])
ax=plt.gca()
ax.tick_params(bottom=True, top=True, left=True, right=True)
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
# font1 = {'family':'serif','weight':'bold','size':20}

plt.ylabel('Counts')
plt.xlabel('Charge [ns.ADC]')
plt.hist(charge, np.arange(-1000, 10000, 50), density=False)
plt.savefig('wave2_ariade_dark_cout_search', dpi=600)
# plt.show()

# TRIGGER RATE
rate = sum(e)/( len(a)*len(a[0])*4e-9*(20*36))  #Hz
print(round(rate,1), 'Hz/mm^2', 'ev', round(sum(e)*1e-3, 1),'k')





