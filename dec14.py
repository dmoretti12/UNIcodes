# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 19:31:28 2022

@author: Davide
"""

from functions import read_long_txt_time, denoising, norm_fit, signal, isto_carica, tot, ampl_pe, mov_av, multi_fit, gauss, multi_gauss, time_res
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit


plt.close("all") #chiude tutte le finestre

def expo(x, a, b):
    return np.exp(a*x + b)

def retta(x, m, q):
    return m*x + q

#  DATA FROM 16th DEC mini arapuca

file = 'C1XArapuca_Efield_0kV_CRPon_A1ch1_250MHz_LEDwith20nsampl30V00001.txt'   #apri file txt e guarda linee
righe_iniz = 504
puntixwvf = 1252

a14, t14, a14_62, t14_62 = read_long_txt_time(file, righe_iniz, puntixwvf)

#%%amp vs car 250
from functions import ampl_pe
# a_filt14 = mov_av(a14, a14)
t14_ns = t14 * 1e9  #s -> ns
carica, peak = ampl_pe(a14, a14, t14_ns, t_min=int((0.23-min(t14[0])*1e6)*1000/4), t_max=int((0.5-min(t14[0])*1e6)*1000/4))


plt.figure()
plt.title('Amplitude vs Charge 14Dec LED 30V, 250 MHz')
plt.ylabel('Amplitude [V]')
plt.xlabel('Charge [V.ns]')
plt.plot(carica, peak, '.')
plt.grid(True)
# plt.savefig('amp_ch250_C1X-14dec', dpi=300)
plt.show()

#%% amp vs car 62

from functions import ampl_pe
# a_filt14 = mov_av(a14, a14)
t14_62_ns = t14_62 * 1e9  #s -> ns
carica, peak = ampl_pe(a14_62, a14_62, t14_62_ns, t_min=int((0.23-min(t14_62[0])*1e6)*1000/16), t_max=int((0.5-min(t14_62[0])*1e6)*1000/16))


plt.figure()
plt.title('Amplitude vs charge 14Dec LED 30V, 62.5 MHz')
plt.ylabel('Amplitude [V]')
plt.xlabel('Charge [V.ns]')
plt.plot(carica, peak, '.')
plt.grid(True)
# plt.savefig('amp_ch62_C1X-14dec', dpi=300)
plt.show()


#%% tot 250

from functions import tot

tre_tot14 = 0.02
t14_ns = t14 * 1e9  #s -> ns
ctot_14, ss_tot14, peaktot_14, peakpos14 = tot(a14, a14, t14_ns, t_min=int((0.23-min(t14[0])*1e6)*1000/4), t_max=int((0.5-min(t14[0])*1e6)*1000/4), threshold=tre_tot14, satur=1.4)

deltaT14 = np.zeros(len(a14))
for i in range (0, len(a14)):
    deltaT14[i] = (ss_tot14[i,1] - ss_tot14[i,0])

satur = 1.4
      
plt.figure()
plt.semilogy()
plt.plot(deltaT14[peaktot_14<satur], peaktot_14[peaktot_14<satur], '.b', label='Non saturated')
plt.plot(deltaT14[peaktot_14>=satur], peaktot_14[peaktot_14>=satur], '.r', label='Saturated')
plt.title('Amplitude vs TOT 14Dec LED 30V, 250 MHz')
plt.ylabel('Amplitude [V]')
plt.xlabel('TOT [ns]')
plt.legend()
plt.grid(True, which='both')
# plt.savefig('amp_tot250_C1X-14dec', dpi=300)
plt.show()

pu=deltaT14[peaktot_14<satur]
po=ctot_14[peaktot_14<satur]
q=pu[pu<250]
qq=po[pu<250]
popt, pcov = curve_fit(expo, q[q>70], qq[q>70], p0=[0.007, 3.8])
fittxt='m= '+str(round(popt[0], 4))+' b= '+str(round(popt[1], 4))

dev = math.sqrt(sum((qq[q>70]-expo(q[q>70], *popt))**2)/(len(q[q>70])-1))
#print('dev', dev)
dev_rel = dev / 300.

ret = deltaT14[peaktot_14>=satur]
yret = expo(ret, *popt)
for i in range(0, len(ret)):
        yret[i] = np.random.normal(yret[i], dev_rel*yret[i])

plt.figure()
plt.semilogy()
plt.plot(deltaT14[peaktot_14<satur], ctot_14[peaktot_14<satur], '.b', label='Non saturated')
plt.plot(deltaT14[deltaT14<300], expo(deltaT14[deltaT14<300], *popt), 'g', label=fittxt)
plt.plot(deltaT14[peaktot_14>=satur], ctot_14[peaktot_14>=satur], '.r', label='Saturated')
plt.plot(deltaT14[peaktot_14>=satur], yret, '.g', label = 'Expected')
plt.title('Charge vs TOT 14Dec LED 30V, 250 MHz')
plt.ylabel('Charge [V.ns]')
plt.xlabel('TOT [ns]')
plt.legend()
plt.grid(True, which='both')
# plt.savefig('ch_tot250_C1X-14dec', dpi=300)
plt.show()    

#%% tot 62

from functions import tot

tre_tot14_62 = 0.02

t14_62_ns = t14_62 * 1e9  #us -> ns
ctot_14_62, ss_tot14_62, peaktot_14_62, peakpos62 = tot(a14_62, a14_62, t14_62_ns, t_min=int((0.23-min(t14_62[0])*1e6)*1000/16), t_max=int((0.5-min(t14[0])*1e6)*1000/16), threshold=tre_tot14_62, satur=0.1) 

deltaT14_62 = np.zeros(len(a14_62))
for i in range (0, len(a14_62)):
    deltaT14_62[i] = (ss_tot14_62[i,1] - ss_tot14_62[i,0])
satur = 1.4
      
plt.figure()
plt.semilogy()
plt.plot(deltaT14_62[peaktot_14_62<satur], peaktot_14_62[peaktot_14_62<satur], '.b', label='Non saturated')
plt.plot(deltaT14_62[peaktot_14_62>=satur], peaktot_14_62[peaktot_14_62>=satur], '.r', label='Saturated')
plt.title('Amplitude vs TOT 14Dec Cathode 62.5 MHz')
plt.ylabel('Amplitude [V]')
plt.xlabel('TOT [ns]')
plt.legend()
plt.grid(True, which='both')
# plt.savefig('test', dpi=300)
plt.show()

pu=deltaT14_62[peaktot_14_62<satur]
po=ctot_14_62[peaktot_14_62<satur]
q=pu[pu<250]
qq=po[pu<250]
popt, pcov = curve_fit(expo, q[q>70], qq[q>70], p0=[0.007, 3.8])
fittxt='m= '+str(round(popt[0], 4))+' b= '+str(round(popt[1], 4))

dev = math.sqrt(sum((qq[q>70]-expo(q[q>70], *popt))**2)/(len(q[q>70])-1))
#print('dev', dev)
dev_rel = dev / 300.

ret = deltaT14_62[peaktot_14_62>=satur]
yret = expo(ret, *popt)
for i in range(0, len(ret)):
        yret[i] = np.random.normal(yret[i], dev_rel*yret[i])

    

plt.figure()
plt.semilogy()
plt.plot(deltaT14_62[peaktot_14_62<satur], ctot_14_62[peaktot_14_62<satur], '.b', label='Non saturated')
plt.plot(deltaT14_62[deltaT14_62<300], expo(deltaT14_62[deltaT14_62<300], *popt), 'g', label=fittxt)
plt.plot(deltaT14_62[peaktot_14_62>=satur], ctot_14_62[peaktot_14_62>=satur], '.r', label='Saturated')
plt.plot(deltaT14_62[peaktot_14_62>=satur], yret, '.g', label = 'Expected')
plt.title('Charge vs TOT 14Dec Cathode 62.5 MHz')
plt.ylabel('Charge [V.ns]')
plt.xlabel('TOT [ns]')
plt.legend()
plt.grid(True, which='both')
# plt.savefig('tes', dpi=300)
plt.show()   

#%% time res 250
from functions import time_res, gauss

tre=0.01
t14_ns = t14 * 1e9
pic, pospic, sigpos = time_res(a14, a14, t14_ns, t_min=int((0.23-min(t14[0])*1e6)*1000/4), t_max=int((0.5-min(t14[0])*1e6)*1000/4), threshold = tre)

    
bb = np.linspace(277, 313, 25)
n, b2 = np.histogram(pospic[pospic>290], bb)
l = (b2[2]-b2[1])/2
bbin = l + b2  # np.resize(b, len(b))
b1 = np.delete(bbin, len(bb)-1)
popt, pcov = curve_fit(gauss, b1, n, maxfev=10000, bounds=(-0.055, 600))
y1 = gauss(np.linspace(min(bb), max(bb), 1000), *popt)

fittxt = 'Dev_std = '+str(round(popt[2], 3))
plt.figure()
plt.title('Time resolution - 14DecLED30V data, 250 MHz')
plt.ylabel('Counts')
plt.xlabel('Peak position [ns]')
plt.hist(pospic, bins=bb)
plt.plot(np.linspace(min(bb), max(bb), 1000), y1, c='r', label=fittxt, linewidth=2, linestyle='-')
plt.legend()
# plt.savefig('trpic250-14dec', dpi=300)
plt.show()

b = np.linspace(236, 245, 19)
l = (b[2]-b[1])/2
bbin2 = l + b  # np.resize(b, len(b))
b2 = np.delete(bbin2, len(b)-1)
n2, b3 = np.histogram(sigpos, b)
popt, pcov = curve_fit(gauss, b2, n2, maxfev=10000, bounds=(-0.055, 600))
y2 = gauss(np.linspace(min(b), max(b), 1000), *popt)

fitxt='Dev_std = '+str(round(popt[2], 3))
plt.figure()
plt.title('Time resolution - 14DecLED30V data, 250 MHz')
plt.ylabel('Counts')
plt.xlabel('Signal over threshold position [ns]')
plt.hist(sigpos, bins=b)
plt.plot(np.linspace(min(b), max(b), 1000), y2, c='r', label=fitxt, linewidth=2, linestyle='-')
plt.legend()
# plt.savefig('trtre250-14dec', dpi=300)
plt.show()

#%% time res 62

from functions import time_res, gauss

tre=0.01
t14_62_ns = t14_62 * 1e9
pic, pospic, sigpos = time_res(a14_62, a14_62, t14_62_ns, t_min=int((0.23-min(t14_62[0])*1e6)*1000/16), t_max=int((0.5-min(t14_62[0])*1e6)*1000/16), threshold = tre)

    
bb = np.linspace(277, 313, 27)
n, b2 = np.histogram(pospic[pospic>290], bb)
l = (b2[2]-b2[1])/2
bbin = l + b2  # np.resize(b, len(b))
b1 = np.delete(bbin, len(bb)-1)
popt, pcov = curve_fit(gauss, b1, n, maxfev=10000, bounds=(-0.055, 600))
y1 = gauss(np.linspace(min(bb), max(bb), 1000), *popt)

fittxt = 'Dev_std = '+str(round(popt[2], 3))
plt.figure()
plt.title('Time resolution - 14DecLED30V data, 62.5 MHz')
plt.ylabel('Counts')
plt.xlabel('Peak position [ns]')
plt.hist(pospic, bins=bb)
plt.plot(np.linspace(min(bb), max(bb), 1000), y1, c='r', label=fittxt, linewidth=2, linestyle='-')
plt.legend()
# plt.savefig('trpic62-14dec', dpi=300)
plt.show()

b = np.linspace(231, 248, 17)
l = (b[2]-b[1])/2
bbin2 = l + b  # np.resize(b, len(b))
b2 = np.delete(bbin2, len(b)-1)
n2, b3 = np.histogram(sigpos, b)
popt, pcov = curve_fit(gauss, b2, n2, maxfev=10000, bounds=(-0.055, 600))
y2 = gauss(np.linspace(min(b), max(b), 1000), *popt)

fitxt='Dev_std = '+str(round(popt[2], 3))
plt.figure()
plt.title('Time resolution - 14DecLED30V data, 62.5 MHz')
plt.ylabel('Counts')
plt.xlabel('Signal over threshold position [ns]')
plt.hist(sigpos, bins=b)
plt.plot(np.linspace(min(b), max(b), 1000), y2, c='r', label=fitxt, linewidth=2, linestyle='-')
plt.legend()
# plt.savefig('trtre62-14dec', dpi=300)
plt.show()

