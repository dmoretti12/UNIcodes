# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 12:04:05 2022

@author: Davide
"""

from functions import read_many_txt, read_long_txt, denoising, norm_fit, signal, isto_carica, tot, ampl_pe, mov_av, multi_fit, gauss, multi_gauss
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
import random


plt.close("all") #chiude tutte le finestre

def expo(x, a, b):
    return np.exp(a*x + b)

def retta(x, m, q):
    return m*x + q

#  DATA FROM 16th DEC mini arapuca

file = 'C2MiniArapuca_A4ch1_250MHz_LED_width20ns_ampl25V_a_00000.txt'
righe_iniz = 3
puntixwvf = 938

aa, tt, aa62, tt62 = read_long_txt(file, righe_iniz, puntixwvf)


#%% charge 250 vs 62 data 16 dec

tt_ns = tt * 1000   #ns
tt62_ns = tt62 * 1000 #ns

carr, peaa = ampl_pe(aa, aa, tt_ns, t_min=int((2.736*1000)/4), t_max=int((2.880*1000)/4))
carr62, peaa62 = ampl_pe(aa62, aa62, tt62_ns, t_min=int((2.736*1000)/16), t_max=int(((2.880*1000)/16)))

popt, pcov = curve_fit(retta, carr62, carr, p0=[1., 0.])   
txt='m= '+str(round(popt[0], 4))+' b= '+str(round(popt[1], 4))
xx = np.linspace(0, 1000, 1000)
plt.figure()
plt.plot(carr62, carr, '.')
plt.plot(xx, retta(xx, *popt), '-', c='r', label = txt)
plt.title('Charge 250 Mhz vs 62.5 Mhz (16Dec-C2MiniArapuca_A4ch1_250MHz_LED_width20ns_ampl25V)')
plt.ylabel('Charge 250 Mhz [ns.V]')
plt.xlabel('Charge 62 Mhz [ns.V]')
plt.grid()
plt.legend()
plt.show()
# plt.savefig('c250vs62_16dec', dpi=300)

#è una merda

#%% CHARGE HISTO 16 dec miniarapuca 250 MHz
from functions import multi_fit

tt_ns = tt * 1000  #us -> ns

t_min = int((2.736*1000)/4)   #4 ns per point
t_max = int((2.880*1000)/4)
mm = np.array(0)
# for j in range (0, 200):
    

_bin=143 #115 #80+j #143  

n, b, charge = isto_carica(aa, aa, tt_ns, t_min, t_max, _bin)

plt.figure()
plt.title('Charge histogram (miniarapuca16thDecLED) 250 MHz')
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
# mm = np.append(mm, SNR)

# mm = np.delete(mm, 0)
    
# print(np.max(mm), np.argmax(mm))

#%%250MHZ
3 gaussian bin=143, 2.736,2.864  SNR=1.32
    2.880  snr=1.28   dif=0.16 bin=143
    2.884   snr =1.31  bin=84  dif 0.17
    2.888 snr =1.27   0.18  bin=80+43 (44)
    2.892 snr=1.23   0.178  bin=80+10  (11)
    2.896 snr = 1.21 dif = 0.18
4 gaussian bin=84, 2.736, 2.884, snr=1.33  dif=0.167
    2.740 bin =80+17   snr=1.32  0.165

a = dif250 optimezed
b= dif62 optimezed
c= dif250 with 62 setup


b/a = 6.1%
c/a = 7.7%
b/c = 1.5%

if 1 of these is more than 10-15% the difference between 250 and 62 is high


#%% CHARGE HISTO 16 dec miniarapuca 62.5 MHz
from functions import multi_fit

tt62_ns = tt62 * 1000 #us -> ns

t_min62 = int((2.720*1000)/16)   #*1000 x ns -> 16 ns per point
t_max62 = int((2.880*1000)/16)
mm = np.array(0)


# for j in range(0, 250):
    
_bin= 170 #60+j

n62, b62, charge62 = isto_carica(aa62, aa62, tt62_ns, t_min62, t_max62, _bin)
# plt.figure()
# plt.title('Charge histogram (miniarapuca16thDecLED) 62.5 MHz')
# plt.ylabel('Counts')
# plt.xlabel('Charge [ns.V]')
# plt.hist(charge62, b62, density=False)
# plt.grid(True)
# #plt.savefig('ch_histo62', dpi=300)
# plt.show()

x = np.linspace(-0.35, 1.9, 1000)
p0 = [2, -0.025, 0.08, 2, 0.14, 0.10, 1.2, 0.3]
#p0 = np.append(p0, [p0[5]*math.sqrt(2)]) 
#p0=[1.8, -0.05, 0.08, 1.6, 0.13, 0.07, 1, 0.3]   #parametri molto giusti
numg62 = 3
for i in range(0, numg62):
    p0 = np.append(p0, p0[6]/(3**(i+1)))#, p0[7+3*i]+p0[7]-p0[4], math.sqrt(2+i)*p0[8+3*i]])
    
p0 = np.array(p0)

pop62, cov62, SNR62, chi2_62 = multi_fit(numg62, p0, n62, b62, charge62, x)

# mm = np.append(mm, SNR62)

print(pop62[7]-pop62[4], SNR62, "chi2 = ", chi2_62)#, j)

# mm = np.delete(mm,0)

# print(np.max(mm), np.argmax(mm))

#%% 62MHZ
4 gaussian SNR=1.36 i=36, 2.736,2.912
2.912 -> 2.836 snr=1.18 i=5 