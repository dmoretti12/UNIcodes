# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 17:02:54 2022

@author: Davide
"""

from functions import read_binary, read_long_txt, norm_fit, ampl_pe, isto_carica, denoising, mov_av, signal
import numpy as np
import matplotlib.pyplot as plt
#import math
import os
from scipy.stats import norm

plt.close('all')
#%% DATI AGOSTO 22 COLDBOX


os.chdir("H:\Il mio Drive\WORK\AUG22COLDBOX")
file="0_wave1_light_leakage_test_led_fibers_connected_dcem_digital_off.dat"
header=12

a, t, _ = read_binary(file, header, adc_res=1, ns=2) #ns tra i punti, 250Mhz=4ns, 500Mhz=2ns
os.chdir("H:\Il mio Drive\WORK")


#%%  BASELINE

a0 = np.zeros((len(a), len(a[0])))

for i in range (0, len(a)):
    
    mu, _ = norm.fit(a[i])
    a0[i] = a[i] - mu
    
#%% FILTER

qff = denoising(a, a0, 1500)
qf = mov_av(a, a0, 13)


#%%

tre = 41
prima = 100
dopo = 400

carica, peak, peakpos, sop, sot, e = signal(a, qf, qff, t, tre, prima, dopo)

#  %time o  %%time per tempo esecuzione riga o cella direttamente sulla console (iPyhton)

#%%

charge = np.zeros(0)
for i in range(0, len(a)):
    c=np.zeros(len(np.where(carica[i]!=0)))
    c = carica[i, np.where(carica[i]!=0)]
    charge = np.append(charge, c)
        
#%%     
plt.figure()
plt.title('Charge histogram 26-07 data 250 MHz')
plt.ylabel('Counts')
plt.xlabel('Charge [ns.V]')
plt.hist(charge, np.arange(-10000, 65000, 30), density=False)
#plt.savefig('ch_histo250', dpi=300)
plt.show()



















#%% DATI AUG22 SPE



#%%   DATILED 17th DEC

os.chdir(r"H:\Il mio Drive\WORK")
filez="C3MiniArapuca_PlateA_250MHz_LED_width20ns_ampl29V_Vbias46V_00000.txt"

a17, t17, a17_62, t17_62 = read_long_txt(filez, 3, 938)

