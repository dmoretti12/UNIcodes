# -*- coding: utf-8 -*-
"""
Created on Wed May  4 11:11:48 2022

@author: Davide
"""
import numpy as np
import matplotlib.pyplot as plt
#import math
#from scipy.optimize import curve_fit

from functions import norm_fit
# import scipy.integrate as integrate
# import scipy.special as special

#from sklearn.metrics import mean_squared_error

import os
import glob
import pandas as pd

#from denoising import Denoiser

plt.close("all") #chiude tutte le finestre

os.chdir(r"C:\Users\Davide\Desktop\20220502_snr_koheron")

filenames = [i for i in glob.glob("C2TraceC*.txt")]

#df = pd.read_csv('C1XArapuca_CRT_250MHz_00000.txt', sep=" ", header=5, index_col = None, encoding = 'unicode_escape')
df = [pd.read_csv(file, sep = " ", header=4, engine = 'python', index_col = None, encoding = 'unicode_escape')
      for file in filenames]

# IL FILE C1X...042 = df[42] HA UNA RIGA IN MENO, L'HO AGGIUNTA A CASO


t = np.empty((len(df),len(df[0])))      #freq 250 Mhz
a = np.empty((len(df),len(df[0])))

for i in range (0, len(df)):
    
    q = df[i].iloc[:, 0].values
    w = df[i].iloc[:, 1].values
    
    if len(df[i]) < len(df[0]):
        
        q = np.append(q, 0)
        w = np.append(w, 0)
    
   # T = np.max(t) + df[i].iloc[:, 0].values #this way (iloc) to access to the values of df visualizing le waveform tutte di seguito
    t[i,:] = q
    a[i,:] = w
    
plt.figure()
plt.plot(t[0], a[0], '-')
plt.plot(t[1], a[1], '.')
plt.show()


#%% MEANS & STD_DEVS

mu0, std0 = norm_fit(df, a, t, -2.5e-07, 0.4e-07)

mu1, std1 = norm_fit(df, a, t, 1.1e-07, 1.6e-07)


#%% SNR

snr = np.empty(len(df))

snr_sum = 0
snr_mean = 0

for i in range(0, len(df)):
    
    snr[i] = (mu1[i]-mu0[i])/std1[i]
    
    snr_sum = snr_sum + snr[i]
    

snr_mean = snr_sum / len(df)

print(snr_mean)
    
#%% PLOT

kk = np.genfromtxt('k.txt', skip_header=5)#, unpack=True)
nr = np.genfromtxt('newr.txt', skip_header=5)

t0=kk[:,0]
a0=kk[:,1]

t1=nr[:,0]
a1=nr[:,1]

t0 = t0*1000000   # time in us
t1 = t1*1000000
a0 = a0 * 5   #koheron scale a 2mV instead new_rec a 10 mV
a1 = a1 + 0.631  #offset new_rec a -611

plt.figure()
plt.plot(t0, a0, 'r', label = 'Koheron: 9.88 mV, 4.50 nV.s, 10.36 RT. SNR: 23.7')
plt.plot(t1, a1, 'c', label = 'New_receiver: 34.4 mV, -306.1 nV.s, 312.70 RT. SNR: 3.6')
plt.title('Koheron and new receiver waveforms')
plt.legend(loc='upper left', fontsize = 'large')
plt.ylabel('Receiver [mV]')
plt.xlabel('Time [us]')
#plt.tight_layout()
plt.savefig('koheron_vs_new_receiver', dpi=300)
#plt.show()


