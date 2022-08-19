# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 14:59:39 2022

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

string = r"C:\Users\Davide\Desktop\2022_06_16_wf_example_new_board_and_sipm\6_3LED"
iniz = "C1*.txt"

df, a, t, ann, tnn = read_many_txt(string, iniz)

#%% pe 250

from functions import multi_fit

t_ns = t * 1000  #us -> ns

t_min = int((2.736*1000)/4)   #4 ns per point
t_max = int((2.880*1000)/4)
mm = np.array(0)
# for j in range (0, 200):
    

_bin=143 #115 #80+j #143  

n, b, charge = isto_carica(df, a, t_ns, t_min, t_max, _bin)

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