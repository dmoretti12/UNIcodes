# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 18:35:19 2022

@author: Davide
"""

from functions import read_many_txt, read_long_txt_time, norm_fit, signal, gauss
import numpy as np
import matplotlib.pyplot as plt
import math, os
from scipy.optimize import curve_fit
from scipy.stats import norm

#%% wvf0
plt.close("all") #chiude tutte le finestre


# string = r"H:\Il mio Drive\WORK\20220621_fibers"
# iniz = 'C1*.txt'
string = r"C:\Users\Davide\Desktop\20220810_dcem36_f9_fiber_test"
iniz = 'C1*.txt'
df, a, t, ann, tnn = read_many_txt(string, iniz)
t = t*1e9

from matplotlib.cm import get_cmap
cmap = get_cmap(name='jet')   # or 'gist_rainbow'
colors = cmap(np.linspace(0, 1, 14))

txt = [50, 25, 20, 15, 11, 10, 9.5, 9, 7, 4, 3, 2, 1, 0.2]

plt.figure()
for i in range(0, len(df)):
    
    if i==9: continue

    tx = ['Input '+str(txt[i])+' mV']
    
    plt.plot(t[i], a[i], color=colors[i], label=tx)
    
plt.title('DCem#36 v_2, black fiber F-09 core 62.5um')
plt.xlabel('Time [ns]')
plt.ylabel('Amplitude [mV]')
plt.legend()
# plt.savefig('plott', dpi=300)
plt.show()

#%%

string = r"C:\Users\Davide\Desktop\20220810_dcem36_orange_fiber_test"
iniz = 'C1*.txt'
df, a, t, ann, tnn = read_many_txt(string, iniz)
t = t*1e9

plt.figure()
for i in range(0, len(df)):
    
    # if i==9: continue

    tx = ['Input '+str(txt[i])+' mV']
    
    plt.plot(t[i], a[i], color=colors[i], label=tx)
    
plt.title('DCem#36 v_2, orange fiber core 200um')
plt.xlabel('Time [ns]')
plt.ylabel('Amplitude [mV]')
plt.legend()
# plt.savefig('plott', dpi=300)
plt.show()