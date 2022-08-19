# -*- coding: utf-8 -*-
"""
Created on Mon May  9 14:38:24 2022

@author: Davide
"""
from functions import read_many_txt, read_long_txt, denoising, norm_fit, signal, isto_carica, tot, ampl_pe, mov_av, multi_fit, gauss, multi_gauss
# from pylab import*
#from scipy.signal import find_peaks
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
#from scipy.optimize import leastsq
import random
# import pandas as pd
# from pathlib import Path  

# filepath = Path(r"C:\Users\Davide\Documents\FAMIGLIA\Davide\UNIV\WORK\test.csv")  

# filepath.parent.mkdir(parents=True, exist_ok=True)  

# import scipy.integrate as integrate
# import scipy.special as special

#from sklearn.metrics import mean_squared_error
#from denoising import Denoiser

plt.close("all") #chiude tutte le finestre

def expo(x, a, b):
    return np.exp(a*x + b)

def retta(x, m, q):
    return m*x + q


#  data 15th dec CRT C1X

string = r"H:\Il mio Drive\WORK\Dec15th_CRT-20220321T101331Z-001\Dec15th_CRT"
iniz = "C1X*.txt"

df, a, t, ann, tnn = read_many_txt(string, iniz)

# plt.figure()
# plt.plot(tnn[17], ann[17], '-')
# plt.plot(tnn[17], sf[17], '-')
# plt.title('C1X 15Dec CRT 62.5 MHz')
# plt.ylabel('Amplitude [V]')
# plt.xlabel('Time [s]')
# plt.savefig('a62_C1X', dpi=300)

#%% TOT (time over threshold)  250 MHz \Dec15th_CRT

t_ns = t * 1e9  #s -> ns
# tre_tot = np.full(len(df), 0.4)
tre_tot = 0.45
#c_tot, ss_tot, peak_tot, peakpos_tot = signal(df, 1, a, t, tre_tot, 15, 40)
#_, _, tot_c, tot_peak = isto_carica(df, a, t, t_min=int((-0.3e-01*1000)/4), t_max=int((1.3*1000)/4), _bin=10) 
#bisogna capire quali wvf sono prese, quelle sopra la threshold
ss_tot = np.zeros((len(df), 2))
c_tot, ss_tot, peak_tot, ppos, wvf = tot(df, a, t_ns, t_min=int(((1.001111-0.03)*1000)/4), t_max=int((1.7*1000)/4), threshold=tre_tot, satur=1.4)
# sstot_ = np.array(0)
# peaktot_= np.array(0)
# deltaT=np.array(0)
# ctot_=np.array(0)

# for i in range(0, len(c_tot)):
#     sstot_ = np.append(sstot_, ss_tot[i])
#     peaktot_ = np.append(peaktot_, peak_tot[i])
#     ctot_ = np.append(ctot_, c_tot[i])
# sstot_ = np.delete(sstot_, 0)
# peaktot_ = np.delete(peaktot_, 0)
# ctot_ = np.delete(ctot_, 0)
# i=0
# while i < len(sstot_)-1:
#         # if sstot_[i]==0 and sstot_[i+1]==0:
#         #     deltaT = np.append(deltaT, 0)
#         # else:
#         deltaT = np.append(deltaT, sstot_[i+1]-sstot_[i])
#         i = i+2
# deltaT = np.delete(deltaT, 0)
# deltaT = deltaT * 4 #  in nanosec               / 1e9   # in secondi
deltaT = np.zeros(len(df))
for i in range (0, len(df)):
    deltaT[i] = (ss_tot[i,1] - ss_tot[i,0])   #già in ns, altrimenti *e-03
    # if deltaT[i]==0:
    #     deltaT = np.delete(deltaT, i)
    #     # ss_tot = np.delete(ss_tot, i)
    #     c_tot = np.delete(c_tot, i)

saturation = 1.4

plt.figure()
plt.semilogy()
plt.plot(deltaT[peak_tot<saturation], peak_tot[peak_tot<saturation], '.b', label='Non saturated')
plt.plot(deltaT[peak_tot>=saturation], peak_tot[peak_tot>=saturation], '.r', label='Saturated')
plt.title('Amplitude vs TOT 15Dec CRT 250 MHz')
plt.ylabel('Amplitude [V]')
plt.xlabel('TOT [ns]')
plt.legend()
plt.grid(True, which='both')
# plt.savefig('amp_tot250_C1X', dpi=300)
plt.show()

pu=deltaT[peak_tot<saturation]
po=c_tot[peak_tot<saturation]
q=pu[pu<250]
qq=po[pu<250]
popt, pcov = curve_fit(expo, q[q>70], qq[q>70], p0=[0.007, 3.8])  
fittxt='m= '+str(round(popt[0], 4))+' b= '+str(round(popt[1], 4))

dev = math.sqrt(sum((qq[q>70]-expo(q[q>70], *popt))**2)/(len(q[q>70])-1))
dev_rel = dev / 300.

ret = deltaT[peak_tot>=saturation]
yret = expo(ret, *popt)
for i in range(0, len(ret)):
        yret[i] = np.random.normal(yret[i], dev_rel*yret[i])

plt.figure()
plt.semilogy()
plt.plot(deltaT[peak_tot<saturation], c_tot[peak_tot<saturation], '.b', label='Non saturated')
plt.plot(deltaT[deltaT<300], expo(deltaT[deltaT<300], *popt), 'g', label=fittxt)
plt.plot(deltaT[peak_tot>=saturation], c_tot[peak_tot>=saturation], '.r', label='Saturated')
plt.plot(deltaT[peak_tot>=saturation], yret, '.g', label = 'Expected')
plt.title('Charge vs TOT 15Dec CRT 250 MHz')
plt.ylabel('Charge [V.ns]')
plt.xlabel('TOT [ns]')
plt.legend()
plt.grid(True, which='both')
# plt.savefig('ch_tot250_C1X', dpi=300)
plt.show()      

#%% TOT (time over threshold)  62.5 MHz  \Dec15th_CRT

tnn_ns = tnn * 1e9    # nanosec
# tre_tot = np.full(len(df), 0.4)
tre_tot = 0.4
ctot_62, ss_tot62, peaktot_62, peakpos_tot62, wvf = tot(df, ann, tnn_ns,  t_min=int(((1.001111-0.03)*1000)/16), t_max=int((1.7*1000)/16), threshold=tre_tot, satur=1.4)

# sstot_62 = np.array(0)
# peaktot_62= np.array(0)
# deltaT62=np.array(0)
# ctot_62=np.array(0)

# for i in range(0, len(df)):
#     sstot_62 = np.append(sstot_62, ss_tot62[i])
#     peaktot_62 = np.append(peaktot_62, peak_tot62[i])
#     ctot_62 = np.append(ctot_62, c_tot62[i])
# sstot_62 = np.delete(sstot_62, 0)
# peaktot_62 = np.delete(peaktot_62, 0)
# ctot_62 = np.delete(ctot_62, 0)
# i=0
# while i < len(sstot_62)-1:
#         deltaT62 = np.append(deltaT62, sstot_62[i+1]-sstot_62[i])
#         i = i+2
# deltaT62 = np.delete(deltaT62, 0)
# deltaT62 = deltaT62 * 16 #in nanosec           / 1e9   # in secondi

deltaT62 = np.zeros(len(df))
for i in range (0, len(df)):
    deltaT62[i] = (ss_tot62[i,1] - ss_tot62[i,0])
    
saturation = 1.4

plt.figure()
plt.semilogy()
plt.plot(deltaT62[peaktot_62<saturation], peaktot_62[peaktot_62<saturation], '.b', label='Non saturated')
plt.plot(deltaT62[peaktot_62>=saturation], peaktot_62[peaktot_62>=saturation], '.r', label='Saturated')
plt.title('Amplitude vs TOT 15Dec CRT 62.5 MHz')
plt.ylabel('Amplitude [V]')
plt.xlabel('TOT [ns]')
plt.legend()
plt.grid(True, which='both')
# plt.savefig('amp_tot62_C1X', dpi=300)
plt.show()

pu=deltaT62[peaktot_62<saturation]
po=ctot_62[peaktot_62<saturation]
q=pu[pu<250]
qq=po[pu<250]
popt, pcov = curve_fit(expo, q[q>70], qq[q>70], p0=[0.007, 3.8])   
fittxt='m= '+str(round(popt[0], 4))+' b= '+str(round(popt[1], 4))

dev = math.sqrt(sum((qq[q>70]-expo(q[q>70], *popt))**2)/(len(q[q>70])-1))
dev_rel = dev / 300.

ret = deltaT62[peaktot_62>=saturation]
yret = expo(ret, *popt)
for i in range(0, len(ret)):
        yret[i] = np.random.normal(yret[i], dev_rel*yret[i])

plt.figure()
plt.semilogy()
plt.plot(deltaT62[peaktot_62<saturation], ctot_62[peaktot_62<saturation], '.b', label='Non saturated')
plt.plot(deltaT62[deltaT62<300], expo(deltaT62[deltaT62<300], *popt), 'g', label=fittxt)
plt.plot(deltaT62[peaktot_62>=saturation], ctot_62[peaktot_62>=saturation], '.r', label='Saturated')
plt.plot(deltaT62[peaktot_62>=saturation], yret, '.g', label = 'Expected')
plt.title('Charge vs TOT 15Dec CRT 62.5 MHz')
plt.ylabel('Charge [V.ns]')
plt.xlabel('TOT [ns]')
plt.legend()
plt.grid(True, which='both')
# plt.savefig('ch_tot62_C1X', dpi=300)
plt.show()

#%%  250 MHz arrays a & t  \Dec15th_CRT

a_filt = denoising(df, a)
t_ns = t * 1e9  #s -> ns
inizio_sig=-1e-07*(1e9)   #nanosec
mu, std = norm_fit(df, a_filt, t_ns, inizio_sig)

sigma = 4 #thresold dev_std sopra baseline x discriminaz segnale
threshold = mu + sigma*abs(std)
# threshold = np.full(len(df), 0.4) 

carica, ss, peak, peakpos = signal(df, sigma, a, a_filt, t_ns, threshold, prima=15, dopo=40)

carica_ = np.asanyarray(carica, dtype=object)
peak_= np.asanyarray(peak, dtype=object)

car=np.array(0)
pe=np.array(0)
for i in range(0, len(df)):
    car = np.append(car,carica_[i])
    pe = np.append(pe,peak_[i])
car = np.delete(car, 0)
pe = np.delete(pe, 0)


sat = 1.4  

plt.figure()
plt.title('Amplitude vs Charge 15Dec CRT 250 MHz')
plt.ylabel('Amplitude [V]')
plt.xlabel('Charge [ns.V]')
#plt.plot(car, pe, '.')
plt.plot(car[pe<sat], pe[pe<sat], '.')
plt.plot(car[pe>=sat], pe[pe>=sat], '.r')
plt.grid(True)
#plt.savefig('amp_ch250_C1X', dpi=300)
plt.show()


#%% 62.5 MHz  \Dec15th_CRT  signals & amp vs ch

tnn_ns = tnn * 1e9   #nanosec
inizio_sig=-1e-07*(1e9)   #nanosec
mu62, std62, p62 = norm_fit(df, ann, tnn_ns, inizio_sig)

sigma62 = 4.2 #thresold dev_std sopra baseline x discriminaz segnale
threshold62 = mu62 + sigma62*abs(std62) 
# threshold62 = np.full(len(df), 0.4)

carica62, ss62, peak62, peakpos62 = signal(df, sigma62, ann, tnn_ns, threshold62, prima=15, dopo=40)


carica_62 = np.asanyarray(carica62, dtype=object)
peak_62= np.asanyarray(peak62, dtype=object)

car62=np.array(0)
pe62=np.array(0)
for i in range(0, len(df)):
    car62 = np.append(car62,carica_62[i])
    pe62 = np.append(pe62,peak_62[i])
car62 = np.delete(car62, 0)
pe62 = np.delete(pe62, 0)

sat = 1.4

plt.figure()
#plt.plot(car62, pe62, '.')
plt.title('Amplitude vs Charge 15Dec CRT 62.5 MHz')
plt.ylabel('Amplitude [V]')
plt.xlabel('Charge [ns.V]')
plt.plot(car62[pe62<sat], pe62[pe62<sat], '.')
plt.plot(car62[pe62>=sat], pe62[pe62>=sat], '.r')
plt.grid(True)
#plt.savefig('amp_ch62_C1X', dpi=300)
plt.show()



#%% charge 250 vs 62 data 15th dec

t_ns = t * 1e9   #ns
tnn_ns = tnn * 1e9 #ns

cari, pea = ampl_pe(df, a, t_ns, t_min=int(((1.001111-0.03)*1000)/4), t_max=int(((1.7*1000)/4)))
cari62, pea62 = ampl_pe(df, ann, tnn_ns, t_min=int(((1.001111-0.03)*1000)/16), t_max=int(((1.7*1000)/16)))

popt, pcov = curve_fit(retta, cari62, cari, p0=[1., 0.])   
txt='m= '+str(round(popt[0], 4))+' b= '+str(round(popt[1], 4))
xx = np.linspace(0, 1000, 1000)
plt.figure()
plt.plot(cari62, cari, '.')
plt.plot(xx, retta(xx, *popt), '-', c='r', label = txt)
plt.title('Charge 250 Mhz vs 62.5 Mhz (15Dec CRT)')
plt.ylabel('Charge 250 Mhz [ns.V]')
plt.xlabel('Charge 62 Mhz [ns.V]')
plt.grid()
plt.legend()
plt.show()
# plt.savefig('c250vs62_15dec', dpi=300)

plt.figure()
plt.title('Amplitude vs Charge 250 Mhz (15Dec CRT)')
plt.ylabel('Amplitude [V]')
plt.xlabel('Charge [ns.V]')
plt.plot(cari[pea<sat], pea[pea<sat], '.', label='Non saturated')
plt.plot(cari[pea>=sat], pea[pea>=sat], '.r', label='Saturated')
plt.grid(True)
plt.legend()
plt.show()
# plt.savefig('ac250_15dec', dpi=300)

plt.figure()
plt.title('Amplitude vs Charge 62.5 Mhz (15Dec CRT)')
plt.ylabel('Amplitude [V]')
plt.xlabel('Charge [ns.V]')
plt.plot(cari62[pea62<sat], pea62[pea62<sat], '.', label='Non saturated')
plt.plot(cari62[pea62>=sat], pea62[pea62>=sat], '.r', label='Saturated')
plt.grid(True)
plt.legend()
plt.show()
# plt.savefig('ac62_15dec', dpi=300)





#%% TO DO

# tune finely the integration interval of 16dec data 16 ns 62.5MHz = 0.016 us
#TOT for CRT DATA: thrshold at least 10 mV (tre done to collect data) -> find index of up&down the tre -> deltaT=TOT
# -> plot ampl vs TOT & plot charge vs TOT
# 15dec fare anche con C2X!!!!!!!!!!
#C1XArapuca_Efield_0kV_A4ch2_250MHz_cosmic00000.txt  ALTRI DATI PER FARE TOT (13DEC)

#find peak with filtered data but analyze them with data normal (also moving average filter with shift)
# create data frames
# use my peak finder code to analyze 16dec out of the time interval of signal
# fit histogram of charge (try also 409 bins other than 182)

List Comprehension for signal function homemade
argparse — Parser for command-line options, arguments and sub-commands per passare da linea di comando il path per funzione che legge i dati

df= pd.DataFrame({'data' : [], again})
df.to_csv('data.csv', sep='; or space', index=False)

 CD = np.where(Re < 0, 0.0, CD)  #+ condizioni per un vettore




