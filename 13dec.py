# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 12:05:33 2022

@author: Davide
"""
from functions import read_many_txt, read_long_txt, denoising, norm_fit, signal, isto_carica, tot, ampl_pe, mov_av, multi_fit, gauss, multi_gauss
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
# import random


plt.close("all") #chiude tutte le finestre

def expo(x, a, b):
    return np.exp(a*x + b)

def retta(x, m, q):
    return m*x + q

#  DATA 13th DEC 

dec13 = 'C1XArapuca_Efield_0kV_A4ch2_250MHz_cosmic00000.txt'
righe_iniz = 3
puntixwvf = 1252

a13, t13, a13_62, t13_62 = read_long_txt(dec13, righe_iniz, puntixwvf)

#%% charge 250 vs 62 data 13 dec

t13_ns = t13 * 1000   #ns
t13_62_ns = t13_62 * 1000 #ns

c13, p13 = ampl_pe(a13, a13, t13_ns, t_min=int((0.8*1000)/4), t_max=int((2*1000)/4))
c13_62, p13_62 = ampl_pe(a13_62, a13_62, t13_62_ns, t_min=int((0.8*1000)/16), t_max=int(((2*1000)/16)))

popt, pcov = curve_fit(retta, c13_62, c13, p0=[1., 0.])   
txt='m= '+str(round(popt[0], 4))+' b= '+str(round(popt[1], 4))
xx = np.linspace(0, 1000/0.3, 1000)

plt.figure()
plt.plot(c13_62/0.3, c13/0.3, '.')
plt.plot(xx, retta(xx, *popt), '-', c='r', label = txt)
plt.title('P.E. 250 Mhz vs 62.5 Mhz (13Dec-C1XArapuca_Efield_0kV_A4ch2_cosmic)')
plt.ylabel('Photoelectrons 250 MHz')
plt.xlabel('Photoelectrons 62.5 MHz')
plt.grid()
plt.legend()
plt.show()
# plt.savefig('c250vs62_13dec_pe', dpi=300)

#%% TOT 13DEC 250 MHz

tre_tot13 = 0.4
#tre_tot = 0.3
#c_tot, ss_tot, peak_tot, peakpos_tot = signal(df, 1, a, t, tre_tot, 15, 40)
#_, _, tot_c, tot_peak = isto_carica(df, a, t, t_min=int((-0.3e-01*1000)/4), t_max=int((1.3*1000)/4), _bin=10) 
#bisogna capire quali wvf sono prese, quelle sopra la threshold
t13_ns = t13 * 1000  #us -> ns
ctot_13, ss_tot13, peaktot_13, peakpos13 = tot(a13, a13, t13_ns, t_min=int(0.8*1000/4), t_max=int(2*1000/4), threshold=tre_tot13, satur=1.4) 

# sstot_13 = np.array(0)
# peaktot_13 = np.array(0)
# deltaT13 = np.array(0)
# ctot_13 = np.array(0)

# for i in range(0, len(c_tot13)):
#     sstot_13 = np.append(sstot_13, ss_tot13[i])
#     peaktot_13 = np.append(peaktot_13, peak13_tot[i])
#     ctot_13 = np.append(ctot_13, c_tot13[i])
# sstot_13 = np.delete(sstot_13, 0)
# peaktot_13 = np.delete(peaktot_13, 0)
# ctot_13 = np.delete(ctot_13, 0)
# i=0
# while i < len(sstot_13)-1:
#         # if sstot_[i]==0 and sstot_[i+1]==0:
#         #     deltaT = np.append(deltaT, 0)
#         # else:
#         deltaT13 = np.append(deltaT13, sstot_13[i+1]-sstot_13[i])
#         i = i+2
# deltaT13 = np.delete(deltaT13, 0)
# deltaT13 = deltaT13 * 4 #in nanosecondi       / 1e3   # in microsecondi

deltaT13 = np.zeros(len(a13))
for i in range (0, len(a13)):
    deltaT13[i] = (ss_tot13[i,1] - ss_tot13[i,0])

satur = 1.4
      
plt.figure()
plt.semilogy()
plt.plot(deltaT13[peaktot_13<satur], peaktot_13[peaktot_13<satur], '.b', label='Non saturated')
plt.plot(deltaT13[peaktot_13>=satur], peaktot_13[peaktot_13>=satur], '.r', label='Saturated')
plt.title('Amplitude vs TOT 13Dec Cathode 250 MHz')
plt.ylabel('Amplitude [V]')
plt.xlabel('TOT [ns]')
plt.legend()
plt.grid(True, which='both')
# plt.savefig('amp_tot250_C1X-13dec', dpi=300)
plt.show()

pu=deltaT13[peaktot_13<satur]
po=ctot_13[peaktot_13<satur]
q=pu[pu<250]
qq=po[pu<250]
popt, pcov = curve_fit(expo, q[q>70], qq[q>70], p0=[0.007, 3.8])
fittxt='m= '+str(round(popt[0], 4))+' b= '+str(round(popt[1], 4))

dev = math.sqrt(sum((qq[q>70]-expo(q[q>70], *popt))**2)/(len(q[q>70])-1))
#print('dev', dev)
dev_rel = dev / 300.

ret = deltaT13[peaktot_13>=satur]
yret = expo(ret, *popt)
for i in range(0, len(ret)):
        yret[i] = np.random.normal(yret[i], dev_rel*yret[i])

plt.figure()
plt.semilogy()
plt.plot(deltaT13[peaktot_13<satur], ctot_13[peaktot_13<satur], '.b', label='Non saturated')
plt.plot(deltaT13[deltaT13<300], expo(deltaT13[deltaT13<300], *popt), 'g', label=fittxt)
plt.plot(deltaT13[peaktot_13>=satur], ctot_13[peaktot_13>=satur], '.r', label='Saturated')
plt.plot(deltaT13[peaktot_13>=satur], yret, '.g', label = 'Expected')
plt.title('Charge vs TOT 13Dec Cathode 250 MHz')
plt.ylabel('Charge [V.ns]')
plt.xlabel('TOT [ns]')
plt.legend()
plt.grid(True, which='both')
# plt.savefig('ch_tot250_C1X-13dec', dpi=300)
plt.show()    


#%% TOT 13DEC 62.5 MHz
from functions import tot

tre_tot13_62 = 0.4
#tre_tot = 0.3
#c_tot, ss_tot, peak_tot, peakpos_tot = signal(df, 1, a, t, tre_tot, 15, 40)
#_, _, tot_c, tot_peak = isto_carica(df, a, t, t_min=int((-0.3e-01*1000)/4), t_max=int((1.3*1000)/4), _bin=10) 
#bisogna capire quali wvf sono prese, quelle sopra la threshold
t13_62_ns = t13_62 * 1000  #us -> ns
ctot_13_62, ss_tot13_62, peaktot_13_62, peakpos62 = tot(a13_62, a13_62, t13_62_ns, t_min=int(0.8*1000/16), t_max=int(2*1000/16), threshold=tre_tot13_62, satur=1.4) 
# sstot_13_62 = np.array(0)
# peaktot_13_62 = np.array(0)
# deltaT13_62 = np.array(0)
# ctot_13_62 = np.array(0)

# for i in range(0, len(c_tot13_62)):
#     sstot_13_62 = np.append(sstot_13_62, ss_tot13_62[i])
#     peaktot_13_62 = np.append(peaktot_13_62, peak13_tot_62[i])
#     ctot_13_62 = np.append(ctot_13_62, c_tot13_62[i])
# sstot_13_62 = np.delete(sstot_13_62, 0)
# peaktot_13_62 = np.delete(peaktot_13_62, 0)
# ctot_13_62 = np.delete(ctot_13_62, 0)
# i=0
# while i < len(sstot_13_62)-1:
#         # if sstot_[i]==0 and sstot_[i+1]==0:
#         #     deltaT = np.append(deltaT, 0)
#         # else:
#         deltaT13_62 = np.append(deltaT13_62, sstot_13_62[i+1]-sstot_13_62[i])
#         i = i+2
# deltaT13_62 = np.delete(deltaT13_62, 0)
# deltaT13_62 = deltaT13_62 * 4 #in nanosecondi       / 1e3   # in microsecondi
deltaT13_62 = np.zeros(len(a13_62))
for i in range (0, len(a13_62)):
    deltaT13_62[i] = (ss_tot13_62[i,1] - ss_tot13_62[i,0])
satur = 1.4
      
plt.figure()
plt.semilogy()
plt.plot(deltaT13_62[peaktot_13_62<satur], peaktot_13_62[peaktot_13_62<satur], '.b', label='Non saturated')
plt.plot(deltaT13_62[peaktot_13_62>=satur], peaktot_13_62[peaktot_13_62>=satur], '.r', label='Saturated')
plt.title('Amplitude vs TOT 13Dec Cathode 62.5 MHz')
plt.ylabel('Amplitude [V]')
plt.xlabel('TOT [ns]')
plt.legend()
plt.grid(True, which='both')
# plt.savefig('test', dpi=300)
plt.show()

pu=deltaT13_62[peaktot_13_62<satur]
po=ctot_13_62[peaktot_13_62<satur]
q=pu[pu<250]
qq=po[pu<250]
popt, pcov = curve_fit(expo, q[q>70], qq[q>70], p0=[0.007, 3.8])
fittxt='m= '+str(round(popt[0], 4))+' b= '+str(round(popt[1], 4))

dev = math.sqrt(sum((qq[q>70]-expo(q[q>70], *popt))**2)/(len(q[q>70])-1))
#print('dev', dev)
dev_rel = dev / 300.

ret = deltaT13_62[peaktot_13_62>=satur]
yret = expo(ret, *popt)
for i in range(0, len(ret)):
        yret[i] = np.random.normal(yret[i], dev_rel*yret[i])

    

plt.figure()
plt.semilogy()
plt.plot(deltaT13_62[peaktot_13_62<satur], ctot_13_62[peaktot_13_62<satur], '.b', label='Non saturated')
plt.plot(deltaT13_62[deltaT13_62<300], expo(deltaT13_62[deltaT13_62<300], *popt), 'g', label=fittxt)
plt.plot(deltaT13_62[peaktot_13_62>=satur], ctot_13_62[peaktot_13_62>=satur], '.r', label='Saturated')
plt.plot(deltaT13_62[peaktot_13_62>=satur], yret, '.g', label = 'Expected')
plt.title('Charge vs TOT 13Dec Cathode 62.5 MHz')
plt.ylabel('Charge [V.ns]')
plt.xlabel('TOT [ns]')
plt.legend()
plt.grid(True, which='both')
# plt.savefig('tes', dpi=300)
plt.show()   

### time res
# pp62 = np.zeros(len(wvf62))
# for i in range(0, len(wvf62)):
#     pp62[i] = peakpos62[int(wvf62[i])]
    
# bb = np.arange(950, 1100, 10)
# plt.figure()
# plt.hist(pp62, bins=bb)
# plt.show()

#%%  250 MHz ampl vs charge  \13DEC

# a_filt13 = mov_av(a13, a13)
t13_ns = t13 * 1e3  #s -> ns
# inizio_sig=0.7*1e3   #nanosec
# mu, std, p = norm_fit(a13, a13, t13_ns, inizio_sig)

# sigma = 4 #thresold dev_std sopra baseline x discriminaz segnale
# threshold = mu + sigma*abs(std)
# threshold = np.full(len(a13), 0.11) 

# carica, ss, peak, peakpos = signal(a13, sigma, a13, a13, t13_ns, threshold, prima=15, dopo=40)
# carica, ss, peak, peakpos = tot(a13, a13, t13_ns, t_min=int(0.8*1000/4), t_max=int(1.5*1000/4), threshold=threshold)
carica, peak = ampl_pe(a13, a13, t13_ns, t_min=int(0.8*1000/4), t_max=int(1.5*1000/4))

# carica_ = np.asanyarray(carica, dtype=object)
# peak_= np.asanyarray(peak, dtype=object)

# car=np.array(0)
# pe=np.array(0)
# for i in range(0, len(a13)):
#     car = np.append(car,carica_[i])
#     pe = np.append(pe,peak_[i])
# car = np.delete(car, 0)
# pe = np.delete(pe, 0)


sat = 1.4  

plt.figure()
plt.title('Amplitude vs P.E. 13Dec 250 MHz')
plt.ylabel('Amplitude [V]')
plt.xlabel('Photo-electrons')   #('Charge [ns.V]')
#plt.plot(car, pe, '.')
plt.plot(carica[peak<sat]/0.3, peak[peak<sat], '.')   #divide axis x /0.3 and is Photo-electrons
plt.plot(carica[peak>=sat]/0.3, peak[peak>=sat], '.r')
plt.grid(True)
# plt.savefig('amp_ch250_C1X', dpi=300)
plt.show()


#%% 62.5 MHz  13DEC  amp vs ch

t13_62_ns = t13_62 * 1e3   #nanosec
# inizio_sig=0.7*1e3   #nanosec
# mu62, std62, p62 = norm_fit(a13, a13_62, t13_62_ns, inizio_sig)

# sigma62 = 4.2 #thresold dev_std sopra baseline x discriminaz segnale
# threshold62 = mu62 + sigma62*abs(std62) 
# threshold62 = np.full(len(a13), 0.11)

#carica62, ss62, peak62, peakpos62 = signal(a13, sigma62, a13_62, t13_62_ns, threshold62, prima=15, dopo=40)
# carica62, ss62, peak62, peakpos62 = tot(a13, a13_62, t13_62_ns, t_min=int(0.8*1000/16), t_max=int(1.8*1000/16), threshold=threshold62)
car62, pe62 = ampl_pe(a13, a13_62, t13_62_ns, t_min=int(0.8*1000/16), t_max=int(1.5*1000/16))

# carica_62 = np.asanyarray(carica62, dtype=object)
# peak_62= np.asanyarray(peak62, dtype=object)

# car62=np.array(0)
# pe62=np.array(0)
# for i in range(0, len(a13)):
#     car62 = np.append(car62,carica_62[i])
#     pe62 = np.append(pe62,peak_62[i])
# car62 = np.delete(car62, 0)
# pe62 = np.delete(pe62, 0)

sat = 1.4

plt.figure()
#plt.plot(car62, pe62, '.')
plt.title('Amplitude vs P.E. 13Dec 62.5 MHz')
plt.ylabel('Amplitude [V]')
plt.xlabel('Photo-electrons')       #('Charge [ns.V]')
plt.plot(car62[pe62<sat]/0.3, pe62[pe62<sat], '.')    #divide axis x /0.3 and is Photo-electrons
plt.plot(car62[pe62>=sat]/0.3, pe62[pe62>=sat], '.r')
plt.grid(True)
#plt.savefig('amp_ch62_C1X', dpi=300)
plt.show()
