# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 16:25:49 2022

@author: Davide
"""

# from pylab import*
from scipy.signal import find_peaks
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit
from scipy.stats import norm

# import scipy.integrate as integrate
# import scipy.special as special

#from sklearn.metrics import mean_squared_error

import os
import glob
import pandas as pd

#from denoising import Denoiser

plt.close("all") #chiude tutte le finestre

os.chdir(r"C:\Users\Davide\Documents\FAMIGLIA\Davide\UNIV\WORK\Dec15th_CRT-20220321T101331Z-001\Dec15th_CRT")

filenames = [i for i in glob.glob("C1X*.txt")]

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
    


#250/62.5
#Out[4]: 4.0

out = int(len(df[0])/4)

tnn = np.empty((len(df),out))  #freq 62.5 Mhz
ann = np.empty((len(df),out))
    
for i in range (0, out):
    
    tnn[:,i] = t[:,4*i]
    ann[:,i] = a[:,4*i]
    
    #t[:,i] = np.full((len(filenames),len(df[0])),T, order='C')
    #a = np.append(a,A)

#t = np.delete(t, 0)
#a = np.delete(a, 0)
# data = np.genfromtxt('data CRT15dec\C1XArapuca_CRT_250MHz_00001.txt', skip_header=6)#, unpack=True)

# t=data[:,0]
# a=data[:,1]

#%% DENOISING


# denoiser = Denoiser()

# denoised = denoiser.denoise(a[0], 1250)

# fig, ax = plt.subplots()
# ax.plot(t[0], a[0])
# ax.plot(t[0], denoised)
# plt.show()



n = len(t[19])
fhat = np.fft.fft(a[19], n) #computes the fft
psd = fhat * np.conj(fhat)/n
#freq = (1/(delta*n)) * np.arange(n) #frequency array
idxs_half = np.arange(1, np.floor(n/2), dtype=np.int32) #first half index
psd_real = np.abs(psd[idxs_half]) #amplitude for first half


## Filter out noise
sort_psd = np.sort(psd_real)[::-1]
# print(len(sort_psd))
threshold = sort_psd[100]
psd_idxs = psd > threshold #array of 0 and 1
psd_clean = psd * psd_idxs #zero out all the unnecessary powers
fhat_clean = psd_idxs * fhat #used to retrieve the signal

signal_filtered = np.fft.ifft(fhat_clean) #inverse fourier transform


## Visualization
# fig, ax = plt.subplots(2,1)
# ax[0].plot(t[0], a[0], color='b', label='Noisy Signal')
# ax[0].set_xlabel('t axis')
# ax[0].set_ylabel('Accn in Gal')
# ax[0].legend()

# ax[1].plot(freq[idxs_half], np.abs(psd[idxs_half]), color='b', lw=0.5, label='PSD noisy')
# ax[1].set_xlabel('Frequencies in Hz')
# ax[1].set_ylabel('Amplitude')
# ax[1].legend()

# ax[2].plot(freq[idxs_half], np.abs(psd_clean[idxs_half]), color='r', lw=1, label='PSD clean')
# ax[2].set_xlabel('Frequencies in Hz')
# ax[2].set_ylabel('Amplitude')
# ax[2].legend()

# ax[1].plot(t[0], signal_filtered, color='r', lw=1, label='Clean Signal Retrieved')
# ax[1].set_ylim([minsignal, maxsignal])
# ax[1].set_xlabel('t axis')
# ax[1].set_ylabel('Accn in Gal')
# ax[1].legend()

# plt.subplots_adjust(hspace=0.6)
# plt.savefig('real-signal-analysis.png', bbox_inches='tight', dpi=300)

fig, ax = plt.subplots()
ax.plot(t[19], a[19])
ax.plot(t[19], signal_filtered)
plt.show()



#%%%   HISTOGRAMMM

def histo_baseline(df, a, t, nb):

    bl = np.empty((len(df),0))
    ni = np.empty((len(df),0))


    
    for i in range(0, len(df)):
        
        
        yy = a[i, t[i] < -1e-07]   #prendo come baseline solo dati prima del segnale
        if i==0:
            bl = np.resize(bl, (len(df),len(yy)))
            tt = np.resize(t, (len(df),len(yy)))
        
        bl[i,:] = yy  #limite sull'ampiezza nel considerare baseline per non avere outliers che spostano la media
        tt[i] = t[i, t[i] < -1e-07]
        
        binn = np.linspace(np.min(bl[i]), np.max(bl[i]), nb)#np.histogram_bin_edges(bl)
        
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
    
        # bx = fig.add_subplot(111)
        # bx.plot(tt[i], bl[i],'.')
        
        #n, b, p = ax.hist(bl[i], binn, density=True, color='b', histtype='step') 
        n, b = np.histogram(bl[i], binn, density=True)#, color='b', histtype='step') 
        
        # if np.any(n == 0) == True:
        #     binn = np.linspace(np.min(bl[i]), np.max(bl[i]), 8)
        #     if np.any(binn == 0) == True:   
        #         binn = np.linspace(np.min(bl[i]), np.max(bl[i]), 7)
        #     n, b = np.histogram(bl[i], binn, density=True)
        
        if i==0: ni = np.resize(ni, (len(df),len(n)))
        ni[i] = n
        

        bbin = np.zeros(len(n))   # x of histogram as center of every bin for the fit
        l = (b[2]-b[1])/2 
        bbin = l + np.resize(b, len(n))
    
    return ni, bbin, bl, binn

#%% RMSE instead of histogram and dev_std


def rmse_function(df, a, t):
    
    rms = np.empty(len(df))
    rmse = np.empty(len(df))
    mse = np.empty(len(df))
    meann = np.empty(len(df))
    rmse_mean = np.empty(len(df))
    
    for i in range(0, len(df)):
        
        yy = a[i, t[i] < -1e-07]   #prendo come baseline solo dati prima del segnale
        
        meann[i] = np.mean(yy)  #è uguale a mu di norm.fit!!!!
        
        yy = yy + 2  #offset x usare rms (intorno a 0 non funziona)
        rms[i] = math.sqrt(np.sum(yy**2)/len(yy))
        rms[i] = rms[i] - 2  #rimuovo l'offset
        yy = yy - 2
        
        mse[i] = np.square(np.subtract(yy,rms[i])).mean()   # e-17 risp a std di norm.fit!!!
        rmse_mean[i] = math.sqrt( np.square(np.subtract(yy,meann[i])).mean() )
        rmse[i] = math.sqrt(mse[i])
        #rmse[i] = mean_squared_error(yy, rms, squared=True) #[rms for _ in yy]
        #mse[i] = mean_squared_error(yy, rms, squared=False)
        
    return meann, rms, rmse, rmse_mean, mse


#%%  FIT CON NORM.FIT

def norm_fit(df, a, t):
    
    mu = np.empty(len(df))
    std = np.empty(len(df))
    pp = []
    
    for i in range(0, len(df)):
        
        yy = a[i, t[i] < -1e-07]   #prendo come baseline solo dati prima del segnale
        
        mu[i], std[i] = norm.fit(yy)    #mu e std sono esattamente np.mean() e suo errore!!!
        
        #xmin, xmax = plt.xlim()
        
        x = np.linspace(-0.03, 0.03, 100)
        
        p = norm.pdf(x, mu[i], std[i]) #best_fit_line
        pp.append(p)
        
        # plt.figure()
        # plt.plot(x, p, 'k', linewidth=2)   #per vedere il fit p=norm.pdf(np.linspace(-0.03, 0.03), mu[i], std[i])
        # plt.show()
        
    return mu, std, p


#%%  HISTOGRAM 1D   

# bl = np.array(len(df[0]))
# bl = bl[t < -1e-07]  #prendo come baseline solo dati prima del segnale
# plt.figure() #apre l'ambiente grafico
# plt.plot(tt[tt<-1e-07], bl, 'o')
# plt.grid(True)  #aggiunge la griglia
# plt.xlabel('Time [s]')
# plt.ylabel('Amplitude [V]')
# plt.title("Baseline histogram (250 Hz)")
# plt.show()
# binn = np.linspace(np.min(bl), np.max(bl), 17)#np.histogram_bin_edges(bl)
# isto = np.histogram(bl, binn, density=True, color='r', histtype='step')




#%%   FIIIITTTTTT con curve_fit

def fit_baseline(df, ni, bbin):
    
    def gauss(x, a, x0, sigma):
        return a*np.exp(-0.5*((x-x0)/sigma)**2)
    
    
    
    popti = np.zeros((len(df),3))
    popti[0] = [100, 0, 0.02]
    #pcovi = np.zeros((len(df),3))
    ym = np.empty((len(df),50))  #why 50???? saw it
    
    for j in range (0, len(df)):
        
        #if j==0:
        popt, pcov = curve_fit(gauss, bbin, ni[j], popti[0], maxfev=10000)
       # else:
          #  popt, pcov = curve_fit(gauss, bbin, ni[j], popti[j-1], maxfev=20000)#, method='dogbox')
        
        
        popti[j] = popt
        #pcov[i] = pcov
    
        ym[j] = gauss(np.linspace(-0.03, 0.03), popt[0], popt[1], popt[2])
    
        #ax.plot(np.linspace(-0.03, 0.03), ym[j], c='r', label='Best fit', linewidth = 2, linestyle='--')
        #ax.legend()
        
    return popti, ym




#%%%   DEV STD TO DISCRIMINIZE THE SIGNALS + INTEGRAL

#plt.close("all")
def signal(df, sigma, a, t, mu, std):
    
    
    dev_std = np.empty(len(df))
    dev_std = std    #popti[:,2]
    mean = np.empty(len(df))
    mean = mu     #popti[:,1]
    
   # info = np.empty((len(df),5))  #array con info x ogni waveform. 1) peak; 2)peak pos; 3)charge;
    
    carica = []
    indsign = []
    peak = []
    peakpos = []
    ss = []
    
    
    for i in range(0, len(df)):
        
        print(i)
        if i==529 or i==542 or i==19 or i==526:   #segnali brutti: CX1 15th Dec
            i=i+1
        tre = mean[i] + sigma*abs(dev_std[i])   #THRESHOLD
        j=0
        c = []   # array integrale segnale (charge) della i-esima waveform
        s = []   # indici segnali
        p = []   # peaks
        pt = []  # peaks position
        soprasotto = []
        
        while j < len(a[i]):
            
            if j>=len(a[i])-2:
                break
            
            if a[i, j] > tre: # and a[i, j+1] > tre:  #and a[i, j+2] > tre:   #if at least 3 points are above the 3-4*dev_std i consider it as a peak
               # s.append(j)
                    
                for k in range(j, len(a[i])-j):

                    if a[i,k] < tre:
                        #s.append(k)
                        break   
                    if k>= len(a[i]):
                        break
              
                if j < k:
                    
                    integral = np.trapz(a[i,j-15:k+40], t[i,j-15:k+40]) #faccio l'integrale prendendo qualche pto prima e qualcuno dopo
                    c.append(integral)
                    s.append([j,k])  #indici sopra e sotto il threshold
    
                    picco = a[i,j:k]
                    p.append(np.max(picco))
                    pt.append(np.argmax(picco))
                
                # plt.figure()
                # plt.plot(t[i,j-2:k+2], a[i,j-2:k+2], '.')
                # plt.show()
                soprasotto.append([j,k])
                
                j = k + int(0.1*len(t)) #in modo che non rientri nel for per i punti appartenenti al picco + evita afterpulse
                
            else: j = j + 1
            
        carica.append(c)
        indsign.append(s)
        peak.append(p)
        peakpos.append(pt)
        ss.append(soprasotto)
        
    return carica, indsign, peak, peakpos, ss
    
    # info[i,0] = *peak[i]
    # info[i,1] = *peakpos[i]
    # info[i,2] = *carica[i]
    # info[i,3] = *indsign[i]
    
    
    
#%% PEAK FINDER

def peak_finder(df, sigma, a, t, mu, std):
    
    dev_std = np.empty(len(df))
    dev_std = std
    mean = np.empty(len(df))
    mean = mu
    
    carica = []
    pospicchi = []
    propr = []
    
    for i in range(0, len(df)):
        
        tre = mean[i] + sigma*abs(dev_std[i])   #THRESHOLD
        
        peaks, properties = find_peaks(a[i], height=None, threshold=tre, distance=None, prominence=None, width=None, wlen=None, rel_height=0.5, plateau_size=None)
    
        pospicchi.append(peaks)
        propr.append(properties)
        
    return pospicchi, propr

#%%  DATA FROM 16th DEC mini arapuca

import numpy as np
import matplotlib.pyplot as plt  

a16 = np.genfromtxt('C2MiniArapuca_A4ch1_250MHz_LED_width20ns_ampl25V_a_00000.txt', skip_header=3)
points = 938
wvf = int(len(a16)/points) # =10000
aa = np.empty((wvf, points))
tt = np.empty((wvf, points))

for i in range (0, wvf):
        
        aa[i] = a16[points*i:points*(i+1)]
        tt[i] = np.arange(0, 4e-03*points, 4e-03)   #250 MHz = 4e-03 us
        

#freq 62.5 Mhz
out = int(len(aa[0])/4)

tt62 = np.empty((len(aa),out))  
aa62 = np.empty((len(aa),out))
    
for i in range (0, out):
    
    tt62[:,i] = tt[:,4*i]
    aa62[:,i] = aa[:,4*i]
    
    
# istogramma carica da 2.7 a 2.9 forse meno   sulle x la carica sulle y i conteggi (dunque devi avere 10000 entrate)
# decidi bene lintervallo
#rifai la stessa cosa, insto carica, per dati a 62.5 MHz
# altezza (y) vs carica (x)  15 dec CRT
# TOT

# tmin = 2.73
# tmax=2.86
# charge = np.empty(len(aa)) 

# for i in range(0, len(aa)):

#     t_int = tt[i, tt[i] > tmin]
#     a_int = aa[i, tt[i] > tmin]
#     yy = a_int[t_int < tmax]
#     zz = t_int[t_int < tmax]

#     integral = np.trapz(yy, zz)
#     charge[i] = integral



# n, b, p = plt.hist(charge, 180, density=True)

# tmin = 2.73
# tmax=2.86
# charge62 = np.empty(len(aa62)) 

# for i in range(0, len(aa62)):

#     t_int = tt62[i, tt62[i] > tmin]
#     a_int = aa62[i, tt62[i] > tmin]
#     yy = a_int[t_int < tmax]
#     zz = t_int[t_int < tmax]

#     integral = np.trapz(yy, zz)
#     charge62[i] = integral



# n62, b62, p = plt.hist(charge62, 220, density=True)
#%% istogramma carica fotopicchi

def isto_carica(df, aa, tt, tmin, tmax, _bin):
    
    charge = np.empty(len(df)) 
    
    for i in range(0, len(df)):
        
        t_int = tt[i, tt[i] > tmin]
        a_int = aa[i, tt[i] > tmin]
        yy = a_int[t_int < tmax]
        zz = t_int[t_int < tmax]
        
        integral = np.trapz(yy, zz)
        charge[i] = integral
        
        
            
    n, b = np.histogram(charge, _bin, density=True)
    
    return n, b, charge
    
#%% dati 16 dec miniarapuca

n, b, charge = isto_carica(aa, aa, tt, 2.73, 2.86, _bin=182)
plt.figure()
plt.hist(charge, b, density=False)
plt.title('Charge histogram (miniarapuca16thDecLED) 250 MHz')
plt.ylabel('Counts')
plt.xlabel('Charge [us.V]')

n62, b62, charge62 = isto_carica(aa62, aa62, tt62, 2.73, 2.86, _bin=220)
plt.figure()
plt.hist(charge62, b62, density=False)
plt.title('Charge histogram (miniarapuca16thDecLED) 62.5 MHz')
plt.ylabel('Counts')
plt.xlabel('Charge [us.V]')
    
#%%  250 MHz arrays a & t

#ni, bbin, bl, binn = histo_baseline(df, a, t, nb=9)  #nb = numero bin (look at data)

#meann, rms, rmse, rmse_mean, mse = rmse_function(df, a, t)

mu, std, p = norm_fit(df, a, t)

#popti, ym = fit_baseline(df, ni, bbin)

sigma = 4 #thresold dev_std sopra baseline x discriminaz segnale

carica, indsign, peak, peakpos, ss = signal(df, sigma, a, t, mu, std)

carica_ = np.asanyarray(carica, dtype=object)
peak_= np.asanyarray(peak, dtype=object)

car=np.array(0)
pe=np.array(0)
for i in range(0, len(df)):
    car = np.append(car,carica_[i])
    pe = np.append(pe,peak_[i])


plt.figure()
plt.plot(car, pe, '.')
plt.title('Amplitude vs Charge 15Dec CRT 250 MHz')
plt.ylabel('Amplitude [V]')
plt.xlabel('Charge [s.V]')


#pospicchi, propr = peak_finder(df, sigma, a, t, mu, std)

#RETURN PROPERTIES PEAK_FINDER

# ‘peak_heights’;; ‘left_thresholds’, ‘right_thresholds’;; ‘prominences’, ‘right_bases’, ‘left_bases’
#‘width_heights’, ‘left_ips’, ‘right_ips’;; ‘plateau_sizes’, left_edges’, ‘right_edges’

#%% 62.5 MHz arrays ann & tnn

ninn, bbinnn, bl, binnn = histo_baseline(df, ann, tnn, nb=6)  #nb = numero bin (look at data)

meannn, rmsnn, rmsenn, rmse_meannn, msenn = rmse_function(df, ann, tnn)

munn, stdnn, pnn = norm_fit(df, ann, tnn)

#poptinn, ymnn = fit_baseline(df, ninn, bbinnn)

sigmann = 4 #thresold dev_std sopra baseline x discriminaz segnale

caricann, indsignnn, peaknn, peakposnn, ssnn = signal(df, sigmann, ann, tnn, munn, stdnn)

#%%

# poptinn[286]
# Out[177]: array([ 3.12697760e+01, -1.02441033e+06, -2.96513564e+07])

# poptinn[285]
# Out[178]: array([ 3.90729755e+01, -9.98613741e+02, -4.59153458e+04])

# poptinn[284]
# Out[179]: array([    39.0739543 ,   -888.30283104, -38630.0083935 ])

# poptinn[283]
# Out[180]: array([5.20857528e+01, 9.57780309e-03, 3.51020838e+00])

# poptinn[282]
# Out[181]: array([ 1.16581779e+02, -2.54819155e-03,  2.80073896e-03])


# RMS DELL'ARRAY BASELINE' 250 e confronta con dev_std isrtoggramma
# poi falla per 62.5m

# plot cariche campionamenti diversi degli stessi picchi (dovrebbe essere lineari) [dipsoersione]

https://www.journaldev.com/16160/python-pip
https://www.journaldev.com/18341/python-scikit-learn-tutorial

NORM.FIT USA COME MU E STD ESATTAMENTE NP.MEAN() E SUO ERRORE, POI FA IL FIT
IN QUALCHE MODO CON SOLI QUESTI DUE PARAMETRI, QUESTO SPIEGHEREBBE PERCHE
L ALTEZZA FA SCHIFO

#scipy.signal.find_peaks(x, height=None, threshold=None, distance=None, prominence=None, width=None, wlen=None, rel_height=0.5, plateau_size=None)

plt.figure()
ss[h]
plt.plot(t[h], a[h], '-')
plt.plot(t[h, (307-15):(308+40)], a[h, (307-15):(40+308)], '.')