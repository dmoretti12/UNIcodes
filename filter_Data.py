# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 16:27:09 2022

@author: Davide
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import lfilter, convolve
from scipy.stats import norm


class DATI:
    
    def __init__(self):
        
        self.a_wout_baseline = 0
        self.a_filt_wout_baseline = 0
        self.a = 0
        self.a_filt = 0
        self.a_dfilt = 0
        self.t = 0
        self.a_buoni = 0
        self.df = 0
        
        self.ns = 0
        
        self.n = 0
        self.b = 0
        self.charge = 0
        
################################################
    
class FILTER_DATA:
    
    def __init__(self, classe_dati):
        

        self.hd = classe_dati
        
        self.l_mov_av = 7
        
        self.threshold_base = 30
        self.bin_base = 100
        self.debug_base = False
        self.exclusion_window = 1000
        self.discard_tre = 800 # number of point minimum for estimating the baseline
        # self.x1 = 0
        # self.x2 = 0
        
        self.inizio_sig = 0
        
        self.do = 'baseline'
        self.filtro = 'mov_av'
        
##################################################################################     
        
    def mov_av(self, a):

        L = self.l_mov_av  # L-point filter
        
        n_wvf = len(a)
        wvf_lenght = len(a[0])
        
        b = (np.ones(L))/L  # numerator co-effs of filter transfer function
        c = np.ones(1)  # denominator co-effs of filter transfer function
        y = np.zeros((n_wvf, wvf_lenght))
        d = a + 2  # shift

        for i in range(0, n_wvf):
            if i==0: print("Doing moving average ..")

            z = convolve(d[i], b)  # filter output using convolution

            z = lfilter(b, c, d[i])  # filter output using lfilter function

            y[i] = z - 2
        
            y[i,0:L-1] = d[i,0:L-1] - 2  # back from shift

        return y  


##################################################################################

    def denoising(self, a, f):

        n_wvf = len(a)
        wvf_lenght = len(a[0])
        signal_filtered = np.zeros((n_wvf, wvf_lenght))

        for i in range(0, n_wvf):

            n = wvf_lenght
            fhat = np.fft.fft(self.hd.a[i], n)  # computes the fft

            psd = fhat * np.conj(fhat)/n
            # freq = (1/(delta*n)) * np.arange(n) #frequency array

            idxs_half = np.arange(1, np.floor(
                n/2), dtype=np.int32)  # first half index
            psd_real = np.abs(psd[idxs_half])  # amplitude for first half

            # Filter out noise
            sort_psd = np.sort(psd_real)[::-1]
            if i==0:
                plt.figure()
                plt.plot(sort_psd)
                plt.show()
            # print(len(sort_psd))
            threshold = sort_psd[f]
            psd_idxs = psd > threshold  # array of 0 and 1
            # psd_clean = psd * psd_idxs #zero out all the unnecessary powers
            fhat_clean = psd_idxs * fhat  # used to retrieve the signal

            signal_filtered[i] = np.fft.ifft(
                fhat_clean)  # inverse fourier transform

        return signal_filtered

#######################################################################################
    
    def center_bins(self, b):
        
        bbin = np.zeros(len(b))
        l = (b[2]-b[1])/2
        bbin = l + b
        b1 = np.delete(bbin, len(b)-1)
        
        return b1    

###################################################

    def sottraggo_baseline(self, a0, a, mu):
        for i in range(len(a)):
            a0[i] = a[i] - mu[i]
        return a0

    def baseline(self, a):#, nwvf, wvf_len):
    
        n_wvf = len(a) 
        wvf_lenght = len(a[0])
        a0 = np.zeros((n_wvf, wvf_lenght))
        mu = np.zeros(n_wvf)
        cont = np.zeros(n_wvf)
        uno = np.ones(n_wvf)
        # valori = a[0:nwvf, 0:wvf_lenght]
        
        
        if self.debug_base==True: 
            plt.figure()                           #PLOT
            n_wvf=2
            
        for i in range (0, n_wvf):
            if i==0: print("Estimating baseline ..")
            if i==int(n_wvf/2): print("Estimated baseline of half of the wvfs ..")
            j=0
            val_acc = np.zeros(wvf_lenght)
            val = a[i] #valori[i]
            # self.x1 = min(val)
            # self.x2 = np.mean(val)+(np.mean(val)-min(val))
            
            if self.debug_base==True:
                plt.hist(val, self.bin_base)#, range=[self.x1, self.x2])         #PLOT
                
            n, b = np.histogram(val, self.bin_base)#,  range=[self.x1, self.x2]) #senza range funziona meglio
            b1 = self.center_bins(b)
            maxx = b1[np.argmax(n)]       # most frequent value
            
            tre = maxx + self.threshold_base
            tre_ = maxx - self.threshold_base

            while j<wvf_lenght:
                if (val[j]<tre and val[j]>tre_):
                    val_acc[j] = val[j]
                    cont[i] += 1
                    j = j + 1
                    
                else: j=j+int(self.exclusion_window/self.hd.ns)

            if cont[i]==0: 
                print(i)
                cont[i]=1
            mu[i] = sum(val_acc)/cont[i]                        #faccio la semplice media perchè è una gaussiana, non serve fare il fit
            if cont[i]<self.discard_tre: uno[i]=0
        # plt.hist(cont, 100, range=[-5, 2500])

        a0 = self.sottraggo_baseline(a0, a, mu)#valori invece di a #salvare in un formato i dati filtrati con baseline
        
        if self.debug_base==True: plt.show()                             #PLOT
        
        return a0, mu, cont, uno#, chi2t
    
###################################################################################

    def norm_fit(self, a, t):
        
        n_wvf = len(a)
        mu = np.empty(n_wvf)
        std = np.empty(self.lenght_wvf)

        for i in range(0, n_wvf):

            # prendo come baseline solo dati prima del segnale
            yy = a[i, t[i] < self.inizio_sig]
            # zz = np.where(np.max(yy) < 0.215 and np.min(yy) > 0.206, yy)
            # print(len(yy), len(zz))

            # mu e std sono esattamente np.mean() e suo errore!!!
            mu[i], std[i] = norm.fit(yy)


        return mu, std

##################################################################################

    def __filter__(self):
        
        if self.do=='norm_fit':
            f = FILTER_DATA.norm_fit(self)
        if self.do=='baseline':
            f = FILTER_DATA.baseline(self)
        if self.do=='mov_av':
            f = FILTER_DATA.mov_av(self)
        if self.do=='denoising':
            f = FILTER_DATA.denoising(self)
        return f
            
    
    
    
    
    
    