# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 17:26:17 2022

@author: Davide
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
plt.rcParams.update({'font.size': 18})
from functions import center_bins

class ANALYZE_DATA:
    
    def __init__(self, classe_dati):
        
        self.hd = classe_dati
        
        self.signal_threshold = 37
        self.prima_t = 30 
        self.dopo_t = 200
        self.time_after_pulse = 400   #o 375
        self.lim_amp = 450
        self.debug_signal = False
        self.windows = 6
        self.number_events_expected = 50
        self.t_min_sig = 0
        self.t_max_sig = 0
        self.sig_plot = True
        
        self.selection=False
        self.isto_plot = True
        self.t_min = 0
        self.t_max = 0
        self.range_x1 = 0
        self.range_x2 = 0
        self.bin = 100
        
        self.p0 = [160, 0, 50, 120, 110, 30, 80, 50]
        self.numg = 2
        self.bound1 = -100
        self.bound2 = 10000
        self.fitdue = False
        self.optimiz = False
        self.start1 = 0
        self.start2 = 0
        self.fin1 = 0
        self.fin2 = 0
        
        self.x1_clean = int(self.t_min/self.hd.ns) - 10
        self.x2_clean = int(self.t_max/self.hd.ns) + 10
        self.lim_inf_clean = -3
        
######################################################################################       
    
    def signal(self, t, a, af):

        n_wvf = len(a)
        wvf_lenght = len(a[0])
        ind_min_sig = int(self.t_min_sig/self.hd.ns)
        ind_max_sig = int(self.t_max_sig/self.hd.ns)
        prima = int(self.prima_t/self.hd.ns)
        dopo = int(self.dopo_t/self.hd.ns)
        
        carica = np.zeros((n_wvf, self.number_events_expected))
        # sop = np.zeros((n_wvf, self.number_events_expected))
        # sot = np.zeros((n_wvf, self.number_events_expected))
        # peak = np.zeros((n_wvf, self.number_events_expected))
        # peakpos = np.zeros((n_wvf, self.number_events_expected))
        e = np.zeros(n_wvf)
        
        if self.debug_signal==True: n_wvf=self.windows
        
        for i in range(0, n_wvf):
            
            if i==0: print("Reading signals ..")
            if i==int(n_wvf/2): print("Read half of the wvfs ..")
            val = a[i]
            val_filt = af[i]
            time = t[i]
            tre = self.signal_threshold # mean[i] + sigma*abs(dev_std[i])   #THRESHOLD
            j = 0
            ev = 0
            trigger = False
            
            if self.debug_signal==True:
                plt.figure()                        #PLOT
                ax=plt.gca()
                ax.tick_params(bottom=True, top=True, left=True, right=True)
                ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
                plt.plot(time, val, '-')           #PLOT
                plt.plot(time, val_filt, '-r')     #PLOT

            j=25     #SKIPPING THE FIRST 100 NS UNTIL I FIXED FILTER
            while j < wvf_lenght:
                
                if j >= wvf_lenght-2:
                    break
                
                if (val_filt[j] > tre and trigger == False):

                    jj = j
                    trigger = True
                    j = j + 15

                elif (val_filt[j] < tre and trigger == True):

                    k = j
                    trigger = False
                    
                    picco = val[jj:k]
                    # peak[i, ev] = np.max(picco)
                    if np.max(picco)>self.lim_amp: continue    #for data random wall leakage

                    ind_picco = np.argmax(picco)+jj
                    if ind_min_sig!=0 and ind_max_sig!=0:
                        if ind_picco>ind_min_sig and ind_picco<ind_max_sig: continue

                    # peakpos[i, ev] = ind_picco
                    # faccio l'integrale prendendo qualche pto prima e qualcuno dopo
                    integral = np.trapz(val[ind_picco-prima:ind_picco+dopo],
                                        time[ind_picco-prima:ind_picco+dopo])
                    carica[i, ev] = integral

                    # indici sopra e sotto il threshold
                    # sop[i, ev] = jj
                    # sot[i, ev] = k
                    
                    if self.debug_signal==True:
                        plt.plot(time[ind_picco], np.max(picco), 'x', markersize=18)              #PLOT
                        plt.plot(time[ind_picco-prima:ind_picco+dopo], val[ind_picco-prima:ind_picco+dopo], '-')      #PLOT
                    
                    # in modo che non rientri nel for per i punti appartenenti al picco + evita afterpulse
                    j = k + int(self.time_after_pulse/self.hd.ns)
                    ev += 1
                    e[i] = int(ev)
                else: j += 1

        if self.debug_signal==True:
            plt.ylabel('Amplitude [ADC]')
            plt.xlabel('Time [ns]')
            plt.show()                 #PLOT
                
        if self.sig_plot==True:
            charge = np.zeros(0)
            for i in range(0, n_wvf):
                c=np.zeros(len(np.where(carica[i]!=0)))
                c = carica[i, np.where(carica[i]!=0)]
                charge = np.append(charge, c)
            del carica
            carica = np.zeros(len(charge))
            carica = charge
            
            plt.figure(figsize=[10, 6])
            ax=plt.gca()
            ax.tick_params(bottom=True, top=True, left=True, right=True)
            ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
            plt.ylabel('Counts')
            plt.xlabel('Charge [ADC.ns]')
            plt.hist(carica, self.bin, range=[self.range_x1, self.range_x2], density=False)
            # plt.savefig('coldbox_aug_rand_trig_darkcounsearch.png', dpi=600)
            plt.show()

        return carica, e  #, peak, peakpos, sop, sot
    
#################################################################################


    def isto_carica(self, t, a, uno):
    
        if self.selection==True:
            n_wvf=int(sum(uno))
        else: n_wvf = len(a)
        
        charge = np.zeros(n_wvf)
        ind_min = int(self.t_min/self.hd.ns)
        ind_max = int(self.t_max/self.hd.ns)
        count=0
        
        for i in range(0, len(a)):
            
            if uno[i]==0: continue
            integral = np.trapz(a[count, ind_min:ind_max], t[count, ind_min:ind_max])
            charge[count] = integral
            count+=1
            if integral==0.:
                print(i)
        
        if self.range_x1==0 and self.range_x2==0: 
            if self.bin==0: n, b = np.histogram(charge)
            else: n, b = np.histogram(charge, self.bin, density=False)
            
        else: n, b = np.histogram(charge, self.bin, range=[self.range_x1, self.range_x2], density=False)
        
        if self.isto_plot==True:
            plt.figure(figsize=[8,6])
            ax=plt.gca()
            ax.tick_params(bottom=True, top=True, left=True, right=True)
            ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
            plt.ylabel('Counts')
            plt.xlabel('Charge [ADC.ns]')
            
            if self.range_x1==0 and self.range_x2==0: 
                if self.bin==0: plt.hist(charge)
                else: plt.hist(charge, self.bin, density=False) 
            else: plt.hist(charge, self.bin, range=[self.range_x1, self.range_x2], density=False)
            #plt.hist(charge, self.bin, density=False)
            
            #plt.savefig('ch_histo250', dpi=300)
            plt.show()
        
        return n, b, charge
    
###################################################################################



    def gauss(self, x, a, x0, sigma):
        return (abs(a)*np.exp(-0.5*((x-x0)/sigma)**2))
    
    
    def multi_gauss(self, x, *par):
    
        q = self.gauss(x, par[0], par[1], par[2]) + self.gauss(x, par[3], par[4],
                                                     par[5]) + self.gauss(x, par[6], par[7], (math.sqrt(2)*par[5]))
        dif = par[7] - par[4]
        media = par[7] + dif
        for i in range(0, int((len(par)-8))):
            q = q + self.gauss(x, par[i+8], media, math.sqrt(3+i)*par[5])
            media = media + dif
        return q
   
    def multi_gauss_n(self, x, *par):
       q = self.gauss(x, par[0], par[1], par[2]) + self.gauss(x, par[3], par[4],
                                                    par[5]) + self.gauss(x, par[6], par[7], par[8])
       dif = par[7] - par[4]
       media = par[7] + dif
       for i in range(0, int((len(par)-9)/2)):
           q = q + self.gauss(x, par[(2*i+9)], media, par[(2*i+10)])
           media = media + dif
       return q
  
###########################
    
    def multi_fit(self):#, bounds):
        # x of histogram as center of every bin for the fit
        b1 = center_bins(self.hd.b)
        p0 = self.p0
        for i in range(0, self.numg):
            p0 = np.append(p0, p0[6]/(3**(i+1)))
        
        popt, pcov = curve_fit(self.multi_gauss, b1, self.hd.n, p0,
                               maxfev=100000)#, bounds=(-10000, 20000))
        yc = self.multi_gauss(b1, *popt)
        x = np.linspace(self.range_x1, self.range_x2, 1000)
        ym = self.multi_gauss(x, *popt)
        if self.fitdue==True:
            p=[]
            j=1
            for i in range(0, len(popt)+self.numg+1):
                if i<8:
                    p.append(popt[i])
                elif i%2==0:
                    p.append(math.sqrt(3+(i-8))*popt[5])
                else:
                    p.append(popt[i-j])
                    j+=1
            len(p)
            popt, pcov = curve_fit(self.multi_gauss_n, b1, self.hd.n, p,
                                    maxfev=100000, bounds=(self.bound1, self.bound2))#, bounds=bounds)     # SECOND FIT
            yc = self.multi_gauss_n(b1, *popt)
            ym = self.multi_gauss_n(x, *popt)
    
        SNR = popt[4]/math.sqrt(popt[2]**2)# + popt[5]**2)
    
        #pop, pc = leastsq(tre_gauss, b1, n, p0, maxfev=10000)
        # popt=p0    #plot delle gaussiane singole e totale con i parametri iniziali che do io
    
        chi2 = 0
        for i in range(0, len(b1)):
            if self.hd.n[i]!=0: chi2 += (yc[i]-self.hd.n[i])**2/(self.hd.n[i]) # FOR CHI2 DENSITY=FALSE IN THIS HIST & IN ISTO_CHARGE HIST
            else: chi2 += (yc[i]-self.hd.n[i])**2
        chi2 = chi2/len(b1)
        
        if self.optimiz==False:
            plt.figure(figsize=[8.5, 6])                                            #PLOT
            plt.hist(self.hd.charge, self.hd.b, density=False)  
            plt.plot(x, ym, c='r', label='Best fit', linewidth=2, linestyle='-')   #PLOT
            
            dif = popt[7]-popt[4]
            media = popt[7]+dif
            c1 = popt[5]
            
            if (self.fitdue==True):
                
                for i in range(0, 3+self.numg):
                    
                    if i<3: 
                        fittxt = 'A%i= '%i+str(round(popt[3*(i)], 2))+'; mu%i= '%i + \
                        str(round(popt[3*i+1], 2))+'; s%i= '%i+str(round(popt[3*i+2], 2))
                        gi=self.gauss(x, popt[3*(i)], popt[3*i+1], popt[3*i+2])
                        
                    else: 
                        fittxt = 'A%i= '%i+str(round(popt[(2*i+3)], 2))
                        gi = self.gauss(x, popt[(2*i+3)], media, popt[(2*i+4)])
                        media = media + dif
    
                    plt.plot(x, gi, label=fittxt, linestyle='-') 
            else:    
                
                for i in range(0, 3+self.numg):
                    
                    if i<2: 
                        fittxt = 'A%i= '%i+str(round(popt[3*(i)], 2))+'; mu%i= '%i + \
                        str(round(popt[3*i+1], 2))+'; s%i= '%i+str(round(popt[3*i+2], 2))
                        gi=self.gauss(x, popt[3*(i)], popt[3*i+1], popt[3*i+2])
                        
                    elif i==2:
                        fittxt = 'A%i= '%i+str(round(popt[3*(i)], 2))+'; mu%i= '%i + \
                        str(round(popt[3*i+1], 2))
                        gi=self.gauss(x, popt[3*(i)], popt[3*i+1], (math.sqrt(i)*c1))

                    else: 
                        gi = self.gauss(x, popt[i+5], media, (math.sqrt(i)*c1))
                        fittxt = 'A%i= '%i+str(round(popt[i+5], 2))
                        media = media + dif
                        
                    plt.plot(x, gi, label=fittxt, linestyle='-') 
            
            ax=plt.gca()
            ax.tick_params(bottom=True, top=True, left=True, right=True)
            ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
            plt.ylabel('Counts')
            plt.xlabel('Charge [ADC.ns]')
            plt.legend(fontsize = 'xx-small', loc='upper right')                      #PLOT
            # plt.savefig('coldbox_aug_spe_spectrum-carica = 5177.131 ,  SNR =  4.57.png', dpi=600)
            plt.show()                                                             #PLOT 
            print('charge =', round(popt[7]-popt[4], 2), ',  SNR = ', round(SNR, 2), ",  chi2 =", round(chi2, 2))#, j)

    
        return popt, pcov, SNR, chi2
        
#####################################################################################


    def clean_events(self):
        
        n_wvf = self.n_wvf
        wvf_lenght = self.wvf_lenght
        puliti=np.ones(n_wvf)
        
        
        for i in range(0, n_wvf):
            if min(self.hd.a_filt[i,0:self.clean_x1])<self.lim_inf_clean or min(self.hd.a_filt[i,self.clean_x2:wvf_lenght])<self.lim_inf_clean:
                puliti[i]=0
                
        n_clean = int(sum(puliti))
        k = np.zeros((n_clean, wvf_lenght))
        kf = np.zeros((n_clean, wvf_lenght))

        j=0
        for i in range(0, n_wvf):
            if puliti[i]==1:
                k[j] = self.hd.a_filt[i]
                kf[j] = self.hd.a_dfilt[i]
                j+=1
        t0 = self.hd.time[0:n_clean]

        return k, kf, t0

#################################################################################


    def optimization(self, integration_time, binning):
    
    
        start1 = self.start1   #divisibili per 4???
        start2 = self.start2
        fin1 = self.fin1
        fin2 = self.fin2
        
        self.optimiz=True
        
        sn = np.zeros(0)
        mu1 = np.zeros(0)
        mu1dif = np.zeros(0)
        var = np.zeros(0)

        if (integration_time==True and binning==False):
            
            while True:
            
                t_min = int((start1)/self.hd.ns)   #4 ns per point
                t_max = int((start2)/self.hd.ns)
                    
                n, b, charge = self.isto_carica()
            
                pop, cov, SNR, chi2 = self.multi_fit()
                
                var = np.append(var, (t_max-t_min)*4)   #indici --> ns
                mu1 = np.append(mu1, pop[4])
                mu1dif = np.append(mu1dif, pop[7]-pop[4])
                sn = np.append(sn, SNR)
                # if SNR>3.026: print(start1, start2)
                
                if start1-self.hd.ns > fin1:
                    start1 = start1 - self.hd.ns
                if start2+self.hd.ns < fin2:
                    start2 = start2 + self.hd.ns
                else: break 
        
        if (integration_time==True and binning==False):

            for i in range(60, 300):
                
                n, b, charge = self.isto_carica()
            
                pop, cov, SNR, chi2 = self.multi_fit()
                
                var = np.append(var, i)   #indici --> ns
                mu1 = np.append(mu1, pop[4])
                mu1dif = np.append(mu1dif, pop[7]-pop[4])
                sn = np.append(sn, SNR)
                
        else: print('Error: put one variable TRUE and the other FALSE!')
        
        
        plt.figure()
        plt.plot(var, sn, '-')
        plt.title('SNR')
        plt.ylabel('SNR')
        plt.xlabel('Integration time [ns]')
        plt.show()
        
        plt.figure()
        plt.plot(var, mu1, '-')
        plt.title('mu1')
        plt.ylabel('mu1 [ADC.ns]')
        plt.xlabel('Integration time [ns]')
        plt.show()
        
        plt.figure()
        plt.plot(var, mu1dif, '-')
        plt.title('mu2-mu1')
        plt.ylabel('mu2-mu1 [ADC.ns]')
        plt.xlabel('Integration time [ns]')
        plt.show()

#################################################################


    # def __analyze__(self):
        
    #     if self.do=='signal':
    #         f = FILTER_DATA.signal(self)
    #     if self.do=='isto_carica':
    #         f = FILTER_DATA.isto_carica(self)
    #     if self.do=='mov_av':
    #         f = FILTER_DATA.mov_av(self)
    #     if self.do=='denoising':
    #         f = FILTER_DATA.denoising(self)
    #     return f
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    