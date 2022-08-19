# -*- coding: utf-8 -*-
"""
Created on Tue May 10 17:47:11 2022

@author: Davide
"""
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

import os
import glob
import pandas as pd
from scipy.signal import lfilter, convolve
from scipy.stats import norm

# %% READ BINARY FILE


def read_binary(filez, header, adc_res):

    file = np.fromfile(filez, dtype=np.int16)
    # file32 = np.fromfile(filez, dtype=np.int32)

# 2500 points x wvf + 12 header = 2512
# len(file)/2512=61674624/2512= 24552 wvf

    pointsxwvf = int(file[0]/2)  # 2 bin per dato
    # pointsxwvf32 = int(pointsxwvf/2)
    wvf = int(len(file)/pointsxwvf)
    points = pointsxwvf-header   # tolgo l'header x ogni wvf

    a = np.empty((wvf, points))
    t = np.empty((wvf, points))
    time = np.empty(wvf)

    for i in range(0, wvf):

        a[i] = file[(i*pointsxwvf)+header:pointsxwvf*(i+1)] * adc_res
        t[i] = np.arange(0, 4e-03*points, 4e-03)  # 250 MHz = 4e-03 us
        # time[i] = file32[pointsxwvf32*i+5]

    # freq 62.5 Mhz
    out = int(len(a[0])/4)

    t62 = np.empty((len(a), out))
    a62 = np.empty((len(a), out))

    for i in range(0, out):

        t62[:, i] = t[:, 4*i]
        a62[:, i] = a[:, 4*i]

    return a, t, a62, t62#, time

# numpy.fromfile(file #string, dtype=float, count=- 1 #quali dati prendere, sep='', offset=0 #skippa primi dati, *, like=None)

# %% READ MANY TXT


def read_many_txt(string, iniz):

    os.chdir(string)

    filenames = [i for i in glob.glob(iniz)]

    #df = pd.read_csv('C1XArapuca_CRT_250MHz_00000.txt', sep=" ", header=5, index_col = None, encoding = 'unicode_escape')
    df = [pd.read_csv(file, sep=",", header=4, engine='python', index_col=None, encoding='unicode_escape')
          for file in filenames]

    t = np.empty((len(df), len(df[0])))  # freq 250 Mhz
    a = np.empty((len(df), len(df[0])))

    for i in range(0, len(df)):

        q = df[i].iloc[:, 0].values
        w = df[i].iloc[:, 1].values

        if len(df[i]) < len(df[0]):

            q = np.append(q, 0)
            w = np.append(w, 0)

       # T = np.max(t) + df[i].iloc[:, 0].values #this way (iloc) to access to the values of df visualizing le waveform tutte di seguito
        t[i, :] = q
        a[i, :] = w

    # 250/62.5
    #Out[4]: 4.0

    out = int(len(df[0])/4)

    tnn = np.empty((len(df), out))  # freq 62.5 Mhz
    ann = np.empty((len(df), out))

    for i in range(0, out):

        tnn[:, i] = t[:, 4*i]
        ann[:, i] = a[:, 4*i]

    return df, a, t, ann, tnn

# %% READ ONE LONG TXT


def read_long_txt(file, righe_iniz, puntixwvf):
    

    a = np.genfromtxt(file, skip_header=righe_iniz)
    points = puntixwvf
    wvf = int(len(a)/points)  # =10000
    aa = np.empty((wvf, points))
    tt = np.empty((wvf, points))

    for i in range(0, wvf):

        aa[i] = a[points*i:points*(i+1)]
        tt[i] = np.arange(0, 4e-03*points, 4e-03)  # 250 MHz = 4e-03 us

    # freq 62.5 Mhz
    out = int(len(aa[0])/4)

    tt62 = np.empty((len(aa), out))
    aa62 = np.empty((len(aa), out))

    for i in range(0, out):

        tt62[:, i] = tt[:, 4*i]
        aa62[:, i] = aa[:, 4*i]

    return aa, tt, aa62, tt62

#%% long txt w/ time

def read_long_txt_time(file, righe_iniz, puntixwvf):

    f = np.genfromtxt(file, skip_header=righe_iniz)
    t = f[:,0]
    a = f[:,1]
    points = puntixwvf
    wvf = int(len(a)/points)  # =10000
    aa = np.empty((wvf, points))
    tt = np.empty((wvf, points))

    for i in range(0, wvf):

        aa[i] = a[points*i:points*(i+1)]
        tt[i] = t[points*i:points*(i+1)]  # 250 MHz = 4e-03 us

    # freq 62.5 Mhz
    out = int(len(aa[0])/4)

    tt62 = np.empty((len(aa), out))
    aa62 = np.empty((len(aa), out))

    for i in range(0, out):

        tt62[:, i] = tt[:, 4*i]
        aa62[:, i] = aa[:, 4*i]

    return aa, tt, aa62, tt62

# %% MOVING AVERAGE + SHIFT FILTER


def mov_av(df, a, l):

    L = l  # L-point filter

    b = (np.ones(L))/L  # numerator co-effs of filter transfer function
    c = np.ones(1)  # denominator co-effs of filter transfer function
    y = np.empty((len(df), len(df[0])))
    a = a + 2  # shift

    for i in range(0, len(df)):

        z = convolve(a[i], b)  # filter output using convolution

        z = lfilter(b, c, a[i])  # filter output using lfilter function

        y[i] = z

    return y - 2  # back from shift


# %% DENOISING

def denoising(df, a, f):

    signal_filtered = np.empty((len(df), len(a[0])))

    for i in range(0, len(df)):

        n = len(a[i])
        fhat = np.fft.fft(a[i], n)  # computes the fft
        psd = fhat * np.conj(fhat)/n
        # freq = (1/(delta*n)) * np.arange(n) #frequency array
        idxs_half = np.arange(1, np.floor(
            n/2), dtype=np.int32)  # first half index
        psd_real = np.abs(psd[idxs_half])  # amplitude for first half

        # Filter out noise
        sort_psd = np.sort(psd_real)[::-1]
        # print(len(sort_psd))
        threshold = sort_psd[f]
        psd_idxs = psd > threshold  # array of 0 and 1
        # psd_clean = psd * psd_idxs #zero out all the unnecessary powers
        fhat_clean = psd_idxs * fhat  # used to retrieve the signal

        signal_filtered[i] = np.fft.ifft(
            fhat_clean)  # inverse fourier transform

    return signal_filtered

# %% FIT BASELINE con norm.pdf (mean and std_dev classiche)


def norm_fit(df, a, t, inizio_sig):

    mu = np.empty(len(df))
    std = np.empty(len(df))
    # pp = []

    for i in range(0, len(df)):

        # prendo come baseline solo dati prima del segnale
        yy = a[i, t[i] < inizio_sig]
        # zz = np.where(np.max(yy) < 0.215 and np.min(yy) > 0.206, yy)
        # print(len(yy), len(zz))

        # mu e std sono esattamente np.mean() e suo errore!!!
        mu[i], std[i] = norm.fit(yy)

        #xmin, xmax = plt.xlim()

        # x = np.linspace(-0.03, 0.03, 100)

        # p = norm.pdf(x, mu[i], std[i])  # best_fit_line
        # pp.append(p)

        # plt.figure()
        # plt.plot(x, p, 'k', linewidth=2)   #per vedere il fit p=norm.pdf(np.linspace(-0.03, 0.03), mu[i], std[i])
        # plt.show()

    return mu, std#, p


# %% MY SIGNAL FINDER (dev_std*sigma to discriminate + charge integral)

def signal(df, a, a_filt, t, threshold, prima, dopo):


    # carica = []
    ss = []
    # peak = []
    # peakpos = []
    carica = np.zeros((len(df), 50))
    sop = np.zeros((len(df), 50))
    sot = np.zeros((len(df), 50))
    peak = np.zeros((len(df), 50))
    peakpos = np.zeros((len(df), 50))
    e = np.zeros(len(df))

    for i in range(0, len(df)):#10, 24):#len(df)):


        tre = threshold[i]  # mean[i] + sigma*abs(dev_std[i])   #THRESHOLD
        j = 0
        # c = []   # array integrale segnale (charge) della i-esima waveform
        s = []   # indici segnali
        # p = []   # peaks
        # pt = []  # peaks position
        
        ev = 0
        trigger = False
        # plt.figure()                        #PLOT
        # plt.plot(t[i], a[i], '-')           #PLOT
        # plt.plot(t[i], a_filt[i], '-r')     #PLOT
        # re=False                            #debug
        j=0
        while j < len(a[i]):
            # if (re==True):   #debug
            #     print(j)
            # re=False
            if j >= len(a[i])-2:
                break

            if (a_filt[i, j] > tre and trigger == False):

                jj = j
                trigger = True
                j = j + 25

            elif (a_filt[i, j] < tre and trigger == True):

                k = j
                trigger = False
                # re=True     #debug
                
                # faccio l'integrale prendendo qualche pto prima e qualcuno dopo
                integral = np.trapz(a_filt[i, jj-prima:k+dopo],
                                    t[i, jj-prima:k+dopo])
                # c.append(integral)
                carica[i, ev] = integral
                s.append([jj, k])  # indici sopra e sotto il threshold
                sop[i, ev] = jj
                sot[i, ev] = k
                
                picco = a_filt[i, jj:k]
                # p.append(np.max(picco))
                peak[i, ev] = np.max(picco)
                # pt.append(np.argmax(picco))
                peakpos[i, ev] = np.argmax(picco)
                
                # plt.plot(t[i, np.argmax(picco)+jj], np.max(picco), 'x')              #PLOT
                # plt.plot(t[i,jj-prima:k+dopo], a_filt[i, jj-prima:k+dopo], '-')      #PLOT
                # print(k)                                                             #debug
                
                # in modo che non rientri nel for per i punti appartenenti al picco + evita afterpulse
                j = k + 90 #int(0.05*len(t[i]))
                ev = ev + 1
                e[i] = ev
            else: j = j + 1


        # plt.show()                 #PLOT
        
        # carica.append(c)
        # ss.append(s)
        # peak.append(p)
        # peakpos.append(pt)

    return carica, peak, peakpos, sop, sot, e

# %% istogramma carica fotopicchi


def isto_carica(df, aa, tt, t_min, t_max, _bin):

    charge = np.empty(len(df))

    for i in range(0, len(df)):

        # for j in range(t_min, t_max):

        # t_int = tt[i, tt[i] > tmin]
        # a_int = aa[i, tt[i] > tmin]
        # yy = a_int[t_int < tmax]
        # zz = t_int[t_int < tmax]

        integral = np.trapz(aa[i, t_min:t_max], tt[i, t_min:t_max])
        charge[i] = integral

    # print(t_min,t_max)

    n, b = np.histogram(charge, _bin, density=False)

    return n, b, charge

# %% TOT
# from pynverse import inversefunc

from scipy.interpolate import Rbf, InterpolatedUnivariateSpline

def tot(df, a, t, t_min, t_max, threshold, satur):

    carica = np.zeros(len(df))
    ss = np.zeros((len(df), 2))
    peak = np.zeros(len(df))
    peakpos = np.zeros(len(df))
    wvf = np.zeros(0)
    
    dig=False
    if (t[0,1]-t[0,0])>5: dig=True

    for i in range(0, len(df)):
        
        # if lowf==True:
            # model = get_natural_cubic_spline_model(t[i], a[i], minval=min(t[i]), maxval=max(t[i]), n_knots=len(a[i])/2)
            # a[i] = model.predict(np.linspace(min(t[i]), max(t[i]), len(df[0])))
            
        
        tre = threshold  # mean[i] + sigma*abs(dev_std[i])   #THRESHOLD

        trigger = False

        for j in range(0, len(a[i])):

            if j >= len(a[i])-2:
                break

            # and a[i, j+1] > tre:  #and a[i, j+2] > tre:   #if at least 3 points are above the 3-4*dev_std i consider it as a peak
            if (a[i, j] > tre and trigger == False):

                jj = j
                trigger = True

            elif (a[i, j] < tre and trigger == True):

                k = j
                # trigger = False

                # faccio l'integrale prendendo qualche pto prima e qualcuno dopo
                integral = np.trapz(a[i, t_min:t_max], t[i, t_min:t_max])

                carica[i] = integral
                # ss[i] = [jj,k]

                picco = a[i, jj:k]
 #               tpicco = t[i, jj:k]

                # m1 = (a[i, jj]-a[i, jj-1])/(t[i, jj]- t[i, jj-1])
                # q1 = (t[i, jj-1]*a[i, jj] - t[i, jj]*a[i, jj-1])/(t[i, jj-1]-t[i, jj])
                # x1 = (tre - q1)/m1
                # m2 = (a[i, k]-a[i, k-1])/(t[i, k]- t[i, k-1])
                # q2 = (t[i, k-1]*a[i, k] - t[i, k]*a[i, k-1])/(t[i, k-1]-t[i, k])
                # x2 = (tre - q2)/m2
                # ss[i] = [x1, x2]
                

                if dig==True and np.max(picco)<satur:
                    ius = Rbf(t[i, t_min:t_max], a[i, t_min:t_max])#InterpolatedUnivariateSpline(t[i], a[i])
                    xi = np.linspace(t[i, t_min], t[i, t_max], len(a[i, t_min:t_max])*4)
                    yi = ius(xi)
                    # if i == 2:
                    #     plt.figure()
                    #     plt.plot(t[i], a[i], '.')
                    #     plt.plot(xi, yi, '.')
                    #     plt.savefig('te', dpi=300)
                    tr=False
                    for k in range(0, len(yi)):
                        if (yi[k]>tre and tr==False):
                            tr = True
                            m1 = (yi[k]-yi[k-1])/(xi[k]- xi[k-1])
                            q1 = (xi[k-1]*yi[k] - xi[k]*yi[k-1])/(xi[k-1]-xi[k])
                            x1 = (tre - q1)/m1
                        elif (yi[k]<tre and tr==True):
                            m2 = (yi[k]-yi[k-1])/(xi[k]- xi[k-1])
                            q2 = (xi[k-1]*yi[k] - xi[k]*yi[k-1])/(xi[k-1]-xi[k])
                            x2 = (tre - q2)/m2
                            break
                    ss[i] = [x1,x2]

                    peak[i] = np.max(yi)
                    peakpos[i] = xi[np.argmax(yi)]
                    # wvf = np.append(wvf, i)
                else:
                # if np.max(picco) > 1.3:
                #     print(i, np.max(picco))
                    m1 = (a[i, jj]-a[i, jj-1])/(t[i, jj]- t[i, jj-1])
                    q1 = (t[i, jj-1]*a[i, jj] - t[i, jj]*a[i, jj-1])/(t[i, jj-1]-t[i, jj])
                    x1 = (tre - q1)/m1
                    m2 = (a[i, k]-a[i, k-1])/(t[i, k]- t[i, k-1])
                    q2 = (t[i, k-1]*a[i, k] - t[i, k]*a[i, k-1])/(t[i, k-1]-t[i, k])
                    x2 = (tre - q2)/m2
                    ss[i] = [x1, x2]
                    peak[i] = np.max(picco)
                    peakpos[i] = t[i, np.argmax(a[i])]

                
                break
                
                # plt.figure()
                # plt.plot(t[i,j-2:k+2], a[i,j-2:k+2], '.')
                # plt.show()

                # j = k + int(0.1*len(t)) #in modo che non rientri nel for per i punti appartenenti al picco + evita afterpulse

    return carica, ss, peak, peakpos#, wvf


#%% time resolution

def time_res(df, a, t, t_min, t_max, threshold):
    
    peak = np.empty(len(df))
    peakpos = np.empty(len(df))
    sigpos = np.empty(len(df))
    
    dig=False
    if (t[0,1]-t[0,0])>5: dig=True

    for i in range(0, len(df)):
        
        tre = threshold

        # trigger = False

        for j in range(t_min, t_max):

            if (a[i, j] > tre):

                jj = j

            # elif (a[i, j] < tre):

                # m1 = (a[i, jj]-a[i, jj-1])/(t[i, jj]- t[i, jj-1])
                # q1 = (t[i, jj-1]*a[i, jj] - t[i, jj]*a[i, jj-1])/(t[i, jj-1]-t[i, jj])
                # x1 = (tre - q1)/m1
                # sigpos[i] = x1
        
                if dig==True:
                    ius = Rbf(t[i, t_min:t_max], a[i, t_min:t_max])#InterpolatedUnivariateSpline(t[i], a[i])
                    xi = np.linspace(t[i, t_min], t[i, t_max], len(a[i, t_min:t_max])*4)
                    yi = ius(xi)
                    for k in range(0, len(yi)):
                        if (yi[k]>tre):
                            m1 = (yi[k]-yi[k-1])/(xi[k]- xi[k-1])
                            q1 = (xi[k-1]*yi[k] - xi[k]*yi[k-1])/(xi[k-1]-xi[k])
                            x1 = (tre - q1)/m1
                            sigpos[i] = x1
                            break
                    # if i == 0:
                    #     plt.figure()
                    #     plt.plot(t[i], a[i], '.')
                    #     plt.plot(xi, yi, '.')
                    peak[i] = np.max(yi)
                    peakpos[i] = xi[np.argmax(yi)]
                else: 
                    m1 = (a[i, jj]-a[i, jj-1])/(t[i, jj]- t[i, jj-1])
                    q1 = (t[i, jj-1]*a[i, jj] - t[i, jj]*a[i, jj-1])/(t[i, jj-1]-t[i, jj])
                    x1 = (tre - q1)/m1
                    sigpos[i] = x1
                    picco = a[i, t_min:t_max]
                    peak[i] = np.max(picco)
                    peakpos[i] = t[i, np.argmax(a[i])]
                    
                break

    return peak, peakpos, sigpos
    
    
    
# %% ampl vs p.e.


def ampl_pe(df, a, t, t_min, t_max):

    carica = np.empty(len(df))
    peak = np.empty(len(df))
    
    dig=False
    if (t[0,1]-t[0,0])>5: dig=True

    for i in range(0, len(df)):

        integral = np.trapz(a[i, t_min:t_max], t[i, t_min:t_max])
        carica[i] = integral
        
        if dig==True:
            ius = Rbf(t[i, t_min:t_max], a[i, t_min:t_max])#InterpolatedUnivariateSpline(t[i], a[i])
            xi = np.linspace(t[i, t_min], t[i, t_max], len(a[i, t_min:t_max])*4)
            yi = ius(xi)
            # plt.figure()
            # plt.plot(t[i], a[i], '.')
            # plt.plot(xi, yi, '.')
            peak[i] = np.max(yi)
        else: 
            picco = a[i, t_min:t_max]
            peak[i] = np.max(picco)

    return carica, peak

# %% FIT PHOTOPEAKS SPECTRUM


def gauss(x, a, x0, sigma):
    return (abs(a)*np.exp(-0.5*((x-x0)/sigma)**2))


def multi_gauss(x, *par):

    q = gauss(x, par[0], par[1], par[2]) + gauss(x, par[3], par[4],
                                                 par[5]) + gauss(x, par[6], par[7], (math.sqrt(2)*par[5]))

    dif = par[7] - par[4]
    media = par[7] + dif

    for i in range(0, int((len(par)-8))):
        q = q + gauss(x, par[i+8], media, math.sqrt(3+i)*par[5])
        media = media + dif

    return q


def multi_fit(numg, p0, n, b, charge, x):

    # x of histogram as center of every bin for the fit
    bbin = np.zeros(len(b))
    l = (b[2]-b[1])/2
    bbin = l + b  # np.resize(b, len(b))
    b1 = np.delete(bbin, len(b)-1)

    popt, pcov = curve_fit(multi_gauss, b1, n, p0,
                           maxfev=10000, bounds=(-0.055, 600))

    SNR = popt[4]/math.sqrt(popt[2]**2 + popt[5]**2)

    #pop, pc = leastsq(tre_gauss, b1, n, p0, maxfev=10000)
    # popt=p0    #plot delle gaussiane singole e totale con i parametri iniziali che do io
    yc = multi_gauss(b1, *popt)
    ym = multi_gauss(x, *popt)
    chi2 = 0
    for i in range(0, len(b1)):

        if n[i] != 0:
            # FOR CHI2 DENSITY=FALSE IN THIS HIST & IN ISTO_CHARGE HIST
            chi2 += (yc[i]-n[i])**2/(n[i])
        else:
            chi2 += (yc[i]-n[i])**2
    chi2 = chi2/len(b1)

    # plt.figure()                                                           #PLOT
    # plt.hist(charge, b, density=False)                                     #PLOT
    # plt.plot(x, ym, c='r', label='Best fit', linewidth=2, linestyle='-')   #PLOT

    g0 = gauss(x, popt[0], popt[1], popt[2])
    g1 = gauss(x, popt[3], popt[4], popt[5])
    g2 = gauss(x, popt[6], popt[7], (math.sqrt(2)*popt[5]))

    dif = popt[7]-popt[4]
    media = popt[7]+dif
    c1 = popt[5]
    fittxt0 = 'A0= '+str(round(popt[0], 5))+'; mu0= ' + \
        str(round(popt[1], 5))+'; s0= '+str(round(popt[2], 5))
    fittxt1 = 'A1= '+str(round(popt[3], 5))+'; mu1= ' + \
        str(round(popt[4], 5))+'; s1= '+str(round(popt[5], 5))
    fittxt2 = 'A2= '+str(round(popt[6], 5))+'; mu2= '+str(round(popt[7], 5))
    # plt.plot(x, g0, label=fittxt0, linestyle='-')                          #PLOT
    # plt.plot(x, g1,  label=fittxt1, linestyle='-')                         #PLOT
    # plt.plot(x, g2,  label=fittxt2, linestyle='-')                         #PLOT

    for i in range(0, numg):
        gi = gauss(x, popt[i+8], media, (math.sqrt(3+i)*c1))
        # plt.plot(x, gi, label='A%i= ' %                                    #PLOT
        #           (i+3)+str(round(popt[i+8], 5)), linestyle='-')            #PLOT
        media = media + dif

    # plt.legend()                                                           #PLOT
    # plt.show()                                                             #PLOT 

    return popt, pcov, SNR, chi2
