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

#%% PLOT

def plot(t, a, qf, i):

    plt.figure(figsize=[8.5, 6])
    ax=plt.gca()
    ax.tick_params(bottom=True, top=True, left=True, right=True)
    ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
    plt.plot(t[i], a[i], label='raw data', )
    plt.plot(t[i], qf[i], label='filtered data')
    plt.legend()
    plt.xlabel('Time [ns]')
    plt.ylabel('Amplitude [ADC]')
    plt.show()
    
def zm():
    plt.close('all')
    
def histo_2d(t, a, x1, x2, step1, y1, y2, step2, draw, mean):
    
    xbi = np.arange(x1, x2, step1)
    ybi = np.arange(y1, y2, step2)
    
    plt.figure()
    h, yb, xb = np.histogram2d(a.flatten(), t.flatten(), bins=(ybi, xbi), range=((y1,y2),(x1,x2)), normed=None, weights=None, density=None)
    # h = np.log(h.T)
    #X, Y = np.meshgrid(xb, yb)
    plt.pcolormesh(xb, yb, np.log(h), cmap='jet')
    plt.xlabel('Time [ns]')
    plt.ylabel('Amplitude [ADC]')
    plt.colorbar()
    if draw==True:
        plt.plot(t[0], mean, c='black')
    plt.show()

    
# %% READ BINARY FILE


def read_binary(filez, header, adc_res, ns):

    file = np.fromfile(filez, dtype=np.int16)
    file32 = np.fromfile(filez, dtype=np.int32)  # header=6 invece di 12
    found = False
    pointsxwvf = int(file32[0]/2)  # 2 bin per dato
    for i in range(1, pointsxwvf):
        if file[i]==file[0]:
            pointsxwvf = int(i)
            # print('trovato', i)      #debug
            found = True
            break
    wvf = int(len(file)/pointsxwvf)
    points = int(pointsxwvf-header)
    # print(wvf, points)               #debug
    if found==False:
    # pointsxwvf32 = int(pointsxwvf/2) 
        wvf = int(len(file)/pointsxwvf)
        points = pointsxwvf-header   # tolgo l'header x ogni wvf

    a = np.empty((wvf, points))
    adc = np.empty((wvf, points))
    t = np.empty((wvf, points))
    time = np.empty(wvf)

    for i in range(0, wvf):

        adc[i] = file[(i*pointsxwvf)+header:pointsxwvf*(i+1)]
        t[i] = np.arange(0, ns*points, ns)  # 250 MHz = 4e-03 us
        # time[i] = file[pointsxwvf32*i+5]
    a = adc * adc_res  #valori in volt o mV
    # freq 62.5 Mhz
    out = int(len(a[0])/4)

    t62 = np.empty((len(a), out))
    a62 = np.empty((len(a), out))

    for i in range(0, out):

        t62[:, i] = t[:, 4*i]
        a62[:, i] = a[:, 4*i]

    return a, t, a62, t62, adc #, time

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


def read_long_txt(file, righe_iniz, puntixwvf, ns):
    

    a = np.genfromtxt(file, skip_header=righe_iniz)   #prende tanto tempo
    points = puntixwvf
    wvf = int(len(a)/points)  # =10000
    aa = np.empty((wvf, points))
    tt = np.empty((wvf, points))

    for i in range(0, wvf):

        aa[i] = a[points*i:points*(i+1)]
        tt[i] = np.arange(0, ns*points, ns)  # 250 MHz = 4e-03 us

    # freq 62.5 Mhz
    # out = int(len(aa[0])/4)

    # tt62 = np.empty((len(aa), out))
    # aa62 = np.empty((len(aa), out))

    # for i in range(0, out):

    #     tt62[:, i] = tt[:, 4*i]
    #     aa62[:, i] = aa[:, 4*i]

    return aa, tt#, aa62, tt62

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
    y = np.zeros((len(df), len(df[0])))
    a = a + 2  # shift

    for i in range(0, len(df)):

        z = convolve(a[i], b)  # filter output using convolution

        z = lfilter(b, c, a[i])  # filter output using lfilter function

        y[i] = z - 2
    
        y[i,0:l-1] = a[i,0:l-1] - 2   # back from shift

    return y  


# %% DENOISING

def denoising(df, a, f):

    signal_filtered = np.zeros((len(df), len(a[0])))

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

#%% BASELINE

def baseline(a, threshold, bins): #p0
    lenght_wvf = len(a[0])
    a0 = np.zeros((len(a), lenght_wvf))
    mu = np.zeros(len(a))
    # chi2t = np.zeros(len(a))
    
    # plt.figure()                           #PLOT
    for i in range (0, len(a)):
        j=0
        v = np.zeros(lenght_wvf)
        
        # plt.hist(a[i], 100)                 #PLOT
        n, b = np.histogram(a[i], bins)
        b1 = center_bins(b)
        # p, _ = curve_fit(gauss, b1, n)#p0
        # yc = gauss(b1, *p)
        # maxx = p[1]                  #media fit
        maxx = b1[np.argmax(n)]       # most frequent value
        # plt.plot(b1, yc, c='r', linewidth=2, linestyle='-')  #PLOT

        # chi2=0
        # for k in range(0, len(b1)):
        #     if n[k] != 0:
        #         # FOR CHI2 DENSITY=FALSE IN THIS HIST & IN ISTO_CHARGE HIST
        #         chi2 += (yc[k]-n[k])**2/(n[k])
        #     else:
        #         chi2 += (yc[k]-n[k])**2
        # chi2t[i] = chi2/len(b1)
        
        tre = maxx + threshold
        tre_ = maxx - threshold
        cont=0
        val = a[i]
        while j<lenght_wvf:
            if (val[j]<tre and val[j]>tre_):
                v[j] = val[j]
                cont += 1
                j = j + 1
            else: j=j+250                    #index at 250Mhz to jump 1us

        # x = np.linspace(3020, 3185, 84)              #PLOT
        # x1 = np.linspace(3020, 3230, 1000)           #PLOT
        mu[i] = sum(v)/cont                         #faccio la semplice media perchè è una gaussiana, non serve fare il fit
        # n, b = np.histogram(v, bins)
        # plt.hist(v, bins=x, alpha=0.5)                #PLOT
        # b1 = center_bins(b)
        # p, _ = curve_fit(gauss, b1, n, p0)
        # plt.plot(x1, gauss(x1, *p), c='r', label='Best fit', linewidth=2, linestyle='-')  #PLOT
        # mu[i] = p[1] 
        a0[i] = a[i] - mu[i]  #p[1]
    # plt.show()                             #PLOT
    return a0, mu#, chi2t

# %% MY SIGNAL FINDER (dev_std*sigma to discriminate + charge integral)

def signal(df, a, a_filt, t, threshold, prima, dopo, time_after_pulse, lim_amp, debug, w, n, t_min, t_max):


    carica = np.zeros((len(df), n))
    # sop = np.zeros((len(df), n))
    # sot = np.zeros((len(df), n))
    # peak = np.zeros((len(df), n))
    # peakpos = np.zeros((len(df), n))
    e = np.zeros(len(df))
    
    if debug==True: leng=w
    else: leng=len(df)
    
    for i in range(0, leng):


        tre = threshold  # mean[i] + sigma*abs(dev_std[i])   #THRESHOLD
        j = 0
        ev = 0
        trigger = False
        
        if debug==True:
            plt.figure()                        #PLOT
            ax=plt.gca()
            ax.tick_params(bottom=True, top=True, left=True, right=True)
            ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
            plt.plot(t[i], a[i], '-')           #PLOT
            plt.plot(t[i], a_filt[i], '-r')     #PLOT

        j=25     #SKIPPING THE FIRST 100 NS UNTIL I FIXED FILTER
        while j < len(a[i]):
            
            if j >= len(a[i])-2:
                break

            if (a_filt[i, j] > tre and trigger == False):

                jj = j
                trigger = True
                j = j + 15

            elif (a_filt[i, j] < tre and trigger == True):

                k = j
                trigger = False
                
                
                picco = a[i, jj:k]
                # peak[i, ev] = np.max(picco)
                if np.max(picco)>lim_amp: continue    #for data random wall leakage

                pk = np.argmax(picco)+jj
                if t_min!=0 and t_max!=0:
                    if pk>t_min and pk<t_max: continue

                # peakpos[i, ev] = pk
                # faccio l'integrale prendendo qualche pto prima e qualcuno dopo
                integral = np.trapz(a[i, pk-prima:pk+dopo],
                                    t[i, pk-prima:pk+dopo])

                carica[i, ev] = integral

                # indici sopra e sotto il threshold
                # sop[i, ev] = jj
                # sot[i, ev] = k
                
                if debug==True:
                    plt.plot(t[i, pk], np.max(picco), 'x', markersize=20)              #PLOT
                    plt.plot(t[i,pk-prima:pk+dopo], a_filt[i, pk-prima:pk+dopo], '-')      #PLOT
                # print(k)                                                             #debug
                
                # in modo che non rientri nel for per i punti appartenenti al picco + evita afterpulse
                j = k + time_after_pulse   #200 ns at 250Mhz
                ev += 1
                e[i] = int(ev)
            else: j += 1

    if debug==True:
        plt.ylabel('Amplitude [ADC]')
        plt.xlabel('Time [ns]')
        plt.show()                 #PLOT

    return carica, e#, peak, peakpos, sop, sot

# %% istogramma carica fotopicchi


def isto_carica(df, aa, tt, t_min, t_max, _bin, x, y):

    charge = np.empty(len(df))

    for i in range(0, len(df)):

        integral = np.trapz(aa[i, t_min:t_max], tt[i, t_min:t_max])
        charge[i] = integral
        if integral==0.:
            print(i)

    # print(t_min,t_max)
    if x==0 and y==0: 
        if _bin==0: n, b = np.histogram(charge)
        else: n, b = np.histogram(charge, _bin, density=False)
    else: n, b = np.histogram(charge, _bin, range=[x, y], density=False)

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

#%% CENTER BIN FOR FIT

def center_bins(b):
    
    bbin = np.zeros(len(b))
    l = (b[2]-b[1])/2
    bbin = l + b  # np.resize(b, len(b))
    b1 = np.delete(bbin, len(b)-1)
    
    return b1


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

def multi_gauss_n(x, *par):
   q = gauss(x, par[0], par[1], par[2]) + gauss(x, par[3], par[4],
                                                par[5]) + gauss(x, par[6], par[7], par[8])
   dif = par[7] - par[4]
   media = par[7] + dif
   for i in range(0, int((len(par)-9)/2)):
       q = q + gauss(x, par[(2*i+9)], media, par[(2*i+10)])
       media = media + dif
   return q


def multi_fit(numg, p0, n, b, charge, x, fitdue, optimiz):

    # x of histogram as center of every bin for the fit
    b1 = center_bins(b)

    popt, pcov = curve_fit(multi_gauss, b1, n, p0,
                           maxfev=100000, bounds=(-10000, 20000))
    yc = multi_gauss(b1, *popt)
    ym = multi_gauss(x, *popt)
    if fitdue==True:
        p=[]
        j=1
        for i in range(0, len(popt)+numg+1):
            if i<8:
                p.append(popt[i])
            elif i%2==0:
                p.append(math.sqrt(3+(i-8))*popt[5])
            else:
                p.append(popt[i-j])
                j+=1
            
        popt, pcov = curve_fit(multi_gauss_n, b1, n, p,
                                maxfev=1000000, bounds=(-100, 50000))     # SECOND FIT
        yc = multi_gauss_n(b1, *popt)
        ym = multi_gauss_n(x, *popt)

    SNR = popt[4]/math.sqrt(popt[2]**2)# + popt[5]**2)

    #pop, pc = leastsq(tre_gauss, b1, n, p0, maxfev=10000)
    # popt=p0    #plot delle gaussiane singole e totale con i parametri iniziali che do io

    chi2 = 0
    for i in range(0, len(b1)):
        if n[i]!=0: chi2 += (yc[i]-n[i])**2/(n[i]) # FOR CHI2 DENSITY=FALSE IN THIS HIST & IN ISTO_CHARGE HIST
        else: chi2 += (yc[i]-n[i])**2
    chi2 = chi2/len(b1)
    
    if optimiz==False:
        plt.figure(figsize=[8.5, 6])                                                         #PLOT
        plt.hist(charge, b, density=False)                                     #PLOT
        plt.plot(x, ym, c='r', label='Best fit', linewidth=2, linestyle='-')   #PLOT
        
        dif = popt[7]-popt[4]
        media = popt[7]+dif
        c1 = popt[5]
        
        if (fitdue==True):
            
            for i in range(0, 3+numg):
                
                if i<3: 
                    fittxt = 'A%i= '%i+str(round(popt[3*(i)], 2))+'; mu%i= '%i + \
                    str(round(popt[3*i+1], 2))+'; s%i= '%i+str(round(popt[3*i+2], 2))
                    gi=gauss(x, popt[3*(i)], popt[3*i+1], popt[3*i+2])
                    
                else: 
                    fittxt = 'A%i= '%i+str(round(popt[(2*i+3)], 2))
                    gi = gauss(x, popt[(2*i+3)], media, popt[(2*i+4)])
                    media = media + dif

                plt.plot(x, gi, label=fittxt, linestyle='-') 
        else:       
            g0 = gauss(x, popt[0], popt[1], popt[2])
            g1 = gauss(x, popt[3], popt[4], popt[5])
            g2 = gauss(x, popt[6], popt[7], (math.sqrt(2)*popt[5]))
            
        
            fittxt0 = 'A0= '+str(round(popt[0], 5))+'; mu0= ' + \
                str(round(popt[1], 5))+'; s0= '+str(round(popt[2], 5))
            fittxt1 = 'A1= '+str(round(popt[3], 5))+'; mu1= ' + \
                str(round(popt[4], 5))+'; s1= '+str(round(popt[5], 5))
            fittxt2 = 'A2= '+str(round(popt[6], 5))+'; mu2= '+str(round(popt[7], 5))
    
            plt.plot(x, g0, label=fittxt0, linestyle='-')                          #PLOT
            plt.plot(x, g1,  label=fittxt1, linestyle='-')                         #PLOT
            plt.plot(x, g2,  label=fittxt2, linestyle='-')                         #PLOT
        
            for i in range(0, numg):
                gi = gauss(x, popt[i+8], media, (math.sqrt(3+i)*c1))
                plt.plot(x, gi, label='A%i= ' %                                    #PLOT
                          (i+3)+str(round(popt[i+8], 5)), linestyle='-')           #PLOT
                media = media + dif
        
        ax=plt.gca()
        ax.tick_params(bottom=True, top=True, left=True, right=True)
        ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
        plt.ylabel('Counts')
        plt.xlabel('Charge [mV.ns]')
        plt.legend(fontsize = 'xx-small', loc='upper right')                      #PLOT
        # plt.savefig('coldbox_aug_spe_spectrum-carica = 5177.131 ,  SNR =  4.57.png', dpi=600)

        plt.show()                                                             #PLOT 

    return popt, pcov, SNR, chi2
