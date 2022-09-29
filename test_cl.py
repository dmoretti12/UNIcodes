# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 03:32:54 2022

@author: Davide
"""

from tutte_le_mie_cazzate import*

hd = [DATI() , DATI()]  #id() per indirizzo memoria, controllare che sia lo stesso altrimenti copia
hd[0].ns = 4
h1 = READ_DATA(hd[0])
h1.type='binary'
#%% AUG LONF WVF
os.chdir("H:\Il mio Drive\WORK\AUG22COLDBOX")
h1.file="0_wave1_light_leakage_test_led_flange_covered.dat"   # 1139.5 Hz/mm^2 ev 410.2 k  - tre=37 15,100, filt = 7 17)
# file="0_wave1_light_leakage_test_led_fibers_connected.dat"   #  1099.0 Hz/mm^2 ev 395.6 k - tre=37 15,100, filt = 7 17)
# file="0_wave1_light_leakage_test_led_fibers_connected_dcem_digital_off.dat"  #  910.6 Hz/mm^2 ev 327.8 k - tre 37 15,100, filt= 7 17)

hd[0].a,  hd[0].t, _, _ = h1.__read__()
os.chdir("H:\Il mio Drive\WORK")
#%% AUG SHOT WVF
hd[1].ns = 4
g1 = READ_DATA(hd[1])
os.chdir("H:\Il mio Drive\WORK\AUG22COLDBOX")

file1="0_wave1_led_7V_50ns_digital_on.dat"  # 633.2 Hz/mm^2  ev 45.6 k     carica = 5101.268 ,  SNR =  4.17 ,  chi2 = 1.58 
file2="1_wave1_led_7V_50ns_digital_on.dat" #  647.7 Hz/mm^2  ev 46.6 k    carica = 5184.381 ,  SNR =  4.13 ,  chi2 = 1.3
file3='2_wave1_led_7V_50ns_digital_on.dat' #  647.0 Hz/mm^2   ev 46.6 k     carica = 5042.97 ,  SNR =  4.13 ,  chi2 = 1.64
files = [file1, file2, file3]
n_file = len(files)

hd[1].a_wout_baseline,  hd[1].t = h1.concatenate(n_file, files)
#%%
hd[1].ns = 4
g1 = READ_DATA(hd[1])

file1="0_wave0_35V20_20ADC_2860mV_50ns.dat"   #14 bit digitizer
file2="1_wave0_35V20_20ADC_2860mV_50ns.dat" 
# file3='2_wave0_34V40_20ADC_2860_50ns.dat' 
# file4='3_wave0_34V40_20ADC_2860_50ns.dat' 
files = [file1, file2]#, file3, file4]
n_file = len(files)

hd[1].a_wout_baseline,  hd[1].t = h1.concatenate(n_file, files)
#%%
hd[1].a = hd[1].a_wout_baseline
#TO DO: save processed data in a dataframe

h2 = FILTER_DATA(hd[1])

hd[1].a_filt_wout_baseline = h2.mov_av(hd[1].a_wout_baseline)
h2.threshold_base = 30
h2.discard_tre = 0
h2.debug_base = False
h2.exclusion_window = 1000 #ns

hd[1].a_filt, mu, cont, uno = h2.baseline(hd[1].a_filt_wout_baseline)   
hd[1].a = h2.sottraggo_baseline(hd[1].a, hd[1].a_wout_baseline, mu)    #sottraggo baseline a dati non filtrati (grezzi)

h2.l_mov_av = 17
hd[1].a_dfilt = h2.mov_av(hd[1].a)     # applico filtro +forte a dati grezzi con baseline
#%% SAVE IN TXT

with open('baseline_LED_22aug.txt', 'w') as f:
    for i in range(len(mu)):
        f.write(str(mu[i]))
        f.write('\n')
f.close()
asd=np.genfromtxt('baseline_LED_22aug.txt')

#%%SALVA

np.save('a_filt', hd[1].a_filt, allow_pickle=True, fix_imports=True)
afilt=np.load('a_filt.npy', mmap_mode=None, allow_pickle=True, fix_imports=True, encoding='bytes')
#%%
histo_2d(hd[1].t, hd[1].a_filt, x1=0, x2=hd[1].t[0,-1], step1=hd[1].ns, y1=-100, y2=200, step2=6)
#%%

h3 = ANALYZE_DATA(hd[1])
h3.selection=False
h3.isto_plot = True
h3.t_min = 5050
h3.t_max = 5450
h3.range_x1 = -3500
h3.range_x2 = 30000
h3.bin = 200
hd[1].n, hd[1].b, hd[1].charge = h3.isto_carica(hd[1].t, hd[1].a_filt, uno)
#%%
h3.p0 = [470, -2000, 800, 240, 2900, 1000, 600, 7000]
h3.numg = 5
h3.bound1 = -500
h3.bound2 = 25000
h3.fitdue = True
# bounds=[(-np.inf, np.inf) for i in range((h3.numg)*2+8)] #costrains on last 2 gaussians
# z = (0, popt[-1]*0.2)
# bounds.append(z)
# del z
popt, pcov, SNR, chi2 = h3.multi_fit(bounds)



#%%
h3 = ANALYZE_DATA(hd[1])
h3.signal_threshold = 18.5 #15
h3.prima_t = 30
h3.dopo_t = 320
h3.debug_signal = False
h3.time_after_pulse = 1000 
h3.sig_plot=True
h3.t_min_sig = 5030
h3.t_max_sig = 5500
charge, events = h3.signal(hd[1].t, hd[1].a_filt, hd[1].a_dfilt)


# TRIGGER RATE
rate = sum(events)/(len(hd[1].a)*len(hd[1].a[0])*4e-9*(20*36))  #Hz
print(round(rate,3), 'Hz/mm^2')
print('ev', round(sum(events)*1e-3,1),'k')


#%%%%
# fai in modo di leggere automaticamente file da diverse cartele di voltage
# e di salvarle in un array, per poi fare un analisi di linearità charge(y)-voltage(x), un voltaggio per volta altrimenti troppi dati
import os, glob
import pandas as pd
n=12
for i in range (n):
    
    os.chdir('%d', d)
    a[i]=0

#filenames = [i for i in glob.glob(self.iniz)]

#WAVEi CORRISPONDE AL CANALE (CH) i !!!!!!!! per cui analizza tutte le wave0 di tutti i voltaggi, fai 1) e 2) poi passa a wave1... distingua sempre per canale!!

# 1) for -> apri file -> leggi dati e fai isto_carica -> prendi media(semplice xkè gaussiana) -> prossimo voltaggio

# 2) fai retta ampl_charge con tutti i voltaggi, sul medesimo plot, differenziando voltaggi per colore









