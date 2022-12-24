# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 18:53:16 2022

@author: davide
"""

from tutte_le_mie_cazzate import*

hd = [DATI()]
#%% DATI LED  argon2x2
j=0
hd[j].ns = 4
h1 = READ_DATA(hd[j])

os.chdir(r"C:\Users\davide\Downloads\run1_all_devices_led_365nm_20ns_3V20")
h1.file="0_wave5_all_devices_led_365nm_20ns_3V15.dat"

hd[j].a_wout_baseline,  hd[j].t = h1.read_binary(0)
os.chdir(r"E:")

#%% DATI random trigger  
j=0
hd[j].ns = 4
h1 = READ_DATA(hd[j])

os.chdir(r"C:\Users\davide\Downloads\Random_trigger\run20_v5_turned_off_random_trigger")
h1.file="0_wave5_v5_turned_off_random_trigger.dat"

hd[j].a_wout_baseline,  hd[j].t = h1.read_binary(0)
os.chdir(r"E:")

#%%

hd[j].a = hd[j].a_wout_baseline
#TO DO: save processed data in a dataframe

h2 = FILTER_DATA(hd[j])
h2.l_filtfilt = 7
hd[j].a_filt_wout_baseline = h2.filtfilt(hd[j].a_wout_baseline)
h2.threshold_base = 30
h2.discard_tre = 0
h2.debug_base = False
h2.exclusion_window = 1000 #ns

hd[j].a_filt, mu, cont, uno = h2.baseline(hd[j].a_filt_wout_baseline)   
hd[j].a = h2.sottraggo_baseline(hd[j].a, hd[j].a_wout_baseline, mu)    #sottraggo baseline a dati non filtrati (grezzi)

h2.l_filtfilt = 13
hd[j].a_dfilt = h2.filtfilt(hd[j].a)




#%%
histo_2d(hd[j].t, hd[j].a_filt, x1=0, x2=hd[j].t[0,-1], step1=hd[j].ns, y1=-60, y2=150, step2=1, draw=False, mean=0)
#%%  MULTIFIT

h3 = ANALYZE_DATA(hd[j])
h3.selection=False
h3.isto_plot = True
h3.t_min = 10350
h3.t_max = 10500
h3.range_x1 = -5000
h3.range_x2 = 55500
h3.bin = 200
hd[j].n, hd[j].b, hd[j].charge = h3.isto_carica(hd[j].t, hd[j].a_filt, uno)
#%%
#h3.p0 = [1300, 0, 1000, 1440, 5700, 1100, 1200, 11100] #spe ch5 14/09
h3.p0 = [250, 0, 200, 50, 6000, 50, 10, 12010]
h3.numg = 2
h3.bound1 = -1000
h3.bound2 = 10000
h3.fitdue = False
# bounds=[(-np.inf, np.inf) for i in range((h3.numg)*2+8)] #costrains on last 2 gaussians
# z = (0, popt[-1]*0.2)
# bounds.append(z)
# del z
popt, pcov, SNR, chi2 = h3.multi_fit()
#%%
m=h3.mean_wvf(popt[4], popt[5], mult1=1, mult2=1)
histo_2d(hd[j].t, hd[j].a_filt, x1=0, x2=hd[j].t[0,-1], step1=hd[j].ns, y1=-60, y2=150, step2=1, draw=True, mean=m)



#%%  SIGNAL
j=0

h3 = ANALYZE_DATA(hd[j])
h3.signal_threshold = 3.5
h3.prima_t = 30
h3.dopo_t = 170                
h3.debug_signal = False
h3.lim_amp = 300  # >10 pe is cosmic
h3.windows = 6
h3.time_after_pulse = 1000 
h3.number_events_expected = 800
h3.sig_plot=True
h3.t_min_sig = 0
h3.t_max_sig = 0
h3.bin = 400
h3.range_x1 = -1000
h3.range_x2 = 6000

h3.time_noise = 250
h3.value_noise = -60
charge, events = h3.signal(hd[j].t, hd[j].a_filt, hd[j].a_dfilt)


# TRIGGER RATE
rate = sum(events)/(len(hd[j].a)*len(hd[j].a[0])*4e-9*(80*36))  #Hz #mini arapuca 20 sipm
time_window = 20*1e-6
rate_f = rate*time_window*80*36
print(round(rate,3), 'Hz/mm^2')
print(f'Events per {int(time_window*1e6)} us: {round(rate_f, 3)}')
print('ev', round(sum(events)*1e-3,1),'k')
print("tre:", h3.signal_threshold)
print("integration:", h3.prima_t, h3.dopo_t)
















