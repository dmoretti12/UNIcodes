# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 12:01:34 2022

@author: Davide
"""

from tutte_le_mie_cazzate import*

hd = [DATI()]
#%% DATI RANDOM TRIGGER 09-11 argon2x2
j=0
hd[j].ns = 4
h1 = READ_DATA(hd[j])

# os.chdir(r"/home/dune/Documents/ADC_data/coldbox_data/November2022run/2022_11_09/random_trigger/argon2x2/run0_miniArapuca_37V_only_random_trigger_PoF_off")
# h1.file="0_wave1_miniArapuca_37V_only_random_trigger_PoF_off.dat"

# os.chdir(r"/home/dune/Documents/ADC_data/coldbox_data/November2022run/2022_11_09/random_trigger/argon2x2/run1_miniArapuca_37V_only_random_trigger_PoF_on")
# h1.file="0_wave1_miniArapuca_37V_only_random_trigger_PoF_on.dat"

os.chdir(r"/home/dune/Documents/ADC_data/coldbox_data/November2022run/2022_11_09/random_trigger/argon2x2/run2_miniArapuca_37V_only_random_trigger_PoF_A_laser_1")
h1.file="0_wave1_miniArapuca_37V_only_random_trigger_PoF_A_laser_1.dat"

hd[j].a_wout_baseline,  hd[j].t = h1.__read__()
os.chdir(r"/home/dune/Documents/apc_python_spyder")
#%% DATI LED 08-11
j=0
hd[j].ns = 4
h1 = READ_DATA(hd[j])

os.chdir(r"/home/dune/Documents/apc_python_spyder/2022_11_07/LED/Argon_2x2/run0_LED_3dot3")

h1.file="0_wave1_LED_3dot3.dat"
hd[j].a_wout_baseline,  hd[j].t = h1.__read__()
os.chdir(r"/home/dune/Documents/apc_python_spyder")

#%% DATI LED 09-11 argon2x2
j=0
hd[j].ns = 4
h1 = READ_DATA(hd[j])

os.chdir(r"/home/dune/Documents/ADC_data/coldbox_data/November2022run/2022_11_09/LED/argon2x2/run0_miniArapuca_37V_only_3V3_20ns_365nm")
h1.file="0_wave1_miniArapuca_37V_only_3V3_20ns_365nm.dat"

hd[j].a_wout_baseline,  hd[j].t = h1.__read__()
os.chdir(r"/home/dune/Documents/apc_python_spyder")

#%% DATI LED 10-11 argon2x2
j=0
hd[j].ns = 4
h1 = READ_DATA(hd[j])

os.chdir(r"/home/dune/Documents/ADC_data/coldbox_data/November2022run/2022_11_10/LED/argon2x2/run0_miniArapuca_37V_only_3V45_20ns_365nm")
h1.file="0_wave2_miniArapuca_37V_only_3V45_20ns_365nm.dat"

hd[j].a_wout_baseline,  hd[j].t = h1.__read__()
os.chdir(r"/home/dune/Documents/apc_python_spyder")



#%% DATI RANDOM TRIGGER 10-11 argon2x2


j=0
hd[j].ns = 4
h1 = READ_DATA(hd[j])

os.chdir(r"/home/dune/Documents/ADC_data/coldbox_data/November2022run/2022_11_10/random_trigger/argon2x2/run8_miniArapuca_37V_random_trigger_PoF_all_on")
h1.file="0_wave2_miniArapuca_37V_random_trigger_PoF_all_on.dat"
# h1.header=0
# os.chdir(r"/home/dune/Documents/apc_python_spyder")
# h1.file="october_hfiltered.dat"
# filez = np.fromfile(h1.file, dtype=np.float64)
hd[j].a_wout_baseline,  hd[j].t = h1.read_binary(0)
os.chdir(r"/home/dune/Documents/apc_python_spyder")

#%% DATI LED 11-11
j=0
hd[j].ns = 4
h1 = READ_DATA(hd[j])

os.chdir(r"C:\Users\davide")
h1.header=0
h1.file="output.dat"
hd[j].a_wout_baseline,  hd[j].t = h1.read_binary(0)
os.chdir(r"E:")
#%%

hd[j].a = hd[j].a_wout_baseline
#TO DO: save processed data in a dataframe

h2 = FILTER_DATA(hd[j])
h2.l_filtfilt = 3
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
h3.t_max = 10600
h3.range_x1 = -3500
h3.range_x2 = 12500
h3.bin = 200
hd[j].n, hd[j].b, hd[j].charge = h3.isto_carica(hd[j].t, hd[j].a_filt, uno)
#%%
#h3.p0 = [1300, 0, 1000, 1440, 5700, 1100, 1200, 11100] #spe ch5 14/09
h3.p0 = [700, 0, 600, 300, 1000, 600, 100, 3010]
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
h3.signal_threshold = 10
h3.prima_t = 50
h3.dopo_t = 200                
h3.debug_signal = False
h3.lim_amp = 300  # >10 pe is cosmic
h3.windows = 6
h3.time_after_pulse = 1000 
h3.number_events_expected = 800
h3.sig_plot=True
h3.t_min_sig = 10330
h3.t_max_sig = 10650
h3.bin = 200
h3.range_x1 = -3000
h3.range_x2 = 30000

h3.time_noise = 250
h3.value_noise = -40
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



























