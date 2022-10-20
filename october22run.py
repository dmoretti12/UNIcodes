# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 15:26:52 2022

@author: Davide
"""
from tutte_le_mie_cazzate import*
plt.rcParams.update({'font.size': 18})

hd = [DATI(), DATI(), DATI(), DATI(), DATI(), DATI()]  #id() per indirizzo memoria, controllare che sia lo stesso altrimenti copia
#%% run0_CATHODE_OFF_48V [500 wvf x 250000 pt = 500 ms x ch5]
j=0
hd[j].ns = 4
h1 = READ_DATA(hd[j])

os.chdir(r"C:\Users\Davide\Downloads\48V_miniXarapuca_DCR\run0_CATHODE_OFF_37V_channelON")
# file1="0_wave5_37V00_20ADC_random_trigger.dat"    #ch5
# file2="1_wave6_37V00_20ADC_random_trigger.dat"   # ch6
# files = [file1, file2]
# n_file = len(files)
# hd[j].a_wout_baseline,  hd[j].t = h1.concatenate(n_file, files)

h1.file="1_wave1_CATHODE_OFF_37V_channelON_LED_OFF.dat"
hd[j].a_wout_baseline,  hd[j].t, _, _ = h1.__read__()
os.chdir(r"C:\Users\Davide\Documents\GitHub\UNIcodes")

#%% SPE
j=1
hd[j].ns = 4
h1 = READ_DATA(hd[j])

os.chdir(r"C:\Users\Davide\Documents\DATI")

h1.file="12_wave2_MiniXArapuca_Wall_48V_MiniXArapuca_37V_Not_Biased_PoF_OFF_LED_4_3V.dat"
hd[j].a_wout_baseline,  hd[j].t, _, _ = h1.__read__()
os.chdir(r"C:\Users\Davide\Documents\GitHub\UNIcodes")

#%%
plot(hd[j].t, hd[j].a_wout_baseline, hd[j].a_wout_baseline, 499)

#%%
hd[j].a = hd[j].a_wout_baseline
#TO DO: save processed data in a dataframe

h2 = FILTER_DATA(hd[j])
h2.l_filtfilt=5
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
plot(hd[j].t, hd[j].a, hd[j].a_filt, 333)
#%%
histo_2d(hd[j].t, hd[j].a_filt, x1=0, x2=hd[j].t[0,-1], step1=hd[j].ns, y1=-50, y2=150, step2=1)
#%%  MULTIFIT

h3 = ANALYZE_DATA(hd[j])
h3.selection=False
h3.isto_plot = False
h3.t_min = 10340
h3.t_max = 10450
h3.range_x1 = -3000
h3.range_x2 = 15000
h3.bin = 200
hd[j].n, hd[j].b, hd[j].charge = h3.isto_carica(hd[j].t, hd[j].a_filt, uno)

#h3.p0 = [1300, 0, 1000, 1440, 5700, 1100, 1200, 11100] #spe ch5 14/09
h3.p0 = [700, 0, 400, 500, 1200, 600, 200, 2400]
h3.numg = 3
h3.bound1 = -1000
h3.bound2 = 20000
h3.fitdue = True
# bounds=[(-np.inf, np.inf) for i in range((h3.numg)*2+8)] #costrains on last 2 gaussians
# z = (0, popt[-1]*0.2)
# bounds.append(z)
# del z
popt, pcov, SNR, chi2 = h3.multi_fit()
print('integrating for', h3.t_max-h3.t_min, 'ns, (', h3.t_min,'-',h3.t_max,')')

#%% signal

h3 = ANALYZE_DATA(hd[j])
h3.signal_threshold = 3
h3.prima_t = 30
h3.dopo_t = 100                
h3.debug_signal = False
h3.lim_amp = 300  # >10 pe is cosmic
h3.windows = 3
h3.time_after_pulse = 1000 
h3.number_events_expected = 800
h3.sig_plot=True
h3.t_min_sig = 0
h3.t_max_sig = 0
h3.bin = 200
h3.range_x1 = -100
h3.range_x2 = 6000
charge, events = h3.signal(hd[j].t, hd[j].a_filt, hd[j].a_dfilt)


# TRIGGER RATE
rate = sum(events)/(len(hd[j].a)*len(hd[j].a[0])*4e-9*(20*36))  #Hz #mini arapuca 20 sipm
print(round(rate,3), 'Hz/mm^2')
print('ev', round(sum(events)*1e-3,1),'k')
print("tre:", h3.signal_threshold)
print("integration:",h3.dopo_t-h3.prima_t, 'ns, (',h3.prima_t, 'ns before peak-', h3.dopo_t,'ns after peak)')

