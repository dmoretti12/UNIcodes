# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 12:01:34 2022

@author: Davide
"""

from allmyimports import*

hd = [DATI(), DATI(), DATI(), DATI()]  #id() per indirizzo memoria, controllare che sia lo stesso altrimenti copia
#%% DATI RANDOM TRIGGER 14/09
j=0
hd[j].ns = 4
h1 = READ_DATA(hd[j])

os.chdir(r"C:\Users\Davide\Documents\DATI\random trigger x e mini")
# file1="0_wave6_37V00_20ADC_random_trigger.dat"
# file2="1_wave6_37V00_20ADC_random_trigger.dat"
# files = [file1, file2]
# n_file = len(files)
# hd[j].a_wout_baseline,  hd[j].t = h1.concatenate(n_file, files)

h1.file="0_wave2_47V00_20ADC_random_trigger.dat"
hd[j].a_wout_baseline,  hd[j].t, _, _ = h1.__read__()
os.chdir(r"C:\Users\Davide\Documents\GitHub\UNIcodes")

#%% DATI RANDOM TRIGGER 16/09
j=0
hd[j].ns = 4
h1 = READ_DATA(hd[j])

# os.chdir(r"C:\Users\Davide\Documents\DATI\run0_47V00_20ADC_random_trigger_digital_on_slow_control_working")
# file1="0_wave5_47V00_20ADC_random_trigger.dat"
# file2="1_wave5_47V00_20ADC_random_trigger.dat"
os.chdir(r"C:\Users\Davide\Documents\DATI\run1_47V00_20ADC_random_trigger_digital_on_slow_control_off")
file1="0_wave6_47V00_20ADC_random_trigger_digital_on_slow_control_off.dat"
file2="1_wave6_47V00_20ADC_random_trigger_digital_on_slow_control_off.dat"

files = [file1, file2]
n_file = len(files)
hd[j].a_wout_baseline,  hd[j].t = h1.concatenate(n_file, files)

# h1.file="0_wave6_47V00_20ADC_random_trigger.dat"
# hd[j].a_wout_baseline,  hd[j].t, _, _ = h1.__read__()
os.chdir(r"C:\Users\Davide\Documents\GitHub\UNIcodes")
#%% DATI LED 14/09
j=0
hd[j].ns = 4
g1 = READ_DATA(hd[j])

os.chdir("H:\Il mio Drive\WORK\\run0_47V00_20ADC_led_7V_50ns")
file1="0_wave6_47V00_20ADC_led_7V_50ns.dat"
file2="1_wave6_47V00_20ADC_led_7V_50ns.dat"
file3="2_wave6_47V00_20ADC_led_7V_50ns.dat"
file4="3_wave6_47V00_20ADC_led_7V_50ns.dat"
files = [file1, file2, file3, file4]
n_file = len(files)

hd[j].a_wout_baseline,  hd[j].t = g1.concatenate(n_file, files)
# h1.file="0_wave5_37V00_20ADC_random_trigger.dat"
# hd[j].a_wout_baseline,  hd[j].t, _, _ = h1.__read__()
os.chdir(r"C:\Users\Davide\Documents\GitHub\UNIcodes")

#%% DATI LED 16/09
j=0
hd[j].ns = 4
h1 = READ_DATA(hd[j])

os.chdir(r"C:\Users\Davide\Documents\DATI\run0_47V00_20ADC_external_trigger_3V3_20ns")
h1.file="0_wave2_47V00_20ADC_external_trigger_3V3_20ns.dat"

hd[j].a_wout_baseline,  hd[j].t, _, _ = h1.__read__()
os.chdir(r"C:\Users\Davide\Documents\GitHub\UNIcodes")

#%%

hd[j].a = hd[j].a_wout_baseline
#TO DO: save processed data in a dataframe

h2 = FILTER_DATA(hd[j])
#h2.l_mov_av = 12
hd[j].a_filt_wout_baseline = h2.mov_av(hd[j].a_wout_baseline)
h2.threshold_base = 30
h2.discard_tre = 0
h2.debug_base = False
h2.exclusion_window = 1000 #ns

hd[j].a_filt, mu, cont, uno = h2.baseline(hd[j].a_filt_wout_baseline)   
hd[j].a = h2.sottraggo_baseline(hd[j].a, hd[j].a_wout_baseline, mu)    #sottraggo baseline a dati non filtrati (grezzi)

h2.l_mov_av = 17
hd[j].a_dfilt = h2.mov_av(hd[j].a)


#%%
histo_2d(hd[j].t, hd[j].a_filt, x1=0, x2=hd[j].t[0,-1], step1=hd[j].ns, y1=-60, y2=150, step2=1)
#%%  MULTIFIT

h3 = ANALYZE_DATA(hd[j])
h3.selection=False
h3.isto_plot = True
h3.t_min = 10380
h3.t_max = 10680
h3.range_x1 = -5000
h3.range_x2 = 55000
h3.bin = 200
hd[j].n, hd[j].b, hd[j].charge = h3.isto_carica(hd[j].t, hd[j].a_filt, uno)
#%%
#h3.p0 = [1300, 0, 1000, 1440, 5700, 1100, 1200, 11100] #spe ch5 14/09
h3.p0 = [1600, 0, 400, 940, 2000, 600, 600, 3510]
h3.numg = 5
h3.bound1 = -100
h3.bound2 = 20000
h3.fitdue = True
# bounds=[(-np.inf, np.inf) for i in range((h3.numg)*2+8)] #costrains on last 2 gaussians
# z = (0, popt[-1]*0.2)
# bounds.append(z)
# del z
popt, pcov, SNR, chi2 = h3.multi_fit()





#%%  SIGNAL
j=0

h3 = ANALYZE_DATA(hd[j])
h3.signal_threshold = 50    
h3.prima_t = 30
h3.dopo_t = 210                
h3.debug_signal = False
h3.lim_amp = 300  # >10 pe is cosmic
h3.windows = 6
h3.time_after_pulse = 1000 
h3.number_events_expected = 800
h3.sig_plot=True
h3.t_min_sig = 0
h3.t_max_sig = 0
h3.bin = 200
h3.range_x1 = -3000
h3.range_x2 = 55000
charge, events = h3.signal(hd[j].t, hd[j].a_filt, hd[j].a_dfilt)


# TRIGGER RATE
rate = sum(events)/(len(hd[j].a)*len(hd[j].a[0])*4e-9*(20*36))  #Hz #mini arapuca 20 sipm
print(round(rate,3), 'Hz/mm^2')
print('ev', round(sum(events)*1e-3,1),'k')
print("tre:", h3.signal_threshold)
print("integration:", h3.prima_t, h3.dopo_t)



























