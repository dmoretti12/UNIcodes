# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 21:36:12 2022

@author: Davide
"""

"""
sept22 only mini
run1_47V00_20ADC_led_5V_50ns 
random trigger only mini

august22
run0_led_7V_50ns
run0_light_leakage_test_led_flange_covered

july22
July22run(minixarapuca)

june22
run0_35V00_100ADC_externa_100ns_1kHz_3V

CONTROLLA CH (5 o 6) e lunghezze e num wvf
"""

from tutte_le_mie_cazzate import*
plt.rcParams.update({'font.size': 18})

hd = [DATI(), DATI(), DATI(), DATI(), DATI()]  #id() per indirizzo memoria, controllare che sia lo stesso altrimenti copia
#%% DATI RANDOM TRIGGER 14/09  [500 wvf x 250000 pt x ch5]
j=0
hd[j].ns = 4
h1 = READ_DATA(hd[j])

os.chdir(r"C:\Users\Davide\Documents\DATI\random trigger only mini")
# file1="0_wave5_37V00_20ADC_random_trigger.dat"    #ch5
# file2="1_wave6_37V00_20ADC_random_trigger.dat"   # ch6
# files = [file1, file2]
# n_file = len(files)
# hd[j].a_wout_baseline,  hd[j].t = h1.concatenate(n_file, files)

h1.file="0_wave5_37V00_20ADC_random_trigger.dat"
hd[j].a_wout_baseline,  hd[j].t, _, _ = h1.__read__()
os.chdir(r"C:\Users\Davide\Documents\GitHub\UNIcodes")

#%% DATI LED 14/09
j=1
hd[j].ns = 4
g1 = READ_DATA(hd[j])

os.chdir(r"C:\Users\Davide\Documents\DATI\run1_47V00_20ADC_led_5V_50ns ")
file1="0_wave5_47V00_20ADC_led_5V_50ns.dat"
file2="0_wave6_47V00_20ADC_led_5V_50ns.dat"
files = [file1, file2]
n_file = len(files)

hd[j].a_wout_baseline,  hd[j].t = g1.concatenate(n_file, files)
# h1.file="0_wave5_37V00_20ADC_random_trigger.dat"
# hd[j].a_wout_baseline,  hd[j].t, _, _ = h1.__read__()
os.chdir(r"C:\Users\Davide\Documents\GitHub\UNIcodes")

#%% DATI RANDOM TRIGGER august
j=2
hd[j].ns = 4
h1 = READ_DATA(hd[j])

os.chdir(r"C:\Users\Davide\Documents\DATI\run0_light_leakage_test_led_flange_covered")
file1="0_wave1_light_leakage_test_led_flange_covered.dat"
file2="0_wave2_light_leakage_test_led_flange_covered.dat"
files = [file1, file2]
n_file = len(files)
hd[j].a_wout_baseline,  hd[j].t = h1.concatenate(n_file, files)

# h1.file="0_wave2_47V00_20ADC_random_trigger.dat"
# hd[j].a_wout_baseline,  hd[j].t, _, _ = h1.__read__()
os.chdir(r"C:\Users\Davide\Documents\GitHub\UNIcodes")

#%% DATI LED august
j=3
hd[j].ns = 4
g1 = READ_DATA(hd[j])

os.chdir(r"C:\Users\Davide\Documents\DATI\run0_led_7V_50ns")
file1="0_wave1_led_7V_50ns.dat"
file2="0_wave2_led_7V_50ns.dat"
files = [file1, file2]
n_file = len(files)

hd[j].a_wout_baseline,  hd[j].t = g1.concatenate(n_file, files)
# h1.file="0_wave5_37V00_20ADC_random_trigger.dat"
# hd[j].a_wout_baseline,  hd[j].t, _, _ = h1.__read__()
os.chdir(r"C:\Users\Davide\Documents\GitHub\UNIcodes")

#%% DATI LED july
j=4
hd[j].ns = 4
g1 = READ_DATA(hd[j])

os.chdir(r"C:\Users\Davide\Documents\DATI\July22run(minixarapuca)")
file1="wave1.dat"
file2="wave2.dat"
file3="wave3.dat"
file4="wave4.dat"
file5="wave5.dat"
file6="wave6.dat"
files = [file1, file2, file3, file4, file5, file6]
n_file = len(files)

hd[j].a_wout_baseline,  hd[j].t = g1.concatenate(n_file, files)
# h1.file="0_wave5_37V00_20ADC_random_trigger.dat"
# hd[j].a_wout_baseline,  hd[j].t, _, _ = h1.__read__()
os.chdir(r"C:\Users\Davide\Documents\GitHub\UNIcodes")

#%% DATI LED june
j=5
hd[j].ns = 4
g1 = READ_DATA(hd[j])

os.chdir(r"C:\Users\Davide\Documents\DATI\run0_35V00_100ADC_externa_100ns_1kHz_3V")
file1="0_wave0_35V00_100ADC_externa_100ns_1kHz_3V.dat"
file2="0_wave1_35V00_100ADC_externa_100ns_1kHz_3V.dat"
file3="0_wave2_35V00_100ADC_externa_100ns_1kHz_3V.dat"
file4="0_wave3_35V00_100ADC_externa_100ns_1kHz_3V.dat"
files = [file1, file2, file3, file4]
n_file = len(files)

hd[j].a_wout_baseline,  hd[j].t = g1.concatenate(n_file, files)
# h1.file="0_wave5_37V00_20ADC_random_trigger.dat"
# hd[j].a_wout_baseline,  hd[j].t, _, _ = h1.__read__()
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

h2.l_mov_av = 35
h2.f_denoising = 100
# hd[j].a_dfilt = h2.mov_av(hd[j].a)
hd[j].a_dfilt = h2.denoising(hd[j].a)

#%%

plot(hd[j].t, hd[j].a, hd[j].a_dfilt, 0)
#%%
from functions import *
n, b = np.histogram(hd[j].a[hd[j].a_dfilt<20]-hd[j].a_dfilt[hd[j].a_dfilt<20], 200, range=[-100, 200])
b1 = h2.center_bins(b)
l = (b[2]-b[1])/2
# mu, std = norm.fit(hd[j].a)
popt, pcov = curve_fit(gauss, b1, n, p0=[4e6, 0, 20], maxfev=100000)
x = np.linspace(min(b1)-5*l, max(b1)+5*l, 1000)
y = gauss(x, *popt)
tx = ['Mean = '+str(round(popt[1], 3))+', std = '+str(round(popt[2], 3))]

#%%
plt.figure(figsize=[10, 6])
ax=plt.gca()
ax.tick_params(bottom=True, top=True, left=True, right=True)
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
plt.plot(b1, n, '-', drawstyle='steps')
plt.plot(x, y, 'r-', label=tx)
plt.ylabel('Counts')
plt.xlabel('Amplitude [ADC]')
plt.legend(loc='upper right')
plt.show()





