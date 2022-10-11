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

hd = [DATI(), DATI(), DATI(), DATI(), DATI()]  #id() per indirizzo memoria, controllare che sia lo stesso altrimenti copia
#%% DATI RANDOM TRIGGER 14/09
j=0
hd[j].ns = 4
h1 = READ_DATA(hd[j])

os.chdir(r"C:\Users\Davide\Documents\DATI\random trigger only mini")
file1="0_wave5_37V00_20ADC_random_trigger.dat"
file2="1_wave6_37V00_20ADC_random_trigger.dat"
files = [file1, file2]
n_file = len(files)
hd[j].a_wout_baseline,  hd[j].t = h1.concatenate(n_file, files)

# h1.file="0_wave2_47V00_20ADC_random_trigger.dat"
# hd[j].a_wout_baseline,  hd[j].t, _, _ = h1.__read__()
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