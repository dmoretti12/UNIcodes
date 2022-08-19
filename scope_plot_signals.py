    #!/usr/bin/env python
# -*- coding: utf-8 -*-

#quick plotter for scope data
# give file names first and plot name last
import matplotlib.pyplot as plt
import sys
import numpy as np
import bisect
import math
import scipy.stats as stats
import math
from scipy.stats import norm

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (12, 3.5),
         'axes.labelsize': '23',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'18',
         'ytick.labelsize':'15'}
plt.rcParams.update(params)

def main():
  files = []
  print(sys.argv)

  i=0
  style = ['r','c','b','g','y','o','m']
  plotname = sys.argv[-1]

  
  #labels = ['input signal','koheron','current']
  #labels=['ch2','ch2 inx10','receiver out']
  #labels = ['cathode x100','anode x100','input signal at LED']
  #labels = ['10 mV', '20 mV', '30 mV', '40 mV', '50 mV', '70 mV', '100 mV']
  #labels = ['Receiver Output']
  #labels = ['Input X200','ACR1 out warm','ACR1 out cold']
  #labels = ['LED pulser']
  labels = ['50 uV, Argon4 board, AMP LDO','50 uV, Argon4 board, old LDO']#, '75 uV, ACR1 board, old LDO', '50 uV, Rev 3.0 board, old LDO']
  #labels = ['A2 out']

  for i, fname in enumerate(sys.argv[1:]):
      if fname.find('txt')<0: continue
      print("Opening file ",fname)
      if len(labels)>= i:
        mylabel = labels[i]
      else: 
        mylabel = '--'
      time1 = []
      volt1 = []
      ff1 = open(fname).readlines()[10:]
      sep = ','
      if ff1[0].find(',')<0:
        sep = ' '
      for l in ff1:
        l=l.strip("\n\r")
        #if file doesn't have time info:
        if len(l.split(sep)) > 1:
          time1.append(float(l.split(sep)[0])/1000000)  #time in us
          if fname.find('input')>0:volt1.append(float(l.split(sep)[1])*1)
          else:volt1.append(float(l.split(sep)[1])*1)
        else:
          if len(time1)<1: 
            time1.append(0)
            volt1.append(float(l.split(',')[0]))
          else: 
            time1.append(0.2*pow(10,-6)+time1[-1])  #time in us
            volt1.append(float(l.split(',')[0]))
        #if i==0 or i==1:
        #  #volt1.append(float(l.split()[1])*10)  
        #  volt1.append(float(l.split(',')[1]))  
        #else : volt1.append(float(l.split()[1]))
    
      print( len(time1))
      print( len(volt1))
      pts = 500
      mean_base = np.sum(volt1[0:pts])
      mean_base /= pts
      baseline = volt1[0:pts]
      mean = np.mean(baseline)
      std = np.std(baseline)
      variance = np.sqrt(std)
      print(mean)
      print(std)
      plt.plot(baseline, stats.norm.pdf(baseline, mu, sigma), color = 'red')
      plt.xlabel('Data points')
      plt.ylabel('Probability Density')  
      #st = 0 #I commented from here down, remove coments to get back the previous code
      #end = len(time1)
      #st = round(len(time1)/5*2) #from here to
      #end = round(len(time1)/2)
      #st =len(time1)/8
      #end= len(time1)/3*2
      #st=bisect.bisect_left(time1, -0.1/1000000)
      #end=bisect.bisect_left(time1, 0.4/1000000) #here it was already comented
      #print(st)
      #print(end)
      #plt.figure(figsize=(12,4)) #also already comented
      #plt.plot(time1[st:end],np.array(volt1[st:end]- mean_base)/np.max(volt1 - mean_base),style[i],label=mylabel)
      
    
    #plt.title('50 uV signal, Argon4 board, new and old LDO comparison')
    #plt.legend(loc='upper right')
    #plt.ylabel('Reciever [mV]')
    #plt.xlabel('Time [us]')
    #plt.tight_layout()
    
    #plt.xlim(time[100],time[len(time)-100])
    #plt.ylim(min(volt1)*0.95, max(volt1)*1.05)
    #plt.ylim(-50, 1000)
    #plt.grid()
      plt.show()
#plt.savefig('AMP_LDO_vs_old_LDO_50uV.png', format='png', dpi= 400)
#plt.savefig('ACR1_waveform_cold',bbox_inches='tight')

################################################  
if __name__ == '__main__':
  main()
