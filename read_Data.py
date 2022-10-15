# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 02:32:49 2022

@author: Davide
"""

import numpy as np
import os
import glob
import pandas as pd


class READ_DATA:
    
    def __init__(self, classe_dati):
        
        self.hd = classe_dati
        
        self.type = 'binary'
        self.file = 'cc.dat'
        
        self.header = 12
        self.adc_res = 1
        self.bit = 16
        
        self.iniz = 'C1X*.txt'
        
        self.puntixwvf = 938
     
 #################################################       
        
    def read_binary(self):
        
        ns = self.hd.ns
        if self.bit==16:
            filez = np.fromfile(self.file, dtype=np.int16)
            print('len file:', len(filez))
        filez32 = np.fromfile(self.file, dtype=np.int32, count=2)  # header=6 invece di 12
        # if self.bit==32:
        #     filez = filez32
        found = False
        pointsxwvf = int(filez32[0]/2)  # 2 bin per dato
        print("points:", pointsxwvf)                   #debug
        # for i in range(1, pointsxwvf):
        #     if filez[i]==filez[0]:
        #         pointsxwvf = int(i)
        #         # print('trovato', i)      #debug
        #         found = True
        #         break
        wvf = int(len(filez)/pointsxwvf)
        print('wvfs:', wvf)
        points = int(pointsxwvf-self.header)
        # print(wvf, points)               #debug
        if found==False:
        # pointsxwvf32 = int(pointsxwvf/2) 
            wvf = int(len(filez)/pointsxwvf)
            points = pointsxwvf-self.header   # tolgo l'header x ogni wvf

        a = np.zeros((wvf, points))
        adc = np.zeros((wvf, points))
        t = np.zeros((wvf, points))
        # time = np.empty(wvf)

        for i in range(0, wvf):

            adc[i] = filez[(i*pointsxwvf)+self.header:pointsxwvf*(i+1)]
            t[i] = np.arange(0, ns*points, ns)  # 250 MHz = 4e-03 us
            # time[i] = file[pointsxwvf32*i+5]
        a = adc * self.adc_res  #valori in volt o mV
        
        # freq 62.5 Mhz se a è 250 Mhz
        out = int(len(a[0])/4)
        
        t62 = np.empty((len(a), out))
        a62 = np.empty((len(a), out))

        for i in range(0, out):

            t62[:, i] = t[:, 4*i]
            a62[:, i] = a[:, 4*i]

        return a, t, a62, t62#, adc , time
        
#########################################################################
        
    def concatenate(self, n_file, files):
        
        for i in range(n_file):
            
            self.file = files[i]
            ai, _, _, _  = self.__read__()
            
            if i==0:
                a = np.zeros((len(ai)*n_file, len(ai[0])))
                a[0:len(ai)] = ai
            
            else: a[len(ai)*i:(i+1)*len(ai)] = ai
                # a = np.concatenate((a, ai), axis=0)#, out=a)
        
        t = np.zeros((len(a),len(a[0])))
        
        for i in range(len(t)):
            t[i] = np.arange(0, self.hd.ns*len(ai[0]), self.hd.ns)
            
        return a, t

#########################################################################

    def read_many_txt(self):

        os.chdir(self.file)

        filenames = [i for i in glob.glob(self.iniz)]
        
        #df = pd.read_csv('C1XArapuca_CRT_250MHz_00000.txt', sep=" ", header=5, index_col = None, encoding = 'unicode_escape')
        df = [pd.read_csv(f, sep=",", header=self.header, engine='python', index_col=None, encoding='unicode_escape')
              for f in filenames]
        
        t = np.empty((len(df), len(df[0]))) 
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
        
        t62 = np.empty((len(df), out))  # freq 62.5 Mhz
        a62 = np.empty((len(df), out))

        for i in range(0, out):

            t62[:, i] = t[:, 4*i]
            a62[:, i] = a[:, 4*i]

        return df, a, t, a62, t62
        
############################################################################  
        

    def read_long_txt(self):
        
        ns = self.hd.ns
        data = np.genfromtxt(self.file, skip_header=self.header)   #PRENDE TANTO TEMPO
        points = self.puntixwvf
        wvf = int(len(data)/points)  

        a = np.empty((wvf, points))
        t = np.empty((wvf, points))

        for i in range(0, wvf):

            a[i] = data[points*i:points*(i+1)]
            t[i] = np.arange(0, ns*points, ns)

        # freq 62.5 Mhz
        out = int(len(a[0])/4)

        t62 = np.empty((len(a), out))
        a62 = np.empty((len(a), out))

        for i in range(0, out):

            t62[:, i] = t[:, 4*i]
            a62[:, i] = a[:, 4*i]

        return a, t, a62, t62
        
        
#######################################################################


    def read_long_txt_time(self):

        f = np.genfromtxt(self.file, skip_header=self.header)
        tt = f[:,0]
        aa = f[:,1]
        points = self.puntixwvf
        wvf = int(len(tt)/points)

        a = np.empty((wvf, points))
        t = np.empty((wvf, points))

        for i in range(0, wvf):

            a[i] = aa[points*i:points*(i+1)]
            t[i] = tt[points*i:points*(i+1)]

        # freq 62.5 Mhz
        out = int(len(a[0])/4)

        t62 = np.empty((len(a), out))
        a62 = np.empty((len(a), out))

        for i in range(0, out):

            t62[:, i] = t[:, 4*i]
            a62[:, i] = a[:, 4*i]

        return a, t, a62, t62

############################################################################à

    def __read__(self):
        if self.type=='binary':
            f = READ_DATA.read_binary(self)
        if self.type=='many_txt':
            f = READ_DATA.read_many_txt(self)
        if self.type=='long_txt':
            f = READ_DATA.read_long_txt(self)
        if self.type=='long_txt_time':
            f = READ_DATA.read_long_txt_time(self)
        return f

     
    