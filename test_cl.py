# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 03:32:54 2022

@author: Davide
"""

from read_Data import read_data


h1 = read_data(file_type='binary', file='wave2.dat', header=12, adc_res=1, ns=4, bit=16, iniz='', puntixwvf=0)
a, t, _, _ = h1.__read__()
#save processed data in a dataframe