# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 16:17:55 2020

@author: clair
"""

# import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.linalg as solver
import decimal as dec
import sympy as symp
import json as json
from scipy.signal import butter, lfilter
from scipy.signal import get_window
import time

#%%Set Model Info

######################################################################################
#model info
with open('modelParameters.json') as json_file:
    data = json.load(json_file)
   
for keys, vals in data.items():
    exec(keys + '=vals')

######################################################################################
    
np.set_printoptions(np.set_printoptions(precision=50))
   



#%% State 2 Function
   
def State2(pin_loc, number_modes, number_nodes, M, K, matrix_size):
    # t0 = time.perf_counter()
    spring_loc= np.int(np.round((pin_loc*(number_nodes-1))))
    
#################################### STATE 2 ########################################
    del_K = np.zeros((matrix_size,matrix_size))
    del_K.itemset(((spring_loc)*2, (spring_loc)*2), data["pin_stiffness_change"])
    
    mM2= M
    kK2= K + del_K
    # t0 = time.perf_counter()
    # Calculation of the natural frequencies. 
    eigvals,eigvects = sp.linalg.eigh(kK2,mM2)
    eigvals=np.expand_dims(np.real(eigvals), axis=0)
    # t1 = time.perf_counter()
    # total = (t1-t0)
    FEA_wn = np.sort(np.real(np.squeeze(np.sqrt(eigvals)))) # Natural frequencies, rad/s
    Frequencies = FEA_wn/(2*np.pi) # Natural freq in Hz
    # t1 = time.perf_counter()
    # total = (t1-t0)
    # print(total)
    # print(Frequencies)
    return(Frequencies[0:number_modes])

