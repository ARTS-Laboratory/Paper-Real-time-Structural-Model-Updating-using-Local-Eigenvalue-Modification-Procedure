# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 09:08:57 2020

@author: clair
"""
import numpy
import numpy as np
# import matplotlib.pyplot as plt
import decimal as dec
import sympy as sp
import scipy.linalg as solve
import scipy as sci
import math
import json as json
from numpy.ctypeslib import ndpointer
import ctypes
from ctypes import *
import time
import secular_equation_functionV4 as solve_sec
#%%Set info from json files
np.set_printoptions(np.set_printoptions(precision=50))

#model info
with open('modelParameters.json') as json_file:
    data = json.load(json_file)
   
for keys, vals in data.items():
    exec(keys + '=vals')

#%%Functions
def state2(spring_node, modes_number, nodes_number):
    
    del_K=np.zeros((nodes_number*2,nodes_number*2))
    del_K.itemset(((spring_node*2), (spring_node*2)), data["pin_stiffness_change"])

    return(del_K)
    
    
###############################################################################
    
#Make this a function!
#DEFORMED STATE
def LEMP(pin_loc, number_modes, number_nodes, a, U1, wn_squared_state1):
    # t0 = time.perf_counter()
    'stage1'
    spring_loc= np.int(np.round((pin_loc*(number_nodes-1))))
    # t1 = time.perf_counter()
    # total = (t1-t0)
    # t0 = time.perf_counter()
    'stage 2'
    #DEFORMED STATE    
    del_K= state2(spring_node=spring_loc, modes_number=number_modes, nodes_number=number_nodes)
    # t1 = time.perf_counter()
    # total = (t1-t0)
    #Spectral Decomposition of del_K
    # t0 = time.perf_counter()
    'stage 3'
    #alpha values and tie vectors (T)
    # alpha, T = np.linalg.eig(del_K)
    alpha = del_K.diagonal()
    
    T = np.zeros((len(alpha), len(alpha)), int) # Create matrix with only 0
    np.fill_diagonal(T, 1) # fill diagonal with 1
    # t1 = time.perf_counter()
    # total = (t1-t0)
    
    # t0 = time.perf_counter()
    'stage 4'
    #only connecting nodes contribute to state changes so modifications necessary...
    # v_matrix= U1.T @ T
    # v_matrix = np.matmul(U1.T, T)
    v_matrix= U1.T
    v_vect_contribute= v_matrix[:,((spring_loc*2))]
    # t1 = time.perf_counter()
    # total = (t1-t0)
    #define how many modes to include i.e. set truncation
    # t0 = time.perf_counter()
    'stage 5'
    v_vect=v_vect_contribute[0:number_modes]
    # wn_squared_state1=wn_squared[:,0:number_modes]
    # t1 = time.perf_counter()
    # total = (t1-t0)
    
    #equation based on Peter's notes
    #wn_squared_state1= np.reshape(wn_squared_state1, (np.size(wn_squared_state1, axis=1), 1))
    
    #creates an equation using summation
    #sp.init_printing(use_unicode=True)
    '''
    To use sympy function solveset
    '''
    # summation = 0
    # vr, wr, o2 = sp.symbols('vr wr o2')
    # expr = ((vr * vr)/ (wr-o2))
            
    # for r in range (number_modes):
    #     v_r= v_vect[r]
    #     wn_r= wn_squared_state1[0,r]
                
    #     part = expr.subs(vr, v_r).subs(wr, wn_r)
    #     print(part)
        
    #     summation= summation + part
    #     print(summation)
    
    #     eq= sp.Eq(summation, -1/data["pin_equivalent_stiffness"])
    #     print(eq)
           
    # sol= sp.solveset(eq, o2)
    # print(sol)
    
    
     
    '''
    To use secular solver written in python
    '''
    
    # t0 = time.perf_counter()
    'stage 6'
    v = v_vect**2
    # print(v)
    mode = number_modes
    d = wn_squared_state1
    p= 1/data["pin_equivalent_stiffness"]
    
    
    # sol= solve_sec.secular_eq_solver(mode,d,v,p)
    
    '''
    To use secular solver written in C
    '''
    d = d.flatten()
    d = list(d)
    v= list(v)
    l_d = len(d)
    l_v = len(v)
    # t1 = time.perf_counter()
    # total = (t1-t0)

    # import C file
    # so_file = "C:/Users/OGUNNIYI/Documents/research/LEMP C version/secular_solver_func_v3.so"
    # my_func = CDLL(so_file)
    # a = my_func
   
    # t0 = time.perf_counter()
    'stage 7'
    a.secular.argtypes=(POINTER(c_float),POINTER(c_float),c_float, c_int, c_int)
    a.secular.restype=ndpointer(dtype=c_double,
                          shape=(l_d,))
    sol = a.secular((c_float * l_d)(*d), (c_float * l_v)(*v), p, mode, l_d)
    # t1 = time.perf_counter()
    # total = (t1-t0)
    
    # print(total)
    modes= list(sol)
    'stage 8'
    # t0 = time.perf_counter()
    new_omegas_squared= np.zeros(len(modes))
    new_omegas= np.zeros(len(modes))
    new_frequencies= np.zeros(len(modes))
    
    for k in range (len(modes)):
        new_omegas_squared[k]= modes[k]
        new_omegas[k]= (np.sqrt(new_omegas_squared[k]))
        new_frequencies[k] = new_omegas[k]/(2*np.pi) # Natural freq in Hz
    new_frequencies.sort()
    # t1 = time.perf_counter()
    # total = (t1-t0)
    # print(total)
    return(new_frequencies)
    
