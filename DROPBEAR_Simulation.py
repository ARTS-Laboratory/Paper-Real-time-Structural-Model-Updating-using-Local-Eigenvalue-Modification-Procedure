# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 09:42:55 2020

@author: clair
"""

#%% import modules
# import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.linalg as solver
import decimal as dec
import sympy as symp
import json
import scipy.linalg as solve
from numpy.ctypeslib import ndpointer
import ctypes
from ctypes import *
from scipy.signal import butter, lfilter
from scipy.signal import get_window
import time
import math

#import LEMP_procedure as LEMP
import LEMP_simplified_procedure as LEMP
import GeneralEigenvalue_stiffness as GE

# plt.rcParams["font.family"] = "Times New Roman"

# plt.close('all')
# cc = plt.rcParams['axes.prop_cycle'].by_key()['color'] 

#%% Model Info



time_data_particle_6 = np.zeros(491)
Trac = np.zeros(491)
solver_time = np.zeros(491)
SNR_particle_6 = np.zeros(491)
MAE_particle_6 = np.zeros(491)
for j in range(10,501,10):
    
    nodes=j
    modes=4
    
    #%% Load Data
    
    # read file
    with open('data_6_with_FFT.json', 'r') as myfile:
        data=myfile.read()
    
    # parse file
    obj = json.loads(data)
    
    #%% function definitions
    
    # Filters for the accelerometer time-series data
    # May not use these filters
    
    #%% Set Parameters
    
    # set the tunning parameters
    
    particle_number = 3 # the number of particles that will be tested each time 
    noise_level = 0.0 # the level of noise in the system. 
    est_pin_loc = 0.5 # set the initial position of the pin at the center
    step_pin_standard_dev = 0.02 # 1/3*particle_number # the standard dev of the pin PDF in terms of pins
    # Previously, 44 nodes with 1 standard dev is equivalent to ~0.02 standard deviation
    freq_mode = 4 # total number of fundamental frequencies to estimate
    
    beam_length = 0.3525 # need this value for the graph
    
    
    #%% Loading saved data and information
    
    measured_pin_location = np.asarray(obj['measured_pin_location'])
    measured_time = np.asarray(obj['measured_pin_location_tt'])
    
    acceleration_time = np.asarray(obj['time_acceleration_data'])
    acceleration = np.asarray(obj['acceleration_data'])
    fft_length = obj['FFT_data_length']
    
    Fs = obj['accelerometer_sample_rate']
    fs=Fs
    Ts = 1/Fs
    datasize = len(acceleration_time)
    
    #K = np.arange(fft_length)
    #freq = K/fft_length*Fs
    #freqs = freq[range(fft_length//2)]
    freqs = np.fft.rfftfreq(fft_length, d=Ts) # or simply use this function
    
    #time_average_cycle = obj['time_average_cycle']
    time_average_cycle = 0.005 # 5 ms is found to be the average cycle from the published paper
    Ff = 1/time_average_cycle # frequency of loop/cycle
    
    FFTstep = int(Fs/Ff)
    startT = obj['time_frequency_calculations_start'][0]
    startFFT = np.argmin(abs(acceleration_time - startT)) 
    delay_time = 0.002 # 2ms delay for FEA computation 
    # (apply delay only for the simulation, not for the actual experiement)
    
    
    #%% Preprocessing - windowing and filtering
    
    # may or may not apply filters
    # if filter is applied, the delay will be larger
    # we can simply apply cut-off filters
    #accelerationF = butter_bandpass_filter(acceleration, 40, 130, Fs, order=3) # with filter
    accelerationF = acceleration # without filter
    
    # cut-off frequencies
    # Range of 1st Freq (35~220)
    min1 = np.argmin(np.abs(freqs-35))
    max1 = np.argmin(np.abs(freqs-220))
    freqR1 = freqs[min1:max1+1]
    # Range of 2nd Freq (240~680)
    min2 = np.argmin(np.abs(freqs-240))
    max2 = np.argmin(np.abs(freqs-680))
    freqR2 = freqs[min2:max2+1]
    
    
    # hanning window is applied before FFT
    win = get_window('hanning', fft_length)
    
    # define empty matrices 
    cnt = 0
    frequencies1 = np.zeros((datasize-startFFT)//FFTstep + 1)
    frequencies2 = np.zeros((datasize-startFFT)//FFTstep + 1)
    estimated_pin_location = np.zeros((frequencies1.size))
    estimated_time = np.zeros((frequencies1.size))
    
    
    # import C file
    
    so_file = "C:/Users/OGUNNIYI/Documents/research/LEMP C version/secular_solver_func_v3.so"
    
    my_func = CDLL(so_file)
    
    
    
    'testing point'
    
    def eig_soln(M,K):
        # Calculation of the natural frequencies using the generalized eigenvalue approach
        #also returns matrix of corresponding eigenvectors U1
        eigvals,eigvects = solve.eigh(K,M)
        eigvals=np.expand_dims(np.real(eigvals), axis=0)
        
        U1=eigvects
        wn_squared= eigvals
        
        FEA_wn= np.sqrt(wn_squared)
        Frequencies = FEA_wn/(2*np.pi) # Natural freq in Hz
        
        return(wn_squared,U1, Frequencies)
    
    
    #%%State 1 info
    ######################################################################################
        
    #np.set_printoptions(np.set_printoptions(precision=50))
    node_new = nodes
    # Generating the uniform grid    
    beam1 = np.linspace(0, 1, node_new , endpoint=True)
    pin_locations_FEA = beam1
        
    
    # Generating the uniform grid    
    # beam1 = np.linspace(0, 1, data["nodes"] , endpoint=True)
    # pin_locations_FEA = beam1
        
    beam_length = 0.3525     # length of the beam in meters
    beam_width = 0.0508   # width of the beam in meters
    beam_height = 0.00635 # thickness of the beam in meters
    beam_E = 209500000000 # Youngs modules of steel in Pa
    beam_density = 7900 # density of steel in kg/m^3
    accelerometer_mass = 0.07 # mass of accelerometer in kg
    beam_I = (beam_width*beam_height**3)/12 # caclulated moment of inertia
    beam_area = beam_width*beam_height
    pin_node_rotation_spring = 0 # set the value of the spring at the pinned connection  
        
    pin_locations_actual = pin_locations_FEA*beam_length
    beam_node_num = pin_locations_FEA.size
    beam_element = beam_node_num-1 # calculate the number of elements in the beam
    beam_el_length = beam_length/beam_element
    
    #%% State 1
    matrix_size = (beam_node_num)*2
    M = np.zeros((matrix_size,matrix_size))
    K = np.zeros((matrix_size,matrix_size))
        
    # for each element, add the element matrix into the global matirx
    for elem_num in range(0,beam_element):
        if elem_num == (beam_element-1):
            beam_density = beam_density + accelerometer_mass/(beam_area*beam_el_length)
          
           
        # define the mass matrix of a Euler-Bernoulli beam
        M_el = (beam_density*beam_area*beam_el_length)/420* \
        np.matrix([[156,22*beam_el_length,54,-13*beam_el_length], \
                        [22*beam_el_length,4*beam_el_length**2,13*beam_el_length,-3*beam_el_length**2], \
                        [54,13*beam_el_length,156,-22*beam_el_length], \
                        [-13*beam_el_length,-3*beam_el_length**2,-22*beam_el_length,4*beam_el_length**2]])
                
        # define the stiffness matrix of a Euler-Bernoulli beam
        K_el = (beam_E*beam_I)/beam_el_length**3* \
        np.matrix([[12,6*beam_el_length,-12,6*beam_el_length], \
                        [6*beam_el_length,4*beam_el_length**2,-6*beam_el_length,2*beam_el_length**2], \
                        [-12,-6*beam_el_length,12,-6*beam_el_length], \
                        [6*beam_el_length,2*beam_el_length**2,-6*beam_el_length,4*beam_el_length**2]])
                    
        n = (elem_num)*2
        M[n:n+4,n:n+4] = np.add(M[n:n+4,n:n+4],M_el)
        K[n:n+4,n:n+4] = np.add(K[n:n+4,n:n+4],K_el)    
                
    # for the fixed end on the left side, u_1 and u_2 = 0, so we can remove these columns
    # and rows form the matrixes. 
    # apply the boundary conditions
    
    K[0,0]=1e10
    K[1,1]=1e10
    
    K_ =K
    M_ = M
    #UNDEFORMED STATE- will always be the same 
    #so no need for function here
    #save as json file for easy reference
    number_modes=modes
    state1_measurements= eig_soln(M,K)
    wn_squared_= state1_measurements[0]
    old_frequencies= state1_measurements[2]
    U1_= state1_measurements[1]
    wn_squared_state1_=wn_squared_[:,0:number_modes]
    
    
    
    
    #%% Model Updating Simulation
    t0_3 = time.perf_counter() # time overall process
    total_2 =0
    for step in range(startFFT, datasize+1, FFTstep):
        t0_2 = time.perf_counter() # iteration time
        # Take FFT from the measured data
        accelerationF = acceleration[step-(fft_length):step]
        w_fft = np.abs(np.fft.rfft(win*accelerationF, fft_length))
        frequency1 = freqR1[np.argmax(w_fft[min1:max1+1])]
        frequencies1[cnt] = frequency1
        frequency2 = freqR2[np.argmax(w_fft[min2:max2+1])]
        frequencies2[cnt] = frequency2
    
        
        #  develop a set of unique system parameters based on the PDF centered around the last pin location
        inputs = np.zeros((particle_number),dtype=float)
        inputs_temp = np.zeros((1),dtype=float)
        ii=0
        
        # selecting particles(models) from the PDF function
        # roller upper limit = 0.9
        # roller lower limit = 0.1
        while ii < particle_number:
            inputs_temp =  np.random.randn(1)*step_pin_standard_dev + est_pin_loc        
            if inputs_temp <= 0.1:
                inputs_temp=0.1
            elif inputs_temp >= 0.9:
                inputs_temp=0.9
            if any((inputs[:]==inputs_temp)) == False:
                inputs[ii] = inputs_temp
                ii += 1
        
        # send the model parameter to the generalized eigenvalue problem
        frequencies_FEA = np.zeros((particle_number,freq_mode))
        for i in range(particle_number):
            t0_1 = time.perf_counter() #time taken to solve for each frequencies of new state
            frequencies_FEA[i,:] = GE.State2(pin_loc= inputs[i], number_modes=modes, number_nodes=nodes, M=M_, K = K_, matrix_size =matrix_size)
            # frequencies_FEA[i,:] = LEMP.LEMP(pin_loc= inputs[i], number_modes=modes, number_nodes=nodes, a= my_func, U1 = U1_, wn_squared_state1 = wn_squared_state1_)
            t1_1 = time.perf_counter()
            total_1 = (t1_1-t0_1)
        # t1 = time.perf_counter()
        # total = (t1-t0)   
        #%% Update the pin location (2 methods) 
        step_estimated_input = np.argmin(np.abs(frequencies_FEA[:,0]-frequency1))
        
        #%% Method 1: Nearest neighbor    
        # est_pin_loc = inputs[step_estimated_input]
        
        #%% Method 2: Linear Regression using Least Squares Fit
        E_freq = np.reshape(frequencies_FEA[:,0] - frequency1, (particle_number,1))
        Hmat = np.concatenate((np.ones((particle_number,1)), np.reshape(inputs,(particle_number,1))), axis = 1)
        a,b = np.linalg.lstsq(Hmat, E_freq, rcond=None)[0]
        est_pin_loc = -a/b
        xmin = np.min(inputs)
        xmax = np.max(inputs)    
        if est_pin_loc <= xmin:
            est_pin_loc = xmin
        elif est_pin_loc >= xmax:
            est_pin_loc = xmax
            
        t1_2 = time.perf_counter()
        t1_3 = time.perf_counter()
        total_2_ = (t1_2-t0_2) 
        total_2 = total_2 + total_2_
        #%% Adaptive Search Space (On or Off)
        step_pin_standard_dev = np.abs(frequencies_FEA[step_estimated_input,0] - frequency1)/frequency1 # On
        if step_pin_standard_dev <= 1e-12:
            step_pin_standard_dev = 1e-12
        # step_pin_standard_dev = 0.02 # Off, fixed search space
         
        #%% save and print
        estimated_pin_location[cnt] = np.copy(est_pin_loc)
        estimated_time[cnt] = np.copy(acceleration_time[step]) + delay_time
        # print('estimated beam location ' + str(np.round(estimated_pin_location[cnt],3)) + '%')
        
        cnt += 1
    
    # t1 = time.perf_counter()
    total_3 = (t1_3-t0_3)
    # print(total_1) # solver time
    # print(total_2/8623) # iteration time
    # print(total_3) # overall time
    time_data_particle_6[j-10] = total_2/8623
    solver_time[j-10] = total_1
    
    
    #%% Results and Plots
    
    estimated_pin_location_len = estimated_pin_location*beam_length
    
    estimated_LEMP = estimated_pin_location_len
    #measured
    measured_pin= np.load('error_data/measured_pin_loc.npy')
    measured_pin_location=measured_pin
    
    #time
    est_time= np.load('error_data/est_time.npy')
    measured_time= np.load('error_data/measured_time.npy')
    
    error_lemp= np.zeros(8623)

    # error_test=np.zeros(8623)
    
    measured=np.zeros(8623)
    
    time_error=np.zeros(18667)
    
    
    
    for i in range (8623):
        time_est= est_time[i]
        time_error= measured_time-time_est
        
        loc=np.argmin(np.abs(time_error))
        
        measured[i]= measured_pin[loc]
        
        error_lemp[i]=(estimated_LEMP[i]-measured[i])*1000
        error_lemp1 = abs(error_lemp)
        
    
    av_error_lemp_new=np.average(error_lemp1[~np.isnan(error_lemp1)])
    # print(av_error_lemp_new)
    MAE_particle_6[j-10] = av_error_lemp_new
    av_measured=np.average(measured[~np.isnan(measured)])*1000
    
    SNR_lemp_new = 10* math.log10((av_measured/av_error_lemp_new))
    # print(SNR_lemp_new)
    SNR_particle_6[j-10] =SNR_lemp_new
    print(j)
  
    TRAC1 = np.matmul(measured, estimated_LEMP)**2 / (np.matmul(measured.T, measured)*np.matmul(estimated_LEMP.T, estimated_LEMP))
    Trac[j-10] = TRAC1
    
time_data_particle_6 = list(time_data_particle_6)
MAE_particle_6 = list(MAE_particle_6)
SNR_particle_6=list(SNR_particle_6)
solver_time = list(solver_time)
Trac = list(Trac)


data_particle_6 = {'time_particle_3':time_data_particle_6,'MAE_particle_3':MAE_particle_6, 'SNR_particle_3':SNR_particle_6, 'solver_time_GE':solver_time, 'Trac':Trac}   

with open('particle_data/500_GE_particle_3_with_solver2.json', 'w') as outfile:
    json.dump(data_particle_6, outfile)
            
    
    # plt.figure(figsize=(7,3))
    # plt.plot(estimated_time, estimated_pin_location_len*1000,'--',label='estimated', linewidth='1')
    # plt.plot(measured_time, measured_pin_location*1000,'-',label='measured', linewidth='1')
    # #plt.title( str(particle_number) +' samples; and '+ str(step_pin_standard_dev*100) + '% standard deviation', fontsize=12) 
    # # plt.title( str(particle_number) +' samples; and ' + 'adaptive standard deviation', fontsize=12) 
    # plt.grid('on')
    # plt.legend(framealpha=1,loc='best', fontsize=12)
    # plt.xlabel('time (s)', fontsize=12)
    # plt.ylabel('roller location (mm)', fontsize=12)
    # plt.tight_layout()
    # plt.show()
    # # plt.savefig('LEMP_26_BR_qual.pdf', dpi=200)
    
    # estimated_pin_location_len = list(estimated_pin_location_len)
    # estimated_time = list(estimated_time)
    # measured_pin_location=list(measured_pin_location)
    # measured_time=list(measured_time)

# save data
# Data = {'EstLoc':estimated_pin_location_len, 'EstTime':estimated_time, 
#         'MeasLoc':measured_pin_location, 'MeasTime':measured_time, 
#         'Samples':particle_number, 'STD':step_pin_standard_dev}



# LEMP_21_BR_P3 = json.dumps(Data)

# with open('Time plots/error_data/data.json', 'w') as outfile:
#     json.dump(Data, outfile)

# np.save('LEMP_21_BR_P3',Data)
#sio.savemat('ExpData_BR_'+str(particle_number)+ 'samp_'+ str(int(step_pin_standard_dev*100)) + 'std.mat', Data)
#sio.savemat('ExpData_BR_'+str(particle_number)+ '_samp_'+ 'adaptive_std.mat', Data)

