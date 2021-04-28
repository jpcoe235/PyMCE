# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 02:00:07 2021

@author: AndresMoreno
"""
import numpy as np
import time 
from matexp import matexpAI
from scipy.linalg import expm
from invmatmine import invertmat
from Sinveigs import Sinv2
t1=np.zeros(1000)
t2=np.zeros(1000)
t3=np.zeros(1000)

for i in range(1000):
    B=np.random.rand(5,5)+np.random.rand(5,5)*1j

    time1=time.perf_counter()
    np.linalg.inv(B)
    time2=time.perf_counter()

    time3=time.perf_counter()

    invertmat(B)

    time4= time.perf_counter()

    time5 = time.perf_counter()

    Sinv2(B)

    time6 = time.perf_counter()

    t1[i]=time2-time1
    t2[i]=time4-time3
    t3[i] = time6 - time5

    
    
print('time my function: ', np.mean(t1))

print('time built-in: ', np.mean(t2))

print('time my function: ', np.mean(t3))

