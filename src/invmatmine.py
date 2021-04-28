# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 02:17:36 2021

@author: AndresMoreno
"""
import scipy.linalg.lapack as lap
import numpy as np 

def invertmat(S):
    
    n=np.size(S[:,0])
    
    thresh=0.00001
    A=S
    
    u,s,vt,info=lap.zgesvd(A) 
    
    sigma=s
    inv_S=np.zeros_like(S, dtype=np.complex128)
    inv_sigma=np.zeros_like(sigma)
    
    for i in range(n): 
        inv_sigma[i]=sigma[i]/(sigma[i]**2+thresh**2)
    
    for i in range(n):
        
        inv_S[:,i]=np.conj(u[i,:])*inv_sigma
        
    inv_S=np.matmul(np.conj(np.transpose(vt)),inv_S)
    
    
    return inv_S
          