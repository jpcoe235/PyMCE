# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 01:09:42 2021

@author: AndresMoreno
"""
import numpy as np


def matexpAI(B):
    Id = np.zeros_like(B)

    for i in range(len(Id[:,0])):
        Id[i, i] = 1.00

    expB_prev = Id
    expB=0.0

    for i in range(4, 21):


        Bn = B / (2.00 ** i)

        Bn2 = np.matmul(Bn, Bn)

        Bn3 = np.matmul(Bn2, Bn)

        Bn4 = np.matmul(Bn2, Bn2)

        tayl = Id + Bn + Bn2 / 2.0 + Bn3 / 6.0 + Bn4 / 24.0

        tayl_prod = tayl

        for j in range(1, i + 1):
            tayl_prod = np.matmul(tayl_prod, tayl_prod)

        expB = tayl_prod

        if np.max(abs(expB_prev - expB)) < 1e-7:
         
            return expB
        else:
            expB_prev = expB

    return expB
