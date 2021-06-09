# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 12:01:32 2021

@author: AndresMoreno
"""

import numpy as np
from src.abinitio import ab_par
import src.mldreader as mld


def explicitconfigs(T):
    ab = ab_par()

    closed = ab.closed_orb
    active = ab.act_orb
    nel = int(ab.n_el)

    count = 0
    count2 = 0
    calpha = 0
    cbeta = 0
    configs = T.getconfigs()

    configsmod = []
    #print(T.getconfigs())
    alpha = np.zeros((np.size(T.configs), int(nel / 2)))
    beta = np.zeros_like(alpha)

    for i in configs:
        configsmod.append('2' * closed + i)
    for i in configsmod:
        count2 = 0
        calpha = 0
        cbeta = 0
        for j in i:
            if j == '2':

                alpha[count, calpha] = count2
                beta[count, cbeta] = count2
                count2 += 1
                calpha += 1
                cbeta += 1

            elif j == 'a':

                alpha[count, calpha] = count2
                count2 += 1
                calpha += 1
            elif j == 'b':

                beta[count, cbeta] = count2
                count2 += 1
                cbeta += 1
            else:

                count2 += 1

        count += 1

    return alpha, beta


def dets(T1, T2, MOovs):
    alpha_1, beta_1 = explicitconfigs(T1)
    alpha_1 = alpha_1.astype(int)
    beta_1 = beta_1.astype(int)

    alpha_2, beta_2 = explicitconfigs(T2)
    alpha_2 = alpha_2.astype(int)
    beta_2 = beta_2.astype(int)

    confign_1 = np.size(alpha_1[:, 0])
    confign_2 = np.size(alpha_2[:, 0])

    detmat_alpha = np.zeros((confign_1, confign_2))
    detmat_beta = np.zeros((confign_1, confign_2))

    detmat_alpha2 = np.zeros_like(detmat_alpha)
    detmat_beta2 = np.zeros_like(detmat_beta)

    alpha_unique_1, u_a_1 = np.unique(alpha_1, return_index=True, axis=0)
    beta_unique_1, u_b_1 = np.unique(beta_1, return_index=True, axis=0)

    alpha_unique_1 = alpha_1[np.sort(u_a_1), :].astype(int)
    beta_unique_1 = beta_1[np.sort(u_b_1), :].astype(int)

    alpha_unique_2, u_a_2 = np.unique(alpha_2, return_index=True, axis=0)
    beta_unique_2, u_b_2 = np.unique(beta_2, return_index=True, axis=0)

    alpha_unique_2 = alpha_2[np.sort(u_a_2), :].astype(int)
    beta_unique_2 = beta_2[np.sort(u_b_2), :].astype(int)

    lists_alpha_1 = np.ones(np.size(alpha_1[:, 0])).astype(int)
    lists_beta_1 = np.ones(np.size(beta_1[:, 0])).astype(int)

    lists_alpha_2 = np.ones(np.size(alpha_2[:, 0])).astype(int)
    lists_beta_2 = np.ones(np.size(beta_2[:, 0])).astype(int)

    norbs_alpha_1 = np.size(alpha_unique_1[0, :])
    norbs_beta_1 = np.size(beta_unique_1[0, :])

    norbs_alpha_2 = np.size(alpha_unique_2[0, :])
    norbs_beta_2 = np.size(beta_unique_2[0, :])

    nconfs_u_alpha_1 = np.size(alpha_unique_1[:, 0])
    nconfs_u_beta_1 = np.size(beta_unique_1[:, 0])

    nconfs_u_alpha_2 = np.size(alpha_unique_2[:, 0])
    nconfs_u_beta_2 = np.size(beta_unique_2[:, 0])

    dets_red_alpha = np.zeros((nconfs_u_alpha_1, nconfs_u_alpha_2))
    dets_red_beta = np.zeros((nconfs_u_beta_1, nconfs_u_beta_2))

    for i in range(confign_1):
        for j in range(nconfs_u_alpha_1):
            if (alpha_1[i, :] == alpha_unique_1[j, :]).all():
                lists_alpha_1[i] = j

    for i in range(confign_2):
        for j in range(nconfs_u_alpha_2):
            if (alpha_2[i, :] == alpha_unique_2[j, :]).all():
                lists_alpha_2[i] = j

    for i in range(confign_1):
        for j in range(nconfs_u_beta_1):
            if (beta_1[i, :] == beta_unique_1[j, :]).all():
                lists_beta_1[i] = j

    for i in range(confign_2):
        for j in range(nconfs_u_beta_2):
            if (beta_2[i, :] == beta_unique_2[j, :]).all():
                lists_beta_2[i] = j

    for i in range(nconfs_u_alpha_1):
        for j in range(nconfs_u_alpha_2):
            mat = np.zeros((norbs_alpha_1, norbs_alpha_2))
            for k in range(norbs_alpha_1):
                for l in range(norbs_alpha_2):
                    ind1 = alpha_unique_1[i, k]
                    ind2 = alpha_unique_2[j, l]
                    mat[k, l] = MOovs[ind1, ind2]
          #  print('matrix: ', i, j, mat)
            dets_red_alpha[i, j] = np.linalg.det(mat)
          #  print('determinant_alpha: ', i, j, dets_red_alpha[i, j])

    for i in range(nconfs_u_beta_1):
        for j in range(nconfs_u_beta_2):
            mat = np.zeros((norbs_beta_1, norbs_beta_2))
            for k in range(norbs_beta_1):
                for l in range(norbs_beta_2):
                    ind1 = alpha_unique_1[i, k]
                    ind2 = alpha_unique_2[j, l]
                    mat[k, l] = MOovs[ind1, ind2]
            dets_red_beta[i, j] = np.linalg.det(mat)

    sign = np.ones_like(detmat_alpha)

    for i in range(np.size(lists_alpha_1)):
        for j in range(np.size(lists_alpha_2)):
            if T1.configs[i].replace('2', '0') != T2.configs[j].replace('2', '0'):
                sign[i, j] = -1.00

    # print('lists_alpha_1', alpha_1)
    # print('lists_alpha_2: ', alpha_unique_1)
    # print('lists_alpha_2: ', lists_alpha_1)
    # print('lists_alpha_2', alpha_2)
    # print('lists_alpha_2: ', alpha_unique_2)
    # print('lists_2: ', lists_alpha_2)

    for i in range(np.size(lists_alpha_1)):
        for j in range(np.size(lists_alpha_2)):
            detmat_alpha[i, j] = sign[i, j] * dets_red_alpha[lists_alpha_1[i], lists_alpha_2[j]]
            # print(i, j, lists_alpha_1[i], lists_alpha_2[j], detmat_alpha[i, j])
    for i in range(np.size(lists_beta_1)):
        for j in range(np.size(lists_beta_2)):
            detmat_beta[i, j] = dets_red_beta[lists_beta_1[i], lists_beta_2[j]]

    # for i in range(np.size(lists_alpha_1)):
    #     for j in range(np.size(lists_alpha_2)):
    #         mat = np.zeros((norbs_alpha_1, norbs_alpha_2))
    #         for k in range(norbs_alpha_1):
    #             for l in range(norbs_alpha_2):
    #                 ind1 = alpha_1[i, k]
    #                 ind2 = alpha_2[j, l]
    #                 mat[k, l] = MOovs[ind1, ind2]
    #         detmat_alpha2[i, j] = np.linalg.det(mat)
    #
    # for i in range(np.size(lists_beta_1)):
    #     for j in range(np.size(lists_beta_2)):
    #         mat = np.zeros((norbs_beta_1, norbs_beta_2))
    #         for k in range(norbs_beta_1):
    #             for l in range(norbs_beta_2):
    #                 ind1 = beta_1[i, k]
    #                 ind2 = beta_2[j, l]
    #                 mat[k, l] = MOovs[ind1, ind2]
    #         detmat_beta2[i, j] = np.linalg.det(mat)

    return detmat_alpha, detmat_beta


def overlap(T1, T2):
    megamatrix1, M1,fff = mld.readingmolden(T1.getfilecalc())
    megamatrix2, M2,fff = mld.readingmolden(T2.getfilecalc())

    MOverlap = mld.MOverlapcalc(megamatrix1, megamatrix2, M1, M2)

    # with np.printoptions(threshold=np.inf):
    #     print('movs :', MOverlap)
    CIS1 = T1.getcivecs()
    CIS2 = T2.getcivecs()

    # CIS2[2:, :] = -CIS2[2:, :]

    detmat_alpha, detmat_beta = dets(T1, T2, MOverlap)

  #  print(detmat_alpha)
    # Ground state overlap calculation
  #  print(np.sum(CIS1[:, 0] * CIS2[:, 0]))
    OverlapG = 0.0
  #  print(np.size(detmat_alpha))
    for i in range(np.size(T1.getconfigs())):
        for j in range(np.size(T2.getconfigs())):
            OverlapG += CIS1[i, 0] * CIS2[j, 0] * detmat_alpha[i, j] * detmat_beta[i, j]
            # print('steps: ', i, j)
            # print('coeffs : ', CIS1[i, 0] * CIS2[j, 0])
            # print('det_alpha ', detmat_alpha[i, j])
            # print('det_beta ', detmat_beta[i, j])

    OverlapExc = 0.0
    for i in range(np.size(T1.getconfigs())):
        for j in range(np.size(T2.getconfigs())):
            OverlapExc += CIS1[i, 1] * CIS2[j, 1] * detmat_alpha[i, j] * detmat_beta[i, j]

    return OverlapG, OverlapExc
