# -*- coding: utf-8 -*-
"""
Created on Sat Jun  5 23:17:22 2021

@author: AndresMoreno
"""

import sys
from sys import argv
import numpy as np
from scipy.special import factorial2 as factd
import math


# print "Specify the mld file u want to use:"
# file = raw_input("> ")


def readingmolden(file):
    with open(file, "r") as txt:
        line = txt.read()

    # Looking for Atoms in molden file
    first1, last1 = line.index("[Atoms]") + 7, line.index("[GTO]")
    first2, last2 = line.index("[GTO]") + 5, line.index("[MO]")
    Atoms = line[first1:last1]
    Atoms = Atoms.split('\n')[:-1]
    matoms = []
    # Defining x,y,z variables for each atom
    del Atoms[0]
    for lines in Atoms:
        coord = lines.split()
        matoms.append([int(coord[1]), int(coord[2]), float(coord[3]) / 0.52917721092, float(coord[4]) / 0.52917721092,
                       float(coord[5]) / 0.52917721092])

    # print x[1],y[1],z[1]

    BASES = []
    GTOs = line[first2:last2]
    GTOs = GTOs.split('\n')[:-1]
    OS = []
    n = 1
    listatoms = []
    listindex = []

    # for i in range(0,len(GTOs)-1):
    for n in range(1, len(Atoms) + 1):
        if n < 10:
            rata = '   ' + str(n) + ' ' + str(0)
            #        if GTOs[i]==rata:
            i = GTOs.index(rata)
            GTOs[i] = 'Atom' + str(n) + ' ' + str(0)
            listatoms.append(n)
            listindex.append(i)



        else:
            rata = '  ' + str(n) + ' ' + str(0)
            i = GTOs.index(rata)
            GTOs[i] = 'Atom' + str(n) + ' ' + str(0)
            listatoms.append(n)
            listindex.append(i)

    finalorder = [x for _, x in sorted(zip(listindex, listatoms))]

    GTOs = ' '.join(GTOs)
    GTOs = GTOs.split(' ')

    GTOs = [i for i in GTOs if i]

    GTOs = list(map(lambda s: s.strip(), GTOs))
    # GTOs.remove('0')
    n = 0
    nn = 0
    Counter = 0

    for k in finalorder:
        An = finalorder[Counter]
        if k == len(Atoms):
            busca1 = 'Atom' + str(len(Atoms))
            first3 = GTOs.index(busca1) + 1
            last3 = len(GTOs)
        else:
            final = finalorder[finalorder.index(k) + 1]
            busca1 = 'Atom' + str(k)
            busca2 = 'Atom' + str(final)
            first3, last3 = GTOs.index(busca1) + 1, GTOs.index(busca2)

        #    OS.append(GTOs[first3-1])

        GTOS1 = GTOs[first3:last3]
        for i in range(1, len(GTOs[first3:last3])):

            if GTOS1[i - 1] == 's':
                n = n + 1
                if (GTOS1[i] == '1'):
                    OS.append(
                        [n, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 0, 0,
                         0])

                elif (GTOS1 != 1):
                    for j in range(2, int(GTOS1[i]) * 2 + 2, 2):
                        OS.append(
                            [n, An, float(GTOS1[i + j].replace('D', 'E')), float(GTOS1[i + j + 1].replace('D', 'E')), 0,
                             0, 0])

            elif GTOS1[i - 1] == 'p':
                n = n + 1
                if GTOS1[i] == '1':

                    OS.append(
                        [n, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 1, 0,
                         0])
                    OS.append(
                        [n + 1, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 0,
                         1, 0])
                    OS.append(
                        [n + 2, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 0,
                         0, 1])

                else:
                    for j in range(2, int(GTOS1[i]) * 2 + 2, 2):
                        OS.append(
                            [n, An, float(GTOS1[i + j].replace('D', 'E')), float(GTOS1[i + j + 1].replace('D', 'E')), 1,
                             0, 0])
                        OS.append([n + 1, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 0, 1, 0])
                        OS.append([n + 2, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 0, 0, 1])
                n = n + 2
            elif GTOS1[i - 1] == 'd':
                n = n + 1
                if GTOS1[i] == '1':

                    OS.append(
                        [n, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 2, 0,
                         0])
                    OS.append(
                        [n + 1, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 0,
                         2, 0])
                    OS.append(
                        [n + 2, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 0,
                         0, 2])
                    OS.append(
                        [n + 3, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 1,
                         1, 0])
                    OS.append(
                        [n + 4, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 1,
                         0, 1])
                    OS.append(
                        [n + 5, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 0,
                         1, 1])

                else:
                    for j in range(2, int(GTOS1[i]) * 2 + 2, 2):
                        OS.append(
                            [n, An, float(GTOS1[i + j].replace('D', 'E')), float(GTOS1[i + j + 1].replace('D', 'E')), 2,
                             0, 0])
                        OS.append([n + 1, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 0, 2, 0])
                        OS.append([n + 2, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 0, 0, 2])
                        OS.append([n + 3, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 1, 1, 0])
                        OS.append([n + 4, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 1, 0, 1])
                        OS.append([n + 5, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 0, 1, 1])
                n = n + 5

            elif GTOS1[i - 1] == 'f':
                n = n + 1
                if GTOS1[i] == '1':

                    OS.append(
                        [n, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 3, 0,
                         0])
                    OS.append(
                        [n + 1, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 0,
                         3, 0])
                    OS.append(
                        [n + 2, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 0,
                         0, 3])
                    OS.append(
                        [n + 3, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 1,
                         2, 0])
                    OS.append(
                        [n + 4, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 2,
                         1, 0])
                    OS.append(
                        [n + 5, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 2,
                         0, 1])
                    OS.append(
                        [n + 6, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 1,
                         0, 2])
                    OS.append(
                        [n + 7, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 0,
                         1, 2])
                    OS.append(
                        [n + 8, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 0,
                         2, 1])
                    OS.append(
                        [n + 9, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 1,
                         1, 1])


                else:
                    for j in range(2, int(GTOS1[i]) * 2 + 2, 2):
                        OS.append(
                            [n, An, float(GTOS1[i + j].replace('D', 'E')), float(GTOS1[i + j + 1].replace('D', 'E')), 3,
                             0, 0])
                        OS.append([n + 1, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 0, 3, 0])
                        OS.append([n + 2, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 0, 0, 3])
                        OS.append([n + 3, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 1, 2, 0])
                        OS.append([n + 4, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 2, 1, 0])
                        OS.append([n + 5, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 2, 0, 1])
                        OS.append([n + 6, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 1, 0, 2])
                        OS.append([n + 7, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 0, 1, 2])
                        OS.append([n + 8, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 0, 2, 1])
                        OS.append([n + 9, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 1, 1, 1])

                n = n + 9

            elif GTOS1[i - 1] == 'g':
                n = n + 1
                if GTOS1[i] == '1':

                    OS.append(
                        [n, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 4, 0,
                         0])
                    OS.append(
                        [n + 1, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 0,
                         4, 0])
                    OS.append(
                        [n + 2, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 0,
                         0, 4])
                    OS.append(
                        [n + 3, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 3,
                         1, 0])
                    OS.append(
                        [n + 4, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 3,
                         0, 1])
                    OS.append(
                        [n + 5, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 1,
                         3, 0])
                    OS.append(
                        [n + 6, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 0,
                         3, 1])
                    OS.append(
                        [n + 7, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 1,
                         0, 3])
                    OS.append(
                        [n + 8, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 0,
                         1, 3])
                    OS.append(
                        [n + 9, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')), 2,
                         2, 0])
                    OS.append(
                        [n + 10, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')),
                         2, 0, 2])
                    OS.append(
                        [n + 11, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')),
                         0, 2, 2])
                    OS.append(
                        [n + 12, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')),
                         2, 1, 1])
                    OS.append(
                        [n + 13, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')),
                         1, 2, 1])
                    OS.append(
                        [n + 14, An, float(GTOS1[i + 2].replace('D', 'E')), float(GTOS1[i + 2 + 1].replace('D', 'E')),
                         1, 1, 2])

                else:
                    for j in range(2, int(GTOS1[i]) * 2 + 2, 2):
                        OS.append(
                            [n, An, float(GTOS1[i + j].replace('D', 'E')), float(GTOS1[i + j + 1].replace('D', 'E')), 4,
                             0, 0])
                        OS.append([n + 1, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 0, 4, 0])
                        OS.append([n + 2, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 0, 0, 4])
                        OS.append([n + 3, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 3, 1, 0])
                        OS.append([n + 4, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 3, 0, 1])
                        OS.append([n + 5, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 1, 3, 0])
                        OS.append([n + 6, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 0, 3, 1])
                        OS.append([n + 7, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 1, 0, 3])
                        OS.append([n + 8, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 0, 1, 3])
                        OS.append([n + 9, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 2, 2, 0])
                        OS.append([n + 10, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 2, 0, 2])
                        OS.append([n + 11, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 0, 2, 2])
                        OS.append([n + 12, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 2, 1, 1])
                        OS.append([n + 13, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 1, 2, 1])
                        OS.append([n + 14, An, float(GTOS1[i + j].replace('D', 'E')),
                                   float(GTOS1[i + j + 1].replace('D', 'E')), 1, 1, 2])

                n = n + 14

        Counter = Counter + 1

    # Begin the best part of the reading
    #
    dim = n + 3
    first4, last4 = line.index("[MO]") + 4, len(line)
    OS = sorted(OS, key=lambda OS: OS[0])

    MO = line[first4:last4]
    MO = MO.split('\n')[:-1]
    init = 0
    end = dim
    del MO[0]

    symm = []
    energy = []
    occ = []
    spin = []
    MM = []
    jotair = 1
    n = 1
    nmo = 1
    mos = []
    while end < len(MO) + 2:

        for i in MO[init:end + 1]:

            data = i.split(' ')
            data = list(filter(None, data))

            if data[0] == 'Sym=':
                symm.append(float(data[1]))

            elif data[0] == 'Ene=':
                energy.append(float(data[1]))

            elif data[0] == 'Occup=':
                occ.append(float(data[1]))
                r = float(data[1])

            elif data[0] == 'Spin=':
                spin.append(data[1])

            elif r != 0.0:

                MM.append([nmo, int(data[0]), r, float(data[1])])

        init = end + 1
        end = dim + init
        if r > 0.0:
            nmo = nmo + 1
            n = n + 1

    megamatrix = []

    for j in OS:

        for i in matoms:

            if j[1] == i[0]:
                megamatrix.append([i[0], i[2], i[3], i[4], j[0], j[2], j[3], j[4], j[5], j[6]])

    megamegamatrix = []

    for j in MM:

        for i in megamatrix:

            if i[4] == j[1]:
                megamegamatrix.append([i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9], j[2], j[3], j[0]])

    MM = np.asarray(MM)

    megamatrix = np.asarray(megamatrix)
    numGa = np.asarray(megamatrix[:, 4]).astype(int)
    MOs = np.zeros((nmo - 1, np.max(numGa)))
    j = 0
    i = 0

    for lines in MM:
        MOs[i, j] = lines[3]
        j = j + 1
        if j >= np.max(numGa):
            j = 0
            i += 1
    symindex=np.argsort(symm).astype(int)
    print(symindex[0:nmo-1])
    MOs=MOs[symindex[0:nmo-1],:]
    return megamatrix, MOs,symindex[0:nmo-1]


def MOverlapcalc(megamatrix1, megamatrix2, MOs1, MOs2):
    numGa_1 = np.asarray(megamatrix1[:, 4]).astype(int)
    ga_1 = np.asarray(megamatrix1[:, 5])
    ci_1 = np.asarray(megamatrix1[:, 6])
    x_1 = np.asarray(megamatrix1[:, 1])
    y_1 = np.asarray(megamatrix1[:, 2])
    z_1 = np.asarray(megamatrix1[:, 3])
    l_1 = np.asarray(megamatrix1[:, 7]).astype(int)
    m_1 = np.asarray(megamatrix1[:, 8]).astype(int)
    n_1 = np.asarray(megamatrix1[:, 9]).astype(int)

    numGa_2 = np.asarray(megamatrix2[:, 4]).astype(int)
    ga_2 = np.asarray(megamatrix2[:, 5])
    ci_2 = np.asarray(megamatrix2[:, 6])
    x_2 = np.asarray(megamatrix2[:, 1])
    y_2 = np.asarray(megamatrix2[:, 2])
    z_2 = np.asarray(megamatrix2[:, 3])
    l_2 = np.asarray(megamatrix2[:, 7]).astype(int)
    m_2 = np.asarray(megamatrix2[:, 8]).astype(int)
    n_2 = np.asarray(megamatrix2[:, 9]).astype(int)

    normA = (2.00 / np.pi) ** 0.75;
    AOoverlap = np.zeros((np.max(numGa_1), np.max(numGa_2)))
    for n1 in range(max(numGa_1)):
        for n2 in range(max(numGa_2)):

            index1 = np.where(numGa_1 == n1 + 1)[0]
            index2 = np.where(numGa_1 == n2 + 1)[0]
            overlap = 0.0

            for i in index1:
                normB = 2.00 ** (l_1[i] + m_1[i] + n_1[i])
                normC = ga_1[i] ** ((2 * l_1[i] + 2 * m_1[i] + 2 * n_1[i] + 3) / 4.0)
                normD = doublefactorial(2 * l_1[i] - 1) * doublefactorial(2 * m_1[i] - 1) * doublefactorial(
                    2 * n_1[i] - 1)
                nrmi = normA * normB * normC / normD ** 0.5

                for j in index2:
                    normB = 2.00 ** (l_2[j] + m_2[j] + n_2[j])
                    normC = ga_2[j] ** ((2 * l_2[j] + 2 * m_2[j] + 2 * n_2[j] + 3) / 4.0)
                    normD = doublefactorial(2 * l_2[j] - 1) * doublefactorial(2 * m_2[j] - 1) * doublefactorial(
                        2 * n_2[j] - 1)
                    nrmj = normA * normB * normC / normD ** 0.5
                    intx = gaussintang(x_1[i], x_2[j], ga_1[i], ga_2[j], l_1[i], l_2[j])
                    inty = gaussintang(y_1[i], y_2[j], ga_1[i], ga_2[j], m_1[i], m_2[j])
                    intz = gaussintang(z_1[i], z_2[j], ga_1[i], ga_2[j], n_1[i], n_2[j])

                    overlap += ci_1[i] * ci_2[j] * nrmi * nrmj * intx * inty * intz

            AOoverlap[n1, n2] = overlap
        # AOoverlap[n2, n1] = overlap
   # print(AOoverlap)
    MOoverlap = np.matmul(np.matmul(MOs1, AOoverlap), np.transpose(MOs2))

    MOverlap2 = np.zeros((np.size(MOs1[:, 0]), np.size(MOs2[:, 0])))

    for m1 in range(np.size(MOs1[:, 0])):
        for m2 in range(np.size(MOs2[:, 0])):
            for i in range(np.size(AOoverlap[:, 0])):
                for j in range(np.size(AOoverlap[0, :])):
                    MOverlap2[m1, m2] += MOs1[m1, i] * MOs2[m2, j] * AOoverlap[i, j]

   # print('second overlap :', MOverlap2)
 #   MOoverlap = MOverlap2
    return MOoverlap


def gaussianint(A, B, a, b):
    p = (a + b)
    mu = (a * b) / p
    xab = A - B
    prexp = np.sqrt(np.pi / p)
    integral = prexp * np.exp(-mu * xab ** 2)
    return integral, xab, p


def gaussintang(A, B, a, b, l1, l2):
    S = np.zeros((20, 20))

    S[0, 0], xab, p = gaussianint(A, B, a, b)
    p = (a + b)
    px = (a * A + b * B) / p
    xpa = -(b / p) * xab
    xpb = (a / p) * xab

    # print(xpa*S[0,0])
    for i in range(1, l1 + 2):
        for j in range(1, l2 + 2):
            S[i, j - 1] = xpa * S[i - 1, j - 1] + 1 / (2 * p) * ((i - 1) * S[i - 2, j - 1] + (j - 1) * S[i - 1, j - 2])
            S[i - 1, j] = xpb * S[i - 1, j - 1] + 1 / (2 * p) * ((i - 1) * S[i - 2, j - 1] + (j - 1) * S[i - 1, j - 2])

    return S[l1, l2]


def doublefactorial(n):
    if n <= 0:
        return 1
    else:
        return factd(n)

        # dimmensions

# Printing important parameters


# data
# f=open('OS.dat','w')
# r=str(len(OS))+'\n'
# f.write(r)
# np.savetxt(f,OS,fmt='%i %i %.18e %.18e %i %i %i')
# f.close()
#
# f=open('MM.dat','w')
# f.write(str(len(MM))+'\n')
# np.savetxt(f,MM,fmt='%i %d %.18e %.18e')
# f.close()
#
# f=open('Atoms.dat','w')
# f.write(str(len(matoms))+'\n')
# np.savetxt(f,matoms,fmt='%i %i %.18f %.18f %.18f')
# f.close()
#
# f=open('Prop.dat','w')
# f.write(str(len(symm))+'\n')
# np.savetxt(f,propMO,fmt='%-.1f %.6e %.18e')
# f.close()
# f=open('TotalMatrix.dat','w')
# f.write(str(len(megamegamatrix))+'\n')
# f.write(str(n-1)+'\n')
# np.savetxt(f,np.asmatrix(megamegamatrix),fmt=' %i %.18f %.18f %.18f %i %.18f %.18f %i %i %i %.18f %.18f %i')
# f.close()
##
#
