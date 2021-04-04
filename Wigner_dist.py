import re
import numpy as np
from input import initgeom as ig
import Constants
import random


def read_freq(file='freq.out', ndf=18, linear=False):
    if linear:
        nc = ndf - 5
        numcut = 5
    else:
        nc = ndf - 6
        numcut = 6

    var = 1
    str2search = 'Mass Weighted 2nd Derivative Matrix Eigenvalues'
    str2search2 = "Mass Weighted 2nd Derivative Matrix Eigenvectors"
    str2search3 = 'Low Vibration'
    eigenval = []
    eigenvec = []
    ndf_var = 0
    ndf_var2 = 0
    ndf_var3 = 0
    with open(file, 'r') as f:
        for line in f:
            if str2search in line:
                var = 2
                continue
            if var == 2:
                n = re.findall("\d+\.\d+", line)
                if n:
                    eigenval.extend(n)

            if str2search2 in line:
                var = 3
                continue
            if var == 3:
                n = re.findall("\d+\.\d+", line)

                if n:

                    if ndf_var <= ndf - 1:
                        eigenvec.append(n)
                        ndf_var += 1
                    elif ndf_var2 <= ndf - 1:
                        eigenvec[ndf_var2].extend(n)
                        ndf_var2 += 1
                    else:
                        eigenvec[ndf_var3].extend(n)
                        ndf_var3 += 1

            if str2search3 in line:
                var = 1
                continue

    eigenval = np.asarray(eigenval)
    eigenvec = np.asarray(eigenvec)

    return eigenval[:], eigenvec[:, :]


def WignerSampling():
    ph = Constants.physconst()

    geo = ig()

    eigval, eigvec = read_freq('freq.out', geo.ndf, False)

    eigval = eigval.astype(np.float)
    eigvec = eigvec.astype(np.float)

    nc = geo.ndf - 6
    inmod = np.asarray(range(6, geo.ndf, 1))

    alpha = np.zeros(nc)
    qmode = np.zeros(nc)
    pmode = np.zeros(nc)
    r_c = np.zeros((3, geo.natoms))
    p_c = np.zeros((3, geo.natoms))

    for i in range(nc):

        index = inmod[i]

        omega = np.sqrt(float(eigval[index]) / ph.amu)

        alpha[i] = 0.5 * omega

        rq = RandomGaussianSample()
        rp = RandomGaussianSample()

        qmode[i] = 0.5 * rq * (np.sqrt(alpha[i])) ** (-1)
        pmode[i] = np.sqrt(alpha[i]) * rp

        ic = 0

        for n in range(geo.natoms):
            for j in range(3):
                r_c[j, n] = r_c[j, n] + qmode[i] * eigvec[ic, index]
                p_c[j, n] = p_c[j, n] + pmode[i] * eigvec[ic, index]
                ic += 1

    p = np.zeros(geo.ndf)
    q = np.zeros(geo.ndf)
    idf = 0
    for i in range(geo.natoms):
        for j in range(3):
            q[idf] = r_c[j, i]
            p[idf] = p_c[j, i]
            idf += 1

    return q, p


class vars:
    iset = 0
    gset = 0


def RandomGaussianSample():
    v = vars()
    if v.iset != 0:
        vars.iset = 0
        gasdev = v.gset

    else:
        while True:
            x = random.uniform(0.0, 1.0)
            v1 = 2.0 * x - 1
            y = random.uniform(0.0, 1.0)
            v2 = 2.0 * y - 1
            rsq = v1 ** 2 + v2 ** 2
            if rsq < 1.0 and rsq != 0.0:
                break
        fac = np.sqrt(-2.0 * np.log(rsq) / rsq)
        vars.gset = v1 * fac
        gasdev = v2 * fac
        vars.iset = 1

    return gasdev / np.sqrt(2.0)


'''Routine to mass weight and normalize the P obtained in every step, including the Wigner sampling'''


def WPrep(p_init):
    ph = Constants.physconst()
    y = p_init
    geo = ig()
    idf = 0
    p = np.zeros((3, geo.natoms))
    for i in range(geo.natoms):
        for j in range(3):
            p[j, i] = y[idf] * np.sqrt(geo.masses[i])
            idf += 1

    totmass = np.sum(geo.masses)

    partial_mass = geo.masses / totmass

    totmom_i = np.sum(p * partial_mass, axis=1)

    p_corr = np.zeros(np.shape(p))
    idf = 0
    y_fin = np.zeros(geo.ndf)

    for i in range(geo.natoms):
        for j in range(3):
            p_corr[j, i] = p[j, i] - totmom_i[j]
            y_fin[idf] = p_corr[j, i] / np.sqrt(geo.masses[i])
            idf += 1

    print(y_fin)
