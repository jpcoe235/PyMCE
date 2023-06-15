import re
import numpy as np
from input import initgeom as ig
import input
import Constants
import random


def read_freq(file='vibs.molden', ndf=21, linear=False):
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
                n = re.findall("(\d*\.\d+|-\d*\.\d+)", line)

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

    eigval, eigvec = read_freq('vibs.molden', geo.ndf, False)

    '''Check that the Hessian is mass weighted '''

    eigval = np.sqrt(eigval.astype(np.float)/ph.amu)
    eigvec = eigvec.astype(np.float)

    N = np.matmul(np.transpose(eigvec[:, 6:]), eigvec[:, 6:])

    nc = geo.ndf - 6
    inmod = np.asarray(range(6, geo.ndf, 1))

    # for i in range(nc):
     #   print(N[i, i])
    alpha = np.zeros(nc)
    qmode = np.zeros(nc)
    pmode = np.zeros(nc)
    r_c = np.zeros((3, geo.natoms))
    p_c = np.zeros((3, geo.natoms))

    ZPE=0
    for i in range(0, nc):
        index = inmod[i]

       # omega = np.sqrt(float(eigval[index]) / ph.amu)

       # alpha[i] = 0.5 * omega
        alpha[i]=float(eigval[index])/(2.0)
        rq = RandomGaussianSample()
        rp = RandomGaussianSample()

        qmode[i] = 0.5 * rq * (np.sqrt(alpha[i])) ** (-1)
        pmode[i] = np.sqrt(alpha[i]) * rp

        ZPE += eigval[index]
    print('ZPE:',ZPE)

    # R = np.zeros((3, geo.natoms))
    # P = np.zeros_like(R)
    # L = 0
    # for i in range(geo.natoms):
    #     for j in range(3):
    #         R[j, i] = 0
    #         P[j, i] = 0
    #         for k in range(nc):
    #
    #             R[j, i] += eigvec[L, inmod[k]] * qmode[k]
    #             P[j, i] += eigvec[L, inmod[k]] * pmode[k]
    #         L = L + 1
    #
    # for i in range(geo.natoms):
    #     for j in range(3):
    #         R[j, i] = R[j, i] / np.sqrt(geo.masses[i])
    #         P[j, i] = P[j, i] * np.sqrt(geo.masses[i])
    # k = 0
    # vel = np.zeros_like(P)
    # for i in range(geo.natoms):
    #     for j in range(3):
    #         R[j, i] = R[j, i]
    #         vel[j, i] = P[j, i] / geo.masses[i]
    #
    # RCM = np.zeros(3)
    # VelCM = np.zeros(3)
    # TotMass = 0
    # for i in range(geo.natoms):
    #     TotMass += geo.masses[i]
    #     for j in range(3):
    #         RCM[j] += R[j, i] * geo.masses[i]
    #         VelCM[j] += geo.masses[i] * vel[j, i]
    #
    # RCM = RCM / TotMass
    # VelCM = VelCM / TotMass
    #
    # for i in range(geo.natoms):
    #     for j in range(3):
    #         # R[j,i]-=RCM[j]
    #         vel[j, i] -= VelCM[j]
    #         P[j, i] = vel[j, i] * geo.masses[i]

    for i in range(nc):
        ic = 0
        index = inmod[i]
        for n in range(0, geo.natoms):
            for j in range(0, 3):
                r_c[j, n] = r_c[j, n] + qmode[i] * eigvec[ic, index]
                p_c[j, n] = p_c[j, n] + pmode[i] * eigvec[ic, index]
                ic += 1

    p = np.zeros(geo.ndf)
    q = np.zeros(geo.ndf)
    idf = 0
    for i in range(0, geo.natoms):
        for j in range(0, 3):
            q[idf] = r_c[j, i] /np.sqrt(geo.masses[i])
            p[idf] = p_c[j, i]*np.sqrt(geo.masses[i])
            idf += 1


    ekin = 0
    k = 0
    for i in range(geo.natoms):
        for j in range(3):
            ekin = ekin + p[k] * p[k] / geo.masses[i]
            k += 1

    print('EKIN:', ekin)


    p_corr = WPrep(p)

    ekin=0
    k=0
    for i in range(geo.natoms):
        for j in range(3):
            ekin=ekin+0.5*p_corr[k]*p_corr[k]/geo.masses[i]
            k+=1

    print('EKIN_corr:',ekin)
    return q, p_corr


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
            p[j, i] = y[idf]/geo.masses[i]
            idf += 1

    totmass = np.sum(geo.masses)




    VelCM=np.zeros(3)
    for i in range(geo.natoms):
        for j in range(3):
            VelCM[j]+=geo.masses[i]*p[j,i]

    VelCM=VelCM/totmass
    print(VelCM)
    p_corr = np.zeros(np.shape(p))
    y_fin=np.zeros(geo.ndf)
    for i in range(geo.natoms):
        for j in range(3):

            p[j,i]=p[j,i]-VelCM[j]
            p_corr[j,i] =p[j,i]*geo.masses[i]

    idf=0
    for i in range(geo.natoms):
        for j in range(3):
            y_fin[idf] = p_corr[j, i]
            idf += 1
    return y_fin
