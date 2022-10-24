import numpy as np
from scipy.linalg import expm
import src.initialize_traj as ip
import src.abinitio as ab
from src.geometry import initgeom
from src.geometry import singlepart
from src.overlaps import coupdotvel
from Wigner_dist import WPrep
from src import bundle
from src import buildhs
from expomol import pade_mol as pm
import decimal as dD
import cmath
from copy import copy
import os
from src.matexp import matexpAI
from src.derivatives import der
from scipy.integrate import solve_ivp
from src.overlaps_wf import overlap as ovwf
from src.mldreader import readingmolden
from src.abinitio import ab_par


def magnus(B, timestep):
    ii = np.clongdouble(0 + 1.00j)
    magnus_slice = 20
    ntraj = B.ntraj
    t0 = B.time
    t1 = t0 + timestep

    Heff_0 = B.Heff

    print('Hamilt_0: ', Heff_0)
    nslice = magnus_slice

    for i in range(ntraj):
        B.Traj[i] = velocityverlet(B.Traj[i], timestep, 1)
        print(np.abs(B.Traj[i].stateAmpE) ** 2)

    B = buildhs.buildsandh(B)

    Heff_1 = B.Heff

    Heff_a = Heff_0

    C_tdt = B.getamps_bundle()

    for n in range(nslice + 2):
        if n == 0 or n == nslice:
            dt = 0.5 * timestep / (nslice + 1)
        else:
            dt = timestep / (nslice + 1)

        f_b = n / np.longdouble(nslice + 1.00)

        Heff_b = (1.00 - f_b) * Heff_0 + f_b * Heff_1

        if ntraj == 1:
            C_tdt[0] = np.exp(-ii * Heff_b[0, 0] * dt) * C_tdt[0]
        else:
            C_tdt = np.matmul(magnus_2(-ii * Heff_b, -ii * Heff_b, dt), C_tdt)

    B.setamps_bundle(C_tdt)
    print('Norm_amps_bundle :', B.get_calc_set_norm())
    print('AMps norm:', np.sum(np.abs(C_tdt) ** 2))
    B.setamps_bundle(C_tdt / cmath.sqrt(B.norm))
    B.settime_bundle(t1)
    return B


def magnus_2(H0, H1, dt):
    ab2 = ab_par()
    ndim = np.size(H0[:, 0])
    Hav = np.clongdouble(0.0)
    for i in range(ndim):
        Hav = Hav + H0[i, i] + H1[i, i]

    Hav = Hav / np.longdouble(2.0 * ndim)

    Htr = np.zeros((ndim, ndim), dtype=np.clongdouble)
    for i in range(ndim):
        Htr[i, i] = Hav
    a0 = (H1 + H0) / 2.0 - Htr
    W1 = dt * a0

    expW2 = pm.exp_pade(W1)

    magH = expW2 * np.exp(Hav * dt, dtype=np.clongdouble)

    return magH


def velocityverlet(T, timestep, NN, calc1, phasewf):
    ab2 = ab_par()
    geo = initgeom()
    geo2 = singlepart()
    ii = np.clongdouble(0 + 1.00j)
    magnus_slice = 10
    nst = T.nstates
    M = T.getmassall_traj()
    R0 = T.getposition_traj()
    P0 = T.getmomentum_traj()
    V0 = T.getvelocity_traj()
    A0 = T.getamplitude_traj()
    phase = T.getphase_traj()
    HE_0 = np.zeros((nst, nst), dtype=np.clongdouble)
    for n1 in range(nst):
        HE_0[n1, n1] = T.getpotential_traj_i(n1)
        for n2 in range(n1 + 1, nst):
            HE_0[n1, n2] = np.clongdouble(-ii * coupdotvel(T, n1, n2))
            HE_0[n2, n1] = -HE_0[n1, n2]

    nslice = magnus_slice

    Ab = A0
    F0 = 0.0
    for i in range(0, nslice):
        dt = timestep / np.longdouble(nslice)
        if T.nstates > 1:
            A1 = np.matmul(magnus_2(-ii * HE_0, -ii * HE_0, dt), Ab, dtype=np.clongdouble)
        else:
            A1 = magnus_2(-ii * HE_0, -ii * HE_0, dt) * Ab

        Ab = A1
        T.setamplitudes_traj(A1)
        F0 += T.get_traj_force() / nslice

    T.setamplitudes_traj(A0)
    es0 = np.zeros(nst)
    fs0 = np.zeros((T.ndim, nst))
    cs0 = np.zeros((T.ndim, nst, nst))

    for i in range(nst):
        es0[i] = T.getpotential_traj_i(i)
        fs0[:, i] = T.getforce_traj(i)
        for j in range(nst):
            cs0[:, i, j] = T.getcoupling_traj(i, j)
    phase += timestep / 2.0 * T.phasedot()

    R1 = R0 + timestep * V0 + timestep ** 2.0 / 2.00 * F0 / M
    print(V0)
    print(T.getmomentum_traj())
    print(T.getmassall_traj())
    P1 = P0 + timestep * F0
    oldcis = T.getcivecs()
    oldcoup = np.zeros((T.ndim, T.nstates - 1))
    T.setoldpos_traj(R0)
    T.setoldmom_traj(P0)
    for i in range(2, T.nstates + 1):
        oldcoup[:, i - 2] = T.getcoupling_traj(0, i - 1)
    T_try = copy(T)
    T.setposition_traj(R1)
    T.setmomentum_traj(P1)

    if not calc1:
        pes, der, cis, configs = ab.inp_out(NN, 0, geo, T)
        print('energies2= ', pes)
    else:
        pes = np.sum(0.5 * geo2.K * T.getposition_traj() ** 2)
        der = np.zeros(3)
        cis = 0
        configs = 0
        for i in range(3):
            der[i] = -geo2.K * T.getposition_traj()[i]

    T.setderivs_traj(der)
    T.setpotential_traj(pes)
    T.setcivecs(cis)
    T.setconfigs(configs)
    phasewf = T.getphasewf()
    print('phasewf= ', phasewf)
    ovs = np.zeros((T.nstates, T.nstates))

    for n1 in range(T.nstates):
        for n2 in range(T.nstates):
            ovs[n1, n2] = (np.dot(T.getcivecs()[:, n1], oldcis[:, n2]))
    print(ovs[0, 0])
    print(ovs[1, 1])
    print(abs(ovs[0, 0]) + abs(ovs[1, 1]), abs(ovs[0, 1]) + abs(ovs[1, 0]))
    # if abs(ovs[0, 0]) + abs(ovs[1, 1]) < abs(ovs[0, 1]) + abs(ovs[1, 0]):
    #     print('Trying to reduce timestep')
    #     T.setposition_traj(T.getoldpos_traj())
    #     T.setmomentum_traj(T.getoldmom_traj())
    #     os.system('cp /home/AndresMoreno/wfu/002.molpro.wfu /home/AndresMoreno/wfu/003.molpro.wfu')
    #     os.system('cp 002.molpro.wfu 003.molpro.wfu')
    #     pes, der, cis, configs = ab.inp_out(NN, 0, geo, T)
    #     T.setderivs_traj(der)
    #     T.setpotential_traj(pes)
    #     T.setcivecs(cis)
    #     print('momentum stored: ', T.getmomentum_traj()[0])
    #     print(T.getkineticlass() + T.getpotential_traj())
    #
    #     timestep = timestep / 2.00
    #     for ts in range(2):
    #         T = velocityverlet(T, timestep, NN+ts, calc1, phasewf)
    #     print('returning to the main routine')
    #     return T

    print(T.getcivecs()[:, 1])
    print(oldcis[:, 1])
    # ovs = np.zeros(T.nstates - 1)
    # for i in range(1,T.nstates):
    #     print(T.getcoupling_traj(0,i))
    #     ovs[i-1] = np.dot(T.getcoupling_traj(0,i), oldcoup[:, i-1])/np.abs(np.dot(T.getcoupling_traj(0,i), oldcoup[:, i-1]))
    #     print('ovs=',ovs)
    # for i in range(T.nstates):
    #     ovs = np.dot(T.getcivecs()[:, i], oldcis[:, i])
    #     print('Ovs= ',ovs)
    #     if abs(ovs) > 0.7:
    #         phasewf = phasewf * np.dot(T.getcivecs()[:, i], oldcis[:, i]) / np.abs(
    #             np.dot(T.getcivecs()[:, i], oldcis[:, i]))
    #     else:
    #         print('STEP TO CHECK; CI OVERLAP IS WRONG')
    #         changingindex = np.ones(np.size(T_try.configs)).astype(int)
    #         megamatrix1, M1, syms1 = readingmolden(T_try.getfilecalc())
    #         megamatrix2, M2, syms2 = readingmolden(T.getfilecalc())
    #
    #         syms1 = (syms1[ab2.closed_orb:] - ab2.closed_orb).astype(int)
    #         syms2 = (syms2[ab2.closed_orb:] - ab2.closed_orb).astype(int)
    #
    #         for nc in range(np.size(T_try.configs)):
    #             cfg = T.configs[nc]
    #             newcfg = ''
    #             for nstring in range(len(cfg)):
    #                 newcfg += cfg[syms2[nstring]]
    #             index1 = T_try.configs.index(newcfg)
    #             changingindex[nc] = int(index1)
    #             print(index1)
    #         newCIs = T.getcivecs()[changingindex, :]
    #         configs = T.getconfigs()
    #         configs = [configs[i] for i in changingindex]
    #         #T.setconfigs(configs)
    #         signs = np.ones_like(newCIs[:, 0])
    #         count = 0
    #         for conf in configs:
    #             if conf.replace('2', '0').replace('0', '') == T_try.configs[count].replace('2', '0').replace('0', ''):
    #                 signs[count] = 1
    #             else:
    #                 signs[count] = -1
    #             count = count + 1
    #         print(oldcis[:, i])
    #         print(newCIs[:, i])
    #         print('changed CIvectors')
    #
    #         ovs_2 = np.dot(newCIs[:, i], signs * oldcis[:, i])
    #         print('newoverlap: ', ovs_2)
    #         if (abs(ovs_2) < 0.9):
    #             ov_11, ov_22 = ovwf(T_try, T)
    #             print('wfoverlaps calculated after CI approach did not work', ov_11, ov_22)
    #             phasewf = phasewf * ov_11 / abs(ov_11) * ov_22 / abs(ov_22)
    #         else:
    #             phasewf = phasewf * ovs_2 / abs(ovs_2)
    #             print('changed phase after reordering CI vectors')

    derivs = np.zeros((T.ndim, T.nstates, T.nstates))
    count = 0
    for n1 in range(T.nstates):
        for n2 in range(T.nstates):
            if n1 != n2:
                derivs[:, n1, n2] = T.getcoupling_traj(n1, n2)

            else:
                derivs[:, n1, n1] = T.getforce_traj(n1)
    T.setphasewf(phasewf)
    T.setderivs_traj(derivs)
    es1 = np.zeros(nst)
    fs1 = np.zeros((T.ndim, nst))
    cs1 = np.zeros((T.ndim, nst, nst))

    for i in range(nst):
        es1[i] = T.getpotential_traj_i(i)
        fs1[:, i] = T.getforce_traj(i)
        for j in range(nst):
            cs1[:, i, j] = T.getcoupling_traj(i, j)

    for i in range(1, T.nstates):
        ovs = np.sign(np.dot(cs1[:, 0, i], cs0[:, 0, i]))
        cs1[:, :, i] = cs1[:, :, i] * ovs
        cs1[:, i, :] = cs1[:, i, :] * ovs

    for n1 in range(T.nstates):
        for n2 in range(T.nstates):
            if n1 != n2:
                derivs[:, n1, n2] = cs1[:, n1, n2]
            else:
                derivs[:, n1, n1] = fs1[:, n1]
    T.setderivs_traj(derivs)
    HE_1 = np.zeros_like(HE_0, dtype=np.clongdouble)

    for n1 in range(nst):
        HE_1[n1, n1] = T.getpotential_traj_i(n1)
        for n2 in range(n1 + 1, nst):
            HE_1[n1, n2] = np.clongdouble(-ii * coupdotvel(T, n1, n2))
            HE_1[n2, n1] = -HE_1[n1, n2]

    nslice = magnus_slice
    Ab = A0
    F1 = 0.0

    print('He_0=', HE_0)
    for n in range(1, nslice + 1):
        dt = timestep / np.longdouble(float(nslice))

        f_b = (n - 0.5) / np.longdouble(float(nslice))

        HE_b_dima = (n * HE_1 + (10 - n) * HE_0) * 0.1
        print(HE_b_dima)
        HE_b = (1.0 - f_b) * HE_0 + f_b * HE_1
        print(HE_b)
        esb = (1.0 - f_b) * es0 + f_b * es1
        fsb = (1.0 - f_b) * fs0 + f_b * fs1
        csb = (1.0 - f_b) * cs0 + f_b * cs1
        if T.nstates > 1:
            A1 = np.matmul(magnus_2(-ii * HE_b, -ii * HE_b, dt), Ab, dtype=np.clongdouble)
        else:
            A1 = magnus_2(-ii * HE_b, -ii * HE_b, dt) * Ab

        Ab = A1

        T.setamplitudes_traj(A1)
        T.HE = HE_1
        fb = T.compforce(A1, fsb, esb, csb)
        F1 += fb / np.longdouble(float(nslice))

    P1 = P0 + timestep * F1
    T.setmomentum_traj(P1)
    phase += timestep / 2.0 * T.phasedot()
    ICycle = np.floor(phase / (2.0000 * np.pi))
    phase = phase - 2.000 * np.pi * ICycle
    T.setphase_traj(phase)
    T.setoldpos_traj(R0)
    T.setoldmom_traj(P0)
    return T


def velocityverlet_dima(T, finaltime,timestep,NN, calc1, phasewf):
    np.set_printoptions(precision=32)
    time=0.00000
    f = open("N_mine.dat", 'w', buffering=1)
    ab2 = ab_par()
    geo = initgeom()
    geo2 = singlepart()
    ii = np.clongdouble(0.00 + 1.00000000000j)
    magnus_slice = 10
    nst = T.nstates
    M = T.getmassall_traj()
    R0 = T.getposition_traj()
    P0 = T.getmomentum_traj()
    V0 = T.getvelocity_traj()
    A0 = T.getamplitude_traj()
    HE_0 = np.zeros((nst, nst), dtype=np.clongdouble)

    for n1 in range(nst):
        HE_0[n1, n1] = T.getpotential_traj_i(n1)
        for n2 in range(n1 + 1, nst):
            HE_0[n1, n2] = np.clongdouble(-ii * np.sum(V0 * T.getcoupling_traj(n1, n2)))
            HE_0[n2, n1] = -HE_0[n1, n2]

    es0 = np.zeros(nst, dtype=np.longdouble)
    fs0 = np.zeros((T.ndim, nst), dtype=np.longdouble)
    cs0 = np.zeros((T.ndim, nst, nst), dtype=np.longdouble)

    for i in range(nst):
        es0[i] = T.getpotential_traj_i(i)
        fs0[:, i] = T.getforce_traj(i)
        for j in range(nst):
            cs0[:, i, j] = T.getcoupling_traj(i, j)

    F0 = T.compforce(A0, fs0, es0, cs0) / np.longdouble(10.000)
    while(time<=finaltime):

        nslice = magnus_slice
        print('He_0=', HE_0)
        Ab = np.matmul(magnus_2(-ii * HE_0, -ii * HE_0, timestep / np.longdouble(20.0000)), A0, dtype=np.clongdouble)

        F0 = T.compforce(A0, fs0, es0, cs0)/10.0000000
        for i in range(1, 10):
            print(i)
            dt = timestep / np.longdouble(nslice)
            if T.nstates > 1:
                A1 = np.matmul(magnus_2(-ii * HE_0, -ii * HE_0, timestep / np.longdouble(10.0000)), Ab,
                               dtype=np.clongdouble)
            else:
                A1 = magnus_2(-ii * HE_0, -ii * HE_0, dt) * Ab
            print('first prop', A1[0])
            Ab = A1
            F0 += T.compforce(A1, fs0, es0, cs0) / np.longdouble(10.000)
        print('first step A and Forces: ', F0)
        T.setamplitudes_traj(A0)

        R1 = R0 + timestep * V0 + timestep ** 2.0 / 2.00 * F0 / M
        P1 = P0 + timestep * F0
        V1 = P1 / M

        T.setoldpos_traj(R0)
        T.setoldmom_traj(P0)

        T_try = copy(T)
        T.setposition_traj(R1)
        T.setmomentum_traj(P1)

        pes, der, cis, configs = ab.inp_out(NN, 0, geo, T)
        print('energies2= ', pes)
        T.setderivs_traj(der)
        T.setpotential_traj(pes)
        T.setcivecs(cis)
        T.setconfigs(configs)
        phasewf = T.getphasewf()
        print('phasewf= ', phasewf)
        ovs = np.zeros((T.nstates, T.nstates))
        es1 = np.zeros(nst, dtype=np.longdouble)
        fs1 = np.zeros((T.ndim, nst), dtype=np.longdouble)
        cs1 = np.zeros((T.ndim, nst, nst), dtype=np.longdouble)

        for i in range(nst):
            es1[i] = T.getpotential_traj_i(i)
            fs1[:, i] = T.getforce_traj(i)
            for j in range(nst):
                cs1[:, i, j] = T.getcoupling_traj(i, j)

        for i in range(1, T.nstates):
            ovi = np.sum(cs1[:, 0, i] * cs0[:, 0, i]) / abs(np.sum(cs1[:, 0, i] * cs0[:, 0, i]))
            cs1[:, :, i] = cs1[:, :, i] * ovi
            cs1[:, i, :] = cs1[:, i, :] * ovi

            print('ovs:', ovi)
        derivs = np.zeros((T.ndim, nst, nst))
        for n1 in range(T.nstates):
            for n2 in range(T.nstates):
                if n1 != n2:
                    derivs[:, n1, n2] = cs1[:, n1, n2]
                else:
                    derivs[:, n1, n1] = fs1[:, n1]
        T.setderivs_traj(derivs)
        HE_1 = np.zeros_like(HE_0, dtype=np.clongdouble)

        for n1 in range(nst):
            HE_1[n1, n1] = T.getpotential_traj_i(n1)
            for n2 in range(n1 + 1, nst):
                HE_1[n1, n2] = np.clongdouble(-ii * np.sum(V1 * cs1[:, n1, n2], dtype=np.longdouble))
                HE_1[n2, n1] = -HE_1[n1, n2]
        print('He_0=' + str(HE_0) + '\n')
        print('He_1=' + str(HE_1) + '\n')
        Ab = np.matmul(magnus_2(-ii * HE_0, -ii * HE_0, timestep / np.longdouble(20.0)), A0, dtype=np.clongdouble)
        esb = np.longdouble(0.05) * es1 + np.longdouble(0.95) * es0
        fsb = np.longdouble(0.05) * fs1 + np.longdouble(0.95) * fs0
        csb = np.longdouble(0.05) * cs1 + np.longdouble(0.95) * cs0
        F1 = T.compforce(Ab, fsb, esb, csb) / np.longdouble(10.0)

        for n in range(1, 10):
            print(n)

            HE_b = (n * HE_1 + (10.0000 - n) * HE_0) * 0.1

            print(Ab)
            esb = (0.1 * n + np.longdouble(0.05)) * es1 + (np.longdouble(0.95) - n * 0.1) * es0
            fsb = (0.1 * n + np.longdouble(0.05)) * fs1 + (np.longdouble(0.95) - n * 0.1) * fs0
            csb = (0.1 * n + np.longdouble(0.05)) * cs1 + (np.longdouble(0.95) - n * 0.1) * cs0
            A1 = np.matmul(magnus_2(-ii * HE_b, -ii * HE_b, timestep / np.longdouble(10.0000000)), Ab, dtype=np.clongdouble)

            Ab = A1
            Fb = T.compforce(A1, fsb, esb, csb)
            F1 = F1 + Fb / np.longdouble(10.00)
            # T.HE = HE_1
        # if T.nstates > 1:
        #         A1 = np.matmul(magnus_2(-ii * HE_b, -ii * HE_b, dt), Ab, dtype=np.clongdouble)
        #     else:
        #         A1 = magnus_2(-ii * HE_b, -ii * HE_b, dt) * Ab
        #
        #     Ab = A1
        #
        #
        #     T.HE = HE_1
        #     fb = T.compforce(A1, fsb, esb, csb)
        #     F1 += fb / np.longdouble(float(nslice))
        A1 = np.matmul(magnus_2(-ii * HE_1, -ii * HE_1, timestep / 20.), Ab, dtype=np.clongdouble)
        print(A1)

        T.setamplitudes_traj(A1)
        P1 = P0 + timestep * F1
        T.setmomentum_traj(P1)
        T.setoldpos_traj(R0)
        T.setoldmom_traj(P0)
        f.write(str(time)+' '+str(abs(A1[0]) ** 2)+' '+str(abs(A1[1])**2)+'\n')
        time=time+timestep
        NN+=1
        R0=R1
        P0=P1
        A0=A1

        fs0 = fs1
        cs0 = cs1
        es0 = es1

        V0 = P1 / M
        F0 = T.compforce(A0, fs0, es0, cs0)
        HE_0 = HE_1

        print('final value of A:', abs(A1) ** 2)
    return T


def rk_method_4(T, dt, time, nstep):
    geo = initgeom()
    d = T.getd_traj()
    q = T.getposition_traj()
    p = T.getmomentum_traj()
    s = T.getphases_traj()
    a = T.getamplitude_traj()
    ii = np.clongdouble(0 + 1j)
    nslice = 20

    nstates = T.nstates
    epot = T.PotEn
    ndf = T.ndim
    grad = np.zeros((nstates, ndf))
    nacs = np.zeros((nstates, nstates, ndf))
    masses = T.getmassall_traj()
    for i in range(nstates):
        grad[i, :] = T.getforce_traj(i)
        for j in range(nstates):
            nacs[i, j, :] = T.getcoupling_traj(i, j)

    '''first define the total energy and the d norm'''
    dnorm = 0.00

    for i in range(nstates):
        dnorm = dnorm + np.abs(np.conj(d[i]) * d[i])

    epot_av = 0.000
    for i in range(nstates):
        epot_av += np.abs(np.conj(d[i]) * d[i]) * epot[i]

    ekin = np.sum(0.5 * p ** 2 / masses)

    tot_en = ekin + epot_av

    print('Total Energy:', tot_en)

    qk1 = dt * p / masses
    pk1 = dt * T.get_traj_force()

    d_dot = np.zeros(nstates, dtype=np.clongdouble)
    for i in range(nstates):
        for j in range(nstates):
            if i != j:
                nac1d = 0.0
                for n in range(ndf):
                    nac1d += p[n] / masses[n] * T.getcoupling_traj(j, i)[n]
                d_dot[i] = d_dot[i] - np.clongdouble(nac1d) * d[j] * np.exp(1.0j * (s[j] - s[i]), dtype=np.clongdouble)
    dk1 = dt * d_dot

    s_dot = np.zeros(nstates, dtype=np.clongdouble)  # Derivative of S (action)
    tmp = 0.0
    for i in range(ndf):
        tmp += 0.5 * (p[i] / masses[i] * p[i] - q[i] * T.get_traj_force()[i])
    for j in range(nstates):
        s_dot[j] = tmp - (0.5 * np.dot(p, p / masses) + epot[j])

    sk1 = dt * s_dot

    T.setposition_traj(q + 0.5 * qk1)
    T.setmomentum_traj(p + 0.5 * pk1)
    T.setd_traj(d + 0.5 * dk1)
    T.setphases_traj(s + 0.5 * sk1)

    amps = T.getd_traj() * np.exp(1j * T.getphases_traj(), dtype=np.clongdouble)
    T.setamplitudes_traj(amps)

    epot, derivs = ab.inp_out(nstep, 0, geo, T)
    T.setderivs_traj(derivs)
    T.setpotential_traj(epot)

    p2 = T.getmomentum_traj()
    q2 = T.getposition_traj()
    # a2 = T.getamplitude_traj()
    d2 = T.getd_traj()
    s2 = T.getphases_traj()

    qk2 = dt * p2 / masses
    pk2 = dt * T.get_traj_force()

    d_dot = np.zeros(nstates, dtype=np.clongdouble)
    for i in range(nstates):
        for j in range(nstates):
            if i != j:
                nac1d = 0.0
                for n in range(ndf):
                    nac1d += p2[n] / masses[n] * T.getcoupling_traj(j, i)[n]
                d_dot[i] = d_dot[i] - np.clongdouble(nac1d) * d2[j] * np.exp(1.0j * (s2[j] - s2[i]),
                                                                             dtype=np.clongdouble)
    dk2 = dt * d_dot

    s_dot = np.zeros(nstates, dtype=np.clongdouble)  # Derivative of S (action)
    tmp = 0.0
    for i in range(ndf):
        tmp += 0.5 * (p2[i] / masses[i] * p2[i] - q2[i] * T.get_traj_force()[i])
    for j in range(nstates):
        s_dot[j] = tmp - (0.5 * np.dot(p2, p2 / masses) + epot[j])

    sk2 = dt * s_dot

    T.setposition_traj(q + 0.5 * qk2)
    T.setmomentum_traj(p + 0.5 * pk2)
    T.setd_traj(d + 0.5 * dk2)
    T.setphases_traj(s + 0.5 * sk2)

    amps = T.getd_traj() * np.exp(1j * T.getphases_traj(), dtype=np.clongdouble)
    T.setamplitudes_traj(amps)

    epot, derivs = ab.inp_out(nstep, 0, geo, T)
    T.setderivs_traj(derivs)
    T.setpotential_traj(epot)

    p3 = T.getmomentum_traj()
    q3 = T.getposition_traj()
    d3 = T.getd_traj()
    s3 = T.getphases_traj()

    qk3 = dt * p3 / masses
    pk3 = dt * T.get_traj_force()

    d_dot = np.zeros(nstates, dtype=np.clongdouble)
    for i in range(nstates):
        for j in range(nstates):
            if i != j:
                nac1d = 0.0
                for n in range(ndf):
                    nac1d += p3[n] / masses[n] * T.getcoupling_traj(j, i)[n]
                d_dot[i] = d_dot[i] - np.clongdouble(nac1d) * d3[j] * np.exp(1.0j * (s3[j] - s3[i]),
                                                                             dtype=np.clongdouble)
    dk3 = dt * d_dot

    s_dot = np.zeros(nstates, dtype=np.clongdouble)  # Derivative of S (action)
    tmp = 0.0
    for i in range(ndf):
        tmp += 0.5 * (p3[i] / masses[i] * p3[i] - q3[i] * T.get_traj_force()[i])
    for j in range(nstates):
        s_dot[j] = tmp - (0.5 * np.dot(p3, p3 / masses) + epot[j])

    sk3 = dt * s_dot

    T.setposition_traj(q + qk3)
    T.setmomentum_traj(p + pk3)
    T.setd_traj(d + dk3)
    T.setphases_traj(s + sk3)

    amps = T.getd_traj() * np.exp(1j * T.getphases_traj(), dtype=np.clongdouble)
    T.setamplitudes_traj(amps)
    epot, derivs = ab.inp_out(nstep, 0, geo, T)
    T.setderivs_traj(derivs)
    T.setpotential_traj(epot)

    p4 = T.getmomentum_traj()
    q4 = T.getposition_traj()
    d4 = T.getd_traj()
    s4 = T.getphases_traj()

    qk4 = dt * p4 / masses
    pk4 = dt * T.get_traj_force()

    d_dot = np.zeros(nstates, dtype=np.clongdouble)
    for i in range(nstates):
        for j in range(nstates):
            if i != j:
                nac1d = 0.0
                for n in range(ndf):
                    nac1d += p4[n] / masses[n] * T.getcoupling_traj(j, i)[n]
                d_dot[i] = d_dot[i] - np.clongdouble(nac1d) * d4[j] * np.exp(1.0j * (s4[j] - s4[i]),
                                                                             dtype=np.clongdouble)
    dk4 = dt * d_dot

    s_dot = np.zeros(nstates, dtype=np.clongdouble)  # Derivative of S (action)
    tmp = 0.0
    for i in range(ndf):
        tmp += 0.5 * (p4[i] / masses[i] * p4[i] - q4[i] * T.get_traj_force()[i])
    for j in range(nstates):
        s_dot[j] = tmp - (0.5 * np.dot(p4, p4 / masses) + epot[j])

    sk4 = dt * s_dot

    T.setposition_traj(q + (qk1 + 2.000 * qk2 + 2.00 * qk3 + qk4) / 6.00)
    T.setmomentum_traj(p + (pk1 + 2.000 * pk2 + 2.00 * pk3 + pk4) / 6.00)
    T.setd_traj(d + (dk1 + 2.000 * dk2 + 2.00 * dk3 + dk4) / 6.00)
    T.setphases_traj(s + (sk1 + 2.000 * sk2 + 2.00 * sk3 + sk4) / 6.00)

    amps = T.getd_traj() * np.exp(1j * T.getphases_traj(), dtype=np.clongdouble)
    T.setamplitudes_traj(amps)

    return T


def qdot(T):
    q_dot = T.getmomentum_traj() / T.getmassall_traj()
    return q_dot


def pdot(T):
    p_dot = T.get_traj_force()
    return p_dot


def ddot(T):
    nstates = T.nstates
    d_dot = np.zeros(nstates, dtype=np.clongdouble)
    d = T.getd_traj()
    s = T.getphases_traj()
    p = T.getmomentum_traj()
    ndf = T.ndim
    masses = T.getmassall_traj()
    for i in range(nstates):
        for j in range(nstates):
            if i != j:
                nac1d = 0.0
                for n in range(ndf):
                    nac1d += p[n] / masses[n] * T.getcoupling_traj(j, i)[n]
                d_dot[i] = d_dot[i] - np.clongdouble(nac1d) * d[j] * np.exp(1.0j * (s[j] - s[i]), dtype=np.clongdouble)

    return d_dot


def sdot(T):
    nstates = T.nstates
    masses = T.getmassall_traj()
    ndf = T.ndim
    p = T.getmomentum_traj()
    q = T.getposition_traj()

    s_dot = np.zeros(nstates, dtype=np.clongdouble)  # Derivative of S (action)
    tmp = 0.0
    for i in range(ndf):
        tmp += 0.5 * (p[i] / masses[i] * p[i] - q[i] * T.get_traj_force()[i])
    for j in range(nstates):
        s_dot[j] = tmp - (0.5 * np.dot(p, p / masses) + T.getpotential_traj_i(j))

    return s_dot


def rkf45(T, dt, time, nstep):
    geo = initgeom()
    d = T.getd_traj()
    q = T.getposition_traj()
    p = T.getmomentum_traj()
    s = T.getphases_traj()
    A0 = T.getamplitude_traj()
    ii = np.clongdouble(0 + 1j)
    nstates = T.nstates
    masses = T.getmassall_traj()

    time = np.linspace(0, dt, 10)
    # Coefficients used to compute the independent variable argument of f

    a2 = 2.500000000000000e-01  # 1/4
    a3 = 3.750000000000000e-01  # 3/8
    a4 = 9.230769230769231e-01  # 12/13
    a5 = 1.000000000000000e+00  # 1
    a6 = 5.000000000000000e-01  # 1/2

    # Coefficients used to compute the dependent variable argument of f

    b21 = 2.500000000000000e-01  # 1/4
    b31 = 9.375000000000000e-02  # 3/32
    b32 = 2.812500000000000e-01  # 9/32
    b41 = 8.793809740555303e-01  # 1932/2197
    b42 = -3.277196176604461e+00  # -7200/2197
    b43 = 3.320892125625853e+00  # 7296/2197
    b51 = 2.032407407407407e+00  # 439/216
    b52 = -8.000000000000000e+00  # -8
    b53 = 7.173489278752436e+00  # 3680/513
    b54 = -2.058966861598441e-01  # -845/4104
    b61 = -2.962962962962963e-01  # -8/27
    b62 = 2.000000000000000e+00  # 2
    b63 = -1.381676413255361e+00  # -3544/2565
    b64 = 4.529727095516569e-01  # 1859/4104
    b65 = -2.750000000000000e-01  # -11/40

    # Coefficients used to compute local truncation error estimate.  These
    # come from subtracting a 4th order RK estimate from a 5th order RK
    # estimate.

    r1 = 2.777777777777778e-03  # 1/360
    r3 = -2.994152046783626e-02  # -128/4275
    r4 = -2.919989367357789e-02  # -2197/75240
    r5 = 2.000000000000000e-02  # 1/50
    r6 = 3.636363636363636e-02  # 2/55

    # Coefficients used to compute 4th order RK estimate

    c1 = 1.157407407407407e-01  # 25/216
    c3 = 5.489278752436647e-01  # 1408/2565
    c4 = 5.353313840155945e-01  # 2197/4104
    c5 = -2.000000000000000e-01  # -1/5

    n = len(time)

    a = 0.0
    b = 2.5
    t = a
    hmax = 1.25
    h = hmax
    hmin = 0.01
    tol = 1e-7
    tol2 = 0.000001
    T2 = T
    while t < b:
        T2 = T
        energy0 = T.getkineticlass() + T.getpotential_traj()
        print(h, energy0)
        if t + h > b:
            h = b - t

        qk1 = h * qdot(T)
        pk1 = h * pdot(T)
        dk1 = h * ddot(T)
        sk1 = h * sdot(T)

        T2.setposition_traj(q + qk1)
        T2.setmomentum_traj(p + pk1)
        T2.setd_traj(d + dk1)
        T2.setphases_traj(s + sk1)
        T2.setamplitudes_traj(T2.getd_traj() * np.exp(T2.getphases_traj() * 1j, dtype=np.clongdouble))

        epot, derivs = ab.inp_out(nstep, 0, geo, T2)
        T.setderivs_traj(derivs)
        T.setpotential_traj(epot)

        qk2 = h * (qdot(T) + b21 * qk1)
        pk2 = h * (pdot(T) + b21 * pk1)
        dk2 = h * (ddot(T) + b21 * dk1)
        sk2 = h * (sdot(T) + b21 * sk1)

        T2.setposition_traj(q + qk2)
        T2.setmomentum_traj(p + pk2)
        T2.setd_traj(d + dk2)
        T2.setphases_traj(s + sk2)
        T2.setamplitudes_traj(T2.getd_traj() * np.exp(T2.getphases_traj() * 1j, dtype=np.clongdouble))

        epot, derivs = ab.inp_out(nstep, 0, geo, T2)
        T.setderivs_traj(derivs)
        T.setpotential_traj(epot)

        qk3 = h * (qdot(T) + b31 * qk1 + b32 * qk2)
        pk3 = h * (pdot(T) + b31 * pk1 + b32 * pk2)
        dk3 = h * (ddot(T) + b31 * dk1 + b32 * dk2)
        sk3 = h * (sdot(T) + b31 * sk1 + b32 * sk2)

        T2.setposition_traj(q + qk3)
        T2.setmomentum_traj(p + pk3)
        T2.setd_traj(d + dk3)
        T2.setphases_traj(s + sk3)
        T2.setamplitudes_traj(T2.getd_traj() * np.exp(T2.getphases_traj() * 1j, dtype=np.clongdouble))

        epot, derivs = ab.inp_out(nstep, 0, geo, T2)
        T.setderivs_traj(derivs)
        T.setpotential_traj(epot)

        qk4 = h * (qdot(T) + b41 * qk1 + b42 * qk2 + b43 * qk3)
        pk4 = h * (pdot(T) + b41 * pk1 + b42 * pk2 + b43 * pk3)
        dk4 = h * (ddot(T) + b41 * dk1 + b42 * dk2 + b43 * dk3)
        sk4 = h * (sdot(T) + b41 * sk1 + b42 * sk2 + b43 * sk3)

        T2.setposition_traj(q + qk4)
        T2.setmomentum_traj(p + pk4)
        T2.setd_traj(d + dk4)
        T2.setphases_traj(s + sk4)
        T2.setamplitudes_traj(T2.getd_traj() * np.exp(T2.getphases_traj() * 1j, dtype=np.clongdouble))

        epot, derivs = ab.inp_out(nstep, 0, geo, T2)
        T.setderivs_traj(derivs)
        T.setpotential_traj(epot)

        qk5 = h * (qdot(T) + b51 * qk1 + b52 * qk2 + b53 * qk3 + b54 * qk4)
        pk5 = h * (pdot(T) + b51 * pk1 + b52 * pk2 + b53 * pk3 + b54 * pk4)
        dk5 = h * (ddot(T) + b51 * dk1 + b52 * dk2 + b53 * dk3 + b54 * dk4)
        sk5 = h * (sdot(T) + b51 * sk1 + b52 * sk2 + b53 * sk3 + b54 * sk4)

        T2.setposition_traj(q + qk5)
        T2.setmomentum_traj(p + pk5)
        T2.setd_traj(d + dk5)
        T2.setphases_traj(s + sk5)
        T2.setamplitudes_traj(T2.getd_traj() * np.exp(T2.getphases_traj() * 1j, dtype=np.clongdouble))

        epot, derivs = ab.inp_out(nstep, 0, geo, T2)
        T.setderivs_traj(derivs)
        T.setpotential_traj(epot)

        qk6 = h * (qdot(T) + b61 * qk1 + b62 * qk2 + b63 * qk3 + b64 * qk4 + b65 * qk5)
        pk6 = h * (pdot(T) + b61 * pk1 + b62 * pk2 + b63 * pk3 + b64 * pk4 + b65 * pk5)
        dk6 = h * (ddot(T) + b61 * dk1 + b62 * dk2 + b63 * dk3 + b64 * dk4 + b65 * dk5)
        sk6 = h * (sdot(T) + b61 * sk1 + b62 * sk2 + b63 * sk3 + b64 * sk4 + b65 * sk5)

        T2.setposition_traj(q + qk6)
        T2.setmomentum_traj(p + pk6)
        T2.setd_traj(d + dk6)
        T2.setphases_traj(s + sk6)
        T2.setamplitudes_traj(T2.getd_traj() * np.exp(T2.getphases_traj() * 1j, dtype=np.clongdouble))

        epot, derivs = ab.inp_out(nstep, 0, geo, T2)
        T.setderivs_traj(derivs)
        T.setpotential_traj(epot)

        rq = abs(r1 * qk1 + r3 * qk3 + r4 * qk4 + r5 * qk5 + r6 * qk6) / h
        rp = abs(r1 * qk1 + r3 * qk3 + r4 * qk4 + r5 * qk5 + r6 * qk6) / h

        energy2 = T2.getkineticlass() + T2.getpotential_traj()

        if len(np.shape(rq)) > 0:
            rq = max(rq)
            rp = max(rp)
        if abs(rp) <= tol2:
            t = t + h
            q = q + c1 * qk1 + c3 * qk3 + c4 * qk4 + c5 * qk5
            p = p + c1 * pk1 + c3 * pk3 + c4 * pk4 + c5 * pk5
            d = d + c1 * dk1 + c3 * dk3 + c4 * dk4 + c5 * dk5
            s = s + c1 * sk1 + c3 * sk3 + c4 * sk4 + c5 * sk5

            T.setposition_traj(q)
            T.setmomentum_traj(p)
            T.setamplitudes_traj(d * np.exp(s * 1j, dtype=np.clongdouble))

        h = h * min(max(0.84 * (tol2 / rp) ** 0.25, 0.1), 4.0)

        if h > hmax:
            h = hmax
        elif h < hmin:
            print("Error: stepsize should be smaller than %e." % hmin)
            break

    print(h, q[0])

    return T
