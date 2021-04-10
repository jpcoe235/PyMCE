from derivatives import der
import input
import numpy as np
from Constants import physconst
import molpro_call_read


def rk_method_4(q, p, s, d, epot, grad, nac, time, dt, ekin_tr, nstep,geo):

    dyn = input.initdyn()
    ph = physconst()

    '''firsty define the total energy and the d norm'''
    dnorm = 0.00

    for i in range(dyn.nstates):
        dnorm = dnorm + np.abs(np.conj(d[i]) * d[i])

    epot_av = 0.000
    for i in range(dyn.nstates):
        epot_av += np.real(np.conj(d[i]) * d[i]) * epot[i] / dnorm

    ekin = np.sum(0.5 * p ** 2)

    tot_en = ekin + epot_av + ekin_tr

    print('Total Energy:', tot_en)

    qrk = np.zeros((geo.ndf, 4))
    dq = np.zeros((geo.ndf, 4))
    prk = np.zeros((geo.ndf, 4))
    dp = np.zeros((geo.ndf, 4))
    trk = np.zeros(4)
    ds = np.zeros((dyn.nstates, 4))
    srk = np.zeros((dyn.nstates, 4))
    dd = np.zeros((dyn.nstates, 4), dtype=np.complex128)
    drk = np.zeros((dyn.nstates, 4), dtype=np.complex128)

    qrk[:, 0] = q
    prk[:, 0] = p
    drk[:, 0] = d
    srk[:, 0] = s
    trk[0] = time

    ind = np.asarray([0.5, 0.5, 1.0, 1.0])

    for j in range(3):

        dq[:, j], dp[:, j], ds[:, j], dd[:, j] = der(qrk[:, j], prk[:, j], srk[:, j], drk[:, j], epot, grad, nac,geo)

        for i in range(geo.ndf):
            qrk[i, j + 1] = qrk[i, 0] + ind[j] * dt * dq[i, j]
            prk[i, j + 1] = prk[i, 0] + ind[j] * dt * dp[i, j]
        for i in range(dyn.nstates):
            srk[i, j + 1] = srk[i, 0] + ind[j] * dt * ds[i, j]
            drk[i, j + 1] = drk[i, 0] + ind[j] * dt * dd[i, j]
        epot, grad, nac = molpro_call_read.inp_out(nstep, j+1, qrk[:, j + 1],geo)
    dq[:, 3], dp[:, 3], ds[:, 3], dd[:, 3] = der(qrk[:, 3], prk[:, 3], srk[:, 3], drk[:, 3], epot, grad, nac,geo)
    q_new = qrk[:, 0] + dt * (0.5 * dq[:, 0] + dq[:, 1] + dq[:, 2] + 0.5 * dq[:, 3]) / 3.00
    p_new = prk[:, 0] + dt * (0.5 * dp[:, 0] + dp[:, 1] + dp[:, 2] + 0.5 * dp[:, 3]) / 3.00
    s_new = srk[:, 0] + dt * (0.5 * ds[:, 0] + ds[:, 1] + ds[:, 2] + 0.5 * ds[:, 3]) / 3.00
    d_new = drk[:, 0] + dt * (0.5 * dd[:, 0] + dd[:, 1] + dd[:, 2] + 0.5 * dd[:, 3]) / 3.00

    totmass = np.sum(geo.masses)

    q_old = qrk[:, 0]
    idf = 0
    q_old_ref = np.zeros((3, geo.natoms))
    q_new_ref = np.zeros((3, geo.natoms))
    p_old_ref = np.zeros((3, geo.natoms))
    p_new_ref = np.zeros((3, geo.natoms))
    for i in range(geo.natoms):
        for j in range(3):
            q_new_ref[j, i] = q_new[idf]
            q_old_ref[j, i] = q_old[idf]
            p_old_ref[j, i] = q_old[idf]
            p_new_ref[j, i] = p_new[idf]
            idf += 1
    dr_com = np.zeros(3)
    v_com = np.zeros_like(dr_com)

    for i in range(geo.natoms):
        for j in range(3):
            dr_com[j] += np.sqrt(geo.masses[i]) * (q_new_ref[j, i] - q_old_ref[j, i])/totmass

    for i in range(3):
        v_com[j] = dr_com[j] / dt

    # v_com = v_com / totmass

    idf = 0
    for i in range(geo.natoms):
        for j in range(3):
            q[idf] = q_new_ref[j, i] - dr_com[j] * np.sqrt(geo.masses[i])
            p[idf] = p_new_ref[j, i] - v_com[j] * np.sqrt(geo.masses[i])
            idf += 1

    for i in range(3):
        ekin_tr = ekin_tr + 0.5 * (np.sqrt(totmass) * v_com[i]) ** 2.0

    a = d * np.exp(1j * s)

    return a, dnorm, q, p, s_new, d_new, ekin_tr
