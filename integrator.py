from derivatives import der
import input
import numpy as np
from Constants import physconst


def rk_method_4(q, p, s, d, epot, grad, nac, time, dt):
    geo = input.initgeom()
    dyn = input.initdyn()
    ph = physconst()
    ekin_tr = 0.0
    print(dt)
    '''firsty define the total energy and the d norm'''
    dnorm = np.sum(np.abs(np.conj(d) * d))
    weight_d = ((d * np.conj(d)) / dnorm)

    epot_av = np.sum(weight_d * epot)
    print('Potential energy before starting: ', epot_av)

    ekin = np.sum(0.5 * p ** 2)
    print('Kinetic energy: ', ekin)

    tot_en = ekin + epot_av
    print('Total energy: ', tot_en)

    qrk = np.zeros((geo.ndf, 4))
    dq = np.zeros((geo.ndf, 4))
    prk = np.zeros((geo.ndf, 4))
    dp = np.zeros((geo.ndf, 4))
    trk = np.zeros(4)
    ds = np.zeros((dyn.nstates, 4))
    srk = np.zeros((dyn.nstates, 4))
    dd = np.zeros((dyn.nstates, 4),dtype=complex)
    drk = np.zeros((dyn.nstates, 4), dtype=complex)

    qrk[:, 0] = q
    prk[:, 0] = p
    drk[:, 0] = d
    srk[:, 0] = s
    trk[0] = time

    ind = np.asarray([0.5, 0.5, 1.0, 1.0])

    for j in range(3):
        dq[:, j], dp[:, j], ds[:, j], dd[:, j] = der(qrk[:, j], prk[:, j], srk[:, j], drk[:, j], epot, grad, nac)
        for i in range(geo.ndf):
            qrk[i, j + 1] = qrk[i, j] + ind[j] * dt * dq[i,j]
            prk[i, j + 1] = prk[i, j] + ind[j] * dt * dp[i,j]
        for i in range(dyn.nstates):
            srk[i, j + 1] = srk[i, j] + ind[j] * dt * ds[i,j]
            drk[i, j + 1] = drk[i, j] + ind[j] * dt * dd[i,j]

    dq[:, 3], dp[:, 3], ds[:, 3], dd[:, 3] = der(qrk[:, 3], prk[:, 3], srk[:, 3], drk[:, 3], epot, grad, nac)
    q_new = qrk[:, 0] + dt * (0.5 * dq[:, 0] + dq[:, 1] + dq[:, 2] + 0.5 * dq[:, 3]) / 3.00
    p_new = prk[:, 0] + dt * (0.5 * dp[:, 0] + dp[:, 1] + dp[:, 2] + 0.5 * dp[:, 3]) / 3.00
    s_new = srk[:, 0] + dt * (0.5 * ds[:, 0] + ds[:, 1] + ds[:, 2] + 0.5 * ds[:, 3]) / 3.00
    d_new = drk[:, 0] + dt * (0.5 * dd[:, 0] + dd[:, 1] + dd[:, 2] + 0.5 * dd[:, 3]) / 3.00

    totmass = np.sum(input.mass2au(geo.masses))

    q_old = qrk[:, 0]


    print('works up to here')
    return np.real(q_new),np.real(p_new),s_new,d_new
