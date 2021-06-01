import numpy as np
import input
from src.overlaps import coupdotvel as  cpdot


def der(T):
    ndf = T.ndim
    nstates = T.nstates
    q = T.getposition_traj()
    p = T.getmomentum_traj()
    s = T.getphases_traj()
    a = T.getamplitude_traj()

    d = T.getd_traj()
    masses = T.getmassall_traj()
    epot = T.getpoten_traj()
    grad = np.zeros((nstates, ndf))
    nac = np.zeros((nstates, nstates, ndf))
    for i in range(nstates):
        grad[i, :] = T.getforce_traj(i)
        for j in range(nstates):
            nac[i, j, :] = T.getcoupling_traj(i, j)

    kinE = 0.5 * np.sum(p ** 2 / masses)  # Kinetic energy definition
    q_dot = p / masses  # Derivative of position in WP frame
    p_dot = np.zeros(ndf)




    p_dot = T.get_traj_force()

    # for n in range(ndf):
    #     count = np.double(0.0)
    #     for i in range(nstates):
    #         count = count - np.double(np.conj(d[i]) * d[i]) * grad[i, n]
    #     for i in range(nstates):
    #         for j in range(nstates):
    #             if i != j:
    #                 count += np.double(
    #                     np.conj(d[i]) * d[j] * np.exp(1.0j * (s[j] - s[i]), dtype=np.complex128)) * (
    #                                  epot[i] - epot[j]) * nac[j, i, n]
    #     p_dot[n] = count
    #
    tmp = 0
    s_dot = np.zeros(nstates)  # Derivative of S (action)
    for i in range(ndf):
        tmp += 0.5 * (q_dot[i] * p[i] - q[i] * p_dot[i])
    for j in range(nstates):
        s_dot[j] = tmp - (kinE + epot[j])

    d_dot = np.zeros(nstates, dtype=np.complex128)

    for i in range(nstates):
        for j in range(nstates):
            if i != j:
                nac1d = 0.0
                for n in range(ndf):
                    nac1d += q_dot[n] * nac[j, i, n]
                d_dot[i] = d_dot[i] - np.complex128(nac1d) * d[j] * np.exp(1.0j * (s[j] - s[i]), dtype=np.complex128)

    a_dot = np.zeros(nstates, dtype=np.complex128)
    term = 0.0
    for i in range(nstates):
        for j in range(nstates):
            if i == j:
                term = T.getpotential_traj_i(i)
            elif i != j:
                term = -1j * cpdot(T, i, j)
            a_dot[i] = a_dot[i] + term * a[j]
    a_dot = -1j * a_dot

    print('a_dot',a_dot)
    print('a_dot_weird:', d_dot*np.exp(1j*s_dot))
    print('sdot:',s_dot)
    print('ddot:',d_dot)
    return q_dot, p_dot



