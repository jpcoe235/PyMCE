import numpy as np
import input


def der(q, p, s, d, epot, grad, nac):

    dyn = input.initdyn()
    geo = input.initgeom()
    kinE = 0.5 * np.sum(p ** 2)  # Kinetic energy definition
    q_dot = p  # Derivative of position in WP frame
    p_dot = np.zeros(geo.ndf)
    count = 0.0  # Derivative of momentum in WP frame

    for n in range(geo.ndf):
        for i in range(dyn.nstates):
            count += -np.real(np.conj(d[i]) * d[i] * grad[i, n])
            for j in range(dyn.nstates):
                if i != j:
                    count += np.real(np.conj(d[i]) * d[j] * np.exp((0+1.0j) * (s[j] - s[i]))) * (
                            epot[j] - epot[i]) * nac[j, i, n]
        p_dot[n] = count
    tmp = 0
    s_dot = np.zeros(dyn.nstates)  # Derivative of S (action)
    for i in range(geo.ndf):
        tmp += 0.5 * (q_dot[i] * p[i] - q[i] * p_dot[i])
    for j in range(dyn.nstates):
        s_dot[j] = tmp - (kinE + epot[j])

    d_dot = np.zeros(dyn.nstates, dtype=complex)
    for i in range(dyn.nstates):
        for j in range(dyn.nstates):
            if i != j:
                nac1d = 0.0
                for n in range(geo.ndf):
                    nac1d += q_dot[n] * nac[j, i, n]
                d_dot[i] = d_dot[i] - np.complex(nac1d) * d[j] * np.exp(1.0j * (s[j] - s[i]))


    return q_dot, np.real(p_dot), s_dot, d_dot
