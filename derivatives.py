import numpy as np
import input


def der(q, p, s, d, epot, grad, nac, geo):
    dyn = input.initdyn()

    kinE = 0.5 * np.sum(p ** 2)  # Kinetic energy definition
    q_dot = p  # Derivative of position in WP frame
    p_dot = np.zeros(geo.ndf)

    for n in range(geo.ndf):
        count = np.double(0.0)
        for i in range(dyn.nstates):
            count = count - np.double(np.conj(d[i]) * d[i]) * grad[i, n]
        for i in range(dyn.nstates):
            for j in range(dyn.nstates):
                if i != j:
                    count += np.double(
                        np.conj(d[i]) * d[j] * np.exp(1.0j * (s[j] - s[i]), dtype=np.complex128)) * (
                                        epot[i] - epot[j]) * nac[j, i, n]
        p_dot[n] = count

    tmp = 0
    s_dot = np.zeros(dyn.nstates)  # Derivative of S (action)
    for i in range(geo.ndf):
        tmp += 0.5 * (q_dot[i] * p[i] - q[i] * p_dot[i])
    for j in range(dyn.nstates):
        s_dot[j] = tmp - (kinE + epot[j])

    d_dot = np.zeros(dyn.nstates, dtype=np.complex128)

    for i in range(dyn.nstates):
        for j in range(dyn.nstates):
            if i != j:
                nac1d = 0.0
                for n in range(geo.ndf):
                    nac1d += q_dot[n] * nac[j,i, n]
                d_dot[i] = d_dot[i] - np.complex128(nac1d) * d[j] * np.exp(1.0j * (s[j] - s[i]), dtype=np.complex128)

    return q_dot, p_dot, s_dot, d_dot
