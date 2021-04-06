import numpy as np
import input


def derivatives(tr):
    dyn = input.initdyn()
    geo = input.initgeom()
    kinE = 0.5 * np.sum(tr.p ** 2)  # Kinetic energy definition
    q_dot = tr.p  # Derivative of position in WP frame
    p_dot = np.zeros(geo.ndf,dtype=complex)
    count = 0.0  # Derivative of momentum in WP frame
    for n in range(geo.ndf):
        for i in range(dyn.nstates):
            count += -np.conj(tr.d[i]) * tr.d[i] * tr.grad[i, n]
            for j in range(dyn.nstates):
                if i != j:
                    count += np.conj(tr.d[i]) * tr.d[j] * tr.nac[i, j, n] * np.exp(1j) * (tr.s[j] - tr.s[i]) * (
                                tr.epot[j] - tr.epot[i])
        p_dot[n] = count
    tmp = 0
    s_dot = np.zeros(dyn.nstates,dtype=complex)  # Derivative of S (action)
    for i in range(geo.ndf):
        tmp += 0.5 * (q_dot[i] * tr.p[i] - tr.q[i] * p_dot[i])
    for j in range(dyn.nstates):
        s_dot[j] = tmp - (kinE + tr.epot[j])

    d_dot = np.zeros(dyn.nstates, dtype=complex)
    for i in range(dyn.nstates):
        for j in range(dyn.nstates):
            if i != j:
                nac1d = 0.0
                for n in range(geo.ndf):
                    nac1d += q_dot[n] * tr.nac[i, j, n]
                d_dot[i] = d_dot[i] - complex(nac1d) * tr.d[j] * np.exp(1j) * (tr.s[j] - tr.s[i])

    tr.dpdt = p_dot
    tr.dqdt = q_dot
    tr.dsdt = s_dot
    tr.dddt = d_dot
    return tr
