import numpy as np
from scipy.linalg import expm
import src.couplings as cp
import src.initialize_traj as ip
import src.abinitio as ab
from src.geometry import initgeom

def magnus_2(H0, H1, dt):
    ndim = np.size(H0[:, 0])
    Hav = np.sum(np.sum(H0 + H1)) / np.double(2.0 * ndim)

    Htr = np.zeros((ndim, ndim),dtype=np.complex128)
    for i in range(ndim):
        Htr[i, i] = Hav
    a0 = (H1 + H0) / 2.0 - Htr
    W1 = dt * a0
    magH = expm(W1) * np.exp(Hav * dt)
    return magH


def velocityverlet(T, timestep,NN):
    geo=initgeom()
    ii = complex(0,1.00)
    magnus_slice = 20
    nst = T.nstates
    M = T.getmassall_traj()
    R0 = T.getposition_traj()
    P0 = T.getmomentum_traj()
    V0 = T.getvelocity_traj()
    FO = T.get_traj_force()
    E0 = T.getpotential_traj()
    A0 = T.getamplitude_traj()

    HE_0 = np.zeros((nst, nst),dtype=np.complex128)
    for n1 in range(nst):
        HE_0[n1, n1] = T.getpotential_traj_i(n1)
        for n2 in range(n1 + 1, nst):
            HE_0[n1, n2] = -ii * cp.coupdotvel(T, n1,n2)
            HE_0[n2, n1] = -HE_0[n1, n2]

    nslice = magnus_slice

    Ab = A0
    F0 = 0.0
    for i in range(1, nslice + 1):
        dt = timestep / nslice
        A1 = np.matmul(magnus_2(-ii * HE_0, -ii * HE_0, dt), Ab)
        Ab = A1
        T.stateAmpE = A1
        F0 += T.get_traj_force() / nslice

    T.stateAmpE = A0
    es0 = np.zeros(nst)
    fs0 = np.zeros((T.ndim, nst))
    cs0 = np.zeros((T.ndim, nst, nst))

    for i in range(nst):
        es0[i] = T.getpotential_traj_i(i)
        fs0[:, i] = T.get_traj_force()
        for j in range(nst):
            cs0[:, i, j] = T.getcoupling_traj(i, j)
    T.phase += timestep / 2.0 * ip.phasedot(T)

    R1 = R0 + timestep * V0 + timestep ** 2 / 2.0 * F0 / M
    P1 = P0 + timestep * F0

    T.setposition_traj(R1)

    T.setmomentum_traj(P1)

    pes, der = ab.inp_out(NN, 0, geo,T)

    T.setderivs_traj(der)
    F1 = T.get_traj_force()
    T.setpotential_traj(pes)
    es1 = np.zeros(nst)
    fs1 = np.zeros((T.ndim, nst))
    cs1 = np.zeros((T.ndim, nst, nst))

    for i in range(nst):
        es1[i] = T.getpotential_traj_i(i)
        fs1[:, i] = T.get_traj_force()
        for j in range(nst):
            cs1[:, i, j] = T.getcoupling_traj(i, j)

    HE_1 = np.zeros_like(HE_0,dtype=np.complex128)

    for n1 in range(nst):
        HE_1[n1, n1] = T.getpotential_traj_i(n1)
        for n2 in range(n1 + 1, nst):
            HE_1[n1, n2] = -ii * cp.coupdotvel(T, n1,n2)
            HE_1[n2, n1] = -HE_1[n1, n2]

    nslice = magnus_slice
    Ab = A0
    F1 = 0.0
    for n in range(1, nslice + 1):
        dt = timestep / nslice

        f_b = (n - 0.5) / np.double(float(nslice))
        HE_b = (1.0 - f_b) * HE_0 + f_b * HE_1
        esb = (1.0 - f_b) * es0 + f_b * es1
        fsb = (1.0 - f_b) * fs0 + f_b * fs1
        csb = (1.0 - f_b) * cs0 + f_b * cs1
        A1 = np.matmul(magnus_2(-ii * HE_b, -ii * HE_b, dt), Ab)
        Ab = A1
        T.stateAmpE = Ab
        T.HE = HE_1
        fb = ip.compforce(T, A1, fsb, esb, csb)
        F1 += fb * dt / timestep

    P1 = P0 + timestep * F1

    T.setamplitudes_traj(Ab)
    T.setmomentum_traj(P1)
    T.phase += timestep / 2.0 * ip.phasedot(T)

    T.setphases_traj(T.phase)
    return T
