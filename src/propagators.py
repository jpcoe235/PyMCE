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
import cmath
from src.matexp import matexpAI


def magnus(B, timestep):
    ii = np.complex128(0 + 1.00j)
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

        f_b = n / np.double(nslice + 1.00)

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
    ndim = np.size(H0[:, 0])
    Hav = np.complex128(0.0)
    for i in range(ndim):
        Hav = Hav + H0[i, i] + H1[i, i]

    Hav = Hav / np.double(2.0 * ndim)

    Htr = np.zeros((ndim, ndim), dtype=np.complex128)
    for i in range(ndim):
        Htr[i, i] = Hav
    a0 = (H1 + H0) / 2.0 - Htr
    W1 = dt * a0

    magH = expm(W1) * np.exp(Hav * dt, dtype=np.complex128)

    return magH


def velocityverlet(T, timestep, NN, calc1):
    geo = initgeom()
    geo2 = singlepart()
    ii = np.complex128(0 + 1.00j)
    magnus_slice = 20
    nst = T.nstates
    M = T.getmassall_traj()
    R0 = T.getposition_traj()
    P0 = T.getmomentum_traj()
    V0 = T.getvelocity_traj()
    A0 = T.getamplitude_traj()

    HE_0 = np.zeros((nst, nst), dtype=np.complex128)
    for n1 in range(nst):
        HE_0[n1, n1] = T.getpotential_traj_i(n1)
        for n2 in range(n1 + 1, nst):
            HE_0[n1, n2] = np.complex128(-ii * coupdotvel(T, n1, n2))
            HE_0[n2, n1] = -HE_0[n1, n2]

    nslice = magnus_slice

    Ab = A0
    F0 = 0.0
    for i in range(0, nslice):
        dt = timestep / np.double(nslice)
        if T.nstates > 1:
            A1 = np.matmul(magnus_2(-ii * HE_0, -ii * HE_0, dt), Ab, dtype=np.complex128)
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
    T.phase += timestep / 2.0 * T.phasedot()

    R1 = R0 + timestep * V0 + timestep ** 2.0 / 2.00 * F0 / M

    P1 = P0 + timestep * F0

    T.setposition_traj(R1)
    T.setmomentum_traj(P1)

    if not calc1:
        pes, der = ab.inp_out(NN, 0, geo, T)
    else:
        pes = np.sum(0.5 * geo2.K * T.getposition_traj() ** 2)
        der = np.zeros(3)
        for i in range(3):
            der[i] = -geo2.K * T.getposition_traj()[i]

    T.setderivs_traj(der)
    T.setpotential_traj(pes)
    es1 = np.zeros(nst)
    fs1 = np.zeros((T.ndim, nst))
    cs1 = np.zeros((T.ndim, nst, nst))

    for i in range(nst):
        es1[i] = T.getpotential_traj_i(i)
        fs1[:, i] = T.getforce_traj(i)
        for j in range(nst):
            cs1[:, i, j] = T.getcoupling_traj(i, j)

    HE_1 = np.zeros_like(HE_0, dtype=np.complex128)

    for n1 in range(nst):
        HE_1[n1, n1] = T.getpotential_traj_i(n1)
        for n2 in range(n1 + 1, nst):
            HE_1[n1, n2] = np.complex128(-ii * coupdotvel(T, n1, n2))
            HE_1[n2, n1] = -HE_1[n1, n2]

    nslice = magnus_slice
    Ab = A0
    F1 = 0.0
    for n in range(1, nslice + 1):
        dt = timestep / np.double(float(nslice))

        f_b = (n - 0.5) / np.double(float(nslice))
        HE_b = (1.0 - f_b) * HE_0 + f_b * HE_1
        esb = (1.0 - f_b) * es0 + f_b * es1
        fsb = (1.0 - f_b) * fs0 + f_b * fs1
        csb = (1.0 - f_b) * cs0 + f_b * cs1
        if T.nstates > 1:
            A1 = np.matmul(magnus_2(-ii * HE_b, -ii * HE_b, dt), Ab, dtype=np.complex128)
        else:
            A1 = magnus_2(-ii * HE_b, -ii * HE_b, dt) * Ab

        Ab = A1

        T.setamplitudes_traj(A1)
        T.HE = HE_1
        fb = T.compforce(A1, fsb, esb, csb)
        F1 += fb / np.double(float(nslice))

    P1 = P0 + timestep * F1
    T.setmomentum_traj(P1)
    T.phase += timestep / 2.0 * T.phasedot()

    T.setphases_traj(T.phase)
    T.setoldpos_traj(R0)
    T.setoldmom_traj(P0)
    return T
