import numpy as np
from src import geometry as g


def coupdotvel(T, i, j):
    if i == j:
        coup = 0.0
    else:
        coup = np.dot(T.getvelocity_traj(), T.getcoupling_traj(i, j))

    return coup


def overlap_CG(x1, x2, p1, p2, a1, a2):
    dx = x1 - x2
    dp = p1 - p2

    real_part = (a1 * a2 * dx ** 2 + 0.25 * dp ** 2) / (a1 + a2)
    pref = np.sqrt(2.0 * np.sqrt(a1 * a2) / (a1 + a2))
    x_cent = (a1 * x1 + a2 * x2) / (a1 + a2)
    imag_part = (p1 * x1 - p2 * x2) - x_cent * dp
    sij = pref * np.exp(-real_part + 1j * imag_part)

    return sij


def overlap_d2x_cg(x1, x2, p1, p2, a1, a2):
    dx = x1 - x2
    p12 = a1 * p2 + a2 * p1
    d2xsij = -(4j * a1 * a2 * dx * p12 + 2 * a1 * a2 * (a1 + a2) - 4 * dx ** 2 * a1 ** 2 * a2 ** 2 + p12 ** 2) / (
            a1 + a2) ** 2
   # d2xsij *= overlap_CG(x1, x2, p1, p2, a1, a2)

    return d2xsij


def overlap_dx_cg(x1, x2, p1, p2, a1, a2):
    dx = x1 - x2
    p12 = a1 * p2 + a2 * p1
    dxCG = -(p12 * 1j - 2 * a1 * a2 * dx) / (a1 + a2)
    dxCG *= overlap_CG(x1, x2, p1, p2, a1, a2)
    return dxCG


def overlap_trajs(T1, T2,t0):
    geo = g.initgeom()

    w1 = T1.get_widths_time(t0)
    w2 = T2.get_widths_time(t0)
    ph1 = T1.get_time_phases(t0)
    ph2 = T2.get_time_phases(t0)

    ndim = np.size(w1)

    S = np.exp(1j * (-ph1 + ph2), dtype=np.complex128)*np.complex128(1.000+0j)
    i = 0

    for n in range(geo.natoms):
        for j in range(3):
            a1 = w1[i]
            a2 = w2[i]
            x1 = T1.get_pos_time(t0)[i]
            x2 = T2.get_pos_time(t0)[i]
            p1 = T1.get_mom_time(t0)[i]
            p2 = T2.get_mom_time(t0)[i]

            S = S * overlap_CG(x1, x2, p1, p2, a1, a2)
            Sdx = overlap_dx_cg(x1,x2,p1,p2,a1,a2)
            Sd2x= overlap_d2x_cg(x1,x2,p1,p2,a1,a2)
            Sdp= overlap_dp_traj(x1,x2,p1,p2,a1,a2)
            i = i + 1

    return S, Sdx, Sd2x,Sdp


def overlap_dp_traj(x1,x2,p1,p2,a1,a2):


    dx = x1 - x2
    dp = p1 - p2

    ndim = np.size(x1)
    pref = 0.5 * (dp + 2j * a1 * dx) / (a1 + a2)
    pref *= overlap_CG(x1, x2, p1, p2, a1, a2)

    dpSij = pref

    return dpSij


def overlap_v_traj(T1, T2):
    ind1 = T1.trajID
    ind2 = T2.trajID

    if ind1 == ind2:
        V = T1.getpotential_traj()
    else:
        nstates1 = T1.nstates
        for i in range(nstates1):
            T1.HE[i, i] = T1.getpotential_traj_i(i)
            for j in range(i + 1, nstates1):
                T1.HE[i, j] = -1j * coupdotvel(T1, i, j)
                T1.HE[j, i] = -T1.HE[i, j]
        nstates2 = T2.nstates
        for i in range(nstates2):
            T2.HE[i, i] = T2.getpotential_traj_i(i)
            for j in range(i + 1, nstates2):
                T2.HE[i, j] = -1j * coupdotvel(T2, i, j)
                T2.HE[j, i] = -T2.HE[i, j]

        Sdp = overlap_dp_traj(T1, T2)

        F1 = T1.get_traj_force()
        F2 = T2.get_traj_force()

        V = np.dot(T1.getamplitude_traj(), np.matmul(0.5 * (T1.HE + T2.HE), T2.getamplitude_traj())) * overlap_trajs(T1,
                                                                                                                     T2) + \
            0.5 * 1j * np.sum(np.real(Sdp) * (F1 + F2)) * np.dot(T1.getamplitude_traj(), T2.getamplitude_traj())


    return V


def overlap_ke_traj(T1, T2):
    geo = g.initgeom()
    ke =complex(0.0+0j)
    sij = overlap_trajs(T1, T2)
    ndim = T1.ndim
    x1 = T1.getposition_traj()
    x2 = T2.getposition_traj()

    p1 = T1.getmomentum_traj()
    p2 = T2.getmomentum_traj()

    a1 = T1.getwidth_traj()
    a2 = T2.getwidth_traj()

    i = 0
    for n in range(geo.natoms):
        for j in range(3):
            ke -= overlap_d2x_cg(x1[i], x2[i], p1[i], p2[i], a1[n], a2[n]) / (2.0 * T1.allmass[i])
            i = i + 1

    ke = ke * sij
    ke = ke * np.dot(T1.getamplitude_traj(), T2.getamplitude_traj())
    return ke


def overlap_dx_traj(T1, T2,t0):
    geo = g.initgeom()


    ndim = T1.ndim
    x1 = T1.get_pos_time(t0)
    x2 = T2.get_pos_time(t0)

    p1 = T1.get_mom_time(t0)
    p2 = T2.get_mom_time(t0)

    a1 = T1.get_widths_time(t0)
    a2 = T2.get_widths_time(t0)

    dxsij = np.zeros(ndim, dtype=np.complex128)
    i = 0
    for n in range(ndim):
            dxsij[n] = overlap_dx_cg(x1[n], x2[n], p1[n], p2[n], a1[n], a2[n])

    #dxsij = dxsij * sij

    return dxsij


def overlap_sdot_traj(T1, T2):
    Sij = overlap_trajs(T1, T2)
    Sdot = np.dot(T2.getvelocity_traj(), overlap_dx_traj(T1, T2)) + np.dot(T2.get_traj_force(), overlap_dp_traj(T1,
                                                                                                                T2)) + 1j * T2.phasedot() * Sij

    nstates = T2.nstates
    HE = T2.getHE_traj()

    if HE.all() == 0.0:
        HE = T2.get_calc_HE()
        T2.setHE_traj(HE)

    SdotE = -1j * np.dot(T1.stateAmpE, np.matmul(T2.HE, T2.stateAmpE))
    Sdot = Sdot * np.dot(T1.stateAmpE, T2.stateAmpE) + SdotE * Sij

    return Sdot
