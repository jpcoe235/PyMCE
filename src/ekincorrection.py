import numpy as np
import src.geometry as geometry
import src.dyn_params as dyn_params
import src.initialize_traj as initialize_traj

''' Not important, makes an energy correction for the COM movement'''
def calc_ekin_tr(T, ekin_tr):
    X = T.getposition_traj()
    X_old = T.getoldpos_traj()
    P = T.getmomentum_traj()
    P_old = T.getoldpos_traj()

    geo = geometry.initgeom()
    dyn = dyn_params.initdyn()
    dr_com = np.zeros(3)
    v_com = np.zeros_like(dr_com)
    k = 0
    totmass = np.sum(geo.masses)
    dt = dyn.dt
    vel = np.zeros_like(P)
    for i in range(geo.natoms):
        for j in range(3):
            dr_com[j] += geo.masses[i] * (X[k] - X_old[k])
            vel[k] = P[k] / geo.masses[i]
            k += 1
    k = 0
    for i in range(geo.natoms):
        for j in range(3):
            v_com[j] += vel[k] * geo.masses[i] / totmass
            k+=1

    idf = 0
    q_corr = np.zeros_like(X)
    p_corr = np.zeros_like(X)
    for i in range(geo.natoms):
        for j in range(3):
    #        q_corr[idf] = X[idf] - dr_com[j]
            vel[idf] = vel[idf] - v_com[j]
            p_corr[idf]=vel[idf]*geo.masses[i]
            idf += 1

    for i in range(3):
        ekin_tr = ekin_tr + 0.5 * (np.sqrt(totmass)* v_com[i]) ** 2.0

    T.setmomentum_traj(p_corr)
   # T.setposition_traj(q_corr)
    return T, ekin_tr
