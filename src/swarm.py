""" Calculate the initial set of trajectories, maybe better starting with a Wigner distribution"""

from src.geometry import initgeom
from src.dyn_params import initdyn
from src import initialize_traj
from src import abinitio
from src.constants import physconst
from src.propagators import velocityverlet
import numpy as np
import Wigner_dist
from src import bundle
from src import branching
from src import buildhs
from src import overlaps as ov
from src import geometry
import cmath

def createswarm(ntraj, npart, ndim, numstates):
    geo = geometry.initgeom()
    trajs = []
    Tinit = inittraj()
    Tinit.settrajid_traj(0)
    # trajs.append(Tinit)
    print(Tinit.PotEn)
    X = Tinit.getposition_traj()
    print(X)
    P = Tinit.getmomentum_traj()
    sigma = 0.1
    p_norm = np.sqrt(np.sum(P ** 2))
    for i in range(ntraj):
        T = initialize_traj.trajectory(npart=npart, ndim=ndim, numstates=numstates)
        n1 = np.random.rand()
        dx = sigma * (2 * n1 - 1.0)
        n2 = np.random.rand()
        dp = sigma / p_norm * (2 * n2 - 1)

        T.setposition_traj(X + dx)

        T.setmomentum_traj(P + dp)
        T.setamplitudes_traj(Tinit.getamplitude_traj())
        T.setphases_traj(Tinit.getphase_traj())
        T.setwidth_traj(Tinit.getwidth_traj())
        print('widths:',T.getwidth_traj())
        T.settrajid_traj(i)
        pec, der = abinitio.inp_out(0, 0, geo, T)  # First ab-initio run
        T.setpotential_traj(pec)  # taking V(R) from ab-initio

        T.setderivs_traj(
            der)  # derivatives matrix mass-weighted (possibly change that), diagonals are forces and off-d are nacmes
        T.setmass_traj(geo.masses)  # mass of every atom in a.u (the dimmension is natoms/nparts)
        T.setmassall_traj(
            geo.massrk)  # mass in every degree of freedom (careful to use it, it can triple the division/multiplication easily)
        trajs.append(T)

    B = bundle.bundle(ntraj, npart, ndim, numstates)
    B.setTraj_bundle(trajs)
    norm = B.get_calc_set_norm()

    for i in range(ntraj):
        B.Traj[i].amp = B.Traj[i].amp / cmath.sqrt(norm)

    B.setamps_bundle(B.getamps_bundle()/cmath.sqrt(norm))

    print('Amplitudes before normalizing: ', B.getamps_bundle())
    B = buildhs.buildsandh(B)

    Sif = np.zeros(ntraj, dtype=np.complex128)
    for i in range(ntraj):
        Sif[i] = ov.overlap_trajs(Tinit, B.Traj[i])

    print('Overlaps: ', Sif)
    print('Sinv :', B.Sinv)

    Sif_2 = np.matmul(np.conj(B.Sinv), Sif)

    print('Overlaps*Sinv: ', Sif_2)
    B.setamps_bundle(Sif_2)

    norm = B.get_calc_set_norm()

    print('norm: ', norm)

    B.setamps_bundle(B.getamps_bundle() / np.sqrt(norm))

    print('Final bundle amplitudes: ', B.getamps_bundle())
    print('Norm Final bundle: ', np.linalg.norm(B.getamps_bundle()))
    print('Norm Final bundle branch: ', B.get_calc_set_norm())
    return B


def inittraj():
    dyn = initdyn()
    geo = initgeom()
    ph = physconst()

    '''First initialize and populate one trajectory'''

    T1 = initialize_traj.trajectory(geo.natoms, 3, dyn.nstates)
    # qin, pin = Wigner_dist.WignerSampling()

    q = np.zeros(geo.ndf)

    p = np.zeros_like(q)
    with open('initialqp.dat', 'r') as f:
        f.readline()
        for i in range(geo.ndf):
            N, M = f.readline().strip().split()
            q[i] = np.double(float(N.replace('D', 'E'))) / np.sqrt(geo.massrk[i])
            p[i] = np.double(float(M.replace('D', 'E'))) * np.sqrt(geo.massrk[i])

    T1.setposition_traj(q + geo.rkinit)
    T1.setmomentum_traj(p)

    pec, der = abinitio.inp_out(0, 0, geo, T1)  # First ab-initio run
    print(pec)
    T1.setpotential_traj(pec)  # taking V(R) from ab-initio
    T1.setderivs_traj(
        der)  # derivatives matrix mass-weighted (possibly change that), diagonals are forces and off-d are nacmes
    T1.setmass_traj(geo.masses)  # mass of every atom in a.u (the dimmension is natoms/nparts)
    T1.setmassall_traj(
        geo.massrk)  # mass in every degree of freedom (careful to use it, it can triple the division/multiplication easily)

    amps = np.zeros(T1.nstates, dtype=np.complex128)
    amps[dyn.inipes - 1] = 1.00  # Amplitudes of Ehrenfest trajectories, they should be defined as a=d *exp(im*S)

    T1.setamplitudes_traj(amps)

    phases = np.zeros(T1.nstates)  # Phase of the wfn, would be S in the previous equation

    T1.setphases_traj(phases)
    T1.setwidth_traj(dyn._gamma)

    return T1


def poptraj(T1):
    dyn = initdyn()
    geo = initgeom()
    ph = physconst()

    pec, der = abinitio.inp_out(0, 0, geo, T1)  # First ab-initio run

    T1.setpotential_traj(pec)  # taking V(R) from ab-initio
    T1.setderivs_traj(
        der)  # derivatives matrix mass-weighted (possibly change that), diagonals are forces and off-d are nacmes
    T1.setmass_traj(geo.masses)  # mass of every atom in a.u (the dimmension is natoms/nparts)
    T1.setmassall_traj(
        geo.massrk)  # mass in every degree of freedom (careful to use it, it can triple the division/multiplication easily)

    amps = np.zeros(T1.nstates, dtype=np.complex128)
    amps[dyn.inipes - 1] = 1.00  # Amplitudes of Ehrenfest trajectories, they should be defined as a=d *exp(im*S)

    T1.setamplitudes_traj(amps)

    phases = np.zeros(T1.nstates)  # Phase of the wfn, would be S in the previous equation

    T1.setphases_traj(phases)
    T1.setwidth_traj(dyn.gamma)

    return T1
