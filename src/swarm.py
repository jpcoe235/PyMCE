""" Calculate the initial set of trajectories, maybe better starting with a Wigner distribution"""

from src.geometry import initgeom
from src.dyn_params import initdyn
from src import initialize_traj
from src import abinitio
from src.constants import physconst
from src.propagators import velocityverlet
import numpy as np
import Wigner_dist


def createswarm(ntraj, npart, ndim, numstates):
    Tinit = inittraj()

    X = Tinit.getposition_traj()
    P = Tinit.getmomentum_traj()
    sigma = 0.1
    p_norm=np.sqrt(np.sum(P**2)
    for i in range(1, ntraj):
        T = initialize_traj.trajectory(npart, ndim, numstates)
        n1 = np.random.rand(1)
        dx = sigma * (2 * n1 - 1.0)
        n2=np.random.rand()
        dp=sigma/p_norm*(2*n2-1)
        T.setposition_traj(X+dx)
        T.setmomentum_traj(P+dp)


def inittraj():
    dyn = initdyn()
    geo = initgeom()
    ph = physconst()

    '''First initialize and populate one trajectory'''

    T1 = initialize_traj.trajectory(geo.natoms, 3, dyn.nstates)
    qin, pin = Wigner_dist.WignerSampling()

    T1.setposition_traj(qin + geo.rkinit)
    T1.setmomentum_traj(pin)

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
    T1.setwidth_traj(dyn._gamma)

    return T1
