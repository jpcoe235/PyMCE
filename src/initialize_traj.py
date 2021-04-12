import numpy as np
from Constants import physconst

''' Class defining the initial geometry mass weighted'''


class trajectory():
    def __init__(self, npart, ndim, numstates):
        self.trajID = 0
        self.centID = -1
        self.zcent = False
        self.nstates = numstates
        self.ndim = ndim * npart
        self.npart = npart
        self.stateID = 0
        self.currentime = 0.0
        self.amp = 0
        self.phase = 0
        self.deadtime = 0
        self.phaseE = np.zeros(self.nstates)
        self.stateAmpE = np.zeros(self.nstates)
        self.HE = np.zeros((self.nstates, self.nstates))
        self.position = np.zeros(ndim)
        self.momentum = np.zeros(ndim)
        self.width = np.zeros(npart)
        self.mass = np.zeros(npart)
        self.allmass = np.zeros(ndim)
        self.transdipole = 0
        self.PotEn = 0
        self.derivmat = 0
        self.dEdx_GA = 0
        self.dipole = 0
        self.modforce = 0
        self.civecs = 0
        self.oldorbs = 0
        self.elecphase = 0
        self.den = 0
        self.SpawnMode = 0
        self.CoupHist = 0.0
        self.SpawnTime = 0
        self.LastSpawn = 0
        self.TunnelSpawnTime = 0.0
        self.SpawnCoupled = 0
        self.old_pos = 0.0
        self.old_mom = 0.0
        self.old_amp = 0.0

    def getwidth_traj(self):
        return self.width

    def getposition_traj(self):
        return self.position

    def setoldpos_traj(self, value):
        self.old_pos = value

    def setoldmom_traj(self, value):
        self.old_mom = value

    def setoldamp_traj(self, value):
        self.old_amp = value

    def getoldamp_traj(self):
        return self.old_amp

    def getoldpos_traj(self):
        return self.old_pos

    def getoldmom_traj(self):
        return self.old_mom

    def getmomentum_traj(self):
        return self.momentum

    def setposition_traj(self, pos):
        self.position = pos

    def setmomentum_traj(self, mom):
        self.momentum = mom

    def setamplitudes_traj(self, value):
        self.stateAmpE = value

    def getamplitude_traj(self):
        return self.stateAmpE

    def getamplitude_traj_i(self, i):
        return self.stateAmpE[i]

    def getforce_traj(self, i):
        return self.derivmat[:, i, i]

    def getcoupling_traj(self, i, j):
        return self.derivmat[:, i, j]

    def setderivs_traj(self, value):
        self.derivmat = value

    def getpotential_traj_i(self, Istate):
        return self.PotEn[Istate]

    def getpotential_traj(self):
        E = 0.

        for i in range(self.nstates):
            E += self.getpotential_traj_i(i) * np.abs(self.stateAmpE[i]) ** 2

        return E

    def setpotential_traj(self, value):
        self.PotEn = value

    def setmass_traj(self, value):
        self.mass = value

    def setmassall_traj(self, value):
        self.allmass = value

    def getmassall_traj(self):
        return self.allmass

    def getmass_traj(self):
        return self.mass

    def setphases_traj(self, value):
        self.phaseE = value

    def setwidth_traj(self, value):
        self.width = value

    def getvelocity_traj(self):
        V = np.zeros(self.ndim)
        k = 0
        for i in range(self.npart):
            for j in range(3):
                V[k] = self.momentum[k] / self.mass[i]
                k += 1
        return V

    def get_traj_force(self):
        nst = self.nstates
        f1 = np.zeros(self.ndim)
        f2 = np.zeros(self.ndim)

        E = np.zeros(nst)
        a = np.zeros(nst, dtype=np.complex128)
        for i in range(nst):
            E[i] = self.getpotential_traj_i(i)
            a[i] = self.getamplitude_traj_i(i)

        for i in range(nst):
            f1 += self.getforce_traj(i) * np.abs(a[i]) ** 2
        for i in range(nst):
            for j in range(i + 1, nst):
                tmp = 2.0 * np.real(np.conj(a[i]) * a[j]) * (E[i] - E[j])
                f2 += tmp * self.getcoupling_traj(i, j)

        fvec = f1 + f2

        return fvec


def compforce(T, A, F, E, C):
    nst = T.nstates
    f1 = np.zeros(T.ndim)
    f2 = np.zeros(T.ndim)

    for i in range(nst):
        f1 = f1 + F[:, i] * np.abs(A[i]) ** 2
    for i in range(nst):
        for j in range(nst):
            tmp = 2.0 * np.real(np.conj(A[i]) * A[j]) * (E[i] - E[j])
            f2 += tmp * C[:, i, j]
    fvec = f1 + f2
    return fvec


def phasedot(T):
    dotph = getkineticlass(T) - 0.5 * np.sum(T.getwidth_traj() / T.getmass_traj())
    return dotph


def getkineticlass(T):
    k = 0
    energy = 0
    p = T.getmomentum_traj()

    for i in range(T.npart):
        for j in range(3):
            energy += 0.5 * p[k] * p[k] / T.mass[i]
            k = k + 1

    return energy
