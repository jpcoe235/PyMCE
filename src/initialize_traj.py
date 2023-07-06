import numpy as np
from Constants import physconst
import overlaps as ov

''' Class defining the initial geometry'''


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
        self.amp = 1.0
        self.d = 1.0
        self.phase = 0
        self.deadtime = 0
        self.phaseE = np.zeros(self.nstates)
        self.stateAmpE = np.zeros(self.nstates)
        self.HE = np.zeros((self.nstates, self.nstates), dtype=np.complex128)
        self.position = np.zeros(ndim)
        self.momentum = np.zeros(ndim)
        self.width = np.zeros(npart)
        self.mass = np.zeros(npart)
        self.allmass = np.zeros(ndim)
        self.transdipole = 0
        self.PotEn = np.zeros(numstates)
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
        self.dotph = 0.0
        self.phasewf = 1
        self.configs = []
        self.file_calc = ''

    def getfilecalc(self):
        return self.file_calc

    def setfilecalc(self, value):
        self.file_calc = value

    def getconfigs(self):
        return self.configs

    def setconfigs(self, value):
        self.configs = value

    def getphasewf(self):
        return self.phasewf

    def setphasewf(self, value):
        self.phasewf = value

    def getcivecs(self):
        return self.civecs

    def setcivecs(self, value):
        self.civecs = value

    def gettrajid_traj(self):
        return self.trajID

    def settrajid_traj(self, value):
        self.trajID = value

    def getpoten_traj(self):
        return self.PotEn

    def getd_traj(self):
        return self.d

    def setd_traj(self, value):
        self.d = value

    def getwidth_traj(self):
        return self.width

    def getposition_traj(self):
        return self.position

    def setoldpos_traj(self, value):
        self.old_pos = value

    def getamp_traj(self):
        return self.amp

    def setamp_traj(self, value):
        self.amp = value

    def getphase_traj(self):
        return self.phase

    def setphase_traj(self, value):
        self.phase = value

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
        if self.nstates > 1:
            return self.stateAmpE[i]
        else:
            return self.stateAmpE

    def getfullforce(self):
        result=np.zeros((self.ndim, self.nstates))
        for i in range(self.nstates):
            result[:,i]=self.derivmat[:,i,i]
        return result
    def getforce_traj(self, i):
        if self.nstates > 1:
            return self.derivmat[:, i, i]
        else:
            return self.derivmat
    def setforce_traj(self,value):
        for i in range(self.nstates):
            self.derivmat[:,i,i]=value[:,i]
        return
    def getcoupling_traj(self, i, j):
        if i == j:
            return np.zeros(self.ndim)
        else:
            return self.derivmat[:, i, j]

    def setderivs_traj(self, value):
        self.derivmat = value

    def getpotential_traj_i(self, Istate):
        if self.nstates > 1:
            return self.PotEn[Istate]
        else:
            return self.PotEn

    def getpotential_traj(self):
        E = np.double(0.0)

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

    def getphases_traj(self):
        return self.phaseE

    def setwidth_traj(self, value):
        self.width = value

    def getvelocity_traj(self):
        V = self.getmomentum_traj() / self.getmassall_traj()
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

    def get_calc_HE_traj(self):
        nstates = self.nstates
        HE = np.zeros((nstates, nstates))
        for i in range(nstates):
            HE[i, i] = self.getpotential_traj_i(i)
            for j in range(i + 1, nstates):
                HE[i, j] = -1j * ov.coupdotvel(self, i, j)
                HE[j, i] = -HE[i, j]
        self.setHE_traj(HE)
        return self.HE

    def getHE_traj(self):
        return self.HE

    def setHE_traj(self, value):
        self.HE = value

    def phasedot(self):
        self.dotph = self.getkineticlass() - 0.50000 * np.sum(self.getwidth_traj() / self.getmassall_traj())
        return np.real(self.dotph)

    def getkineticlass(self):

        p = self.getmomentum_traj()
        energy = 0.5 * np.dot(p, p / self.getmassall_traj())

        return energy

    def compforce(self, A, F, E, C):
        nst = self.nstates
        f1 = np.zeros(self.ndim)
        f2 = np.zeros(self.ndim)

        for i in range(nst):
            f1 = f1 + F[:, i] * abs(A[i]) ** 2
        for i in range(nst):
            for j in range(i + 1, nst):
                tmp = np.cdouble(2.0) * np.real(np.conj(A[i]) * A[j]) * (E[i] - E[j])
                f2 = f2 + tmp * C[:, i, j]
        fvec = f1 + f2
        return fvec
