import initialize_traj
import numpy as np
from src import overlaps as ov
from src import branching


class bundle():

    def __init__(self, ntraj, numstates):

        self.ntraj = ntraj #Number of trajectories
        self.nstates = numstates # numberofstates
        self.Traj = [] #Dict for the trajectories
        self.pops = np.zeros(numstates) #population of each state
        self.S = np.zeros((ntraj, ntraj), dtype=np.complex128) #Overlap matrices
        self.H = np.zeros((ntraj, ntraj), dtype=np.complex128) #Hamiltonian
        self.Sdot = np.zeros_like(self.S, dtype=np.complex128) #Derivative of the S matrix
        self.Sinv = np.zeros_like(self.S, dtype=np.complex128) #inverse of S
        self.T = np.zeros_like(self.H, dtype=np.complex128) # Kinetic energy of the bundle
        self.V = np.zeros_like(self.H, dtype=np.complex128) # Potential energy of the bundle
        self.SE = np.zeros_like(self.H, dtype=np.complex128) #Ehrenfest overlap
        self.Heff = np.zeros((ntraj, ntraj), dtype=np.complex128) # Effective hamiltonian (H-iS)
        self.Heff1 = np.zeros((ntraj, ntraj), dtype=np.complex128) #Effective hamiltonian 1 (matmul Sinv,H-iS)
        self.time = 0.0 #Time
        self.norm = 1.0 #Nrom
        self.amps = np.ones(ntraj) # Amplitudes C

    def getTraj_bundle(self):
        return self.Traj

    def getS_bundle(self):
        return self.S

    def getpops_bundle(self):
        return self.pops

    def getSinv_bundle(self):
        return self.Sinv

    def getH_bundle(self):
        return self.H

    def getT_bundle(self):
        return self.T

    def getV_bundle(self):
        return self.V

    def getSE_bundle(self):
        return self.SE

    def getHeff_bundle(self):
        return self.Heff

    def getHeff1_bundle(self):
        return self.Heff1

    def gettime_bundle(self):
        return self.time

    def getnorm_bundle(self):
        return self.norm

    def getamps_bundle(self):
        return self.amps

    def setamps_bundle(self, value):
        self.amps = value

    def setTraj_bundle(self, value):
        self.Traj = value

    def setS_bundle(self, value):
        self.S = value

    def setpops_bundle(self, value):
        self.pops = value

    def setSinv_bundle(self, value):
        self.Sinv = value

    def setH_bundle(self, value):
        self.H = value

    def setT_bundle(self, value):
        self.T = value

    def setV_bundle(self, value):
        self.V = value

    def setSE_bundle(self, value):
        self.SE = value

    def setHeff_bundle(self, value):
        self.Heff = value

    def setHeff1_bundle(self, value):
        self.Heff1 = value

    def settime_bundle(self, value):
        self.time = value

    def setnorm_bundle(self,value):
        self.norm = value

    def addtraj(self, newtraj):
        '''All set to 0, call the appropiate routines to recalculate the couplings and Hamilt'''
        self.ntraj = self.ntraj + 1
        ntraj = self.ntraj
        self.Traj[self.ntraj - 1] = newtraj
        self.pops = np.zeros(self.nstates)
        self.S = np.zeros((ntraj, ntraj))
        self.Sinv = np.zeros((ntraj, ntraj))
        self.H = np.zeros((ntraj, ntraj))
        self.Heff = np.zeros((ntraj, ntraj))
        self.Heff1 = np.zeros((ntraj, ntraj))

    def killtraj(self, ind):
        self.ntraj = self.ntraj - 1
        ntraj = self.ntraj
        self.T[ind] = []
        self.pops = np.zeros(self.nstates)
        self.S = np.zeros((ntraj, ntraj))
        self.Sinv = np.zeros((ntraj, ntraj))
        self.H = np.zeros((ntraj, ntraj))
        self.Heff = np.zeros((ntraj, ntraj))
        self.Heff1 = np.zeros((ntraj, ntraj))

    def get_calc_set_overlaps(self):
        for i in range(self.ntraj):
            for j in range(self.ntraj):
                self.S[i, j] = ov.overlap_trajs(self.Traj[i], self.Traj[j])

        return self.S

    def get_calc_set_norm(self):
        ratio = branching.branching_ratios(self)
        print('ratios :',ratio)
        self.norm = sum(ratio)

        return self.norm
