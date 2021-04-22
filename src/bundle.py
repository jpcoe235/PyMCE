import initialize_traj
import numpy as np
import overlaps as ov
import branching


class bundle():

    def __init__(self, ntraj, npart, ndim, numstates):

        self.ntraj = ntraj
        self.nstates = numstates
        self.Traj = []
        self.pops = np.zeros(numstates)
        self.S = np.zeros((ntraj, ntraj))
        self.H = np.zeros((ntraj, ntraj))
        self.Sdot = np.zeros_like(self.S)
        self.Sinv = np.zeros_line(self.S)
        self.T = np.zeros_line(self.H)
        self.V = np.zeros_like(self.H)
        self.SE = 0.0
        self.Heff = np.zeros((ntraj, ntraj))
        self.Heff1 = np.zeros((ntraj, ntraj))
        self.time = 0.0
        self.norm = 1.0
        self.amps=np.ones(ntraj)

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

    def setamps_bundle(self,value):
        self.amps=value

    def setTraj_bundle(self,value):
        self.Traj=value

    def setS_bundle(self,value):
        self.S=value

    def setpops_bundle(self,value):
        self.pops=value

    def setSinv_bundle(self,value):
        self.Sinv=value

    def setH_bundle(self,value):
        self.H=value

    def setT_bundle(self,value):
        self.T=value

    def setV_bundle(self,value):
        self.V=value

    def setSE_bundle(self,value):
        self.SE=value

    def setHeff_bundle(self,value):
        self.Heff=value

    def setHeff1_bundle(self,value):
        self.Heff1=value

    def settime_bundle(self,value):
        self.time=value

    def setnorm_bundle(self):
        self.norm=value


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
                self.S[i, j] = ov.overlap_trajs(self.T[i], self.T[j])

        return self.S

    def get_calc_set_norm(self):
        ratio = branching.branching_ratios(self)

        self.norm = sum(ratio)

        return self.norm
