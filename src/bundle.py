import initialize_traj
import numpy as np
import overlaps as ov
import branching

class bundle():

    def __init__(self, ntraj, npart, ndim, numstates):

        self.ntraj = ntraj
        self.nstates = numstates
        self.T = swarmsampling(ntrajs, npart, ndim, numstates)
        self.pops = np.zeros(nstates)
        self.S = np.zeros((ntraj, ntraj))
        self.Sinv = np.zeros((ntraj, ntraj))
        self.H = np.zeros((ntraj, ntraj))
        self.Heff = np.zeros((ntraj, ntraj))
        self.Heff1 = np.zeros((ntraj, ntraj))
        self.time = 0.0
        self.norm=1.0

    def addtraj(self, newtraj):
        '''All set to 0, call the appropiate routines to recalculate the couplings and Hamilt'''
        self.ntraj = self.ntraj + 1
        ntraj = self.ntraj
        self.T[self.ntrajs - 1] = newtraj
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
        ratio=branching.branching_ratios(self)

        self.norm=sum(ratio)

        return self.norm

