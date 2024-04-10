from src.geometry import initgeom
import numpy as np
from src.constants import physconst

class initdyn():

    def __init__(self):
        ph = physconst()
        self._ntraj = 50  # Number of trajectories
        self._gamma = np.zeros((14*3))  # Gamma var, width of the gaussians
        self._nstep = 100  # Time-steps by trajectory
        self._dt = (0.1* 1e-15 / (ph.au2sec))
        self._nstates = 2
        self._state = range(0, self._nstates)  # array of states
        self._inipes = 2  # Initial state
        self._e_ref = np.double(-70) #Reference energy for the ground state


    @property
    def ntraj(self):
        return self._ntraj

    @ntraj.setter
    def ntraj(self, value):
        self._ntraj = value

    @property
    def nstep(self):
        return self._nstep

    @nstep.setter
    def nstep(self, value):
        self._nstep = value

    @property
    def dt(self):
        return self._dt

    @dt.setter
    def dt(self, value):
        self._dt = value

    @property
    def inipes(self):
        return self._inipes

    @inipes.setter
    def inipes(self, value):
        self._inipes = value

    @property
    def nstates(self):
        return self._nstates

    @nstates.setter
    def nstates(self, value):
        self._nstates = value

    @property
    def e_ref(self):
        return self._e_ref

    @e_ref.setter
    def e_ref(self, value):
        self._e_ref = value

    @property
    def compress(self):
        return self._compress

    @compress.setter
    def compress(self, value):
        self._compress = value

    @property
    def gamma(self):
        return self._gamma

    @gamma.setter
    def gamma(self, value):
        self._gamma= value

