from src.geometry import initgeom
import numpy as np
from src.constants import physconst

class initdyn():

    def __init__(self):
        ph = physconst()
        self._ntraj = 50  # Number of trajectories
        self._ndiffbasis = 1  # Number of basis sets
        self._gamma = np.asarray([22.7,22.7,4.7,4.7,4.7,4.7])  # Gamma var, width of the gaussians
        self._nstep = 100  # Time-steps by trajectory
        self._dt = (0.1* 1e-15 / (ph.au2sec)) # femtoseconds
        self._nstates = 2
        self._state = range(0, self._nstates)  # array of states, Kenichiro made the eqs. up to 3
        self._inipes = 2  # Initial state
        self._e_ref = np.double(-7.80549751000000e+1)

    @property
    def ndiffbasis(self):
        return self._ndiffbasis

    @ndiffbasis.setter
    def ndiffbasis(self, value):
        self._ndiffbasis = value

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

