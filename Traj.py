import input
import Constants
import Wigner_dist
import numpy as np


class traj:
    def __init__(self):
        geo = input.initgeom()
        dyn = input.initdyn()
        ph = Constants.physconst()
        self._n = 0
        self._totN = dyn.ntraj
        self._ind = 0
        self._q = np.zeros(geo.ndf)
        self._p = np.zeros(geo.ndf)
        self._dc = 0
        self._time = 0
        self._dt = (dyn.dt / (ph.au2sec/ 1e-15))
        self._nprop=1500
        self._nstep=1
        self._iprop = 1
        self._istep = 0
        self._isubstep = 0
        self._a=0
        self._d=np.zeros(dyn.nstates)
        self._s=np.zeros(dyn.nstates,dtype=complex)
        self._dqdt=0
        self._dpdt=0
        self._dsdt=0
        self._dddt=0
        self._epot=0
        self._grad=0
        self._nac=0
        self._ekin=0


    @property
    def dddt(self):
        return self._dddt
    @dddt.setter
    def dddt(self,value):
        self._dddt=value

    @property
    def nstep(self):
        return self._nstep

    @nstep.setter
    def nstep(self, value):
        self._nstep = value



    @property
    def n(self):
        return self._n

    @n.setter
    def n(self, value):
        self._n = value

    @property
    def q(self):
        return self._q

    @q.setter
    def q(self, value):
        self._q = value

    @property
    def p(self):
        return self._p

    @p.setter
    def p(self, value):
        self._p = value

    @property
    def ind(self):
        return self._ind

    @ind.setter
    def ind(self, value):
        self._ind = value

    @property
    def time(self):
        return self._time

    @time.setter
    def time(self, value):
        self._time = value

    @property
    def iprop(self):
        return self._iprop

    @iprop.setter
    def iprop(self, value):
        self._iprop = value

    @property
    def istep(self):
        return self._istep

    @istep.setter
    def istep(self, value):
        self._iprop = value

    @property
    def isubstep(self):
        return self._isubstep

    @isubstep.setter
    def isubstep(self, value):
        self._isubstep = value


    @property
    def a(self):
        return self._a

    @a.setter
    def a(self, value):
        self._a = value


    @property
    def s(self):
        return self._s

    @s.setter
    def s(self, value):
        self._s = value

    @property
    def d(self):
        return self._d

    @d.setter
    def d(self, value):
        self._d = value

    @property
    def dqdt(self):
        return self._dqdt

    @dqdt.setter
    def dqdt(self, value):
        self._dqdt = value

    @property
    def dpdt(self):
        return self._dpdt

    @dpdt.setter
    def dpdt(self, value):
        self._dpdt = value

    @property
    def dsdt(self):
        return self._dsdt

    @dsdt.setter
    def dsdt(self, value):
        self._dsdt = value

    @property
    def epot(self):
        return self._epot

    @epot.setter
    def epot(self, value):
        self._epot = value

    @property
    def grad(self):
        return self._grad

    @grad.setter
    def grad(self, value):
        self._grad = value

    @property
    def nac(self):
        return self._nac

    @nac.setter
    def nac(self, value):
        self._nac = value

    @property
    def ekin(self):
        return self._ekin

    @ekin.setter
    def ekin(self, value):
        self._ekin = value

    @property
    def epot(self):
        return self._epot
    @epot.setter
    def epot(self,value):
        self._epot=value

