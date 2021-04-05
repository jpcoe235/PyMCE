import input
import Wigner_dist
import numpy as np
import Traj

'''AIMCE propagation'''

''' First, common variables are initialized using the class Traj, based on the dynamics and geometry inputs'''
'''q,p initial values, taken from the wigner dist'''

geo = input.initgeom()
dyn = input.initdyn()
tr = Traj.traj()

qin, pin = Wigner_dist.WignerSampling()

tr.n=1
tr.
tr.q = qin
tr.p = pin
tr.d = np.zeros(dyn.nstates, dtype=complex)
tr.d[dyn.inipes - 1] = 1.0 + 0j
tr.s=np.zeros(dyn.nstates)
tr.a=tr.d*np.exp(1j*tr.s)
