from src.geometry import singlepart
from src.dyn_params import initdyn
from src import initialize_traj
from src import abinitio
from src.constants import physconst
from src.propagators import velocityverlet
import numpy as np
import src.Wigner_dist as Wigner_dist
import os
from matplotlib import pyplot as plt
import multiprocessing as mp
from src import swarm
from src import branching as br

from src import buildhs
import src.ekincorrection as ek

print("Number of processors: ", mp.cpu_count())
'''AIMCE propagation'''

''' First, common variables are initialized using the class Traj, based on the dynamics and geometry inputs'''
'''q,p initial values, taken from the wigner dist and d,s,a'''

'''Remove previous files if needed'''
os.system('rm molpro.pun')
os.system('rm molpro_traj*')

''' call initial geometry and dynamic parameters along with pyhisical constants'''

dyn = initdyn()
geo = singlepart()
ph = physconst()

T1 = initialize_traj.trajectory(1, 3, 1)

q = np.zeros(3)
q[0] = 4.00
p = np.zeros_like(q)

T1.setposition_traj(q)
T1.setmomentum_traj(p)

T1.setpotential_traj(np.sum(geo.K * q ** 2 / 2))

T1.setderivs_traj(geo.K * (T1.getposition_traj()))

T1.setmass_traj(geo.mass)
T1.setmassall_traj(np.ones(3) * geo.mass)

amps = np.complex128(1.00 + 0.00j)
T1.setamplitudes_traj(amps)

phases = np.zeros(T1.nstates)  # Phase of the wfn, would be S in the previous equation
T1.setphases_traj(phases)

T1.setwidth_traj(dyn.gamma[3])

dt = np.double(2.0)
print(dt)
time = np.linspace(0, 150, 1500)
amps = np.zeros((1500, 2))
ekin_tr = 0
t = 0
calc1 = True

for i in range(1500):
    #    print('coups1:', np.sqrt(np.sum(T1.getcoupling_traj(0, 0) ** 2)), np.sqrt(np.sum(T1.getcoupling_traj(0, 1) ** 2)))
    print('forces:', T1.get_traj_force())
    t = t + 0.05

    T1 = velocityverlet(T1, dt, i, calc1)

    #   T1,ekin_tr=ek.calc_ekin_tr(T1,ekin_tr)
    #  print(ekin_tr)
    # B = magnus(B, dt)

    print('step ', i)
    print('norm1', np.sum(np.abs(T1.stateAmpE) ** 2))

    print('pop s01: ', np.abs(T1.stateAmpE[0]) ** 2)
    # print('pop s11: ', np.abs(T1.stateAmpE[1]) ** 2)

    # ratios=br.branching_ratios(B)

    energy1 = T1.getpotential_traj() + T1.getkineticlass() - ekin_tr
    print('Energy1: ', energy1)
    # energy2 = B.Traj[1].getpotential_traj() + B.Traj[1].getkineticlass() - ekin_tr
    # print('Energy1: ', energy2)
    plt.scatter(t, np.double(T1.getpotential_traj()-2.8), c='blue')
    plt.scatter(t, np.double(T1.getkineticlass()-2.8), c='red')
    plt.scatter(t, np.double(energy1-2.8), c='black')
    # plt.scatter(t, np.double(T1.getposition_traj()[0]), c='red')
    # plt.scatter(t, np.double(T1.getposition_traj()[1]), c='blue')
    # plt.scatter(t, np.double(T1.getposition_traj()[2]), c='black')
    plt.pause(0.001)
    amps[i, 0] = np.abs(T1.stateAmpE[0]) ** 2
# amps[i, 1] = np.abs(T1.stateAmpE[1]) ** 2
# plt.plot(time, ), c='blue')
# plt.plot(time, np.abs(amps[:, 1]), c='red')

plt.show()
