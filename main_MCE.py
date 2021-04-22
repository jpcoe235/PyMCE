from src.geometry import initgeom
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
import src.ekincorrection as ek

print("Number of processors: ", mp.cpu_count())
'''AIMCE propagation'''

''' First, common variables are initialized using the class Traj, based on the dynamics and geometry inputs'''
'''q,p initial values, taken from the wigner dist and d,s,a'''

'''Remove previous files if needed'''
os.system('rm molpro.pun')
os.system('rm molpro_traj*')

''' call initial geometry and dynamic parameters along with pyhisical constants'''
# As I am doing everything mass-weighted, also applied to the widths
dyn = initdyn()
geo = initgeom()
ph = physconst()

'''First initialize and populate one trajectory'''

T1 = initialize_traj.trajectory(geo.natoms, 3, dyn.nstates)
qin, pin = Wigner_dist.WignerSampling()
q = np.zeros_like(qin)
p = np.zeros_like(pin)
with open('initialqp.dat', 'r') as f:
    f.readline()
    for i in range(geo.ndf):
        N, M = f.readline().strip().split()
        q[i] = np.double(float(N.replace('D', 'E'))) / np.sqrt(geo.massrk[i])
        p[i] = np.double(float(M.replace('D', 'E'))) * np.sqrt(geo.massrk[i])

T1.setposition_traj(q + geo.rkinit)
T1.setmomentum_traj(p)

pec, der = abinitio.inp_out(0, 0, geo, T1)  # First ab-initio run

T1.setpotential_traj(pec)  # taking V(R) from ab-initio
T1.setderivs_traj(
    der)  # derivatives matrix not mass-weighted (possibly change that), diagonals are forces and off-d are nacmes
T1.setmass_traj(geo.masses)  # mass of every atom in a.u (the dimmension is natoms/nparts)
T1.setmassall_traj(
    geo.massrk)  # mass in every degree of freedom (careful to use it, it can triple the division/multiplication easily)

amps = np.zeros(T1.nstates, dtype=np.complex128)
amps[dyn.inipes - 1] = np.complex128(
    1.00 + 0.00j)  # Amplitudes of Ehrenfest trajectories, they should be defined as a=d *exp(im*S)

T1.setamplitudes_traj(amps)

phases = np.zeros(T1.nstates)  # Phase of the wfn, would be S in the previous equation

T1.setphases_traj(phases)
T1.setwidth_traj(dyn._gamma)

dt = np.double(dyn.dt)
print(dt)
time = np.linspace(0, 150, 1500)
amps = np.zeros((1500, 2))
ekin_tr = 0
t = 0
for i in range(1500):
    t = t + 0.05
    T1 = velocityverlet(T1, dt, i)

    #   T1,ekin_tr=ek.calc_ekin_tr(T1,ekin_tr)
    #  print(ekin_tr)

    print('step ', i)
    print('norm', np.sum(np.abs(T1.stateAmpE) ** 2))
    print('pop s0: ', np.abs(T1.stateAmpE[0]) ** 2)
    print('pop s1: ', np.abs(T1.stateAmpE[1]) ** 2)
    print('Amplitude 0: ', T1.stateAmpE[0])
    print('Amplitude 1: ', T1.stateAmpE[1])
    energy = T1.getpotential_traj() + T1.getkineticlass() - ekin_tr

    print('potential energy:', T1.getpotential_traj_i(0), T1.getpotential_traj_i(1))
    print('kinetic energy:', T1.getkineticlass())
    print('Energy: ', energy)

    print('nacs', np.sum(T1.getcoupling_traj(0, 1)))

    # plt.scatter(t, np.double(energy), c='blue')
    # plt.pause(0.1)
    amps[i, 0] = np.abs(T1.stateAmpE[0]) ** 2
    amps[i, 1] = np.abs(T1.stateAmpE[1]) ** 2
plt.plot(time, np.abs(amps[:, 0]), c='blue')
plt.plot(time, np.abs(amps[:, 1]), c='red')

plt.show()
