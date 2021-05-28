from src.geometry import initgeom
from src.dyn_params import initdyn
from src import initialize_traj
from src import abinitio
from src.constants import physconst
from src.propagators import velocityverlet
import src.overlaps as ovs
import numpy as np
import src.Wigner_dist as Wigner_dist
import os
from matplotlib import pyplot as plt
import multiprocessing as mp
from src.outputs import output_traj as wrtout
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
os.system('rm *.wfu')
os.system('rm /home/AndresMoreno/wfu/*')

''' call initial geometry and dynamic parameters along with pyhisical constants'''

dyn = initdyn()
geo = initgeom()
ph = physconst()

'''First initialize and populate one trajectory'''

T1 = initialize_traj.trajectory(geo.natoms, 3, dyn.nstates)
# qin, pin = Wigner_dist.WignerSampling()
q = np.zeros(geo.ndf)
p = np.zeros_like(q)
with open('initialqp_2.dat', 'r') as f:
    f.readline()
    for i in range(geo.ndf):
        N, M = f.readline().strip().split()
        q[i] = np.double(float(N.replace('D', 'E')))
        p[i] = np.double(float(M.replace('D', 'E')))

T1.setposition_traj(q)
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
T1.setwidth_traj(dyn.gamma)
# B = swarm.createswarm(2, geo.natoms, 3, dyn.nstates)
# B = buildhs.buildsandh(B)

dt = 2.5000000
print(dt)
time = np.linspace(0, 150, 15000)
amps = np.zeros((15000, 2))
ekin_tr = 0
t = 0
calc1 = False
wrtout(True,T1,0.000)
for i in range(15000):

    # for n1 in range(T1.nstates):
    #     for n2 in range(n1,T1.nstates):
    #         print('coupsdotvel:', ovs.coupdotvel(T1, n1, n2))
    #         if np.abs(ovs.coupdotvel(T1, n1, n2)) > 0.005:
    #             dt = dt / 4
    #             print('going to smaller time-steps')


    if i==12:
        dt=dt/2.0000
    if i==62:
        dt=dt/2.00
    if i==185:
        dt=dt/5.00
    print('step ', i)
    print('time:', t)

    print('coups1:', np.sqrt(np.sum(T1.getcoupling_traj(0, 0) ** 2)), np.sqrt(np.sum(T1.getcoupling_traj(0, 1) ** 2)),
          np.sqrt(np.sum(ovs.coupdotvel(T1, 0, 1) ** 2)))

    print('forces:', T1.get_traj_force())

    #   T1,ekin_tr=ek.calc_ekin_tr(T1,ekin_tr)
    #  print(ekin_tr)
    # B = magnus(B, dt)

    print('norm1', np.sum(np.abs(T1.stateAmpE) ** 2))

    print('pop s01: ', np.abs(T1.stateAmpE[0]) ** 2)
    print('pop_branching', np.real(np.dot(T1.stateAmpE[0], T1.stateAmpE[0])))
    print('pop s11: ', np.abs(T1.stateAmpE[1]) ** 2)
    print('traj momentum:', T1.getmomentum_traj()[0])
    print('traj position:', T1.getposition_traj()[0])

    # ratios=br.branching_ratios(B)

    energy1 = T1.getpotential_traj() + T1.getkineticlass() - ekin_tr
    print('Energy1: ', energy1)
    print('Epotential:', T1.getpotential_traj())
    print('Ekinetic:', T1.getkineticlass())
    print('Ekinetic_ov:', np.abs(ovs.overlap_ke_traj(T1, T1)))
    print('overlap:', ovs.overlap_trajs(T1, T1))
    print('phase:',T1.phase)
    print('constant_term_phase:,',0.500* np.sum(T1.getwidth_traj()/T1.getmass_traj()))
    # energy2 = B.Traj[1].getpotential_traj() + B.Traj[1].getkineticlass() - ekin_tr
    # print('Energy1: ', energy2)
    # plt.scatter(t, np.double(energy), c='blue')
    # plt.pause(0.1)

    amps[i, 0] = np.abs(T1.stateAmpE[0]) ** 2
    amps[i, 1] = np.abs(T1.stateAmpE[1]) ** 2
    Told = T1
    T2 = velocityverlet(Told, dt, i, calc1)
    energy2 = T2.getpotential_traj() + T2.getkineticlass() - ekin_tr
    print('Endif:', np.abs(energy1 - energy2))
    # if np.abs(energy1 - energy2) > 9E-7:
    #     os.system('cp /home/AndresMoreno/wfu/' + str(i - 1) + '.check.wfu 003.molpro.wfu')
    #     dt_2 = dt / 2.0000
    #     for j in range(2):
    #         Told = velocityverlet(Told, dt_2,j, calc1)
    #     print('Adaptative timestep')
    #     T1=Told
    #     energy2 = T1.getpotential_traj() + T1.getkineticlass() - ekin_tr
    #     print('Adapif:', np.abs(energy1 - energy2))
    # else:
    T1 = T2

    t = t + dt
    wrtout(False, T1,t)
plt.plot(time, np.abs(amps[:, 0]), c='blue')
plt.plot(time, np.abs(amps[:, 1]), c='red')

plt.show()
