from src.geometry import initgeom
from src.dyn_params import initdyn
from src import initialize_traj
from src import abinitio
from src.constants import physconst
from src.propagators import velocityverlet
from src.propagators import rk_method_4
from src.propagators import rkf45
import src.overlaps as ovs
import numpy as np
from copy import copy
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

pec, der,cis = abinitio.inp_out(0, 0, geo, T1)  # First ab-initio run

T1.setpotential_traj(pec)  # taking V(R) from ab-initio
T1.setderivs_traj(
    der)  # derivatives matrix not mass-weighted (possibly change that), diagonals are forces and off-d are nacmes
T1.setmass_traj(geo.masses)  # mass of every atom in a.u (the dimmension is natoms/nparts)
T1.setmassall_traj(
    geo.massrk)  # mass in every degree of freedom (careful to use it, it can triple the division/multiplication easily)
T1.setcivecs(cis)
amps = np.zeros(T1.nstates, dtype=np.complex128)
amps[dyn.inipes - 1] = np.complex128(
    1.00 + 0.00j)  # Amplitudes of Ehrenfest trajectories, they should be defined as a=d *exp(im*S)

T1.setamplitudes_traj(amps)

T1.setd_traj(amps)

phases = np.zeros(T1.nstates)  # Phase of the wfn, would be S in the previous equation

T1.setphases_traj(phases)

T1.setwidth_traj(dyn.gamma)
# B = swarm.createswarm(2, geo.natoms, 3, dyn.nstates)
# B = buildhs.buildsandh(B)

dt = 2.5
print(dt)
time = np.linspace(0, 150, 1100)
amps = np.zeros((15000, 2))
ekin_tr = 0
t = 0
calc1 = False
wrtout(True, T1, 0.000)
repeat=True
T0=copy(T1)
for i in range(1100):
    if i==0:
        T1=copy(T0)

    # for n1 in range(T1.nstates):
    #     for n2 in range(n1,T1.nstates):
    #         print('coupsdotvel:', ovs.coupdotvel(T1, n1, n2))
    #         if np.abs(ovs.coupdotvel(T1, n1, n2)) > 0.005:
    #             dt = dt / 4
    #             print('going to smaller time-steps')

    # if i == 10:
    #     dt = dt / 4.0000
    # if i == 62:
    #     dt = dt / 2.00
    # if i == 185:
    #     dt = dt / 5.00
    print('step ', i)
    print('time:', t)

    print('coups1:', np.sqrt(np.sum(T1.getcoupling_traj(0, 0) ** 2)), np.sqrt(np.sum(T1.getcoupling_traj(0, 1) ** 2)),
          (np.sum(ovs.coupdotvel(T1, 0, 1))))

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
    print('phase:', T1.phase)
    print('constant_term_phase:,', 0.500 * np.sum(T1.getwidth_traj() / T1.getmass_traj()))
    # energy2 = B.Traj[1].getpotential_traj() + B.Traj[1].getkineticlass() - ekin_tr
    # print('Energy1: ', energy2)
    # plt.scatter(t, np.double(energy), c='blue')
    # plt.pause(0.1)


    amps[i, 0] = np.abs(T1.stateAmpE[0]) ** 2
    amps[i, 1] = np.abs(T1.stateAmpE[1]) ** 2
    Told = copy(T1)
    if Told.getcoupling_traj(0,1)[0]>0:
        phasewf=1
    else:
        phasewf=-1

    os.system('cp /home/AndresMoreno/wfu/003.molpro.wfu /home/AndresMoreno/wfu/002.molpro.wfu')
    os.system('cp 003.molpro.wfu 002.molpro.wfu')
    T2 = velocityverlet(Told, dt, i, calc1,phasewf)

    energy2 = T2.getpotential_traj() + T2.getkineticlass() - ekin_tr
    print('coupling',T2.getcoupling_traj(0,1)[0])
    # if T0.getcoupling_traj(0, 1)[0] / T2.getcoupling_traj(0, 1)[0] < 0 and abs(T0.getcoupling_traj(0, 1)[0]-T2.getcoupling_traj(0, 1)[0])>1.5:
    #     print('phase changed')
    #     print(T2.getcoupling_traj(0, 1)[0])
    #     derivs = np.zeros((T2.ndim, T2.nstates, T2.nstates))
    #     for n1 in range(T2.nstates):
    #         for n2 in range(T2.nstates):
    #             if n1 != n2:
    #                 derivs[:, n1, n2] = -T2.getcoupling_traj(n1, n2)
    #             else:
    #                 derivs[:, n1, n1] = T2.getforce_traj(n1)
    #
    #     T2.setderivs_traj(derivs)
    #     print(T2.getcoupling_traj(0, 1)[0])
    #     print(T2.getcoupling_traj(1, 0)[0])
    #     print(T2.get_traj_force())
    #     print(np.abs(T2.getamplitude_traj())**2)


    print('Endif:', np.abs(energy1 - energy2))
    # if np.abs(energy1 - energy2) > 9E-7:
    #     print('Trying to adapt time-step')
    #     print(T1.getmomentum_traj()[0])
    #     os.system('cp /home/AndresMoreno/wfu/002.molpro.wfu /home/AndresMoreno/wfu/003.molpro.wfu')
    #     os.system('cp 002.molpro.wfu 003.molpro.wfu')
    #     #     os.system('cp /home/AndresMoreno/wfu/' + str(i - 1) + '.check.wfu 003.molpro.wfu')
    #     dt_2 = dt /5.000
    #     for j in range(5):
    #         #T1 = velocityverlet(T1, 1.25, i, calc1)
    #         T1= rk_method_4(T1, 0.5, i, calc1)
    #         print(T1.getmomentum_traj()[0])
    # #     print('Adaptative timestep')
    # #     T1=Told
    # #     energy2 = T1.getpotential_traj() + T1.getkineticlass() - ekin_tr
    # #     print('Adapif:', np.abs(energy1 - energy2))
    # else:
    T1 = copy(T2)

    t = t + dt
    wrtout(False, T1, t)
plt.plot(time, np.abs(amps[:, 0]), c='blue')
plt.plot(time, np.abs(amps[:, 1]), c='red')

plt.show()
