from src.geometry import initgeom
from src.dyn_params import initdyn
from src import initialize_traj
from src import abinitio
from src.constants import physconst
from src.propagators import velocityverlet
from src.propagators import velocityverlet_dima
from src.propagators import velocityverlet_restart
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
from src.overlaps_wf import overlap as ovwf
from src.Full_traj import full_trajectory

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
# The 3 here referes to the degrees of freedom per atom, in case you want a 1D system or 2D
# qin, pin = Wigner_dist.WignerSampling()
# print(qin,pin)
# exit()
# q=qin
# p=pin
q = np.zeros(geo.ndf, dtype=np.longdouble)
p = np.zeros_like(q, dtype=np.longdouble)

sharc = True
numtraj = 4
count = 0

if sharc:
    index1 = False
    index2 = False
    atomc = 1
    #   os.system('python src/wigner.py -n 1 vibs.molden')
    with open('initconds', 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('Index'):
                index1 = line.strip().split()
                print('reading initial conditions ', index1[1])
                index1 = True
            elif index1:
                index2 = True
                index1 = False

            elif index2:

                index1 = False
                line1 = line.strip().split()
                q[count:count + 3] = line1[2:5]
                p[count:count + 3] = line1[6:9]
                count = count + 3
                atomc += 1
                if atomc > geo.natoms:
                    index1 = False
                    index2 = False
                    atomc = 1

else:
    with open('initial_pyrazine.dat', 'r') as f:
        f.readline()
        for i in range(geo.ndf):
            N, M = f.readline().strip().split()
            q[i] = np.longdouble(float(N.replace('D', 'E')))
            p[i] = np.longdouble(float(M.replace('D', 'E')))

T1.setmassall_traj(geo.massrk)
T1.setposition_traj(q)
T1.setmomentum_traj(p * T1.getmassall_traj())
# Sharc calculates velocities, we need to multiply by the mass to get the momentum

pec, der, cis, configs = abinitio.inp_out(0, 0, geo, T1)  # First ab-initio run

T1.setpotential_traj(pec)  # taking V(R) from ab-initio
T1.setderivs_traj(
    der)  # derivatives matrix not mass-weighted (possibly change that), diagonals are forces and off-d are nacmes
T1.setcivecs(cis)

f = open("coups_forces_first.dat", 'w', buffering=1)
for i in range(T1.nstates):
    for j in range(T1.nstates):
        for k in range(geo.ndf):
            f.write(str(der[k, i, j]) + '\n')
        f.write('\n')
f.close()

print('dermatrix=', der)
print('Energies= ', pec)

T1.setmass_traj(geo.masses)  # mass of every atom in a.u (the dimmension is natoms/nparts)
T1.setmassall_traj(
    geo.massrk)  # mass in every degree of freedom (careful to use it, it can triple the division/multiplication easily)
T1.setcivecs(cis)
T1.setconfigs(configs)
amps = np.zeros(T1.nstates, dtype=np.clongdouble)
amps[dyn.inipes - 1] = np.clongdouble(
    1.00000000 + 0.000000000j)  # Amplitudes of Ehrenfest trajectories, they should be defined as a=d *exp(im*S)

T1.setamplitudes_traj(amps)
T1.setd_traj(amps)
phases = np.zeros(T1.nstates)  # Phase of the wfn, would be S in the previous equation
T1.setphases_traj(phases)
T1.setwidth_traj(dyn.gamma)

# We need to change gamma if we plan to couple the trajectories afterwards


dt = 0.1 * 1e-15 / ph.au2sec
# dt = np.longdouble(4.130000) #0.1 fs
print('dt:', dt * ph.au2sec / 1e-15)

amps = np.zeros((15000, 2))
ekin_tr = 0
t = 0
calc1 = False
# wrtout(True, T1, 0.000)
repeat = True
T0 = copy(T1)
# fig1, ax1 = plt.subplots()

# ax1.axes.set_xlabel('Time(fs)')
# ax1.axes.set_ylabel('E difference')
f = open("N_mine.dat", 'w', buffering=1)
f2 = open("CS_mine.dat", 'w', buffering=1)
f3 = open("F_mine.dat", 'w', buffering=1)
colorvec = ['blue', 'green', 'black', 'yellow']

for i in range(1):
    if i == 0:
        T1 = copy(T0)

    # ratios=br.branching_ratios(B)

    # energy2 = B.Traj[1].getpotential_traj() + B.Traj[1].getkineticlass() - ekin_tr
    # print('Energy1: ', energy2)

    # ax = fig.add_subplot(111)

    # ax2.scatter(t * ph.au2sec / 1e-15, np.abs(T1.stateAmpE[0]) ** 2, c='red')
    # ax2.scatter(t * ph.au2sec / 1e-15, np.abs(T1.stateAmpE[1]) ** 2, c='blue')
    # plt.pause(0.1)
    f.write(str(t) + ' ' + str(np.abs(T1.stateAmpE[0]) ** 2) + ' ' + str(np.abs(T1.stateAmpE[1]) ** 2) + '\n')
    f2.write(str(t) + ' ' + str(T1.getcoupling_traj(0, 1)[0]) + '\n')

    # amps[i, 0] = np.abs(T1.stateAmpE[0]) ** 2
    # amps[i, 1] = np.abs(T1.stateAmpE[1]) ** 2
    Told = copy(T1)
    # if Told.getcoupling_traj(0, 1)[0] > 0:
    #     phasewf = 1
    # else:
    #     phasewf = -1
    phasewf = 1
    finaltime=600
    finaltime = finaltime / (ph.au2sec / 1E-15)
    time_vec = np.linspace(0.000, finaltime,int(finaltime/dt))
    # os.system('cp /home/AndresMoreno/wfu/003.molpro.wfu /home/AndresMoreno/wfu/002.molpro.wfu')
    # os.system('cp 003.molpro.wfu 002.molpro.wfu')
    restarting = False
    if not restarting:
        FT, T2, Bundle = velocityverlet_dima(Told, finaltime, dt, i + 1, numtraj, calc1, phasewf, 0.0000)
    else:

        NN = 90
        time = time_vec[NN-1]
        FT, T2, Bundle = velocityverlet_restart(finaltime, dt, NN, numtraj, calc1, phasewf, time)

    # exit
    # for ns in range(T2.nstates):
    #     ovs_ci = np.dot(T2.getcivecs()[:, ns], T1.getcivecs()[:, ns])
    #     if ovs_ci < 0.0:
    #         ccoo = 'red'
    #     else:
    #         ccoo = colorvec[ns]
    # #   ax3.scatter(t*ph.au2sec/1e-15, abs(ovs_ci),c=ccoo)
    # ov_11, ov_22 = ovwf(T1, T2)
    # if ov_11 < 0.0:
    #     ccoo = 'red'
    # else:
    #     ccoo = colorvec[2]
    # # ax3.scatter(t * ph.au2sec / 1e-15, abs(ov_11), c=ccoo)
    # if ov_22 < 0.0:
    #     ccoo = 'red'
    # else:
    #     ccoo = colorvec[3]
    # # ax3.scatter(t * ph.au2sec / 1e-15, abs(ov_22), c=ccoo)
    # # plt.pause(0.01)
    # energy2 = T2.getpotential_traj() + T2.getkineticlass()
    # f3.write(str(t) + ' ' + str(energy2) + '\n')
    # print('coupling', T2.getcoupling_traj(0, 1)[0])
    #
    # T1 = copy(T2)
    #
    # t = t + dt
# wrtout(False, T1, t)
f.close()
f2.close()
f3.close()
# fig1.savefig('energydif.png')
# fig2.savefig('amplitudes.png')
# fig3.savefig('overlaps.png')
