import input
import Wigner_dist
import numpy as np
import Traj
import molpro_call_read
import os
import time

'''AIMCE propagation'''

''' First, common variables are initialized using the class Traj, based on the dynamics and geometry inputs'''
'''q,p initial values, taken from the wigner dist and d,s,a'''

geo = input.initgeom()
dyn = input.initdyn()
tr = Traj.traj()

qin, pin = Wigner_dist.WignerSampling()

tr.n = 1
tr.iprop = 0

tr.q = qin
tr.p = pin
tr.d = np.zeros(dyn.nstates, dtype=complex)
tr.d[dyn.inipes - 1] = 1.0 + 0j
tr.s = np.zeros(dyn.nstates)
tr.a = tr.d * np.exp(1j * tr.s)
# molpro_call_read.create_input(tr.n,tr.iprop,tr.q)
# if not os.path.isfile('scratch_files/'):
#     os.system('mkdir scratch_files')
# os.system('E:/Molpro/bin/molpro.exe -d scratch_files/ -s molpro_traj_1_0.inp')
#
#
# time_to_wait = 1000
# time_counter = 0
# while not os.path.exists('molpro.pun'):
#     time.sleep(1)
#     time_counter += 1
#     if time_counter > time_to_wait:
#         print('more than 1000 secs')
#         break
#
# print(time_counter)
N,L,M=molpro_call_read.readpun()
print(M[0,2,0,1])

