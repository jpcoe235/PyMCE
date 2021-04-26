import numpy as np

from src import overlaps as ov
from src import bundle
from src import MatrixInv
from scipy import linalg as sla
from src import GJ_matrix_inv
from src import Sinveigs as la
def buildsandh(B):

    overlap_tresh = 1e-8
    ntraj = B.ntraj
    nstates = B.nstates
    S = np.zeros((ntraj, ntraj),dtype=np.complex128)
    H = np.zeros((ntraj, ntraj),dtype=np.complex128)
    T = np.zeros((ntraj, ntraj),dtype=np.complex128)
    V = np.zeros((ntraj, ntraj),dtype=np.complex128)
    V0 = np.zeros((ntraj, ntraj),dtype=np.complex128)
    V1 = np.zeros((nstates,ntraj, ntraj),dtype=np.complex128)
    SE = np.zeros((ntraj, ntraj),dtype=np.complex128)
    Sdot = np.zeros((ntraj, ntraj),dtype=np.complex128)
    # Sdot_fb = np.zeros((ntraj, ntraj))
    # Sinv = np.zeros((ntraj, ntraj))
    # Sp5i = np.zeros((ntraj, ntraj))
    # Sp5 = np.zeros((ntraj, ntraj))
    # Heff = np.zeros((ntraj, ntraj))
    # Heff1 = np.zeros((ntraj, ntraj))
    Ampdot = np.zeros(ntraj)
    D0 = np.zeros((nstates, nstates, ntraj, ntraj))
    D1 = np.zeros((nstates, nstates, ntraj, ntraj))

    for i in range(ntraj):
        for j in range(i,ntraj):
            T1 = B.Traj[i]
            T2 = B.Traj[j]

            S[i, j] = ov.overlap_trajs(T1, T2)

            print('i,j: ',S[i,j],i,j)
            print('AMPSE: ', T1.stateAmpE)
            SE[i, j] = np.dot(T1.stateAmpE, T2.stateAmpE)
            if abs(S[i, j]) < overlap_tresh:
                S[i, j] = 0.0
                H[i, j] = 0.0
                V1[:,i, j] = 0.0
                V[i, j] = 0.0
                T[i, j] = 0.0
                D0[:,:,i, j] = 0.0
                D1[:,:,i, j] = 0.0

            else:
                if i==j:
                    V0[i,j]=0.
                    V1[:,i,j]=ov.overlap_v_traj(T1,T2)
                    V[i,j]=ov.overlap_v_traj(T1,T2)
                    T[i,j]=ov.overlap_ke_traj(T1,T2)
                    H[i,j]=T[i,j]+V[i,j]
                else:
                    V1[:,i,j]=0.0
                    V[i,j]=ov.overlap_v_traj(T1,T2)
                    V0[i,j]=ov.overlap_v_traj(T1,T2)
                    T[i,j]=ov.overlap_ke_traj(T1,T2)

                    H[i,j]=T[i,j]+V[i,j]

                    H[j,i]=np.conj(H[i,j])
                    T[j, i] = np.conj(T[i, j])
                    V[j, i] = np.conj(V[i, j])
                    V0[j,i]=V0[i,j]
                    V1[:,j,i]=V1[:,i,j]
                    S[j,i]=np.conj(S[i,j])
                    SE[j,i]=np.conj(SE[i,j])

                    Sdot[i,j]=ov.overlap_sdot_traj(T1,T2)
                    Sdot[j,i]=ov.overlap_sdot_traj(T2,T1)

    # Oimg=np.real(-1j*S*SE)
    # Oreal=np.real(S*SE)
    # wimag,evimag=np.linalg.eig(Oimg)
    # wreal, evreal = np.linalg.eig(Oreal)
    # CEvec=evreal+1j*evimag
    #
    print(' in HS :',T)
    print('S in Hamilt :', S)
    print('SE in hamilt :', SE)

    fullmat=S*SE

    Sinv=la.Sinv2(fullmat)
    #Sinv=np.linalg.inv(fullmat)

    #Sinv2=MatrixInv.getMatrixInverse(fullmat)

    #Sinv3=GJ_matrix_inv.MatrixInversion.solve(-1j*S*SE)

    print('Fullmat in hamilt: ', fullmat)
    print('Sinv in hamilt: ', Sinv)
    # print('Sinv2 in hamilt: ', Sinv3)

    tmpmat=H-1j*Sdot
    Heff=np.matmul(Sinv,tmpmat)

    Heff1=H-1j*Sdot

    '''Set variables'''
    B.setH_bundle(H)
    B.setHeff1_bundle(Heff1)
    B.setHeff_bundle(Heff)
    B.setS_bundle(S)
    B.setSE_bundle(SE)
    B.setSinv_bundle(Sinv)
    B.setT_bundle(T)
    B.setV_bundle(V)
    return B









