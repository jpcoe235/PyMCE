import numpy as np
import initialize_traj


def coupdotvel(T, i, j):
    if i == j:
        coup = 0.0
    else:
        coup = np.dot(T.getvelocity_traj(), T.getcoupling_traj( i, j))

    return coup

def overlap_CG(x1, x2, p1, p2, a1, a2):
    dx = x1 - x2
    dp = p1 - p2
    real_part = (a1 * a2 * dx ** 2 + 0.25 * dp ** 2) / (a1 + a2)
    pref = np.sqrt(2 * np.sqrt(a1 * a2) / (a1 * a2))
    x_cent = (a1 * x1 + a2 * x2) / (a1 + a2)
    imag_part = (p1 * x1 + p2 * x2) - x_cent * dp
    sij = pref * np.exp(-real_part + 1j * imag_part)

    return sij


def overlap_d2x_cg(x1, x2, p1, p2, a1, a2):
    dx=x1-x2
    p12=a1*p2+a2*p1
    d2xsij=(4j*a1*a2*dx*p12+2*a1*a2*(a1+a2)-4*dx**2*a1**2*a2**2+p12**2)/(a1+a2)**2
    d2xsij*=overlap_CG(x1, x2, p1, p2, a1, a2)

    return d2xsij


def overlap_dx_cg(x1, x2, p1, p2, a1, a2):
    dx = x1 - x2
    p12 = a1 * p2 + a2 * p1
    dxCG=(-p12*1j-2*a1*a2*dx)/(a1+a2)
    dxCG*=overlap_CG(x1, x2, p1, p2, a1, a2)
    return dxCG

def overlap_trajs(T1, T2):

    amps1 = T1.getamp_traj()
    amps2 = T2.getamp_traj()
    ph1 = T1.getphase_traj()
    ph2 = T2.getphase_traj()

    ndim = T1.ndim

    S = np.exp(1j * (-ph1 + ph2), dtype=np.complex128)
    for i in range(ndim):
        a1 = amps1[i]
        a2 = amps2[i]
        x1 = T1.getposition_traj()[i]
        x2 = T2.getposition_traj()[i]
        p1 = T1.getposition_traj()[i]
        p2 = T2.getposition_traj()[i]

        S = S * overlap_CG(x1, x2, p1, p2, a1, a2)

    return S


def overlap_dp_traj(T1,T2):

    x1=T1.getposition_traj()
    x2=T2.getposition_traj()

    p1=T1.getmomentum_traj()
    p2=T2.getmomentum_traj()

    a1=T1.getphase_traj()
    a2=T2.getphase_traj()

    dx = x1 - x2
    dp = p1 - p2

    ndim=T1.ndim
    pref=np.zeros(ndim)

    for i in range(ndim):
        pref[i]=0.5*(dp[i]+2j*a1[i]*dx[i])/(a1[i]+a2[i])
        pref[i]*=overlap_CG(x1[i],x2[i],p1[i],p2[i],a1[i],a2[i])

    dpSij=pref*overlap_trajs(T1,T2)

    return dpSij





def overlap_v_traj(T1,T2):
    ind1=T1.ntraj()
    ind2=T2.ntraj()

    if ind1==ind2:
        V=T1.getpotential_traj()
    else
        nstates1=T1.nstates
        for i in range(nstates1):
            T1.HE[i,i]=T1.getpotential_traj_i(i)
            for j in range(i+1,nstates1):
                T1.HE[i,j]=-1j*coupdotvel(T1,i,j)
                T1.HE[j,i]=-T1.HE[i,j]
        nstates2 = T2.nstates
        for i in range(nstates2):
            T2.HE[i, i] = T2.getpotential_traj_i(i)
            for j in range(i + 1, nstates2):
                T2.HE[i, j] = -1j * coupdotvel(T2, i, j)
                T2.HE[j, i] = -T2.HE[i, j]

        Sdp=overlap_dp_traj(T1,T2)

        F1=T1.get_traj_force()
        F2=T2.get_traj_force()

        V=np.dot(T1.getamplitude_traj(),np.matmul(0.5*(T1.HE+T2.HE),T2.getamplitude_traj()))*overlap_trajs(T1,T2)+\
    0.5*1j*np.sum(np.real(Sdp))*(F1+F2)*np.dot(T1.getamplitude_traj(),T2.getamplitude_traj())


def overlap_ke_traj(T1,T2):
    ke=0.0
    sij=overlap_trajs(T1,T2)
    ndim=T1.ndim
    x1 = T1.getposition_traj()
    x2 = T2.getposition_traj()

    p1 = T1.getmomentum_traj()
    p2 = T2.getmomentum_traj()

    a1 = T1.getphase_traj()
    a2 = T2.getphase_traj()

    for i in range(ndim):

        ke-=overlap_d2x_cg(x1[i],x2[i],p1[i],p2[i],a1[i],a2[i])

    ke=ke*sij
    ke=ke*np.dot(T1.getamplitude_traj(),T2.getamplitude_traj())
    return ke

def overlap_dx_traj(T1,T2):
    sij = overlap_trajs(T1, T2)
    ndim = T1.ndim
    x1 = T1.getposition_traj()
    x2 = T2.getposition_traj()

    p1 = T1.getmomentum_traj()
    p2 = T2.getmomentum_traj()

    a1 = T1.getphase_traj()
    a2 = T2.getphase_traj()

    dxsij=np.zeros(ndim)
    for i in range(ndim):
        dxsij[i]=overlap_dx_cg(x1[i],x2[i],p1[i],p2[i],a1[i],a2[i])

    dxsij=dxsij*sij

    return dxsij

def overlap_sdot_traj(T1,T2):

    Sij=overlap_trajs(T1,T2)
    Sdot= np.dot(T2.getvelocity_traj(),overlap_dx_traj(T1,T2))+np.dot(T2.get_traj_force(), overlap_dp_traj(T1,T2))+1j*T2.phasedot()*Sij


    nstates=T2.nstates
    HE=T2.getHE_traj()

    if HE==0.0:
        HE=T2.get_calc_HE()
        T2.setHE_traj(HE)

    SdotE=-1j *np.dot(T1.stateAmpE,np.matmul(T2.HE,T2.stateAmpE))
    Sdot=Sdot*np.dot(T1.stateAmpE,T2.stateAmpE)+SdotE*Sij

    return Sdot