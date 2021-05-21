import numpy as np



def branching_ratios(B):
    '''Branching of a bundle with ntraj trajectories, this function needs to be refactored'''

    ntraj = B.ntraj
    nstates = B.nstates

    # c_l=np.zeros(ntraj,dtype=np.complex128)
    # for i in range(ntraj):
    #     c_l[i] = B.Traj[i].amp

    c_l=B.getamps_bundle()

    print('c_l: ',c_l)
    cse = np.zeros((ntraj, nstates),dtype=np.complex128)

    for i in range(ntraj):
        cse[i,:] = B.Traj[i].stateAmpE

    print('cse: ',cse)
    s_l = B.get_calc_set_overlaps()
    print('Overlaps branching: ', s_l)
    ratio = np.zeros(nstates)

    for i in range(nstates):
        ratio[i] = np.real(np.dot(c_l * cse[:,i], np.matmul(s_l, c_l * cse[:, i])))


    norm=np.sum(ratio)
    return ratio


