import numpy as np
import os



def MCCI_reader(q,geo,dyn):
   # os.chdir('../')
    os.chdir('bathCI_1/files_required_1')
    os.system('rm -rf *')
    with open('InputGeom.txt', 'w') as f:
        for i in q:
            f.write(str(float(i)) + '\n')



    # os.system('cd filesRequired_v4_CSFs')


    os.system('cp ../run/config.json .')
    os.system('cp ../run/SHCIconfig.json .')
    os.system('cp ../run/geom_ini.xyz .')
    os.system('python3 ../run/SHCI_nacsandgrads.py >log.txt')
    pes = np.zeros(dyn.nstates, dtype=np.double)
    grads = np.zeros((geo.ndf, dyn.nstates), dtype=np.double)
    nacs = np.zeros((geo.ndf, dyn.nstates, dyn.nstates), dtype=np.double)
    with open('Singlet_FinalE', 'r') as f:
        for i in range(dyn.nstates):
            pes[i] = f.readline()


    for i in range(dyn.nstates):
        for j in range(i + 1, dyn.nstates):
            with open("NAC_singlets_state" + str(i + 1) +"state"+ str(j + 1)+'.txt') as f:
                print(i, j)
                for n in range(geo.natoms):
                    M = f.readline()
                    M = M.strip().split()
                    print(M)
                    nacs[3 * n:3 * n + 3, i, j] = M

    for i in range(dyn.nstates):
        with open(str("SCIgradient_Singlet_state" + str(i + 1) + '.txt'), 'r') as f:
            for n in range(geo.natoms):
                M = f.readline()
                M = M.strip().split()
                print(M)
                grads[3 * n:3 * n + 3, i] = M[0:]

    pes = np.asarray(pes)
    nacs = np.asarray(nacs)

    der = np.zeros(
        (geo.ndf, dyn.nstates, dyn.nstates))  # Define variable der, grads in the diag and nacmes in the offdiags
    for i in range(dyn.nstates):
        der[:, i, i] = -grads[:, i]
        for j in range(i + 1, dyn.nstates):
            der[:, i, j] = nacs[:, i, j]
            der[:, j, i] = -nacs[:, i, j]
    os.chdir("../../")
    return pes, der