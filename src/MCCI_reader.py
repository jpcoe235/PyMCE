import numpy as np
import os



def MCCI_reader(q,geo,dyn):
   # os.chdir('../')
    with open('InputGeom.txt', 'w') as f:
        for i in q:
            f.write(str(float(i)) + '\n')

    os.chdir('filesneeded')
    # os.system('cd filesRequired_v4_CSFs')
    os.system('cp ../InputGeom.txt .')
    os.system('python3 MCCI_output_v2.py')
    pes = np.zeros(dyn.nstates, dtype=np.double)
    grads = np.zeros((geo.ndf, dyn.nstates), dtype=np.double)
    nacs = np.zeros((geo.ndf, dyn.nstates, dyn.nstates), dtype=np.double)
    with open('FinalE', 'r') as f:
        for i in range(dyn.nstates):
            pes[i] = f.readline()

    with open('NACresults.txt', 'r') as f:
        for i in range(dyn.nstates):
            for j in range(i + 1, dyn.nstates):
                f.readline()
                print(i, j)
                for n in range(geo.natoms):
                    M = f.readline()
                    M = M.strip().split()
                    print(M)
                    nacs[3 * n:3 * n + 3, i, j] = M[1:]

    for i in range(dyn.nstates):
        with open(str("SCIgradient_state" + str(i + 1) + '.txt'), 'r') as f:
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
    os.chdir("../")
    return pes, der