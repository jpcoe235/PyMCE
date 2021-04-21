import input
import Constants
import Traj
import numpy as np
import os

''' Routine to call and read parameters from molpro, it is recommended to check the convergence of the calculation before starting
a dynamics calculation'''

'''First we create the input, no symmetry, no oriented and using the parameters specified by the user'''

'''ab_par is a class for the abinitio calculation parameters, NPI should be included'''


class ab_par():
    def __init__(self):
        self.molecule = 'Ethylene'
        self.act_orb = 9
        self.closed_orb = 7
        self.basis = 'avdz'
        self.civec = False
        self.first = True
        self.molden = False
        self.n_el = 16


def inp_out(i, substep, q,geo):
    os.system('rm molpro.pun')
    os.system('rm molpro_traj*')
    create_input(i, substep, q,geo)
    os.system('E:/Molpro/bin/molpro.exe -d . -s molpro_traj_' + str(i) + '_' + str(substep) + '.inp')
    time_counter = 0

    while not os.path.exists('molpro.pun'):
        time.sleep(1)
        time_counter += 1
        if time_counter > time_to_wait:
            print('more than 1000 secs')
            break
    N, L, M = readpun()

    pes, grads, nacs = update_vars(N, L, M,geo)

    return pes, grads, nacs


def transform_q_to_r(q,geo):

    rkinit_1 = geo.rkinit

    ph = Constants.physconst()
    k = 0
    r = np.zeros((3, geo.natoms))

    diff=np.double(q + rkinit_1)
    for i in range(geo.natoms):
        for j in range(3):
            mass_inv = 1.0 / np.sqrt(geo.masses[i])
            r[j, i] = mass_inv * (diff[k]) * ph.bohr
            k += 1

    return r


def create_input(trajN, istep, q,geo):
    ab = ab_par()

    dyn = input.initdyn()
    if trajN == 0 and istep == 0:
        first = True
    else:
        first = False

    '''routine to create molpro input, valid for molpro2012'''

    file = 'molpro_traj_' + str(trajN) + '_' + str(istep) + '.inp'

    r = transform_q_to_r(q,geo)

    # if istep==0 and trajN==0:
    #     r=np.transpose(np.asarray([[7.0544201503E-01,-8.8768894715E-03,-6.6808940143E-03],[-6.8165975936E-01,1.7948934206E-02,4.0230273972E-03],[ 1.2640943417E+00,9.0767471618E-01,-3.0960211126E-02],[-1.4483835697E+00 ,8.7539956319E-01 ,4.6789959947E-02],[1.0430033100E+00,-9.0677165411E-01 ,8.1418967247E-02],[-1.1419988770E+00,-9.8436525752E-01,-6.5589218426E-02]]))
    with open(file, 'w') as f:
        f.write(
            '***,' + ab.molecule + ' ' + 'calculation of ' + str(istep) + ' step in trejectory ' + str(trajN) + '\n')
        if first:
            line_wr = 'file,3,molpro.wfu,new\n'
        else:
            line_wr = 'file,3,molpro.wfu\n'
        f.write(line_wr)
        f.write('memory,100,m\n')
        f.write('punch,molpro.pun,new\n')
        f.write('basis={\n')
        f.write('''!
! HYDROGEN       (4s,1p) -> [2s,1p]
 s, H ,   13.0100000,    1.9620000,   0.4446000
 c, 1.3,   0.0196850,    0.1379770,   0.4781480
 s, H ,    0.1220000
 c, 1.1,   1.0000000
 p, H ,    0.7270000
 c, 1.1,   1.0000000
!
!
! CARBON       (9s,4p,1d) -> [3s,2p,1d]
 s, C , 6665.0000000, 1000.0000000, 228.0000000, 64.7100000, 21.0600000,  7.4950000,  2.7970000, 0.5215000
 c, 1.8,   0.0006920,    0.0053290,   0.0270770,  0.1017180,  0.2747400,  0.4485640,  0.2850740, 0.0152040
 s, C , 6665.0000000, 1000.0000000, 228.0000000, 64.7100000, 21.0600000,  7.4950000,  2.7970000, 0.5215000
 c, 1.8,  -0.0001460,   -0.0011540,  -0.0057250, -0.0233120, -0.0639550, -0.1499810, -0.1272620, 0.5445290
 s, C ,    0.1596000
 c, 1.1,   1.0000000
 p, C ,    9.4390000,    2.0020000,   0.5456000
 c, 1.3,   0.0381090,    0.2094800,   0.5085570
 p, C ,    0.1517000
 c, 1.1,   1.0000000
 d, C ,    0.5500000
 c, 1.1,   1.0000000
}
''')
        f.write('''symmetry,nosym;
orient,noorient;
angstrom;
geomtype=xyz;
geom={
        ''')
        f.write(str(geo.natoms) + '\n')
        f.write('\n')
        for i in range(geo.natoms):
            line_wr = geo.atnames[i] + ' ' + str(r[0, i]) + ' ' + str(r[1, i]) + ' ' + str(r[2, i]) + '\n'
            # print(file)
            # print(line_wr)
            f.write(line_wr)
        f.write('}\n')
        f.write('''{multi,failsafe;
maxiter,40;
''')
        line_wr = 'occ,' + str(ab.act_orb) + ';\n'
        line_wr2 = 'closed,' + str(ab.closed_orb) + ';\n'
        line_wr3 = 'wf,' + str(ab.n_el) + ',1,0;'
        line_wr4 = 'state,' + str(dyn.nstates) + ';\n'

        f.write(line_wr)
        f.write(line_wr2)
        f.write(line_wr3)
        f.write(line_wr4)
        f.write('''weight,1,1;
orbital,2140.3;
canonical,2140.2,ci;
ORBITAL,IGNORE_ERROR;
''')

        record = 5100.1
        for i in range(dyn.nstates):
            for j in range(i, dyn.nstates):

                if i == j:

                    f.write('CPMCSCF,GRAD,' + str(i + 1) + '.1,record=' + str(record) + ';\n')
                    record += 1
                else:
                    f.write('CPMCSCF,NACM,' + str(i + 1) + '.1,' + str(j + 1) + '.1,' + 'record=' + str(record) + ';\n')
                    record += 1
        f.write('}\n')
        record = 5100.1

        for i in range(dyn.nstates):
            for j in range(i, dyn.nstates):

                if i == j:

                    f.write('{FORCE;SAMC,' + str(record) + '};\n')
                    record += 1
                else:
                    f.write('{FORCE;SAMC,' + str(record) + '};\n')
                    record += 1
        f.write('''{OPTG,maxit=1;
coord,3n,norot;
}
---''')


def readpun():
    dyn = input.initdyn()
    geo = input.initgeom()
    v_c = np.zeros(dyn.nstates)
    grad = np.zeros((3, geo.natoms, dyn.nstates))
    nacmes = np.zeros((3, geo.natoms, dyn.nstates, dyn.nstates))
    with open('molpro.pun', 'r') as f:
        cV = 0
        cPes = 0
        cNacs1 = 0
        cNacs2 = cNacs1 + 1
        for lines in f:

            if 'MCSCF STATE ' in lines:
                string = lines.strip().split()
                if string[3] == 'Energy':
                    v_c[cV] = np.double(string[4])
                    cV += 1
            if 'SA-MC GRADIENT' in lines:
                string = lines.strip().split()
                atom = int(string[3].replace(':', '')) - 1
                grad[0, atom, cPes] = float(string[4])
                grad[1, atom, cPes] = float(string[5])
                grad[2, atom, cPes] = float(string[6])
                if atom == geo.natoms - 1:
                    cPes += 1
            if 'SA-MC NACME' in lines:
                string = lines.strip().split()
                atom = int(string[3].replace(':', '')) - 1
                nacmes[0, atom, cNacs1, cNacs2] = float(string[4])
                nacmes[1, atom, cNacs1, cNacs2] = float(string[5])
                nacmes[2, atom, cNacs1, cNacs2] = float(string[6])
                if atom == geo.natoms - 1:
                    if cNacs2 == dyn.nstates - 1:
                        cNacs1 += 1
                    else:
                        cNacs2 += 1
    # for i in range(2):
    #     for j in range(geo.natoms):
    #         print(grad[0, j, i], grad[1, j, i], grad[2, j, i])
    return v_c, grad, nacmes


def update_vars(v_c, grad, nacmes,geo):

    dyn = input.initdyn()
    '''Update the variables and store them using vectorization'''
    pes = np.zeros(dyn.nstates)
    grads = np.zeros((dyn.nstates, geo.ndf))
    nacs = np.zeros((dyn.nstates, dyn.nstates, geo.ndf))
    for i in range(
            dyn.nstates):  # Updates pes energy with the reference value (value corresponding to the initial geometry GS)
        pes[i] = v_c[i]-dyn.e_ref

    for n in range(dyn.nstates):
        idf = 0
        for i in range(geo.natoms):
            for j in range(3):
                grads[n, idf] = grad[j, i, n] / np.sqrt(geo.masses[i])
                idf += 1

    for n in range(dyn.nstates):
        for k in range(n + 1, dyn.nstates):
            idf = 0
            for i in range(geo.natoms):
                for j in range(3):
                    nacs[n, k, idf] = nacmes[j, i, n, k] / np.sqrt(geo.masses[i])
                    nacs[k, n, idf] = -nacmes[j, i, n, k] / np.sqrt(geo.masses[i])
                    idf += 1

    return pes, grads, nacs
