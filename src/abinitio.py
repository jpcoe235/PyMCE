from src.constants import physconst
from src.dyn_params import initdyn
from src.geometry import initgeom
import time
import numpy as np
import os
from decimal import Decimal

''' Routine to call and read parameters from molpro, it is recommended to check the convergence of the calculation before starting
a dynamics calculation'''

'''First we create the input, no symmetry, no oriented and using the parameters specified by the user'''

'''ab_par is a class for the abinitio calculation parameters, NPI should be included'''


class ab_par():
    def __init__(self):
        self.molecule = 'Pyrazyne'  # Name of the molecule, it is not really used for anything
        self.act_orb = 21  # Active orbitals
        self.closed_orb = 13  # Closed orbitals
        self.basis = 'avdz'  # Basis, not really used as the basis should be partitioned (pople or copying them from exchange)
        self.civec = False  # prints CIs
        self.first = True  # First calculation creates the wf file
        self.molden = False  # create a molden
        self.n_el = 36  # number of electrons to define the wf


def inp_out(i, substep, geo, T1):
    os.system('rm molpro.pun')
   # os.system('rm molpro_traj*')
    q = T1.getposition_traj()

    file2 = create_input(i, substep, q, geo)  # call to create_inp function

    T1.setfilecalc(file2)
    os.system('/usr/bin/molpro -d . -s molpro_traj_' + str(i) + '_' + str(
        substep) + '.inp')  # running molpro, change it for any run in a different computer
    time_counter = 0  #
    time_to_wait = 10000000
    while not os.path.exists(
            'molpro.pun'):  # Dodgy way to wait for the punch file to be created, must be other way more elegant
        time.sleep(1)
        time_counter += 1
        if time_counter > time_to_wait:
            print('more than 1000 secs')
            break
    N, L, M, cis, configs = readpun()  # Read the punch file, only PES, GRADS and NACMES are read in a (natoms,3) format

    pes, grads, nacs, = update_vars(N, L, M,
                                    geo)  # Update the vars to be vectors, it gives the option to mass-weigth the params

    der = np.zeros(
        (geo.ndf, T1.nstates, T1.nstates))  # Define variable der, grads in the diag and nacmes in the offdiags
    for i in range(T1.nstates):
        der[:, i, i] = grads[:, i]
        for j in range(T1.nstates):
            if i != j:
                der[:, i, j] = nacs[:, i, j]

    return pes, der, cis, configs  # Return both variables, potential energies and derivatives


def transform_q_to_r(q, geo, first):
    first = False
    ph = physconst()
    k = 0
    r = np.zeros((3, geo.natoms))

    diff = np.longdouble(q)
    for i in range(geo.natoms):
        for j in range(3):
            if first:
                mass_inv = 1.0/np.sqrt(geo.masses[i])
                r[j, i] = geo.rk[k] - (mass_inv * (diff[k]) * ph.bohr)
            else:
                mass_inv = 1.0

                r[j, i] = mass_inv * (diff[k]) * np.longdouble(5.2917721092e-1)
               # r[j, i] = geo.rk[k]-(mass_inv * (diff[k]) * ph.bohr)

            k += 1

    return r  # vector->matrix, if first mass weights


def create_input(trajN, istep, q, geo):
    ab = ab_par()
    ph = physconst()
    dyn = initdyn()
    if trajN == 0 and istep == 0:
        first = True
    else:
        first = False
    # Molpro input, it must be changed for other codes
    '''routine to create molpro input, valid for molpro2012'''

    file = 'molpro_traj_' + str(trajN) + '_' + str(istep) + '.inp'
    file2 = 'molpro_traj_' + str(trajN) + '_' + str(istep) + '.mld'
    # if first:
    r = transform_q_to_r(q, geo, first)

    # if istep==0 and trajN==0:
    #     r=np.transpose(np.asarray([[7.0544201503E-01,-8.8768894715E-03,-6.6808940143E-03],[-6.8165975936E-01,1.7948934206E-02,4.0230273972E-03],[ 1.2640943417E+00,9.0767471618E-01,-3.0960211126E-02],[-1.4483835697E+00 ,8.7539956319E-01 ,4.6789959947E-02],[1.0430033100E+00,-9.0677165411E-01 ,8.1418967247E-02],[-1.1419988770E+00,-9.8436525752E-01,-6.5589218426E-02]]))
    with open(file, 'w') as f:
        f.write(
            '***,' + ab.molecule + ' ' + 'calculation of ' + str(istep) + ' step in trejectory ' + str(trajN) + '\n')
        line_wr_2 = 'file,2,' + str(trajN) + '.check.wfu,new\n'
        # f.write(line_wr_2)
        if first:
            line_wr = 'file,3,004.molpro.wfu,new\n'
        else:
            line_wr = 'file,3,004.molpro.wfu\n'
        
        f.write(line_wr)

        f.write('memory,1600,m\n')

        f.write('''gprint,orbitals,civector,angles=-1,distance=-1
 gthresh,twoint=1.0d-13
 gthresh,energy=1.0d-7,gradient=1.0d-2
 gthresh,thrpun=0.001\n''')

        f.write('punch,molpro.pun,new\n')

        if first:
           f.write('''
  SYMMETRY,NOSYM
  GEOMTYPE=XYZ
            ''')
           f.write('geom={' + '\n')
           f.write(str(geo.natoms) + '\n')
           f.write('\n')

           for i in range(geo.natoms):
               line_wr = geo.atnames[i] + ' ' + str("{:.16E}".format(Decimal(r[0, i]))) + ' ' + \
                         str("{:.16E}".format(Decimal(r[1, i]))) + ' ' + str("{:.16E}".format(Decimal(r[2, i]))) + '\n'
               # print(file)
               # print(line_wr)
               f.write(line_wr)
           f.write('}\n')

           f.write('''
            BASIS=sto-3g
  {RHF
  ORBITAL,2100.2
  ORBPRINT,12
  }

  put,molden,in.molden

  {MULTI
  OCC,20
  START,2100.2
  ORBITAL,2070.2
  CANONICAL,2070.2
  CLOSED,14
  ORBPRINT,19
  WF,ELEC=36,SYM=1,SPIN=0
  STATE,1
  }

  put,molden,sto3gorbs.molden

  basis=aug-cc-pvdz

  {MULTI
  OCC,22
  START,2070.2
  ORBITAL,2071.2
  CLOSED,12
  ORBPRINT,19
  WF,ELEC=36,SYM=1,SPIN=0
  STATE,3
  }

  {MULTI
  OCC,21
  START,2071.2
  ORBITAL,2072.2
  CLOSED,13
  ORBPRINT,19
  WF,ELEC=36,SYM=1,SPIN=0
  STATE,3
  CPMCSCF,GRAD,1.1,RECORD=5101.1
  CPMCSCF,GRAD,2.1,RECORD=5102.1
  CPMCSCF,GRAD,3.1,RECORD=5103.1
  CPMCSCF,NACM,1.1,2.1,SAVE=5104.1
  CPMCSCF,NACM,1.1,3.1,SAVE=5105.1
  CPMCSCF,NACM,2.1,3.1,SAVE=5106.1
  }

  force;samc,5101.1;
  force;samc,5102.1; 
  force;samc,5103.1;
  force;samc,5104.1;
  force;samc,5105.1;
  force;samc,5106.1;
''')
           f.write('put,molden, ' + file2 + '\n')
           f.write('''---''')
        else:
            f.write('basis=avdz\n')
        # f.write('basis={\n')
        # f.write('''                                                                               !
        #                                                                                 ! HYDROGEN       (4s,1p) -&gt; [2s,1p]
        #  s, H ,   13.0100000,    1.9620000,   0.4446000
        #  c, 1.3,   0.0196850,    0.1379770,   0.4781480
        #  s, H ,    0.1220000
        #  c, 1.1,   1.0000000
        #  p, H ,    0.7270000
        #  c, 1.1,   1.0000000
        #                                                                                 !
        #                                                                                 !
        #                                                                                 ! CARBON       (9s,4p,1d) -&gt; [3s,2p,1d]
        #  s, C , 6665.0000000, 1000.0000000, 228.0000000, 64.7100000, 21.0600000,  7.4950000,  2.7970000, 0.5215000
        #  c, 1.8,   0.0006920,    0.0053290,   0.0270770,  0.1017180,  0.2747400,  0.4485640,  0.2850740, 0.0152040
        #  s, C , 6665.0000000, 1000.0000000, 228.0000000, 64.7100000, 21.0600000,  7.4950000,  2.7970000, 0.5215000
        #  c, 1.8,  -0.0001460,   -0.0011540,  -0.0057250, -0.0233120, -0.0639550, -0.1499810, -0.1272620, 0.5445290
        #  s, C ,    0.1596000
        #  c, 1.1,   1.0000000
        #  p, C ,    9.4390000,    2.0020000,   0.5456000
        #  c, 1.3,   0.0381090,    0.2094800,   0.5085570
        #  p, C ,    0.1517000
        #  c, 1.1,   1.0000000
        #  d, C ,    0.5500000
        #  c, 1.1,   1.0000000
        # }
        # ''')
            f.write('''symmetry,nosym;
orient,noorient;
angstrom;
geomtype=xyz;
''')
            if 2>3:
                f.write('geom=geom_0.in \n')
            else:
                f.write('geom={'+'\n')
                f.write(str(geo.natoms) + '\n')
                f.write('\n')

                for i in range(geo.natoms):
                    line_wr = geo.atnames[i] + ' ' + str("{:.16E}".format(Decimal(r[0, i]))) + ' ' + \
                              str("{:.16E}".format(Decimal(r[1, i]))) + ' ' + str("{:.16E}".format(Decimal(r[2, i]))) + '\n'
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
            line_wr4 = 'state,' + str(3) + ';\n'

            f.write(line_wr)
            f.write(line_wr2)
            f.write(line_wr3)
            f.write(line_wr4)

            f.write('''pspace,10.0        
orbital,2140.3;
ORBITAL,IGNORE_ERROR;
ciguess,2501.2 
save,ci=2501.2}
''')
            f.write('''data,copy,2140.3,3000.2
''')
            f.write('''{multi,failsafe;
maxiter,40;
''')
            line_wr = 'occ,' + str(ab.act_orb) + ';\n'
            line_wr2 = 'closed,' + str(ab.closed_orb) + ';\n'
            line_wr3 = 'wf,' + str(ab.n_el) + ',1,0;'
            line_wr4 = 'state,' + str(3) + ';\n'

            f.write(line_wr)
            f.write(line_wr2)
            f.write(line_wr3)
            f.write(line_wr4)

            f.write('''pspace,10.0
orbital,2140.3;
dm,2140.3
save,ci=2501.2
diab,3000.2,save=2140.3}
            ''')

            f.write('''{multi,failsafe;
        maxiter,40;
            ''')
            line_wr = 'occ,' + str(ab.act_orb) + ';\n'
            line_wr2 = 'closed,' + str(ab.closed_orb) + ';\n'
            line_wr3 = 'wf,' + str(ab.n_el) + ',1,0;'
            line_wr4 = 'state,' + str(3) + ';\n'

            f.write(line_wr)
            f.write(line_wr2)
            f.write(line_wr3)
            f.write(line_wr4)

            f.write('''pspace,10.0
            orbital,2140.3;
            dm,2140.3
            save,ci=2501.2
            diab,3000.2,save=2140.3
        
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

                        f.write('{FORCES;SAMC,' + str(record) + '};\n')
                        record += 1
                    else:
                        f.write('{FORCES;SAMC,' + str(record) + '};\n')
                        record += 1
            f.write('put,molden, ' + file2 + '\n')
            f.write('''---''')
    return file2


def readpun():
    '''This routine reads a punch molpro file, it outputs as matrices the potential energies, gradients
    and non-adiabatic coupling matrices'''
    dyn = initdyn()
    geo = initgeom()
    v_c = np.zeros(dyn.nstates)
    grad = np.zeros((3, geo.natoms, dyn.nstates))
    nacmes = np.zeros((3, geo.natoms, dyn.nstates, dyn.nstates))
    pos = np.zeros((3, geo.natoms))
    cis = np.zeros((40, 3))
    CIV = True
    j = 0
    configs = []
    with open('molpro.pun', 'r') as f:
        cV = 0
        cPes = 0
        cNacs1 = 0
        cNacs2 = cNacs1 + 1
        readV = True
        for lines in f:

            if lines.startswith('ATOM'):
                string = lines.strip().split()
                atom = int(string[1]) - 1
                pos[0, atom] = np.double(float(string[4]))
                pos[1, atom] = np.double(float(string[5]))
                pos[2, atom] = np.double(float(string[6]))

            if 'MCSCF STATE ' in lines and readV:
                string = lines.strip().split()
                if string[3] == 'Energy':
                    v_c[cV] = np.double(float(string[4]))
                    cV += 1
                if cV == dyn.nstates:
                    readV = False
            if 'SA-MC GRADIENT' in lines:
                string = lines.strip().split()
                atom = int(string[3].replace(':', '')) - 1
                grad[0, atom, cPes] = np.double(float(string[4]))
                grad[1, atom, cPes] = np.double(float(string[5]))
                grad[2, atom, cPes] = np.double(float(string[6]))
                if atom == geo.natoms - 1:
                    cPes += 1
            if 'SA-MC NACME' in lines:
                string = lines.strip().split()
                atom = int(string[3].replace(':', '')) - 1
                nacmes[0, atom, cNacs1, cNacs2] = np.double(float(string[4]))
                nacmes[1, atom, cNacs1, cNacs2] = np.double(float(string[5]))
                nacmes[2, atom, cNacs1, cNacs2] = np.double(float(string[6]))
                if atom == geo.natoms - 1:
                    if cNacs2 == dyn.nstates - 1:
                        cNacs1 += 1
                    else:
                        cNacs2 += 1
            if lines.startswith(' ') and CIV:
                ff = lines.strip().split()

                for i in range(3):
                    cis[j, i] = float(ff[i + 1])
                configs.append(ff[0])
                j = j + 1

    print(configs)
    total = len(configs)
    print('total CIS ',total)
    oneciv = int(total / 3)
    cis = cis[0:oneciv, :]
    configs = configs[0:oneciv]
    # for i in range(2):
    #     for j in range(geo.natoms):
    #         print(grad[0, j, i], grad[1, j, i], grad[2, j, i])
    return v_c, grad, nacmes, cis, configs


def update_vars(v_c, grad, nacmes, geo):
    dyn = initdyn()
    '''Update the variables and store them using vectorization'''
    pes = np.zeros(dyn.nstates, dtype=np.double)
    grads = np.zeros((geo.ndf, dyn.nstates), dtype=np.double)
    nacs = np.zeros((geo.ndf, dyn.nstates, dyn.nstates), dtype=np.double)
    for i in range(
            dyn.nstates):  # Updates pes energy with the reference value (value corresponding to the initial geometry GS)
        pes[i] = v_c[i] - dyn.e_ref

    for n in range(dyn.nstates):
        idf = 0
        for i in range(geo.natoms):
            for j in range(3):
                grads[idf, n] = (grad[j, i, n])
                idf += 1

    for n in range(dyn.nstates):
        for k in range(n + 1, dyn.nstates):
            idf = 0
            for i in range(geo.natoms):
                for j in range(3):
                    nacs[idf, n, k] = nacmes[j, i, n, k]  # / np.sqrt(geo.masses[i])
                    nacs[idf, k, n] = -nacmes[j, i, n, k]  # /np.sqrt(geo.masses[i])
                    idf += 1
    '''Gradients are negative as we use them as forces'''
    return pes, -grads, nacs
