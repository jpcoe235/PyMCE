import input
import Constants
import traj
import numpy as np

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


def transform_q_to_r(q):
    geo = input.initgeom()
    rkinit_1 = geo.rkinit

    ph = Constants.physconst()
    k = 0

    for i in range(geo.natoms):
        mass_inv = 1.0 / np.sqrt(geo.masses[i])
        for j in range(3):
            r[j, i] = (mass_inv * q[k] * rkinit_1[k]) * ph.bohr
            k += 1

    return r


def create_input(trajN, istep, q):
    ab = ab_par()
    geo = input.initgeom()
    dyn=input.initdyn()
    '''routine to create molpro input, valid for molpro2012'''

    file = 'molpro_traj_' + str(trajN) + '_' + str(istep) + '.inp'
    r =transform_q_to_r(q)
    with open(file,'w') as f:
        f.write('***'+ab.molecule+' '+'calculation of '+str(istep)+' step in trejectory '+str(trajN)+'\n')
        if ab.first:
            line_wr='file,2,molpro.wfu,new\n'
        else:
            line_wr = 'file,2,molpro.wfu\n'
        f.write(line_wr)
        f.write('punch,molpro.pun\n')
        f.write('basis='+ab.basis+'\n')
        f.write('''symmetry,nosym;
        orient,noorient;
        angstrom;
        geomtype=xyz;
        geom={
        ''')
        f.write(str(geo.natoms)+'\n')
        f.write('\n')
        for i in range(geo.natoms):
            line_wr=geo.atnames[i]+' '+str(r[0,i])+' '+str(r[1,i])+' '+str(r[2,i])+'\n'
            f.write(line_wr)
        f.write('}\n')
        f.write('''{multi,failsafe;
        maxiter,40;
        ''')
        line_wr='occ,'+str(ab.act_orb)+';\n'
        line_wr2='closed,'+str(ab.closed_orb)+';\n'
        line_wr3='wf.'+str(ab.n_el)+',1,1;'
        line_wr4='state,'+str(dyn.nstates)+';\n'

        f.write(line_wr)
        f.write(line_wr2)
        f.write(line_wr3)
        f.write(line_wr4)
        f.write('''weight,1,1,1;
        orbital,2140.3;
        canonical,2140.2,ci;
        ORBITAL,IGNORE_ERROR;
        ''')







