#!/usr/bin/env python
# must be RHF ref and spin 0 at the moment
# input InputGeom.txt  
# output SCIgradient_stateI.txt, NAC_stateIstateJ.txt  (J>I) and states must be singlets (spin<0.1). All state spins (FinalSpinValues) and all state energies (FinalE)  
# modifies SHCIconfig.json to create config.json with new values for states, electrons and symmetry but not cutoffs
# Edit 315 to 333 for different molecules and installations

import os
import time
import pyscf
from pyscf import gto,symm,numpy

def mySelectedCI_Energy(x):
	global createIntsLocation
	global SHCIcommand
	global n_double,states

	mol.build() #this is to get a geometry array mol._atom which we will then populate using x vector and rebuild molecule
	
	newmol=mol._atom
	

	for i in range(mol.natm):
		for j in range(3):
			newmol[i][1][j]=x[3*i+j]
		
		
			
	mol.atom=newmol
	#print(mol.atom)
	mol.build()
	
	nbft=mol.nao_nr()
	
	myfile=open("Overlap.txt","w")
	
	Enuc=mol.energy_nuc()
	#print('Enuc=',Enuc)
	
	overlap = mol.intor('int1e_ovlp_sph',aosym='s1')
	
	
	
	print(nbft,file=myfile)
	print(Enuc,file=myfile)
	
	#print('AO overlap')
	for i in range(nbft):
		for j in range(i+1):
			print(overlap[i][j],i+1,j+1,file=myfile)
		
	
	myfile.close()

	
	myfile=open("OneEints.txt","w")
	#print('1e ints')
	kin = mol.intor('int1e_kin_sph',aosym='s1')
	vnuc = mol.intor('int1e_nuc_sph',aosym='s1')
	
	for i in range(nbft):
		for j in range(i+1):
			print(kin[i][j]+vnuc[i][j],i+1,j+1,file=myfile)
	
	myfile.close
	
	myfile=open("TwoEints.txt","w")
	#modify later so that it uses full permutation symmetry see write out of MOs in my HF program
	#print('2e ints')
	eri = mol.intor('int2e_sph',aosym='s1')
	
	for i in range(nbft):
		for j in range(i+1):
			for k in range(nbft):
				for l in range(k+1):
					print(eri[i][j][k][l],i+1,j+1,k+1,l+1,file=myfile)
	
	myfile.close
	
	#write out derivative AO integrals here also
	
	overlapGrad = mol.intor('int1e_ipovlp_sph',aosym='s1')
	
	myfile=open("OverlapGrad.txt","w")
	
	#also need to write number of atoms and AOs range for each atom  
	print(mol.natm,file=myfile)
	aoslices = mol.aoslice_by_atom()
	
	for i in range(mol.natm):
		print(aoslices[i,0],aoslices[i,1],aoslices[i,2],aoslices[i,3],file=myfile) # aoslices[i,3]-1 is end basis function
	
	#print('AO x,y,z  overlap grad')
	for k in range(3):
		for i in range(nbft):
			for j in range(nbft): #not symmetric for grad of overlap
				if abs(overlapGrad[k][i][j])>0.0:   #just write non zero values
					print(overlapGrad[k][i][j],k+1,i+1,j+1,file=myfile)
		
	myfile.close()
	
	#1e gradients
	myfile=open("OneEGrad.txt","w")
	KinGrad = mol.intor('int1e_ipkin_sph',aosym='s1')
	vnucGrad= mol.intor('int1e_ipnuc_sph',aosym='s1')
	for k in range(3):
		for i in range(nbft):
			for j in range(nbft):
					if abs(KinGrad[k][i][j]+vnucGrad[k][i][j])>0.0:   #just write non zero values
						print(KinGrad[k][i][j]+vnucGrad[k][i][j],k+1,i+1,j+1,file=myfile)
	
	
	
	myfile.close()
	
	#<\nabla|1/r|> for each atom so we can use to calculate < |  (r-R)/|r-R|^3 | > contribution from derivative of H^{core}
	myfile=open("rinvGrad.txt","w")
	for myatom in range(mol.natm):
		with mol.with_rinv_at_nucleus(myatom):
			rinvGrad = mol.intor('int1e_iprinv_sph',aosym='s1')
			for k in range(3):
				for i in range(nbft):
					for j in range(nbft): 
						if abs(rinvGrad[k][i][j])>0.0:   #just write non zero values
							print(rinvGrad[k][i][j],myatom+1,k+1,i+1,j+1,file=myfile)
	
		
	myfile.close()
	
	#charges and geometry for nuclear derivatives
	myfile=open("geometry.txt","w")
	for i in range(mol.natm):
		print(mol.atom_charge(i),file=myfile)
		mycoord=mol.atom_coord(i)
		print(mycoord[0],mycoord[1],mycoord[2],file=myfile)
	myfile.close()
	
	
	#2e gradients
	myfile=open("TwoEGrad.txt","w")
	eriGrad = mol.intor('int2e_ip1_sph',aosym='s1')
	
	for coord in range(3):
		for i in range(nbft):
			for j in range(nbft):
				for k in range(nbft):
					for l in range(nbft):
						if abs(eriGrad[coord][i][j][k][l])>0.0:   #just write non zero values
							print(eriGrad[coord][i][j][k][l],coord+1,i+1,j+1,k+1,l+1,file=myfile)
	
	
	myfile.close()

	


	
	
	myhf = mol.HF()
	myhf.kernel()
	
	
	orbsym = symm.label_orb_symm(myhf.mol, myhf.mol.irrep_id, myhf.mol.symm_orb, myhf.mo_coeff)
	
	#print(orbsym)
	
	
	myfile=open("MOenergies.txt","w")
	for i in range(nbft):
	  print(myhf.mo_energy[i],orbsym[i],i+1,file=myfile)
	
	myfile.close()
	
	myfile=open("MOcoeffs.txt","w")
	for i in range(nbft):
	  for j in range(nbft):
	  	print(myhf.mo_coeff[i][j],i+1,j+1,file=myfile)
	
	myfile.close()
	
	#os.system('cp Ints_selectedCIgradIn.txt selectedCIgradIn.txt')
	myfile=open("selectedCIgradIn.txt","w")
	print('double occupied',file=myfile)
	print(n_double,file=myfile)
	print('write MO integrals',file=myfile)
	print('.TRUE.',file=myfile)
	print('FCI gradient',file=myfile)
	print('.FALSE.',file=myfile)
	print('Selected CI gradient',file=myfile)
	print('.FALSE.',file=myfile)
	print('Selected CI NACs',file=myfile)
	print('.FALSE.',file=myfile)
	print('States',file=myfile)
	print(states,file=myfile)
	myfile.close()
	
	t1=time.perf_counter()	
	os.system(createIntsLocation) #create FCIDUMP
	t2=time.perf_counter()
	print("Time in create ints",t2-t1)
	
	#need to get point group and if infinite group choose largest abelian
	
	my_pointgroup=mol.groupname
	if (my_pointgroup=='Coov'):
		my_pointgroup='C2v'	
	
	if (my_pointgroup=='Doov'):
		my_pointgroup='D2h'	
	
	
	#print(my_pointgroup)
	myTemplate=open("SHCIconfig.json","r")
	
	myfile=open("config.json","w")
	
	
	
	for line in myTemplate:
		if "\"n_up\"" in line:
			line="\"n_up\": "+str(n_double)+',\n'
		if "\"n_dn\"" in line:
			line="\"n_dn\": "+str(n_double)+',\n'
		if "\"n_states\"" in line:
			line="\"n_states\": "+str(states)+',\n'
		if "\"point_group\"" in line:
			line="\"point_group\": \""+my_pointgroup+'\"\n'	
			
		myfile.write(line)
	

	
	myfile.close()
	myTemplate.close()
	

	
	#end of testing for writing orbs and restart to mcci.in

	t1=time.perf_counter()
	os.system(SHCIcommand) # run to give civ_out files and energy 
	t2=time.perf_counter()
	print("Time in SHCI",t2-t1)
	print(SHCIcommand)	
	

	
	
	
	if os.path.isfile('FinalSpins'):  #Check calculation has finished and delete wf_*.dat
		print("Selected CI finished correctly")
		os.system('mv FinalSpins FinalSpinValues')
		os.system('rm wf_*.dat')
	else:
		print("Selected CI did not finish")
		quit()
	
	
	return 


def mySCI_Grad_Nac(x): 
	global 	SCIgradLocation
	global NACcommand
	global n_double,states
	
	#os.system('cp Grads_selectedCIgradIn.txt selectedCIgradIn.txt')
	myfile=open("selectedCIgradIn.txt","w")
	print('double occupied',file=myfile)
	print(n_double,file=myfile)
	print('write MO integrals',file=myfile)
	print('.FALSE.',file=myfile)
	print('FCI gradient',file=myfile)
	print('.FALSE.',file=myfile)
	print('Selected CI gradient',file=myfile)
	print('.TRUE.',file=myfile)
	print('Selected CI NACs',file=myfile)
	print('.FALSE.',file=myfile)
	print('states',file=myfile)
	print(states,file=myfile)
	myfile.close()


	t1=time.perf_counter()
	os.system(SCIgradLocation)
	t2=time.perf_counter()
	print("Time in grad",t2-t1)
	
	
	#os.system('cp NACs_selectedCIgradIn.txt selectedCIgradIn.txt')
	myfile=open("selectedCIgradIn.txt","w")
	print('double occupied',file=myfile)
	print(n_double,file=myfile)
	print('write MO integrals',file=myfile)
	print('.FALSE.',file=myfile)
	print('FCI gradient',file=myfile)
	print('.FALSE.',file=myfile)
	print('Selected CI gradient',file=myfile)
	print('.FALSE.',file=myfile)
	print('Selected CI NACs',file=myfile)
	print('.TRUE.',file=myfile)
	print('NACstates',file=myfile)
	print(states,file=myfile)
	myfile.close()
	t1=time.perf_counter()
	os.system(NACcommand) 
	t2=time.perf_counter()	
	print("Time in NAC",t2-t1)
	return


		
	
	
#this part  is the only place that needs changing for different molecules
mol = gto.Mole()
mol.atom='''C  -0.05692030  -0.01970668  -1.30651553 ;C 0.00730436  -0.00165497   1.27257832  ; H  0.24302756   1.88072730  -2.08957572 ;H 0.25564585  -1.71591723  -2.23392567   ; H 0.05235700   1.69199743   2.31779601 ; H 0.03973809  -1.60245793   2.40978991  '''
mol.basis='STO-3G'
mol.unit='B' #Use bohr
mol.symmetry = True
#mol.symmetry_subgroup = 'C1'  #uncomment to always use C1

mybin = '../bin/'
SHCIcommand = 'mpirun -n 1 /home/andres/PycharmProjects/PyMCE/PyMCE/ModifiedArrowSHCI/shci_t2rdm/shci'

SCIgradLocation = mybin + 'OutputSHCIGrad_singlet'
createIntsLocation = SCIgradLocation

NACcommand = mybin + 'OutputSHCINac_singlet'

states = 6
#change for different molecules


mol.build()

#print(mol.symmetry)
#print(mol.groupname) #need to move this to SCI energy after HF orbs at read in geom calculated  use if to print C2v if Coov and D2h if Dooh

print(mol.basis)


n_double=mol.nelec[0]  # number of alphas

if(n_double != mol.nelec[1]):
	print('Need even number of electrons')
	quit()



x0=numpy.zeros([mol.natm*3]) #create the 1d array x0


#copy the pyscf coords to the 1d array
for i in range(mol.natm):
	for j in range(3):
		x0[3*i+j]=mol._atom[i][1][j]


myfile=open("InitialGeomFromPython.txt","w")

for i in range(mol.natm):
	print(mol._atom[i][1][0],mol._atom[i][1][1],mol._atom[i][1][2],file=myfile)
	
myfile.close()

myfile=open("InputGeom.txt","r")
x0= numpy.loadtxt(myfile) 
myfile.close()

#x0=x0.flatten()

#print(x0)


mol.unit='B' #as using the mol._atom in myHF_Energy, and myHF_Grad and don't want to convert twice if we used angstrom


mySelectedCI_Energy(x0)


mySCI_Grad_Nac(x0)  






