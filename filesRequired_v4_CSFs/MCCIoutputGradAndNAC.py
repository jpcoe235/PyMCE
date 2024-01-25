#!/usr/bin/env python
# must be RHF ref and spin 0 at the moment
# input InputGeom.txt  whichstate.txt
# output SCIgradient.txt, NAC_stateIstateJ.txt  (J>I one is state of interest)  and FinalE
#now writes orbital labels to mcci.in and number of alpha etc. Also modifies restart flag if civ_in exists.
#also write input for CSFtoSD

import os
import pyscf
from pyscf import gto,symm,numpy

def mySelectedCI_Energy(x):
	global createIntsLocation
	global twoRDMcommand
	global MCCIcommand
	global CSFtoSDcommand
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
	myfile.close()
	
	
	os.system(createIntsLocation) #create FCIDUMP
	
	#use orbsym to get orblabels for mcci.in
	
	mysymcount=numpy.zeros([8]) #create an array to count number of each irrrep
	for i in range(nbft):
		mysymcount[orbsym[i]]=mysymcount[orbsym[i]]+1
	
	
	
	mysumsym=numpy.zeros([8]) 
	for i in range(7):
		mysumsym[i+1]=mysumsym[i]+mysymcount[i]
	
	

	
	mysymcount=numpy.zeros([8]) # zero this again
	mydoubleocc=numpy.zeros([n_double]) 
	for i in range(n_double):
		mydoubleocc[i]=mysymcount[orbsym[i]]+mysumsym[orbsym[i]]+1
		mysymcount[orbsym[i]]=mysymcount[orbsym[i]]+1
	
	
	mysorted_doubleocc=(sorted(mydoubleocc))
	
	my_doubleoccString=""
	for i in range(n_double-1):
		my_doubleoccString=my_doubleoccString+str(int(mysorted_doubleocc[i]))+","
	my_doubleoccString=my_doubleoccString+str(int(mysorted_doubleocc[n_double-1]))
	
		
	
	myTemplate=open("mcciCalc.in","r")
	
	myfile=open("mcci.in","w")
	
	myrestartFlag=".false."
	if os.path.isfile('civ_in'):
		myrestartFlag=".true."
	
	for line in myTemplate:
		if "restart" in line:
			line="restart = "+myrestartFlag+'\n'
		if "n_up" in line:
			line="n_up= "+str(n_double)+'\n'
		if "mo_up" in line:
			line="mo_up= "+my_doubleoccString+'\n'
		if "n_dn" in line:
			line="n_dn= "+str(n_double)+'\n'
		if "mo_dn" in line:
			line="mo_dn= "+my_doubleoccString+'\n'
		if "ieig" in line:
			line="ieig = "+str(states)+'\n'
			
		myfile.write(line)
	

	
	myfile.close()
	myTemplate.close()
	
	#end of testing for writing orbs and restart to mcci.in

	
	os.system(MCCIcommand) # run to give civ_out files and energy 
	
	if os.path.isfile('done.txt'):
		print("Selected CI finished correctly")
		os.system('rm done.txt')
		
		#create CSFtoSDinput.txt
		myfile=open("CSFtoSDinput.txt","w")
		print('nbft',file=myfile)
		print(nbft,file=myfile)
		print('ntotal',file=myfile)
		print(2*n_double,file=myfile)
		print('n_beta',file=myfile)
		print(n_double,file=myfile)
		myfile.close()
		
		#then loop over states to convert to SDs
		for i in range(states):
			myCommandString="cp civ_out_state"+str(states-i)+" civ_in"
			#print(myCommandString)
			os.system(myCommandString) #convert state  to CSFs
			os.system(CSFtoSDcommand) # gives civ_out_SDs
			myCommandString="cp civ_out_SDs civ_out_state"+str(states-i)
			#print(myCommandString)
			os.system(myCommandString)
	else:
		print("Selected CI did not finish")
		quit()


	myTemplate=open("mcciTRDM.in","r")
	
	myfile=open("mcci.in","w")
	
	
	
	for line in myTemplate:
		if "n_up" in line:
			line="n_up= "+str(n_double)+'\n'
		if "mo_up" in line:
			line="mo_up= "+my_doubleoccString+'\n'
		if "n_dn" in line:
			line="n_dn= "+str(n_double)+'\n'
		if "mo_dn" in line:
			line="mo_dn= "+my_doubleoccString+'\n'
		if "ieig" in line:
			line="ieig = "+str(states)+'\n'
			
		myfile.write(line)
	

	
	myfile.close()
	myTemplate.close()

	#os.system('cp mcciTRDM.in mcci.in')
	
	os.system(twoRDMcommand) # run to get 2rdm for whichstate.txt and transition 2rdms
	
	
	return 


def mySCI_Grad_Nac(x): 
	global 	SCIgradLocation
	global NACcommand
	global twoRDMcommand
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
	myfile.close()
	os.system(twoRDMcommand)
	os.system(SCIgradLocation)

	
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
	os.system(NACcommand) 
	
	return


		
	
	
#this part  is the only place that needs changing for different molecules
mol = gto.Mole()
mol.atom='''C  1.269698160   -0.000000000    0.000000000 ; C, -1.269698160   -0.000000000    0.000000000; H 2.337457645    1.745780932    0.000000000
             ;H -2.337457645    1.745780932    0.000000000; H 2.337457645    -1.745780932    0.000000000; H -2.337457645    -1.745780932    0.000000000 '''  # what if we change system to O2 then it still works as it overwrites mcci.in etc with new orb labels etc

mol.basis='sto-3g'
mol.unit='B' #Use bohr
mol.symmetry = True
mol.symmetry_subgroup = 'C1'  #uncomment to always use C1

mybin='../bin/'


SCIgradLocation=mybin+'OutputSCIgrad'
createIntsLocation=SCIgradLocation
MCCIcommand='mpirun -n 1 '+mybin+'mcci/mcci' 
twoRDMcommand='mpirun -n 1 '+mybin+'transition2RDMuseWhichstate/mcci'
NACcommand=mybin+'Output_NAC'   # this only calculates NAC between state in whichstate.txt and other states
CSFtoSDcommand=mybin+'CSFsToSDs_UpTo96orbs'

states=3
#change for different molecules


mol.build()


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
for i in range(3*mol.natm):
	print(x0[i],file=myfile)
	
myfile.close()

myfile=open("InputGeom.txt","r")
x0= numpy.loadtxt(myfile) 
myfile.close()

#x0=x0.flatten()

print(x0)


mol.unit='B' #as using the mol._atom in myHF_Energy, and myHF_Grad and don't want to convert twice if we used angstrom


mySelectedCI_Energy(x0)
mySCI_Grad_Nac(x0)






