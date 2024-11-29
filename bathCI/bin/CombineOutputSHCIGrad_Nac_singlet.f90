
! Analytic gradients and non-adiabatic couplings for selected configuration interaction.

! Jeremy P. Coe  (2024)

! Please cite the following papers in any publications resulting from this code:

! [1] Analytic Gradients for Selected Configuration Interaction. J. P. Coe, J. Chem. Theory Comput. 19, 874 (2023).
! [2] Analytic Non-adiabatic Couplings for Selected Configuration Interaction via Approximate Degenerate Coupled Perturbed Hartree–Fock. J. P. Coe, J. Chem. Theory Comput. 19, 8053 (2023)



! This version only prints singlet states labelled in order of singlets not in order of total states

! Requires blas and lapack.





module CI_NAC_arrays
double precision, allocatable ::  TtwoRDM(:,:,:,:)  ! transition 2RDM
double precision, allocatable ::  ToneRDM(:,:)   !transition 1RDM
double precision, allocatable ::  backT1RDM(:,:)
double precision   ,allocatable::   NACs(:,:)   !(1,2) is x on atom 2
end module CI_NAC_arrays


module CIgrad_arrays
double precision, allocatable ::  TwoRDM(:,:,:,:)
double precision, allocatable ::  OneRDM(:,:)
double precision   ,allocatable     ::   Sgrad_atom(:,:) 
double precision   ,allocatable     ::   O1(:,:) 
double precision   ,allocatable     ::   H1(:,:)
double precision   ,allocatable     ::   Hgrad_atom(:,:)  
double precision   ,allocatable     ::   TwoEgrad_atom(:,:,:,:)  
double precision   ,allocatable     ::   TwoEgrad_MO(:,:,:,:)  
double precision   ,allocatable     ::   X(:,:)  
double precision   ,allocatable     ::   Zvec(:,:)  
double precision deltaXthresh,deltaOrbEthresh
end module CIgrad_arrays

module RHFpolari_arrays
double precision   ,allocatable     ::   dipoleints(:,:,:)
double precision   ,allocatable     ::   MOdipoleints(:,:,:)
double precision   ,allocatable     ::   u1(:,:)
end module RHFpolari_arrays

module RHFgrad_arrays
double precision   ,allocatable     ::   dmE(:,:)  ! Fock eigenvalue weighted density matrix
integer   ,allocatable     ::   firstAO(:) !first AO for each atom 
integer   ,allocatable     ::   lastAO(:) !last AO for each atom 
integer natoms
double precision   ,allocatable     ::   Sgrad(:,:,:) !(1,:,:) is x  (3,:,:) is z   overlap gradient integrals
double precision   ,allocatable     ::   e1grad(:,:,:) !(1,:,:) is x  (3,:,:) is z   one electron gradient integrals
double precision   ,allocatable     ::   gradients(:,:)   !(1,2) is x on atom 2
double precision   ,allocatable     ::   geometry(:,:)   
double precision   ,allocatable     ::   charge(:)   
double precision   ,allocatable     ::   rinvgrad(:,:,:,:)  !  (1,2,:,:) for atom 1 y coord etc 
double precision   ,allocatable     ::   e2grad(:,:,:,:,:) !(1,:,:,:,:) is x  (3,:,:,:,:) is z   two electron gradient integrals
end module RHFgrad_arrays

module AOintegral_arrays
integer nbft
double precision   ,allocatable     ::   overlap(:,:)
double precision   ,allocatable     ::   e1ints(:,:)
double precision   ,allocatable     ::   e2ints(:,:,:,:)
double precision   ,allocatable     ::   SMH(:,:)
double precision   ,allocatable     ::   Fock(:,:)
double precision   ,allocatable     ::   OrthoFock(:,:)
double precision   ,allocatable     ::   eigenfn(:,:)
double precision   ,allocatable     ::   FockE(:)
double precision   ,allocatable     ::   denmat(:,:)
double precision   ,allocatable     ::   DIISerrorMatrix(:,:,:)
double precision   ,allocatable     ::   DIIScoeffs(:)
double precision   ,allocatable     ::   FockHist(:,:,:)
double precision ecore,SeigenTol
integer maxDIIS
double precision , allocatable:: MOe2ints(:,:,:,:)
double precision , allocatable:: MOe1ints(:,:)
integer   ,allocatable     ::   MOsym(:)
end module AOintegral_arrays

PROGRAM main
use AOintegral_arrays
implicit none
integer i,j,k,docc,maxloops,loop,mat_i,nDIIS
double precision energy,denmaterror,DIISconvTol,DIISerror
logical converged,calcMOints,writeMOints,FCIgrad,calcRHFpolar
logical CISDgrad,selectedCIgrad,readinRHF
logical selectedCI_NACs
integer NACstates,NACstate1,NACstate2
CHARACTER (LEN=30)  :: tempchar

docc=1
maxloops=30
DIISconvTol=1.0D-8 ! tolerance for largest value of DIIS error matrix for convergence
SeigenTol=0.0001D0 ! tolerance for eigenvalues of overlap matrix to check for linear dependence
maxDIIS=8
calcMOints=.TRUE.
writeMOints=.FALSE. ! this now writes energy not symmetry ordered FCIDUMP file - better for SHCI but will not work with standard MCCI
FCIgrad=.FALSE.
calcRHFpolar=.FALSE. ! dipole and polarizability for RHF
CISDgrad=.FALSE.
selectedCIgrad=.FALSE.
selectedCI_NACs=.FALSE.
readinRHF=.TRUE. !if true this reads in MO coeffs,energies and symms from pySCF output and does not do my own no sym RHF calc

!READ in settings

OPEN(UNIT=14, FILE='selectedCIgradIn.txt')
READ(14,*) tempchar
READ(14,*) docc
READ(14,*) tempchar
READ (14,*) writeMOints
READ(14,*) tempchar
READ(14,*) FCIgrad
READ(14,*) tempchar
READ(14,*) selectedCIgrad
READ(14,*) tempchar
READ(14,*) selectedCI_NACs
READ(14,*) tempchar
READ(14,*) NACstates
CLOSE(14)



!!!!above variables can be modifed
converged=.FALSE.
loop=1





call ReadinAOintsPySCF() !allocates and reads in overlap,e1ints,e2ints. reads in nbft (basis size), ecore




allocate(SMH(nbft,nbft))
allocate(FOCK(nbft,nbft))
allocate(orthoFOCK(nbft,nbft))
allocate(eigenFn(nbft,nbft))
allocate(denmat(nbft,nbft))
allocate(DIISerrorMatrix(nbft,nbft,maxDIIS))
allocate(DIIScoeffs(maxDIIS))
allocate(FockHist(nbft,nbft,maxDIIS))
allocate(FockE(nbft))
allocate(MOsym(nbft))

!set to no symmetry intially but may be overwritten by ReadinRHFdata with symmetry values
do i=1,nbft
MOsym(i)=1  
end do


if(readinRHF) THEN



call ReadinRHFdata()

! Create denmat as we'll use it later and check energy is right
call createDenmat(docc)
call FormFock()  !Form Fock Matrix from denmat
call HFenergy(energy) 

!PRINT *,'Energy=',energy+ecore
converged=.TRUE. !as we just read in the converged values

ELSE

call createSMH() ! create SMH=overlap^-0.5




!guess FOCK as 1e only  
FOCK=e1ints




call TRAN2orthoFOCK() !transform FOCK matrix to orthoFOCK= SMH^dagger FOCK SMH



call DIAGorthoFOCK(converged) !diagonalize orthoFOCK matrix and transform eigenfunctions back to non orthogonal AO basis . Prints eigenvalues if converged is true




do loop=1,maxloops



call createDenmat(docc) ! given number of double occupied docc create density matrix by setting k to docc in dgemm





call FormFock()  !Form Fock Matrix from denmat



call HFenergy(energy) !Calc HF energy as 0.5D0*Trace[(e1ints+F)*D]  

PRINT *, loop,'Energy=',energy+ecore



mat_i=mod(loop-1,maxDIIS)+1
call CalcDIISerrorMatrix(mat_i)  

call GetMaxDIISErrorMatrix(DIISerror,mat_i)

PRINT *, 'DIIS error matrix max value', DIISerror




if(DIISerror.lt.DIISconvTol) THEN
converged=.TRUE.
EXIT
END IF





!Add fock matrix to the history
do j=1,nbft
do i=1,nbft
FockHist(i,j,mat_i)=Fock(i,j)
end do
end do



!if we have more than two error matrices stored then calculate DIIS coeffs
if(loop.gt.1) THEN ! have at least 2 DIISerrorMatrices so run DIIS update on Fock matrix
nDIIS=min(loop,maxDIIS) !we only have maxDIIS error Matrices stored
call CalcDIIScoeffs(nDIIS) 

!use FockHist and DIIScoeffs to replace Fock matrix w
Fock=0.0D0
do k=1,nDIIS
do j=1,nbft
do i=1,nbft
Fock(i,j)=Fock(i,j)+DIIScoeffs(k)*FockHist(i,j,k)
end do
end do
end do

END IF






call TRAN2orthoFocK() !transform FOCK matrix to orthoFOCK= SMH^dagger FOCK SMH
call DIAGorthoFocK(converged) !diagonalize orthoFOCK matrix and transform eigenfunctions back to non orthogonal AO basis . Prints eigenvalues if converged is TRUE


end do

END IF !end of  ELSE part of if(readinRHF)

if(converged) THEN

if(.NOT.readinRHF) then !only do this if we did not just read in RHF orb energies, syms and coeffs
PRINT *, 'Converged to',DIISconvTol,'as max DIIS error matrix is',DIISerror
! 1 last calculation using actual FOCK matrix not DIIS approx
call TRAN2orthoFocK() !transform FOCK matrix to orthoFOCK= SMH^dagger FOCK SMH
call DIAGorthoFocK(converged) !diagonalize orthoFOCK matrix and transform eigenfunctions back to non orthogonal AO basis . Prints eigenvalues if converged is TRUE
end if



call RHFgradients(docc)




! call calcMOoverlap() ! calcs and prints MO overlap matrix for testing that this is identity
if(calcMOints.OR.calcRHFpolar) then
call calcMOoneElec(writeMOints)  ! one electron integrals in MO basis  check that these are the same as Columbus  - yes to 9 decimal places
call calcMOtwoElec(writeMOints,docc)

!PRINT *, 'MO integrals stored in MOe1ints and MOe2ints'

if(calcRHFpolar) call RHFpolarizability(docc)  

!!!!!!!!!assuming 2RDM has been calculated and stored 
if(FCIGrad) call FCIgradients(docc) ! energy invariant to any orbital rotation if no frozen so don't need CPHF

!Call two routines to make sure Zvec approach reproduces original CPHF
if(CISDgrad) call CISDgradients(docc)

if(CISDgrad) call CISDgradientsZvec(docc)

if(selectedCIgrad) call selectedCIgradZvec(docc,NACstates)

if(selectedCI_NACs) call selectedCI_NACsZvec(docc,NACstates)  ! uses Z vector approach

!!!!!!!!!!!end of CIgrad 
deallocate(MOe2ints)
deallocate(MOe1ints)
end if
END IF !end of converged 

if(loop.gt.maxloops) PRINT *,'No convergence in',maxloops,'SCF iterations'


DEALLOCATE(denmat)
DEALLOCATE(FOCK)
DEALLOCATE(orthoFOCK)
DEALLOCATE(eigenFn)
DEALLOCATE(e1ints)
DEALLOCATE(e2ints)
DEALLOCATE(overlap)
DEALLOCATE(SMH)
DEALLOCATE(DIISerrorMatrix)
DEALLOCATE(DIIScoeffs)
DEALLOCATE(FockHist)
DEALLOCATE(FockE)
DEALLOCATE(MOsym)

END PROGRAM



subroutine selectedCI_NACsZvec(docc,NACstates)   
use CI_NAC_arrays
use CIgrad_arrays
use AOintegral_arrays, only: nbft,FockE
use RHFgrad_arrays, only:e2grad,gradients,Sgrad,natoms,firstAO,lastAO,&
geometry,charge,rinvgrad,e1grad
implicit none
integer, intent(in) :: docc,NACstates
integer i,coord,atom,j,whichstate
double precision overDeltaE
double precision , allocatable:: stateE(:)
integer NACstateA,NACstateB
CHARACTER (LEN=30)  :: f1 
CHARACTER (LEN=30)  :: f2
double precision , allocatable:: Spins(:)
integer  singletcountA,singletcountB

deltaXthresh=1.0D-8
deltaOrbEthresh=1.0D-6
!deltaOrbEthresh=1.0D-2

PRINT *, 'NACstates',NACstates
allocate(stateE(NACstates))

!Read in state energies
OPEN(UNIT=14,FILE='FinalE')
do i=1,NACstates
READ (14,*) stateE(i)
end do
CLOSE(14)


allocate(Spins(NACstates))

OPEN(UNIT=14,FILE='FinalSpinValues')
do i=1,NACstates
READ (14,*) Spins(i)
end do
PRINT *, 'If state has spin>0.1 then excluding'
CLOSE(14)










 allocate(Zvec(nbft,nbft))  !we use map to 1D to get vector when solving linear equations then map back to 2D for ease of use 
 
 allocate(Sgrad(3,nbft,nbft))

call readinSgrad() !reads in natoms, firstAO array, last AO array and Sgrad array
allocate(NACs(3,natoms))  ! first index is x,y,z

 
allocate(e2grad(3,nbft,nbft,nbft,nbft))
call readinTwoEgrad()  !this may have already been called by RHF gradients and CI gradients but they are self contained so e2grad array has been deallocated
 
 

 
 !need following to get things for oneE
allocate(geometry(3,natoms))
allocate(charge(natoms))

allocate(rinvgrad(natoms,3,nbft,nbft))



call readinGeometry()


allocate(e1grad(3,nbft,nbft))
call readinOneEgrad() 
call readin_rinvgrad() 
 
  allocate(X(nbft,nbft))
 
  allocate(TtwoRDM(nbft,nbft,nbft,nbft))
  allocate(ToneRDM(nbft,nbft)) !calculate from 2RDM

allocate(BackT1RDM(nbft,nbft))
  

 
 OPEN(UNIT=14,FILE='Singlet_FinalE')
 do NACstateA=1,NACstates
  if(ABS(Spins(NACstateA)).gt.0.1D0) cycle
  WRITE(14,*) stateE(NACstateA)
 
 end do
 CLOSE(14)
 
 
 singletcountA=0
 singletcountB=0
 do NACstateA=1,NACstates-1
  if(ABS(Spins(NACstateA)).gt.0.1D0) cycle
  singletcountA=singletcountA+1
  singletcountB=singletcountA
 do NACstateB=NACstateA+1,NACstates
 

 if(ABS(Spins(NACstateB)).gt.0.1D0) cycle
 singletcountB=singletcountB+1
 
 NACs=0.0D0
 call readinT2RDM_Arrow(docc,NACstateA,NACstateB) !   reads in T2RDM between state1 and state2

 call pureAOtwoE_NAC()   ! creates then uses back transform of TtwoRDM  

 call pureAOoneE_NAC() !pure AO derivative 1e contribution also stores back transformation of T1RDM (backT1RDM) as we need it for the "DET" part
 

 

 


 call calcXlagrangian_NAC()
 
 


 
 
 call MO_CIcontributionOverlap_NAC() 
 


    
  
    

 overdeltaE=1.0D0/(stateE(NACstateB)-stateE(NACstateA))
 NACs=overdeltaE*NACs   
 

 




 call DET_OverlapContribution_NAC() 











 
 X=X*overdeltaE ! include scaling by energy 
X=X+0.5D0*ToneRDM  !combine terms from CI and SD z vector parts 





Zvec=0.0D0

 call calcZvecSelectedCI(docc) ! this is using the combined X

allocate(gradients(3,natoms)) 
gradients=0.0D0

call contractZwithBSelectedCIdegen(docc) !this writes to gradients but rest does not need to be changed for NACs so let's just just add its gradient contribution to nacs

 do i=1,natoms
 do j=1,3
 NACs(j,i)=NACs(j,i)+gradients(j,i)
end do
end do
deallocate(gradients)


X=0.5D0*ToneRDM





 call MO_CIcontributionOverlap_NAC()   






write(f1,'(i2)') singletcountA
 write(f2,'(i2)') singletcountB
            f1='state'//trim(adjustl(f1))
            f2='state'//trim(adjustl(f2))
            f1=trim(adjustl(f1))//trim(adjustl(f2))
            f1='NAC_singlets_'//trim(adjustl(f1))
            f2='.txt'
            f1=trim(adjustl(f1))//trim(adjustl(f2))

OPEN(UNIT=16,FILE=f1)
do i=1,natoms

WRITE (16,*) NACs(1,i),NACs(2,i),NACs(3,i)

end do
CLOSE(16)





end do
end do


deallocate(StateE)
deallocate(e2grad)
deallocate(geometry)
deallocate(charge)
deallocate(rinvgrad)
deallocate(e1grad)
deallocate(Zvec)
deallocate(X)
deallocate(Sgrad)  
deallocate(firstAO)
deallocate(lastAO)
deallocate(TtwoRDM)
deallocate(ToneRDM)
deallocate(NACs)
deallocate(backT1RDM)

deallocate(Spins)
end subroutine selectedCI_NACsZvec



subroutine contractZwithBSelectedCIdegen(docc)   
use AOintegral_arrays, only: nbft,MOe2ints,FockE,eigenFn,denmat,&
FockE
use RHFgrad_arrays, only:natoms,firstAO,lastAO,Sgrad,e1grad&
,rinvgrad,charge,e2grad,gradients
use CIgrad_arrays, only:  Sgrad_atom,O1,H1,Hgrad_atom,&
TwoEgrad_atom,TwoEgrad_MO,Zvec,deltaOrbEthresh,X
implicit none
integer, intent(in) :: docc
integer atom,coord
integer j,sigma,k,r ,l! use notation of Gerratt and Mills JCP 49, 1719 (1968) 
integer n,m
integer i,p,q,s
double precision dtemp,deltaOrbE,deltaX
double precision , allocatable:: B(:,:)
double precision , allocatable:: Zback(:,:)
double precision , allocatable:: D(:,:)
integer N_degen
integer , allocatable:: degen_pairs(:,:)


allocate(degen_pairs(nbft,2))
!!!!!!Get degen pairs
!PRINT *, 'Degenerate pairs for threshold',deltaOrbEthresh
N_degen=0
do q=1,nbft-1
do p=q+1,nbft

deltaOrbE=FockE(p)-FockE(q)

if(ABS(deltaOrbE).lt.deltaOrbEthresh) THEN
N_degen=N_degen+1
degen_pairs(N_degen,1)=p
degen_pairs(N_degen,2)=q


END IF

end do
end do

!!!!!!!




allocate(B(nbft,nbft))
allocate(Zback(nbft,nbft))
allocate(D(nbft,nbft))


n=docc
m=nbft



Zback=0.0D0



!use D for temporary storage in back transform first
D=0.0D0


do j=1,nbft-1  
do i=j+1,nbft
do p=1,nbft

D(p,j)=D(p,j)+&
(eigenFn(p,i)*Zvec(i,j))

end do
end do
end do


do j=1,nbft-1 
do q=1,nbft
do p=1,nbft

Zback(p,q)=Zback(p,q)+&
(eigenFn(q,j)*D(p,j))

end do
end do
end do







D=0.5D0*Denmat  



allocate(Sgrad_atom(nbft,nbft))
allocate(O1(nbft,nbft))
allocate(H1(nbft,nbft))
allocate(Hgrad_atom(nbft,nbft))  
allocate(TwoEgrad_atom(nbft,nbft,nbft,nbft))  
  allocate(TwoEgrad_MO(nbft,nbft,nbft,nbft)) 

do atom=1,natoms
do coord=1,3


!Transform deriv integrals to MO. 

Sgrad_atom=0.0D0

do j=1,nbft
do i=firstAO(atom),lastAO(atom) 
Sgrad_atom(i,j)=-1.0D0*Sgrad(coord,i,j) 
end do
end do



do i=1,nbft
do j=firstAO(atom),lastAO(atom)
Sgrad_atom(i,j)=Sgrad_atom(i,j)-1.0D0*Sgrad(coord,j,i) 
end do
end do


!transform Sgrad_atom to MO basis to give O1
call calcO1() 

!now H1   Gerrat and Mills use H for one-e and G for two-e terms




Hgrad_atom=0.0D0
!copy relevant part of e1grad to Hgrad_atom to get one e deriv w.r.t nuclear
do j=1,nbft
do i=firstAO(atom),lastAO(atom) 
Hgrad_atom(i,j)=-1.0D0*e1grad(coord,i,j) 
end do
end do


do i=1,nbft
do j=firstAO(atom),lastAO(atom)
Hgrad_atom(i,j)=Hgrad_atom(i,j)-1.0D0*e1grad(coord,j,i) 
end do
end do

!also need derivative of One e Hamiltonian term
do i=1,nbft
do j=1,nbft ! not orbital deriv so don't limit this
Hgrad_atom(i,j)=Hgrad_atom(i,j)-charge(atom)*rinvgrad(atom,coord,i,j)
Hgrad_atom(i,j)=Hgrad_atom(i,j)-charge(atom)*rinvgrad(atom,coord,j,i) 
end do
end do

! transform Hgrad_atom to MO basis to give H1
call calcH1() 

!now Two electron gradients


TwoEgrad_atom=0.0D0
!
do i=firstAO(atom),lastAO(atom)
do j=1,nbft
do k=1,nbft
do r=1,nbft
TwoEgrad_atom(i,j,k,r)=-1.0D0*e2grad(coord,i,j,k,r) 
end do
end do
end do
end do

!second index
do i=1,nbft
do j=firstAO(atom),lastAO(atom)
do k=1,nbft
do r=1,nbft
TwoEgrad_atom(i,j,k,r)=TwoEgrad_atom(i,j,k,r)-1.0D0*e2grad(coord,j,i,r,k) 
end do
end do
end do
end do


do i=1,nbft
do j=1,nbft
do k=firstAO(atom),lastAO(atom)
do r=1,nbft
TwoEgrad_atom(i,j,k,r)=TwoEgrad_atom(i,j,k,r)-1.0D0*e2grad(coord,k,j,i,r) 
end do
end do
end do
end do

do i=1,nbft
do j=1,nbft
do k=1,nbft
do r=firstAO(atom),lastAO(atom)
TwoEgrad_atom(i,j,k,r)=TwoEgrad_atom(i,j,k,r)-1.0D0*e2grad(coord,r,k,j,i) 
end do
end do
end do
end do









!create B which is all terms independent of u in CPHF


B=0.0D0
do j=1,m-1  ! as selected CI go over all orbital pairs when sigma>j
do sigma=j+1,m

B(sigma,j)=H1(sigma,j)-FockE(j)*O1(sigma,j)


do k=1,n

do l=1,n
B(sigma,j)=B(sigma,j)-O1(k,l)*(2.0D0*MOe2ints(sigma,l,j,k)-&
MOe2ints(sigma,l,k,j))
end do
end do ! end of k loop

end do
end do


!contract with Zvec  
dtemp=0.0D0

do j=1,nbft-1
do i=j+1,nbft ! as selected CI all pairs where i>j

dtemp=dtemp+(Zvec(i,j)*B(i,j))
end do
end do



!add in two electron deriv terms using back transformed Z O(N^4) not O(N^5) each time like transforming AO 2e derivs to MO

do s=1,nbft
do r=1,nbft
do q=1,nbft
do p=1,nbft
dtemp=dtemp+(TwoEgrad_atom(p,q,r,s)*&
(2.0D0*Zback(p,r)*D(q,s)-Zback(p,s)*D(q,r)))

end do 
end do
end do
end do


gradients(coord,atom)=gradients(coord,atom)+2.0D0*dtemp 


!!!!!!!!!!!!!!!degeneracy contribution as we have O1 in MO basis here


dtemp=0.0D0
do i=1,N_degen
p=degen_pairs(i,1)
q=degen_pairs(i,2)
deltaX=X(p,q)-X(q,p)
dtemp=dtemp-deltaX*O1(q,p)
end do

gradients(coord,atom)=gradients(coord,atom)+dtemp 

!!!!!!!!!!!!!!!!!!end of degenerate contribution


end do !loop over coords
end do !loop over atoms


deallocate(Sgrad_atom)
deallocate(O1)
deallocate(H1)
deallocate(Hgrad_atom)  
deallocate(TwoEgrad_atom)
deallocate(TwoEgrad_MO)
deallocate(B)
deallocate(Zback)
deallocate(D)

deallocate(degen_pairs)
end subroutine contractZwithBSelectedCIdegen





subroutine DET_OverlapContribution_NAC() 
use AOintegral_arrays, only: nbft
use CI_NAC_arrays, only: backT1RDM,NACs
use RHFgrad_arrays, only:natoms,Sgrad,firstAO,lastAO
implicit none
integer i,j
double precision dtemp
integer coord,atom




!overlap gradients
do atom=1,natoms
do coord=1,3

dtemp=0.0D0 
do i=1,nbft   
do j=firstAO(atom),lastAO(atom) 

dtemp=dtemp-(backT1RDM(i,j)*Sgrad(coord,j,i)) 

end do
end do


NACs(coord,atom)=NACs(coord,atom)+dtemp
end do
end do




end subroutine DET_OverlapContribution_NAC









subroutine MO_CIcontributionOverlap_NAC()  
use AOintegral_arrays, only: nbft,MOe2ints,MOe1ints,eigenFn
use CI_NAC_arrays, only: NACs
use CIgrad_arrays, only: X
use RHFgrad_arrays, only:natoms,Sgrad,firstAO,lastAO
implicit none
double precision , allocatable:: backX(:,:)
double precision , allocatable:: Tran1(:,:)
integer t,p,r,s,q,i,j
double precision dtemp
integer coord,atom


allocate(backX(nbft,nbft))
allocate(Tran1(nbft,nbft))





!Back transform X noting that now we have Xij when i.le. j  and Xji otherwise 
BackX=0.0D0
Tran1=0.0D0


do p=1,nbft
do j=1,nbft
do i=1,nbft
if(i.le.j) THEN
Tran1(p,j)=Tran1(p,j)+&
(eigenFn(p,i)*X(i,j))
ELSE
Tran1(p,j)=Tran1(p,j)+&
(eigenFn(p,i)*X(j,i))
END IF

end do
end do
end do

do p=1,nbft
do q=1,nbft
do j=1,nbft
BackX(p,q)=BackX(p,q)+&
(eigenFn(q,j)*Tran1(p,j))
end do
end do
end do




!!!!!!!!!!!!!!!!!!!!!!!



!overlap gradients
do atom=1,natoms
do coord=1,3

dtemp=0.0D0
do j=1,nbft
do i=firstAO(atom),lastAO(atom) 

dtemp=dtemp+(backX(i,j)*Sgrad(coord,i,j)) 

end do
end do

! derivative on ket terms
do i=1,nbft
do j=firstAO(atom),lastAO(atom) 

dtemp=dtemp+(backX(i,j)*Sgrad(coord,j,i))

end do
end do


NACs(coord,atom)=NACs(coord,atom)+dtemp
end do
end do



deallocate(backX)
deallocate(Tran1)
end subroutine MO_CIcontributionOverlap_NAC




subroutine calcXlagrangian_NAC() 
use AOintegral_arrays, only: nbft,MOe2ints,MOe1ints
use CI_NAC_arrays, only: ToneRDM,TtwoRDM
use CIgrad_arrays, only:X
implicit none
integer t,p,r,s,q
double precision dtemp




X=0.0D0


!Store X
do t=1,nbft
do p=1,nbft

dtemp=0.0D0
do q=1,nbft
dtemp=dtemp+(0.5D0*(ToneRDM(p,q)+&
ToneRDM(q,p))*MOe1ints(t,q))
end do

do q=1,nbft
do s=1,nbft
do r=1,nbft
dtemp=dtemp+(0.5D0*(TtwoRDM(p,r,s,q)+TtwoRDM(q,r,s,p))*& 
MOe2ints(t,r,q,s))
end do
end do
end do

X(t,p)=dtemp
end do
end do


end subroutine calcXlagrangian_NAC




subroutine pureAOoneE_NAC()
use CI_NAC_arrays, only: ToneRDM,NACs,backT1RDM
use AOintegral_arrays, only: nbft,eigenFn
use RHFgrad_arrays, only:natoms,e1grad,firstAO,lastAO,&
charge,rinvgrad
implicit none
double precision , allocatable:: Tran1(:,:)
integer p,q,i,j
integer atom,coord
double precision dtemp,dtemp2







allocate(Tran1(nbft,nbft))




BackT1RDM=0.0D0
Tran1=0.0D0


do p=1,nbft
do j=1,nbft
do i=1,nbft
Tran1(p,j)=Tran1(p,j)+&
(eigenFn(p,i)*TOneRDM(i,j))
end do
end do
end do

do p=1,nbft
do q=1,nbft
do j=1,nbft
BackT1RDM(p,q)=BackT1RDM(p,q)+&
(eigenFn(q,j)*Tran1(p,j))
end do
end do
end do


do atom=1,natoms
do coord=1,3

dtemp=0.0D0
do j=1,nbft
do i=firstAO(atom),lastAO(atom) 

dtemp=dtemp-(e1grad(coord,i,j)*(BackT1RDM(i,j)+&
BackT1RDM(j,i))) 
end do
end do



!operator deriv
dtemp2=0.0D0
do j=1,nbft 
do i=1,nbft

dtemp2=dtemp2+(rinvgrad(atom,coord,i,j)*(BackT1RDM(i,j)+&
BackT1RDM(j,i)))
end do
end do



dtemp=dtemp-charge(atom)*dtemp2

NACs(coord,atom)=NACs(coord,atom)+dtemp
end do
end do


deallocate(Tran1)

end subroutine pureAOoneE_NAC





subroutine pureAOtwoE_NAC()  
use CI_NAC_arrays, only: TtwoRDM,NACs
use RHFgrad_arrays, only:natoms,e2grad,firstAO,lastAO
use AOintegral_arrays, only: nbft,eigenFn
implicit none
double precision , allocatable:: Tran1(:,:,:,:)
integer p,q,r,s,i,j,k,m
double precision dtemp
integer atom,coord
double precision, allocatable ::  backT2RDM(:,:,:,:)

allocate(backT2RDM(nbft,nbft,nbft,nbft))

allocate(Tran1(nbft,nbft,nbft,nbft))








Tran1=0.0D0
do m=1,nbft
do k=1,nbft
do j=1,nbft
do p=1,nbft
do i=1,nbft
Tran1(p,j,k,m)=Tran1(p,j,k,m)+&
(eigenFn(p,i)*TtwoRDM(i,j,k,m))
end do
end do
end do
end do 
end do

BackT2RDM=0.0D0 !this is used as temporary storage for transformation at this stage
do m=1,nbft
do k=1,nbft
do q=1,nbft
do j=1,nbft
do p=1,nbft
BackT2RDM(p,q,k,m)=BackT2RDM(p,q,k,m)+&
(eigenFn(q,j)*Tran1(p,j,k,m))
end do
end do
end do
end do 
end do


Tran1=0.0D0 !reuse this storage space
do m=1,nbft
do r=1,nbft
do k=1,nbft
do q=1,nbft
do p=1,nbft
Tran1(p,q,r,m)=Tran1(p,q,r,m)+&
(eigenFn(r,k)*BackT2RDM(p,q,k,m)) 
end do
end do
end do
end do 
end do

!final step to get back2RDM in Tran2

backT2RDM=0.0D0 ! reuse this space for final values
do s=1,nbft
do m=1,nbft
do r=1,nbft
do q=1,nbft
do p=1,nbft
backT2RDM(p,q,r,s)=backT2RDM(p,q,r,s)+&
(eigenFn(s,m)*Tran1(p,q,r,m))
end do
end do
end do
end do 
end do



do atom=1,natoms
do coord=1,3



dtemp=0.0D0
do s=1,nbft
do q=1,nbft
do r=1,nbft
do p=firstAO(atom),lastAO(atom)  
dtemp=dtemp+((backT2RDM(p,r,s,q)+&
backT2RDM(q,r,s,p))*e2grad(coord,p,r,q,s))  
end do
end do
end do
end do



NACs(coord,atom)=NACs(coord,atom)-1.0D0*dtemp

end do
end do





deallocate(Tran1)
deallocate(BackT2RDM)

end subroutine pureAOtwoE_NAC


subroutine readinT2RDM_Arrow(docc,NACstate1,NACstate2) 
use AOintegral_arrays, only: nbft
use CI_NAC_arrays, only: TtwoRDM,ToneRDM
implicit none
integer, intent(in) :: docc,NACstate1,NACstate2
integer info,p,q,r,s
integer i,j,k,m
double precision dtemp
CHARACTER (LEN=30)  :: f1 
CHARACTER (LEN=30)  :: f2
integer whichline

ToneRDM=0.0D0
  TtwoRDM=0.0D0




write(f1,'(i2)') NACstate1
 write(f2,'(i2)') NACstate2
            f1='state'//trim(adjustl(f1))
            f2='state'//trim(adjustl(f2))//'.txt'
            f1=trim(adjustl(f1))//trim(adjustl(f2))
            f1='Tran2RDM_'//trim(adjustl(f1))

PRINT *,'Read in ',f1

OPEN(UNIT=14,FILE=f1)  !now this is from symmetry ordered orbitals and we need to map back to energy order
READ(14,*) i  ! first entry is nbft
info=0
!whichline=1
do while (info==0)  !loop until end of file
READ (14,*,IOSTAT=info)  p,q,r,s,dtemp




i=p+1
j=q+1
k=r+1
m=s+1



TtwoRDM(i,j,m,k)=dtemp ! arrow order

!if(ABS(TtwoRDM(i,j,m,k)).gt.0.1) PRINT *,i,j,k,m,TtwoRDM(i,j,m,k)

end do
CLOSE(14)

!get 1RDM from 2RDM and check we get 1e energy and total energy

do r=1,nbft
do q=1,nbft
do p=1,nbft
ToneRDM(p,q)=ToneRDM(p,q)+TtwoRDM(p,r,r,q)
end do
end do
end do

ToneRDM=ToneRDM/((2.0D0*docc)-1.0D0)








end subroutine readinT2RDM_Arrow












subroutine readinT2RDM(docc,NACstate1,NACstate2) 
use AOintegral_arrays, only: nbft
use CI_NAC_arrays, only: TtwoRDM,ToneRDM
implicit none
integer, intent(in) :: docc,NACstate1,NACstate2
integer info,p,q,r,s
integer i,j,k,m
double precision dtemp
integer , allocatable:: mapToSymOrder(:)
CHARACTER (LEN=30)  :: f1 
CHARACTER (LEN=30)  :: f2


ToneRDM=0.0D0
  TtwoRDM=0.0D0


allocate(mapToSymorder(nbft))
OPEN(UNIT=11,FILE='MapToSymOrder.txt')
do i=1,nbft
READ(11,*) mapToSymOrder(i)
end do
CLOSE(11)

write(f1,'(i2)') NACstate1
 write(f2,'(i2)') NACstate2
            f1='state'//trim(adjustl(f1))
            f2='state'//trim(adjustl(f2))
            f1=trim(adjustl(f1))//trim(adjustl(f2))
            f1='Tran2RDM_'//trim(adjustl(f1))

PRINT *,'Read in',f1

OPEN(UNIT=14,FILE=f1)  !now this is from symmetry ordered orbitals and we need to map back to energy order

info=0
do while (info==0)  !loop until end of file
READ (14,*,IOSTAT=info)  p,q,r,s,dtemp

i=mapToSymOrder(p)
j=mapToSymOrder(q)
k=mapToSymOrder(r)
m=mapToSymOrder(s)

!TtwoRDM(p,r,s,q)=dtemp
TtwoRDM(i,k,m,j)=dtemp
end do
CLOSE(14)

!get 1RDM from 2RDM and check we get 1e energy and total energy

do r=1,nbft
do q=1,nbft
do p=1,nbft
ToneRDM(p,q)=ToneRDM(p,q)+TtwoRDM(p,r,r,q)
end do
end do
end do

ToneRDM=ToneRDM/((2.0D0*docc)-1.0D0)




deallocate(mapToSymorder)



end subroutine readinT2RDM





subroutine ReadinRHFdata() 
use AOintegral_arrays, only: nbft,FockE,eigenFn,MOsym
implicit none
integer i,isym,mo,j,i1,j1
double precision dtemp 

OPEN(UNIT=10,FILE='MOenergies.txt')
do i=1,nbft
READ (10,*) dtemp,isym,mo
FockE(mo)=dtemp
MOsym(mo)=isym+1
end do
CLOSE(10)



OPEN(UNIT=10,FILE='MOcoeffs.txt')
do i=1,nbft
do j=1,nbft

READ(10,*) dtemp,i1,j1
eigenFn(i1,j1)=dtemp
end do
end do

CLOSE(10)


end subroutine ReadinRHFdata




subroutine selectedCIgradZvec(docc,NACstates) !Seems to be working and reproduces CISD result but test with high cutoff NH3 , low cutoff larger basis NH3 as well then CO that should have symmetry but we won't use it
use CIgrad_arrays
use AOintegral_arrays, only: nbft,FockE
use RHFgrad_arrays, only:e2grad,gradients,Sgrad,natoms,firstAO,lastAO,&
geometry,charge,rinvgrad,e1grad
implicit none
integer, intent(in) :: docc
integer, intent(in) :: NACstates
integer i,coord,atom,j,mystate
double precision dtemp
double precision , allocatable:: Spins(:)
CHARACTER (LEN=30)  :: f1 
integer singletcount

deltaXthresh=1.0D-8
deltaOrbEthresh=1.0D-5
!deltaOrbEthresh=1.0D-2


allocate(Spins(NACstates))


OPEN(UNIT=14,FILE='FinalSpinValues')
do i=1,NACstates
READ (14,*) Spins(i)
end do
PRINT *, 'If state has spin>0.1 then excluding'
CLOSE(14)

allocate(Zvec(nbft,nbft))  !we use map to 1D to get vector when solving linear equations then map back to 2D for ease of use 


allocate(Sgrad(3,nbft,nbft))

call readinSgrad() !reads in natoms, firstAO array, last AO array and Sgrad array


allocate(gradients(3,natoms))  ! first index is x,y,z



allocate(e2grad(3,nbft,nbft,nbft,nbft))
call readinTwoEgrad()  !this may have already been called by RHF gradients but that is self contained so e2grad array has been deallocated




!need following to get things for oneE
allocate(geometry(3,natoms))
allocate(charge(natoms))

allocate(rinvgrad(natoms,3,nbft,nbft))



call readinGeometry()


allocate(e1grad(3,nbft,nbft))
call readinOneEgrad() 
call readin_rinvgrad() 

allocate(X(nbft,nbft))


allocate(TwoRDM(nbft,nbft,nbft,nbft))


allocate(OneRDM(nbft,nbft)) !calculate from 2RDM

singletcount=0
do mystate=1,NACstates

if(ABS(Spins(mystate)).gt.0.1D0) cycle
PRINT *, 'state',mystate
singletcount=singletcount+1
PRINT *, 'singlet state',singletcount


call readin2RDM_Arrow(docc,mystate) !  reads in 2rdm and calculates energy using this 2rdm as check -- does not allocate 2rdm now we do that above

gradients=0.0D0

call pureAOtwoE() ! pure AO derivative 2e contribution creates and stores back transformation of 2RDM (back2RDM) rather than transforming all AO derivative integrals to MO basis


call pureAOoneE() !pure AO derivative 1e contribution also creates and store back transformation of 1RDM (back1RDM)





 call Vnngrad()  ! add in nuclear nuclear



call calcXlagrangian()

!PRINT *, 'ABS(delta X) >= deltaXthresh and orbital energy difference for MO pairs'
!do j=1,nbft-1
!do i=j+1,nbft
!if(ABS(X(i,j)-X(j,i)).ge.deltaXthresh) then
!PRINT *,i,j,X(i,j)-X(j,i),FockE(i)-FockE(j)
!end if
!end do
!end do

call MO_CIcontributionOverlap() 

!PRINT *,'Selected CI gradients without CPHF contrib'
!do i=1,natoms
!PRINT *,i,gradients(1,i),gradients(2,i),gradients(3,i) 
!end do

!above this point CISD and selected CI are the same


!call calcZvec(docc)
call calcZvecSelectedCI(docc)  ! Zvec different to CISD case as occ-occ RHF pairs and unocc-unocc RHF pairs might have non zero delta X so need to be included 

!here
call contractZwithBSelectedCI(docc) !this loops over atoms and coords in subroutine


PRINT *,'Selected CI gradients with CPHF contrib using Z vector approach'
do i=1,natoms
PRINT *,i,gradients(1,i),gradients(2,i),gradients(3,i) 
end do

write(f1,'(i2)') singletcount
f1='SCIgradient_Singlet_state'//trim(adjustl(f1))//trim(adjustl('.txt'))
!!!!!!!!!!!!!! write out
OPEN(UNIT=14,FILE=f1)
do i=1,natoms
WRITE (14,*) gradients(1,i),gradients(2,i),gradients(3,i)
end do
CLOSE(14)

end do !loop over states

!!!!!!!!!!!!!!!!!!!!


deallocate(TwoRDM)
deallocate(OneRDM)
deallocate(e2grad)
deallocate(gradients)
deallocate(Sgrad)
deallocate(firstAO)
deallocate(lastAO)
deallocate(geometry)
deallocate(charge)
deallocate(rinvgrad)
deallocate(e1grad)

deallocate(Zvec)
deallocate(X)
deallocate(Spins)
end subroutine SelectedCIgradZvec


subroutine contractZwithBSelectedCI(docc)   ! Only small changes to CISD version to  make sure it goes over all orbital pairs for Zback, B, and Z with B contraction -- done  ! see about using RHF density matrix instead of D and removing factor of 2
use AOintegral_arrays, only: nbft,MOe2ints,FockE,eigenFn,denmat
use RHFgrad_arrays, only:natoms,firstAO,lastAO,Sgrad,e1grad&
,rinvgrad,charge,e2grad,gradients
use CIgrad_arrays, only:  Sgrad_atom,O1,H1,Hgrad_atom,&
TwoEgrad_atom,TwoEgrad_MO,Zvec
implicit none
integer, intent(in) :: docc
integer atom,coord
integer j,sigma,k,r ,l! use notation of Gerratt and Mills JCP 49, 1719 (1968) 
integer n,m
integer i,p,q,s
double precision dtemp
double precision , allocatable:: B(:,:)
double precision , allocatable:: Zback(:,:)
double precision , allocatable:: D(:,:)


allocate(B(nbft,nbft))
allocate(Zback(nbft,nbft))
allocate(D(nbft,nbft))


n=docc
m=nbft



Zback=0.0D0

!make this N^3 instead of N^4 later ! this is only over occ unocc pairs for CISD
!do p=1,nbft
!do j=1,docc
!do i=docc+1,nbft
!!do j=1,nbft-1
!!do i=j+1,nbft
!do q=1,nbft

!Zback(p,q)=Zback(p,q)+&
!(eigenFn(p,i)*eigenFn(q,j)*Zvec(i,j))

!end do
!end do
!end do
!end do

!made O(N^3)
!use D for temporary storage in back transform first
D=0.0D0

!do i=docc+1,nbft
!do j=1,docc
do j=1,nbft-1  ! swap for docc loop etc as using selected CI so all orbital pairs might contribute
do i=j+1,nbft
do p=1,nbft

D(p,j)=D(p,j)+&
(eigenFn(p,i)*Zvec(i,j))

end do
end do
end do

!do j=1,docc
do j=1,nbft-1 !as for selected CI
do q=1,nbft
do p=1,nbft

Zback(p,q)=Zback(p,q)+&
(eigenFn(q,j)*D(p,j))

end do
end do
end do






!also pre contract MO coeffs to allow O(N^4) two electron derivive terms - see my z vector notes page 5 note k only goes to docc=n
!D=0.0D0
!do k=1,n
!do q=1,nbft
!do p=1,nbft


!D(p,q)=D(p,q)+&
!(eigenFn(p,k)*eigenFn(q,k))

!end do
!end do
!end do

D=0.5D0*Denmat  !actually we have already calculated D as it is 0.5*RHF density matrix



allocate(Sgrad_atom(nbft,nbft))
allocate(O1(nbft,nbft))
allocate(H1(nbft,nbft))
allocate(Hgrad_atom(nbft,nbft))  
allocate(TwoEgrad_atom(nbft,nbft,nbft,nbft))  
  allocate(TwoEgrad_MO(nbft,nbft,nbft,nbft)) 

do atom=1,natoms
do coord=1,3


!Transform deriv integrals to MO. Later use back transform of Z for two electron derivs

Sgrad_atom=0.0D0
!copy relevant part of Sgrad to Sgrad atom to get S deriv w.r.t nuclear
do j=1,nbft
do i=firstAO(atom),lastAO(atom) ! only AOs on this atom have non zero
Sgrad_atom(i,j)=-1.0D0*Sgrad(coord,i,j) !sign change as nuclear derivative
end do
end do


! this is only derivative on first index so have to add value with indices swapped as (i, deriv j) = (deriv j, i) for real AOs
do i=1,nbft
do j=firstAO(atom),lastAO(atom)
Sgrad_atom(i,j)=Sgrad_atom(i,j)-1.0D0*Sgrad(coord,j,i) !sign change as nuclear derivative
end do
end do


!transform Sgrad_atom to MO basis to give O1
call calcO1() 

!now H1   Gerrat and Mills use H for one-e and G for two-e terms


!Copy parts on atom of e1grad for orbital derivs

Hgrad_atom=0.0D0
!copy relevant part of e1grad to Hgrad_atom to get one e deriv w.r.t nuclear
do j=1,nbft
do i=firstAO(atom),lastAO(atom) ! only AOs on this atom have non zero
Hgrad_atom(i,j)=-1.0D0*e1grad(coord,i,j) !sign change as nuclear derivative
end do
end do

! this is only derivative on first index so have to add value with indices swapped as (i, deriv j) = (deriv j, i) for real AOs
do i=1,nbft
do j=firstAO(atom),lastAO(atom)
Hgrad_atom(i,j)=Hgrad_atom(i,j)-1.0D0*e1grad(coord,j,i) !sign change as nuclear derivative
end do
end do

!also need derivative of One e Hamiltonian term
do i=1,nbft
do j=1,nbft ! not orbital deriv so don't limit this
Hgrad_atom(i,j)=Hgrad_atom(i,j)-charge(atom)*rinvgrad(atom,coord,i,j)
Hgrad_atom(i,j)=Hgrad_atom(i,j)-charge(atom)*rinvgrad(atom,coord,j,i) ! need deriv on other index see RHF deriv notes
end do
end do

! transform Hgrad_atom to MO basis to give H1
call calcH1() 

!now Two electron gradients


TwoEgrad_atom=0.0D0
!
do i=firstAO(atom),lastAO(atom)
do j=1,nbft
do k=1,nbft
do r=1,nbft
TwoEgrad_atom(i,j,k,r)=-1.0D0*e2grad(coord,i,j,k,r) !e2grad only has derivative on first index after coord
end do
end do
end do
end do

!second index
do i=1,nbft
do j=firstAO(atom),lastAO(atom)
do k=1,nbft
do r=1,nbft
TwoEgrad_atom(i,j,k,r)=TwoEgrad_atom(i,j,k,r)-1.0D0*e2grad(coord,j,i,r,k) ! have to do two swaps to make j first index
end do
end do
end do
end do


do i=1,nbft
do j=1,nbft
do k=firstAO(atom),lastAO(atom)
do r=1,nbft
TwoEgrad_atom(i,j,k,r)=TwoEgrad_atom(i,j,k,r)-1.0D0*e2grad(coord,k,j,i,r) ! only one swap as both x1
end do
end do
end do
end do

do i=1,nbft
do j=1,nbft
do k=1,nbft
do r=firstAO(atom),lastAO(atom)
TwoEgrad_atom(i,j,k,r)=TwoEgrad_atom(i,j,k,r)-1.0D0*e2grad(coord,r,k,j,i) ! two swaps
end do
end do
end do
end do

!transform 2e gradients This will slow it down as 4 index transformation every time costs Nbft^5 and there are 3Natoms calls to solveCPHF- this will be why z vector approach used actually solving linear system will be O(Adim^3) so this could becomr like Nbft^6 as Adim=(Nbft-nelec)*nelec and if nelec approachs Nbft/2
!call calcTwoEGrad_MO()   !we don't need to do this now we instead use Zvector back transfrom to do this with O(Nbft^4)


!now have the H1, O1 , TwoEGrad_MO and already had Moe2ints so put cphf together







!create B which is all terms independent of u in CPHF

!for selected CI this will loop over all but only need unocc occ for CISD
B=0.0D0
!do sigma=n+1,m
!do j=1,n
do j=1,m-1  ! as selected CI go over all orbital pairs when sigma>j
do sigma=j+1,m

B(sigma,j)=H1(sigma,j)-FockE(j)*O1(sigma,j)

!comment out two e_grad terms as we will now inlcude them more efficiently in contraction with B
do k=1,n
!B(sigma,j)=B(sigma,j)+(2.0D0*TwoEgrad_MO(sigma,k,j,k))-&  !note we need to swap middle two compared with Gerrat and Mills CPHF for our notation
!TwoEgrad_MO(sigma,k,k,j)
do l=1,n
B(sigma,j)=B(sigma,j)-O1(k,l)*(2.0D0*MOe2ints(sigma,l,j,k)-&
MOe2ints(sigma,l,k,j))
end do
end do ! end of k loop

end do
end do


!contract with Zvec  !i gt j but only occ unocc then double see notes.  This will become i gt j generally not just occ and unocc for selected CI
dtemp=0.0D0
!do j=1,docc
!do i=docc+1,nbft
do j=1,nbft-1
do i=j+1,nbft ! as selected CI all pairs where i>j

dtemp=dtemp+(Zvec(i,j)*B(i,j))
end do
end do



!add in two electron deriv terms using back transformed Z O(N^4) not O(N^5) each time like transforming AO 2e derivs to MO

do s=1,nbft
do r=1,nbft
do q=1,nbft
do p=1,nbft
dtemp=dtemp+(TwoEgrad_atom(p,q,r,s)*&
(2.0D0*Zback(p,r)*D(q,s)-Zback(p,s)*D(q,r)))

end do 
end do
end do
end do


gradients(coord,atom)=gradients(coord,atom)+2.0D0*dtemp 


end do !loop over coords
end do !loop over atoms


deallocate(Sgrad_atom)
deallocate(O1)
deallocate(H1)
deallocate(Hgrad_atom)  
deallocate(TwoEgrad_atom)
deallocate(TwoEgrad_MO)
deallocate(B)
deallocate(Zback)
deallocate(D)

end subroutine contractZwithBSelectedCI






subroutine calcZvecSelectedCI(docc)  !called by selectedCIgradZvec ! Think this is ok  !modifying to see if can cope with degenerate orbital energies
use AOintegral_arrays, only: nbft,FockE,MOe2ints
use CIgrad_arrays, only:  Zvec,X,deltaXthresh,deltaOrbEthresh
implicit none
integer, intent(in) :: docc
integer j,sigma,k,r ! use notation of Gerratt and Mills JCP 49, 1719 (1968) 
integer n,m
integer i,i2
double precision , allocatable:: A(:,:)
double precision , allocatable:: y(:)
double precision dtemp,deltaX,deltaOrbE
integer Adim
integer, allocatable:: IPIV(:)
double precision , allocatable:: work(:)
integer Lwork,info
double precision , allocatable:: deltaXtilde(:,:)



Zvec=0.0D0
n=docc
m=nbft

do k=1,n-1
do r=k+1,n

deltaX=X(r,k)-X(k,r)
if(ABS(deltaX).lt.deltaXthresh) THEN
cycle 
END IF

deltaOrbE=FockE(k)-FockE(r)
 ! set degenerate orbital terms to zero even if deltaX not small and give warning
if(ABS(deltaOrbE).lt.deltaOrbEthresh) THEN
PRINT *, 'In calcZvecSelectedCI'
PRINT *, 'problem orb',k,'too close in energy with',r
PRINT *, 'deltaOrbE',deltaOrbE,'but deltaX',deltaX
cycle
END IF


Zvec(r,k)=deltaX/deltaOrbE

end do
end do

!unocc-unocc
do k=n+1,m-1
do r=k+1,m

deltaX=X(r,k)-X(k,r)
if(ABS(deltaX).lt.deltaXthresh) THEN
!PRINT *,k,r,deltaX
cycle ! ignore this contribution so it is set to zero
END IF

deltaOrbE=FockE(k)-FockE(r)
 ! set degenerate orbital terms to zero even if deltaX not small and give warning
if(ABS(deltaOrbE).lt.deltaOrbEthresh) THEN
PRINT *, 'In calcZvecSelectedCI'
PRINT *, 'problem orb',k,'too close in energy with',r
PRINT *, 'deltaOrbE',deltaOrbE,'but deltaX',deltaX
cycle
END IF


Zvec(r,k)=deltaX/deltaOrbE

end do
end do






allocate(deltaXtilde(nbft,nbft)) 

deltaXtilde=0.0D0


do sigma=n+1,m
do j=1,n
deltaXtilde(sigma,j)=X(sigma,j)-X(j,sigma)
end do
end do


! need values for RHF indendenpent pairs
do sigma=n+1,m
do j=1,n

!occ-occ terms
do k=1,n-1
do r=k+1,n

deltaX=X(r,k)-X(k,r)
if(ABS(deltaX).lt.deltaXthresh) cycle 

deltaOrbE=FockE(k)-FockE(r)
if(ABS(deltaOrbE).lt.deltaOrbEthresh) cycle 

! remember  have to change order of coulomb integrals to Gerrat and Mills (swap middle two indices) 

dtemp=4.0D0*MOe2ints(sigma,k,j,r)-&
MOe2ints(sigma,r,k,j)-MOe2ints(sigma,k,r,j)

deltaXtilde(sigma,j)=deltaXtilde(sigma,j)&
+(dtemp*deltaX/deltaOrbE)

end do
end do

!unocc-unocc pairs
do k=n+1,m-1
do r=k+1,m

deltaX=X(r,k)-X(k,r)
if(ABS(deltaX).lt.deltaXthresh) cycle 

deltaOrbE=FockE(k)-FockE(r)
if(ABS(deltaOrbE).lt.deltaOrbEthresh) cycle 

! remember  have to change order of coulomb integrals to Gerrat and Mills (swap middle two indices) 

dtemp=4.0D0*MOe2ints(sigma,k,j,r)-&
MOe2ints(sigma,r,k,j)-MOe2ints(sigma,k,r,j)

deltaXtilde(sigma,j)=deltaXtilde(sigma,j)&
+(dtemp*deltaX/deltaOrbE)

end do
end do



end do ! end of loop over r
end do ! end of loop over k






!map  to occ with unocc pairs (m-n)*n  




!create y
Adim=n*(m-n)
allocate(y(Adim))
y=0.0D0
do sigma=n+1,m
do j=1,n
!map to 1D
i=j+(sigma-n-1)*n
y(i)=deltaXtilde(sigma,j) !deltaXtilde now

end do
end do




!create A which is now transposed note Gerrat and Mills Coulomb MO integrals are [x1,x1,x2,x2] but ours are [x1,x2,x1,x2]
allocate(A(Adim,Adim))
A=0.0D0
do sigma=n+1,m
do j=1,n

!map to 1D
i=j+(sigma-n-1)*n

A(i,i)=FockE(j)-FockE(sigma)

do r=n+1,m
do k=1,n

!map to 1D
i2=k+(r-n-1)*n
! remember  have to change order of coulomb integrals to Gerrat and Mills (swap middle two indices) and sign as we have put on left side of CPHF

dtemp=4.0D0*MOe2ints(sigma,k,j,r)-&
MOe2ints(sigma,r,k,j)-MOe2ints(sigma,k,r,j)

A(i2,i)=A(i2,i)-dtemp ! transpose and change sign here

end do
end do


end do
end do




!solve symmetric linear system then map back to sigma and j
allocate(IPIV(Adim))
Lwork=64*Adim  
ALLOCATE(work(Lwork))

!PRINT *,y
 Call dsysv('U',Adim,1,A,Adim,IPIV,&
 y,Adim,work,Lwork,info)  ! on output y is now the solution z
 
if(info.ne.0) STOP 'Problem with dsysv in calcZvecSelectedCI'




!need to map to sigma and j then store in 2D array Zvec where we alread have the occ -occ and unocc-unocc values

do sigma=n+1,m
do j=1,n

i=j+(sigma-n-1)*n

Zvec(sigma,j)=y(i)

end do
end do





deallocate(work)
deallocate(IPIV)
deallocate(A)
deallocate(y)

deallocate(deltaXtilde)

end subroutine calcZvecSelectedCI








subroutine CISDgradientsZvec(docc)
use CIgrad_arrays
use AOintegral_arrays, only: nbft
use RHFgrad_arrays, only:e2grad,gradients,Sgrad,natoms,firstAO,lastAO,&
geometry,charge,rinvgrad,e1grad
implicit none
integer, intent(in) :: docc
integer i,coord,atom,j

call readin2RDM(docc) ! allocates TwoRDM and OneRDM and calculates energy using this as check



allocate(Zvec(nbft,nbft))  !we use map to 1D to get vector when solving linear equations then map back to 2D for ease of use 


allocate(Sgrad(3,nbft,nbft))

call readinSgrad() !reads in natoms, firstAO array, last AO array and Sgrad array


allocate(gradients(3,natoms))  ! first index is x,y,z
gradients=0.0D0


allocate(e2grad(3,nbft,nbft,nbft,nbft))
call readinTwoEgrad()  !this may have already been called by RHF gradients but that is self contained so e2grad array has been deallocated

call pureAOtwoE() ! pure AO derivative 2e contribution creates and stores back transformation of 2RDM (back2RDM) rather than transforming all AO derivative integrals to MO basis


!need following to get things for oneE
allocate(geometry(3,natoms))
allocate(charge(natoms))

allocate(rinvgrad(natoms,3,nbft,nbft))



call readinGeometry()


allocate(e1grad(3,nbft,nbft))
call readinOneEgrad() 
call readin_rinvgrad() 


call pureAOoneE() !pure AO derivative 1e contribution also creates and store back transformation of 1RDM (back1RDM)





 call Vnngrad()  ! add in nuclear nuclear


allocate(X(nbft,nbft))
call calcXlagrangian()



call MO_CIcontributionOverlap() 

PRINT *,'CISD gradients without CPHF contrib'
do i=1,natoms
PRINT *,i,gradients(1,i),gradients(2,i),gradients(3,i) 
end do




call calcZvec(docc)




PRINT *,'CISD gradients with CPHF contrib using Z vector approach'
do i=1,natoms
PRINT *,i,gradients(1,i),gradients(2,i),gradients(3,i) 
end do

deallocate(TwoRDM)
deallocate(OneRDM)
deallocate(e2grad)
deallocate(gradients)
deallocate(Sgrad)
deallocate(firstAO)
deallocate(lastAO)
deallocate(geometry)
deallocate(charge)
deallocate(rinvgrad)
deallocate(e1grad)

deallocate(Zvec)
deallocate(X)
end subroutine CISDgradientsZvec



subroutine contractZwithB(docc)  
use AOintegral_arrays, only: nbft,MOe2ints,FockE,eigenFn
use RHFgrad_arrays, only:natoms,firstAO,lastAO,Sgrad,e1grad&
,rinvgrad,charge,e2grad,gradients
use CIgrad_arrays, only:  Sgrad_atom,O1,H1,Hgrad_atom,&
TwoEgrad_atom,TwoEgrad_MO,Zvec
implicit none
integer, intent(in) :: docc
integer atom,coord
integer j,sigma,k,r ,l! use notation of Gerratt and Mills JCP 49, 1719 (1968) 
integer n,m
integer i,p,q,s
double precision dtemp
double precision , allocatable:: B(:,:)
double precision , allocatable:: Zback(:,:)
double precision , allocatable:: D(:,:)


allocate(B(nbft,nbft))
allocate(Zback(nbft,nbft))
allocate(D(nbft,nbft))


n=docc
m=nbft



Zback=0.0D0



!made O(N^3)
!use D for temporary storage in back transform first
D=0.0D0

do i=docc+1,nbft
do j=1,docc
!do j=1,nbft-1  ! swap for docc loop etc if using selected CI
!do i=j+1,nbft
do p=1,nbft

D(p,j)=D(p,j)+&
(eigenFn(p,i)*Zvec(i,j))

end do
end do
end do

do j=1,docc
!do j=1,nbft-1 !swap for above if using selected CI
do q=1,nbft
do p=1,nbft

Zback(p,q)=Zback(p,q)+&
(eigenFn(q,j)*D(p,j))

end do
end do
end do






!also pre contract MO coeffs to allow O(N^4) two electron derivative terms 
D=0.0D0
do k=1,n
do q=1,nbft
do p=1,nbft


D(p,q)=D(p,q)+&
(eigenFn(p,k)*eigenFn(q,k))

end do
end do
end do





allocate(Sgrad_atom(nbft,nbft))
allocate(O1(nbft,nbft))
allocate(H1(nbft,nbft))
allocate(Hgrad_atom(nbft,nbft))  
allocate(TwoEgrad_atom(nbft,nbft,nbft,nbft))  
  allocate(TwoEgrad_MO(nbft,nbft,nbft,nbft)) 

do atom=1,natoms
do coord=1,3


!Transform deriv integrals to MO. 

Sgrad_atom=0.0D0
!copy relevant part of Sgrad to Sgrad atom to get S deriv w.r.t nuclear
do j=1,nbft
do i=firstAO(atom),lastAO(atom) 
Sgrad_atom(i,j)=-1.0D0*Sgrad(coord,i,j) 
end do
end do


! this is only derivative on first index so have to add value with indices swapped as (i, deriv j) = (deriv j, i) for real AOs
do i=1,nbft
do j=firstAO(atom),lastAO(atom)
Sgrad_atom(i,j)=Sgrad_atom(i,j)-1.0D0*Sgrad(coord,j,i) 
end do
end do


!transform Sgrad_atom to MO basis to give O1
call calcO1() 

!now H1   


!Copy parts on atom of e1grad for orbital derivs

Hgrad_atom=0.0D0
!copy relevant part of e1grad to Hgrad_atom to get one e deriv w.r.t nuclear
do j=1,nbft
do i=firstAO(atom),lastAO(atom) 
Hgrad_atom(i,j)=-1.0D0*e1grad(coord,i,j) 
end do
end do

! this is only derivative on first index so have to add value with indices swapped as (i, deriv j) = (deriv j, i) for real AOs
do i=1,nbft
do j=firstAO(atom),lastAO(atom)
Hgrad_atom(i,j)=Hgrad_atom(i,j)-1.0D0*e1grad(coord,j,i) 
end do
end do

!also need derivative of One e Hamiltonian term
do i=1,nbft
do j=1,nbft 
Hgrad_atom(i,j)=Hgrad_atom(i,j)-charge(atom)*rinvgrad(atom,coord,i,j)
Hgrad_atom(i,j)=Hgrad_atom(i,j)-charge(atom)*rinvgrad(atom,coord,j,i) ! need deriv on other index 
end do
end do

! transform Hgrad_atom to MO basis to give H1
call calcH1() 

!now Two electron gradients


TwoEgrad_atom=0.0D0
!
do i=firstAO(atom),lastAO(atom)
do j=1,nbft
do k=1,nbft
do r=1,nbft
TwoEgrad_atom(i,j,k,r)=-1.0D0*e2grad(coord,i,j,k,r) !e2grad only has derivative on first index after coord
end do
end do
end do
end do

!second index
do i=1,nbft
do j=firstAO(atom),lastAO(atom)
do k=1,nbft
do r=1,nbft
TwoEgrad_atom(i,j,k,r)=TwoEgrad_atom(i,j,k,r)-1.0D0*e2grad(coord,j,i,r,k) ! have to do two swaps to make j first index
end do
end do
end do
end do


do i=1,nbft
do j=1,nbft
do k=firstAO(atom),lastAO(atom)
do r=1,nbft
TwoEgrad_atom(i,j,k,r)=TwoEgrad_atom(i,j,k,r)-1.0D0*e2grad(coord,k,j,i,r) ! only one swap as both x1
end do
end do
end do
end do

do i=1,nbft
do j=1,nbft
do k=1,nbft
do r=firstAO(atom),lastAO(atom)
TwoEgrad_atom(i,j,k,r)=TwoEgrad_atom(i,j,k,r)-1.0D0*e2grad(coord,r,k,j,i) ! two swaps
end do
end do
end do
end do









!create B which is all terms independent of u in CPHF

!for selected CI this will loop over all but only need unocc occ for CISD
B=0.0D0
do sigma=n+1,m
do j=1,n


B(sigma,j)=H1(sigma,j)-FockE(j)*O1(sigma,j)


do k=1,n

do l=1,n
B(sigma,j)=B(sigma,j)-O1(k,l)*(2.0D0*MOe2ints(sigma,l,j,k)-&
MOe2ints(sigma,l,k,j))
end do
end do ! end of k loop

end do
end do


!contract with Zvec  
dtemp=0.0D0
do j=1,docc
do i=docc+1,nbft
dtemp=dtemp+(Zvec(i,j)*B(i,j))
end do
end do



!add in two electron deriv terms using back transformed Z O(N^4) not O(N^5) each time like transforming AO 2e derivs to MO

do s=1,nbft
do r=1,nbft
do q=1,nbft
do p=1,nbft
dtemp=dtemp+(TwoEgrad_atom(p,q,r,s)*&
(2.0D0*Zback(p,r)*D(q,s)-Zback(p,s)*D(q,r)))

end do 
end do
end do
end do


gradients(coord,atom)=gradients(coord,atom)+2.0D0*dtemp 


end do !loop over coords
end do !loop over atoms


deallocate(Sgrad_atom)
deallocate(O1)
deallocate(H1)
deallocate(Hgrad_atom)  
deallocate(TwoEgrad_atom)
deallocate(TwoEgrad_MO)
deallocate(B)
deallocate(Zback)
deallocate(D)

end subroutine contractZwithB


subroutine calcZvec(docc)
use AOintegral_arrays, only: nbft,FockE,MOe2ints
use CIgrad_arrays, only:  Zvec,X
implicit none
integer, intent(in) :: docc
integer j,sigma,k,r ! use notation of Gerratt and Mills JCP 49, 1719 (1968) 
integer n,m
integer i,i2
double precision , allocatable:: A(:,:)
double precision , allocatable:: y(:)
double precision dtemp
integer Adim
integer, allocatable:: IPIV(:)
double precision , allocatable:: work(:)
integer Lwork,info




n=docc
m=nbft


!map  to occ with unocc pairs (m-n)*n  



!create y
Adim=n*(m-n)
allocate(y(Adim))
y=0.0D0
do sigma=n+1,m
do j=1,n
!map to 1D
i=j+(sigma-n-1)*n
y(i)=X(sigma,j)-X(j,sigma)

end do
end do




!create A which is now transposed note Gerrat and Mills Coulomb MO integrals are [x1,x1,x2,x2] but ours are [x1,x2,x1,x2]
allocate(A(Adim,Adim))
A=0.0D0
do sigma=n+1,m
do j=1,n

!map to 1D
i=j+(sigma-n-1)*n

A(i,i)=FockE(j)-FockE(sigma)

do r=n+1,m
do k=1,n

!map to 1D
i2=k+(r-n-1)*n
! remember  have to change order of coulomb integrals to Gerrat and Mills (swap middle two indices) and sign as we have put on left side of CPHF

dtemp=4.0D0*MOe2ints(sigma,k,j,r)-&
MOe2ints(sigma,r,k,j)-MOe2ints(sigma,k,r,j)

A(i2,i)=A(i2,i)-dtemp ! transpose and change sign here

end do
end do


end do
end do




!solve symmetric linear system then map back to sigma and j
allocate(IPIV(Adim))
Lwork=64*Adim 
ALLOCATE(work(Lwork))

!PRINT *,y
 Call dsysv('U',Adim,1,A,Adim,IPIV,&
 y,Adim,work,Lwork,info)  ! on output y is now the solution z
 
if(info.ne.0) STOP 'Problem with dsysv in calcZvec'




!need to map to sigma and j then store in 2D array Zvec
Zvec=0.0D0
do sigma=n+1,m
do j=1,n

i=j+(sigma-n-1)*n

Zvec(sigma,j)=y(i)

end do
end do





deallocate(work)
deallocate(IPIV)
deallocate(A)
deallocate(y)


end subroutine calcZvec






subroutine CISDgradients(docc) 
use CIgrad_arrays
use AOintegral_arrays, only: nbft
use RHFgrad_arrays, only:e2grad,gradients,Sgrad,natoms,firstAO,lastAO,&
geometry,charge,rinvgrad,e1grad
use RHFpolari_arrays, only:u1 
implicit none
integer, intent(in) :: docc
integer i,coord,atom,j

call readin2RDM(docc) ! allocates TwoRDM and OneRDM and calculates energy using this as check


allocate(u1(nbft,nbft))

allocate(Sgrad(3,nbft,nbft))

call readinSgrad() !reads in natoms, firstAO array, last AO array and Sgrad array


allocate(gradients(3,natoms))  ! first index is x,y,z
gradients=0.0D0


allocate(e2grad(3,nbft,nbft,nbft,nbft))
call readinTwoEgrad()  !this may have already been called by RHF gradients but that is self contained so e2grad array has been deallocated

call pureAOtwoE() ! pure AO derivative 2e contribution creates and stores back transformation of 2RDM (back2RDM) rather than transforming all AO derivative integrals to MO basis


!need following to get things for oneE
allocate(geometry(3,natoms))
allocate(charge(natoms))

allocate(rinvgrad(natoms,3,nbft,nbft))



call readinGeometry()


allocate(e1grad(3,nbft,nbft))
call readinOneEgrad() 
call readin_rinvgrad() 


call pureAOoneE() !pure AO derivative 1e contribution also creates and store back transformation of 1RDM (back1RDM)





 call Vnngrad()  ! add in nuclear nuclear


allocate(X(nbft,nbft))
call calcXlagrangian()



call MO_CIcontributionOverlap() 

PRINT *,'CISD gradients without CPHF contrib'
do i=1,natoms
PRINT *,i,gradients(1,i),gradients(2,i),gradients(3,i) 
end do




do atom=1,natoms
do coord=1,3
call solveCPHF(docc,atom,coord)  
call MO_CIcontributionCPHF(docc,atom,coord)
end do
end do

PRINT *,'CISD gradients with CPHF contrib'
do i=1,natoms
PRINT *,i,gradients(1,i),gradients(2,i),gradients(3,i) 
end do

deallocate(TwoRDM)
deallocate(OneRDM)
deallocate(e2grad)
deallocate(gradients)
deallocate(Sgrad)
deallocate(firstAO)
deallocate(lastAO)
deallocate(geometry)
deallocate(charge)
deallocate(rinvgrad)
deallocate(e1grad)

deallocate(u1)
deallocate(X)
end subroutine CISDgradients


subroutine MO_CIcontributionCPHF(docc,atom,coord)   
use AOintegral_arrays, only: eigenFn,nbft
use RHFpolari_arrays, only: u1
use CIgrad_arrays, only: X
use RHFgrad_arrays, only:gradients
implicit none
integer, intent(in) :: docc
integer, intent(in) :: atom
integer, intent(in) :: coord
integer i,j
double precision dtemp





dtemp=0.0D0
!do i=1,docc
!do j=docc+1,nbft
!dtemp=dtemp+((X(i,j)-X(j,i))*u1(i,j))
!end do
!end do

!  i>j 
do j=1,docc
do i=docc+1,nbft
dtemp=dtemp+((X(i,j)-X(j,i))*u1(i,j))
end do
end do



gradients(coord,atom)=gradients(coord,atom)+2.0D0*dtemp 


end subroutine MO_CIcontributionCPHF



subroutine MO_CIcontributionOverlap() 
use AOintegral_arrays, only: nbft,MOe2ints,MOe1ints,eigenFn
use CIgrad_arrays, only: OneRDM,TwoRDM,X
use RHFgrad_arrays, only:natoms,Sgrad,firstAO,lastAO,gradients
implicit none
double precision , allocatable:: backX(:,:)
double precision , allocatable:: Tran1(:,:)
integer t,p,r,s,q,i,j
double precision dtemp
integer coord,atom


allocate(backX(nbft,nbft))
allocate(Tran1(nbft,nbft))





!Back transform X noting that now we have Xij when i.le. j  and Xji otherwise 
BackX=0.0D0
Tran1=0.0D0


do p=1,nbft
do j=1,nbft
do i=1,nbft
if(i.le.j) THEN
Tran1(p,j)=Tran1(p,j)+&
(eigenFn(p,i)*X(i,j))
ELSE
Tran1(p,j)=Tran1(p,j)+&
(eigenFn(p,i)*X(j,i))
END IF

end do
end do
end do

do p=1,nbft
do q=1,nbft
do j=1,nbft
BackX(p,q)=BackX(p,q)+&
(eigenFn(q,j)*Tran1(p,j))
end do
end do
end do




!!!!!!!!!!!!!!!!!!!!!!!



!overlap gradients
do atom=1,natoms
do coord=1,3

dtemp=0.0D0
do j=1,nbft
do i=firstAO(atom),lastAO(atom) 

dtemp=dtemp+(backX(i,j)*Sgrad(coord,i,j)) 

end do
end do

! derivative on ket terms
do i=1,nbft
do j=firstAO(atom),lastAO(atom) 

dtemp=dtemp+(backX(i,j)*Sgrad(coord,j,i))

end do
end do


gradients(coord,atom)=gradients(coord,atom)+dtemp
end do
end do



deallocate(backX)
deallocate(Tran1)
end subroutine MO_CIcontributionOverlap



subroutine calcXlagrangian()
use AOintegral_arrays, only: nbft,MOe2ints,MOe1ints
use CIgrad_arrays, only: OneRDM,TwoRDM,X
implicit none
integer t,p,r,s,q
double precision dtemp




X=0.0D0


!Store X
do t=1,nbft
do p=1,nbft

dtemp=0.0D0
do q=1,nbft
dtemp=dtemp+OneRDM(p,q)*MOe1ints(t,q)
end do

do q=1,nbft
do s=1,nbft
do r=1,nbft
dtemp=dtemp+(0.5D0*(TwoRDM(p,r,s,q)+TwoRDM(q,r,s,p))*& 
MOe2ints(t,r,q,s))
end do
end do
end do

X(t,p)=dtemp
end do
end do


end subroutine calcXlagrangian




subroutine solveCPHF(docc,atom,coord) 
use AOintegral_arrays, only: nbft,FockE,MOe2ints
use RHFpolari_arrays, only: u1  !
use RHFgrad_arrays, only:natoms,firstAO,lastAO,Sgrad,e1grad&
,rinvgrad,charge,e2grad
use CIgrad_arrays, only:  Sgrad_atom,O1,H1,Hgrad_atom,&
TwoEgrad_atom,TwoEgrad_MO
implicit none
integer, intent(in) :: docc
integer, intent(in) :: atom
integer, intent(in) :: coord
integer j,sigma,k,r ,l! use notation of Gerratt and Mills JCP 49, 1719 (1968) 
integer n,m
integer i,tempj,tempsigma,i2
double precision , allocatable:: A(:,:)
double precision , allocatable:: y(:)
double precision dtemp
integer Adim
integer, allocatable:: IPIV(:)
double precision , allocatable:: work(:)
integer Lwork,info
double precision , allocatable:: F1(:,:)


allocate(Sgrad_atom(nbft,nbft))
allocate(O1(nbft,nbft))

Sgrad_atom=0.0D0
!copy relevant part of Sgrad to Sgrad atom to get S deriv w.r.t nuclear
do j=1,nbft
do i=firstAO(atom),lastAO(atom) 
Sgrad_atom(i,j)=-1.0D0*Sgrad(coord,i,j) 
end do
end do


! this is only derivative on first index so have to add value with indices swapped as (i, deriv j) = (deriv j, i) for real AOs
do i=1,nbft
do j=firstAO(atom),lastAO(atom)
Sgrad_atom(i,j)=Sgrad_atom(i,j)-1.0D0*Sgrad(coord,j,i) 
end do
end do


!transform Sgrad_atom to MO basis to give O1
call calcO1() 

!now H1   Gerrat and Mills use H for one-e and G for two-e terms
allocate(H1(nbft,nbft))
allocate(Hgrad_atom(nbft,nbft))  

!Copy parts on atom of e1grad for orbital derivs

Hgrad_atom=0.0D0
!copy relevant part of e1grad to Hgrad_atom to get one e deriv w.r.t nuclear
do j=1,nbft
do i=firstAO(atom),lastAO(atom) 
Hgrad_atom(i,j)=-1.0D0*e1grad(coord,i,j) 
end do
end do

! this is only derivative on first index so have to add value with indices swapped as (i, deriv j) = (deriv j, i) for real AOs
do i=1,nbft
do j=firstAO(atom),lastAO(atom)
Hgrad_atom(i,j)=Hgrad_atom(i,j)-1.0D0*e1grad(coord,j,i) 
end do
end do

!also need derivative of One e Hamiltonian term
do i=1,nbft
do j=1,nbft ! not orbital deriv so don't limit this
Hgrad_atom(i,j)=Hgrad_atom(i,j)-charge(atom)*rinvgrad(atom,coord,i,j)
Hgrad_atom(i,j)=Hgrad_atom(i,j)-charge(atom)*rinvgrad(atom,coord,j,i) 
end do
end do

! transform Hgrad_atom to MO basis to give H1
call calcH1() 

!now Two electron gradients
allocate(TwoEgrad_atom(nbft,nbft,nbft,nbft))  
  allocate(TwoEgrad_MO(nbft,nbft,nbft,nbft)) 

TwoEgrad_atom=0.0D0
!
do i=firstAO(atom),lastAO(atom)
do j=1,nbft
do k=1,nbft
do r=1,nbft
TwoEgrad_atom(i,j,k,r)=-1.0D0*e2grad(coord,i,j,k,r) !e2grad only has derivative on first index after coord
end do
end do
end do
end do

!second index
do i=1,nbft
do j=firstAO(atom),lastAO(atom)
do k=1,nbft
do r=1,nbft
TwoEgrad_atom(i,j,k,r)=TwoEgrad_atom(i,j,k,r)-1.0D0*e2grad(coord,j,i,r,k) ! have to do two swaps to make j first index
end do
end do
end do
end do


do i=1,nbft
do j=1,nbft
do k=firstAO(atom),lastAO(atom)
do r=1,nbft
TwoEgrad_atom(i,j,k,r)=TwoEgrad_atom(i,j,k,r)-1.0D0*e2grad(coord,k,j,i,r) ! only one swap as both x1
end do
end do
end do
end do

do i=1,nbft
do j=1,nbft
do k=1,nbft
do r=firstAO(atom),lastAO(atom)
TwoEgrad_atom(i,j,k,r)=TwoEgrad_atom(i,j,k,r)-1.0D0*e2grad(coord,r,k,j,i) ! two swaps
end do
end do
end do
end do


call calcTwoEGrad_MO() 


!now have the H1, O1 , TwoEGrad_MO and already had Moe2ints so put cphf together

n=docc
m=nbft



!map CPHF eqn 50 from JCP 49, 1719 (1968) to (m-n)*n  


!Want to solve Ax=y from rearranging CPHF eqn 50
! y is H^{1}_sigma,j where we map sigma and j to i
!create y
Adim=n*(m-n)
allocate(y(Adim))
y=0.0D0
do sigma=n+1,m
do j=1,n
!map to 1D
i=j+(sigma-n-1)*n
y(i)=H1(sigma,j)-FockE(j)*O1(sigma,j)

do k=1,n
y(i)=y(i)+(2.0D0*TwoEgrad_MO(sigma,k,j,k))-&  !note we need to swap middle two compared with Gerrat and Mills CPHF for our notation
TwoEgrad_MO(sigma,k,k,j)
do l=1,n
y(i)=y(i)-O1(k,l)*(2.0D0*MOe2ints(sigma,l,j,k)-&
MOe2ints(sigma,l,k,j))
end do
end do ! end of k loop

end do
end do




!create A note Gerrat and Mills Coulomb MO integrals are [x1,x1,x2,x2] but ours are [x1,x2,x1,x2]
allocate(A(Adim,Adim))
A=0.0D0
!PRINT *,FockE
do sigma=n+1,m
do j=1,n

!map to 1D
i=j+(sigma-n-1)*n

A(i,i)=FockE(j)-FockE(sigma)

do r=n+1,m
do k=1,n

!map to 1D
i2=k+(r-n-1)*n
! remember  have to change order of coulomb integrals to Gerrat and Mills (swap middle two indices) and sign as we have put on left side of CPHF

dtemp=4.0D0*MOe2ints(sigma,k,j,r)-&
MOe2ints(sigma,r,k,j)-MOe2ints(sigma,k,r,j)
A(i,i2)=A(i,i2)-dtemp 

end do
end do


end do
end do



!solve symmetric linear system then map back to sigma and j
allocate(IPIV(Adim))
Lwork=64*Adim 
ALLOCATE(work(Lwork))

!PRINT *,y
 Call dsysv('U',Adim,1,A,Adim,IPIV,&
 y,Adim,work,Lwork,info)  ! on output y is now the solution x
 
if(info.ne.0) STOP 'Problem with dsysv in solveCPHF'

!PRINT *,y

!seems ok to here


!need to map to sigma and j then store in u1
u1=0.0D0
do sigma=n+1,m
do j=1,n

i=j+(sigma-n-1)*n

u1(sigma,j)=y(i)

end do
end do

!PRINT *,u1



!get remaining u1 by calculating F1
allocate(F1(m,m))
F1=0.0D0
!one electron contribution
do j=1,m
do sigma=1,m
F1(sigma,j)=H1(sigma,j)
end do
end do
!two electron contribution these depend on the u1 calculated by solving the previous Ax=y
do j=1,m
do sigma=1,m

do r=n+1,m
do k=1,n

dtemp=4.0D0*MOe2ints(sigma,k,j,r)-&
MOe2ints(sigma,r,k,j)-MOe2ints(sigma,k,r,j)

dtemp=dtemp*u1(r,k)

F1(sigma,j)=F1(sigma,j)+dtemp
end do
end do

dtemp=0.0D0
do k=1,n
dtemp=dtemp+(2.0D0*TwoEgrad_MO(sigma,k,j,k))-&  !note we need to swap middle two compared with Gerrat and Mills CPHF for our notation
TwoEgrad_MO(sigma,k,k,j)
do l=1,n
dtemp=dtemp-O1(k,l)*(2.0D0*MOe2ints(sigma,l,j,k)-&
MOe2ints(sigma,l,k,j))
end do
end do ! end of k loop

F1(sigma,j)=F1(sigma,j)+dtemp




end do
end do


do j=n+1,m  ! swap j and sigma to get remaining non zero
do sigma=1,n


if(j.eq.sigma) cycle


u1(sigma,j)=(F1(sigma,j)-FockE(j)*O1(sigma,j))/(FockE(j)-FockE(sigma))


end do
end do




deallocate(work)
deallocate(IPIV)
deallocate(A)
deallocate(y)
deallocate(F1)

deallocate(Sgrad_atom)
deallocate(O1)
deallocate(H1)
deallocate(Hgrad_atom)  
deallocate(TwoEgrad_atom)
deallocate(TwoEgrad_MO)
end subroutine solveCPHF

subroutine calcTwoEGrad_MO() 
use AOintegral_arrays, only: eigenFn,nbft
use CIgrad_arrays, only:  TwoEgrad_atom,TwoEgrad_MO
implicit none
integer i,j,k,m,p,q,r,s
double precision , allocatable:: Tran1(:,:,:,:)
integer info,pr,qs



allocate(Tran1(nbft,nbft,nbft,nbft),stat=info)
if(info.ne.0) STOP "Error allocating 4 dim array Tran1 in calcTwoEGrad_MO"



Tran1=0.0D0
do m=1,nbft
do k=1,nbft
do j=1,nbft
do p=1,nbft
do i=1,nbft
Tran1(p,j,k,m)=Tran1(p,j,k,m)+&
(eigenFn(i,p)*TwoEgrad_atom(i,j,k,m))
end do
end do
end do
end do 
end do

TwoEgrad_MO=0.0D0 !this is used as temporary storage for transformation at this stage
do m=1,nbft
do k=1,nbft
do q=1,nbft
do j=1,nbft
do p=1,nbft
TwoEgrad_MO(p,q,k,m)=TwoEgrad_MO(p,q,k,m)+&
(eigenFn(j,q)*Tran1(p,j,k,m))
end do
end do
end do
end do 
end do


Tran1=0.0D0 !reuse this storage space
do m=1,nbft
do r=1,nbft
do k=1,nbft
do q=1,nbft
do p=1,nbft
Tran1(p,q,r,m)=Tran1(p,q,r,m)+&
(eigenFn(k,r)*TwoEgrad_MO(p,q,k,m)) 
end do
end do
end do
end do 
end do

!final step to get MO 2e integrals in Tran2

TwoEgrad_MO=0.0D0 ! reuse this space for final values
do s=1,nbft
do m=1,nbft
do r=1,nbft
do q=1,nbft
do p=1,nbft
TwoEgrad_MO(p,q,r,s)=TwoEgrad_MO(p,q,r,s)+&
(eigenFn(m,s)*Tran1(p,q,r,m))
end do
end do
end do
end do 
end do




deallocate(Tran1)

end subroutine calcTwoEGrad_MO









subroutine calcH1() 
use AOintegral_arrays, only: eigenFn,nbft
use CIgrad_arrays, only:  Hgrad_atom,H1
implicit none
integer i,j,p,q
double precision , allocatable:: Tran1(:,:)


allocate(Tran1(nbft,nbft))



H1=0.0D0

Tran1=0.0D0


do p=1,nbft
do j=1,nbft
do i=1,nbft
Tran1(p,j)=Tran1(p,j)+&
(eigenFn(i,p)*Hgrad_atom(i,j))
end do
end do
end do

do p=1,nbft
do q=1,nbft
do j=1,nbft
H1(p,q)=H1(p,q)+&
(eigenFn(j,q)*Tran1(p,j))
end do
end do
end do



deallocate(Tran1)
end subroutine calcH1


subroutine calcO1() 
use AOintegral_arrays, only: eigenFn,nbft
use CIgrad_arrays, only:  Sgrad_atom,O1
implicit none
integer i,j,p,q
double precision , allocatable:: Tran1(:,:)


allocate(Tran1(nbft,nbft))



O1=0.0D0

Tran1=0.0D0


do p=1,nbft
do j=1,nbft
do i=1,nbft
Tran1(p,j)=Tran1(p,j)+&
(eigenFn(i,p)*Sgrad_atom(i,j))
end do
end do
end do

do p=1,nbft
do q=1,nbft
do j=1,nbft
O1(p,q)=O1(p,q)+&
(eigenFn(j,q)*Tran1(p,j))
end do
end do
end do






deallocate(Tran1)
end subroutine calcO1










subroutine RHFpolarizability(docc)
use RHFpolari_arrays
use AOintegral_arrays, only: nbft
use RHFgrad_arrays, only: geometry,charge,natoms
implicit none
integer, intent(in) :: docc
integer i,j
double precision polariz(3)

allocate(dipoleints(3,nbft,nbft))
allocate(u1(nbft,nbft))

!to get nuclear contribution to dipole
OPEN(UNIT=14,FILE='OverlapGrad.txt') 
READ (14,*) natoms
CLOSE(14)
allocate(geometry(3,natoms))
allocate(charge(natoms))

call readinGeometry()  ! get charge and coords for nuclear dipole

call readinDipoleInts() 

call calcRHFdipole()




!call checkFock0Tran() !this works fine

!transform the dipoleints to the MO basis
allocate(MOdipoleints(3,nbft,nbft))

call calcMOdipoleInts()



!store polarizability array  
PRINT *, 'Polarizability'
do j=1,3
call solveCPHFpolariz(docc,j) ! this coord j is for derivative of density matrix


call calcPolariz(docc,polariz) 

PRINT *,polariz
end do





deallocate(dipoleints)
deallocate(geometry)
deallocate(charge)
deallocate(MOdipoleints)
deallocate(u1)

end subroutine RHFpolarizability


subroutine calcPolariz(docc,polariz) 
use AOintegral_arrays, only: eigenFn,nbft
use RHFpolari_arrays, only: u1,dipoleints
implicit none
integer, intent(in) :: docc
double precision, intent(out) :: polariz(3)
double precision , allocatable:: c1(:,:)
double precision alpha,beta
double precision , allocatable:: D1(:,:)
double precision , allocatable:: Dtemp(:,:)
integer i,j,coord



!transform u1 to c1  c1=C0 u1    [C0 is eigenFn matrix]
allocate(c1(nbft,nbft))
c1=0.0D0

!DGEMM to perform C0*u1=c1
alpha=1.0D0
beta=0.0D0
Call dgemm('N','N',nbft,nbft,nbft,alpha,&
eigenFn,nbft,u1,nbft,beta,c1,nbft)

!Form D1 derivative of density matrix 
allocate(D1(nbft,nbft))
allocate(Dtemp(nbft,nbft))
D1=0.0D0
Dtemp=0.0D0
!C1*C0(Tran) term
alpha=2.0D0  ! as all double occupied -- every spin up also has a spin down
beta=0.0D0
Call dgemm('N','T',nbft,nbft,docc,alpha,&   
c1,nbft,eigenFn,nbft,beta,Dtemp,nbft)  ! note docc used as only want these columns of eigenFn

!C0*C1(Tran) term
alpha=2.0D0  ! as all double occupied -- every spin up also has a spin down
beta=0.0D0
Call dgemm('N','T',nbft,nbft,docc,alpha,&   
eigenFn,nbft,c1,nbft,beta,D1,nbft)  ! note docc used as only want these columns of eigenFn

D1=D1+Dtemp

!calculate polarization chosen by  combining coord here and coord in solveCPHFpolariz. Here we have the coord for derivative of one-electron operator
polariz=0.0D0
do coord=1,3

do j=1,nbft
do i=1,nbft
polariz(coord)=polariz(coord)-(D1(i,j)*dipoleints(coord,i,j)) 
end do
end do

end do


deallocate(c1)
deallocate(D1)
deallocate(Dtemp)
end subroutine calcPolariz


subroutine solveCPHFpolariz(docc,coord)  
use AOintegral_arrays, only: nbft,FockE,MOe2ints
use RHFpolari_arrays, only:MOdipoleints,u1
implicit none
integer, intent(in) :: docc
integer, intent(in) :: coord
integer j,sigma,k,r ! use notation of Gerratt and Mills JCP 49, 1719 (1968) 
integer n,m
integer i,tempj,tempsigma,i2
double precision , allocatable:: A(:,:)
double precision , allocatable:: y(:)
double precision dtemp
integer Adim
integer, allocatable:: IPIV(:)
double precision , allocatable:: work(:)
integer Lwork,info
double precision , allocatable:: F1(:,:)



n=docc
m=nbft






! y is H^{1}_sigma,j where we map sigma and j to i
!create y
Adim=n*(m-n)
allocate(y(Adim))
y=0.0D0
do sigma=n+1,m
do j=1,n
!map to 1D
i=j+(sigma-n-1)*n
y(i)=MOdipoleints(coord,sigma,j)
end do
end do

!create A note Gerrat and Mills Coulomb MO integrals are [x1,x1,x2,x2] but ours are [x1,x2,x1,x2]
allocate(A(Adim,Adim))
A=0.0D0
!PRINT *,FockE
do sigma=n+1,m
do j=1,n

!map to 1D
i=j+(sigma-n-1)*n

A(i,i)=FockE(j)-FockE(sigma)

do r=n+1,m
do k=1,n

!map to 1D
i2=k+(r-n-1)*n
! remember  have to change order of coulomb integrals to Gerrat and Mills (swap middle two indices) and sign as we have put on left side of CPHF

dtemp=4.0D0*MOe2ints(sigma,k,j,r)-&
MOe2ints(sigma,r,k,j)-MOe2ints(sigma,k,r,j)
A(i,i2)=A(i,i2)-dtemp !we change sign here

end do
end do


end do
end do



!solve symmetric linear system then map back to sigma and j
allocate(IPIV(Adim))
Lwork=64*Adim  
ALLOCATE(work(Lwork))

!PRINT *,y
 Call dsysv('U',Adim,1,A,Adim,IPIV,&
 y,Adim,work,Lwork,info)  ! on output y is now the solution x
 
if(info.ne.0) STOP 'Problem with dsysv in solveCPHFpolariz'

!PRINT *,y

!need to map to sigma and j then store in u1
u1=0.0D0
do sigma=n+1,m
do j=1,n

i=j+(sigma-n-1)*n

u1(sigma,j)=y(i)

end do
end do

!PRINT *,u1

!get remaining u1 by calculating F1
allocate(F1(m,m))
F1=0.0D0
!one electron contribution
do j=1,m
do sigma=1,m
F1(sigma,j)=MOdipoleints(coord,sigma,j)
end do
end do
!two electron contribution these depend on the u1 calculated by solving the previous Ax=y
do j=1,m
do sigma=1,m

do r=n+1,m
do k=1,n

dtemp=4.0D0*MOe2ints(sigma,k,j,r)-&
MOe2ints(sigma,r,k,j)-MOe2ints(sigma,k,r,j)

dtemp=dtemp*u1(r,k)

F1(sigma,j)=F1(sigma,j)+dtemp

end do
end do


end do
end do





deallocate(work)
deallocate(IPIV)
deallocate(A)
deallocate(y)
deallocate(F1)
end subroutine solveCPHFpolariz




subroutine calcMOdipoleInts() 
use AOintegral_arrays, only: eigenFn,nbft
use RHFpolari_arrays, only:dipoleints,MOdipoleints
implicit none
integer i,j,p,q,coord
double precision , allocatable:: Tran1(:,:)


allocate(Tran1(nbft,nbft))



MOdipoleints=0.0D0

do coord=1,3
Tran1=0.0D0


do p=1,nbft
do j=1,nbft
do i=1,nbft
Tran1(p,j)=Tran1(p,j)+&
(eigenFn(i,p)*dipoleints(coord,i,j))
end do
end do
end do

do p=1,nbft
do q=1,nbft
do j=1,nbft
MOdipoleints(coord,p,q)=MOdipoleints(coord,p,q)+&
(eigenFn(j,q)*Tran1(p,j))
end do
end do
end do

end do !end of loop over coords




deallocate(Tran1)
end subroutine calcMOdipoleInts





subroutine checkFock0Tran() 
use AOintegral_arrays, only: Fock,eigenFn,nbft
implicit none
integer i,j,p,q
double precision , allocatable:: Fock0Tran(:,:)
double precision , allocatable:: Tran1(:,:)
allocate(Fock0Tran(nbft,nbft))

allocate(Tran1(nbft,nbft))




Fock0Tran=0.0D0
Tran1=0.0D0


do p=1,nbft
do j=1,nbft
do i=1,nbft
Tran1(p,j)=Tran1(p,j)+&
(eigenFn(i,p)*Fock(i,j))
end do
end do
end do

do p=1,nbft
do q=1,nbft
do j=1,nbft
Fock0Tran(p,q)=Fock0Tran(p,q)+&
(eigenFn(j,q)*Tran1(p,j))
end do
end do
end do


do i=1,nbft
PRINT *,(Fock0Tran(i,j),j=1,nbft)   
end do

deallocate(Fock0Tran)
deallocate(Tran1)


end subroutine checkFock0Tran



subroutine calcRHFdipole()  
use AOintegral_arrays, only: nbft,denmat
use RHFpolari_arrays, only:dipoleints
use RHFgrad_arrays, only: geometry,charge,natoms
implicit none
integer i,j,coord
double precision dipole(3),nucdipole(3)

!electronic contribution
do coord=1,3

dipole(coord)=0.0D0

do j=1,nbft
do i=1,nbft
dipole(coord)=dipole(coord)+(dipoleints(coord,i,j)*denmat(i,j))
end do
end do

end do


!nuclear contribution
do coord=1,3
nucdipole(coord)=0.0D0
do i=1,natoms
nucdipole(coord)=nucdipole(coord)+charge(i)*geometry(coord,i)
end do

dipole(coord)=nucdipole(coord)-dipole(coord)
end do

Print *,' Total dipole'
PRINT *,dipole


end subroutine calcRHFdipole




subroutine readinDipoleInts() 
use RHFpolari_arrays, only:dipoleints
implicit none
integer i,j,k,l
integer info
double precision dtemp

dipoleints=0.0D0

OPEN(UNIT=14,FILE='DipoleInts.txt') 


info=0
do while (info==0)  !loop until end of file
READ (14,*,IOSTAT=info)  dtemp,i,j,k  !i is x y or z, j and k label AOs
dipoleints(i,j,k)=dtemp
end do
CLOSE(14)



end subroutine readinDipoleInts






subroutine FCIgradients(docc) 
use CIgrad_arrays,only: oneRDM,twoRDM
use AOintegral_arrays, only: nbft
use RHFgrad_arrays, only:e2grad,gradients,Sgrad,natoms,firstAO,lastAO,&
geometry,charge,rinvgrad,e1grad
implicit none
integer, intent(in) :: docc
integer i

call readin2RDM(docc) ! allocates TwoRDM and OneRDM and calculates energy using this as check

allocate(Sgrad(3,nbft,nbft))

call readinSgrad() !reads in natoms, firstAO array, last AO array and Sgrad array


allocate(gradients(3,natoms))  ! first index is x,y,z
gradients=0.0D0


allocate(e2grad(3,nbft,nbft,nbft,nbft))
call readinTwoEgrad()  !this may have already been called by RHF gradients but that is self contained so e2grad array has been deallocated

call pureAOtwoE() ! pure AO derivative 2e contribution creates and stores back transformation of 2RDM (back2RDM) rather than transforming all AO derivative integrals to MO basis


!need following to get things for oneE
allocate(geometry(3,natoms))
allocate(charge(natoms))

allocate(rinvgrad(natoms,3,nbft,nbft))



call readinGeometry()


allocate(e1grad(3,nbft,nbft))
call readinOneEgrad() 
call readin_rinvgrad() 


call pureAOoneE() !pure AO derivative 1e contribution also creates and store back transformation of 1RDM (back1RDM)





 call Vnngrad()  ! add in nuclear nuclear


call MO_FCIcontribution() ! add in MO contribution for FCI only at the moment
PRINT *,'FCI gradients'
do i=1,natoms
PRINT *,i,gradients(1,i),gradients(2,i),gradients(3,i) 
end do

deallocate(TwoRDM)
deallocate(OneRDM)
deallocate(e2grad)
deallocate(gradients)
deallocate(Sgrad)
deallocate(firstAO)
deallocate(lastAO)
deallocate(geometry)
deallocate(charge)
deallocate(rinvgrad)
deallocate(e1grad)
end subroutine FCIgradients


subroutine MO_FCIcontribution() 
use AOintegral_arrays, only: nbft,MOe2ints,MOe1ints,eigenFn
use CIgrad_arrays, only: OneRDM,TwoRDM
use RHFgrad_arrays, only:natoms,Sgrad,firstAO,lastAO,gradients
implicit none
double precision , allocatable:: X(:,:)
double precision , allocatable:: backX(:,:)
double precision , allocatable:: Tran1(:,:)
integer t,p,r,s,q,i,j
double precision dtemp
integer coord,atom

allocate(X(nbft,nbft))
allocate(backX(nbft,nbft))
allocate(Tran1(nbft,nbft))

X=0.0D0


!Store X
do t=1,nbft
do p=1,nbft

dtemp=0.0D0
do q=1,nbft
dtemp=dtemp+OneRDM(p,q)*MOe1ints(t,q)
end do

do q=1,nbft
do s=1,nbft
do r=1,nbft
dtemp=dtemp+(0.5D0*(TwoRDM(p,r,s,q)+TwoRDM(q,r,s,p))*&
MOe2ints(t,r,q,s))
end do
end do
end do

X(t,p)=dtemp
end do
end do


!Back transform X
BackX=0.0D0
Tran1=0.0D0


do p=1,nbft
do j=1,nbft
do i=1,nbft
Tran1(p,j)=Tran1(p,j)+&
(eigenFn(p,i)*X(i,j))
end do
end do
end do

do p=1,nbft
do q=1,nbft
do j=1,nbft
BackX(p,q)=BackX(p,q)+&
(eigenFn(q,j)*Tran1(p,j))
end do
end do
end do

!overlap gradients
do atom=1,natoms
do coord=1,3

dtemp=0.0D0
do j=1,nbft
do i=firstAO(atom),lastAO(atom) 
dtemp=dtemp+(backX(i,j)*Sgrad(coord,i,j))
end do
end do


do i=1,nbft
do j=firstAO(atom),lastAO(atom) 
dtemp=dtemp+(backX(i,j)*Sgrad(coord,j,i))
end do
end do


gradients(coord,atom)=gradients(coord,atom)+dtemp
end do
end do


deallocate(X)
deallocate(backX)
deallocate(Tran1)
end subroutine MO_FCIcontribution

subroutine pureAOoneE()
use CIgrad_arrays, only: OneRDM
use AOintegral_arrays, only: nbft,eigenFn
use RHFgrad_arrays, only:natoms,e1grad,firstAO,lastAO,gradients,&
charge,rinvgrad
implicit none
double precision , allocatable:: Tran1(:,:)
integer p,q,i,j
integer atom,coord
double precision dtemp,dtemp2
double precision, allocatable ::  back1RDM(:,:)




allocate(Back1RDM(nbft,nbft))

allocate(Tran1(nbft,nbft))



Back1RDM=0.0D0
Tran1=0.0D0


do p=1,nbft
do j=1,nbft
do i=1,nbft
Tran1(p,j)=Tran1(p,j)+&
(eigenFn(p,i)*OneRDM(i,j))
end do
end do
end do

do p=1,nbft
do q=1,nbft
do j=1,nbft
Back1RDM(p,q)=Back1RDM(p,q)+&
(eigenFn(q,j)*Tran1(p,j))
end do
end do
end do


do atom=1,natoms
do coord=1,3

dtemp=0.0D0
do j=1,nbft
do i=firstAO(atom),lastAO(atom) 

dtemp=dtemp-e1grad(coord,i,j)*Back1RDM(i,j) 
end do
end do

dtemp=2.0D0*dtemp

!operator deriv
dtemp2=0.0D0
do j=1,nbft 
do i=1,nbft

dtemp2=dtemp2+rinvgrad(atom,coord,i,j)*Back1RDM(i,j)
end do
end do

dtemp2=2.0D0*dtemp2

dtemp=dtemp-charge(atom)*dtemp2

gradients(coord,atom)=gradients(coord,atom)+dtemp
end do
end do


deallocate(Tran1)
deallocate(Back1RDM)
end subroutine pureAOoneE


subroutine pureAOtwoE()
use CIgrad_arrays, only: TwoRDM
use RHFgrad_arrays, only:natoms,e2grad,firstAO,lastAO,gradients
use AOintegral_arrays, only: nbft,eigenFn
implicit none
double precision , allocatable:: Tran1(:,:,:,:)
integer p,q,r,s,i,j,k,m
double precision dtemp
integer atom,coord
double precision, allocatable ::  back2RDM(:,:,:,:)
!calculate back transformed 2RDM so we don't need to transform integrals for every degree of freedom
allocate(Back2RDM(nbft,nbft,nbft,nbft))

allocate(Tran1(nbft,nbft,nbft,nbft))


Tran1=0.0D0
do m=1,nbft
do k=1,nbft
do j=1,nbft
do p=1,nbft
do i=1,nbft
Tran1(p,j,k,m)=Tran1(p,j,k,m)+&
(eigenFn(p,i)*TwoRDM(i,j,k,m))
end do
end do
end do
end do 
end do

Back2RDM=0.0D0 !this is used as temporary storage for transformation at this stage
do m=1,nbft
do k=1,nbft
do q=1,nbft
do j=1,nbft
do p=1,nbft
Back2RDM(p,q,k,m)=Back2RDM(p,q,k,m)+&
(eigenFn(q,j)*Tran1(p,j,k,m))
end do
end do
end do
end do 
end do


Tran1=0.0D0 !reuse this storage space
do m=1,nbft
do r=1,nbft
do k=1,nbft
do q=1,nbft
do p=1,nbft
Tran1(p,q,r,m)=Tran1(p,q,r,m)+&
(eigenFn(r,k)*Back2RDM(p,q,k,m)) 
end do
end do
end do
end do 
end do

!final step to get back2RDM in Tran2

back2RDM=0.0D0 ! reuse this space for final values
do s=1,nbft
do m=1,nbft
do r=1,nbft
do q=1,nbft
do p=1,nbft
back2RDM(p,q,r,s)=back2RDM(p,q,r,s)+&
(eigenFn(s,m)*Tran1(p,q,r,m))
end do
end do
end do
end do 
end do

!now have back2RDM calculate pure AO 2electron part

do atom=1,natoms
do coord=1,3



dtemp=0.0D0
do s=1,nbft
do q=1,nbft
do r=1,nbft
do p=firstAO(atom),lastAO(atom)  ! if p orbital is not on atom then nuclear derivative is zero
dtemp=dtemp+back2RDM(p,r,s,q)*e2grad(coord,p,r,q,s)  
end do
end do
end do
end do



gradients(coord,atom)=gradients(coord,atom)-2.0D0*dtemp

end do
end do





deallocate(Tran1)
deallocate(Back2RDM)
end subroutine pureAOtwoE

subroutine readin2RDM_Arrow(docc,mystate) !in MO basis energy ordered
use AOintegral_arrays, only: nbft,MOe2ints,MOe1ints,ecore
use CIgrad_arrays, only: TwoRDM,OneRDM
implicit none
integer, intent(in) :: docc
integer, intent(in) :: mystate
integer info,p,q,r,s
integer i,j,k,m
double precision dtemp,energy
CHARACTER (LEN=30)  :: f1 
CHARACTER (LEN=30)  :: f2



TwoRDM=0.0D0
OneRDM=0.0D0



write(f1,'(i2)') mystate
 write(f2,'(i2)') mystate
            f1='state'//trim(adjustl(f1))
            f2='state'//trim(adjustl(f2))//'.txt'
            f1=trim(adjustl(f1))//trim(adjustl(f2))
            f1='Tran2RDM_'//trim(adjustl(f1))

PRINT *,'Read in ',f1


OPEN(UNIT=14,FILE=f1)  





READ(14,*) i  ! first entry is nbft
info=0
do while (info==0)  !loop until end of file
READ (14,*,IOSTAT=info)  p,q,r,s,dtemp

!need to add one to orbitals as from 0 to nbft-1 in SHCI


i=p+1
j=q+1
k=r+1
m=s+1



TwoRDM(i,j,m,k)=dtemp





end do
CLOSE(14)

!check we get 2electron energy
dtemp=0.0D0

do s=1,nbft
do r=1,nbft
do q=1,nbft
do p=1,nbft

dtemp=dtemp+TwoRDM(p,r,s,q)*MOe2ints(p,r,q,s) 
end do
end do
end do
end do

PRINT *, '2RDM 2e contribution',0.5D0*dtemp 
energy=0.5D0*dtemp

!get 1RDM from 2RDM and check we get 1e energy and total energy

do r=1,nbft
do q=1,nbft
do p=1,nbft
OneRDM(p,q)=OneRDM(p,q)+TwoRDM(p,r,r,q)
end do
end do
end do

OneRDM=OneRDM/((2.0D0*docc)-1.0D0)

dtemp=0.0D0

do q=1,nbft
do p=1,nbft
dtemp=dtemp+OneRDM(p,q)*MOe1ints(p,q)
end do
end do
PRINT *, '2RDM 1e contribution',dtemp

energy=energy+dtemp+ecore
PRINT *, 'total energy from 2RDM',energy





end subroutine readin2RDM_Arrow


subroutine readin2RDM(docc) !in MO basis This version reads in symmetry ordered and maps to energy ordered 
use AOintegral_arrays, only: nbft,MOe2ints,MOe1ints,ecore
use CIgrad_arrays, only: TwoRDM,OneRDM
implicit none
integer, intent(in) :: docc
integer info,p,q,r,s
integer i,j,k,m
double precision dtemp,energy
integer , allocatable:: mapToSymOrder(:)

allocate(TwoRDM(nbft,nbft,nbft,nbft))

TwoRDM=0.0D0
allocate(OneRDM(nbft,nbft)) !calculate from 2RDM
OneRDM=0.0D0


allocate(mapToSymorder(nbft))
OPEN(UNIT=11,FILE='MapToSymOrder.txt')
do i=1,nbft
READ(11,*) mapToSymOrder(i)
end do
CLOSE(11)




OPEN(UNIT=14,FILE='NonZero2RDM_MO_Coeffs')  !now this is from symmetry ordered orbitals and we need to map back to energy order

info=0
do while (info==0)  !loop until end of file
READ (14,*,IOSTAT=info)  p,q,r,s,dtemp

i=mapToSymOrder(p)
j=mapToSymOrder(q)
k=mapToSymOrder(r)
m=mapToSymOrder(s)


TwoRDM(i,k,m,j)=dtemp
end do
CLOSE(14)

!check we get 2electron energy
dtemp=0.0D0

do s=1,nbft
do r=1,nbft
do q=1,nbft
do p=1,nbft

dtemp=dtemp+TwoRDM(p,r,s,q)*MOe2ints(p,r,q,s) 
end do
end do
end do
end do

PRINT *, '2RDM 2e contribution',0.5D0*dtemp 
energy=0.5D0*dtemp

!get 1RDM from 2RDM and check we get 1e energy and total energy

do r=1,nbft
do q=1,nbft
do p=1,nbft
OneRDM(p,q)=OneRDM(p,q)+TwoRDM(p,r,r,q)
end do
end do
end do

OneRDM=OneRDM/((2.0D0*docc)-1.0D0)

dtemp=0.0D0

do q=1,nbft
do p=1,nbft
dtemp=dtemp+OneRDM(p,q)*MOe1ints(p,q)
end do
end do
PRINT *, '2RDM 1e contribution',dtemp

energy=energy+dtemp+ecore
PRINT *, 'total energy from 2RDM',energy


deallocate(mapToSymorder)

end subroutine




subroutine RHFgradients(docc)
use AOintegral_arrays, only: nbft
use RHFgrad_arrays
implicit none
integer, intent(in) :: docc
integer i,j
!PRINT *, 'in HF gradients subroutine',docc
allocate(dmE(nbft,nbft))

call createEweightDenmat(docc)  !creates dmE
allocate(Sgrad(3,nbft,nbft))
allocate(e1grad(3,nbft,nbft))
allocate(e2grad(3,nbft,nbft,nbft,nbft))


call readinSgrad() !reads in natoms, firstAO array, last AO array and Sgrad array

allocate(gradients(3,natoms))  ! first index is x,y,z
gradients=0.0D0

allocate(geometry(3,natoms))
allocate(charge(natoms))

allocate(rinvgrad(natoms,3,nbft,nbft))




call readinGeometry()

call VnnGrad()


call readinOneEgrad() 
call readin_rinvgrad() 



call contractSgrad_dmE() ! subtracts result from gradients array



call contractOneEgrad_denmat()



call readinTwoEgrad() 




call contractTwoEgrad_denmats()





deallocate(rinvgrad)
deallocate(dmE)
deallocate(Sgrad)
deallocate(firstAO)
deallocate(lastAO)
deallocate(gradients)
deallocate(e1grad)
deallocate(geometry)
deallocate(charge)
deallocate(e2grad)
end subroutine RHFgradients


subroutine contractTwoEgrad_denmats()
use RHFgrad_arrays, only:natoms,e2grad,firstAO,lastAO,gradients
use AOintegral_arrays, only: nbft,denmat
implicit none
integer i,j,q,p,coord,atom
double precision dtempJ,dtempK



do atom=1,natoms
do coord=1,3


!J contribution 
dtempJ=0.0D0
do j=1,nbft
do q=1,nbft
do i=1,nbft
do p=firstAO(atom),lastAO(atom)  
dtempJ=dtempJ-denmat(q,p)*denmat(i,j)*e2grad(coord,p,i,q,j) 
end do
end do
end do
end do



!K contribution
dtempK=0.0D0
do j=1,nbft
do q=1,nbft
do i=1,nbft
do p=firstAO(atom),lastAO(atom)  
dtempK=dtempK-denmat(q,p)*denmat(i,j)*e2grad(coord,p,q,i,j)  
end do
end do
end do
end do


gradients(coord,atom)=gradients(coord,atom)+(2.0D0*dtempJ)-dtempK

end do
end do



end subroutine contractTwoEgrad_denmats



subroutine contractOneEgrad_denmat()
use RHFgrad_arrays, only:natoms,e1grad,firstAO,lastAO,gradients,charge,rinvgrad
use AOintegral_arrays, only: nbft,denmat
implicit none
integer coord,atom,i,j
double precision dtemp,dtemp2


do atom=1,natoms
do coord=1,3

dtemp=0.0D0
do j=1,nbft
do i=firstAO(atom),lastAO(atom) 

dtemp=dtemp-e1grad(coord,i,j)*denmat(j,i) 
end do
end do

dtemp=2.0D0*dtemp ! transpose contribution as H^core hermitiian and denmat symmetric



!rinv contribution from deriv of H^{core}
dtemp2=0.0D0
do j=1,nbft 
do i=1,nbft

dtemp2=dtemp2+rinvgrad(atom,coord,i,j)*denmat(j,i)
end do
end do

dtemp2=2.0D0*dtemp2

dtemp=dtemp-charge(atom)*dtemp2

gradients(coord,atom)=gradients(coord,atom)+dtemp
end do
end do


end subroutine contractOneEgrad_denmat


subroutine Vnngrad() !nuclear nuclear contribution 
use RHFgrad_arrays, only:natoms,gradients,geometry,charge
implicit none
integer coord,atom,j
double precision dtemp

do coord=1,3
do atom=1,natoms

dtemp=0.0d0
do j=1,natoms
if(j.eq.atom) cycle
dtemp=dtemp+(charge(j)*(geometry(coord,j)-geometry(coord,atom))*&
( (geometry(1,j)-geometry(1,atom))**2 +(geometry(2,j)-geometry(2,atom))**2&
  +(geometry(3,j)-geometry(3,atom))**2  ) **(-1.5D0)  )

end do

gradients(coord,atom)=gradients(coord,atom)+charge(atom)*dtemp
end do
end do


end subroutine Vnngrad



subroutine contractSgrad_dmE()
use RHFgrad_arrays, only:natoms,Sgrad,firstAO,lastAO,dmE,gradients
use AOintegral_arrays, only: nbft
implicit none
integer coord,atom,i,j
double precision dtemp



do atom=1,natoms
do coord=1,3

dtemp=0.0D0
do j=1,nbft
do i=firstAO(atom),lastAO(atom) 



dtemp=dtemp-Sgrad(coord,i,j)*dmE(i,j)  
end do
end do

gradients(coord,atom)=gradients(coord,atom)-2.0D0*dtemp
end do
end do


!PRINT *,'in contractSgrad_dmE',-2.0D0*dtemp 

end subroutine contractSgrad_dmE


subroutine readin_rinvgrad() 
use RHFgrad_arrays, only:rinvgrad
implicit none
integer i,j,k,atom,info
double precision dtemp
rinvgrad=0.0D0

OPEN(UNIT=14,FILE='rinvGrad.txt') 


info=0
do while (info==0)  !loop until end of file
READ (14,*,IOSTAT=info)  dtemp,atom,i,j,k !atom is atom, i is x y or z, j and k label AOs
rinvgrad(atom,i,j,k)=dtemp
end do
CLOSE(14)




end subroutine readin_rinvgrad


subroutine readinTwoEgrad() 
use RHFgrad_arrays, only:e2grad
implicit none
integer i,j,k,l,coord
integer info
double precision dtemp

!no permutation symmetry is being used
e2grad=0.0D0

OPEN(UNIT=14,FILE='TwoEGrad.txt') 


info=0
do while (info==0)  !loop until end of file
READ (14,*,IOSTAT=info)  dtemp,coord,i,j,k,l  !coord is x y or z, j and k label AOs

e2grad(coord,i,k,j,l)=dtemp  
end do
CLOSE(14)



end subroutine readinTwoEgrad



subroutine readinOneEgrad() 
use RHFgrad_arrays, only:e1grad
implicit none
integer i,j,k,l
integer info
double precision dtemp

e1grad=0.0D0

OPEN(UNIT=14,FILE='OneEGrad.txt') 


info=0
do while (info==0)  !loop until end of file
READ (14,*,IOSTAT=info)  dtemp,i,j,k  !i is x y or z, j and k label AOs
e1grad(i,j,k)=dtemp
end do
CLOSE(14)



end subroutine readinOneEgrad



subroutine readinGeometry()
use RHFgrad_arrays, only: geometry,charge,natoms
implicit none
integer i

OPEN(unit=14,file='geometry.txt')

do i=1,natoms
READ (14,*) charge(i)
READ (14,*) geometry(1,i), geometry(2,i), geometry(3,i)
end do

CLOSE(14)

end subroutine readinGeometry

subroutine readinSgrad() !read in number of atoms, last AO for each atom, overlap integrals derivatives
use RHFgrad_arrays, only:natoms,Sgrad,firstAO,lastAO
implicit none
integer i,j,k,l
integer info
double precision dtemp

Sgrad=0.0D0

OPEN(UNIT=14,FILE='OverlapGrad.txt') 
READ (14,*) natoms
allocate(firstAO(natoms))
allocate(lastAO(natoms))

do i=1,natoms !read in first and last AO for each atom 
READ (14,*) j,k,firstAO(i),lastAO(i)

firstAO(i)=firstAO(i)+1 !fortran labels start for 1 not 0
end do
info=0
do while (info==0)  !loop until end of file
READ (14,*,IOSTAT=info)  dtemp,i,j,k  !i is x y or z, j and k label AOs
Sgrad(i,j,k)=dtemp
end do
CLOSE(14)



end subroutine readinSgrad




subroutine createEweightDenmat(docc) ! for HF gradients dmE is Fock eigenvalue weighted density matrix
use AOintegral_arrays, only: nbft,eigenFn,FockE
use RHFgrad_arrays, only:dmE
implicit none
integer, intent(in) :: docc
integer p,q,i

dmE=0.0D0

do i=1,docc
do q=1,nbft
do p=1,nbft
dmE(p,q)=dmE(p,q)+FockE(i)*eigenFn(p,i)*eigenFn(q,i)
end do
end do
end do

dmE=2.0D0*dmE


end subroutine createEweightDenmat


!!!!!!!!!!!!!!!!!! RHF subroutines below

subroutine calcMOtwoElec(writeMOints,docc) !two electron integrals in MO basis. Transformation cost nbft^5 
use AOintegral_arrays, only: e2ints,eigenFn,nbft,MOe2ints,MOe1ints,ecore,MOsym
implicit none
logical, intent(in) :: writeMOints
integer,intent(in) ::docc
integer i,j,k,m,p,q,r,s
double precision , allocatable:: Tran1(:,:,:,:)
integer info,pr,qs
 integer norb,nelec,ms2,isym
 integer , allocatable:: orbsym(:)
 NAMELIST /FCI/ norb,nelec,ms2,orbsym,isym

allocate(Tran1(nbft,nbft,nbft,nbft),stat=info)
if(info.ne.0) STOP "Error allocating 4 dim array Tran1 in calcMOtwoElec"
allocate(MOe2ints(nbft,nbft,nbft,nbft),stat=info)
if(info.ne.0) STOP "Error allocating 4 dim array MOe2ints in calcMOtwoElec"


Tran1=0.0D0
do m=1,nbft
do k=1,nbft
do j=1,nbft
do p=1,nbft
do i=1,nbft
Tran1(p,j,k,m)=Tran1(p,j,k,m)+&
(eigenFn(i,p)*e2ints(i,j,k,m))
end do
end do
end do
end do 
end do

MOe2ints=0.0D0 !this is used as temporary storage for transformation at this stage
do m=1,nbft
do k=1,nbft
do q=1,nbft
do j=1,nbft
do p=1,nbft
MOe2ints(p,q,k,m)=MOe2ints(p,q,k,m)+&
(eigenFn(j,q)*Tran1(p,j,k,m))
end do
end do
end do
end do 
end do


Tran1=0.0D0 !reuse this storage space
do m=1,nbft
do r=1,nbft
do k=1,nbft
do q=1,nbft
do p=1,nbft
Tran1(p,q,r,m)=Tran1(p,q,r,m)+&
(eigenFn(k,r)*MOe2ints(p,q,k,m)) 
end do
end do
end do
end do 
end do

!final step to get MO 2e integrals in Tran2

MOe2ints=0.0D0 ! reuse this space for final values
do s=1,nbft
do m=1,nbft
do r=1,nbft
do q=1,nbft
do p=1,nbft
MOe2ints(p,q,r,s)=MOe2ints(p,q,r,s)+&
(eigenFn(m,s)*Tran1(p,q,r,m))
end do
end do
end do
end do 
end do



if(writeMOints) THEN
PRINT *, 'Writing MO integrals to FCIDUMP' 


allocate(orbsym(nbft))






do i=1,nbft

orbsym(i)=MOsym(i)
end do





!!!!!!!!!!!!!!!!!!
OPEN(UNIT=14,FILE='FCIDUMP')
norb=nbft
nelec=2*docc 
ms2=0
isym=1
        write(14,'(1X,''&FCI NORB='',I3,'',NELEC='',I2,'',MS2='',I2,'','')') norb,nelec,0
          write (14,'(2X,''ORBSYM='',30(I1,'',''),6(/1X,40(I1,'','')))') (orbsym(i),i=1,norb)
          write (14,'(2X,''ISYM='',I1)') isym
          write (14,'(1X,''&END'')')


do p=1,nbft
do r=1,p
do q=1,nbft 
do s=1,q



i=p
k=r
j=q
m=s


pr=r+(p*(p-1)/2)
qs=s+(q*(q-1)/2)
if(pr.lt.qs) cycle

WRITE(14,'(E28.20,4I4)')  MOe2ints(i,j,k,m),p,r,q,s  
end do
end do
end do
end do

!1eints
do p=1,nbft
do q=1,p  


i=p
j=q


WRITE(14,'(E28.20,4I4)')  MOe1ints(i,j),p,q,0,0 
end do
end do

WRITE(14,'(E28.20,4I4)')  ecore,0,0,0,0
CLOSE(14)

deallocate(orbsym)

END IF




deallocate(Tran1)

end subroutine calcMOtwoElec


subroutine calcMOoneElec(writeMOints) !one electron integrals in MO basis 
use AOintegral_arrays, only: e1ints,eigenFn,nbft,MOe1ints
implicit none
logical, intent(in) :: writeMOints
integer i,j,p,q
double precision , allocatable:: Tran1(:,:)

allocate(MOe1ints(nbft,nbft))
allocate(Tran1(nbft,nbft))



MOe1ints=0.0D0
Tran1=0.0D0


do p=1,nbft
do j=1,nbft
do i=1,nbft
Tran1(p,j)=Tran1(p,j)+&
(eigenFn(i,p)*e1ints(i,j))
end do
end do
end do

do p=1,nbft
do q=1,nbft
do j=1,nbft
MOe1ints(p,q)=MOe1ints(p,q)+&
(eigenFn(j,q)*Tran1(p,j))
end do
end do
end do





deallocate(Tran1)
end subroutine calcMOoneElec



subroutine calcMOoverlap() 
use AOintegral_arrays, only: overlap,eigenFn,nbft
implicit none
integer i,j,p,q
double precision , allocatable:: MOoverlap(:,:)
double precision , allocatable:: Tran1(:,:)
allocate(MOoverlap(nbft,nbft))

allocate(Tran1(nbft,nbft))


MOoverlap=0.0D0

!naive nbft^4 cost transform works
do q=1,nbft
do p=1,nbft
do i=1,nbft
do j=1,nbft
MOoverlap(p,q)=MOoverlap(p,q)+&
(eigenFn(i,p)*eigenFn(j,q)*overlap(i,j))
end do
end do
end do
end do


MOoverlap=0.0D0
Tran1=0.0D0


do p=1,nbft
do j=1,nbft
do i=1,nbft
Tran1(p,j)=Tran1(p,j)+&
(eigenFn(i,p)*overlap(i,j))
end do
end do
end do

do p=1,nbft
do q=1,nbft
do j=1,nbft
MOoverlap(p,q)=MOoverlap(p,q)+&
(eigenFn(j,q)*Tran1(p,j))
end do
end do
end do


do i=1,nbft
PRINT *,(MOoverlap(i,j),j=1,nbft)   
end do

deallocate(MOoverlap)
deallocate(Tran1)
end subroutine calcMOoverlap


subroutine GetMaxDIISErrorMatrix(DIISerror,mat_i) ! get largest value of DIISErrorMatrix and return as DIISerror for convergence checking see P. Pulay J. Comp. Chem. 3(4), 556 (1982) 
use AOintegral_arrays, only: DIISerrorMatrix,nbft
implicit none
integer, intent(in) :: mat_i
double precision, intent(out) :: DIISerror
integer i,j

DIISerror=0.0D0
do j=1,nbft
do i=1,nbft
if(ABS(DIISerrorMatrix(i,j,mat_i)).gt.DIISerror) THEN
DIISerror=ABS(DIISerrorMatrix(i,j,mat_i))
END IF 
end do
end do

end subroutine GetMaxDIISErrorMatrix




subroutine CalcDIIScoeffs(Bdim)  !solve DIIS linear equations using B_ij=Tr(E_i E_j^{Transpose}) see P. Pulay J. Comp. Chem. 3(4), 556 (1982) 
use AOintegral_arrays, only: DIISerrorMatrix,nbft,DIIScoeffs,maxDIIS
implicit none
integer, intent(in) :: Bdim
double precision , allocatable:: A(:,:)
double precision , allocatable:: y(:)
integer, allocatable:: IPIV(:)
double precision , allocatable:: work(:)
integer i,j,k,r,Adim,Lwork,info

if(Bdim.gt.maxDIIS)  STOP 'Bdim.gt.maxDIIS in CalcDIIScoeffs'

Adim=Bdim+1

allocate(A(Adim,Adim))

allocate(y(Adim))

A=0.0D0

y=0.0D0

! Construct A and y
y(1)=-1.0D0

do i=2,Adim
A(1,i)=-1.0D0
A(i,1)=-1.0D0
end do

do i=2,Adim
do j=i,Adim ! only store upper triangle of A

do k=1,nbft
do r=1,nbft
A(i,j)=A(i,j)+DIISerrorMatrix(r,k,i-1)*DIISerrorMatrix(r,k,j-1)   ! A(i+1,j+1)=B_ij=Tr(E_i E_j^{Transpose})
end do
end do

end do
end do


!Solve Ax=y 

allocate(IPIV(Adim))
Lwork=64*Adim  
ALLOCATE(work(Lwork))

 Call dsysv('U',Adim,1,A,Adim,IPIV,&
 y,Adim,work,Lwork,info)  ! on output y is now the solution x
 
if(info.ne.0) STOP 'Problem with dsysv in CalcDIIScoeffs'
!PRINT *, 'optimum LWORK',work(1)

!copy solution in y to DIIScoeffs

do i=2,Adim
DIIScoeffs(i-1)=y(i)
end do


deallocate(work)
deallocate(IPIV)
deallocate(A)
deallocate(y)
end subroutine CalcDIIScoeffs




subroutine CalcDIISerrorMatrix(mat_i) !stores i'th error matrix of FDS-SDF  where F=F(D_i) D=D_i
use AOintegral_arrays, only: denmat,Fock,nbft,DIISerrorMatrix,overlap,maxDIIS
implicit none
integer, intent(in) :: mat_i
integer i,j
double precision alpha,beta
double precision , allocatable:: dtemp(:,:)
double precision , allocatable:: FDS(:,:)
double precision , allocatable:: SDF(:,:)

if(mat_i.gt.maxDIIS)  STOP 'mat_t.gt.maxDIIS in CalcDIISerrorMatrix'

allocate(dtemp(nbft,nbft))
allocate(FDS(nbft,nbft))
allocate(SDF(nbft,nbft))


!denmat*overlap=dtemp
dtemp=0.0D0
alpha=1.0D0  
beta=0.0D0
Call dgemm('N','N',nbft,nbft,nbft,alpha,&   
denmat,nbft,overlap,nbft,beta,dtemp,nbft) 

!Fock*dtemp=FDS
FDS=0.0D0
alpha=1.0D0  
beta=0.0D0
Call dgemm('N','N',nbft,nbft,nbft,alpha,&   
Fock,nbft,dtemp,nbft,beta,FDS,nbft)

!denmat*Fock=dtemp
dtemp=0.0D0
alpha=1.0D0  
beta=0.0D0
Call dgemm('N','N',nbft,nbft,nbft,alpha,&   
denmat,nbft,Fock,nbft,beta,dtemp,nbft) 

! overlap*dtemp=SDF
SDF=0.0D0
alpha=1.0D0  
beta=0.0D0
Call dgemm('N','N',nbft,nbft,nbft,alpha,&   
overlap,nbft,dtemp,nbft,beta,SDF,nbft) 

!store FDS-SDF
do j=1,nbft
do i=1,nbft
DIISerrorMatrix(i,j,mat_i)=FDS(i,j)-SDF(i,j)
end do
end do


deallocate(dtemp)
deallocate(FDS)
deallocate(SDF)
end subroutine CalcDIISerrorMatrix




subroutine HFenergy(energy) !0.5D0*Trace[(e1ints+F)*D] 
use AOintegral_arrays, only: e1ints,nbft,denmat,Fock
implicit none
double precision, intent(out) :: energy
integer i,j

energy=0.0D0

do j=1,nbft
do i=1,nbft
energy=energy+((e1ints(i,j)+Fock(i,j))*denmat(j,i))
end do
end do

energy=0.5D0*energy

end subroutine HFenergy


subroutine FormFock()
use AOintegral_arrays, only: e1ints,e2ints,nbft,denmat,Fock
implicit none
integer i,j,p,q
double precision dtempJ,dtempK

Fock=0.0D0


do j=1,nbft
do i=j,nbft 

dtempJ=0.0D0
do q=1,nbft
do p=1,nbft
dtempJ=dtempJ+&
(denmat(p,q)*(e2ints(p,i,q,j))) 
end do
end do

dtempK=0.0D0
do q=1,nbft
do p=1,nbft
dtempK=dtempK+&
(denmat(p,q)*(e2ints(p,q,i,j)))
end do
end do

Fock(i,j)=Fock(i,j)+dtempJ-0.5D0*dtempK
Fock(j,i)=Fock(i,j)
end do
end do

!Add 1e
Fock=Fock+e1ints

end subroutine FormFock








subroutine createDenmat(docc)
use AOintegral_arrays, only: nbft,eigenFn,denmat
implicit none
integer, intent(in) :: docc
double precision alpha,beta

denmat=0.0D0

alpha=2.0D0  
beta=0.0D0
Call dgemm('N','T',nbft,nbft,docc,alpha,&   
eigenFn,nbft,eigenFn,nbft,beta,denmat,nbft)  


end subroutine createDenmat




subroutine DIAGorthoFOCK(print_eigen) !diagonalize orthoFOCK matrix and transform eigenfunctions back to AO basis. Prints eigenvalues  if input is true
use AOintegral_arrays, only: OrthoFock,nbft,SMH,eigenFn,FockE
implicit none
integer info,i,j
double precision , allocatable:: work(:,:)
double precision , allocatable:: U(:,:)
double precision alpha,beta
logical, intent(in) :: print_eigen


ALLOCATE(U(nbft,nbft))
ALLOCATE(work(4*nbft,4*nbft))
!diagonalize orthoFOCK matrix
U=orthoFOCK ! as this will be overwritten with eigenvectors
 call dsyev('V','U',nbft,U,nbft,FockE,work,4*nbft,info) 
if(info.ne.0) STOP 'Problem with dsyev in DIAGorthoFOCK'

if(print_eigen) PRINT *,'eigenvalues of orthoFock',FockE
!Transform Orthoeigenfunctions U back to non-ortho AO  use dgemm for SMH*U=eigenFn

alpha=1.0D0
beta=0.0D0
eigenFn=0.0D0
Call dgemm('N','N',nbft,nbft,nbft,alpha,&
SMH,nbft,U,nbft,beta,eigenFn,nbft)



DEALLOCATE(U)
DEALLOCATE(work)
end subroutine DIAGorthoFOCK


subroutine TRAN2orthoFOCK !calculates  orthoFOCK= SMH^dagger FOCK SMH
use AOintegral_arrays, only: FOCK,OrthoFock,nbft,SMH
implicit none
double precision , allocatable:: dtemp(:,:)
double precision alpha,beta

ALLOCATE(dtemp(nbft,nbft))
OrthoFock=0.0D0
!DGEMM to perform SMH^dagger FOCK=dtemp
alpha=1.0D0
beta=0.0D0
Call dgemm('T','N',nbft,nbft,nbft,alpha,&
SMH,nbft,FOCK,nbft,beta,dtemp,nbft)
!DGEMM to perform dtemp*SMH=orthoFOCK
alpha=1.0D0
beta=0.0D0
Call dgemm('N','N',nbft,nbft,nbft,alpha,&
dtemp,nbft,SMH,nbft,beta,orthoFOCK,nbft)




DEALLOCATE(dtemp)
end subroutine TRAN2orthoFOCK





subroutine createSMH ! create S^-0.5
use AOintegral_arrays, only: SMH,overlap,nbft,SeigenTol
implicit none
integer info,i,j
double precision, allocatable:: diag(:)
double precision , allocatable:: work(:,:)
double precision , allocatable:: U(:,:)
double precision , allocatable:: dtemp(:,:)
double precision alpha,beta


SMH=0.0D0
ALLOCATE(diag(nbft))
ALLOCATE(U(nbft,nbft))
ALLOCATE(dtemp(nbft,nbft))
ALLOCATE(work(4*nbft,4*nbft))
!diagonalize overlap matrix
U=overlap ! as this will be overwritten with eigenvectors
 call dsyev('V','U',nbft,U,nbft,diag,work,4*nbft,info) 
if(info.ne.0) STOP 'Problem with dsyev in createSMH'

!first diagonal element is smallest so check if this is less than SeigenTol
if(diag(1).lt.SeigenTol) STOP 'Eigenvalues of S less than SeigenTol in createSMH'


!now turn diag into 1/DSQRT(diag) to create SMH=U 1/SQRT(diag) U^dagger

do i=1,nbft
diag(i)=1.0D0/DSQRT(diag(i))
end do

!Do  U*diag=dtemp myself

do j=1,nbft
do i=1,nbft
dtemp(i,j)=U(i,j)*diag(j)
end do
end do

!DGEMM to perform dtemp*U^dagger=SMH
alpha=1.0D0
beta=0.0D0
Call dgemm('N','T',nbft,nbft,nbft,alpha,&
dtemp,nbft,U,nbft,beta,SMH,nbft)



 


DEALLOCATE(diag)
DEALLOCATE(U)
DEALLOCATE(work)

end subroutine createSMH 








subroutine ReadInAOintsPySCF
use AOintegral_arrays, only: ecore,overlap,nbft,e1ints,e2ints
implicit none
integer  ntitle, nsym, ninfo, nenrgy, nmap
integer nskip,nbpsy(8) ! max 8 symmetries
integer num, lrecl, idum1, itypea, itypeb, idum2, last
integer i,i1,j1,i2,j2,k2,l2,j,k,l
integer iskip


double precision eint

OPEN(UNIT=10,FILE='Overlap.txt')

READ (10,*) nbft
READ (10,*) ecore  
  
  ALLOCATE(overlap(nbft,nbft))
  overlap=0.0D0
  
  do i=1,nbft
  do j=i,nbft
  READ(10,*) eint,i1,j1
           if(j1.gt.i1) STOP 'Error get_int: overlap ints '
             overlap(i1,j1) = eint
             overlap(j1,i1)=eint
           end do
           end do
  
  CLOSE(10)
  

  
  
  
  OPEN(UNIT=10,FILE='OneEints.txt')
  ALLOCATE(e1ints(nbft,nbft))
  e1ints=0.0D0
  
    do i=1,nbft
  do j=i,nbft
  READ(10,*) eint,i1,j1
           if(j1.gt.i1) STOP 'Error get_int:  OneEints '
            e1ints(i1,j1) = eint
                  e1ints(j1,i1) = e1ints(i1,j1) 
           end do
           end do
  
  CLOSE(10)
  
  
    OPEN(UNIT=10,FILE='TwoEints.txt')
    ALLOCATE(e2ints(nbft,nbft,nbft,nbft))
    e2ints=0.0D0
  
  do i=1,nbft
  do j=i,nbft
  do k=1,nbft
  do l=k,nbft
   read(10,*) eint, i2, j2, k2, l2
           if(j2.gt.i2 .or. l2.gt.k2) STOP 'Error get_int: 2e ints'
      
        	 e2ints(i2,k2,j2,l2)=eint  
     e2ints(j2,k2,i2,l2)=eint        
     e2ints(i2,l2,j2,k2)=eint   
     e2ints(j2,l2,i2,k2)=eint       
          
     e2ints(k2,i2,l2,j2)=eint  
     e2ints(l2,i2,k2,j2)=eint                    
     e2ints(k2,j2,l2,i2)=eint     
     e2ints(l2,j2,k2,i2)=eint                   
end do
end do
end do
end do
  
  
  CLOSE(10)
   



END subroutine ReadInAOintsPySCF



