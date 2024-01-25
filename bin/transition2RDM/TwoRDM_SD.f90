subroutine TwoRDM_SD(length,TwoRDM_e) !length of wavefunction, energy from 2RDM (without nuclear contribution)
use commonarrays, only: nbft,ntotal,c,list,icij,nword  ! basis functions, electrons, coefficients, list of orbitals in SD
implicit none
integer length
double precision, allocatable :: SpinFree2RDM(:,:,:,:)
integer porb,rorb,sorb,qorb
integer spins,spinq,spinp,spinr
double precision dtemp,ep,eg
integer ici,jci,ndiff,idiff1,idiff2
integer i,i1,i2,l,l2,k,k2
double precision Eone_e,Etwo_e,TwoRDM_e
integer newdiff,mytemp,n

ALLOCATE(SpinFree2RDM(nbft,nbft,nbft,nbft))

  !set to zero (a while ago Array=0.0D0 did not work for some compilers)
 do porb=1,nbft
 do rorb=1,nbft 
 do sorb=1,nbft
 do qorb=1,nbft
 SpinFree2RDM(porb,rorb,sorb,qorb)=0.0D0
 end do
 end do
 end do
 end do

dtemp=0.0D0
do ici=1,length
do jci=1,length

!!!!!!!!!!!!!
newdiff=0
do n=1,nword
mytemp=IEOR(icij(1,n,ici),icij(1,n,jci))

newdiff=newdiff+POPCNT(mytemp) ! calcs number of bits set to 1 as we used xor (IEOR) bitwise that must be - note one difference adds two to newdiff as there are two places where the bits will be set to 1 

mytemp=IEOR(icij(2,n,ici),icij(2,n,jci))

newdiff=newdiff+POPCNT(mytemp)

end do

 if (newdiff.gt.4) cycle !more than two differences so matrix element is zero

!!!!!!!!!!!!!!!!!!!!

          
         call reorder_sd(ici,jci,ndiff,idiff1,idiff2,ep)
!         if(ndiff.gt.2) cycle !can't be true as we have already verified not gt 2

if(ndiff.eq.2) THEN
   !Case 1

     sorb=list(2,idiff1)
     qorb=list(2,idiff2)
     porb=list(1,idiff2)
     rorb=list(1,idiff1)
      eg=1.0D0

 call getSpinAndSpatial4(sorb,qorb,porb,rorb,spins,spinq,spinp,spinr)
if(spins.eq.spinr) THEN ! check spins are allowed
if(spinq.eq.spinp) THEN

dtemp=SpinFree2RDM(porb,rorb,sorb,qorb)
  dtemp=dtemp+(c(ici)*c(jci)*ep*eg)
SpinFree2RDM(porb,rorb,sorb,qorb)=dtemp
end if
end if
    
      !Case 2 Swap p and r introduces negative sign
 sorb=list(2,idiff1)
     qorb=list(2,idiff2)
     rorb=list(1,idiff2)
     porb=list(1,idiff1)
      eg=-1.0D0
 call getSpinAndSpatial4(sorb,qorb,porb,rorb,spins,spinq,spinp,spinr)
if(spins.eq.spinr) THEN ! check spins are allowed
if(spinq.eq.spinp) THEN

dtemp=SpinFree2RDM(porb,rorb,sorb,qorb)
  dtemp=dtemp+(c(ici)*c(jci)*ep*eg)
SpinFree2RDM(porb,rorb,sorb,qorb)=dtemp
end if
end if
!Case 3 from Case 1 swap s and q to give negative sign
     qorb=list(2,idiff1)
     sorb=list(2,idiff2)
     porb=list(1,idiff2)
     rorb=list(1,idiff1)
      eg=-1.0D0
 call getSpinAndSpatial4(sorb,qorb,porb,rorb,spins,spinq,spinp,spinr)
if(spins.eq.spinr) THEN ! check spins are allowed
if(spinq.eq.spinp) THEN

dtemp=SpinFree2RDM(porb,rorb,sorb,qorb)
  dtemp=dtemp+(c(ici)*c(jci)*ep*eg)
SpinFree2RDM(porb,rorb,sorb,qorb)=dtemp
end if
end if

!Case 4 from Case 1 swap s and q then swap p and r so no sign change
     qorb=list(2,idiff1)
     sorb=list(2,idiff2)
     rorb=list(1,idiff2)
     porb=list(1,idiff1)
      eg=1.0D0
 call getSpinAndSpatial4(sorb,qorb,porb,rorb,spins,spinq,spinp,spinr)
if(spins.eq.spinr) THEN ! check spins are allowed
if(spinq.eq.spinp) THEN

dtemp=SpinFree2RDM(porb,rorb,sorb,qorb)
  dtemp=dtemp+(c(ici)*c(jci)*ep*eg)
SpinFree2RDM(porb,rorb,sorb,qorb)=dtemp
end if
end if

cycle
end if ! end of ndiff.eq.2


if(ndiff.eq.1) THEN ! 

!     Case1 
     qorb=list(2,idiff1)
     porb=list(1,idiff1)
      eg=1.0D0
        !get spins and spatial for these two
        call getSpinAndSpatial(qorb,porb,spinq,spinp)

       do i=1,ntotal !loop over other occupied
       if(i.eq.idiff1) cycle !can't be one of the single differences
             
       sorb=list(2,i)
       rorb=list(1,i)

       call getSpinAndSpatial(sorb,rorb,spins,spinr)

       if(spins.eq.spinr) THEN ! check spins are allowed
if(spinq.eq.spinp) THEN

dtemp=SpinFree2RDM(porb,rorb,sorb,qorb)
  dtemp=dtemp+(c(ici)*c(jci)*ep*eg)
SpinFree2RDM(porb,rorb,sorb,qorb)=dtemp
end if
end if
         
        end do  ! end of loop over occupied Case 1

 !Case 2 now q and r have the different orbitals - causes sign change to make both come before other occupied
qorb=list(2,idiff1)
     rorb=list(1,idiff1)
      eg=-1.0D0

      call getSpinAndSpatial(qorb,rorb,spinq,spinr)

       do i=1,ntotal !loop over other occupied
       if(i.eq.idiff1) cycle !can't be one of the single differences
             
       sorb=list(2,i)
       porb=list(1,i)

       call getSpinAndSpatial(sorb,porb,spins,spinp)

 

       if(spins.eq.spinr) THEN ! check spins are allowed
if(spinq.eq.spinp) THEN

dtemp=SpinFree2RDM(porb,rorb,sorb,qorb)
  dtemp=dtemp+(c(ici)*c(jci)*ep*eg)
SpinFree2RDM(porb,rorb,sorb,qorb)=dtemp
end if
end if
         
        end do  ! end of loop over occupied Case 2

!Case 3 s and p are the different ones - sign change so both operators to come before occupied when acting on their det (left det for p and r, right det for s and q)
sorb=list(2,idiff1)
     porb=list(1,idiff1)
      eg=-1.0D0

      call getSpinAndSpatial(sorb,porb,spins,spinp)

       do i=1,ntotal !loop over other occupied
       if(i.eq.idiff1) cycle !can't be one of the single differences
             
       qorb=list(2,i)
       rorb=list(1,i)
       call getSpinAndSpatial(qorb,rorb,spinq,spinr)

       if(spins.eq.spinr) THEN ! check spins are allowed
if(spinq.eq.spinp) THEN

dtemp=SpinFree2RDM(porb,rorb,sorb,qorb)
  dtemp=dtemp+(c(ici)*c(jci)*ep*eg)
SpinFree2RDM(porb,rorb,sorb,qorb)=dtemp
end if
end if
         
        end do  ! end of loop over occupied Case 3  
!     Case4  s and r are the differences so signs of moving through occupied will cancel
     sorb=list(2,idiff1)
     rorb=list(1,idiff1)
      eg=1.0D0
         call getSpinAndSpatial(sorb,rorb,spins,spinr)


       do i=1,ntotal !loop over other occupied
       if(i.eq.idiff1) cycle !can't be one of the single differences
             
       qorb=list(2,i)
       porb=list(1,i)
          call getSpinAndSpatial(qorb,porb,spinq,spinp)

    

       if(spins.eq.spinr) THEN ! check spins are allowed
if(spinq.eq.spinp) THEN

dtemp=SpinFree2RDM(porb,rorb,sorb,qorb)
  dtemp=dtemp+(c(ici)*c(jci)*ep*eg)
SpinFree2RDM(porb,rorb,sorb,qorb)=dtemp
end if
end if
         
        end do  ! end of loop over occupied Case 4


cycle
end if  ! end of ndiff.eq.1

if(ndiff.eq.0) THEN !Now spin free derivation see notes

      do i1=1,ntotal  !   all possible pairs from orbital list. 
      do i2=1,ntotal  !was i1  
      if(i1.eq.i2) cycle
         
        sorb=list(2,i1)
        qorb=list(2,i2)
         call getSpinAndSpatial(sorb,qorb,spins,spinq)

  !Case 1 p is same as s and r is q
  if(spinq.eq.spins)  THEN
 porb=sorb
 rorb=qorb
 eg=-1.0D0
dtemp=SpinFree2RDM(porb,rorb,sorb,qorb)
  dtemp=dtemp+(c(ici)*c(jci)*ep*eg)
SpinFree2RDM(porb,rorb,sorb,qorb)=dtemp
  END IF  !case 1        

!Case 2 p is same as q and r is same as s
!All s and q spins valid for this contribution
 porb=qorb
 rorb=sorb
 eg=1.0D0
dtemp=SpinFree2RDM(porb,rorb,sorb,qorb)
  dtemp=dtemp+(c(ici)*c(jci)*ep*eg)
SpinFree2RDM(porb,rorb,sorb,qorb)=dtemp


      end do
      end do
    end if ! ndiff=0



end do !end of ici loop
end do !end of jci loop



!Write out 2RDM

OPEN(UNIT=15,FILE='NonZero2RDM_MO_Coeffs')
do l=1,nbft
do l2=1,nbft
do k=1,nbft
do k2=1,nbft


if (ABS(SpinFree2RDM(l,l2,k,k2)).eq.0.0D0) cycle
!WRITE(15,*) l,l2,k,k2,SpinFree2RDM(l,l2,k,k2)

WRITE(15,*) l,k2,l2,k,SpinFree2RDM(l,l2,k,k2) ! notation from arxiv 1809.09058
end do
end do
end do
end do
CLOSE(15)

 call One_eEnergy(SpinFree2RDM,Eone_e)
 PRINT *, '1e Energy from 2RDM subroutine',Eone_e
 call Two_eEnergy(SpinFree2RDM,Etwo_e)
 PRINT *, '2e Energy from 2RDM subroutine',Etwo_e

TwoRDM_E=Eone_e+Etwo_e



DEALLOCATE(SpinFree2RDM)
end subroutine


subroutine Two_eEnergy(SpinFree2RDM,Etwo_e)
use commonarrays, only: nbft,e2ints,ipoint 
implicit none
double precision SpinFree2RDM(nbft,nbft,nbft,nbft),Etwo_e
integer iorb,jorb,korb,lorb,ik,jl,ikjl
!2RDM energy calc
Etwo_e=0.0D0
do iorb=1,nbft
do jorb=1,nbft
do korb=1,nbft
do lorb=1,nbft

if(korb.gt.iorb) THEN
          ik = ipoint(korb) + iorb  
          ELSE
          ik = ipoint(iorb) + korb  
          END IF
           

           if(lorb.gt.jorb) THEN
          jl = ipoint(lorb) + jorb  
          ELSE
          jl = ipoint(jorb) + lorb
          END IF
         
         if(jl.gt.ik) THEN
         ikjl = ipoint(jl) + ik 
         ELSE
         ikjl = ipoint(ik) + jl 
         END IF

         

Etwo_e=Etwo_e+(SpinFree2RDM(iorb,jorb,lorb,korb)*e2ints(ikjl)) ! note swapped last two terms as 2RDM_p,r,s,q, goes with <pr|qs> integral
end do
end do
end do
end do


Etwo_e=0.5D0*Etwo_e!remember we have 0.5 in the energy term not density matrix
end subroutine



subroutine One_eEnergy(SpinFree2RDM,Eone_e)
use commonarrays, only: nbft,ntotal,e1ints,ipoint 
implicit none
double precision SpinFree2RDM(nbft,nbft,nbft,nbft),Eone_e
integer korb,jorb,iorb,kk

!Energy test get 1RDM from 2RDM contribution
Eone_e=0.0D0

do korb=1,nbft
do iorb=1,nbft
do jorb=1,nbft

 if(iorb.gt.jorb) THEN
 kk = ipoint(iorb) + jorb
 ELSE
 kk = ipoint(jorb) + iorb
  END IF
!PRINT *,iorb,jorb,kk,e1ints(kk)
Eone_e=Eone_e+(SpinFree2RDM(iorb,korb,korb,jorb)*e1ints(kk))
      
end do
end do
end do

Eone_e=Eone_e/(1.0D0*(ntotal-1))



end subroutine



subroutine getSpinAndSpatial(sorb,qorb,spins,spinq)
use commonarrays, only: nbft
implicit none
integer sorb,qorb,porb,rorb
integer spins,spinq,spinp,spinr
! get spins and spatial orbs
        spins=1
        if(sorb.gt.nbft) THEN
         sorb = sorb - nbft
         spins=-1
         END IF

        spinq=1
        if(qorb.gt.nbft) THEN
        qorb = qorb - nbft
        spinq=-1
        END IF

     
end subroutine



subroutine getSpinAndSpatial4(sorb,qorb,porb,rorb,spins,spinq,spinp,spinr)
use commonarrays, only: nbft
implicit none
integer sorb,qorb,porb,rorb
integer spins,spinq,spinp,spinr
! get spins and spatial orbs
        spins=1
        if(sorb.gt.nbft) THEN
         sorb = sorb - nbft
         spins=-1
         END IF

        spinq=1
        if(qorb.gt.nbft) THEN
        qorb = qorb - nbft
        spinq=-1
        END IF

       spinp=1
        if(porb.gt.nbft) THEN
        porb = porb - nbft
        spinp=-1
        END IF

           spinr=1
        if(rorb.gt.nbft) THEN
        rorb = rorb - nbft
        spinr=-1
        END IF
end subroutine
