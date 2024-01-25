!only terms with c(ici)*c(jci) will have to change so should be able to do at single 2RDM cost but more memory
! see if we can store 2RDMs for all pairs of interest not including same state
!seems reasonable but let's calc T1RDM from T2RDM and see if it agrees with T1RDM computed directly  -- this is fine for state 1 to state 2 calculated directly  so comment out T1RDM calc now
!also ok with 3 states giving 3 pairs 1,2; 1,3; 2,3   - if we run out of memory will have to create slower version where we write out T2RDM before calculating next



subroutine TransitionTwoRDM_SD(length,ieig) !length of wavefunction and number of wavefunction states that have been calculated  ! calls some subroutines in TwoRDM_SD
use commonarrays, only: nbft,ntotal,ctemp,list,icij,nword  ! basis functions, electrons, coefficients, list of orbitals in SD
implicit none
integer length,ieig
double precision, allocatable :: SpinFreeT2RDM(:,:,:,:,:) !SpinFreeT2RDM   last variable is state pair label
integer porb,rorb,sorb,qorb
integer spins,spinq,spinp,spinr
double precision dtemp,ep,eg
integer ici,jci,ndiff,idiff1,idiff2
integer i,i1,i2,l,l2,k,k2
double precision Eone_e,Etwo_e,TwoRDM_e
integer newdiff,mytemp,n
integer is,is2,slabel,statepairs
CHARACTER (LEN=30)  :: f1 
CHARACTER (LEN=30)  :: f2

PRINT *, 'in TransitionTwoRDM',ieig
statepairs=((ieig-1)*ieig)/2
PRINT  *,statepairs,'pairs'




ALLOCATE(SpinFreeT2RDM(nbft,nbft,nbft,nbft,statepairs))

  !set to zero (a while ago Array=0.0D0 did not work for some compilers)
 do porb=1,nbft
 do rorb=1,nbft 
 do sorb=1,nbft
 do qorb=1,nbft
 do slabel=1,statepairs
 
 SpinFreeT2RDM(porb,rorb,sorb,qorb,slabel)=0.0D0
 end do
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

slabel=0
do is=1,ieig
do is2=is+1,ieig
slabel=slabel+1

dtemp=SpinFreeT2RDM(porb,rorb,sorb,qorb,slabel)
  dtemp=dtemp+(ctemp(ici,is)*ctemp(jci,is2)*ep*eg)
SpinFreeT2RDM(porb,rorb,sorb,qorb,slabel)=dtemp

end do
end do

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

slabel=0
do is=1,ieig
do is2=is+1,ieig
slabel=slabel+1

dtemp=SpinFreeT2RDM(porb,rorb,sorb,qorb,slabel)
  dtemp=dtemp+(ctemp(ici,is)*ctemp(jci,is2)*ep*eg)
SpinFreeT2RDM(porb,rorb,sorb,qorb,slabel)=dtemp

end do
end do

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

slabel=0
do is=1,ieig
do is2=is+1,ieig
slabel=slabel+1

dtemp=SpinFreeT2RDM(porb,rorb,sorb,qorb,slabel)
  dtemp=dtemp+(ctemp(ici,is)*ctemp(jci,is2)*ep*eg)
SpinFreeT2RDM(porb,rorb,sorb,qorb,slabel)=dtemp

end do
end do

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

slabel=0
do is=1,ieig
do is2=is+1,ieig
slabel=slabel+1

dtemp=SpinFreeT2RDM(porb,rorb,sorb,qorb,slabel)
  dtemp=dtemp+(ctemp(ici,is)*ctemp(jci,is2)*ep*eg)
SpinFreeT2RDM(porb,rorb,sorb,qorb,slabel)=dtemp

end do
end do

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

slabel=0
do is=1,ieig
do is2=is+1,ieig
slabel=slabel+1

dtemp=SpinFreeT2RDM(porb,rorb,sorb,qorb,slabel)
  dtemp=dtemp+(ctemp(ici,is)*ctemp(jci,is2)*ep*eg)
SpinFreeT2RDM(porb,rorb,sorb,qorb,slabel)=dtemp

end do
end do

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

slabel=0
do is=1,ieig
do is2=is+1,ieig
slabel=slabel+1

dtemp=SpinFreeT2RDM(porb,rorb,sorb,qorb,slabel)
  dtemp=dtemp+(ctemp(ici,is)*ctemp(jci,is2)*ep*eg)
SpinFreeT2RDM(porb,rorb,sorb,qorb,slabel)=dtemp

end do
end do

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

slabel=0
do is=1,ieig
do is2=is+1,ieig
slabel=slabel+1

dtemp=SpinFreeT2RDM(porb,rorb,sorb,qorb,slabel)
  dtemp=dtemp+(ctemp(ici,is)*ctemp(jci,is2)*ep*eg)
SpinFreeT2RDM(porb,rorb,sorb,qorb,slabel)=dtemp

end do
end do

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

slabel=0
do is=1,ieig
do is2=is+1,ieig
slabel=slabel+1

dtemp=SpinFreeT2RDM(porb,rorb,sorb,qorb,slabel)
  dtemp=dtemp+(ctemp(ici,is)*ctemp(jci,is2)*ep*eg)
SpinFreeT2RDM(porb,rorb,sorb,qorb,slabel)=dtemp

end do
end do

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

slabel=0
do is=1,ieig
do is2=is+1,ieig
slabel=slabel+1

dtemp=SpinFreeT2RDM(porb,rorb,sorb,qorb,slabel)
  dtemp=dtemp+(ctemp(ici,is)*ctemp(jci,is2)*ep*eg)
SpinFreeT2RDM(porb,rorb,sorb,qorb,slabel)=dtemp

end do
end do

  END IF  !case 1        

!Case 2 p is same as q and r is same as s
!All s and q spins valid for this contribution
 porb=qorb
 rorb=sorb
 eg=1.0D0

slabel=0
do is=1,ieig
do is2=is+1,ieig
slabel=slabel+1

dtemp=SpinFreeT2RDM(porb,rorb,sorb,qorb,slabel)
  dtemp=dtemp+(ctemp(ici,is)*ctemp(jci,is2)*ep*eg)
SpinFreeT2RDM(porb,rorb,sorb,qorb,slabel)=dtemp

end do
end do

      end do
      end do
    end if ! ndiff=0



end do !end of ici loop
end do !end of jci loop



!Write out T2RDMs
slabel=0
do is=1,ieig
do is2=is+1,ieig

slabel=slabel+1
write(f1,'(i2)') is
 write(f2,'(i2)') is2
            f1='state'//trim(adjustl(f1))
            f2='state'//trim(adjustl(f2))
            f1=trim(adjustl(f1))//trim(adjustl(f2))
            f1='Tran2RDM_'//trim(adjustl(f1))



OPEN(UNIT=15,FILE=f1)
do l=1,nbft
do l2=1,nbft
do k=1,nbft
do k2=1,nbft


if (ABS(SpinFreeT2RDM(l,l2,k,k2,slabel)).eq.0.0D0) cycle
!WRITE(15,*) l,l2,k,k2,SpinFree2RDM(l,l2,k,k2)

WRITE(15,*) l,k2,l2,k,SpinFreeT2RDM(l,l2,k,k2,slabel) ! notation from arxiv 1809.09058
end do
end do
end do
end do
CLOSE(15)


! call calc_T1RDM(SpinFreeT2RDM,ieig)

end do
end do



DEALLOCATE(SpinFreeT2RDM)
end subroutine TransitionTwoRDM_SD


subroutine calc_T1RDM(SpinFreeT2RDM,ieig)
use commonarrays, only: nbft,ntotal
implicit none
double precision SpinFreeT2RDM(nbft,nbft,nbft,nbft,ieig)
integer korb,jorb,iorb,kk,ieig,slabel
double precision, allocatable :: TOneRDM(:,:) 
allocate(TOneRDM(nbft,nbft)) 


TOneRDM=0.0D0
do korb=1,nbft
do iorb=1,nbft
do jorb=1,nbft

TOneRDM(iorb,jorb)=TOneRDM(iorb,jorb)+SpinFreeT2RDM(iorb,korb,korb,jorb,1)   ! this is state 1 to 2 have to change last value to get other transitions
!PRINT *,iorb,jorb,kk,e1ints(kk)
!Eone_e=Eone_e+(SpinFree2RDM(iorb,korb,korb,jorb)*e1ints(kk))
      
end do
end do
end do

TOneRDM=TOneRDM/(1.0D0*(ntotal-1))

PRINT *, 'Transition 1RDM for state 1 to 2 from Transition 2RDM'
do iorb=1,nbft
do jorb=1,nbft
if(ABS(TOneRDM(iorb,jorb)).eq.0.0D0) cycle

PRINT *, iorb,jorb,TOneRDM(iorb,jorb)

end do
end do

deallocate(TOneRDM)
 
end subroutine calc_T1RDM






