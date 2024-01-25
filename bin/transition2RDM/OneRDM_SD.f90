
subroutine OneRDM_SD(length) !length of wavefunction
use commonarrays, only: nbft,ntotal,c,list,e1ints,ipoint  ! basis functions, electrons, coefficients, list of orbitals in SD
implicit none
integer length
integer ici,jci,ndiff,idiff1,idiff2
 integer l,k
 double precision ep 
double precision, allocatable :: OneRDMCoeffs(:,:)
integer iorb,jorb,kk
double precision e

ALLOCATE(OneRDMCoeffs(nbft,nbft))
  
   !set to zero (think OneRDMCoeffs=0.0D0 did not work for some compilers)
 do l=1,nbft
 do k=1,nbft 
 OneRDMCoeffs(l,k)=0.0D0
 end do
 end do


do ici=1,length
do jci=1,length  

 
 call reorder_sd(ici,jci,ndiff,idiff1,idiff2,ep) ! Returns list of orbitals, sign from putting in maximal coincidence and differences 
 if(ndiff.gt.1) cycle

    if(ndiff.eq.1) THEN
       l=list(1,idiff1)
       k=list(2,idiff1)
       !alpha and beta fixed so spins must be the same for one difference
       if(l.gt.nbft) l=l-nbft

       if(k.gt.nbft) k=k-nbft
     
 OneRDMCoeffs(l,k)=OneRDMCoeffs(l,k)+(ep*c(ici)*c(jci)) 

    ELSE !ndiff=0

      do l=1,ntotal 
      k=list(1,l)
      if (k.gt.nbft) k=k-nbft
 OneRDMCoeffs(k,k)=OneRDMCoeffs(k,k)+(c(ici)*c(jci))
      end do

    END IF


end do
end do


!output
OPEN(UNIT=15,FILE='OneRDM.txt')
do l=1,nbft
do k=l,nbft  !OneRDM is symmetric (transition 1RDM may not be)
if (OneRDMCoeffs(l,k).eq.0.0D0) cycle ! only store non zero
WRITE(15,*) l,k, OneRDMCoeffs(l,k)
end do
end do
CLOSE(15)



!Energy test 1RDM contribution
e=0.0D0

do iorb=1,nbft
do jorb=1,nbft

 if(iorb.gt.jorb) THEN
 kk = ipoint(iorb) + jorb
 ELSE
 kk = ipoint(jorb) + iorb
  END IF
!PRINT *,iorb,jorb,kk,e1ints(kk)
e=e+(OneRDMCoeffs(iorb,jorb)*e1ints(kk))
      
end do
end do
PRINT *, 'Energy from 1RDM',e
DEALLOCATE(OneRDMCoeffs)
end subroutine

