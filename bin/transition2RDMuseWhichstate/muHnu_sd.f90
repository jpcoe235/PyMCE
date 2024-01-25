subroutine muHnu_sd(ici,jci,e,n_2p,s_overlap)
  use commonarrays, only: nbft, icij, nword, i_sx2
  use dyn_par
  implicit real*8    (a-h,o-z)
  integer           ijsingles(iword), ijdoubles(iword),ij_sd(iword)
  !common  /aodat/   nsym, nbft, nbpsy(irmax)
  !common  /config/  icij(2,iword,maxc), nword
  !common  /occupy/  ntotal, n_alpha, n_beta, i_sx2
  integer newdiff,mytemp
  e = 0.0d0
  s_overlap = 0.0d0

!!!!!!!!!!!!!
newdiff=0
do n=1,nword
mytemp=IEOR(icij(1,n,ici),icij(1,n,jci))

newdiff=newdiff+POPCNT(mytemp) ! calcs number of bits set to 1 as we used xor (IEOR) bitwise that must be - note one difference adds two to newdiff as there are two places where the bits will be set to 1 

mytemp=IEOR(icij(2,n,ici),icij(2,n,jci))

newdiff=newdiff+POPCNT(mytemp)

end do

 if (newdiff.gt.4) return !more than two differences so matrix element is zero

 call reorder_sd(ici,jci,ndiff,idiff1,idiff2,ep)

  if(ndiff.gt.2) return

  !     Slater Condon Harris rules
  if( ndiff .eq. 0) then
     s_overlap = 1.0D0  !ep=1 if ndiff=0 as duplicate dictionary ordered SDs ahve been removed
     call slater0(e,i_am_mu,i_am_nu,nu_doubly,kck,n_2p,ep)       ! no orbital differences
  elseif( ndiff .eq. 1) then
     call slater1(e,i_am_mu,i_am_nu,nu_doubly,idiff1,kck,n_2p,ep) 		! one orbital difference
  elseif( ndiff .eq. 2) then
     call slater2(e,i_am_mu,i_am_nu,nu_doubly,idiff1,idiff2,kck,n_2p,ep) 		! two orbital differences
  else
     write(*,*) 'ndiff', ndiff
     STOP 'muHnu: no conditions were met'
  endif

      return
      end
