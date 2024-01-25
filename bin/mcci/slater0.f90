
      subroutine slater0(e,i_am_mu,i_am_nu,nu_doubly,kck,n_2p,ep)
      use commonarrays, only: nbft, ntotal, i_sx2, e1ints, e2ints, ipoint, list
      use dyn_par
      implicit real*8    (a-h,o-z)
      integer ispin,jspin
   

!     diagonal Slater rules

!     zero energy
      e = 0.0d0
!     sum one particle contributions for doubly and singly occupied orbitals
      do n = 1, ntotal
	 korb = list(1,n) 
	 if( korb .gt. nbft) korb = korb - nbft
         kk = ipoint(korb) + korb
         e = e + e1ints(kk)
      enddo 
      
!   two particle contributions  - set exc to zero if i and j (i is korb and j lorb) spins don't match - check i j loop is ok as not orbitals but if we swap i and j should be ok - just avoid double counting
       do i=1,ntotal
        do j=1,i-1
        korb = list(1,i)
        lorb = list(1,j)
! check spins are the same for exchange contribution        
         ispin=1
         jspin=1
       
         if(korb.gt.nbft) THEN
         korb = korb - nbft
         ispin=-1
         END IF

         if(lorb.gt.nbft) THEN
         lorb=lorb-nbft   
         jspin=-1
         END IF


         kk = ipoint(korb) + korb
         ll = ipoint(lorb) + lorb
         if(ll.gt.kk) THEN
         kkll = ipoint(ll) + kk
         ELSE
         kkll = ipoint(kk) + ll
         END IF
         e = e + e2ints(kkll)
! exchange term
         if(ispin.EQ.jspin) THEN
        
         if(korb.gt.lorb) THEN
            kl = ipoint(korb) + lorb 
          ELSE
           kl = ipoint(lorb) + korb 
          END IF
          
         
!         if(lk.gt.kl) THEN
!         klkl = ipoint(lk) + kl
!         ELSE
          klkl = ipoint(kl) + kl
!         END IF 
         e = e - e2ints(klkl)
         END IF
         end do
         end do


      
     
       e = e*ep
  
      return
      end
