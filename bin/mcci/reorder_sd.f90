! rewrote to just work with SDs. JPC

      subroutine reorder_sd(ici,jci,ndiff,idiff1,idiff2,ep)
     use commonarrays, only: icij, nword, nbft, ntotal, n_alpha, n_beta, list, my_pair
     use dyn_par
     implicit none
      integer           idoubles(iword), jdoubles(iword)
      logical           spin1, spin2, swapped1, swapped2
     integer l,k,itemp,j,jshift,n,i
     integer idiff1,idiff2,ndiff,ici,jci
     double precision ep
     

!     convert from "dictionary ordering" to "maximal coincidence ordering" 


      ep=1.0D0
!     initialize list  to 0
      do n=1,ntotal
         list(1,n) = 0
         list(2,n) = 0

      enddo

! convert string to orbital labels

      k = 0
      l = 0

      
      do j=0, nbft-1

 	 n = j/int_bits + 1
	 jshift = j - (n-1)*int_bits
     
	 spin1   = btest(icij(1,n,ici), jshift)
	 spin2   = btest(icij(2,n,ici), jshift)
	 
	    
	    if(spin1) THEN
            k=k+1
            list(1,k) = j+1
            END IF
	    if(spin2) THEN
            k=k+1
            list(1,k) = j+1 + nbft
            END IF

	 spin1   = btest(icij(1,n,jci), jshift)
	 spin2   = btest(icij(2,n,jci), jshift)
	
	    
	    if(spin1) THEN
            l = l+1
            list(2,l) = j+1
            END IF
	    
            if(spin2) THEN
            l = l+1
            list(2,l) = j+1 + nbft
            END IF

            enddo

!     these should simply be the number of occupied orbitals
      if(k.ne.ntotal .or. l.ne.ntotal) then
!          write(*,*)'nttl nalpha nbeta',ntotal,n_alpha,n_beta
!          write(*,*)'kdoubly ldoubly',kdoubly,ldoubly
!          write(*,*)'k l',k,l
!         call ldump(1,2,0)
         STOP 'error reorder: wrong particle number'
      endif
! find out where they share and swap - if not equal then look through remaining to make equal
! no duplicates so can consider full list each time

!      ndiff=0
      do i=1,ntotal

      if (list(1,i).NE.list(2,i)) THEN
!      ndiff=ndiff+1
      do j=1,ntotal
      if (list(1,i).EQ.list(2,j)) THEN
!      ndiff=ndiff-1
!     SWAP
      itemp=list(2,j)
      list(2,j)=list(2,i)
      list(2,i)=itemp
      ep=-ep
     EXIT
      END IF
      
     
       
      end do      

      END IF
      end do
! locate differences
       ndiff=0
       do i=1,ntotal
       if((list(1,i).NE.list(2,i))) THEN
       ndiff=ndiff+1
        if(ndiff.EQ.1) THEN
        idiff1=i
        END IF
         if(ndiff.EQ.2) THEN
        idiff2=i
        END IF

        END IF
       end do

     

   
  return 
end subroutine reorder_sd
