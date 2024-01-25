subroutine muSnu_sd(ici,jci,exS2)
  use commonarrays, only: nbft, icij, nword, i_sx2,list,n_beta,n_alpha,ntotal
  use dyn_par
  implicit none
  integer           ijsingles(iword), ijdoubles(iword),ij_sd(iword)
  integer iSame,sumdoubles,n,ndiff,ispin,ispin1,ispin2
  integer idoubles,iorb1,iorb2,idiff1,idiff2,i,ici,jci
  double precision exS2,ep
  !common  /aodat/   nsym, nbft, nbpsy(irmax)
  !common  /config/  icij(2,iword,maxc), nword
  !common  /occupy/  ntotal, n_alpha, n_beta, i_sx2

  exS2 = 0.0d0
  
!        PRINT *,'*****'
!        PRINT *,ici,jci

        call reorder_sd(ici,jci,ndiff,idiff1,idiff2,ep)

!         PRINT *,(list(1,n), n=1,ntotal)
!         PRINT *,(list(2,n), n=1,ntotal)
!         PRINT *,ep,ndiff
!         READ *,ispin

     if(ndiff.eq.0) THEN 

!   iSame=1
!      do n=1,nword

!      do ispin=1,2
!      if(icij(ispin,n,ici).ne.icij(ispin,n,jci)) THEN

!      iSame=0
!      end if
!      end do
!      end do
 
!      do n=1,nword
!      do ispin=1,2
!      PRINT *,iSame,icij(ispin,n,ici),icij(ispin,n,jci)
!      end do
!      end do


! If they are the same state then work out how many doubles
!       if(iSame.eq.1) THEN
       sumdoubles=0
       do  n=1,nword
       idoubles = iand(icij(1,n,ici),icij(2,n,ici))

       do i=0,int_bits-1
       if(BTEST(idoubles,i)) sumdoubles=sumdoubles+1
       end do
        

       end do
!     PRINT *,sumdoubles
!       PRINT *,idoubles
       exS2=0.5D0*(n_alpha-n_beta)
       exS2=(exS2*(exS2+1.0D0))+1.0D0*n_beta-1.0D0*sumdoubles !n_beta from i=j -doubles from i.ne.j terms as we then need to swap orbital and spin to give max coincidence.
       exS2=exS2*ep
!       PRINT *,'exS2',exS2
!       READ *,ispin
       RETURN
       END IF

   


!         PRINT *,'ndiff=',ndiff,idiff1,idiff2
!         PRINT *,(list(1,i),i=1,ntotal)
!         PRINT *,(list(2,i),i=1,ntotal)
      
            if(ndiff.ne.2) THEN
            RETURN
! if two differences can we swap two spins to remove difference?
            ELSE
            
            !check that they are alpha and beta
            ispin1=1
            if(list(1,idiff1).gt.nbft) ispin1=2
            ispin2=1
            if(list(1,idiff2).gt.nbft) ispin2=2
            
!            PRINT *,ispin1,ispin2
            if(ispin1.eq.ispin2) THEN
!            PRINT *, 'same spins'
            RETURN
            end if

! ok the spins are different - can we swap
! turn beta to alpha and vice versa from two diff in list 1
! and check if this is now the same as the two differences in list 2

            if(list(1,idiff1).gt.nbft) THEN
            iorb1=list(1,idiff1)-nbft
            ELSE
            iorb1=list(1,idiff1)+nbft
            endif
            
          if(list(1,idiff2).gt.nbft) THEN
            iorb2=list(1,idiff2)-nbft
            ELSE
            iorb2=list(1,idiff2)+nbft
            endif
!            PRINT *,iorb1,iorb2


! do we have a match now without changing order of list2

            if(list(2,idiff1).eq.iorb1) THEN
              if(list(2,idiff2).eq.iorb2) THEN
!              PRINT *,'spin swap gives match no change to ep'
                exS2=1.0D0*ep
                 
            
            end if
            end if 

                 if(list(2,idiff1).eq.iorb2) THEN
              if(list(2,idiff2).eq.iorb1) THEN
!              PRINT *,'spin swap then orbital and spin swap'
              
               exS2=-1.0D0*ep
          

                   
             end if
            end if     


             end if



        
          Return
      end
