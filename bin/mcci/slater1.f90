!determinants differ by one orbital

     subroutine slater1(e,i_am_mu,i_am_nu,nu_doubly,idiff1,kck,n_2p,ep)
      use commonarrays, only: nbft, ntotal, i_sx2, e1ints, e2ints, ipoint, list
      use dyn_par
      implicit real*8    (a-h,o-z)
      integer spini,spinj,spink

       

        iorb = list(1,idiff1)
        jorb = list(2,idiff1)

!     zero energy
        e = 0.0d0
!     one particle contribution if spins are the same

         spini=1
         spinj=1
        
         if(iorb.gt.nbft) THEN
         iorb = iorb - nbft
          spini=-1
          END IF
        
  
          if(jorb.gt.nbft) THEN
          jorb=jorb-nbft
          spinj=-1  
          END IF
  
	 if(spini.EQ.spinj) THEN 
         if(iorb.gt.jorb) THEN
         kk = ipoint(iorb) + jorb
         ELSE
         kk = ipoint(jorb) + iorb
         END IF
         e = e + e1ints(kk)
         END IF




!   two particle contributions 
       do k=1,ntotal
       
        if (k.EQ.idiff1) cycle
        korb = list(1,k)

        
        spink=1
        if(korb.gt.nbft) THEN
         korb = korb - nbft
         spink=-1
         END IF
! check i j spins are the same for Coulomb contribution        
         if(spini.EQ.spinj) THEN 
          
          if (jorb.gt.iorb) THEN
          ij = ipoint(jorb) + iorb
          ELSE
          ij = ipoint(iorb) + jorb 
          END IF 
         kk = ipoint(korb) + korb
         
         if(kk.gt.ij) THEN
          ijkk = ipoint(kk) + ij 
         ELSE
         ijkk = ipoint(ij) + kk 
         END IF
         e = e + e2ints(ijkk) 
         END IF


! check if i and k and k and j are the same for Exchange      
        if(spini.EQ.spink)  THEN
        if (spinj.EQ.spink) THEN

        if(korb.gt.iorb) THEN
        ik=ipoint(korb)+iorb
        ELSE
        ik=ipoint(iorb)+korb
        END IF

         if(jorb.gt.korb) THEN
        kj=ipoint(jorb)+korb
        ELSE
        kj=ipoint(korb)+jorb
        END IF

        if(kj.gt.ik) THEN
         ikkj=ipoint(kj)+ik
        ELSE
         ikkj=ipoint(ik)+kj
        END IF
        

        e = e - e2ints(ikkj) 
        END IF
        END IF
         
        end do
         



     
      e = e*ep
  
      return
      end
