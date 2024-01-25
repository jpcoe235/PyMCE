!determinants differ by two orbitals
      subroutine slater2(e,i_am_mu,i_am_nu,nu_doubly,idiff1,idiff2,kck,n_2p,ep)
      use commonarrays, only: nbft, ntotal, i_sx2, e1ints, e2ints, ipoint, list
      use dyn_par
      implicit real*8    (a-h,o-z)
      integer spini,spinj,spink,spinl
 
       
! different orbitals using notation of advances in highly correlated approaches
        iorb = list(1,idiff1)
        jorb = list(1,idiff2)
        korb = list(2,idiff1)
        lorb = list(2,idiff2)

!     zero energy
      e = 0.0d0

        
        spini=1
        spinj=1
        spink=1
        spinl=1

        if(iorb.gt.nbft) THEN
         iorb = iorb - nbft
         spini=-1
         END IF
         
          if(jorb.gt.nbft) THEN
         jorb = jorb - nbft
         spinj=-1
         END IF
         
          if(korb.gt.nbft) THEN
         korb = korb - nbft
         spink=-1
         END IF 
          
           if(lorb.gt.nbft) THEN
         lorb = lorb - nbft
         spinl=-1
         END IF


! check  for Coulomb contribution        
         if(spini.EQ.spink) THEN 
         if(spinj.EQ.spinl) THEN
! access correct coloumb integral
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

         

         e = e + e2ints(ikjl) 
         END IF
         END IF

! check if i and k and k and j are the same for Exchange      
        if(spini.EQ.spinl)  THEN
        if (spinj.EQ.spink) THEN
         
        if(lorb.gt.iorb) THEN
        il=ipoint(lorb)+iorb
        ELSE
        il=ipoint(iorb)+lorb
        END IF
        
        if(korb.gt.jorb) THEN
        jk=ipoint(korb)+jorb
        ELSE
        jk=ipoint(jorb)+korb
        END IF

        if(jk.gt.il) THEN
        iljk=ipoint(jk)+il
        ELSE
        iljk=ipoint(il)+jk
        END IF
        e = e - e2ints(iljk) 
        END IF
        END IF
         
       
         



     
      e = e*ep
  
      return
      end
