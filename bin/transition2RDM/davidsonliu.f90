  subroutine davidsonliu(length,ieig,idiag,ieigmax)
  use commonarrays, only: hf, sf, b, h, ijh, s, ijs, c, e,ctemp
  use dyn_par
  !implicit none
  implicit   real*8  (a-h,o-z)
  parameter           (eps=1.0d-15,maxit=300)
  !common  /reduced/   hf(kmax*(kmax+1)/2), sf(kmax*(kmax+1)/2)
  !common  /bk/        b(maxc,kmax)
  !common  /hands/     h(maxh), ijh(maxh), s(maxs), ijs(maxs)
  !common  /coef/      c(maxc)
  !common  /eigen/     e(kmax)
  integer             length, info, idiag,is
  integer             ieigmax
  integer             i, ici, k, kk, kkk, klim
  real*8              z(kmax,kmax), btemp(maxc,kmax)
  real*8              vnorm, rnorm, vnormA(ieig)
  real*8              a(length,kmax),d(length,kmax),r(length,ieig)
  real*8              work(3*kmax)

  !     write(34,*)
  !     write(34,*) 'entering davidson...' 

        cstop=davidson_stop**2
idiag = 1
 111  continue

           klim = min(kmax,length)
            !need at least ieig steps as we step by ieig but we always have 100
             
            kk=ieig
            do kloops=ieig, klim,ieig !step by ieig as we create ieig new bks each iteration
           
!           generate reduced h and s matrices need to fix lower limit
            call SA_h_s_reduced(length,1,kk,a,d)
!           solve nonorthogonal eigenvalue problem using routine from LAPACK
!           hred * z = e * sred * z
            call dspgv(1 , 'V' ,'U', kk, hf, sf, e, z, kmax,work,info)
            if(info.ne.0) THEN
            PRINT *, 'error in dspgv',info,kmax,idiag,kk,length
            STOP
            ENDIF
            r=0.0D0
!            PRINT *,r(1,ieig), 'r'
            do k=1, kk
!             We have stored a and d during the call to SA_h_s_reduced
!             

              
            do is=1,ieig
            do ici=1, length

              
            r(ici,is) = r(ici,is)+z(k,is)*(a(ici,k)-e(is)*d(ici,k))
            end do

               enddo
            enddo



!           residual norm
             vnormA = 0.0D0  !Array
            do is=1,ieig
            do ici = 1, length
               
            vnormA(is)=vnormA(is)+r(ici,is)*r(ici,is)
!               vnorm1 = vnorm1 + r1(ici)*r1(ici)
            enddo
            enddo

! remove sqrt for rnorm expression and set cstop to davidson_stop**2 initially instead.
           
             rnorm = maxval(vnormA)
       
            
             
             if(rnorm.lt.cstop) goto 222          ! finished?
!           form correction vector and orthogonalize to b_ks
            if(kk.le.(klim-ieig)) then       ! Can we add ieig values to kk without passing kmax array size
               

! other bk ground
                 do is=1,ieig  
                 do ici=1, length
                 if( dabs(r(ici,is)) .lt. eps) then
                     b(ici,kk+1) = 0.0D0
                  else
                     b(ici,kk+1) = r(ici,is)/(e(is)*s(ici)-h(ici))
                  endif




                 enddo


               vnorm = 0.0D0
               do ici = 1, length
                  vnorm = vnorm + b(ici,kk+1)*b(ici,kk+1)
               enddo
               vnorm = dsqrt(vnorm)
               do ici= 1,length
                  b(ici,kk+1) = b(ici,kk+1)/vnorm
               enddo
!              orthogonalize b_k vectors               
               do k=kk,1, -1
                  dot = 0.0D0
                  do ici = 1, length
                     dot = dot + b(ici,k)*b(ici,kk+1)
                  enddo
                  do ici = 1, length
                     b(ici,kk+1) = b(ici,kk+1) - dot*b(ici,k)
                  enddo
	       enddo
               vnorm = 0.0D0
               do ici = 1, length
                  vnorm = vnorm + b(ici,kk+1)*b(ici,kk+1)
               enddo
               vnorm = dsqrt(vnorm)
               do ici= 1,length
                  b(ici,kk+1) = b(ici,kk+1)/vnorm
               enddo
               


                 kk=kk+1    !We have added a new bk
                enddo
                               endif

          enddo
  
  !!!EDITED moved continue from here

!      generate new b_k vectors from eigenvectors of hred
       do k=1, kmax
          do i=1, maxc
             btemp(i,k) = 0.0D0
          enddo
       enddo
       do i= 1,ieigmax
          do ici=1, length
             do k=1, kk           ! last value of kk when jumped above
                btemp(ici,i) = btemp(ici,i) + z(k,i)*b(ici,k)
             enddo
          enddo
       enddo
!      normalize b_k vectors
       do k=1,ieigmax
          vnorm = 0.0D0
          do ici=1, length
             b(ici,k) = btemp(ici,k)
             vnorm  = vnorm + b(ici,k)*b(ici,k)
          enddo
          vnorm = dsqrt(vnorm)
         
          do ici=1, length
             b(ici,k) = b(ici,k)/ vnorm
          enddo
          do ici=length+1, maxc
             b(ici,k) = 0.0D0
          enddo
       enddo
!      orthogonalize the new b_k vectors
       do kkk=2,ieigmax
          do k=kkk-1, 1, -1
             dot = 0.0D0
             do ici = 1, length
                dot = dot + b(ici,k)*b(ici,kkk)
             enddo
             do ici = 1, length
                b(ici,kkk) = b(ici,kkk) - dot*b(ici,k)
             enddo
          enddo
          vnorm = 0.0D0
          do ici=1, length
             vnorm = vnorm + b(ici,kkk)*b(ici,kkk)
          enddo
          vnorm= dsqrt(vnorm)
          do ici = 1, length
             b(ici,kkk) = b(ici,kkk)/vnorm
          enddo
       enddo

!      if not converged, jump back and try again or STOP
       if(rnorm.gt.cstop) then
         idiag = idiag + 1 
         if(idiag.gt.maxit) STOP 'TOO MANY DAVIDSON ITERATIONS'
         goto 111
       endif
  
 222   continue  

!!!!!!! From v4 mcci davidson - with moved 222 continue.  We do not orthonormalise on final step
! so the higher eigenvalues agree with the expectation values as the eigenfunctions have fewer numerical issues.
! might lose orthogonality of states


!     generate new b_k vectors from eigenvectors of hred
      do i=1, maxc
         do k=1, kmax
            btemp(i,k) = 0.0
         enddo
      enddo
      do i= 1, ieig
         do ici=1, length
            do k=1, kk        ! we've just expanded in this many b_k's
               btemp(ici,i) = btemp(ici,i) + z(k,i)*b(ici,k)
            enddo
         enddo
      enddo
!     We guess these eigenvectors are linearly independent: true if they
!     come from different eigenvalues, and likely to be true otherwise
!     if the LAPACK routine is doing its job.

      do k=1, ieig
         do ici=1, maxc         ! better than length, for branching
            b(ici,k) = btemp(ici,k)
         enddo
      enddo
!!!!!!!!!!
  
!      update c vectors            
       
           do ici=1,length  
           c(ici) = 0.0D0
           end do


            do is=1,ieig-1
            do ici=1,length
            ctemp(ici,is)=0.0D0
            end do
            end do
 

       do ici= 1, length
!          do k=1, kk
!             c(ici) = c(ici) + z(k,ieig)*b(ici,k)
              c(ici) =  b(ici,ieig)
             do is=1,ieig-1
!             ctemp(ici,is) = ctemp(ici,is) + z(k,is)*b(ici,k)
             ctemp(ici,is) = b(ici,is)
             end do


 !         enddo
       enddo   

      return
      end

