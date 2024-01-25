subroutine spin_sparse_sd(length,llast)
  use commonarrays!, only: h, ijh, s, ijs, cnorm, i_sx2
  use dyn_par
  implicit   real*8    (a-h,o-z)



  !     diagonals
  do ici = llast+1, length
     call muSnu_sd(ici,ici,exS2)
!     cnorm(ici) = dsqrt( dabs(ck(i_sx2,n_2p,0)) )
     h(ici) = exS2
!     s(ici) = ck(i_sx2,n_2p,0)
  enddo
  ijh(1) = length + 2
  ijs(1) = length + 2

  if(llast.gt.0) then 
     !        number of new configurations
     lnew = length - llast
     !        set index to last location to resume generating h
     k = ijh(llast+1) - 1
     l = ijs(llast+1) - 1
  else
     k = length + 1
     l = length + 1
  endif

  do ici = llast+1, length
     do jci = 1, ici-1

        call muSnu_sd(ici,jci,exS2)
        
     
       
           if( dabs(exS2).ge.hmin ) then
           k = k+1
           if(k.gt.maxh) STOP 'gen_h_s: exceeded matrix storage h'
           h(k) = exS2
           ijh(k) = jci
           end if


     enddo

     ijh(ici+1) = k+1

  enddo
  return
end subroutine spin_sparse_sd
