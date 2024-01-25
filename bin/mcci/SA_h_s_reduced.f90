subroutine SA_h_s_reduced(length,kl,ku,a,d)
  use commonarrays, only: h, ijh, s, ijs, hf, sf, b
  use dyn_par
  implicit   real*8  (a-h,o-z)
  !common  /hands/     h(maxh), ijh(maxh), s(maxs), ijs(maxs)
  !common  /reduced/   hf(kmax*(kmax+1)/2), sf(kmax*(kmax+1)/2) 
  !common  /bk/        b(maxc,kmax)
  integer             i, j, k, kl, ku, length
  real*8              d(length,kmax),a(length,kmax)

!     matrices are UPPER packed
      if(kl.le.0) kl=l 
      do j = kl, ku
! moved from inner loop
!           call mxv_sparse(length,h,ijh,b(1,j),ai)
!  mxv_sparse moved here and a and d are stored


        do ici=1, length
!        diagonal
        a(ici,j) = h(ici)*b(ici,j)
      enddo
                   
      do ici=1, length
         do ki=ijh(ici), ijh(ici+1)-1
            jci =ijh(ki)
!           lower off diagonals
	    a(ici,j) = a(ici,j) + h(ki)*b(jci,j)
!           upper off diagonals
 	    a(jci,j) = a(jci,j) + h(ki)*b(ici,j)
         enddo
      enddo

      do ici=1, length
!        diagonal
         d(ici,j) = s(ici)*b(ici,j)
      enddo
!            mv=m*v          
      do ici=1, length
         do ki=ijs(ici), ijs(ici+1)-1
            jci =ijs(ki)
!           lower off diagonals
	    d(ici,j) = d(ici,j) + s(ki)*b(jci,j)
!           upper off diagonals
 	    d(jci,j) = d(jci,j) + s(ki)*b(ici,j)
         enddo
      enddo
      






         do i=1, j

            

            hf(i+(j-1)*j/2) = 0.0d0
            sf(i+(j-1)*j/2) = 0.0d0
            do k = 1, length   
               hf(i+(j-1)*j/2) = hf(i+(j-1)*j/2) + b(k,i)*a(k,j)
               sf(i+(j-1)*j/2) = sf(i+(j-1)*j/2) + b(k,i)*d(k,j)
            enddo

         enddo
      enddo

      return
end subroutine SA_h_s_reduced
