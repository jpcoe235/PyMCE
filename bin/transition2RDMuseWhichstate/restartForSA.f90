subroutine restartForSA(length,ieig)
  use commonarrays, only: ctemp, icij, nword
  use dyn_par
  integer           length, idummy !icij !declaration moved to commonarrays.f90 
  integer*8  bigtemp1,bigtemp2
  integer ieig,is
  CHARACTER (LEN=30)  :: f1 

PRINT *,'ieig in restartForSA',ieig

do is=1,ieig

write(f1,'(i2)') is

            f1='state'//trim(adjustl(f1))
      
 
            f1='civ_out_'//trim(adjustl(f1))

print *,f1
  open(unit=60,file=f1)
  i = 0
11 continue
  i = i + 1	
  do n=1, nword
     !        read(60,*,end=22) idummy,icij(1,n,i),icij(2,n,i)
     if(n.eq.1) then
        !           read(60,*,end=22) c(i),icij(1,n,i),icij(2,n,i)
        read(60,*,end=22) idummy,ctemp(i,is),bigtemp1,bigtemp2
        icij(1,n,i)=bigtemp1
        icij(2,n,i)=bigtemp2 
    else
        read(60,*,end=22) bigtemp1,bigtemp2
         icij(1,n,i)=bigtemp1
        icij(2,n,i)=bigtemp2 
     endif
  enddo
  goto 11
22 continue
  length = i-1

PRINT *,'state',is,length

end do ! end of is loop over states


  return
end subroutine restartForSA
