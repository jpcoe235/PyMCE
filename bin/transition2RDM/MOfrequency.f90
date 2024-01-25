




PROGRAM MOcount
implicit none
integer i,csize,j,cmaxi,nword,k,itemp
double precision, allocatable:: carray(:)
integer, allocatable:: alphas(:,:)
integer, allocatable:: betas(:,:)
integer, allocatable:: MOcountAlpha(:)
integer, allocatable:: MOcountBeta(:)
double precision, allocatable:: WeightMOAlpha(:)
double precision, allocatable:: WeightMOBeta(:)
integer*8 bigtemp1,bigtemp2
double precision cmax,norm,weightcut
integer pmax,bintemp,alphae,subalpha,betae,subbeta
integer int_bits,nbft,orbital,subtotal
integer n,jshift
logical spin1,spin2




int_bits=32 !This is set in params at compile


 
PRINT *, 'number of basis functions'
READ *, nbft


PRINT *, 'Print weighted percentage larger than ?%'
READ *,weightcut

ALLOCATE(MOcountAlpha(nbft))
ALLOCATE(MOcountBeta(nbft))
ALLOCATE(WeightMOBeta(nbft))
ALLOCATE(WeightMOAlpha(nbft))

! nword = nbft/int_bits + 1 int_bits=8*nbyte =8*4=32
nword=nbft/int_bits+1
PRINT *, 'nword'
READ *, nword
OPEN(UNIT=14, FILE='civ_out')

j=0
 csize=-1
do while (j==0)
READ (14,*,IOSTAT=j)
 csize=csize+1
end do
PRINT *,csize
 csize=csize/nword


ALLOCATE(carray(csize))
ALLOCATE(alphas(nword,csize))
ALLOCATE(betas(nword,csize))
 CLOSE(14)

OPEN(UNIT=14, FILE='civ_out')
do i=1,csize     

READ(14,*) itemp,carray(i),bigtemp1,bigtemp2

alphas(1,i)=bigtemp1
betas(1,i)=bigtemp2



do k=2,nword
READ (14,*)  bigtemp1,bigtemp2
alphas(k,i)=bigtemp1
betas(k,i)=bigtemp2
end do

end do
 CLOSE(14)


!Approximate normalization
norm=0.0D0
do i=1,csize
norm=norm+carray(i)**2
end do

 carray=carray/DSQRT(norm)

PRINT *,'norm^2=',norm




! Find MO labels

MOcountAlpha=0
MOcountBeta=0
do i=1,csize

!alpha
do j=0, nbft-1

     n = j/int_bits + 1
     jshift = j - (n-1)*int_bits
     spin1   = btest(alphas(n,i), jshift)
!     spin2   = btest(icij(2,n,ici), jshift)
   
        if(spin1) THEN 


MOcountAlpha(j+1)=MOcountAlpha(j+1)+1
WeightMOAlpha(j+1)=WeightMOAlpha(j+1)+carray(i)**2
END IF
end do

!beta
do j=0, nbft-1

     n = j/int_bits + 1
     jshift = j - (n-1)*int_bits
     spin1   = btest(betas(n,i), jshift)
   
 if(spin1) THEN 



MOcountBeta(j+1)=MOcountBeta(j+1)+1
WeightMOBeta(j+1)=WeightMOBeta(j+1)+carray(i)**2
END IF

end do

end do


PRINT *, 'MOs in percentage of configurations'

PRINT *,'csize',csize
do j=1,nbft

PRINT *,j,100.0D0*MOcountAlpha(j)/csize,100.0D0*MOcountBeta(j)/csize

end do


PRINT *, 'Larger than 10%'
do j=1,nbft

if (100.0D0*MOcountAlpha(j)/csize.gt.10.0D0.OR.100.0D0*MOcountBeta(j)/csize.gt.10.0D0) THEN
PRINT *,j,100.0D0*MOcountAlpha(j)/csize,100.0D0*MOcountBeta(j)/csize
END IF



end do


PRINT *,'Weighted by coefficient squared approximately normalised'
do j=1,nbft

PRINT *,j,100.0D0*WeightMOAlpha(j),100.0D0*WeightMOBeta(j)

end do


PRINT *, 'Larger than',weightcut,'%'
do j=1,nbft

if (100.0D0*WeightMOAlpha(j).gt.weightcut.OR.100.0D0*WeightMOBeta(j).gt.weightcut) THEN
PRINT *,j,100.0D0*WeightMOAlpha(j),100.0D0*WeightMOBeta(j)
END IF



end do



DEALLOCATE(carray)
DEALLOCATE(alphas)
DEALLOCATE(betas)
DEALLOCATE(MOcountAlpha)
DEALLOCATE(MOcountBeta)
DEALLOCATE(WeightMOAlpha)
DEALLOCATE(WeightMOBeta)
END PROGRAM





