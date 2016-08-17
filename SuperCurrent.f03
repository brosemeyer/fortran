!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine current_y(t,b,J,Ky,N,Ny,i_t,i_b,rr, &
				             Nc,Neigs,NNc,Ep,Fp,Uxp,Vxp)
implicit none
integer :: i_t, i_b, NNc, N, Ny, i
integer, dimension(Ny) :: Nc, Neigs
real(8), dimension(Ny) :: Ky
real(8) :: t, b, D0
REAL(8), dimension(Ny,NNc,2) :: Ep, Fp
REAL(8), dimension(Ny,N) :: Jsum
REAL(8), dimension(N) :: J, rr
COMPLEX(8), DIMENSION(NNc,Ny,N) :: Uxp, Vxp
CHARACTER(len=2048) :: filename

Jsum(:,:) = 0.0
!$OMP PARALLEL DO SCHEDULE(DYNAMIC)
do i = 1, Ny
	if (Neigs(i)>0) then
	Jsum(i,:) = sum(abs(Uxp(1:Neigs(i),i,:))**2 &
	      			  - abs(Vxp(1:Neigs(i),i,:))**2,Dim=1)*Ky(i)
	end if
end do
!$OMP END PARALLEL DO

J = sum(Jsum, Dim=1)

write(filename,'(A,I2,I2,A)') "./SELF_CONSISTENT/current",i_t,i_b, ".dat"
call save_file_2(rr,N,1,filename,J)

end subroutine current_y
