!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine Free_Energy(D0,Delta,t,b,EE,Ux,Vx,NN,Nc,NNc,Nci,Ny,kx,ky,Consistent0,Free_E,rr,it,ib)
implicit none
INTEGER :: NN, Ny, NNc, jj, it, ib
INTEGER, DIMENSION(Ny) :: Nc
INTEGER, DIMENSION(Ny,NNc) :: Nci
DOUBLE PRECISION, parameter :: PI=3.141592
DOUBLE PRECISION, DIMENSION(Ny,NNc,2) :: EE, E0
COMPLEX(8), DIMENSION(NNc,Ny,NN) :: Ux, Vx
DOUBLE PRECISION, DIMENSION(Ny,NNc) :: E_presort
DOUBLE PRECISION :: t, b, D0
DOUBLE PRECISION, DIMENSION(NN) :: Delta, kx, Free_x, rr
DOUBLE PRECISION, DIMENSION(Ny) :: ky, F2s
DOUBLE PRECISION, DIMENSION(Ny,NN) :: fs
DOUBLE PRECISION, DIMENSION(3) :: Free_E
DOUBLE PRECISION ::  Consistent0, term
LOGICAL, DIMENSION(Ny,NNc) :: msk
INTEGER :: i, j
CHARACTER(len=2048) :: filename

msk = .TRUE.
!$OMP PARALLEL DO PRIVATE(jj)
do i=1,Ny
	do j=1,Nc(i)
		E_presort(i,j) = kx(Nci(i,j))**2+ky(i)**2-1
	end do
	do j=1,Nc(i)
		jj = minloc(abs(E_presort(i,1:Nc(i))),MASK=msk(i,1:Nc(i)),DIM=1)
		msk(i,jj) = .FALSE.
		E0(i,j,1) = E_presort(i,jj) - b*D0
		E0(i,j,2) = E_presort(i,jj) + b*D0
	end do
	do j=1,NN
	fs(i,j) =-sum((abs(Ux(1:Nc(i),i,j))**2+abs(Vx(1:Nc(i),i,j))**2) &
				*log(cosh(EE(i,1:Nc(i),1)/(2*t*D0))/cosh(E0(i,1:Nc(i),1)/(2*t*D0))) &
				+ (abs(Ux(1:Nc(i),i,j))**2+abs(Vx(1:Nc(i),i,j))**2) &
				*log(cosh(EE(i,1:Nc(i),2)/(2*t*D0))/cosh(E0(i,1:Nc(i),2)/(2*t*D0))))
	end do
	F2s(i) = sum(abs(E0(i,1:Nc(i),:)) - EE(i,1:Nc(i),:))
end do
!$OMP END PARALLEL DO
fs = fs*t/(2*Consistent0)
Free_E(1) = sum(fs)/dble(NN)
Free_E(2) = sum(Delta**2)/dble(2*NN)
Free_E(3) = sum(F2s)/(dble(NN)*Consistent0)
Free_x = sum(fs,dim=1)/dble(NN) + Delta**2/dble(NN)

write(filename,'(A,I2,I2,A)') "./SELF_CONSISTENT/Free_E",it,ib, ".dat"
call save_file_2(rr,NN,1,filename,Free_x)
end subroutine Free_Energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DOS(Ux, Vx, EE, Ky, Kx, rr, D0, M, NN,NNc,Nci,Ny,NE,E_cut,eta,it,ib)

implicit none
REAL(8) :: eta
INTEGER :: NN, Ny, NNc, NE
INTEGER ::  i, j, k, ii,it,ib
INTEGER, DIMENSION(Ny,NNc) :: Nci
INTEGER, dimension(Ny) :: M
REAL(8), dimension(Ny) :: Ky
REAL(8), DIMENSION(NE) :: E_DOS
REAL(8), DIMENSION(NN) :: rr, Kx
REAL(8) :: D0, E_cut
CHARACTER(len=2048) :: filename
REAL(8), DIMENSION(Ny,NNc,2) :: EE
COMPLEX(8), dimension(NNc,Ny,NN) :: Ux, Vx
REAL(8), dimension(NE,NN,2) :: G_x

do i = 1, NE
	E_DOS(i) = D0*E_cut*(2*dble(i-1)/dble(NE-1) - 1)
end do
G_x(:,:,:) = cmplx(0.0,0.0)
! TODO make the delta function like relaxation rate
!$OMP PARALLEL DO SCHEDULE(DYNAMIC)
do i = 1, NE
	do k=1,2
		do j = 1, NN
			do ii=1,Ny
				G_x(i,j,k) = G_x(i,j,k)+3.141592*sum(abs(Ux(1:M(ii),ii,j))**2*exp(-(E_DOS(i)-EE(ii,1:M(ii),k))**2/eta**2)  &
										           + abs(Vx(1:M(ii),ii,j))**2*exp(-(E_DOS(i)+EE(ii,1:M(ii),3-k))**2/eta**2))/dble(Ny)
			end do
		end do
	end do
end do
!$OMP END PARALLEL DO

E_DOS = E_DOS/D0
write(filename,'(A,I2,I2,I2,A)') "./DOS/LDOSx.dat"
call save_file_2(rr,NN,2,filename, sum(G_x,dim=1)/dble(NE))
write(filename,'(A,I2,I2,A)') "./DOS/LDOS", it, ib, ".dat"
call save_file_3(E_dos,NE,rr,NN,2,filename, G_x)
end subroutine DOS
