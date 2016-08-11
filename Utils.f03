!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Fourier_q(Fout,FinS,N,Ns,Nci,dir)
use, intrinsic :: iso_c_binding
       include 'fftw3.f03'
!~ implicit none
integer :: j, N, dir, Ns
INTEGER, DIMENSION(Ns) :: Nci
COMPLEX(8), DIMENSION(N) :: Fout, Fin
COMPLEX(8), DIMENSION(Ns) :: FinS
type(C_PTR) :: plan

if (Ns < N) then
	Fin = cmplx(0.0,0.0)
	do j=1,Ns
		Fin(Nci(j)) = FinS(j)
	end do
elseif (Ns == N) then
	Fin = FinS
else
	write(*,*) 'BAD LENGTHS IN FFT'
end if

	if (dir == 1) then
		plan = fftw_plan_dft_1d(N, Fin, Fout, FFTW_FORWARD, FFTW_ESTIMATE)
		call fftw_execute_dft(plan,Fin, Fout)
		call fftw_destroy_plan(plan)
	else if (dir == 2) then
		plan = fftw_plan_dft_1d(N, Fin, Fout, FFTW_BACKWARD, FFTW_ESTIMATE)
		call fftw_execute_dft(plan,Fin, Fout)
		call fftw_destroy_plan(plan)
	else
		write(*,*) 'fftw direction not supported'
	end if

end subroutine Fourier_q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE save_file(X,N,filename, Y1)
implicit none
integer :: N, i, j, k, l
REAL(8), dimension(N) :: X, Y1
character(*) :: filename
101 format (E12.6)
102 format (A)
110 format (E12.6,A)

open(unit = 1, file = trim(filename))
do i=1,N
	write(1,110,advance="no") X(i), ', '
	write(1,101,advance="yes") Y1(i)
end do
close(1)
end SUBROUTINE save_file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE save_file_2(X,Nx,Ny,filename, Y1)
implicit none
integer :: Nx, Ny, i, j, k, l
REAL(8), dimension(Nx,Ny) :: Y1
REAL(8), dimension(Nx) :: X
character(*) :: filename
101 format (E12.6)
102 format (A)
110 format (E12.6,A)

open(unit = 1, file = trim(filename))
do i=1,Nx
	write(1,110,advance="no") X(i), ', '
	do j=1,Ny
		write(1,110,advance="no") Y1(i,j), ', '
	end do
	write(1,102,advance="yes") ''
end do
close(1)
end SUBROUTINE save_file_2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE save_file_3(X,Nx,Y,Ny,Nz,filename, Z)
implicit none
integer :: Nx, Ny, i, j, k, Nz
REAL(8), dimension(Nx,Ny,Nz) :: Z
REAL(8), dimension(Nx) :: X
REAL(8), dimension(Ny) :: Y
character(*) :: filename
101 format (E12.6)
102 format (A)
110 format (E12.6,A)
111 format (A, E12.6)

open(unit = 1, file = trim(filename))
do j=1,Ny
	do i=1,Nx
		write(1,101,advance="no") X(i)
		write(1,111,advance="no") ', ', Y(j)
		do k=1,Nz
			write(1,111,advance="no") ', ', Z(i,j,k)
		end do
		write(1,102,advance="yes") ''
		if (i == Nx .and. Nx /= 1) then
			write(1,102,advance="yes") ''
		end if
	end do
end do
close(1)
end SUBROUTINE save_file_3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE save_file_3c(X,Nx,Y,Ny,Nz,filename, Z)
implicit none
integer :: Nx, Ny, i, j, k, Nz
COMPLEX(8), dimension(Nx,Ny,Nz) :: Z
REAL(8), dimension(Nx) :: X
REAL(8), dimension(Ny) :: Y
character(*) :: filename
101 format (E12.6)
102 format (A)
110 format (E12.6,A)
111 format (A, E12.6)

open(unit = 1, file = trim(filename))
do j=1,Ny
	do i=1,Nx
		write(1,101,advance="no") X(i)
		write(1,111,advance="no") ', ', Y(j)
		do k=1,Nz
			write(1,111,advance="no") ', ', real(Z(i,j,k))
			write(1,111,advance="no") ', ', aimag(Z(i,j,k))
		end do
		write(1,102,advance="yes") ''
		if (i == Nx .and. Nx /= 1) then
			write(1,102,advance="yes") ''
		end if
	end do
end do
close(1)
end SUBROUTINE save_file_3c
