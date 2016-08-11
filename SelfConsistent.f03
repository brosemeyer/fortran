!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine self_consistent(Fr,Kx,Ky,rr,NN,Ny,t,b,D0,flags,phase, &
						   VECS, EIGS, sc_const00,it, ib, &
						   success,w_len,wr_len,wi_len,Nc,Nci,NNc,tol,max_it,E_cut)
implicit none
REAL(8), parameter :: PI=3.141592654
integer :: ii, i, j, k, jj, Ny,j_it,NNc, max_it, badness, it, ib
integer :: NN, st, last_it, success
INTEGER, DIMENSION(Ny) :: w_len, wr_len, wi_len, Nc, Neigs
INTEGER, DIMENSION(Ny,NNc) :: Nci
integer, dimension(9) :: flags
REAL(8) :: t, b, D0, sc_const00, phase, tol, E_cut
REAL(8), DIMENSION(Ny) :: Ky
REAL(8), DIMENSION(NN) :: Kx, Q, rr, Fr, Fr2, Fr_b
REAL(8), DIMENSION(Ny,2*NNc) :: EIGS
COMPLEX(8), DIMENSION(NN) :: Fk2, Fk, Frc, Pk
COMPLEX(8), DIMENSION(2*NNc,Ny,2*NNc) :: VECS
CHARACTER(len=2048) :: filename
Pk(:) = cmplx(0.0,0.0)
Fr_b(:) = cmplx(1.0,0.0)
Neigs = Nc
if (flags(5)<0) then
 	Neigs = Nc-2
 	Fr_b(:) = 0.0
	! Fr_b(1:2) = 1e3
 	call Fourier_q(Pk,cmplx(Fr_b,0.0,kind=8),NN,NN,Nci(1,:),1)
 	Pk = Pk/dble(NN)
	Fr_b(:) = 1.0
	Fr_b(1:10) = 0.0
end if
call Fourier_q(Fk,cmplx(Fr*Fr_b,0.0,KIND=8),NN,NN,Nci(1,:),1)
Fk = Fk/dble(NN)
last_it = 0
st = 0
j_it = 0
success = 0
do while(st == 0)
	write(filename,'(A,I2,I2,A)') "./SELF_CONSISTENT/sc_delta.dat"
	call save_file_2(rr,NN,1,filename,Fr)
	j_it = j_it+1
	badness = 0
	call get_new_profile(Fk2(1:NN/2+1),VECS,t,b,EIGS,Fk,Neigs,Pk, &
						 Kx,Ky,D0,NN,Ny,flags,phase,w_len,wr_len,wi_len,Nc,Nci,NNc,badness,E_cut)
	Fk2(NN/2+2:NN) = conjg(Fk2(NN/2:2:-1))
	Fk2 = Fk2/sc_const00
	call Fourier_q(Frc,Fk2,NN,NN,Nci(1,:),2)
	Fr2 = real(Frc)*Fr_b
	call Fourier_q(Fk2,cmplx(Fr2,0.0,kind=8),NN,NN,Nci(1,:),1)
	Fk2 = Fk2/dble(NN)

	if (last_it == 0) then
		if (sum(abs(Fr - Fr2))/(NN) < tol .and. badness == 0) then  !!!!CONVERGENCE
			success = 1
		else if (j_it > max_it) then
			write(*,*) 'ABORT SELF CONSISTENCY:', t, b, 'MAX ITERATIONS'
			success = -1
		else if (maxval(abs(real(Frc)),dim=1) > 1.1) then
			write(*,*) 'ABORT SELF CONSISTENCY:', t, b, 'TOO BIG'
			success = -1
		else if (maxval(abs(real(Frc)),dim=1) < 0.05) then
			write(*,*) 'ABORT SELF CONSISTENCY:', t, b, 'TOO SMALL'
			success = -1
		end if !! END CONVERGENCE CHECK
		if (success == 1) then
			st = 1
		else if (success == -1) then
			last_it = 1
			Fk = cmplx(0.0,0.0)
		else
			Fk = Fk2
		end if
	else
		st = 1
	end if
	Fr = Fr2
end do

write(filename,'(A,I2,I2,A)') "./SELF_CONSISTENT/sc_delta",it,ib, ".dat"
call save_file_2(rr,NN,1,filename,Fr)

end subroutine self_consistent
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine self_consistent_00(Kx,Ky,NN,Ny,D0, flags,phase,sc_const00, &
							  w_len,wr_len,wi_len,Nc,Nci,NNc,E_cut)

implicit none
REAL(8), parameter :: PI=3.141592654
integer :: Ny, NN, i, NNc, badness
integer, dimension(9) :: flags
REAL(8) :: D0, sc_const00, phase, E_cut
REAL(8), DIMENSION(Ny) :: Ky
COMPLEX(8), DIMENSION(NN) :: Fk, Pk
COMPLEX(8), DIMENSION(NN/2+1) :: Fk2
REAL(8), DIMENSION(NN) :: Kx
COMPLEX(8), DIMENSION(2*NNc,Ny,2*NNc) :: VECS
REAL(8), DIMENSION(Ny,2*NNc) :: EIGS
INTEGER, DIMENSION(Ny,NNc) :: Nci
INTEGER, DIMENSION(Ny) :: Nc
!!!!!!!!!MUST BE TYPE 8 FOR LAPACK ROUTINE!!!!!!!!!!
INTEGER,DIMENSION(Ny) :: w_len, wr_len, wi_len
INTEGER :: INFO, M, ISUPPZ
COMPLEX(8), DIMENSION(2*NNc,2*NNc) :: MAT
INTEGER :: IWORK
REAL(8) :: RWORK
COMPLEX(8) :: WORK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,Ny
	call zheevr('V', 'I', 'U', 2*Nc(i), MAT(1:2*Nc(i),1:2*Nc(i)), 2*Nc(i), dble(0), dble(0), Nc(i)+1, 2*Nc(i), dble(-1), &
				M, MAT(1:2*Nc(i),1), MAT(1:2*Nc(i),1:Nc(i)), 2*Nc(i), ISUPPZ, WORK, -1, RWORK, -1, IWORK, -1, INFO)
	w_len(i) = int(dble(WORK))
	wr_len(i) = int(RWORK)
	wi_len(i) = IWORK
end do
Fk(:) = cmplx(0.0,0.0)
Pk = Fk
Fk(1) = cmplx(1.0,0.0)
call get_new_profile(Fk2,VECS,dble(0),dble(0),EIGS,Fk,Nc,Pk, &
					 Kx,Ky,D0,NN,Ny,flags,phase,w_len,wr_len,wi_len,Nc,Nci,NNc,badness,E_cut)

sc_const00 = abs(Fk2(1))

end subroutine self_consistent_00
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_new_profile(Fk2,VECS,t,b,EIGS_pass,Fk,Neigs,Pk, &
						   Kx,Ky,D0,NN,Ny,flags,phase,w_len,wr_len,wi_len,Nc,Nci,NNc,badness,E_cut)
IMPLICIT NONE
REAL(8), PARAMETER :: PI=3.141592654
INTEGER :: Ny, NN, i, j, jj, NNc, jjm, j_x, badness
INTEGER, DIMENSION(Ny) :: w_len, wr_len, wi_len, Nc, Neigs
INTEGER, DIMENSION(Ny,NNc) :: Nci
REAL(8) :: t, b, D0, phase, E_cut
REAL(8), DIMENSION(NN) :: Kx
REAL(8), DIMENSION(Ny) :: Ky
COMPLEX(8), DIMENSION(NN/2+1) :: Fk2
COMPLEX(8), DIMENSION(NN) :: Fk, Pk
COMPLEX(8), DIMENSION(NN/2+1,Ny) :: Fk_sum
REAL(8), DIMENSION(Ny,NNc) :: F_sc
integer, dimension(9) :: flags
INTEGER, dimension(Ny) :: M
REAL(8), DIMENSION(Ny,NNc) :: EIGS_pass
REAL(8), DIMENSION(Ny,NNc,NNc) :: theta
!!!!!!!!!MUST BE TYPE 8 FOR LAPACK ROUTINE!!!!!!!!!!
INTEGER, DIMENSION(Ny) :: INFO, ISUPPZ
REAL(8), DIMENSION(Ny,2*NNc) :: EIGS
COMPLEX(8), DIMENSION(2*NNc,Ny,2*NNc) :: MAT
COMPLEX(8), DIMENSION(2*NNc,Ny,2*NNc) :: VECS
INTEGER, DIMENSION(maxval(wi_len,dim=1),Ny) :: IWORK
REAL(8), DIMENSION(maxval(wr_len,dim=1),Ny) :: RWORK
COMPLEX(8), DIMENSION(maxval(w_len,dim=1),Ny) :: WORK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Fk_sum(:,:) = cmplx(0.0,0.0)
MAT(:,:,:) = cmplx(0.0,0.0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC)
do i=1,Ny
	do j=1,Nc(i)
		MAT(j,i,j) = (Kx(Nci(i,j))**2 + Ky(i)**2 - 1)/D0
		do jj=1,Nc(i)
			if ((Kx(Nci(i,j))-Kx(Nci(i,jj)))<=maxval(Kx) .and. (Kx(Nci(i,j))-Kx(Nci(i,jj)))>=minval(Kx)) then
				M(i) = minloc(abs(Kx - (Kx(Nci(i,j))-Kx(Nci(i,jj)))),dim=1)
				MAT(j,i,Nc(i)+jj) = Fk(M(i))
				MAT(j,i,jj) = MAT(j,i,jj)+ Pk(M(i))
				if (flags(1) == 2) then
					theta(i,j,jj) = sign(dble(1),Ky(i))* &
									acos(((Kx(Nci(i,j)) + Kx(Nci(i,jj)))/2)/sqrt(((Kx(Nci(i,j)) + Kx(Nci(i,jj)))/2)**2 + Ky(i)**2))
					MAT(j,i,Nc(i)+jj) = MAT(j,i,Nc(i)+jj)*sin(2*(theta(i,j,jj)+phase))
				end if
			end if
		end do
	end do
	MAT(Nc(i)+1:2*Nc(i),i,Nc(i)+1:2*Nc(i)) = - conjg(MAT(1:Nc(i),i,1:Nc(i)))
	call zheevr('V', 'I', 'U', 2*Nc(i), MAT(1:2*Nc(i),i,1:2*Nc(i)), 2*Nc(i), dble(0), dble(0), Nc(i)+1, Nc(i)+Neigs(i), dble(-1), &
		M(i), EIGS(i,1:2*Nc(i)), VECS(1:2*Nc(i),i,1:2*Nc(i)), 2*Nc(i), ISUPPZ(i), &
		WORK(1:w_len(i),i), w_len(i), RWORK(1:wr_len(i),i), wr_len(i), IWORK(1:wi_len(i),i), wi_len(i), INFO(i))

  do j = 1,Nc(i)
    if (EIGS(i,j) < 0 .or. EIGS(i,j) > E_cut+0.1) then
      write(*,*) "BAD EIGENVALUE!!!!!!!!! ", EIGS(i,j)
    end if
  end do
	EIGS(i,:) = abs(EIGS(i,:))
	F_sc(i,1:Neigs(i)) = (1 - (exp((EIGS(i,1:Neigs(i)) + b)/t)+1)**(-1) - (exp((EIGS(i,1:Neigs(i)) - b)/t)+1)**(-1))
	do j=1,NN/2+1
		do jj=1,Nc(i)
		do jjm=1,Nc(i)
			if (Kx(Nci(i,jjm)) == Kx(Nci(i,jj))-Kx(j)) then
			if (flags(1) == 1) then
				Fk_sum(j,i) = Fk_sum(j,i) &
							+ sum(F_sc(i,1:Neigs(i))*VECS(jj,i,1:Neigs(i))*conjg(VECS(Nc(i)+jjm,i,1:Neigs(i))))
			elseif (flags(1) == 2) then
				Fk_sum(j,i) = Fk_sum(j,i) &
							+ sum(F_sc(i,1:Neigs(i))*VECS(jj,i,1:Neigs(i))*conjg(VECS(Nc(i)+jjm,i,1:Neigs(i)))) &
							 *sin(2*(theta(i,jj,jjm)+ phase))
			end if
			end if
		end do
		end do
	end do
end do
!$OMP END PARALLEL DO
Fk2 = sum(Fk_sum,dim=2)
EIGS_pass = EIGS(:,1:NNc)
END SUBROUTINE get_new_profile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE delta_guess(j,i, DD, rr, D00,N,Nt,Nb,flags)
implicit none
integer :: i, j, N, Nt, Nb, jj
integer, dimension(9) :: flags
REAL(8), dimension(Nt,Nb,N) :: DD
REAL(8), dimension(N) :: rr
REAL(8), dimension(2*N) :: guess_file
REAL(8) :: D00
CHARACTER(len=2048) :: filename

write(filename,'(A,I2,I2,A)') "./SELF_CONSISTENT/sc_delta",j,i, ".dat"
open(1, file=trim(filename), status='old', iostat=jj)
if (flags(6) == 0 .and. jj == 0) then
	read(1,*) guess_file
	do jj=1,N
		DD(j,i,jj) = guess_file(jj*2)
	end do
	close(1)
else
	if (i == 1 .and. j == 1) then
		if (flags(5) < 0) then
			DD(j,i,:) = D00
		elseif (flags(5) == 0) then
			DD(j,i,:) = D00
		elseif (flags(5) > 0) then
			DD(j,i,:) = D00*sin(flags(5)*2*3.141592*rr/dble(N))
		end if
	elseif (i == 1) then
		if (j > 2) then
			DD(j,i,:) = 2*DD(j-1,1,:) - DD(j-2,1,:)
		else
			DD(j,i,:) = DD(j-1,1,:)
		end if
		if (maxval(abs(DD(j-1,1,:))) == 0.0) then
			DD(j,i,:) = 0.0
		end if
	else
		if (i > 2) then
			DD(j,i,:) = 2*DD(j,i-1,:) - DD(j,i-2,:)
		else
			DD(j,i,:) = DD(j,i-1,:)
		end if
		if (maxval(abs(DD(j,i-1,:))) == 0.0) then
			DD(j,i,:) = 0.0
		end if
	endif
end if
end subroutine delta_guess
