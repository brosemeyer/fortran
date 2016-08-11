program MAIN
implicit none

! TODO clean up PI usage
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(8), parameter :: PI=3.141592654
!!!!!!!!!!!!!!!INTEGERS
INTEGER, parameter :: Nb=1, Nt=10, NL=1, NN = 512, N_norm = 1000, NE = 101
!!!!!!!!!!!!!!!BOUNDS
INTEGER, parameter :: Lmin = 0, Lmax = 400
REAL(8), parameter :: t_min = 0.05, t_max = 0.45
REAL(8), parameter :: b_min = 0.3, b_max = 0.8
REAL(8), parameter :: E_cut = 5, r_range = 150
!!!!!!!!!!!!!!!SELF CONSISTENCY
INTEGER, parameter :: max_it = 200
REAL(8), parameter :: tol = 1e-3, bzones = 1.0
!!!!!!!!!!!!!!!ORDER PARAMETER
REAL(8), parameter :: phase = 0*PI/8, D0 = 0.05, D00 = 1.0 !First guess
!!!!!!!!!!!!!!!FREQUENCY
REAL(8), parameter :: eta = 2.5e-3
!!!!!!!!!!!!!!!COMPUTED
INTEGER, parameter :: Ny = 2*int(NN*sqrt(1+E_cut*D0)/(bzones*2*PI))
!!!!!!!!!!!!!!!FLAGS
INTEGER, DIMENSION(9) :: flags = (/1, & !	1=S wave, 2=D wave
								   0, & !	1=calculate real part, 0=no
								   0, & !	1=calculate relaxation rate, 0=no
								   0, & !	1=calculate DOS, 0=no save
								   -1, & !	0=0 DW, 1=2 DW, 2=4 DW, 3=6 DW, 4=8 DW, 5=10 DW
								   1, & !	1=Calculate self consistent, 0=use data file guess
								   0, & !	1=Calculate Homogeneous, 0=no homogeneous
								   0, &	!	1=Calculate Free energy, 0=no free energy
								   0  &	!	1=Calculate Magnetization, 0=no Magnetization
								   /)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER, DIMENSION(NL) :: L_sn
REAL(8), DIMENSION(NN,2) :: chi_norm
REAL(8), DIMENSION(NN) :: Kx, q0, rr
REAL(8), DIMENSION(Ny) :: Ky
INTEGER, DIMENSION(Ny) :: Nc, w_len, wr_len, wi_len
INTEGER, DIMENSION(Ny,NN) :: Nci
LOGICAL, DIMENSION(NN/2+1) :: msk_space
REAL(8), DIMENSION(Nt) :: TT
REAL(8), DIMENSION(Nb) :: HH
REAL(8), DIMENSION(Nt,Nb,NN) :: DD
REAL(8), DIMENSION(Nt,Nb,3) :: Free_E
REAL(8), DIMENSION(Nt,Nb,NL,5) :: SUS0
REAL(8), DIMENSION(Nt,Nb,NL) :: rel_time
REAL(8) :: sc_const00
INTEGER :: it, ib, i, j, NNc
character(8)  :: date
character(10) :: time
INTEGER, DIMENSION(3) :: time_s, time_f, date_s, date_f
character(5)  :: zone
integer,DIMENSION(8) :: values
CHARACTER(len=2048) :: filename

call date_and_time(date,time,zone,values)
read(time(1:2),'(I2)') time_s(1)
read(time(3:4),'(I2)') time_s(2)
read(time(5:6),'(I2)') time_s(3)
read(date(1:4),'(I2)') date_s(1)
read(date(5:6),'(I2)') date_s(2)
read(date(7:8),'(I2)') date_s(3)

call build_arrays()
! temp and field loop
do it=1,Nt
 TT(it) = t_min + (t_max-t_min)*dble(it-1)/dble(Nt-1)
 if (Nt == 1) then
	 TT(1) = t_min
 end if
 do ib=1,Nb
	 HH(ib) = b_min + (b_max-b_min)*dble(ib-1)/dble(Nb-1)
	 if (Nb == 1) then
		 HH(1) = b_min
	 end if
	 call calculations()
 end do
 if (Nb /= 1) then
	 write(filename,'(A,I2,A)') "./RELAXATION/rel_time_fixedT",it, ".dat"
	 call save_file_3(HH,Nb,dble(L_sn),NL,1,filename,rel_time(it,:,:))
	 write(filename,'(A,I2,A)') "./SUS/SUS0_fixedT",it, ".dat"
	 call save_file_3(HH,Nb,dble(L_sn),NL,5,filename,SUS0(it,:,:,:))
 end if
end do

call save_files()

call date_and_time(date,time,zone,values)
read(time(1:2),'(I2)') time_f(1)
read(time(3:4),'(I2)') time_f(2)
read(time(5:6),'(I2)') time_f(3)
read(date(1:4),'(I2)') date_f(1)
read(date(5:6),'(I2)') date_f(2)
read(date(7:8),'(I2)') date_f(3)

write(*,*) ((time_f(1)-time_s(1))*3600 &
		 + (time_f(2)-time_s(2))*60 &
		 + (time_f(3)-time_s(3)) &
		 + (date_f(3)-date_s(3))*24*60*60), ' seconds'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
	!CALCULATIONS
	SUBROUTINE calculations()
		implicit none
		INTEGER :: success
		REAL(8), DIMENSION(Ny,2*NNc) :: EIGS
		COMPLEX(8), DIMENSION(2*NNc,Ny,2*NNc) :: VECS
		REAL(8), dimension(Ny,NNc,2) :: Ep, Fp
		COMPLEX(8), DIMENSION(NNc,Ny,NN) :: Uxp, Vxp
		REAL(8), DIMENSION(NN,2) :: Mag

		if (flags(7) == 1) then
 			call Normal_sus(TT(it),HH(ib),D0,q0/4,NN,N_norm,chi_norm,E_cut,flags)
 	  end if

 	  call delta_guess(it,ib,DD,rr,D00,NN,Nt,Nb,flags)
 	  if (maxval(abs(DD(it,ib,:))) == 0.0) then
			DD(j,i,:) = 0.0
 		  rel_time(j,i,:) = 0.0
		  SUS0(j,i,:,:) = 0.0
		else
		call self_consistent(DD(it,ib,:), Kx, Ky, rr, NN, Ny, TT(it), HH(ib), D0, flags,phase,VECS, EIGS, sc_const00, &
													it, ib, success,w_len,wr_len,wi_len,Nc,Nci,NNc,tol,max_it,E_cut)
		Ep(:,:,2) = D0*(EIGS(:,1:NNc) + HH(ib))
		Ep(:,:,1) = D0*(EIGS(:,1:NNc) - HH(ib))
		Fp = (exp(Ep/(TT(it)*D0))+1)**(-1)
		do i=1,Ny
			do j=1,Nc(i)
				call Fourier_q(Uxp(j,i,:),VECS(1:Nc(i),i,j),NN,Nc(i),Nci(i,1:Nc(i)),2)
				call Fourier_q(Vxp(j,i,:),VECS(Nc(i)+1:2*Nc(i),i,j),NN,Nc(i),Nci(i,1:Nc(i)),2)
			end do
		end do
		if (flags(2) == 1) then
			call Suscept(TT(it),HH(ib),SUS0(it,ib,:,1:2),D0,rr,Kx,Ky,q0/4,NN,Ny, &
									 L_sn,NL,it,ib,chi_norm,Nc,Nci(:,1:NNc),NNc,msk_space, &
									 Ep, Fp, Uxp, Vxp)
		end if
		if (flags(9) == 1) then
			call Magnetization(D0,TT(it),HH(ib),Fp,Uxp,Vxp,NN,Nc,NNc,Nci,Ny,kx,ky,Mag,it,ib,rr)
			SUS0(it,ib,:,3) = Mag(1,1)
			SUS0(it,ib,:,4) = Mag(NN/4,1)
			SUS0(it,ib,:,5) = Mag(1,2)
		end if
		if (flags(3) == 1) then
			do i=1,NL
				call get_relaxation(Fp,Ep,Uxp(:,:,L_sn(i)+1),Vxp(:,:,L_sn(i)+1),Ny,NNc,NN,Nc,rel_time(it,ib,i),TT(it),eta)
			end do
		end if
		if (flags(8) == 1) then
			call Free_Energy(D0,DD,TT(it),HH(ib),Ep,Uxp,Vxp,NN,Nc,NNc,Nci,Ny,kx,ky, &
											sc_const00,Free_E(it,ib,:),rr,it,ib)
		end if
		if (flags(4) == 1) then
			call DOS(Uxp, Vxp, Ep, Ky, Kx, rr, D0, Nc, NN,NNc,Nci,Ny,NE,E_cut,eta,it,ib)
		end if

		end if
	END SUBROUTINE calculations

  ! MAKE ARRAYS
	SUBROUTINE build_arrays()
	 do j=1,NN
		 rr(j) = dble(j-1)/bzones
		 q0(j) = bzones*(2*PI)*dble(j-1)/dble(NN)
	 end do
	 msk_space = (2*abs(rr(1:NN/2+1)-rr(NN/4+1)) <= r_range)
	 do i=1,Ny/2
		 Ky(Ny/2+i) = sqrt(1+E_cut*D0)*dble(i-0.5)/dble(Ny/2)
		 Ky(Ny/2+1-i) = -Ky(Ny/2+i)
	 end do
	 do i=1,NL
		 if (NL /= 1) then
			 L_sn(i) = Lmin + int((Lmax-Lmin)*dble(i-1)/dble(NL-1))
		 else
			 L_sn(i) = Lmin
		 end if
	 end do
	 Kx(1:NN/2+1) = q0(1:NN/2+1)
	 kx(NN/2+2:NN) = -q0(NN/2:2:-1) ! CAN DO FOR FFT
	 Nc = 0
	 do i=1,Ny
	 do j=1,NN
		 if (E_cut*D0-abs(kx(j)**2+Ky(i)**2-1)>=0) then
			 Nc(i) = Nc(i) + 1
			 Nci(i,Nc(i)) = j
		 end if
	 end do
	 if (Nc(i)<=2) then
		 write(*,*) 'small number of momenta in energy range'
	 end if
	 end do
	 NNc = maxval(Nc,dim=1)
	 !initialize sc
	 call self_consistent_00(kx,Ky,NN,Ny,D0,flags,phase,sc_const00, &
	 							w_len,wr_len,wi_len,Nc,Nci(:,1:NNc),NNc,E_cut)
	call system('rm -rf DOS')
	call system('rm -rf SUS')
	call system('rm -rf RELAXATION')
	call system('mkdir DOS')
	call system('mkdir SUS')
	call system('mkdir RELAXATION')
	if (flags(6) == 1) then
		call system('rm -rf SELF_CONSISTENT')
		call system('mkdir SELF_CONSISTENT')
	end if
 END SUBROUTINE build_arrays

  !SAVE FILES
 SUBROUTINE save_files()
	 if (Nt /= 1 .or. (Nt == 1 .and. Nb == 1)) then
	 	do i=1,Nb
	 		write(filename,'(A,I2,A)') "./RELAXATION/rel_time_fixedB",i, ".dat"
	 		call save_file_3(TT,Nt,dble(L_sn),NL,1,filename,rel_time(:,i,:))
	 		write(filename,'(A,I2,A)') "./SUS/SUS0_fixedB",i, ".dat"
	 		call save_file_3(TT,Nt,dble(L_sn),NL,5,filename,SUS0(:,i,:,:))
	 	end do
	 end if
	 if (Nt /= 1 .and. Nb /= 1) then
	 	do i=1,NL
	 		write(filename,'(A,I2,A)') "./RELAXATION/rel_time_fixedL",i, ".dat"
	 		call save_file_3(TT,Nt,HH,Nb,1,filename,rel_time(:,:,i))
	 		write(filename,'(A,I2,A)') "./SUS/SUS0_fixedL",i, ".dat"
	 		call save_file_3(TT,Nt,HH,Nb,5,filename,SUS0(:,:,i,:))
	 	end do
	 end if
	 call save_file_3(TT,Nt,HH,Nb,1,'./SELF_CONSISTENT/delta_max.dat',maxval(DD,dim=3))
	 if (flags(8) == 1) then
	 	call save_file_3(TT,Nt,HH,Nb,3,'./SELF_CONSISTENT/Free_E.dat',Free_E)
	 end if
 END SUBROUTINE save_files
END PROGRAM MAIN
