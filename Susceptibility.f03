!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Suscept(t,b,SUS0,D0,rr,Kx,Ky,q_norm,NN,Ny,L_sn, &
				   NL,i_t,i_b, &
				   chi_norm,Nc,Nci,NNc,msk_space, &
					 Ep, Fp, Uxp, Vxp)

implicit none
INTEGER :: NN, NL, Ny, NNc
REAL(8), PARAMETER :: PI=3.141592654
INTEGER :: i, j, ii, i_b, i_t
INTEGER, DIMENSION(Ny) :: Nc
INTEGER, DIMENSION(Ny,NNc) :: Nci
LOGICAL, DIMENSION(NN/2+1) :: msk_space
INTEGER, DIMENSION(NL) :: L_sn
REAL(8), DIMENSION(NL,2) :: SUS0
REAL(8), DIMENSION(Ny/2) :: qy
REAL(8), DIMENSION(NN) ::  Kx, q_norm, rr, Kx_p
REAL(8), DIMENSION(Ny) :: Ky
REAL(8), DIMENSION(NN,2) :: chi_norm
REAL(8), DIMENSION(NN,Ny/2,2) :: Chi_xxq, chi_zzq
COMPLEX(8), DIMENSION(NN/2+1,Ny/2) :: Chi_xx, chi_zz
REAL(8) :: t, b, D0
CHARACTER(len=2048) :: filename
REAL(8), dimension(Ny,NNc,2) :: Ep, Fp
COMPLEX(8), DIMENSION(NNc,Ny,NN) :: Uxp, Vxp

	qy = ky(Ny/2+1:Ny)-abs(ky(1)-ky(2))/2
	qy(1) = 0.0
	Kx_p = (/Kx(NN/2+2:NN),Kx(1:NN/2+1)/)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,NL
	Chi_xx(:,:) = cmplx(0.0,0.0)
	Chi_xxq(:,:,:) = 0.0
	Chi_zz(:,:) = cmplx(0.0,0.0)
	Chi_zzq(:,:,:) = 0.0
!$OMP PARALLEL DO SCHEDULE(DYNAMIC)
do j=1,Ny/2
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call get_suscept_xx(j,Fp,Ep,Uxp,Vxp,Ny,NNc,NN,Nc,L_sn(i)+1,msk_space,chi_xx(:,j),chi_zz(:,j))
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	chi_xx(:,j) = chi_xx(:,j)*4*abs(ky(1)-ky(2))/dble(NN)
	chi_zz(:,j) = chi_zz(:,j)*4*abs(ky(1)-ky(2))/dble(NN)
	do ii=1,NN
		if (sqrt(Kx_p(ii)**2 + (2*qy(j))**2) <= 2*qy(Ny/2)) then
			chi_xxq(ii,j,1) = sum(real(chi_xx(:,j)*exp(-cmplx(0,2*(rr(1:NN/2+1)-rr(NN/4+1))*kx_p(ii)))))/dble(NN)
			chi_xxq(ii,j,2) = chi_norm(minloc(abs(sqrt(Kx_p(ii)**2 + (2*qy(j))**2) - q_norm),dim=1),1)
			chi_zzq(ii,j,1) = sum(real(chi_zz(:,j)*exp(-cmplx(0,2*(rr(1:NN/2+1)-rr(NN/4+1))*kx_p(ii)))))/dble(NN)
			chi_zzq(ii,j,2) = chi_norm(minloc(abs(sqrt(Kx_p(ii)**2 + (2*qy(j))**2) - q_norm),dim=1),2)
		end if
	end do
end do !!!!!!!!!!!!!!!!!!!!!! END Qy j
!$OMP END PARALLEL DO
	SUS0(i,1) = chi_zzq(NN/2,2,1)
	SUS0(i,2) = chi_xxq(NN/2,2,1)
	write(filename,'(A,I2,I2,I2,A)') "./SUS/SUS_qx_LHT", i, i_b, i_t, ".dat"
	call save_file_3(Kx_p,NN,2*qy,Ny/2,2,filename,chi_xxq)
	write(filename,'(A,I2,I2,I2,A)') "./SUS/SUS_qz_LHT", i, i_b, i_t, ".dat"
	call save_file_3(Kx_p,NN,2*qy,Ny/2,2,filename,chi_zzq)
end do  !!!!!!!!!!!!!!!!!!!!!!! END POSITIONS i

end subroutine Suscept
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_suscept_xx(j,Fp,Ep,Uxp,Vxp,Ny,NNc,NN,Nc,i_R,msk_space,chix,chiz)
implicit none
INTEGER :: Ny, NNc, NN, i_R
INTEGER :: jj, k, ip, im, j
INTEGER :: i1, i2, i1m, i2m
INTEGER, DIMENSION(Ny) :: Nc
COMPLEX(8), DIMENSION(NN/2+1) :: chix, chiz
REAL(8), DIMENSION(Ny,NNc,2) :: Ep, Fp
LOGICAL, DIMENSION(NN/2+1) :: msk_space
COMPLEX(8), DIMENSION(NNc,Ny,NN) :: Uxp, Vxp
REAL(8), DIMENSION(NNc,NNc,3,2) :: F_terms
COMPLEX(8), DIMENSION(NNc,NNc,3) :: AA

do i1 = 1, Ny-2*(j-1)
	i2 = i1 + 2*(j-1)
	i1m = Ny+1-i1
	i2m = Ny+1-i2
	do jj = 1,Nc(i1)
		F_terms(jj,1:Nc(i2m),1,1) = (Fp(i1,jj,1)-1+Fp(i2m,1:Nc(i2m),1)) &
						       /((Ep(i1,jj,1) + Ep(i2m,1:Nc(i2m),1))) &
						      + (Fp(i1,jj,2)-1+Fp(i2m,1:Nc(i2m),2)) &
						       /((Ep(i1,jj,2) + Ep(i2m,1:Nc(i2m),2)))
		F_terms(jj,1:Nc(i2),2,1) = (Fp(i1m,jj,1)-1+Fp(i2,1:Nc(i2),1)) &
						       /((Ep(i1m,jj,1) + Ep(i2,1:Nc(i2),1))) &
						      + (Fp(i1m,jj,2)-1+Fp(i2,1:Nc(i2),2)) &
						       /((Ep(i1m,jj,2) + Ep(i2,1:Nc(i2),2)))
		F_terms(jj,1:Nc(i2),3,1) = (Fp(i1,jj,1)-Fp(i2,1:Nc(i2),2)) &
						       /((Ep(i1,jj,1) - Ep(i2,1:Nc(i2),2))) &
						      + (Fp(i1,jj,2)-Fp(i2,1:Nc(i2),1)) &
						       /((Ep(i1,jj,2) - Ep(i2,1:Nc(i2),1)))
		F_terms(jj,1:Nc(i2m),1,2) = (Fp(i1,jj,1)-1+Fp(i2m,1:Nc(i2m),2)) &
						       /((Ep(i1,jj,1) + Ep(i2m,1:Nc(i2m),2))) &
						      + (Fp(i1,jj,2)-1+Fp(i2m,1:Nc(i2m),1)) &
						       /((Ep(i1,jj,2) + Ep(i2m,1:Nc(i2m),1)))
		F_terms(jj,1:Nc(i2),2,2) = (Fp(i1m,jj,1)-1+Fp(i2,1:Nc(i2),2)) &
						       /((Ep(i1m,jj,1) + Ep(i2,1:Nc(i2),2))) &
						      + (Fp(i1m,jj,2)-1+Fp(i2,1:Nc(i2),1)) &
						       /((Ep(i1m,jj,2) + Ep(i2,1:Nc(i2),1)))
		F_terms(jj,1:Nc(i2),3,2) = (Fp(i1,jj,1)-Fp(i2,1:Nc(i2),1)) &
						       /((Ep(i1,jj,1) - Ep(i2,1:Nc(i2),1))) &
						      + (Fp(i1,jj,2)-Fp(i2,1:Nc(i2),2)) &
						       /((Ep(i1,jj,2) - Ep(i2,1:Nc(i2),2)))
	end do
	do k=0,NN/4
		if (msk_space(k+NN/4+1)) then
			ip = i_R + k
			im = i_R - k
			if (ip > NN) then
				ip = ip - NN
			elseif (ip < 1) then
				ip = ip + NN
			end if
			if (im < 1) then
				im = im + NN
			elseif (im > NN) then
				im = im - NN
			end if
			do jj = 1,Nc(i1)
AA(jj,1:Nc(i2m),1) = conjg(Uxp(jj,i1,ip)*Vxp(1:Nc(i2m),i2m,ip) - Vxp(jj,i1,ip)*Uxp(1:Nc(i2m),i2m,ip)) * &
					      (Uxp(jj,i1,im)*Vxp(1:Nc(i2m),i2m,im) - Vxp(jj,i1,im)*Uxp(1:Nc(i2m),i2m,im))/2
AA(jj,1:Nc(i2),2) =      (Vxp(jj,i1m,ip)*Uxp(1:Nc(i2),i2,ip) - Uxp(jj,i1m,ip)*Vxp(1:Nc(i2),i2,ip)) * &
					conjg(Vxp(jj,i1m,im)*Uxp(1:Nc(i2),i2,im) - Uxp(jj,i1m,im)*Vxp(1:Nc(i2),i2,im))/2
AA(jj,1:Nc(i2),3) = (conjg(Uxp(jj,i1,ip))*Uxp(1:Nc(i2),i2,ip) + conjg(Vxp(jj,i1,ip))*Vxp(1:Nc(i2),i2,ip)) * &
					(Uxp(jj,i1,im)*conjg(Uxp(1:Nc(i2),i2,im)) + Vxp(jj,i1,im)*conjg(Vxp(1:Nc(i2),i2,im)))
			end do
			chix(k+NN/4+1) = chix(k+NN/4+1) - sum(AA(1:Nc(i1),1:Nc(i2m),1)*F_terms(1:Nc(i1),1:Nc(i2m),1,1) &
													 +AA(1:Nc(i1m),1:Nc(i2),2)*F_terms(1:Nc(i1m),1:Nc(i2),2,1) &
													 +AA(1:Nc(i1),1:Nc(i2),3)*F_terms(1:Nc(i1),1:Nc(i2),3,1) &
													 )
			chiz(k+NN/4+1) = chiz(k+NN/4+1) - sum(AA(1:Nc(i1),1:Nc(i2m),1)*F_terms(1:Nc(i1),1:Nc(i2m),1,2) &
													 +AA(1:Nc(i1m),1:Nc(i2),2)*F_terms(1:Nc(i1m),1:Nc(i2),2,2) &
													 +AA(1:Nc(i1),1:Nc(i2),3)*F_terms(1:Nc(i1),1:Nc(i2),3,2) &
													 )
		end if
	end do! end k position
end do ! end ii y momentum
chix(NN/4:1:-1) = conjg(chix(NN/4+2:NN/2+1))
chiz(NN/4:1:-1) = conjg(chiz(NN/4+2:NN/2+1))
end subroutine get_suscept_xx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Magnetization(D0,t,b,Fp,Uxp,Vxp,NN,Nc,NNc,Nci,Ny,kx,ky,Mag,it,ib,rr)
	implicit none
INTEGER :: NN, Ny, NNc, i, j, it, ib
INTEGER, DIMENSION(Ny) :: Nc
REAL(8) :: D0, t, b
INTEGER, DIMENSION(Ny,NNc) :: Nci
REAL(8), dimension(Ny,NNc,2) :: Fp
REAL(8), DIMENSION(Ny,2*NNc) :: EIGS
COMPLEX(8), DIMENSION(NNc,Ny,NN) :: Uxp, Vxp
REAL(8), DIMENSION(NN) :: Kx, rr
REAL(8), DIMENSION(Ny) :: Ky, Eny
REAL(8), DIMENSION(NN,Ny,2) :: Fn
REAL(8), DIMENSION(NN,2) :: Mag
CHARACTER(len=2048) :: filename

Mag(:,:) = 0.0
do i=1,Ny
do j=1,NN
	Mag(j,1) = sum(abs(Uxp(1:Nc(i),i,j))**2*Fp(i,1:Nc(i),1) - abs(Uxp(1:Nc(i),i,j))**2*Fp(i,1:Nc(i),2) &
			     + abs(Vxp(1:Nc(i),i,j))**2*(1-Fp(i,1:Nc(i),2)) - abs(Vxp(1:Nc(i),i,j))**2*(1-Fp(i,1:Nc(i),1))) &
			 + Mag(j,1)
end do
end do
do i=1,NN
	Eny = Kx(i)**2 + Ky(:)**2 - 1
	Fn(i,:,1) = (exp((Eny-b*D0)/(t*D0))+1)**(-1)
	Fn(i,:,2) = (exp((Eny+b*D0)/(t*D0))+1)**(-1)
end do
Mag(:,2) = sum(Fn(:,:,1)-Fn(:,:,2))
write(filename,'(A,I2,I2,A)') "./SELF_CONSISTENT/Mag",it,ib, ".dat"
call save_file_2(rr,NN,2,filename,Mag)
end subroutine Magnetization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_relaxation(Fp,Ep,Uxp,Vxp,Ny,NNc,NN,Nc,rel_time,t,eta)
implicit none
INTEGER :: Ny, NNc, NN
INTEGER :: jj
INTEGER :: i1, i2
INTEGER, DIMENSION(Ny) :: Nc
REAL(8) :: t, eta, rel_time
REAL(8), DIMENSION(Ny,NNc,2) :: Ep, Fp
COMPLEX(8), DIMENSION(NNc,Ny) :: Uxp, Vxp
REAL(8), DIMENSION(Ny,NNc,2) :: DF
DF(:,:,1) = Fp(:,:,1)*(1-Fp(:,:,1))/t
DF(:,:,2) = Fp(:,:,2)*(1-Fp(:,:,2))/t
rel_time = dble(0.0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) REDUCTION(+:rel_time)
do i1=1,Ny
	do jj = 1,Nc(i1)
	do i2=1,Ny
rel_time = rel_time + sum(abs(conjg(Uxp(jj,i1))*Uxp(1:Nc(i2),i2) + conjg(Vxp(jj,i1))*Vxp(1:Nc(i2),i2))**2 &
				         *(DF(i1,jj,1)/((Ep(i1,jj,1) - Ep(i2,1:Nc(i2),2))**2/eta**2+1) &
				          +DF(i1,jj,2)/((Ep(i1,jj,2) - Ep(i2,1:Nc(i2),1))**2/eta**2+1)) &
						+abs(Uxp(jj,i1)*Vxp(1:Nc(i2),i2) - Vxp(jj,i1)*Uxp(1:Nc(i2),i2))**2 &
				         *(DF(i1,jj,1)/((Ep(i1,jj,1) + Ep(i2,1:Nc(i2),1))**2/eta**2+1) &
				          +DF(i1,jj,2)/((Ep(i1,jj,2) + Ep(i2,1:Nc(i2),2))**2/eta**2+1)) &
				          )
	end do
	end do
end do
!$OMP END PARALLEL DO
end subroutine get_relaxation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Normal_sus(t,b,D0,q,NN,N,chi,E_cut,flags)

	implicit none
REAL(8), PARAMETER :: PI=3.141592654
INTEGER :: i,j, ii, N, NN, k
integer, dimension(9) :: flags
REAL(8), DIMENSION(NN) :: q
REAL(8), DIMENSION(N,N) :: Kx, Ky, em, ep
LOGICAL, DIMENSION(N,N) :: msk
REAL(8), DIMENSION(NN,2) :: chi
REAL(8) :: E_cut, t, b, D0, DF!, em, ep

do ii=1,N
	Kx(ii,:)= PI*dble(ii-1)/dble(N)
	Ky(:,ii)= PI*dble(ii-1)/dble(N)
end do
!$OMP PARALLEL DO SCHEDULE(DYNAMIC)
do k=1,NN
  em = (Kx-q(k)/2)**2 + Ky**2 - 1
  ep = (Kx+q(k)/2)**2 + Ky**2 - 1
  msk = (abs(em) <  E_cut*D0 .and. abs(ep) <  E_cut*D0)
  chi(k,1) = sum(((exp((em/D0-b)/t)+1)**(-1)-(exp((ep/D0+b)/t)+1)**(-1))/(em-ep-2*b*D0) &
  		   + ((exp((em/D0+b)/t)+1)**(-1)-(exp((ep/D0-b)/t)+1)**(-1))/(em-ep+2*b*D0), MASK=msk)
  msk = (abs(em) <  E_cut*D0 .and. abs(ep) <  E_cut*D0 .and. em-ep /= 0.0)
  chi(k,2) = sum(((exp((em/D0-b)/t)+1)**(-1)-(exp((ep/D0-b)/t)+1)**(-1))/(em-ep) &
  		   + ((exp((em/D0+b)/t)+1)**(-1)-(exp((ep/D0+b)/t)+1)**(-1))/(em-ep), MASK=msk)
  end do
!$OMP END PARALLEL DO
chi = -2*(Kx(1,1)-Kx(2,1))**2*chi/PI
call save_file_2(q,NN,2, "./SUS/sus_norm.dat",chi)

end subroutine Normal_sus
