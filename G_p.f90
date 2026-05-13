Program G   ! Fotran codes for integration
implicit none
  integer, parameter :: n=3*10074   ! N=10074 Number of monomer
  integer, parameter :: m=3*10074
  integer, parameter :: T=300  ! Temperature
  integer, parameter :: samp=2
  real(8), parameter :: tu=1e+12 ! frequency in THz 
  real(8), dimension(m) :: v, rk, Gx,Gy,Gz
  real(8), dimension(m) :: swp,wp
  real(8), dimension(m) :: gamma
  integer i, j, k, l, p, tej
  real(8) :: sum,  constant, nu, f, w, GA
  CHARACTER(50) F0,F1,F2,F3,F4,F5,F6


rk=0
   do tej=1,4  !trejectories started 1 to 4 runs

   write(F0,'(a,I1,a,I1,a)')'../T300/runs/run',samp,'/data/',tej,'/GA.dat'
open(unit=8,file=F0,status="old")
write(F1,'(a,I1,a,I1,a)')'../T300/runs/run',samp,'/data/',tej,'/volume.dat'
open(unit=9,file=F1,status="old")
write(F2,'(a,I1,a,I1,a)')'../T300/runs/run',samp,'/data/',tej,'/eigenvalues.data'
open(unit=10,file=F2,status="old")
write(F4,'(a,I1,a,I1,a,a)')'../T300/runs/run',samp,'/data/',tej,'/gamma.data'
open(unit=111,file=F4,access='sequential',status='old',action='read')
write(F6,'(a,I1,a,I1,a,I3,a)')'../T300/runs/run',samp,'/data/',tej,'/G_prime_T',T,'.data'
open(unit=13,file=F6,status="unknown")

do k=1,1
read(9,*) v(k)  ! volume
read(8,*) Gx(k),Gy(k),Gz(k)  !Affine modulus
end do


do i=1,m
read(10,*) swp(i)    !eigenfrequency in term of omega^2
read(111,*) gamma(i) !Affine force field correlator
end do
!++++++++++++++++++++++++++++

do i=1,m
   if(swp(i).le.0.d0) then
   wp(i)=-(sqrt(-swp(i)))
else
   wp(i)=  sqrt(swp(i))
end if
end do
!++++++++++++++++++++++++++++
do p=1,235
w=(tu)*exp(-15+0.1d0*p) !External frequency in term of real unit

nu=5e+13 ! Kg/s  ! Kg/s friction kernal
constant=(1e+30)*v(1)**(-1)
GA= (101325*(Gx(1)+Gy(1)+Gz(1))/3)  ! Affine modulus

sum=0.d0
do i=1,m
if(abs(swp(i)).lt.+0.00001d0) then  !Avoid Goldstone modes
gamma(i)=0.d0
elseif((swp(i)).gt.-2.570d0.and.(swp(i)).lt.+2.570d0) then  ! Discart Unphysical modes belwo cutoff frequency or qudratic fit
gamma(i)=0!0.1157d0*abs(swp(i))
else
gamma(i)=gamma(i)
endif
f=-2.906*(1e+6)*((tu**2)*swp(i)-w*w)*gamma(i)*(((tu**2)*swp(i)-w*w)**(+2)+1*nu*nu*w*w)**(-1)  !Integration perform in term of sum and in proper SI unit
sum=sum+f
end do
write(13,*)w*(tu**(-1)), (1*(GA)+sum*constant)*(1e-9)  !modulus results as a function of external frequency
end do
write(*,15) (101325*(Gx(1)+Gy(1)+Gz(1))/3) ! write on screen to check
end do
14 format(2I7,15f10.2)
15 format(2e15.7)

!end do
end Program G

