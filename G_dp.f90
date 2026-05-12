Program G
implicit none
  integer, parameter :: n=3*10074
  integer, parameter :: m=3*10074
  integer, parameter :: T=300
  integer, parameter :: samp=1
  real(8), parameter :: tu=1e+12
  real(8), dimension(m) :: v, rk
  real(8), dimension(m) :: swp,wp
  real(8), dimension(m) :: gamma
  integer i, j, k, l, p, tej
  real(8) :: sum,  constant, nu, f, w
  CHARACTER(50) F1,F2,F3,F4,F5,F6


rk=0
   do tej=1,4

write(F1,'(a,I1,a,I1,a)')'../T300/runs/run',samp,'/data/',tej,'/volume.dat'
open(unit=9,file=F1,status="old")
write(F2,'(a,I1,a,I1,a)')'../T300/runs/run',samp,'/data/',tej,'/eigenvalues.data'
open(unit=10,file=F2,status="old")
write(F4,'(a,I1,a,I1,a,a)')'../T300/runs/run',samp,'/data/',tej,'/gamma.data'
open(unit=111,file=F4,access='sequential',status='old',action='read')
write(F6,'(a,I1,a,I1,a,I3,a)')'../T300/runs/run',samp,'/data/',tej,'/G_dprime_T',T,'.data'
open(unit=13,file=F6,status="unknown")

do k=1,1
read(9,*) v(k)
end do


do i=1,m
read(10,*) swp(i)
read(111,*) gamma(i)
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
w=(tu)*exp(-15+0.1d0*p)

nu=5e+13 ! Kg/s
constant=(1e+30)*v(1)**(-1)

sum=0.d0
do i=1,m-1
if(abs(swp(i)).lt.+0.00001d0) then
gamma(i)=0.d0
elseif((swp(i)).gt.-2.570d0.and.(swp(i)).lt.+2.570d0) then 
gamma(i)=0!0.1157d0*abs(swp(i))
else
gamma(i)=gamma(i)
endif
f=2.906*(1e+6)*w*gamma(i)*nu*(((tu**2)*swp(i)-w*w)**(+2)+1*nu*nu*w*w)**(-1)
sum=sum+f
end do
write(13,*)w*(tu**(-1)), sum*constant*(1e-9)

end do

end do
14 format(2I7,15f10.2)
15 format(15f10.2)

!end do
end Program G

