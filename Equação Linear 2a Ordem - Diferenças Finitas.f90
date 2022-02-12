!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!***************************************************************!!
!! Objetivo: escrever um programa que resolva uma equação linear !!
!! de segunda ordem pelo método direto das diferenças finitas    !!
!! Aluna: Laura Simonassi Raso de Paiva                          !!
!! Matrícula: 180021885                                          !!
!! Data: 05/11/2019                                              !!
!! Última modificação: 12/11/2019                                !!
!!***************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module vg

real(kind=8),allocatable :: y(:),x(:)
real(kind=8)             :: h
integer                  :: n

end module vg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!***************************************************************!!
!!  Programa que resolve a equação diferencial de segunda ordem  !! 
!!***************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program edosegunda
use vg
implicit none
integer :: i

h=0.1d0
n=int(1.0d0/h)

allocate(x(n),y(n))

x(0)=0.0d0
x(n)=1.0d0

y(0)=0.0d0
y(n)=100.0d0


do i=1,n-1
call resolve_sistema
  if (abs(y(i+1)-y(i))/y(i) <= 1d-10) exit
  write(*,*)'y=',y(i),'x=',x(i) 
  x(i+1)=x(i)+h
  x(i)=x(i+1) 
end do

end program edosegunda
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!****************************************************************!!
!!    Subrotina que resolve o sistema de equações lineares        !!
!!****************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine resolve_sistema
use vg
implicit none
real(kind=8),dimension(n)  :: beta, p
real(kind=8),dimension(n)  :: b, a, c, r
integer:: i,j

 b(j)= 10.0d0 - 2.0d0/(h**2.0d0) 
 a(j)= 1.0d0/(h**2.0d0) + 5.0d0/(2.0d0*h)
 c(j)= 1.0d0/(h**2.0d0) - 5.0d0/(2.0d0*h)

do j=2,n-2
r(j)= 10*x(j)
enddo

r(1)=10*x(1)-(1.0d0/(h**2.0d0) + 5.0d0/(2.0d0*h))*y(0)
r(n-1)=10*x(n-1)-(1.0d0/(h**2.0d0) - 5.0d0/(2.0d0*h))*y(n)

beta(1)=b(1)
p(1)=r(1)

do j=2,n-1
beta(j)=b(j)-a(j)/beta(j-1)*c(j-1)
end do

 do j=2,n-1
p(j)=r(j)-a(j)*p(j-1)/beta(j-1)
end do

y(n)=p(n)/beta(n)

do j=1,n-1
y(n-j)=(p(n-j)-c(n-j)*y(n-j+1))/beta(n-j)
end do

end subroutine resolve_sistema

