!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!*********************************************************!!
!! Objetivo: escrever um sistema que resolva a equação     !!
!! diferencial y''-5y'+10y=10x pelo método de Gauss Seidel !!
!!                                                         !!
!! Aluna: Laura Simonassi Raso de Paiva                    !!
!! Matrícula: 180021885                                    !!
!! Data: 10/11/2019                                        !!
!! Última modificação: 26/11/2019                          !! 
!!*********************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program gauss_seidel

implicit none
real(kind=8), allocatable :: yini(:), yfin(:),x(:)
real(kind=8):: h,erro
integer:: n
integer:: i

open(unit=1, file='resultados.txt')

write(*,*) 'Escolha um valor para n'
read(*,*) n

allocate(x(n),yini(n),yfin(n))

x(1)=0.0d0
h=1.0d0/dfloat(n-1)
do i=1,n-1
  x(i+1)=x(i)+h
end do
yini=100.0d0*x

yfin=yini

do 
  erro=0.0d0
  do i=2,n-1
     yfin(i)=((1.0d0-(5.0d0*h)/2.0d0)*yini(i+1) + (1.0d0+(5.0d0*h)/2.0d0)*yfin(i-1) - 10.0d0*(h**2.0d0)*x(i))/(2.0d0-10.0d0*(h**2.0d0))
     erro=erro+(yfin(i)-yini(i))**2
     yini(i)=yfin(i)
  end do
  erro=erro/sum(yini**2)
    if(erro <= 1.0d-10) exit
end do
 
do i =1,n
   write(1,*) x(i),yfin(i)
end do

close(1)
end program gauss_seidel
