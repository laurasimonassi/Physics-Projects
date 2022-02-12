!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!****************************************************************!!
!! Objetivo: resolver uma equação diferencial ordinária           !!
!! usando os métodos de Euler, Euler modificado e série de Taylor !!
!!                                                                !!
!! Aluna: Laura Simonassi Raso de Paiva                           !!
!! Matrícula: 180021885                                           !!
!! Data: 30/10/2019                                               !!
!! Ultima modificação: 31/10/2019                                 !!
!!****************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!****************************************************************!!
!!            Programa que calcula a solução da edo               !!
!!****************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
program edo
implicit none

open(unit=1, file='resultados.txt', status='new')
call euler
call euler_modificado
call taylor
close(1)

open(unit=10, file='tabelaeuler.txt', status='new')
call euler
close(10)

open(unit=20, file='tabelaeulermodificado.txt', status='new')
call euler_modificado
close(20)

open(unit=30, file='tabelataylor1.txt', status='new')
call taylor
close(30)
end program edo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!****************************************************************!!
!!                      Edo a ser resolvida                       !!
!!****************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
function derivada(x,y)
implicit none
real(kind=8):: x,y, derivada
derivada=-(x*y)
return 
end function derivada
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!****************************************************************!!
!!        Subrotina que calcula a edo pelo método de Euler        !!
!!****************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
subroutine euler
implicit none
real(kind=8):: b,a,h, n
real(kind=8):: x,y, derivada
integer:: i
external derivada

a=0.0d0
b=1.0d0
y=1.0d0
h=0.1d0     !Variamos o valor
n=((b-a)/h)
x=a
do i=0,n-1
 y=y+h*derivada(x,y)
 x=x+h
write(1,*)'Os valores obtidos pelo método de Euler foram x=',x,'e y=',y
write(10,*)x,y
end do

return
end subroutine euler
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!****************************************************************!!
!!  Subrotina que calcula a edo pelo método de Euler modificado   !!
!!****************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
subroutine euler_modificado
implicit none
real(kind=8):: b,a,h, n
real(kind=8):: x,y, derivada
integer:: i
external derivada

a=0.0d0
b=1.0d0
y=1.0d0
h=0.1d0     !Variamos o valor
n=((b-a)/h)
x=a
do i=0,n-1
 y=y+h*derivada((x+h/2),y+(h/2)*derivada(x,y))
 x=x+h
write(1,*)'Os valores obtidos pelo método de Euler modificado foram x=',x,'e y=',y
write(20,*)x,y
end do

return
end subroutine euler_modificado
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!****************************************************************!!
!!  Subrotina que calcula a edo pelo método de série de Taylor    !!
!!****************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
subroutine taylor
implicit none
real(kind=8):: b,a,h, n
real(kind=8):: x,y, derivada
integer:: i
external derivada

a=0.0d0
b=1.0d0
y=1.0d0
h=0.1d0     !Variamos o valor
n=((b-a)/h)
x=a
do i=0,n-1
 y=y+h*derivada(x,y)+((h**2)/2)*(-y-(x*derivada(x,y)))
 x=x+h
write(1,*)'Os valores obtidos pelo método da serie de taylor foram x=',x,'e y=',y
write(30,*)x,y
end do
return
end subroutine taylor
