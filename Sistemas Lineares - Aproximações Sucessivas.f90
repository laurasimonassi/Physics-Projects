!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!**************************************************!!
!! Objetivo: escrever um programa que calcule a     !!
!!   a solução de sistemas lineares pelo método de  !!
!!             aproximação sucessivas               !!
!!                                                  !!
!! Aluna: Laura Simonassi Raso de Paiva             !!
!! Matrícula: 180021885                             !!
!! Data: 22/10/2019                                 !!
!!**************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module vg
real(kind=8):: xk1                    !x1(k+1)- Solução do sistema
real(kind=8):: xk2                    !x2(k+1)- Solução do sistema
real(kind=8):: xkant1,xkant2          !X(k) - solução anterior como passada em sala
end module vg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!**************************************************!!
!!       Programa que resolve o sistema linear      !!
!!**************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program sistemalinear
use vg
implicit none 

!Usamos valores iniciais para fazermos iterações até que o sistema seja resolvido
xkant1 = 0.01d0
xkant2 = 0.01d0

open (unit=1, file='resultados.txt', status='new')
call calculo_solucao
write(1,*)'Como solução para o sistema linear, temos x1=',xk1,'e x2=',xk2
close(1)
end program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!**************************************************!!
!!    Subrotina que calcula a solução do sistema    !!
!!**************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculo_solucao
use vg
implicit none
  
do
  xk1=sqrt((1.0d0/3.0d0)*((7.0d0/2.0d0)-xkant2))
  xk2=((13.0d0/8.0d0)-xkant1)**(1.0d0/3.0d0)
  
  if(abs((xk1-xkant1)/xkant1)<1d-12 .and. abs((xk2-xkant2)/xkant2)<1d-12) exit
   
  xkant1=xk1
  xkant2=xk2

end do
return
end subroutine
