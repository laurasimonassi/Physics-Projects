!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!*******************************************************!!
!! Objetivo: escrever um programa que resolva um sistema !!
!!      não linear pelo método de Newton-Raphson         !!        
!! Aluna: Laura Simonassi Raso de Paiva                  !!
!! Matrícula: 180021885                                  !!
!! Data: 28/10/2019                                      !!
!! Ultima modificação: 31/10/2019                        !!
!!*******************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module vg
real(kind=8), dimension(3,3) :: Mat_Int    
real(kind=8), dimension(3,3) :: L, U       
real(kind=8), dimension(3)   :: b,y
real(kind=8), dimension(3)   :: s                     !Delta
real(kind=8), dimension(3,3) :: F_linha
real(kind=8), dimension(3)   :: F_negativa
real(kind=8)                 :: x1,x2,x3
end module vg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!*******************************************************!!
!!      Programa que resolve sistemas não lineares       !!
!!*******************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Temos que o sistema linear a ser resolvido é:
!x1²+x2²+x3²=9
!x1x2x3=1
!x1+x2-x3²=0

program sistema
use vg
implicit none
real(kind=8):: x1_ant,x2_ant,x3_ant

!Chute inicial para as raízes
x1=1.0d0
x2=2.0d0
x3=3.0d0

call derivadaf
call negativaf

!Nomeação para decomposição LU
Mat_Int=F_linha
b=F_negativa

do 
  x1_ant=x1
  x2_ant=x2
  x3_ant=x3
 
  call lu
  call calcula_solucao
   x1=x1+s(1)
   x2=x2+s(2)
   x3=x3+s(3)
 
  if ((abs((x1-x1_ant)/x1_ant)<=1d-10) .and.(abs((x2-x2_ant)/x2_ant)<=1d-10).and. (abs((x3-x3_ant)/x3_ant)<=1d-10))  exit

   call derivadaf
   call negativaf
   Mat_Int=F_linha
   b=F_negativa
end do

open(unit=1, file='resultados.txt', status='new')
write(1,*)'A resolucao do sistema nao linear resulta nas raizes x1=',x1,'x2=',x2,'x3=',x3
close(1)

end program sistema

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!*******************************************************!!
!!        Definição das funções por subrotinas           !!
!!*******************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine derivadaf
use vg
implicit none
F_linha(1,1)=2.0d0*x1
F_linha(1,2)=2.0d0*x2
F_linha(1,3)=2.0d0*x3
F_linha(2,1)=x2*x3
F_linha(2,2)=x1*x3
F_linha(2,3)=x1*x2
F_linha(3,1)=1.0d0
F_linha(3,2)=1.0d0
F_linha(3,3)=-2.0d0*x3
return
end subroutine derivadaf

subroutine negativaf
use vg
implicit none
F_negativa(1)=-(x1**2.0d0)-(x2**2.0d0)-(x3**2.0d0) + 9.0d0
F_negativa(2)=-(x1*x2*x3) + 1.0d0
F_negativa(3)=-x1-x2+(x3**2.0d0)
return
end subroutine negativaf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!******************************************************************************************!!
!!  Vamos usar o exercício 3 para resolver o sistema linear, mas agora para dimensão 3      !!
!!******************************************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!*********************************************!!
!! Subrotina que realiza a decomposição Lu     !!
!!*********************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine lu
use vg
implicit none
real(kind=8):: coef
integer :: i,j,k

!Definição dos elementos de L e U como nulos, para iniciar a decomposição.
L=0.0
U=0.0

do k=1,2
 do i=k+1,3
  coef=Mat_Int(i,k)/Mat_Int(k,k)
  L(i,k)=coef
  Mat_Int(i,k)=0
   do j=k+1,3
    Mat_Int(i,j)=Mat_Int(i,j)-coef*Mat_Int(k,j)
   end do
  end do
end do

!Diagonal de L=1
do i=1,3
 L(i,i)=1.0
end do

do j=1,3
 do i=1,j
  U(i,j)=Mat_Int(i,j)
 end do
end do

return
end subroutine lu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!*****************************************************!!
!!   Subrotina que realiza os cálculos na resolução    !!
!!                  do sistema linear                  !!
!!*****************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calcula_solucao
use vg
implicit none
integer :: i,j

!Primeiramente, fazemos em que Ly=b e utilizamos a substituição progressiva
!Primeiro termo de y
 y(1)=b(1)/L(1,1)
!Definição vetor y
 do i=2,3
  y(i)=b(i)
  do j=1, i-1
   y(i)=y(i)-L(i,j)*y(j)
  end do
  y(i)=y(i)/L(i,i)
 end do

!Em seguida, usamos a definição Us=y 
s(3)=y(3)/U(3,3)
 do i=2,1,-1!Devido a substituição regressiva
  s(i)=y(i)
  do j=i+1,3
   s(i)=s(i)-U(i,j)*s(j) 
  end do
  s(i)=s(i)/U(i,i)
 end do
return
end subroutine calcula_solucao

