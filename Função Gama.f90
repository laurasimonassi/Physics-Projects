!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!*****************************************************************!!
!!     Objetivo: fazer um programa para calcular a função gamma    !!
!!                                                                 !!
!!  Aluna: Laura Simonassi Raso de Paiva                           !!
!!  Matrícula: 180021885                                           !!
!!  Data: 19/10/2019                                               !!
!!  Última modificação: 24/10/2019                                 !!
!!*****************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module vg
real(kind=8):: gama
real(kind=8):: x
real(kind=8):: a(0:20)
integer:: n

end module vg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!*****************************************************************!!
!!           Programa que avalia a função gamma no ponto           !!
!!*****************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program valor_gamma 
use vg
implicit none

open(unit=1, file='coef.txt', status='old')

do n=0,20
 read(1,*)a(n)
end do
close(1)

open(unit=10, file='resultados.txt', status='new')
!No arquivo de saída teremos o valor de x utilizado, o valor da função gamma no ponto pelo método aplicado, o valor no ponto pela função gamma do fortran e o erro.

x=3.8d0
call calculo_gamma
write(10,*)x,gama,gamma(x),(gama-gamma(x))/gamma(x)

x=0.6d0
call calculo_gamma
write(10,*)x,gama,gamma(x),(gama-gamma(x))/gamma(x)

x=5.4d0
call calculo_gamma
write(10,*)x,gama,gamma(x),(gama-gamma(x))/gamma(x)

close(10)
end program valor_gamma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!*****************************************************************!!
!!           Subrotina que faz o cálculo da função gamma           !!
!!*****************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculo_gamma
use vg
implicit none
real(kind=8):: b,bn,bn_ant, bsoma !Termos para encontrar os coeficientes b
real(kind=8):: x_s,x_ss    !X auxiliares para não modificarmos o valor de x
real(kind=8):: mult        !Produtório
integer:: i


gama=0.0d0
mult=1.0d0 !Porque se não iremos zerar o produtório

x_s=x

!Verificação que x se enquadra nas condições iniciais necessária

!X maior que 2
 if (x_s > 2.0d0) then
  do while (x_s > 2.0d0)
   mult=mult*(x_s-1.0d0)
   x_s=x_s-1.0d0
  end do
 end if

!X menor que 1
 if (x_s < 1.0d0) then
  do while (x_s < 1.0d0)
   mult=mult/x_s
   x_s=x_s+1.0d0
  end do
 end if

!Calculo do valor de gamma - Expansão em Polinômios de Chebyshev

 do n=0,20 !Cálculo para gamma com os 20 coeficientes a(n)
  
  b=0.0d0        
  bn=0.0d0
  bn_ant=1.0d0

x_ss=x_s-1.0d0 !Para começarmos do B(i-1)

  do i=1,n
   bn=bn_ant*((-n+i-1)*(n+i-1)*(1-x_ss))/((i-0.5d0)*i)
   b=b+bn
   bn_ant=bn
  end do
  bsoma=1.0d0+b      !Pela definição do polinômio Tn*(x)
 
  gama=gama+a(n)*bsoma
end do

!Caso x não esteja no intervalo necessário, correção para o valor de gamma com o produtório
if (x>2.0d0) then
 gama=gama*mult
end if

if (x<1.0d0) then
 gama=gama*mult
end if

return
end subroutine calculo_gamma

