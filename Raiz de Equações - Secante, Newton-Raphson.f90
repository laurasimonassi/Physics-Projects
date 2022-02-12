!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!************************************************************!!
!! Objetivo: escrever um programa que avalie as raízes de uma !!
!! equação pelos métodos da Secante de de Newton-Raphson      !!
!!                                                            !!  
!! Aluna: Laura Simonassi Raso de Paiva  		      !!
!! Matrícula: 180021885 				      !!
!! Data: 19/09/2019					      !!
!! Ultima modificação: 03/10/2019                             !!
!!************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!************************************************************!!
!! Programa que efetua a avaliação das raízes obtidas pelos   !! 
!!                          dois métodos                      !!
!!************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program raizes
implicit none
integer:: g,k
external f
external derivada 

open(unit=1,file="resultados.txt",status='new')
write(1,*)'Pelo metodo da secante, obtemos os resultados:'
call secante(g)
write(1,*)'Pelo metodo Newton Raphson, obtemos:'
call newtonraphson(k)


close(1)
end program raizes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!************************************************************!!
!!                        Função                              !!
!!************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Function f(x)
implicit none
real(kind=8):: f,x
f= x**(4.0d0)-(26.0d0)*x**(3.0d0)+(131.0d0)*x**(2.0d0)-(226.0d0)*x+120.0d0
end function f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!************************************************************!!
!!                    Função Derivada                         !!
!!************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Function derivada(x)
implicit none
real(kind=8):: derivada,x
derivada= 4.0d0*x**(3.0d0)-(78.0d0)*x**(2.0d0)+(262.0d0)*x-226.0d0
end function derivada


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!************************************************************!!
!! Subrotina que calcula raízes através do método da secante  !!
!!************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine secante(raiz)
implicit none
real(kind=8):: f,raiz, p      !Função usada, raiz desta função e precisão
real(kind=8):: xa,xb,v        !Variáveis de auxílio
real(kind=8):: df             !Derivada da função f
integer:: i,j                 !Variável de auxílio
integer:: n                   !Número de interações
p = 1.0e-12                   !Precisão desejada

do j=1,4
write(*,*)'Escolha um valor para uma das raizes'
read(*,*) xa
write(*,*)'Escolha um segundo valor possível para esta mesma raiz'
read(*,*) xb


n=0
write(1,*)'Os chutes iniciais para a raiz sao', xa,'e',xb

!Investigação da raiz
 do      
  n=n+1 
  df=(f(xb)-f(xa))/(xb-xa)
  raiz=xb-(f(xb)/df)
 
!Condição para parar a investigação
  v=raiz-xb
  if (abs(v/xb)<=p) exit
  xa=xb
  xb=raiz


 end do


write(1,*) 'A raiz encontrada para essa funcao tem o valor',raiz,', que e o ponto no qual a funcao assume valor', f(raiz), 'apos',n,'iteracoes'
end do
return
end subroutine secante

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!*******************************************************************!!
!! Subrotina que calcula raízes através do método de newton raphson  !!
!!*******************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine newtonraphson(raiz)
implicit none
real(kind=8):: f, raiz, p, derivada             !Função usada, raiz desta função, precisão e derivada da função
real(kind=8):: f_val, f_der,x0                  !Variáveis de auxílio
real(kind=8):: v                                !Variável de auxílio para verificar a precisão
integer:: i,j                                   !Variável de auxílio
integer:: n                                     !Número de interações
p = 10e-12                                      !Precisão desejada

do j=1,4                                        !Cálculo realizado 4 vezes para obter as quatro raizes
write(*,*)'Escolha um valor possível para a raiz'
read(*,*) x0                                    !"Chute" inicial

n=0
write(1,*)'O chute inicial para a raiz é de', x0

!Investigação da raiz
 do                                             !Quantas iterações forem precisas
  n=n+1
  f_val=f(x0)
  f_der=derivada(x0)

 raiz= x0-(f_val/f_der)                         !Cálculo da raiz
 v=raiz-x0 

 if (abs(v/x0)<=p) exit                         !Verificação da precisão desejada

  x0=raiz
 end do

write(1,*) 'A raiz encontrada para essa funcao tem o valor',raiz,' que e o ponto no qual a funcao assume o valor', f(raiz), 'apos',n, 'iteracoes'
end do

return
end subroutine newtonraphson
