!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!************************************************************!!
!! Objetivo: escrever um programa que calcule os valores da   !!
!!         funcao de Bessel J0,J1,Y0 e Y1 para x=5            !!  
!! Aluna: Laura Simonassi Raso de Paiva  		      !!
!! Matrícula: 180021885 				      !!
!! Data: 03/10/2019                                           !!
!! Ultima modificação: 18/10/2019                             !!
!!************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program funcao_bessel
implicit none
real(kind=8):: x,j,a
real(kind=8):: y
real(kind=8):: e1,e2
real, parameter:: valor_jz= -0.1775967713d0
real, parameter:: valor_ju= -0.3275791376d0
real, parameter:: valor_yz= -0.3085176252d0
real, parameter:: valor_yu=  0.1478631434d0
external j
external y

open(unit=1,file='resultados.txt',status='new')

write(1,*)'O resultado obtido ao aplicar x=5 na funcao de bessel J0 foi=',j(0.d0,5.0d0)
write(1,*)'com erro de', j(0.d0,5.0d0)-valor_jz
write(1,*)'O resultado obtido ao aplicar x=5 na funcao de bessel J1 foi=',j(1.d0,5.0d0)
write(1,*)'com erro de', j(1.d0,5.0d0)-valor_ju

e1=0.01d0
do 
if ((abs(y(0.0d0,e1,5.0d0))-abs(y(0.0d0,e1*10.0d0,5.0d0))/abs(y(0.0d0,e1*10.0d0,5.0d0)))<=1d-10) exit
e1=e1/10.0d0
end do
write(1,*)'O resultado obtido ao aplicar x=5 na funcao de bessel Y0 foi=',y(0.0d0,e1,5.0d0)
write(1,*)'com erro de', y(0.0d0,e1,5.0d0)-valor_yz

e2=0.01d0
do
if ((abs(y(1.0d0,e2,5.0d0))-abs(y(1.0d0,e2*10.0d0,5.0d0))/abs(y(1.0d0,e2*10.0d0,5.0d0)))<=1d-10) exit
e2=e2/10.0d0
end do
write(1,*) 'O resultado obtido ao aplicar x=5 na funcao de bessel Y1 foi=',y(1.0d0,e2,5.0d0)
write(1,*)'com erro de', y(1.0d0,e2,5.0d0)-valor_yu


end program funcao_bessel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!************************************************************!!
!!           Função que efetua o calculo para J(5)            !!
!!************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Function j(v,x)

implicit none
real(kind=8):: x               !Valor da aplicação
real(kind=8):: v               !Parametro da função
real(kind=8):: s              
real(kind=8):: soma            !Somatório da função Bessel J0
real(kind=8):: j               !Resultado
integer:: k                    !Variável de auxílio


soma=0.0d0 !Valor inicial do somatório

s=1d0/GAMMA(v+1.0d0)
k=1

do
 soma = ((-(x**2.0d0)/4.0d0)**k)/(GAMMA(k+1.0d0)*GAMMA(v+k+1.0d0))
 if(abs(soma/s)<1d-12)exit
 s=s+soma
 k=k+1
end do
 j=((x/2.0d0)**v)*s


end function j
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!************************************************************!!
!!           Função que efetua o calculo para Y(5)            !!
!!************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Function y(n,e,x)

implicit none
real(kind=8):: x               !Valor da aplicação
real(kind=8):: v,n,e           !Parametro da função
real(kind=8):: y, termo
real(kind=8):: j               !Resultado
integer:: k                    !Variável de auxílio
real(kind=8), parameter:: pi=3.1415926535897932d0
external j
termo=0.0d0

 termo= ((j(n+e,x)*cos((n+e)*pi)-j(-n-e,x))/(sin((n+e)*pi)) + (j(n-e,x)*cos((n-e)*pi)-j(-n+e,x))/(sin((n-e)*pi)))
 y=(0.5)*termo

end function y
