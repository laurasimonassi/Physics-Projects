!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!**************************************************************************!!
!!   Objetivo: realizar um programa que calcule numericamente integrais     !!
!!  utilizando o método de retângulo, trapézio, regra de Simpson e de Gauss !!
!!                                                                          !!
!! ALuna: Laura Simonassi Raso de Paiva                                     !!
!! Matrícula: 180021885                                                     !!
!! Data: 18/09/2019                                                         !!
!! Utilma modificação: 03/10/2019                                           !!
!!**************************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module vg

real(kind=8):: a,b     !Limites inferior e superior da integração

end module vg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!**************************************************************************!!
!!        Programa que efetua a integração por todos os métodos             !!
!!**************************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program integral

use vg
implicit none
real(kind=8):: f                                                                                  !Função a ser integrada
real(kind=8):: gauss_8, resultado_gauss8,gauss_12, resultado_gauss12,gauss_20, resultado_gauss20  !Definição de valores obtidos pela quadratura de gauss
real, parameter:: pi=3.14159265358979d0                                                           !Valor de pi para comparação 
integer :: i,j,k  
external f

a=-1.0d0
b=1.0d0

open(unit=1, file='resultados.txt', status='new')
write(1,*) 'Regra do Retangulo'
do i=500,15000,500
call retangulo(i)
end do
write(1,*) 'Regra do Trapezio'
do j=500,15000,500
call trapezio(j)
end do
write(1,*) 'Regra de Simpson'
do k=500,15000,500
call simpson(k)
end do

!Impressão dos resultados do método de quadratura de Gauss
write(1,*)'Para quadratura de Gauss com n=8, obtemos o resultado para a integral da funcao como='
resultado_gauss8=gauss_8(f,a,b)
write(1,*) resultado_gauss8, 'e o erro e', resultado_gauss8-(pi/2.0d0)
write(1,*)
write(1,*)'Para quadratura de Gauss com n=12, obtemos o resultado para a integral da funcao como='
resultado_gauss12=gauss_12(f,a,b)
write(1,*) resultado_gauss12, 'e o erro e', resultado_gauss12-(pi/2.0d0)
write(1,*)
write(1,*)'Para quadratura de Gauss com n=20, obtemos o resultado para a integral da funcao como='
resultado_gauss20=gauss_20(f,a,b)
write(1,*) resultado_gauss20, 'e o erro e', resultado_gauss20-(pi/2.0d0)


close(1)
end program integral

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!**************************************************************************!!
!!                            Função a ser integrada                        !!
!!**************************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Function f(x)
real(kind=8) :: f,x

f = sqrt(1.0d0-(x**2)) !Função utilizada  para o exercício
return
end 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!**************************************************************************!!
!!          Subrotina que integra usando a regra do retângulo               !!
!!**************************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine retangulo(n)
use vg
implicit none
real(kind=8) :: resultado_ret                                          !Resultado da integral
real(kind=8) :: f,x                                                    !Função a ser integrada e valor em x
real(kind=8) :: h                                                      !Variação em x
real, parameter:: pi=3.14159265358979d0                                !Valor de pi para comparação 
integer :: i                                                           !Variável de auxilio
integer:: n                                                            !Número de subintervalos


h=(b-a)/dfloat(n)
resultado_ret=0.0d0
 do i=0, n-1
  x=a+dfloat(i)*h
  resultado_ret = resultado_ret + h*f((x+x+h)/2.0d0)
 end do
    
!Impressão de resultados

write(1,*)'A integracao da funcao usando a regra do retangulo para',n,'intervalos resulta em=', resultado_ret, 'e o erro e', resultado_ret-(pi/2.0d0)
write(1,*)
return
end subroutine retangulo      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!**************************************************************************!!
!!           Subrotina que integra usando a regra do trapézio               !!
!!**************************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine trapezio(n)
use vg

implicit none
real(kind=8):: resultado_trap                     !Resultado da integral
real(kind=8):: f,x                                !Função a ser integrada e valor em x
real(kind=8):: h                                  !Variação em x
real, parameter:: pi=3.14159265358979d0           !Valor de pi para comparação 
integer :: i                                      !Variável de auxilio
integer:: n                                       !Número de subintervalos
    



   h=(b-a)/dfloat(n)
   resultado_trap=0.0d0   
   do i=1, n-1
     x=a+dfloat(i)*h
     resultado_trap = resultado_trap + f(x)
   end do
   
   resultado_trap = h*(resultado_trap +(f(a)+f(b))/2.0d0)
  
!Impressão de resultados

write(1,*)'A integracao da funcao usando a regra do trapezio para',n,'intervalos resulta em=', resultado_trap, 'e o erro e', resultado_trap-(pi/2.0d0)
write(1,*)
return
end subroutine trapezio
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!**************************************************************************!!
!!            Subrotina que integra usando a regra de Simpson               !!
!!**************************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine simpson(n)
use vg
implicit none
real(kind=8):: resultado_simp                     !Resultado da integral
real(kind=8):: f,x                                !Função a ser integrada e valor em x
real(kind=8) :: h                                 !Variação em x
real, parameter:: pi=3.14159265358979d0           !Valor de pi para comparação 
integer :: i                                      !Variável de auxilio
integer:: n                                       !Número de subintervalos

if ((n/2)*2 .ne. n) then
     n=n+1
end if

   h=(b-a)/dfloat(n)

   resultado_simp=0.0d0
   do i=0, n-4,2
     x=a+dfloat(i+1)*h
     resultado_simp = resultado_simp + 4.0d0*f(x) + 2.0d0*f(x+h)
   end do
     resultado_simp= (h*(resultado_simp + f(a) + f(b) + 4.0d0*f(b-h))/3.0d0)

!Impressão de resultados

write(1,*)'A integracao da funcao usando a regra de simpson para',n,'intervalos resulta em=', resultado_simp, 'e o erro e', resultado_simp-(pi/2.0d0)
write(1,*)
return
end subroutine simpson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!**************************************************************************!!
!!                       Integração por função Gauss                        !!
!!**************************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!
!!***********!!
!! Gauss n=8 !!
!!***********!!
!!!!!!!!!!!!!!!
function gauss_8 (f,a,b)
implicit none
integer,parameter:: n=4                                                 !n dos pontos da quadratura de Gauss/2
integer:: i                                                             !Variável de auxílio
real(kind=8)::gauss_8, a, b, f                                          !Definição das funções e limites
real(kind=8):: x(n), w(n)
real(kind=8):: res, pos, neg                                            !Variáveis de auxílio
DATA x/0.183434642495650d0, 0.525532409916329d0, 0.796666477413627d0, 0.960289856497536d0/          !Valores obtidos pela tabela
DATA w/0.362683783378362d0, 0.313706645877887d0, 0.222381034453374d0, 0.101228536290376d0/ 


res=0.0d0
pos=(b+a)/2.0d0
neg=(b-a)/2.0d0

do i=1,n
 res=res+(w(i)*(f(neg*(-x(i))+pos)+f(neg*x(i)+pos)))
end do
gauss_8=res*neg


return
end function
!!!!!!!!!!!!!!!!
!!************!!
!! Gauss n=12 !!
!!************!!
!!!!!!!!!!!!!!!!
function gauss_12 (f,a,b)
implicit none
integer,parameter:: n=6                                                  !n dos pontos da quadratura de Gauss/2
integer:: i                                                              !Variáveis de auxílio
real(kind=8)::gauss_12, a, b, f                                          !Definição das funções e limites
real(kind=8):: x(n) , w(n)
real(kind=8):: res, pos, neg                                             !Variáveis de auxílio
DATA x/0.125233408511469d0, 0.367831498998180d0, 0.587317954286617d0, 0.769902674194305d0, &         !Valores obtidos pela tabela
    0.904117256370475d0, 0.981560634246719d0/
  DATA w/0.249147045813403d0, 0.233492536538355d0, 0.203167426723066d0, 0.160078328543346d0, &
    0.106939325995318d0, 0.047175336386512d0/

res=0.0d0
pos=(b+a)/2.0d0
neg=(b-a)/2.0d0

do i=1,n
 res=res+(w(i)*(f(neg*(-x(i))+pos)+f(neg*x(i)+pos)))
end do
gauss_12=res*neg


return
end function
!!!!!!!!!!!!!!!!
!!************!!
!! Gauss n=20 !!
!!************!!
!!!!!!!!!!!!!!!!
function gauss_20 (f,a,b)
implicit none
integer,parameter:: n=10                                           !n dos pontos da quadratura de Gauss/2
integer:: i                                                        !Variáveis de auxílio
real(kind=8)::gauss_20, a, b, f                                     !Definição das funções e limites
real(kind=8):: x(n) , w(n)
real(kind=8):: res, pos, neg                                        !Variáveis de auxílio
DATA x/0.076526521133497333755d0, 0.227785851141645078080d0, 0.373706088715419560673d0, 0.510867001950827098004d0, &    !Valores obtidos pela tabela
    0.636053680726515025453d0, 0.746331906460150792614d0, 0.839116971822218823395d0, 0.912234428251325905868d0, &
    0.963971927277913791268d0, 0.993128599185094924786d0/
  DATA w/0.152753387130725850698d0, 0.149172986472603746788d0, 0.142096109318382051329d0, 0.131688638449176626898d0, &
    0.118194531961518417312d0, 0.101930119817240435037d0, 0.083276741576704748725d0, 0.062672048334109063570, &
    0.040601429800386941331d0, 0.017614007139152119312d0/ 

res=0.0d0
pos=(b+a)/2.0d0
neg=(b-a)/2.0d0

do i=1,n
 res=res+w(i)*(f(neg*(-x(i))+pos)+f(neg*x(i)+pos))
end do
gauss_20=res*neg


return
end function


















