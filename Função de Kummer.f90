!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!***********************************************************!!
!!   Objetivo: Fazer um programa que calcule os valores da   !!
!!                       função de Kummer                    !!
!!  Aluna: Laura Simonassi Raso de Paiva                     !!
!!  Matrícula: 180021885                                     !!
!!  Data:04/10/2019                                          !!
!!  Ultima modificação: 15/10/2019                           !!
!!***********************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module globais

real(kind=8):: a,b,x,p
real(kind=8):: res, res_in, an, an_um
integer:: n

end module globais
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!***********************************************************!!
!!    Programa que calcula os valores da função de Kummer    !!
!!                   e seus respectivos erros                !!  
!!***********************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program kummer
use globais

 open(unit=1,file='resultados.txt',status='new')

a=1.0d0
b=0.5d0
x=0.3d0
call coef
  write(1,*) 'Para M(',a,',',b,',',x,'), temos o valor', res, 'com erro de', res-1.7357213


a=0.5d0
b=1.0d0
x=5.0d0
call coef
  write(1,*) 'Para M(',a,',',b,',',x,'), temos o valor', res, 'com erro de', res-40.078446

a=1.0d0
b=1.0d0
x=7.0d0
call coef
  write(1,*) 'Para M(',a,',',b,',',x,'), temos o valor', res, 'com erro de', res-1096.6332

 close(1)

end program kummer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!*****************************************************************!!
!! Subrotina que calcula os coeficientes An e efetua seu somatório !!
!!         para obtemos o resultado da função de Kummer            !!
!!*****************************************************************!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine coef
use globais

p=1.0e-10                 !Precisão desejada que deve ser atingida para que haja convergencia no somatório
an_um=1.0d0               !Definição de A0=1
res_in=1.0d0              
n=0.0d0
an=0.0d0

res=res_in

  do 
  n=dfloat(n)+1.0d0
  an = ((a+dfloat(n)-1.0d0)/(b+dfloat(n)-1.0d0))*(x/dfloat(n))*an_um
  res=res+an
  if (abs(an/res)<=p) exit
  an_um=an
  end do

return
end subroutine coef

