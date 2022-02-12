!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!******************************************************************************!!
!! Objetivo: escrever um programa que calcule os valores da função geométrica   !!
!! F(a,b,c,x) e testá-lo para ver se satisfaz algumas relações para |x|<1       !!
!!                                                                              !!
!! Aluna: Laura Simonassi Raso de Paiva                                         !!
!! Matrícula:180021885                                                          !!
!! Data: 17/10/2019                                                             !!
!!******************************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!******************************************************************************!!
!!   Programa que calcula o valor das funções hipergeomérticas. Ele terá como   !!
!! resultado o valor calculado pela função e os valores para as funções a serem !!
!!   comparadas, de forma que possamos avaliar se as relações são satisfeitas   !!
!!******************************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program geometrica

implicit none
real(kind=8):: a,b,c,x,y,res
integer:: n,i

open(unit=1, file='resultadosln.txt', status='new')
 
 a=1.0d0
 b=1.0d0
 c=2.0d0
 
do x=-0.9d0,0.9d0,0.01d0                                                       !Variar o valor de x 
  call hiper(a,b,c,x,res)
  write(1,*) x,res,(-1.0d0/x)*log(1.0d0-x), res-(-1.0d0/x)*log(1.0d0-x)  !Valor de x, resultado da função para esse x, a função a ser comparada e o erro
end do 

close(1)

open(unit=10, file='resultadosarctan.txt', status='new')

 a=0.5d0
 b=1.0d0
 c=1.5d0
 
do x=-0.9d0,0.9d0,0.01d0                                                        !Variar o valor de x 
  y=-(x**2)
  call hiper(a,b,c,y,res)                                 
  write(10,*) x,res,(1.0d0/x)*atan(x), res-(1.0d0/x)*atan(x)              !Valor de x, resultado da função para esse x, a função a ser comparada e o erro
end do 

close(10)

open(unit=11, file='resultadosa.txt', status='new')

do i=1,19

a=1.0d0
b=2.0d0
 do x=-0.9d0,0.9d0,0.1d0                                                        !Variar o valor de x 
  call hiper(a,b,b,x,res)                                 
  write(11,*) a,b,x,res,(1.0d0/((1.0d0-x)**a)), res-(1.0d0/((1.0d0-x)**a)) !Valor de x, resultado da função para esse x, a função a ser comparada e o erro

a=a+1.0d0
b=b+1.0d0
end do
end do 

close(11)

end program geometrica

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!*****************************************************************!!
!!          Subrotina que calcula a função hipergeométrica         !!
!!*****************************************************************!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine hiper(a,b,c,x,res)

implicit none
real(kind=8):: a,b,c,x,p
real(kind=8):: res, res_in, an, an_um
integer:: n

p=1.0e-10                 !Precisão desejada que deve ser atingida para que haja convergencia no somatório
an_um=1.0d0               !Definição de A0=1
res_in=1.0d0              
n=0.0d0
an=0.0d0

res=res_in

  do 
   n=dfloat(n)+1.0d0
   an = an_um*(((a+dfloat(n)-1.0d0)*(b+dfloat(n)-1.0d0))/(c+n-1.0d0))*(x/dfloat(n))
   res=res+an
    if (abs(an/res)<=p) exit
   an_um=an
  end do
 
return
end subroutine hiper

