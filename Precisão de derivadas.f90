!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!***********************************************************!!
!! Objetivo: Fazer um programa que investigue a precisão de  !! 
!!          fórmulas aproximadas de derivadas                !!
!!                                                           !!
!! Aluna: Laura Simonassi Raso de Paiva                      !!   
!! Matrícula: 180021885                                      !!
!! Data: 17/09/2019                                          !! 
!!***********************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module vg

real(kind=8) :: primeira, segunda, terceira, quarta !Derivadas exatas da função escolhida
real(kind=8) :: fqpa, fqpb, fqpc                    !Aproximação para quatro pontos
real(kind=8) :: fcpa, fcpb, fcpc, fcpq              !Aproximação para cinco pontos
real(kind=8) :: x,h                                 !Parâmetros da função

end module vg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!************************************************************!!
!! Programa que efetua a verificação das fórmulas aproximadas !!
!!                         de derivadas                       !!
!!************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program verificacao_derivada

call efetuacao_derivada
call impressao_resultados

end program verificacao_derivada

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!************************************************************!!
!!      Subrotina que calcula as derivadas aproximadas        !!
!!************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine efetuacao_derivada
use vg

!Escolha do valor para o parâmetro h
Write(*,*) 'Escolha um valor para h'
read(*,*) h

!A função escolhida é cos(x)
Write(*,*) 'Escolha um valor para x'
read(*,*) x


!Definição da das derivadas exatas
primeira=-sin(x)
segunda=-cos(x)
terceira=sin(x)
quarta=cos(x)

!Obtenção dos resultados por fórmulas de derivadas aproximadas para diferentes pontos

!Para quatro pontos

!Primeira derivada
fqpa=((-2)*cos(x-h)-(3)*cos(x)+(6)*cos(x+h)-cos(x+2*h))/(6*h)
!Segunda derivada
fqpb=(cos(x-h)-(2)*cos(x)+cos(x+h))/(h**2)
!Terceira derivada
fqpc=(-cos(x-h)+3*cos(x)-3*cos(x+h)+cos(x+2*h))/(h**3)

!Para cinco pontos

!Primeira derivada
fcpa=(cos(x-2*h)-(8)*cos(x-h)+(8)*cos(x+h)-cos(x+2*h))/(12*h)
!Segunda derivada
fcpb=(-cos(x-2*h)+(16)*cos(x-h)-(30)*cos(x)+(16)*cos(x+h)-cos(x+2*h))/(12*(h**2))
!Terceira derivada
fcpc=(-cos(x-2*h)+(2)*cos(x-h)-(2)*cos(x+h)+cos(x+2*h))/(2*(h**3))
!Quarta derivada
fcpq=((cos(x-2*h)-(4)*cos(x-h)+(6)*cos(x)-(4)*cos(x+h)+cos(x+2*h)))/(h**4)

return
end subroutine efetuacao_derivada

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!************************************************************!!
!!       Subrotina que imprime os resultados obtidos          !!
!!************************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine impressao_resultados
use vg

write(*,*) 'Verifique seus resultados no novo arquivo gerado'
 open(1,file='resultados.txt',status='new')

write(1,*) 'Derivadas exatas para a funcao cos(x) para o x escolhido=',x
write(1,*) 'df/dx=', primeira
write(1,*)
write(1,*) 'd^2f/dx^2=', segunda
write(1,*)
write(1,*) 'd^3f/dx^3=', terceira
write(1,*)
write(1,*) 'd^4f/dx^4=', quarta
write(1,*)
write(1,*) 'Aproximacao das derivadas com 4 pontos'
write(1,*) 'f**(1)=', fqpa
write(1,*) 'Erro=', primeira-fqpa
write(1,*)
write(1,*) 'f**(2)=', fqpb
write(1,*) 'Erro=', segunda-fqpb
write(1,*)
write(1,*) 'f**(3)=', fqpc
write(1,*) 'Erro=', terceira-fqpc
write(1,*)
write(1,*)
write(1,*) 'Aproximacao das derivadas com 5 pontos'
write(1,*) 'f**(1)=', fcpa
write(1,*) 'Erro=', primeira-fcpa
write(1,*)
write(1,*) 'f**(2)=', fcpb
write(1,*) 'Erro=', segunda-fcpb
write(1,*)
write(1,*) 'f**(3)=', fcpc
write(1,*) 'Erro=', terceira-fcpc
write(1,*)
write(1,*) 'f**(4)=', fcpq
write(1,*) 'Erro=', quarta-fcpq

close(1)
return
end subroutine impressao_resultados
