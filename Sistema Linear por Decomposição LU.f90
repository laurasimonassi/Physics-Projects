!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!*********************************************************!!
!! Objetivo: Fazer um programa que resolva um sistema      !!
!! linear de equações a partir da decomposição LU          !!
!!                                                         !!
!! Aluna: Laura Simonassi Raso de Paiva                    !!
!! Matrícula: 180021885                                    !!
!! Data: 03/09/2019                                        !!
!!*********************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module vg

integer                  :: N              !Variável que representa a dimensão da matriz
integer :: i,j                             !Variáveis de auxílio que descrevem elementos das matrizes e do vetor
real(kind=8),allocatable :: Mat(:,:)       !Matriz a ser lida
real(kind=8),allocatable :: Mat_Int(:,:)   !Matriz Intermediária
real(kind=8),allocatable :: L(:,:), U(:,:) !Respectivamente: Matriz L, Matriz U e Matriz A(Resultante da Multiplicação das matrizes L e U)
real(kind=8),allocatable :: b(:)
real(kind=8),allocatable :: X(:),y(:)

end module vg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!***********************************************!!
!!   Programa que resolve o sistema linear       !!
!!***********************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program sistemalinear

  call leitura_matriz
  call lu
  call calcula_solucao
  call escreve_resultados

end program sistemalinear

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!****************************************************!!
!! Subrotina que lê a matriz e o vetor em um arquivo  !! 
!!                    de entrada                      !!
!!****************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine leitura_matriz
use vg

!No arquivo de entrada devem estar escritas a dimensão da matriz e a matriz, respectivamente.
!Aqui, ao invés de usar a matriz original do exercício 1, realiza-se uma troca de linhas nessa matriz, como foi dito em sala que era
!possível fazer. A matriz utilizada com as linhas já trocadas está descrita no arquivo de entrada.

open(unit=1,file="dados.in",status="old")
   read(unit=1,fmt=*) N
   allocate(Mat(N,N), Mat_Int(N,N),L(N,N),U(N,N))
   allocate(b(N),x(N),y(N))
 
   do i = 1, N
    read(1,*)(Mat(i,j),j=1,N)
   end do   
 
   do i=1,N
    read(1,*)(b(i))
   end do
close(unit=1)
Mat_Int=Mat
return
end subroutine leitura_matriz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!*********************************************!!
!! Subrotina que realiza a decomposição Lu     !!
!!*********************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine lu
use vg
implicit none
real(kind=8):: coef
integer:: k

!Definição dos elementos de L e U como nulos, para iniciar a decomposição.
L=0.0
U=0.0

do k=1, n-1
 do i=k+1,n
  coef=Mat_Int(i,k)/Mat_Int(k,k)
  L(i,k)=coef
  Mat_Int(i,k)=0
   do j=k+1,n
    Mat_Int(i,j)=Mat_Int(i,j)-coef*Mat_Int(k,j)
   end do
  end do
end do

!Diagonal de L=1
do i=1,n
 L(i,i)=1.0
end do

do j=1,n
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

!Primeiramente, fazemos em que Ly=b e utilizamos a substituição progressiva
!Primeiro termo de y
 y(1)=b(1)/L(1,1)
!Definição vetor y
 do i=2,N
  y(i)=b(i)
  do j=1, i-1
   y(i)=y(i)-L(i,j)*y(j)
  end do
  y(i)=y(i)/L(i,i)
 end do

!Em seguida, usamos a definição Ux=y 
x(N)=y(N)/U(N,N)
 do i=N-1,1,-1!Devido a substituição regressiva
  x(i)=y(i)
  do j=i+1,N
   x(i)=x(i)-U(i,j)*x(j) 
  end do
  x(i)=x(i)/U(i,i)
 end do

end subroutine calcula_solucao
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!****************************************************!!
!!    Subrotina que imprimirá os resultados obtidos   !!
!!       para cada Matriz e para os vetores           !!
!!****************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine escreve_resultados
use vg

 open(1,file='resultados.txt',status='new')
   write(1,*) "Matriz lida"
   do i=1,N
    write(1,10)(Mat(i,j), j=1,N)
   enddo
 
   write(1,*) "Matriz L"
   do i=1,N
     write(1,10)(L(i,j), j=1,N)
   enddo
  
  write(1,*) "Matriz U"
   do i=1,N
    write(1,10)(U(i,j), j=1,N)
  enddo

  write(1,*) "Vetor b lido"
   do i=1,N
    write(1,10)(b(i))
  enddo

  write(1,*) "Vetor y"
   do i=1,N
    write(1,10)(y(i))
  enddo
  
  write(1,*) "Vetor x solucao do sistema linear"
   do i=1,N
    write(1,10)(x(i))
  enddo

close(1)

 10 format( <N>(f6.2,1x) )
return
end subroutine escreve_resultados 

