!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!*********************************************************!!
!! Objetivo: Fazer um programa que realize decomposição LU !!
!!                                                         !!
!! Aluna: Laura Simonassi Raso de Paiva                    !!
!! Matrícula: 180021885                                    !!
!! Data: 03/09/2019                                        !!
!!*********************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module vg

integer                  :: N              !Variável que representa a dimensão da matriz
real(kind=8),allocatable :: Mat(:,:)       !Matriz a ser lida
real(kind=8),allocatable :: Mat_Int(:,:)
real(kind=8),allocatable :: L(:,:), U(:,:), A(:,:) !Respectivamente: Matriz L, Matriz U e Matriz A(Resultante da Multiplicação das
!matrizes L e U)

end module vg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!***********************************************!!
!! Programa que efetua o cálculo do determinante !!
!!***********************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program decomposicao

  call leitura_matriz
  call lu
  call escreve_resultados

end program decomposicao
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!****************************************************!!
!! Subrotina que lê a matriz em um arquivo de entrada !!
!!****************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine leitura_matriz
use vg
implicit none
integer :: i,j !Variáveis de auxílio que descrevem elementos da matriz

!No arquivo de entrada devem estar escritas a dimensão da matriz e a matriz, respectivamente.
!Aqui, ao invés de usar a matriz original do exercício 1, realiza-se uma troca de linhas nessa matriz, como foi dito em sala que era
!possível fazer. A matriz utilizada com as linhas já trocadas está descrita no arquivo de entrada.

open(unit=1,file="dados.in",status="old")
   read(unit=1,fmt=*) N
   allocate(Mat(N,N), Mat_Int(N,N),L(N,N),U(N,N), A(N,N))
     do i = 1, N
      read(1,*)(Mat(i,j),j=1,N)
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
integer:: i, j, k

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

A=matmul(L,U)

return
end subroutine lu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!****************************************************!!
!!       Subrotina que testa a leitura da matriz,     !! 
!!            como uma forma de verificação.          !!
!!     Esse teste pode ser realizado antes e após a   !!
!!                 eliminacao de Gauss                !!
!!****************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine escreve_resultados
use vg
implicit none
integer :: i,j

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
write(1,*) "Matriz A que e verifica a decomposicao LU ao multiplicar as matrizes L e U. Deve ser a mesma que a matriz lida."
  do i=1,N
    write(1,10)(A(i,j), j=1,N)
 enddo

close(1)

10 format( <N>(f6.2,1x) )

return
end subroutine escreve_resultados
