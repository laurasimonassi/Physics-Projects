!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!******************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Objetivo: Calcular o determinante de uma matriz NxN !!
!!            pelo método de eliminação de Gauss.       !!
!!                                                      !!
!!     Aluna: Laura Simonassi Raso de Paiva             !!
!!     Matrícula: 180021885                             !!
!!     Data: 20/08/2019                                 !! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!******************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module vg

integer                  :: N         !Variável que representa a dimensão da matriz
real(kind=8),allocatable :: Mat(:,:)  !Matriz a ser lida
real(kind=8)             :: det       !Determinante
real(kind=8)             :: t         !Variável que armazena o número de trocas entre linhas
end module vg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!***********************************************!!
!! Programa que efetua o cálculo do determinante !!
!!***********************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program determinante

  call leitura_matriz
  call el_gauss
  call escreve_resultados

end program determinante
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
open(unit=1,file="dados.in",status="old")
   read(unit=1,fmt=*) N
   allocate(Mat(N,N))
   do i = 1, N
      read(1,*)(Mat(i,j),j=1,N)
   end do
close(unit=1)

return
end subroutine leitura_matriz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!*********************************************!!
!! Subrotina que realiza a eliminação de Gauss !!
!!*********************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine el_gauss
use vg
implicit none
integer:: i, j, k, l            !Variáveis de auxílo
real(kind=8) :: gauss, pivot    !Variáveis que auxiliam na eliminação de gauss
                                ! e descrevem o pivot, respectivamente


do k=1, n-1

 !Identificação da fileira com maior elemento para realizar a troca do pivot

 pivot=Abs(Mat(k,k))
 l=k

 do j=k+1,n
  if (Abs(Mat(j,k) > pivot)) then 
    pivot = Abs(Mat(j,k))
    l=j
  end if
 end do
  
  
 !Troca de fileiras l e k, caso necessário
 
 if (l /= k) then
   do j=k,n
     t= Mat(k,j)
     Mat(k,j)=Mat(l,j)
     Mat(l,j)=t
   end do
 end if

 !Verificação da singularidade da matriz ao analisar elementos nulos
    if(pivot == 0.0) then
      write(*,*) ' A matriz é singular '
    return
    end if

 !Eliminação e substituição 
 
 do i=k+1,n
  gauss=Mat(i,k)/Mat(k,k)
  Mat(i,k)=0.0
  do j=k+1, n
   Mat(i,j)=Mat(i,j) - gauss*Mat(k,j)
   end do
 end do
end do
 
 !Com a matriz triangularizada, retornamos para o programa principal que efetua o calculo do determinante.
 
return
end subroutine el_gauss
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!**************************************************!!
!! Subrotina que escreve os resultados em uma saída !!
!!**************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine escreve_resultados
use vg
implicit none
integer :: i,j

open(unit=10, file='resultadodet.txt',status='unknown')
  !Cálculo do determinante caso a matriz não seja singular 
  
  det=1
  
 write(10,*) 'A matriz apos a eliminacao de Gauss e dada por:'
  do i=1,N
   write(10,10)(Mat(i,j), j=1,N)
  
  10 format( <N>(f4.1,1x) )
 end do
  
  do i=1,n
   det=Mat(i,i)*det
   det=det*((-1)**t)
  end do
 
 write(10,*) 'O determinante desta matriz é dado por ', det, '.'
  
  close(10)

return
end subroutine escreve_resultados
