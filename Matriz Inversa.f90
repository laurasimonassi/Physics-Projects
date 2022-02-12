!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!*********************************************************!!
!! Objetivo: Fazer um programa que calcule a inversa de    !!
!!                    uma matriz NxN                       !!
!! Aluna: Laura Simonassi Raso de Paiva                    !!
!! Matrícula: 180021885                                    !!
!! Data: 03/09/2019                                        !!
!!*********************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
module vg
        
integer:: N !Variável que representa a dimensão da matriz
integer :: i,j !Variáveis de auxílio que descrevem elementos da matriz
real(kind=8),allocatable :: A(:,:) !Matriz lida no arquivo de entrada
real(kind=8),allocatable :: B(:,:) !Matriz Intermediára
real(kind=8),allocatable :: C(:,:) !Matriz inversa da matriz lida
real(kind=8),allocatable :: D(:,:) !Matriz Identidade
        
      end module vg
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!***********************************************!!
!! Programa que efetua a inversão da matriz      !!
!!***********************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      program matrizinversa
      
      call leitura_matriz
      call calcula_inversa
      call escreve_resultados

      end program matrizinversa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!****************************************************!!
!! Subrotina que lê a matriz em um arquivo de entrada !!
!!****************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine leitura_matriz
use vg

!No arquivo de entrada devem estar escritas a dimensão da matriz e a matriz, respectivamente.
 open(unit=1,file="dados.in",status="old")
   read(unit=1,fmt=*) N
   allocate(A(N,N))
   allocate(B(N,N))
   allocate(C(N,N))
   allocate(D(N,N))
     do i = 1, N
       read(1,*)(A(i,j),j=1,N)
     end do
 close(unit=1)

   B=A !Escrevemos a matriz intermediária igual a matriz lida, para que ela auxilie na formação da matriz inversa sem ser destruída.

   return
end subroutine leitura_matriz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!****************************************************!!
!! Subrotina que efetua o cálculo da matriz inversa   !!
!!****************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calcula_inversa
use vg
implicit none
integer:: linha, k, m
real(kind=8):: pivot, l

!Definimos inicialmente uma matriz identidade, a ser transformada:
 do i=1,N
   do j=1,N
    C(i,j)=0.0
   end do
   C(i,i)=1.0
 end do

!Definindo que o melhor pivot. 
     
 do i=1,N
   pivot=B(i,i)
   do j=1,N
   if (B(j,i)>pivot) then
    pivot=B(j,i)
    linha=j
   end if
 end do
      
!Efetuação das trocas das linhas e alocação na matriz que será a inversa.
 if (pivot>B(i,i)) then
   do k=1,N
      l=B(i,k)
      B(i,k)=B(linha,k)
      B(linha,k)=l
        
      l=C(i,k)
      C(i,k)=C(linha,k)
      C(linha,k)=l
    end do
 end if
      
!Método de Gauss-Jordan
!Divisão pelo variável de auxílio para obter as matrizes que, multiplicadas, resultarão na matriz identidade.    

  l=B(i,i)
  do j=1,N
   B(i,j)=B(i,j)/l
   C(i,j)=C(i,j)/l
  end do

  do j=i+1,N
  l=B(j,i)
   do k=1,N
    B(j,k)=B(j,k)-l*B(i,k)
    C(j,k)=C(j,k)-l*C(i,k)
   end do
  end do
end do

  do i=1, N-1
   do j=i+1,N
    l=B(i,j)
     do m=1,N
      B(i,m)=B(i,m)-l*B(j,m)
      C(i,m)=C(i,m)-l*C(j,m)
     end do
    end do
  end do

   D=matmul(C,A) !Matriz Identidade resultará da matriz A lida multiplicada por sua inversa.
return
end subroutine calcula_inversa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!****************************************************!!
!!    Subrotina que imprimirá os resultados obtidos   !!
!!                  para cada Matriz                  !!
!!****************************************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine escreve_resultados
use vg

 open(1,file='resultados.txt',status='new')
   write(1,*) "Matriz lida"
   do i=1,N
    write(1,10)(A(i,j), j=1,N)
   enddo

   write(1,*) "Matriz Inversa"
    do i=1,N
     write(1,10)(C(i,j), j=1,N)
    enddo
     
   write(1,*) "Matriz Identidade, obtida atraves da multiplicacao de A e de sua inversa C"
     do i=1,N
      write(1,10)(D(i,j), j=1,N)
     enddo

close(1)

 10 format( <N>(f6.2,1x) )
return
end subroutine escreve_resultados 























