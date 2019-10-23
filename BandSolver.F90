program BandSolver

implicit none

integer, parameter :: n = 4

real(8) A(4, n), B(n), At(n, 3), S(n)
integer ipiv(n)

integer info

integer ku, kl, j, i

real(8), allocatable :: C(:, :), D(:)



write(*, '(a, 4f5.0)') 'B = ', B

! Sistema Linear
!     |  2 -1 0   0 | | x1 |    | 1 |
!     | -1  2 -1  0 | | x2 |    | 0 |
!     |  0 -1  2 -1 | | x3 | =  | 0 |
!     |  0  0 -1  2 | | x4 |    | 1 |

A = 0

A(3, :) = 2
A(2, 2:) = -1
A(4, :n-1) = -1

B = 0
B(1) = 1
B(n) = 1

S = 1

write(*,*) "      A"
write(*, '(4f5.0)') A(1,:)
write(*, '(4f5.0)') A(2,:)
write(*, '(4f5.0)') A(3,:)
write(*, '(4f5.0)') A(4,:)

call dgbsv( n, 1, 1, 1, A, 4, ipiv, B, n, info )

write(*, '(a, i0)') ' info = ', info
write(*, '(a, 4f5.0)') 'B = ', B


if(norm2(B - S) < 1.0d-15) then
    write(*, *) "Correct Solution!"
else
    write(*,*) "Wrong Solution !"
end if


! Sistema Linear
!     |  1  1  0  0 | | x1 |    | 2 |
!     |  1  1  1  0 | | x2 |    | 3 |
!     |  0  1  1  1 | | x3 | =  | 3 |
!     |  0  0  1  1 | | x4 |    | 2 |

A = 0

A(3, :) = 1
A(2, 2:) = 1
A(4, :n-1) = 1

write(*,*) "      A"
write(*, '(4f5.0)') A(1,:)
write(*, '(4f5.0)') A(2,:)
write(*, '(4f5.0)') A(3,:)
write(*, '(4f5.0)') A(4,:)

B(1) = 2
B(2) = 3
B(3) = 3
B(4) = 2

write(*, *) 'B = ', B
call dgbsv( n, 1, 1, 1, A, 4, ipiv, B, n, info )

write(*, '(a, i0)') ' info = ', info
write(*, '(a, 4f5.0)') 'B = ', B


! Sistema Linear - Setando matriz com coordenadas de acordo com site intel
! (https://software.intel.com/en-us/mkl-developer-reference-fortran-matrix-arguments#BAND)
!     |  2 -1  0  0 | | x1 |    | 1 |
!     | -1  2 -1  0 | | x2 |    | 0 |
!     |  0 -1  2 -1 | | x3 | =  | 0 |
!     |  0  0 -1  2 | | x4 |    | 1 |


ku = 1
kl = 1

allocate(C(2*kl+ku+1, n), D(n))

C = 0

do i=1, n
    j = i
    C(kl+ku+1+i-j, j) = 2

    if(i > 0) then
        j = i-1
        C(kl+ku+1+i-j, j) = -1
    end if

    if( i < n) then
        j = i+1
        C(kl+ku+1+i-j, j) = -1
    end if
end do

do i=1, size(C, dim=1)
    write(*, '(4f5.0)') C(i, :)
end do


D = 0
D(1) = 1
D(n) = 1

call dgbsv( n, 1, 1, 1, C, 2*kl+ku+1, ipiv, D, n, info )
write(*, *) info
write(*, '(4f5.0)') D


deallocate(D)


end program