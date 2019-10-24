program BandSolver

implicit none

integer :: n = 200

integer, allocatable :: ipiv(:)

integer info

integer ku, kl, j, i

real(8), allocatable :: A(:, :), B(:), sol(:), x(:)

real(8) :: PI = acos(-1.0d0), dx, normb

character(len=32) :: format


ku = 1

kl = 1

write(*,'(a15, a15, a15)') 'N', 'Sol Diff', 'info'


do n=20, 1000, 40


    allocate(A(2*kl + ku + 1, n))
    allocate(x(n), sol(n), B(n), ipiv(n))
    A = 0.0d0

    dx = PI / (n-1)

    do i = 1,n
        x(i) = (i-1) * dx

        sol(i) = sin(x(i))
        B(i) = sin(x(i))
    end do

    if(n < 10) then
        write(*, *) " Coordinates Vector"
        write(*, *) x

        write(*, *) " Rhs Vector"
        write(*, *) B
    end if

    do i=2,n-1

        j = i
        A(kl+ku+1+i-j, j) = 2 / (dx*dx)

        j = i-1
        A(kl+ku+1+i-j, j) = -1 / (dx*dx)

        j = i+1
        A(kl+ku+1+i-j, j) = -1 / (dx*dx)

    end do

    i = 1 ; j = 1
    A(kl+ku+1+i-j, j) = 1

    i = n ; j = n
    A(kl+ku+1+i-j, j) = 1

    if(n < 10) then
        write(*, *) " A matrix"
        write(format, '(a,i0,a)') '(', n, 'f5.0)'
        do i=kl+1, 2*kl + ku + 1
            write(*,format) A(i, :)
        end do
    end if

    normb = norm2(B)

    call dgbsv( n, 1, 1, 1, A, 2*kl+ku+1, ipiv, B, n, info )


    write(*, '(i15, e15.5, i15)') n, norm2(B - sol)/normb, info


    open(10, file="out.txt")
    do i=1,n
        write(10, *) B(i)
    end do
    close(10)
    deallocate(A, B, sol, x, ipiv)
end do

end program


