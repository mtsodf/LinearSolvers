      subroutine read_matrix(c_ia, c_ja, c_a, c_loads, c_neq, c_nnz) bind(C)
            use, intrinsic :: iso_c_binding
            implicit none
            character(len=1024)  :: filename
            integer              :: i, cg_limit, neq, nnz
            real(8)              :: cg_tol
            real(8), allocatable, save :: diag_precon(:), a(:), b(:), auxRealDb(:)
            integer(c_int), allocatable, save :: ia(:), ja(:)


            !C Variables
            type(c_ptr), intent(out) :: c_ia, c_ja, c_a, c_loads
            integer, intent(out)     :: c_neq, c_nnz


            write(*,*) "Entrou no fortran"

            OPEN(99, file="convert.aux", action='read')
            READ(99, '(A)') filename
            CLOSE(99)

            write(*,*) "Lendo Arquivo", trim(filename)

            OPEN(10, file=filename, action='read', form='UNFORMATTED')

            READ(10)i
            ALLOCATE(diag_precon(i))
            READ(10)i
            ALLOCATE(a(i))
            READ(10)i
            ALLOCATE(ia(i))
            READ(10)i
            ALLOCATE(ja(i))
            READ(10)i
            ALLOCATE(b(i))
          
            READ(10)cg_limit
            READ(10)neq; WRITE(*,*)'Numero de equacoes:',neq
            ALLOCATE(auxRealDb(neq))
            READ(10)cg_tol
          
            READ(10)diag_precon
            READ(10)a
            READ(10)ia
            READ(10)ja
            READ(10)b

            nnz = size(ja, dim=1)
            CLOSE(10)

            write(*,*) "Fortran"
            write(*,'(A6,3I10,A5,3I10)') "ia = ", ia(1), ia(2), ia(3), "...", ia(neq-1), ia(neq), ia(neq+1)
            write(*,'(A6,3I10,A5,3I10)') "ja = ", ja(1), ja(2), ja(3), "...", ja(nnz-2), ja(nnz-1), ja(nnz)
            write(*,'(A6,3D10.2,A5,3D10.2)') " a = ", a(1), a(2), a(3), "...", a(nnz-2), a(nnz-1), a(nnz)
            write(*,'(A6,3D10.2,A5,3D10.2)') " b = ", b(1), b(2), b(3), "...", b(neq-2), b(neq-1), b(neq+1)

            c_ia    = c_loc(ia)
            c_ja    = c_loc(ja)
            c_a     = c_loc(a)
            c_loads = c_loc(b)
            c_neq   = neq
            c_nnz   = nnz           

            write(*,*) "Fortran", c_neq, c_nnz

      end subroutine 


      subroutine write_test() bind(C)
            use, intrinsic :: iso_c_binding
            character(len=1024) :: num
            integer n

            

            num = ""

            read(99, *) num
            write(*,*) "O numero eh ", num

            close(99)
      end subroutine 