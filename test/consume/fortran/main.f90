! Pulls one object out of each installed archive so the link is real and
! not an empty-archive no-op. Both routines are pure arithmetic and need
! no MPI_Init, so this runs outside mpirun.
program consume_fortran
   implicit none
   double precision, external :: ddot
   integer, external :: numroc
   double precision :: x(3), y(3), d
   integer :: n

   x = [1.0d0, 2.0d0, 3.0d0]
   y = [4.0d0, 5.0d0, 6.0d0]

   ! ddot from eplinalg::blas
   d = ddot(3, x, 1, y, 1)
   if (abs(d - 32.0d0) > 1.0d-12) then
      write (*, *) 'FAIL: ddot returned ', d, ' expected 32'
      stop 1
   end if

   ! numroc from eplinalg::scalapack: rows of a 100-row matrix in
   ! 10-row blocks landing on process 0 of a 1-process grid.
   n = numroc(100, 10, 0, 0, 1)
   if (n /= 100) then
      write (*, *) 'FAIL: numroc returned ', n, ' expected 100'
      stop 1
   end if

   write (*, *) 'OK: Fortran-only consumer linked and ran'
end program consume_fortran
