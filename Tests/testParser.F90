program TestParser
#include <petsc/finclude/petsc.h>
   use m_MEF90
   use petsc
   use symengine

   implicit none(type, external)

   PetscErrorCode                                      :: ierr
   type(Basic)                                         :: f, feval
   integer, parameter                                   :: nvars = 2
   type(symbol), dimension(nvars)                       :: vars
   type(realdouble), dimension(nvars)                   :: vals

   PetscReal                                           :: theta
   PetscInt                                            :: i, j, n = 11
   character(len=256)                                  :: func, IOBuffer
   PetscBool                                           :: set

   PetscCallA(PetscInitialize(ierr))
   PetscCallA(MEF90Initialize(PETSC_COMM_WORLD, ierr))

   PetscCallA(PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-s', func, set, ierr))

   write (IOBuffer, "('Parsing: ',A,'\n')") trim(func)
   PetscCallA(PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr))

   f = parse(func)
   write (IOBuffer, '("Parsed:  ",A,"\n")') f % str()
   PetscCallA(PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr))

   vars = [Symbol("x"), Symbol("y")]
   do i = 1, N
      theta = -4.0_kr * atan(1.0_kr) + 8.0_kr * atan(1.0_kr) * (i - 1.0_kr) / (n - 1.0_kr)
      vals = [RealDouble(cos(theta)), RealDouble(sin(theta))]
      feval = f
      do j = 1, nvars
         feval = feval % subs(vars(j), vals(j))
      end do
      feval = feval % evalf()
      write (IOBuffer, '("Evaluated as a PetscReal:  ",ES12.4, ES12.4,"\n")') theta, feval % dbl()
      PetscCallA(PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr))
   end do

   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
end program testParser
