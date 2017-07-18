module TestSNLPF90_mod
#include "petsc/finclude/petsc.h"
   use m_MEF90
   implicit NONE
   !!! note that this type is NOT C interoperable, which is not an issue, since we only
   !!! need SNLP to carry its address
   type :: ctx 
      Type(MatS3D)           :: A
      Type(Vect3D)           :: b,p
      real(Kind = Kr)        :: sigmac
      Type(Tens4OS2D)        :: C
   end type ctx
contains
   subroutine fhg(x,f,h,g,myctx) bind(c)
      use,intrinsic :: iso_c_binding
      use m_MEF90
      real(kind=c_double)           :: x(*)
      real(kind=c_double)           :: f(*)
      real(kind=c_double)           :: h(*)
      real(kind=c_double)           :: g(*)
      type(c_ptr),intent(in),value  :: myctx
      type(ctx),pointer             :: myctx_ptr
      type(Vect3D)                  :: x3D

      x3D = x(1:3)
      !!! This is the fortran equivalent of casting ctx into a c_ptr
      call c_f_pointer(myctx,myctx_ptr)      
      write(*,*) 'A:      ',myctx_ptr%A
      write(*,*) 'C:      ',myctx_ptr%C
      write(*,*) 'b:      ',myctx_ptr%b
      write(*,*) 'p:      ',myctx_ptr%p
      write(*,*) 'sigmac: ',myctx_ptr%sigmac

      
      f(1) = ((myctx_ptr%A * (x3D-myctx_ptr%p)) .DotP. (x3D-myctx_ptr%p))/2. - (x3D .DotP. myctx_ptr%b)

      h(1) = sum(x(1:3))

      g(1) =  x(1) - x(2) - myctx_ptr%sigmac
      g(2) = -x(1) + x(2) - myctx_ptr%sigmac
      g(3) =  x(1) - x(3) - myctx_ptr%sigmac
      g(4) = -x(1) + x(3) - myctx_ptr%sigmac
      g(5) =  x(2) - x(3) - myctx_ptr%sigmac
      g(6) = -x(2) + x(3) - myctx_ptr%sigmac
   end subroutine fhg

   subroutine Dfhg(x,Df,Dh,Dg,myctx) bind(c)
      use,intrinsic :: iso_c_binding
   use m_MEF90
      real(kind=c_double)           :: x(*)
      real(kind=c_double)           :: Df(*)
      type(c_ptr)                   :: Dh
      type(c_ptr)                   :: Dg
      type(c_ptr),value             :: myctx
      type(ctx),pointer             :: myctx_ptr

      real(kind=c_double),dimension(:), pointer    :: Dhptr      
      real(kind=c_double),dimension(:,:), pointer  :: Dgptr      
      type(Vect3D)                  :: x3D,Df3D

      x3D = x(1:3)

      !!! This is the fortran equivalent of casting ctx into a c_ptr
      call c_f_pointer(myctx,myctx_ptr)      
      Df3D = myctx_ptr%A * (x3D-myctx_ptr%p) - myctx_ptr%b
      Df(1) = Df3D%X
      Df(2) = Df3D%Y
      Df(3) = Df3D%Z
      
      call c_f_pointer(Dh,Dhptr,[3])
      Dhptr = 1.
      
      !!! remember that fortran uses column-major ordering whereas C uses row-major ordering
      !!! so that we need to compute Dg^t instead of Dg
      call c_f_pointer(Dg,Dgptr,[3,6])
      Dgptr(1:3,1) = [1., -1.,  0.]
      Dgptr(1:3,2) = [-1., 1.,  0.]
      Dgptr(1:3,3) = [1.,  0., -1.]
      Dgptr(1:3,4) = [-1., 0.,  1.]
      Dgptr(1:3,5) = [0.,  1., -1.]
      Dgptr(1:3,6) = [0., -1.,  1.]

   end subroutine Dfhg
end module TestSNLPF90_mod

program testSNLP
#include "petsc/finclude/petsc.h"
   use,intrinsic :: iso_c_binding
   use TestSNLPF90_mod
   use m_MEF90
#ifdef MEF90_HAVE_SNLP   
   use SNLPF90
   implicit NONE
   
   integer(kind=c_int)  :: n = 3
   integer(kind=c_int)  :: m = 1
   integer(kind=c_int)  :: p = 6
   type(SNLP),pointer   :: s
   integer(kind=c_int)  :: exit_code
   real(kind=c_double),dimension(:),pointer  ::x
   type(ctx),target     :: ctx_ptr
   
   ctx_ptr%A = 0.0_Kr
   ctx_ptr%A%XX = 1.0_Kr
   ctx_ptr%A%YY = 2.0_Kr
   ctx_ptr%A%ZZ = 5.0_Kr

   ctx_ptr%b = 0.0_Kr

   ctx_ptr%p = [-1.0_Kr,2.0_Kr,-5.0_Kr]
   ctx_ptr%sigmac = 1.0_Kr
   
   ctx_ptr%C = -1.23_Kr
   ctx_ptr%C%YYYY = 999.9_Kr

   allocate(x(n))
   x = 0.
   
   call SNLPNew(s,n,m,p,c_funloc(fhg),c_funloc(Dfhg),c_loc(ctx_ptr))
   s%show_progress = 1
   
   exit_code = SNLPL1SQP(s,x)
   write(*,*) 'exit_code: ',exit_code
   write(*,*) 'x:         ',x
   call SNLPDelete(s)
   deallocate(x)
#else
   write(*,*) 'This example needs SNLP'
#endif
end program testSNLP