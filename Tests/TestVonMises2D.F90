module TestVonMises2D_mod
#include "petsc/finclude/petsc.h"
   use m_MEF90
   implicit NONE
   !!! note that this type is NOT C interoperable, which is not an issue, since we only
   !!! need SNLP to carry its address

   type :: ctx 
      Type(MEF90HookesLaw2D)  :: HookesLaw
      real(Kind = Kr)         :: YieldStress
      Type(MatS2D)            :: Strain
      Type(MatS2D)            :: OldPlasticStrain
      Type(MatS2D)            :: PlasticStrain
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
      type(MatS2D)                  :: x3D

      x3D%XX = x(1)
      x3D%YY = x(2)
      x3D%XY = x(3)

      !!! This is the fortran equivalent of casting ctx into a c_ptr
      call c_f_pointer(myctx,myctx_ptr)      
      
      !write(*,*) 'x:         ',x3D
      !write(*,*) 'HookesLaw:         ',myctx_ptr%HookesLaw
      !write(*,*) 'sigma_D:         ',ctx_ptr%sigma_D

      f(1) = ( (myctx_ptr%HookesLaw * (x3D-myctx_ptr%OldPlasticStrain)) .DotP. (x3D-myctx_ptr%OldPlasticStrain) ) /2.
      h(1) = Trace(x3D)
      g(1) = sqrt( 2.0*trace(  deviatoricPart(myctx_ptr%HookesLaw*(myctx_ptr%Strain-x3D))  *  deviatoricPart(myctx_ptr%HookesLaw*(myctx_ptr%Strain-x3D)) )) - myctx_ptr%YieldStress


   end subroutine fhg

!]   subroutine Dfhg(x,Df,Dh,Dg,myctx) bind(c)
!]      use,intrinsic :: iso_c_binding
!]   use m_MEF90
!]      real(kind=c_double)           :: x(*)
!]      real(kind=c_double)           :: Df(*)
!]      type(c_ptr)                   :: Dh
!]      type(c_ptr)                   :: Dg
!]      type(c_ptr),value             :: myctx
!]      type(ctx),pointer             :: myctx_ptr
!]
!]      real(kind=c_double),dimension(:), pointer    :: Dhptr      
!]      real(kind=c_double),dimension(:,:), pointer  :: Dgptr      
!]      type(Vect3D)                  :: x3D,Df3D
!]
!]      x3D = x(1:3)
!]
!]      !!! This is the fortran equivalent of casting ctx into a c_ptr
!]      call c_f_pointer(myctx,myctx_ptr)      
!]      Df3D = myctx_ptr%A * (x3D-myctx_ptr%p) - myctx_ptr%b
!]      Df(1) = Df3D%X
!]      Df(2) = Df3D%Y
!]      Df(3) = Df3D%Z
!]
!]      call c_f_pointer(Dh,Dhptr,[3])
!]      Dhptr = 1.
!]
!]      !!! remember that fortran uses column-major ordering whereas C uses row-major ordering
!]      !!! so that we need to compute Dg^t instead of Dg
!]      call c_f_pointer(Dg,Dgptr,[3,6])
!]      Dgptr(1:3,1) = [1., -1.,  0.]
!]      Dgptr(1:3,2) = [-1., 1.,  0.]
!]      Dgptr(1:3,3) = [1.,  0., -1.]
!]      Dgptr(1:3,4) = [-1., 0.,  1.]
!]      Dgptr(1:3,5) = [0.,  1., -1.]
!]      Dgptr(1:3,6) = [0., -1.,  1.]
!]
!]   end subroutine Dfhg
end module TestVonMises2D_mod

program testVonMises2D
#include "petsc/finclude/petsc.h"
   use,intrinsic :: iso_c_binding
   use TestVonMises2D_mod
   use m_MEF90
#ifdef MEF90_HAVE_SNLP
   use SNLPF90
   implicit NONE

   integer(kind=c_int)  :: n = 3
   integer(kind=c_int)  :: m = 1
   integer(kind=c_int)  :: p = 1
   type(SNLP),pointer   :: s
   integer              :: i,j
   integer(kind=c_int)  :: exit_code
   real(kind=c_double),dimension(:),pointer  ::x
   type(ctx),target     :: ctx_ptr
   
   ctx_ptr%HookesLaw%type = MEF90HookesLawTypeIsotropic
   ctx_ptr%HookesLaw%YoungsModulus    = 1.0_Kr
   ctx_ptr%HookesLaw%PoissonRatio     = 0.0_Kr
   ctx_ptr%HookesLaw%lambda = ctx_ptr%HookesLaw%YoungsModulus * ctx_ptr%HookesLaw%PoissonRatio / (1.0_Kr + ctx_ptr%HookesLaw%PoissonRatio) / (1.0_Kr - 2.0_Kr * ctx_ptr%HookesLaw%PoissonRatio)
   ctx_ptr%HookesLaw%mu     = ctx_ptr%HookesLaw%YoungsModulus / (1.0_Kr + ctx_ptr%HookesLaw%PoissonRatio) * .5_Kr   

   ctx_ptr%Strain = 0.0_Kr
   ctx_ptr%Strain%XX = 2.0_Kr

   ctx_ptr%OldPlasticStrain = 0.0_Kr
   ctx_ptr%PlasticStrain = 0.0_Kr
   ctx_ptr%YieldStress = 1.0_Kr
   
   allocate(x(n))
   x = ctx_ptr%PlasticStrain
   
   call SNLPNew(s,n,m,p,c_funloc(fhg),c_null_funptr,c_loc(ctx_ptr))
   call SNLPSolve(s,n,m,p,c_funloc(fhg),c_null_funptr,c_loc(ctx_ptr))
   s%show_progress = 1
   
   exit_code = SNLPL1SQP(s,x)
   write(*,*) 'exit_code: ',exit_code
   write(*,*) 'x:         ',x
   call SNLPDelete(s)
   deallocate(x)


#else
   write(*,*) 'This example needs SNLP'
#endif
end program testVonmises2D