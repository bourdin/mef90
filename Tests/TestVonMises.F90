module TestVonMises_mod
#include "finclude/petscdef.h"
   use m_MEF90
   implicit NONE
   !!! note that this type is NOT C interoperable, which is not an issue, since we only
   !!! need SNLP to carry its address
   type :: ctx 
      Type(MEF90HookesLaw3D) :: HookesLaw
      Type(MatS3D)           :: sigma_D,p
      real(Kind = Kr)        :: sigmac
   end type ctx
contains
      !subroutine Hookelaw(Contex,lambda,mu)
      !real,intent(in)           :: lambda
      !real,intent(in)           :: mu
      !type(ctx)  :: Contex
      !Contex%lambda=lambda
      !Contex%mu=mu
      !Contex%XXXX = Contex%lambda + Contex%mu * 2.0
      !Contex%XXYY = Contex%lambda
      !Contex%XXZZ = Contex%lambda
      !Contex%XYXY = Contex%mu
      !Contex%XZXZ = Contex%mu
      !Contex%YYYY = Contex%lambda + Contex%mu * 2.0
      !Contex%YYZZ = Contex%lambda
      !Contex%YZYZ = Contex%mu
      !Contex%ZZZZ = Contex%lambda + Contex%mu * 2.0
      !Contex%XXYZ = 0.0
      !Contex%XXXZ = 0.0
      !Contex%XXXY = 0.0
      !Contex%YYXY = 0.0
      !Contex%YZXY = 0.0
      !Contex%ZZXY = 0.0
      !Contex%YYYZ = 0.0
      !Contex%YYXZ = 0.0
      !Contex%ZZYZ = 0.0
      !Contex%ZZXZ = 0.0
      !Contex%YZXZ = 0.0
      !Contex%XZXY = 0.0
      !Contex%A(:,1) = [Contex%XXXX,Contex%XXYY,Contex%XXZZ,Contex%XXYZ,Contex%XXXZ,Contex%XXXY]       ! XXXX,XXYY,XXZZ,XXYZ,XXXZ,XXXY 
      !Contex%A(:,2) = [Contex%XXYY,Contex%YYYY,Contex%YYZZ,Contex%YYYZ,Contex%YYXZ,Contex%YYXY]       !      YYYY,YYZZ,YYYZ,YYXZ,YYXY 
      !Contex%A(:,3) = [Contex%XXZZ,Contex%YYZZ,Contex%ZZZZ,Contex%ZZYZ,Contex%ZZXZ,Contex%ZZXY]       !           ZZZZ,ZZYZ,ZZXZ,ZZXY 
      !Contex%A(:,4) = [Contex%XXYZ,Contex%YYYZ,Contex%ZZYZ,Contex%YZYZ,Contex%YZXZ,Contex%YZXY]       !                YZYZ,YZXZ,YZXY 
      !Contex%A(:,6) = [Contex%XXXY,Contex%YYXY,Contex%ZZXY,Contex%YZXY,Contex%XZXY,Contex%XYXY]       !                         ,XYXY 
      !end subroutine Hookelaw

   subroutine fhg(x,f,h,g,myctx) bind(c)
      use,intrinsic :: iso_c_binding
    use m_MEF90
      real(kind=c_double)           :: x(*)
      real(kind=c_double)           :: f(*)
      real(kind=c_double)           :: h(*)
      real(kind=c_double)           :: g(*)
      type(c_ptr),intent(in),value  :: myctx
      type(ctx),pointer             :: myctx_ptr
      type(MatS3D)                  :: x3D

      !x3D = x(1:6)

      x3D%XX = x(1)
      x3D%YY = x(2)
      x3D%ZZ = x(3)
      x3D%XY = x(4)
      x3D%XZ = x(5)
      x3D%YZ = x(6)
      !!! This is the fortran equivalent of casting ctx into a c_ptr
      call c_f_pointer(myctx,myctx_ptr)      
      
      write(*,*) 'x:         ',x3D
      write(*,*) 'A:         ',myctx_ptr%HookesLaw
      !write(*,*) 'sigma_D:         ',ctx_ptr%sigma_D

      !f(1) = Trace(Tens4OS3DXMatS3D(myctx_ptr%A,(x3D-myctx_ptr%p)))

      f(1) = Trace(myctx_ptr%HookesLaw * x3D)

      !((myctx_ptr%A * (x3D-myctx_ptr%p)) .DotP. (x3D-myctx_ptr%p) )/2.

      !!f(1) = dot_product(matmul(myctx_ptr%A,x3D-myctx_ptr%p), x3D-myctx_ptr%p)/2.0
      h(1) = Trace(x3D)

      g(1) = sqrt(3.0*trace((x3D-myctx_ptr%p)*(x3D-myctx_ptr%p))/2.0) - myctx_ptr%sigmac

      !!g(1) =  x(1) - x(2) - myctx_ptr%sigmac
      !!g(2) = -x(1) + x(2) - myctx_ptr%sigmac
      !!g(3) =  x(1) - x(3) - myctx_ptr%sigmac
      !!g(4) = -x(1) + x(3) - myctx_ptr%sigmac
      !!g(5) =  x(2) - x(3) - myctx_ptr%sigmac
      !!g(6) = -x(2) + x(3) - myctx_ptr%sigmac
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

!]      real(kind=c_double),dimension(:), pointer    :: Dhptr      
!]      real(kind=c_double),dimension(:,:), pointer  :: Dgptr      
!]      type(Vect3D)                  :: x3D,Df3D

!]      x3D = x(1:3)

!]      !!! This is the fortran equivalent of casting ctx into a c_ptr
!]      call c_f_pointer(myctx,myctx_ptr)      
!]      Df3D = myctx_ptr%A * (x3D-myctx_ptr%p) - myctx_ptr%b
!]      Df(1) = Df3D%X
!]      Df(2) = Df3D%Y
!]      Df(3) = Df3D%Z
     
!]      call c_f_pointer(Dh,Dhptr,[3])
!]      Dhptr = 1.
     
!]      !!! remember that fortran uses column-major ordering whereas C uses row-major ordering
!]      !!! so that we need to compute Dg^t instead of Dg
!]      call c_f_pointer(Dg,Dgptr,[3,6])
!]      Dgptr(1:3,1) = [1., -1.,  0.]
!]      Dgptr(1:3,2) = [-1., 1.,  0.]
!]      Dgptr(1:3,3) = [1.,  0., -1.]
!]      Dgptr(1:3,4) = [-1., 0.,  1.]
!]      Dgptr(1:3,5) = [0.,  1., -1.]
!]      Dgptr(1:3,6) = [0., -1.,  1.]

!]   end subroutine Dfhg
end module TestVonMises_mod

program testVonMises
#include "finclude/petscdef.h"
   use,intrinsic :: iso_c_binding
   use TestVonMises_mod
   use m_MEF90
#ifdef MEF90_HAVE_SNLP   
   use SNLPF90
   implicit NONE
   
   integer(kind=c_int)  :: n = 6
   integer(kind=c_int)  :: m = 1
   integer(kind=c_int)  :: p = 1
   type(SNLP),pointer   :: s
   integer              :: i,j
   integer(kind=c_int)  :: exit_code
   real(kind=c_double),dimension(:),pointer  ::x
   type(ctx),target     :: ctx_ptr
   


   !ctx_ptr%A = 0.0_Kr
   !ctx_ptr%A%XX = 1.0_Kr
   !ctx_ptr%A%YY = 2.0_Kr
   !ctx_ptr%A%ZZ = 5.0_Kr

   !allocate(ctx_ptr%A(6,6))
   !call Hookelaw(ctx_ptr,1.0,1.0)

   ctx_ptr%HookesLaw%type = MEF90HookesLawTypeIsotropic
   ctx_ptr%HookesLaw%YoungsModulus    = 1.0_Kr
   ctx_ptr%HookesLaw%PoissonRatio     = 0.22_Kr
   ctx_ptr%HookesLaw%lambda = ctx_ptr%HookesLaw%YoungsModulus * ctx_ptr%HookesLaw%PoissonRatio / (1.0_Kr + ctx_ptr%HookesLaw%PoissonRatio) / (1.0_Kr - 2.0_Kr * ctx_ptr%HookesLaw%PoissonRatio)
   ctx_ptr%HookesLaw%mu     = ctx_ptr%HookesLaw%YoungsModulus / (1.0_Kr + ctx_ptr%HookesLaw%PoissonRatio) * .5_Kr      
   
   ctx_ptr%sigma_D = 0.0_Kr
   ctx_ptr%sigma_D%XX = 3.0_Kr

   ctx_ptr%p = 0.0_Kr
   ctx_ptr%p%XX = 1.0_Kr


   !ctx_ptr%p = [0.0_Kr, 0.0_Kr,0.0_Kr,0.0_Kr ,0.0_Kr, 0.0_Kr]
   ctx_ptr%sigmac = 1.0_Kr
   
   allocate(x(n))
   x = [0.0_Kr ,  0.0_Kr , 0.0_Kr , 0.0_Kr , 0.0_Kr , 0.0_Kr ]
   
   call SNLPNew(s,n,m,p,c_funloc(fhg),c_null_funptr,c_loc(ctx_ptr))
   s%show_progress = 1
   
   exit_code = SNLPL1SQP(s,x)
   write(*,*) 'exit_code: ',exit_code
   write(*,*) 'x:         ',x
   call SNLPDelete(s)
   deallocate(x)
#else
   write(*,*) 'This example needs SNLP'
#endif
end program testVonmises