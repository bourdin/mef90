#define PASS(s) s
#define MEF90_APPEND(a,b)       PASS(a)b

#if MEF90_DIM == 3
#define MEF90_ELEMENT_SCAL  MEF90Element3DScal
#define MEF90_ELEMENT_VECT  MEF90Element3DVect
#define MEF90_ELEMENT_ELAST MEF90Element3DVect
#define MEF90_VECT          Vect3D
#define MEF90_MAT           Mat3D
#define MEF90_MATS          MatS3D
#define MEF90_TENS4OS       Tens4OS3D
#define MEF90_HOOKESLAW     MEF90HookesLaw3D
#define MEF90_MATPROP       MEF90MatProp3D_Type
#define SIZEOFMEF90_VECT    3_Ki
#define SIZEOFMEF90_MAT     9_Ki
#define SIZEOFMEF90_MATS    6_Ki
#define SIZEOFMEF90_TENS4OS 21_Ki
#define SIZEOFMATPROP       sizeofMEF90MatProp3D
#define MEF90_MAT_IDENTITY  MEF90Mat3DIdentity
#define MEF90_MATS_IDENTITY MEF90MatS3DIdentity
#else
#define MEF90_ELEMENT_SCAL  MEF90Element2DScal
#define MEF90_ELEMENT_VECT  MEF90Element2DVect
#define MEF90_ELEMENT_ELAST MEF90Element2DVect
#define MEF90_VECT          Vect2D
#define MEF90_MAT           Mat2D
#define MEF90_MATS          MatS2D
#define MEF90_TENS4OS       Tens4OS2D
#define MEF90_HOOKESLAW     MEF90HookesLaw2D
#define MEF90_MATPROP       MEF90MatProp2D_Type
#define SIZEOFMEF90_VECT    2_Ki
#define SIZEOFMEF90_MAT     4_Ki
#define SIZEOFMEF90_MATS    3_Ki
#define SIZEOFMEF90_TENS4OS 6_Ki
#define SIZEOFMATPROP       sizeofMEF90MatProp2D
#define MEF90_MAT_IDENTITY  MEF90Mat2DIdentity
#define MEF90_MATS_IDENTITY MEF90MatS2DIdentity
#define MEF90_DEFMECHSPLIT  MEF90DefMechSplit2D
#endif
