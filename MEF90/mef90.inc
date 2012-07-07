#if MEF90_DIM == 3
#define MEF90_ELEMENT_SCAL  Element3D_Scal
#define MEF90_ELEMENT_VECT  Element3D_Vect
#define MEF90_ELEMENT_ELAST Element3D_Elast
#define MEF90_VECT          Vect3D
#define MEF90_MAT           Mat3D
#define MEF90_MATS          MatS3D
#define MEF90_TENS4OS       Tens4OS3D
#define MEF90_MATPROP       MEF90_MatProp3D_Type
#define SIZEOFMEF90_VECT    3
#define SIZEOFMEF90_MAT     9
#define SIZEOFMEF90_MATS    6
#define SIZEOFMEF90_TENS4OS 21
#else
#define MEF90_ELEMENT_SCAL  Element2D_Scal
#define MEF90_ELEMENT_VECT  Element2D_Vect
#define MEF90_ELEMENT_ELAST Element2D_Elast
#define MEF90_VECT          Vect2D
#define MEF90_MAT           Mat2D
#define MEF90_MATS          MatS2D
#define MEF90_TENS4OS       Tens4OS2D
#define MEF90_MATPROP       MEF90_MatProp2D_Type
#define SIZEOFMEF90_VECT    2
#define SIZEOFMEF90_MAT     4
#define SIZEOFMEF90_MATS    3
#define SIZEOFMEF90_TENS4OS 6
#endif
