#if MEF90_DIM == "3
#define MEF90_ELEMENT_SCAL  Element3D_Scal
#define MEF90_ELEMENT_VECT  Element3D_Vect
#define MEF90_ELEMENT_ELAST Element3D_Elast
#define MEF90_VECT          Vect3D
#define MEF90_MAT           Mat3D
#define MEF90_MATS          MatS3D
#define MEF90_TENS4OS       Tens4OS3D
#else
#define MEF90_ELEMENT_SCAL  Element2D_Scal
#define MEF90_ELEMENT_VECT  Element2D_Vect
#define MEF90_ELEMENT_ELAST Element2D_Elast
#define MEF90_VECT          Vect2D
#define MEF90_MAT           Mat2D
#define MEF90_MATS          MatS2D
#define MEF90_TENS4OS       Tens4OS2D
#endif
