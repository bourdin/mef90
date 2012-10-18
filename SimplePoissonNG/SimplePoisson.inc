#include "../MEF90/mef90.inc"

#if MEF90_DIM == 3
#define M_POISSON m_Poisson3D
#define M_POISSONASSEMBLY m_PoissonAssembly3D
#define DEFAULT_ELEMENT_SHORTID MEF90_P1_Lagrange_3D_Scal_ShortID 
#define DEFAULT_ELEMENT MEF90_P1_Lagrange_3D_Scal
#define DEFAULT_MATERIAL MEF90_Mathium3D
#else
#define M_POISSON m_Poisson2D
#define M_POISSONASSEMBLY m_PoissonAssembly2D
#define DEFAULT_ELEMENT_SHORTID MEF90_P1_Lagrange_2D_Scal_ShortID 
#define DEFAULT_ELEMENT MEF90_P1_Lagrange_2D_Scal
#define DEFAULT_MATERIAL MEF90_Mathium2D
#endif


