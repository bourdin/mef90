#include "../MEF90/mef90.inc"

#if MEF90_DIM == 3
#define M_POISSON_TYPES m_Poisson3D_Types
#define M_POISSONCELLSETPROPERTY_INTERFACE m_PoissonCellSetProperty3D_Interface
#define M_POISSONVERTEXSETPROPERTY_INTERFACE m_PoissonVertexSetProperty3D_Interface
#define M_POISSON m_Poisson3D
#else
#define M_POISSON_TYPES m_Poisson2D_Types
#define M_POISSONCELLSETPROPERTY_INTERFACE m_PoissonCellSetProperty2D_Interface
#define M_POISSONVERTEXSETPROPERTY_INTERFACE m_PoissonVertexSetProperty2D_Interface
#define M_POISSON m_Poisson2D
#endif


