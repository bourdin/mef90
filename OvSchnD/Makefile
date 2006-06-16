PROJ=OvSch
VER=1

all: $(PROJ)$(DIM)-$(VER) 

MEF90VER        = 0.12.0
MEF90DIR        = $(HOME)/mef90/MEF90-$(MEF90VER)
MEF90OBJS       = $(MEF90DIR)/m_Constantes.o    $(MEF90DIR)/m_AlgebLin.o  \
                  $(MEF90DIR)/m_MEF_Types.o     $(MEF90DIR)/m_MEF_Tools.o \
                  $(MEF90DIR)/m_MEF_IO.o        $(MEF90DIR)/m_MEF_Gauss.o \
                  $(MEF90DIR)/m_MEF90.o         $(MEF90DIR)/m_Utils.o     \
                  $(MEF90DIR)/m_MEF_Ensight.o   $(MEF90DIR)/m_MEF_EXO.o   \
                  $(MEF90DIR)/m_MEF_SD.o        $(MEF90DIR)/m_MEF_PU.o    \
                  $(MEF90DIR)/m_MEF_Communication.o

M_Poisson_STRUCT_DIR=../m_Poisson_Struct-0.0.1
M_Poisson_STRUCT_OBJS=${M_Poisson_STRUCT_DIR}/m_Poisson.o

EXO_CFLAGS      = -I$(HPCDIST_DIR)/include
EXO_LFLAGS      = -L$(HPCDIST_DIR)/lib -lexoIIv2for -lexoIIv2c -lnetcdf

all: ${PROJ}${DIM}-${VER}

include ${PETSC_DIR}/bmake/common/base

OBJS=  $(PROJ)$(DIM).o m_OvSch$(DIM)_Vars.o m_OvSch$(DIM)_Proc.o

FFLAGS = -I${MEF90DIR} -I${M_Poisson_STRUCT_DIR} ${EXO_CFLAGS} -nothreads     \
         -DPB_${DIM} 


$(PROJ)$(DIM)-${VER}: ${PROJ}${DIM}.o
	-${FLINKER} -o ${PROJ}${DIM}-${VER} ${OBJS} ${MEF90OBJS}   \
          ${M_Poisson_STRUCT_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}          \
          ${EXO_LFLAGS} 

$(PROJ)-$(VER):
	-@ echo use 'make DIM=2D' or 'make DIM=3D'
# Dependences
m_OvSch$(DIM)_Vars.o: $(MEF90OBJS)
m_OvSch$(DIM)_Proc.o: $(MEF90OBJS) m_OvSch$(DIM)_Vars.o
OvSch$(DIM).o:  $(MEF90OBJS) m_OvSch$(DIM)_Vars.o m_OvSch$(DIM)_Proc.o 
