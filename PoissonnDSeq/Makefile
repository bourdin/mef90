PROJ      = PoissonSeq
MEF90DIR  = ../MEF90
MEF90OBJS = $(MEF90DIR)/m_Constantes.o    $(MEF90DIR)/m_AlgebLin.o  \
            $(MEF90DIR)/m_MEF_Types.o     $(MEF90DIR)/m_MEF_Gauss.o \
            $(MEF90DIR)/m_MEF90.o         $(MEF90DIR)/m_Utils.o     \
            $(MEF90DIR)/m_MEF_Ensight.o   $(MEF90DIR)/m_MEF_EXO.o   \
            $(MEF90DIR)/m_MEF_SD.o        $(MEF90DIR)/m_MEF_PU.o    \
            $(MEF90DIR)/m_MEF_Communication.o

M_POISSON_STRUCT_DIR=../m_Poisson_Struct
M_POISSON_STRUCT_OBJS=${M_POISSON_STRUCT_DIR}/m_Poisson.o

EXO_CFLAGS = -I$(HPCDIST_DIR)/include
EXO_LFLAGS = -L$(HPCDIST_DIR)/lib -lexoIIv2for -lexoIIv2c -lnetcdf -lmetis -lparmetis

all: ${PROJ}${DIM}

include ${PETSC_DIR}/bmake/common/base

OBJS = ${PROJ}${DIM}.o m_Poisson${DIM}_Vars.o m_Poisson${DIM}_Proc.o 

FFLAGS = -I${MEF90DIR} -I${M_POISSON_STRUCT_DIR} ${EXO_CFLAGS} -DPB_${DIM} 

${PROJ}${DIM}: ${PROJ}${DIM}.o
	-${FLINKER} -o ${PROJ}${DIM} ${OBJS} ${MEF90OBJS} ${M_POISSON_STRUCT_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${EXO_LFLAGS}


${PROJ}:
        -@ echo use 'make DIM=2D' or 'make DIM=3D'

# Dependencies
m_Poisson${DIM}_Vars.o:  ${MEF90OBJS}
m_Poisson${DIM}_Proc.o:  ${MEF90OBJS} m_Poisson${DIM}_Vars.o
PoissonSeq${DIM}.o:      ${MEF90OBJS} m_Poisson${DIM}_Vars.o m_Poisson${DIM}_Proc.o






