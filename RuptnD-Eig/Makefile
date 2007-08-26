PROJ      = Rupt
MEF90DIR  = ../MEF90
MEF90OBJS = $(MEF90DIR)/m_Constantes.o    $(MEF90DIR)/m_AlgebLin.o  \
            $(MEF90DIR)/m_MEF_Types.o     $(MEF90DIR)/m_MEF_Gauss.o \
            $(MEF90DIR)/m_MEF90.o         $(MEF90DIR)/m_Utils.o     \
            $(MEF90DIR)/m_MEF_Ensight.o   $(MEF90DIR)/m_MEF_EXO.o   \
            $(MEF90DIR)/m_MEF_SD.o        $(MEF90DIR)/m_MEF_PU.o    \
            $(MEF90DIR)/m_MEF_Communication.o

M_RUPT_STRUCT_DIR  = ../m_Rupt_Struct
M_RUPT_STRUCT_OBJS = $(M_RUPT_STRUCT_DIR)/m_Rupt_Struct.o

EXO_CFLAGS = -I$(HPCDIST_DIR)/include
EXO_LFLAGS = -L$(HPCDIST_DIR)/lib -lexoIIv2for -lexoIIv2c -lnetcdf -lmetis -lparmetis

all: ${PROJ}${DIM}

#include ${PETSC_DIR}/bmake/common/base
include ${SLEPC_DIR}/bmake/slepc_common

OBJS=    ${PROJ}${DIM}.o m_Rupt${DIM}_Vars.o m_Rupt${DIM}_Proc.o            \
         m_Rupt${DIM}_U.o m_Rupt${DIM}_V.o m_Rupt${DIM}_Ener.o 

FFLAGS = -I${MEF90DIR} -I${M_RUPT_STRUCT_DIR} ${EXO_CFLAGS} -DPB_${DIM} 

${PROJ}${DIM}: ${PROJ}${DIM}.o
	-${FLINKER} -o ${PROJ}${DIM} ${OBJS} ${MEF90OBJS} ${M_RUPT_STRUCT_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SLEPC_LIB} ${EXO_LFLAGS}

${PROJ}:
        -@ echo use 'make DIM=2D' or 'make DIM=3D' or 'make DIM=2DA'

# Dependences
m_Rupt${DIM}_Vars.o:  ${MEF90OBJS}
m_Rupt${DIM}_U.o:     ${MEF90OBJS} 
m_Rupt${DIM}_V.o:     ${MEF90OBJS} 
m_Rupt${DIM}_Ener.o:  ${MEF90OBJS} 
m_Rupt${DIM}_Proc.o:  ${MEF90OBJS} m_Rupt${DIM}_Vars.o m_Rupt${DIM}_U.o m_Rupt${DIM}_V.o m_Rupt${DIM}_Ener.o
Rupt${DIM}.o:         ${MEF90OBJS} m_Rupt${DIM}_Proc.o

