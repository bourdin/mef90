PROJ            = Elast
VER             = 0.12.trunk
MEF90VER        = 0.12.0
MEF90DIR        = $(MEF90ROOT)/MEF90-$(MEF90VER)
MEF90OBJS	= $(MEF90DIR)/m_Constantes.o    $(MEF90DIR)/m_AlgebLin.o  \
                  $(MEF90DIR)/m_MEF_Types.o     $(MEF90DIR)/m_MEF_Tools.o \
                  $(MEF90DIR)/m_MEF_IO.o        $(MEF90DIR)/m_MEF_Gauss.o \
                  $(MEF90DIR)/m_MEF90.o         $(MEF90DIR)/m_Utils.o     \
                  $(MEF90DIR)/m_MEF_Ensight.o   $(MEF90DIR)/m_MEF_EXO.o   \
                  $(MEF90DIR)/m_MEF_SD.o        $(MEF90DIR)/m_MEF_PU.o    \
                  $(MEF90DIR)/m_MEF_Communication.o
M_RUPT_STRUCT_DIR = $(MEF90ROOT)/m_Rupt_Struct-0.0.6
M_RUPT_STRUCT_OBJS = $(M_RUPT_STRUCT_DIR)/m_Rupt_Struct.o


EXO_CFLAGS      = -I$(MATSCI_DIR)/SEAMS.SRC/inc
EXO_LFLAGS      = -L$(MATSCI_DIR)/SEAMS.SRC/lib -lexoIIv2for -lexoIIv2c   \
                  -lnetcdf -lexoIIv2for -lsuplib -lsupes

all: ${PROJ}${DIM}-${VER}

#include ../Makefile.${PETSC_ARCH}
include ${PETSC_DIR}/bmake/common/base

OBJS=    ${PROJ}${DIM}.o m_Elast${DIM}_Vars.o m_Elast${DIM}_Proc.o            \
         m_Elast${DIM}_Debug.o
FFLAGS = -I${MEF90DIR} -I${M_RUPT_STRUCT_DIR} ${EXO_CFLAGS} -nothreads        \
         -DMEF90_TIMING -DPB_${DIM}
#  -CB -traceback

${PROJ}${DIM}-${VER}: ${PROJ}${DIM}.o
	-${FLINKER} -nothreads -o ${PROJ}${DIM}-${VER} ${OBJS} ${MEF90OBJS}   \
          ${M_RUPT_STRUCT_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${EXO_LFLAGS}

${PROJ}-${VER}:
        -@ echo use 'make DIM=2D' or 'make DIM=3D'

# Dependences
m_Elast${DIM}_Vars.o:  ${MEF90OBJS}
m_Elast${DIM}_Proc.o:  ${MEF90OBJS} m_Elast${DIM}_Vars.o m_Elast${DIM}_Debug.o
m_Elast${DIM}_Debug.o: ${MEF90OBJS} m_Elast${DIM}_Vars.o
Elast${DIM}.o:         ${MEF90OBJS} m_Elast${DIM}_Vars.o m_Elast${DIM}_Proc.o \
                        m_Elast${DIM}_Debug.o
