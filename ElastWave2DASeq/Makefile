PROJ            = ElastWave2DASeq
MEF90DIR        = ../MEF90
MEF90OBJS	    = $(MEF90DIR)/m_Constantes.o    $(MEF90DIR)/m_AlgebLin.o  \
                  $(MEF90DIR)/m_MEF_Types.o     $(MEF90DIR)/m_MEF_Gauss.o \
                  $(MEF90DIR)/m_MEF90.o         $(MEF90DIR)/m_Utils.o     \
                  $(MEF90DIR)/m_MEF_Ensight.o   $(MEF90DIR)/m_MEF_EXO.o   \
                  $(MEF90DIR)/m_MEF_SD.o        $(MEF90DIR)/m_MEF_PU.o    \
                  $(MEF90DIR)/m_MEF_Communication.o
M_RUPT_STRUCT_DIR = ../m_Rupt_Struct
M_RUPT_STRUCT_OBJS = $(M_RUPT_STRUCT_DIR)/m_Rupt_Struct.o


EXO_CFLAGS = -I$(HPCDIST_DIR)/include
EXO_LFLAGS = -L$(HPCDIST_DIR)/lib -lexoIIv2for -lexoIIv2c -lnetcdf -lmetis -lparmetis

all: ${PROJ}

include ${PETSC_DIR}/bmake/common/base

OBJS=    ${PROJ}.o m_Elast2DA_Vars.o m_Elast2DA_Proc.o
FFLAGS = -I${MEF90DIR} -I${M_RUPT_STRUCT_DIR} ${EXO_CFLAGS} -DPB_${DIM} 

${PROJ}: ${PROJ}.o
	-${FLINKER} -o ${PROJ} ${OBJS} ${MEF90OBJS} ${M_RUPT_STRUCT_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${EXO_LFLAGS}

# Dependences
m_Elast2DA_Vars.o:  ${MEF90OBJS}
m_Elast2DA_Proc.o:  ${MEF90OBJS} m_Elast2DA_Vars.o
ElastWave2DASeq.o:      ${MEF90OBJS} m_Elast2DA_Vars.o m_Elast2DA_Proc.o
