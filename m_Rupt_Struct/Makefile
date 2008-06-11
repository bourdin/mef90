PROJ=m_Rupt_Struct
VER=0.0.6
MEF90DIR        = ../MEF90
MEF90OBJS = $(MEF90DIR)/m_Constantes.o    $(MEF90DIR)/m_AlgebLin.o  \
            $(MEF90DIR)/m_MEF_Types.o     $(MEF90DIR)/m_MEF_Gauss.o \
            $(MEF90DIR)/m_MEF90.o         $(MEF90DIR)/m_Utils.o     \
            $(MEF90DIR)/m_MEF_EXO.o       $(MEF90DIR)/m_MEF_MPI.o   \
            $(MEF90DIR)/m_MEF_Layout.o

EXO_CFLAGS      = -I$(HPCDIST_DIR)/include

all: ${PROJ}.o
include ${PETSC_DIR}/conf/base

OBJS=    $(PROJ).o

FFLAGS = -I${MEF90DIR} ${EXO_CFLAGS} 

# Dependences
$(PROJ).o: $(MEF90OBJS)

really_clean: clean
	- rm -f ${PROJ}-${VER} *~ core


