PROJ      = m_Poisson
MEF90DIR  = ../MEF90
MEF90OBJS = $(MEF90DIR)/m_Constantes.o    $(MEF90DIR)/m_AlgebLin.o      \
            $(MEF90DIR)/m_MEF_Types.o     $(MEF90DIR)/m_MEF_InitElem.o  \
            $(MEF90DIR)/m_MEF90.o         $(MEF90DIR)/m_Utils.o         \
            $(MEF90DIR)/m_MEF_MPI.o

all: ${PROJ}.o
include ${PETSC_DIR}/conf/base

OBJS=    $(PROJ).o

FFLAGS = -I${MEF90DIR}

# Dependences
$(PROJ).o: $(MEF90OBJS)

really_clean: clean
	- rm -f ${PROJ}-${VER} *~ core


