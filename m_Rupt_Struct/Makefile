PROJ      = m_Rupt_Struct
MEF90DIR  = ../MEF90
MEF90OBJS = $(MEF90DIR)/m_MEF_Parameters.o $(MEF90DIR)/m_MEF_LinAlg.o    \
            $(MEF90DIR)/m_MEF_Types.o      $(MEF90DIR)/m_MEF_Elements.o  \
            $(MEF90DIR)/m_MEF90.o          $(MEF90DIR)/m_MEF_Utils.o         \
            $(MEF90DIR)/m_MEF_EXO.o        $(MEF90DIR)/m_MEF_MPI.o       \
            $(MEF90DIR)/m_MEF_Sieve.o 

all: ${PROJ}.o
include ${PETSC_DIR}/conf/base

FFLAGS = -I${MEF90DIR}

# Dependences
$(PROJ).o: $(MEF90OBJS)


