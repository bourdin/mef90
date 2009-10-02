PROJ      = Elast
MEF90DIR  = ../MEF90/${PETSC_ARCH}
MEF90OBJS = $(MEF90DIR)/m_MEF_Parameters.o $(MEF90DIR)/m_MEF_LinAlg.o    \
            $(MEF90DIR)/m_MEF_Types.o      $(MEF90DIR)/m_MEF_Elements.o  \
            $(MEF90DIR)/m_MEF90.o          $(MEF90DIR)/m_MEF_Utils.o         \
            $(MEF90DIR)/m_MEF_EXO.o        $(MEF90DIR)/m_MEF_MPI.o       \
            $(MEF90DIR)/m_MEF_Sieve.o 
M_VARFRACDIR  = ../m_VarFrac_Struct/${PETSC_ARCH}
M_VARFRACOBJS = ${M_VARFRACDIR}/m_VarFrac_Struct.o

all: ${PROJ}${DIM}

#include ${PETSC_DIR}/conf/base
include ${TAO_DIR}/bmake/tao_common

FFLAGS = -DPB_${DIM} ${FC_MODULE_FLAG}${PETSC_ARCH} -I. -I${MEF90DIR} -I${M_VARFRACDIR} --free-line-length-none
# -ftrace=full

.F90.o: 
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${FFLAGS} ${FC_FLAGS} ${FCPPFLAGS} -o ${PETSC_ARCH}/${DIM}$@ $<

${PETSC_ARCH}:
	@mkdir ${PETSC_ARCH}

${PROJ}${DIM}: ${PETSC_ARCH} ${PROJ}.o   ${PETSC_ARCH/${DIM}${PROJ}.o \
                             m_${PROJ}.o ${PETSC_ARCH}/${DIM}m_${PROJ}.o
	${FLINKER} -o ${PETSC_ARCH}/${PROJ}${DIM} ${PETSC_ARCH}/${DIM}${PROJ}.o  ${PETSC_ARCH}/${DIM}m_${PROJ}.o  ${MEF90OBJS} ${M_VARFRACOBJS} ${TAO_FORTRAN_LIB} ${TAO_LIB} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 

${PROJ}.o:   ${PETSC_ARCH} m_${PROJ}.o ${PETSC_ARCH}/${DIM}m_${PROJ}.o

m_${PROJ}.o: ${PETSC_ARCH}

clean::
	@rm -Rf ${PETSC_ARCH}/*o ${PETSC_ARCH}/*mod
	@rmdir ${PETSC_ARCH}
