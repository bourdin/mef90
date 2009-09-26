PROJ      = VarFracQS
MEF90DIR  = ../MEF90/${PETSC_ARCH}
MEF90OBJS = $(MEF90DIR)/m_MEF_Parameters.o $(MEF90DIR)/m_MEF_LinAlg.o    \
            $(MEF90DIR)/m_MEF_Types.o      $(MEF90DIR)/m_MEF_Elements.o  \
            $(MEF90DIR)/m_MEF90.o          $(MEF90DIR)/m_MEF_Utils.o         \
            $(MEF90DIR)/m_MEF_EXO.o        $(MEF90DIR)/m_MEF_MPI.o       \
            $(MEF90DIR)/m_MEF_Sieve.o 
M_VARFRACDIR  = ../m_VarFrac_Struct/${PETSC_ARCH}
M_VARFRACOBJS = ${M_VARFRACDIR}/m_VarFrac_Struct.o

all: VarFracQS${DIM}

#include ${PETSC_DIR}/conf/base
include ${TAO_DIR}/bmake/tao_common

FFLAGS = -DWITH_TAO -DPB_${DIM} -module ${PETSC_ARCH} -I. -I${MEF90DIR} -I${M_VARFRACDIR}
# -ftrace=full -ffree-line-length-huge

.F90.o: 
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${FFLAGS} ${FC_FLAGS} ${FCPPFLAGS} -o ${PETSC_ARCH}/${DIM}$@ $<

${PETSC_ARCH}:
	@mkdir ${PETSC_ARCH}

VarFracQS${DIM}: ${PETSC_ARCH} VarFracQS.o   ${PETSC_ARCH}/${DIM}VarFracQS.o \
                                      m_VarFracQS_Types.o ${PETSC_ARCH}/${DIM}m_VarFracQS_Types.o \
                                      m_VarFracQS_U.o ${PETSC_ARCH}/${DIM}m_VarFracQS_U.o \
                                      m_VarFracQS_V.o ${PETSC_ARCH}/${DIM}m_VarFracQS_V.o \
                                      m_VarFracQS.o ${PETSC_ARCH}/${DIM}m_VarFracQS.o \
                                      m_VarFracQS_Post.o ${PETSC_ARCH}/${DIM}m_VarFracQS_Post.o
	${FLINKER} -o ${PETSC_ARCH}/VarFracQS${DIM} ${PETSC_ARCH}/${DIM}VarFracQS.o ${PETSC_ARCH}/${DIM}m_VarFracQS_Types.o ${PETSC_ARCH}/${DIM}m_VarFracQS_U.o ${PETSC_ARCH}/${DIM}m_VarFracQS_V.o ${PETSC_ARCH}/${DIM}m_VarFracQS.o  ${PETSC_ARCH}/${DIM}m_VarFracQS_Post.o ${MEF90OBJS} ${M_VARFRACOBJS} ${TAO_FORTRAN_LIB} ${TAO_LIB} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

VarFracQS.o:   ${PETSC_ARCH} m_VarFracQS.o ${PETSC_ARCH}/${DIM}m_VarFracQS.o \
                                      m_VarFracQS_Types.o ${PETSC_ARCH}/${DIM}m_VarFracQS_Types.o \
                                      m_VarFracQS_U.o ${PETSC_ARCH}/${DIM}m_VarFracQS_U.o \
                                      m_VarFracQS_V.o ${PETSC_ARCH}/${DIM}m_VarFracQS_V.o \
                                      m_VarFracQS_Post.o ${PETSC_ARCH}/${DIM}m_VarFracQS_Post.o

m_VarFracQS.o: ${PETSC_ARCH} m_VarFracQS_Types.o ${PETSC_ARCH}/${DIM}m_VarFracQS_Types.o \
                                      m_VarFracQS_U.o ${PETSC_ARCH}/${DIM}m_VarFracQS_U.o \
                                      m_VarFracQS_V.o ${PETSC_ARCH}/${DIM}m_VarFracQS_V.o \
                                      m_VarFracQS_Post.o ${PETSC_ARCH}/${DIM}m_VarFracQS_Post.o
m_VarFracQS_U.o:  ${PETSC_ARCH} m_VarFracQS_Types.o ${PETSC_ARCH}/${DIM}m_VarFracQS_Types.o \
                                m_VarFracQS_Post.o ${PETSC_ARCH}/${DIM}m_VarFracQS_Post.o
m_VarFracQS_V.o:  ${PETSC_ARCH} m_VarFracQS_Types.o ${PETSC_ARCH}/${DIM}m_VarFracQS_Types.o \
                                m_VarFracQS_Post.o ${PETSC_ARCH}/${DIM}m_VarFracQS_Post.o
m_VarFracQS_Post.o:  ${PETSC_ARCH} m_VarFracQS_Types.o ${PETSC_ARCH}/${DIM}m_VarFracQS_Types.o
m_VarFracQS_Types.o: ${PETSC_ARCH}


clean::
	@rm -Rf ${PETSC_ARCH}/*o ${PETSC_ARCH}/*mod
	@rmdir ${PETSC_ARCH}
