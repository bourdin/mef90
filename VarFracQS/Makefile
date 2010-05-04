PROJ      = VarFracQS
MEF90DIR  = ../MEF90/${PETSC_ARCH}
MEF90OBJS = $(MEF90DIR)/m_MEF_Parameters.o $(MEF90DIR)/m_MEF_LinAlg.o    \
            $(MEF90DIR)/m_MEF_Types.o      $(MEF90DIR)/m_MEF_Elements.o  \
            $(MEF90DIR)/m_MEF90.o          $(MEF90DIR)/m_MEF_Utils.o         \
            $(MEF90DIR)/m_MEF_EXO.o        $(MEF90DIR)/m_MEF_MPI.o       \
            $(MEF90DIR)/m_MEF_Sieve.o 
M_VARFRACSTRUCTDIR  = ../m_VarFrac_Struct/${PETSC_ARCH}
M_VARFRACSTRUCTOBJS = ${M_VARFRACSTRUCTDIR}/m_VarFrac_Struct.o

all: VarFracQS${DIM}

#include ${PETSC_DIR}/conf/rules
#include ${PETSC_DIR}/conf/variables
include ${TAO_DIR}/bmake/tao_common

ifeq (${FC},gfortran)
	MEF90_FC_FLAGS=-DWITH_TAO -DPB_${DIM} ${FC_MODULE_OUTPUT_FLAG} ${PETSC_ARCH} -I./${PETSC_ARCH} -I../MEF90/${PETSC_ARCH} -I${M_VARFRACSTRUCTDIR} ${FC_FLAGS} --free-line-length-none
else
	MEF90_FC_FLAGS=-DWITH_TAO -DPB_${DIM} ${FC_MODULE_OUTPUT_FLAG} ${PETSC_ARCH} -I./${PETSC_ARCH} -I../MEF90/${PETSC_ARCH} -I${M_VARFRACSTRUCTDIR} ${FC_FLAGS}
endif

.F90.o:
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${MEF90_FC_FLAGS} ${FCPPFLAGS} -o ${PETSC_ARCH}/${DIM}$@ $<

${PETSC_ARCH}:
	@mkdir ${PETSC_ARCH}

VarFracQS${DIM}: ${PETSC_ARCH} VarFracQS.o   ${PETSC_ARCH}/${DIM}VarFracQS.o \
                                      m_VarFracQS_Types.o ${PETSC_ARCH}/${DIM}m_VarFracQS_Types.o \
                                      m_VarFracQS_U.o ${PETSC_ARCH}/${DIM}m_VarFracQS_U.o \
                                      m_VarFracQS_V.o ${PETSC_ARCH}/${DIM}m_VarFracQS_V.o \
                                      m_VarFracQS.o ${PETSC_ARCH}/${DIM}m_VarFracQS.o \
                                      m_VarFracQS_Post.o ${PETSC_ARCH}/${DIM}m_VarFracQS_Post.o
	${FLINKER} -o ${PETSC_ARCH}/VarFracQS${DIM} ${PETSC_ARCH}/${DIM}VarFracQS.o ${PETSC_ARCH}/${DIM}m_VarFracQS_Types.o ${PETSC_ARCH}/${DIM}m_VarFracQS_U.o ${PETSC_ARCH}/${DIM}m_VarFracQS_V.o ${PETSC_ARCH}/${DIM}m_VarFracQS.o  ${PETSC_ARCH}/${DIM}m_VarFracQS_Post.o ${MEF90OBJS} ${M_VARFRACSTRUCTOBJS} ${TAO_FORTRAN_LIB} ${TAO_LIB} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

VarFracQS.o:   ${PETSC_ARCH} m_VarFracQS.o ${PETSC_ARCH}/${DIM}m_VarFracQS.o \
                                      m_VarFracQS_Types.o ${PETSC_ARCH}/${DIM}m_VarFracQS_Types.o \
                                      m_VarFracQS_U.o ${PETSC_ARCH}/${DIM}m_VarFracQS_U.o \
                                      m_VarFracQS_V.o ${PETSC_ARCH}/${DIM}m_VarFracQS_V.o \
                                      m_VarFracQS_Post.o ${PETSC_ARCH}/${DIM}m_VarFracQS_Post.o

VarFracQSRecomputeEnergy${DIM}: ${PETSC_ARCH} VarFracQSRecomputeEnergy.o   ${PETSC_ARCH}/${DIM}VarFracQSRecomputeEnergy.o \
                                          m_VarFracQS_Types.o ${PETSC_ARCH}/${DIM}m_VarFracQS_Types.o \
                                          m_VarFracQS_U.o ${PETSC_ARCH}/${DIM}m_VarFracQS_U.o \
                                          m_VarFracQS_V.o ${PETSC_ARCH}/${DIM}m_VarFracQS_V.o \
                                          m_VarFracQS.o ${PETSC_ARCH}/${DIM}m_VarFracQS.o \
                                          m_VarFracQS_Post.o ${PETSC_ARCH}/${DIM}m_VarFracQS_Post.o
	${FLINKER} -o ${PETSC_ARCH}/VarFracQSRecomputeEnergy${DIM} ${PETSC_ARCH}/${DIM}VarFracQSRecomputeEnergy.o ${PETSC_ARCH}/${DIM}m_VarFracQS_Types.o ${PETSC_ARCH}/${DIM}m_VarFracQS_U.o ${PETSC_ARCH}/${DIM}m_VarFracQS_V.o ${PETSC_ARCH}/${DIM}m_VarFracQS.o  ${PETSC_ARCH}/${DIM}m_VarFracQS_Post.o ${MEF90OBJS} ${M_VARFRACSTRUCTOBJS} ${TAO_FORTRAN_LIB} ${TAO_LIB} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

VarFracQSRecomputeEnergy.o:   ${PETSC_ARCH} m_VarFracQS.o ${PETSC_ARCH}/${DIM}m_VarFracQS.o \
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
