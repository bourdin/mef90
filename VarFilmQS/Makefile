PROJ      = VarFilmQS
MEF90DIR  = ../MEF90/${PETSC_ARCH}
MEF90OBJS = $(MEF90DIR)/m_MEF_Parameters.o $(MEF90DIR)/m_MEF_LinAlg.o    \
            $(MEF90DIR)/m_MEF_Types.o      $(MEF90DIR)/m_MEF_Elements.o  \
            $(MEF90DIR)/m_MEF90.o          $(MEF90DIR)/m_MEF_Utils.o         \
            $(MEF90DIR)/m_MEF_EXO.o        $(MEF90DIR)/m_MEF_MPI.o       \
            $(MEF90DIR)/m_MEF_Sieve.o 
M_FILMSTRUCTDIR  = ../m_Film_Struct/${PETSC_ARCH}
M_FILMSTRUCTOBJS = ${M_FILMSTRUCTDIR}/m_Film_Struct.o

all: VarFilmQS

#include ${PETSC_DIR}/conf/rules
#include ${PETSC_DIR}/conf/variables
# include ${TAO_DIR}/bmake/tao_common       # for tao-lite
include ${TAO_DIR}/conf/tao_base        # for tao-2.0beta7

ifneq (,$(findstring gcc,$(PETSC_ARCH)))
	MEF90_FC_FLAGS=-DWITH_TAO ${FC_MODULE_OUTPUT_FLAG}${PETSC_ARCH} -I./${PETSC_ARCH} -I../MEF90/${PETSC_ARCH} -I${M_FILMSTRUCTDIR} ${FC_FLAGS} -ffree-line-length-none -ffixed-line-length-none  -ffree-form -Wtabs	
else
	MEF90_FC_FLAGS=-DWITH_TAO ${FC_MODULE_OUTPUT_FLAG} ${PETSC_ARCH} -I./${PETSC_ARCH} -I../MEF90/${PETSC_ARCH} -I${M_FILMSTRUCTDIR} ${FC_FLAGS}
endif

.F90.o:
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${MEF90_FC_FLAGS} ${FCPPFLAGS} -o ${PETSC_ARCH}/$@ $<

${PETSC_ARCH}:
	@mkdir ${PETSC_ARCH}

VarFilmQS: ${PETSC_ARCH} VarFilmQS.o   ${PETSC_ARCH}/VarFilmQS.o \
                                      m_VarFilmQS_Types.o ${PETSC_ARCH}/m_VarFilmQS_Types.o \
                                      m_VarFilmQS_U.o ${PETSC_ARCH}/m_VarFilmQS_U.o \
                                      m_VarFilmQS_V.o ${PETSC_ARCH}/m_VarFilmQS_V.o \
                                      m_VarFilmQS_W.o ${PETSC_ARCH}/m_VarFilmQS_W.o \
                                      m_VarFilmQS_NL_U.o ${PETSC_ARCH}/m_VarFilmQS_NL_U.o \
                                      m_VarFilmQS.o ${PETSC_ARCH}/m_VarFilmQS.o \
                                      m_VarFilmQS_Post.o ${PETSC_ARCH}/m_VarFilmQS_Post.o
				${FLINKER} -o ${PETSC_ARCH}/VarFilmQS ${PETSC_ARCH}/VarFilmQS.o ${PETSC_ARCH}/m_VarFilmQS_Types.o ${PETSC_ARCH}/m_VarFilmQS_U.o ${PETSC_ARCH}/m_VarFilmQS_V.o ${PETSC_ARCH}/m_VarFilmQS_W.o ${PETSC_ARCH}/m_VarFilmQS_NL_U.o ${PETSC_ARCH}/m_VarFilmQS.o  ${PETSC_ARCH}/m_VarFilmQS_Post.o ${MEF90OBJS} ${M_FILMSTRUCTOBJS} ${TAO_FORTRAN_LIB} ${TAO_LIB} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

VarFilmQS.o:   ${PETSC_ARCH} m_VarFilmQS.o ${PETSC_ARCH}/m_VarFilmQS.o \
                                      m_VarFilmQS_Types.o ${PETSC_ARCH}/m_VarFilmQS_Types.o \
                                      m_VarFilmQS_U.o ${PETSC_ARCH}/m_VarFilmQS_U.o \
                                      m_VarFilmQS_V.o ${PETSC_ARCH}/m_VarFilmQS_V.o \
                                      m_VarFilmQS_W.o ${PETSC_ARCH}/m_VarFilmQS_W.o \
                                      m_VarFilmQS_NL_U.o ${PETSC_ARCH}/m_VarFilmQS_NL_U.o \
                                      m_VarFilmQS_Post.o ${PETSC_ARCH}/m_VarFilmQS_Post.o


m_VarFilmQS.o: ${PETSC_ARCH} m_VarFilmQS_Types.o ${PETSC_ARCH}/m_VarFilmQS_Types.o \
                                      m_VarFilmQS_U.o ${PETSC_ARCH}/m_VarFilmQS_U.o \
                                      m_VarFilmQS_V.o ${PETSC_ARCH}/m_VarFilmQS_V.o \
                                      m_VarFilmQS_W.o ${PETSC_ARCH}/m_VarFilmQS_W.o \
                                      m_VarFilmQS_NL_U.o ${PETSC_ARCH}/m_VarFilmQS_NL_U.o \
                                      m_VarFilmQS_Post.o ${PETSC_ARCH}/m_VarFilmQS_Post.o
m_VarFilmQS_U.o:  ${PETSC_ARCH} m_VarFilmQS_Types.o ${PETSC_ARCH}/m_VarFilmQS_Types.o \
                                m_VarFilmQS_Post.o ${PETSC_ARCH}/m_VarFilmQS_Post.o \
                                m_VarFilmQS_W.o ${PETSC_ARCH}/m_VarFilmQS_W.o 
m_VarFilmQS_V.o:  ${PETSC_ARCH} m_VarFilmQS_Types.o ${PETSC_ARCH}/m_VarFilmQS_Types.o \
                                m_VarFilmQS_Post.o ${PETSC_ARCH}/m_VarFilmQS_Post.o
m_VarFilmQS_W.o:  ${PETSC_ARCH} m_VarFilmQS_Types.o ${PETSC_ARCH}/m_VarFilmQS_Types.o \
                                m_VarFilmQS_Post.o ${PETSC_ARCH}/m_VarFilmQS_Post.o
m_VarFilmQS_NL_U.o:  ${PETSC_ARCH} m_VarFilmQS_Types.o ${PETSC_ARCH}/m_VarFilmQS_Types.o \
                                   m_VarFilmQS_Post.o ${PETSC_ARCH}/m_VarFilmQS_Post.o

m_VarFilmQS_Post.o:  ${PETSC_ARCH} m_VarFilmQS_Types.o ${PETSC_ARCH}/m_VarFilmQS_Types.o
m_VarFilmQS_Types.o: ${PETSC_ARCH}


clean::
	@rm -Rf ${PETSC_ARCH}
debug::
	@echo MEF90_FC_FLAGS is ${MEF90_FC_FLAGS}
	@echo FC is ${FC}
	@echo FC_FLAGS is ${FC_FLAGS}
	@echo PETSC_MAKE_STOP_ON_ERROR is ${PETSC_MAKE_STOP_ON_ERROR}
	@echo FCPPFLAGS is ${FCPPFLAGS}
	@echo FC_MODULE_FLAG is ${FC_MODULE_FLAG}
	@echo FC_MODULE_OUTPUT_FLAG is ${FC_MODULE_OUTPUT_FLAG}
