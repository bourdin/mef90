MODULE_OPTS = ${FC_MODULE_FLAG}${PETSC_ARCH}
FFLAGS = --free-line-length-none

all: m_MEF90.o ${PETSC_ARCH}

include ${PETSC_DIR}/conf/base

.F90.o:
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${MODULE_OPTS} ${FFLAGS} ${FC_FLAGS} ${FCPPFLAGS} -o ${PETSC_ARCH}/$@ $<

debug:
	@echo FC is ${FC}
	@echo FFLAGS is ${FFLAGS}
	@echo FC_FLAGS is ${FC_FLAGS}
	@echo PETSC_MAKE_STOP_ON_ERROR is ${PETSC_MAKE_STOP_ON_ERROR}
	@echo FCPPFLAGS is $FCPPFLAGS}

${PETSC_ARCH}:
	@mkdir ${PETSC_ARCH}

m_MEF_Parameters.o: ${PETSC_ARCH}   

m_MEF_LinAlg.o:     m_MEF_Parameters.o ${PETSC_ARCH}/m_MEF_Parameters.o

m_MEF_Utils.o:      m_MEF_Parameters.o ${PETSC_ARCH}/m_MEF_Parameters.o

m_MEF_Types.o:      m_MEF_LinAlg.o     ${PETSC_ARCH}/m_MEF_LinAlg.o

m_MEF_MPI.o:        m_MEF_Types.o      ${PETSC_ARCH}/m_MEF_Types.o        \
                    m_MEF_Parameters.o ${PETSC_ARCH}/m_MEF_Parameters.o 

m_MEF_Elements.o:   m_MEF_Types.o      ${PETSC_ARCH}/m_MEF_Types.o        \
                    m_MEF_Parameters.o ${PETSC_ARCH}/m_MEF_Parameters.o

m_MEF_Sieve.o:      m_MEF_Types.o      ${PETSC_ARCH}/m_MEF_Types.o

m_MEF_Ensight.o:    m_MEF_Types.o      ${PETSC_ARCH}/m_MEF_Types.o

m_MEF_EXO.o:        m_MEF_Types.o      ${PETSC_ARCH}/m_MEF_Types.o        \
                    m_MEF_MPI.o        ${PETSC_ARCH}/m_MEF_MPI.o          \
                    m_MEF_Elements.o   ${PETSC_ARCH}/m_MEF_Elements.o

m_MEF90.o:          m_MEF_Types.o      ${PETSC_ARCH}/m_MEF_Types.o        \
                    m_MEF_LinAlg.o     ${PETSC_ARCH}/m_MEF_LinAlg.o       \
                    m_MEF_Parameters.o ${PETSC_ARCH}/m_MEF_Parameters.o   \
                    m_MEF_Utils.o      ${PETSC_ARCH}/m_MEF_Utils.o        \
                    m_MEF_MPI.o        ${PETSC_ARCH}/m_MEF_MPI.o          \
                    m_MEF_Elements.o   ${PETSC_ARCH}/m_MEF_Elements.o     \
                    m_MEF_EXO.o        ${PETSC_ARCH}/m_MEF_EXO.o          \
                    m_MEF_Sieve.o      ${PETSC_ARCH}/m_MEF_Sieve.o

clean::
	@rm -Rf ${PETSC_ARCH}/*o ${PETSC_ARCH}/*mod
	@rmdir ${PETSC_ARCH}

restore:
	@hg revert *.F90
	@rm *.orig
