all: m_MEF90.o ${PETSC_ARCH}

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include ../Makefile.inc

m_MEF_Parameters.o: ${PETSC_ARCH}   

m_MEF_LinAlg.o:     m_MEF_Parameters.o ${PETSC_ARCH}/m_MEF_Parameters.o   \
		    m_MEF_Utils.o      ${PETSC_ARCH}/m_MEF_Utils.o

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
