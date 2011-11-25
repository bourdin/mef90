PROJ = MEF90 
OBJDIR_PROJ = ${MEF90DIR}

all: m_MEF90.o 

include ../Makefile.inc 
include ../Makefile.mk

.F90.o:
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${MEF90_FC_FLAGS} ${FCPPFLAGS} -o ${MEF90DIR}/$@ $<

m_MEF_Parameters.o: dir   

m_MEF_LinAlg.o:     m_MEF_Parameters.o ${MEF90DIR}/m_MEF_Parameters.o   \
		    m_MEF_Utils.o      ${MEF90DIR}/m_MEF_Utils.o

m_MEF_Utils.o:      m_MEF_Parameters.o ${MEF90DIR}/m_MEF_Parameters.o 

m_MEF_IO.o:         m_MEF_Types.o      ${MEF90DIR}/m_MEF_Types.o        \
                    m_MEF_Parameters.o ${MEF90DIR}/m_MEF_Parameters.o 

m_MEF_Types.o:      m_MEF_LinAlg.o     ${MEF90DIR}/m_MEF_LinAlg.o

m_MEF_MPI.o:        m_MEF_Types.o      ${MEF90DIR}/m_MEF_Types.o        \
                    m_MEF_Parameters.o ${MEF90DIR}/m_MEF_Parameters.o 

m_MEF_Elements.o:   m_MEF_Types.o      ${MEF90DIR}/m_MEF_Types.o        \
                    m_MEF_Parameters.o ${MEF90DIR}/m_MEF_Parameters.o

m_MEF_BoundaryElements.o:   m_MEF_Elements.o      ${MEF90DIR}/m_MEF_Elements.o

m_MEF_Sieve.o:      m_MEF_Types.o      ${MEF90DIR}/m_MEF_Types.o

m_MEF_Ensight.o:    m_MEF_Types.o      ${MEF90DIR}/m_MEF_Types.o

m_MEF_EXO.o:        m_MEF_Types.o      ${MEF90DIR}/m_MEF_Types.o        \
                    m_MEF_MPI.o        ${MEF90DIR}/m_MEF_MPI.o          \
                    m_MEF_Elements.o   ${MEF90DIR}/m_MEF_Elements.o

m_MEF90.o:          m_MEF_Types.o              ${MEF90DIR}/m_MEF_Types.o            \
                    m_MEF_LinAlg.o             ${MEF90DIR}/m_MEF_LinAlg.o           \
                    m_MEF_Parameters.o         ${MEF90DIR}/m_MEF_Parameters.o       \
                    m_MEF_Utils.o              ${MEF90DIR}/m_MEF_Utils.o            \
                    m_MEF_MPI.o                ${MEF90DIR}/m_MEF_MPI.o              \
                    m_MEF_Elements.o           ${MEF90DIR}/m_MEF_Elements.o         \
                    m_MEF_BoundaryElements.o   ${MEF90DIR}/m_MEF_BoundaryElements.o \
                    m_MEF_EXO.o                ${MEF90DIR}/m_MEF_EXO.o              \
                    m_MEF_Sieve.o              ${MEF90DIR}/m_MEF_Sieve.o            \
                    m_MEF_IO.o                 ${MEF90DIR}/m_MEF_IO.o


