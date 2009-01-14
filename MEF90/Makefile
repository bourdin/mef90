PROJ=MEF90

MEF90OBJS	=     m_MEF_Parameters.o    m_MEF_LinAlg.o      \
                  m_MEF_Types.o                             \
                  m_MEF90.o         m_MEF_Utils.o           \
                  m_MEF_MPI.o       m_MEF_Elements.o        \
                  m_MEF_EXO.o       m_MEF_Sieve.o

all: m_MEF90.o

include ${PETSC_DIR}/conf/base

FFLAGS = ${EXO_CFLAGS}

#Dependencies:
m_MEF_Parameters.o:    
m_MEF_LinAlg.o:     m_MEF_Parameters.o
m_MEF_Utils.o:      m_MEF_Parameters.o
m_MEF_Types.o:      m_MEF_LinAlg.o
m_MEF_MPI.o:        m_MEF_Parameters.o m_MEF_Types.o
m_MEF_Elements.o:   m_MEF_Types.o m_MEF_Parameters.o
m_MEF_Sieve.o:      m_MEF_Types.o
m_MEF_Ensight.o:    m_MEF_Types.o
m_MEF_EXO.o:        m_MEF_Types.o m_MEF_MPI.o m_MEF_Elements.o
m_MEF90.o:          m_MEF_Types.o m_MEF_LinAlg.o m_MEF_Parameters.o m_MEF_Utils.o m_MEF_MPI.o m_MEF_Elements.o m_MEF_EXO.o m_MEF_Sieve.o

really_clean: clean
	- rm -f *~ core* *mod

