PROJ=MEF90

MEF90OBJS	=     m_Constantes.o    m_AlgebLin.o      \
                  m_MEF_Types.o                       \
                  m_MEF90.o         m_Utils.o         \
                  m_MEF_MPI.o       m_MEF_Elements.o  \
                  m_MEF_EXO.o       m_MEF_Sieve.o

all: m_MEF90.o

include ${PETSC_DIR}/conf/base

FFLAGS = ${EXO_CFLAGS}

#Dependencies:
m_Constantes:    
m_AlgebLin.o:       m_Constantes.o
m_Utils.o:          m_Constantes.o
m_MEF_Types.o:      m_AlgebLin.o
m_MEF_MPI.o:        m_Constantes.o m_MEF_Types.o
m_MEF_Elements.o:   m_MEF_Types.o m_Utils.o m_Constantes.o
m_MEF_Sieve.o:      m_MEF_Types.o
m_MEF_Ensight.o:    m_MEF_Types.o
m_MEF_EXO.o:        m_MEF_Types.o m_MEF_MPI.o m_MEF_Elements.o
m_MEF90.o:          m_MEF_Types.o m_AlgebLin.o m_Constantes.o m_Utils.o m_MEF_MPI.o m_MEF_Elements.o m_MEF_EXO.o m_MEF_Sieve.o

really_clean: clean
	- rm -f *~ core* *mod

