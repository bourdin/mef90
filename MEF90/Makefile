all: m_MEF90.o ${PETSC_ARCH}

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

ifneq (,$(findstring xlf,$(FC)))
	FC_MODULE_OUTPUT_FLAG=-qmoddir=
endif

ifneq (,$(findstring gcc,$(PETSC_ARCH)))
	MEF90_FC_FLAGS=${FC_MODULE_OUTPUT_FLAG}${PETSC_ARCH} ${FC_FLAGS} -ffixed-line-length-none  -ffree-line-length-none -ffree-form -I../MEF90/${PETSC_ARCH}
else
	MEF90_FC_FLAGS=${FC_MODULE_OUTPUT_FLAG}${PETSC_ARCH} ${FC_FLAGS} -I../MEF90/${PETSC_ARCH}
endif

.F90.o:
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${MEF90_FC_FLAGS} ${FCPPFLAGS} -o ${PETSC_ARCH}/$@ $<

${PETSC_ARCH}:
	@mkdir ${PETSC_ARCH}

m_MEF_Parameters.o: ${PETSC_ARCH}   

m_MEF_LinAlg.o:     m_MEF_Parameters.o ${PETSC_ARCH}/m_MEF_Parameters.o   \
		    m_MEF_Utils.o      ${PETSC_ARCH}/m_MEF_Utils.o

m_MEF_Utils.o:      m_MEF_Parameters.o ${PETSC_ARCH}/m_MEF_Parameters.o 

m_MEF_IO.o:         m_MEF_Types.o      ${PETSC_ARCH}/m_MEF_Types.o        \
                    m_MEF_Parameters.o ${PETSC_ARCH}/m_MEF_Parameters.o 

m_MEF_Types.o:      m_MEF_LinAlg.o     ${PETSC_ARCH}/m_MEF_LinAlg.o

m_MEF_MPI.o:        m_MEF_Types.o      ${PETSC_ARCH}/m_MEF_Types.o        \
                    m_MEF_Parameters.o ${PETSC_ARCH}/m_MEF_Parameters.o 

m_MEF_Elements.o:   m_MEF_Types.o      ${PETSC_ARCH}/m_MEF_Types.o        \
                    m_MEF_Parameters.o ${PETSC_ARCH}/m_MEF_Parameters.o

m_MEF_BoundaryElements.o:   m_MEF_Elements.o      ${PETSC_ARCH}/m_MEF_Elements.o

m_MEF_Sieve.o:      m_MEF_Types.o      ${PETSC_ARCH}/m_MEF_Types.o

m_MEF_Ensight.o:    m_MEF_Types.o      ${PETSC_ARCH}/m_MEF_Types.o

m_MEF_EXO.o:        m_MEF_Types.o      ${PETSC_ARCH}/m_MEF_Types.o        \
                    m_MEF_MPI.o        ${PETSC_ARCH}/m_MEF_MPI.o          \
                    m_MEF_Elements.o   ${PETSC_ARCH}/m_MEF_Elements.o

m_MEF90.o:          m_MEF_Types.o              ${PETSC_ARCH}/m_MEF_Types.o            \
                    m_MEF_LinAlg.o             ${PETSC_ARCH}/m_MEF_LinAlg.o           \
                    m_MEF_Parameters.o         ${PETSC_ARCH}/m_MEF_Parameters.o       \
                    m_MEF_Utils.o              ${PETSC_ARCH}/m_MEF_Utils.o            \
                    m_MEF_MPI.o                ${PETSC_ARCH}/m_MEF_MPI.o              \
                    m_MEF_Elements.o           ${PETSC_ARCH}/m_MEF_Elements.o         \
                    m_MEF_BoundaryElements.o   ${PETSC_ARCH}/m_MEF_BoundaryElements.o \
                    m_MEF_EXO.o                ${PETSC_ARCH}/m_MEF_EXO.o              \
                    m_MEF_Sieve.o              ${PETSC_ARCH}/m_MEF_Sieve.o            \
                    m_MEF_IO.o                 ${PETSC_ARCH}/m_MEF_IO.o

clean::
	@rm -Rf ${PETSC_ARCH}/*o ${PETSC_ARCH}/*mod
	@rmdir ${PETSC_ARCH}

debug::
	@echo MEF90_FC_FLAGS is ${MEF90_FC_FLAGS}
	@echo FC is ${FC}
	@echo FC_FLAGS is ${FC_FLAGS}
	@echo PETSC_MAKE_STOP_ON_ERROR is ${PETSC_MAKE_STOP_ON_ERROR}
	@echo FCPPFLAGS is ${FCPPFLAGS}
	@echo FC_MODULE_FLAG is ${FC_MODULE_FLAG}
	@echo FC_MODULE_OUTPUT_FLAG is ${FC_MODULE_OUTPUT_FLAG}

