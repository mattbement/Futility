
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                              Copyright (C) 2012                              !
!                   The Regents of the University of Michigan                  !
!              MPACT Development Group and Prof. Thomas J. Downar              !
!                             All rights reserved.                             !
!                                                                              !
! Copyright is reserved to the University of Michigan for purposes of          !
! controlled dissemination, commercialization through formal licensing, or     !
! other disposition. The University of Michigan nor any of their employees,    !
! makes any warranty, express or implied, or assumes any liability or          !
! responsibility for the accuracy, completeness, or usefulness of any          !
! information, apparatus, product, or process disclosed, or represents that    !
! its use would not infringe privately owned rights. Reference herein to any   !
! specific commercial products, process, or service by trade name, trademark,  !
! manufacturer, or otherwise, does not necessarily constitute or imply its     !
! endorsement, recommendation, or favoring by the University of Michigan.      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Module provides a preconditioner type and methods to set up different
!> preconditioners
!>
!> The precontioner type owns its own matrix. It is initialized with a matrix. 
!> In general any third party library (TPL) that may be interfaced with this 
!> module should be implemented such that it's optional.
!>
!> Currently supported TPLs include:
!>  - 
!>
!>
!> @par Module Dependencies
!>  - @ref IntrType "IntrType": @copybrief IntrType
!>  - @ref BLAS "BLAS": @copybrief BLAS
!>  - @ref Times "Times": @copybrief Times
!>  - @ref ExceptionHandler "ExceptionHandler": @copybrief ExceptionHandler
!>  - @ref Allocs "Allocs": @copybrief Allocs
!>  - @ref ParameterLists "ParameterLists": @copybrief ParameterLists
!>  - @ref ParallelEnv "ParallelEnv": @copybrief ParallelEnv
!>  - @ref VectorTypes "VectorTypes": @copybrief VectorTypes
!>  - @ref MatrixTypes "MatrixTypes": @copybrief MatrixTypes
!>
!> @par EXAMPLES
!> @code
!> 
!> @endcode
!>
!> @author Aaron Graham, Andrew Gerlach
!>   @date 11/01/2013
!>
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE PreconditionerTypes

  USE IntrType
  USE BLAS
  USE ExceptionHandler
  USE Allocs
  USE ParameterLists
  USE ParallelEnv
  USE VectorTypes
  USE MatrixTypes

  IMPLICIT NONE
  PRIVATE

  !
  ! List of Public members
  PUBLIC :: PreconditionerType
  PUBLIC :: LU_PreCondType
  PUBLIC :: ILU_PreCondType
  PUBLIC :: BILU_PreCondType
  PUBLIC :: ePreCondType


  TYPE,ABSTRACT :: PreConditionerType
    LOGICAL(SBK) :: isInit=.FALSE.
    CLASS(MatrixType),POINTER :: A=>NULL()

    CONTAINS
      PROCEDURE(precond_init_absintfc),DEFERRED,PASS :: init
      PROCEDURE(precond_absintfc),DEFERRED,PASS :: clear
      PROCEDURE(precond_absintfc),DEFERRED,PASS :: setup
      PROCEDURE(precond_apply_absintfc),DEFERRED,PASS :: apply
  ENDTYPE PreConditionerType

  TYPE,ABSTRACT,EXTENDS(PreConditionerType) :: LU_PreCondType
    CLASS(MatrixType),POINTER :: L
    CLASS(MatrixType),POINTER :: U

    CONTAINS
      PROCEDURE,PASS :: init => init_LU_PreCondType
      PROCEDURE,PASS :: clear => clear_LU_PreCondType
      PROCEDURE(precond_LU_absintfc),DEFERRED,PASS :: setup
      PROCEDURE,PASS :: apply => apply_LU_PreCondType
  ENDTYPE LU_PreCondType

  TYPE,EXTENDS(LU_PreCondType) :: ILU_PreCondType
    CONTAINS
      PROCEDURE,PASS :: setup => setup_ILU_PreCondType
  ENDTYPE ILU_PreCondType
  
  TYPE,EXTENDS(LU_PreCondType) :: BILU_PreCondType
    INTEGER(SIK) :: nPlane
    INTEGER(SIK) :: nPin
    INTEGER(SIK) :: nGrp
    INTEGER(SIK) :: BILUType  !for 2D, 3D and other variants
    CONTAINS
      PROCEDURE,PASS :: setup => setup_BILU_PreCondType
  ENDTYPE BILU_PreCondType

  ABSTRACT INTERFACE
    SUBROUTINE precond_init_absintfc(PC,A)
      IMPORT :: PreconditionerType,Matrixtype
      CLASS(PreconditionerType),INTENT(INOUT) :: PC
      CLASS(MatrixType),TARGET,INTENT(IN) :: A
    ENDSUBROUTINE precond_init_absintfc
  ENDINTERFACE

  ABSTRACT INTERFACE
    SUBROUTINE precond_apply_absintfc(PC,v)
      IMPORT :: PreconditionerType,VectorType
      CLASS(PreconditionerType),INTENT(INOUT) :: PC
      CLASS(VectorType),INTENT(INOUT) :: v
    ENDSUBROUTINE precond_apply_absintfc
  ENDINTERFACE

  ABSTRACT INTERFACE
    SUBROUTINE precond_absintfc(PC)
      IMPORT :: PreconditionerType
      CLASS(PreconditionerType),INTENT(INOUT) :: PC
    ENDSUBROUTINE precond_absintfc

    SUBROUTINE precond_LU_absintfc(PC)
      IMPORT :: LU_PreCondType
      CLASS(LU_PreCondType),INTENT(INOUT) :: PC
    ENDSUBROUTINE precond_LU_absintfc
  ENDINTERFACE

  CHARACTER(LEN=*),PARAMETER :: modName='PreconditionerTypes'

  TYPE(ExceptionHandlerType),POINTER,SAVE :: ePreCondType => NULL()
!
!===============================================================================
  CONTAINS
!
!-------------------------------------------------------------------------------
!> @brief Initializes the Linear Solver Type with a parameter list
!> @param pList the parameter list
!>
!> @param solver The linear solver to act on
    SUBROUTINE init_LU_PreCondtype(PC,A)
      CLASS(LU_PrecondType),INTENT(INOUT) :: PC
      CLASS(MatrixType),TARGET,INTENT(IN) :: A

      CHARACTER(LEN=*),PARAMETER :: myName='init_LU_PreCondType'
      LOGICAL(SBK) :: localalloc

      localalloc=.FALSE.
      IF(.NOT.ASSOCIATED(ePreCondType)) THEN
        ALLOCATE(ePreCondType)
        localalloc=.TRUE.
      ENDIF

      IF(PC%isinit) THEN
        CALL ePreCondType%raiseError('Incorrect input to '//modName//'::'//myName// &
          ' - Preconditioner is already initialized!')
      ELSE
        IF(.NOT.(A%isInit)) THEN
          CALL ePreCondType%raiseError('Incorrect input to '//modName//'::'//myName// &
            ' - Matrix being used for LU Preconditioner is not initialized!')
        ELSE
          PC%A => A

          ! This might not be necessary here, but not sure
          SELECTTYPE(mat => PC%A)
            CLASS IS(SparseMatrixType)
              ALLOCATE(SparseMatrixType :: PC%L)
              ALLOCATE(SparseMatrixType :: PC%U)
              PC%isInit=.TRUE.
            CLASS IS(PETScMatrixType)
              ALLOCATE(PETScMatrixType :: PC%L)
              ALLOCATE(PETScMatrixType :: PC%U)
              PC%isInit=.TRUE.
            CLASS DEFAULT
              CALL ePreCondType%raiseError('Incorrect input to '//modName//'::'//myName// &
                ' - LU Preconditioners are not supported for input matrix type!')
          ENDSELECT 

        ENDIF
      ENDIF
    ENDSUBROUTINE init_LU_PreCondtype
!
!-------------------------------------------------------------------------------
!> @brief Initializes the Linear Solver Type with a parameter list
!> @param pList the parameter list
!>
!> @param solver The linear solver to act on
    SUBROUTINE clear_LU_PreCondtype(PC)
      CLASS(LU_PrecondType),INTENT(INOUT) :: PC

      IF(ASSOCIATED(PC%A)) NULLIFY(PC%A)
      CALL PC%U%clear()
      CALL PC%L%clear()
      DEALLOCATE(PC%U)
      DEALLOCATE(PC%L)
      PC%isInit=.FALSE.

    ENDSUBROUTINE clear_LU_PreCondtype
!
!-------------------------------------------------------------------------------
!> @brief Initializes the Linear Solver Type with a parameter list
!> @param pList the parameter list
!>
!> @param solver The linear solver to act on
    SUBROUTINE apply_LU_PreCondType(PC,v)
      CLASS(LU_PrecondType),INTENT(INOUT) :: PC
      CLASS(Vectortype),INTENT(INOUT) :: v

      CHARACTER(LEN=*),PARAMETER :: myName='apply_LU_PreCondType'
      LOGICAL(SBK) :: localalloc

      localalloc=.FALSE.
      IF(.NOT.ASSOCIATED(ePreCondType)) THEN
        ALLOCATE(ePreCondType)
        localalloc=.TRUE.
      ENDIF

      IF(.NOT.(PC%isInit)) THEN
        CALL ePreCondType%raiseError('Incorrect input to '//modName//'::'//myName// &
          ' - Preconditioner is not initialized.')
      ELSEIF(.NOT.(v%isInit)) THEN
        CALL ePreCondType%raiseError('Incorrect input to '//modName//'::'//myName// &
          ' - VectorType is not initialized.')
      ELSE
        ! Perform U^-1 * L^-1 * v
      ENDIF
    ENDSUBROUTINE apply_LU_PreCondType
!
!-------------------------------------------------------------------------------
!> @brief Initializes the Linear Solver Type with a parameter list
!> @param pList the parameter list
!>
!> @param solver The linear solver to act on
    SUBROUTINE setup_ILU_PreCondtype(PC)
      CLASS(ILU_PrecondType),INTENT(INOUT) :: PC

      CHARACTER(LEN=*),PARAMETER :: myName='setup_ILU_PreCondType'
      INTEGER(SIK) :: row,col,col2,i,j,k,nL,nU,nnzL,nnzU
      REAL(SRK) :: val1,val2,val3
      LOGICAL(SBK) :: localalloc
      TYPE(ParamType) :: PL

      localalloc=.FALSE.
      IF(.NOT.ASSOCIATED(ePreCondType)) THEN
        ALLOCATE(ePreCondType)
        localalloc=.TRUE.
      ENDIF

      IF(.NOT.(PC%isinit)) THEN
        CALL ePreCondType%raiseError('Incorrect input to '//modName//'::'//myName// &
          ' - Preconditioner is not initialized!')
      ELSE
        IF(.NOT.(PC%A%isInit)) THEN
          CALL ePreCondType%raiseError('Incorrect input to '//modName//'::'//myName// &
            ' - Matrix being used for LU Preconditioner is not initialized!')
        ELSE
          ! This might not be necessary here, but not sure
          SELECTTYPE(mat => PC%A)
            CLASS IS(SparseMatrixType)
              SELECTTYPE(U => PC%U); TYPE IS(SparseMatrixType)
                SELECTTYPE(L => PC%L); TYPE IS(SparseMatrixType)
                  ! Loop over A to get initialization data for L and U
                  j=0
                  nU=0; nL=0  !number of rows
                  nnzU=0; nnzL=0 !number of non-zero elements
                  DO row=1,SIZE(mat%ia)-1
                    DO i=1,mat%ia(row+1)-mat%ia(row)
                      j=j+1
                      col=mat%ja(j)
                      ! This may be redundant since mat is sparse, but be safe for now
                      CALL mat%get(row,col,val1)
                      IF(.NOT.(val1 .APPROXEQA. 0.0_SRK)) THEN
                        IF(col == row) THEN
                          nnzU=nnzU+1 ! This location is in U
                        ELSEIF(col > row) THEN
                          nnzU=nnzU+1 !This location is in U
                        ELSE
                          nnzL=nnzL+1 !This location is in L
                        ENDIF
                      ENDIF
                    ENDDO
                    nnzL=nnzL+1 ! Account for 1's on diagonal of L
                  ENDDO
                  
                  ! Initialize L and U
                  nU=mat%n
                  nL=mat%n
                  CALL PL%add('MatrixType->n',nU)
                  CALL PL%add('MatrixType->nnz',nnzU)
                  CALL U%init(PL)
                  CALL PL%set('MatrixType->n',nL)
                  CALL PL%set('MatrixType->nnz',nnzL)
                  CALL L%init(PL)
                  CALL PL%clear()
                    
                  ! Make sure initialization worked
                  IF(.NOT.(U%isInit)) THEN
                    CALL ePreCondType%raiseError('Incorrect input to '//modName//'::'//myName// &
                      ' - In LU decomposition, U was not properly initialized!')
                  ELSEIF(.NOT.(L%isInit)) THEN
                    CALL ePreCondType%raiseError('Incorrect input to '//modName//'::'//myName// &
                      ' - In LU decomposition, L was not properly initialized!')
                  ! Now loop through A again and set values of L and U
                  ELSE
                    ! Set the shape
                    j=0
                    DO row=1,SIZE(mat%ia)-1
                      j=mat%ia(row)
                      DO i=1,mat%ia(row+1)-mat%ia(row)
                        col=mat%ja(j)
                        ! This may be redundant since mat is sparse, but be safe for now
                        CALL mat%get(row,col,val1)
                        IF(.NOT.(val1 .APPROXEQA. 0.0_SRK)) THEN
                          IF(col >= row) THEN
                            CALL U%setShape(row,col,0.0_SRK)
                          ELSE
                            CALL L%setShape(row,col,0.0_SRK)
                          ENDIF
                        ENDIF
                        j=j+1
                      ENDDO
                      CALL L%setShape(row,row,0.0_SRK)
                    ENDDO

                    ! Set the first row of U and L
                    DO i=1,mat%ia(2)-mat%ia(1)
                      col=mat%ja(i)
                      CALL mat%get(1,col,val1)
                      CALL U%set(1,col,val1)
WRITE(*,*) 'set:',1,col,val1
                    ENDDO
                    CALL L%set(1,1,1.0_SRK)
                    ! Now complete LU Decomposition
                    DO row=2,SIZE(mat%ia)-1
                      j=mat%ia(row)
                      DO i=1,mat%ia(row+1)-mat%ia(row)
                        col=mat%ja(j)
                        IF(col > row-1) EXIT
                        CALL mat%get(row,col,val1)
                        CALL mat%get(col,col,val2)
                        val2=val1/val2
                        CALL L%set(row,col,val2)
WRITE(*,*) 'set:',row,col,val2
                        DO k=i+1,mat%ia(row+1)-mat%ia(row)
                          col2=mat%ja(j-i+k)
                          CALL mat%get(row,col2,val1)
WRITE(*,*) 'get:',row,col2,val1
                          IF(.NOT.(val1 .APPROXEQA. 0.0_SRK)) THEN
                            CALL U%get(col,col2,val3)
WRITE(*,*) 'get:',col,col2,val3
                            IF(col2 < row) THEN
                              CALL L%set(row,col2,val1-val2*val3)
                            ELSE
                              CALL U%set(row,col2,val1-val2*val3)
                            ENDIF
                          ENDIF
                        ENDDO
                        j=j+1
                      ENDDO
                      CALL L%set(row,row,1.0_SRK)
                    ENDDO
                  ENDIF
                ENDSELECT
              ENDSELECT
            CLASS DEFAULT
              CALL ePreCondType%raiseError('Incorrect input to '//modName//'::'//myName// &
                ' - LU Preconditioners are not supported by input matrix type!')
          ENDSELECT 
        ENDIF
      ENDIF

      IF(localalloc) DEALLOCATE(ePreCondType)
    ENDSUBROUTINE setup_ILU_PreCondtype
!
!-------------------------------------------------------------------------------
!> @brief Sets up the BILU Preconditioner
!> @param pList the parameter list
!>
!> @param solver The linear solver to act on
    SUBROUTINE setup_BILU_PreCondtype(PC)
      CHARACTER(LEN=*),PARAMETER :: myName='setup_BILU_PreCondType'
      CLASS(BILU_PrecondType),INTENT(INOUT) :: PC
      
      INTEGER(SIK) :: row,col,col2,i,j,k,nL,nU,nnzL,nnzU
      REAL(SRK) :: val1
      LOGICAL(SBK) :: localalloc
      TYPE(ParamType) :: PL

      localalloc=.FALSE.
      IF(.NOT.ASSOCIATED(ePreCondType)) THEN
        ALLOCATE(ePreCondType)
        localalloc=.TRUE.
      ENDIF

      IF(.NOT.(PC%isinit)) THEN
        CALL ePreCondType%raiseError('Incorrect input to '//modName//'::'//myName// &
          ' - Preconditioner is not initialized!')
      ELSE
        IF(.NOT.(PC%A%isInit)) THEN
          CALL ePreCondType%raiseError('Incorrect input to '//modName//'::'//myName// &
            ' - Matrix being used for LU Preconditioner is not initialized!')
        ELSE
          ! This might not be necessary here, but not sure
          SELECTTYPE(mat => PC%A)
            CLASS IS(SparseMatrixType)
              SELECTTYPE(U => PC%U); TYPE IS(SparseMatrixType)
                SELECTTYPE(L => PC%L); TYPE IS(SparseMatrixType)
                  ! Loop over A to get initialization data for L and U
                  ! may need to change for block matrices
                  j=0
                  nU=0; nL=0  !number of rows
                  nnzU=0; nnzL=0 !number of non-zero elements
                  DO row=1,SIZE(mat%ia)-1
                    DO i=1,mat%ia(row+1)-mat%ia(row)
                      j=j+1
                      col=mat%ja(j)
                      ! This may be redundant since mat is sparse, but be safe for now
                      CALL mat%get(row,col,val1)
                      IF(.NOT.(val1 .APPROXEQA. 0.0_SRK)) THEN
                        IF(col == row) THEN
                          nnzU=nnzU+1 ! This location is in U
                        ELSEIF(col > row) THEN
                          nnzU=nnzU+1 !This location is in U
                        ELSE
                          nnzL=nnzL+1 !This location is in L
                        ENDIF
                      ENDIF
                    ENDDO
                    nnzL=nnzL+1 ! Account for 1's on diagonal of L
                  ENDDO
                  
                  ! Initialize L and U
                  nU=mat%n
                  nL=mat%n
                  CALL PL%add('MatrixType->n',nU)
                  CALL PL%add('MatrixType->nnz',nnzU)
                  CALL U%init(PL)
                  CALL PL%set('MatrixType->n',nL)
                  CALL PL%set('MatrixType->nnz',nnzL)
                  CALL L%init(PL)
                  CALL PL%clear()
                    
                ENDSELECT
              ENDSELECT
            CLASS IS(PETScMatrixType)
              SELECTTYPE(U => PC%U); TYPE IS(SparseMatrixType)
                SELECTTYPE(L => PC%L); TYPE IS(SparseMatrixType)
                  !to be supported at some point soon
                  
                  ! Initialize L and U (add preallocation eventually)
                  nU=mat%n
                  nL=mat%n
                  CALL PL%add('MatrixType->n',nU)
                  CALL U%init(PL)
                  CALL PL%set('MatrixType->n',nL)
                  CALL L%init(PL)
                  CALL PL%clear()
                  
                ENDSELECT
              ENDSELECT
              
            CLASS DEFAULT
              CALL ePreCondType%raiseError('Incorrect input to '//modName//'::'//myName// &
                ' - LU Preconditioners are not supported by input matrix type!')
          ENDSELECT 
        ENDIF
      ENDIF

      IF(localalloc) DEALLOCATE(ePreCondType)
    ENDSUBROUTINE
!
!-------------------------------------------------------------------------------
END MODULE
