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
PROGRAM testLinearSolver
  
  USE IntrType
  USE ExceptionHandler
  USE LinearSolverTypes
  USE MatrixTypes
  USE VectorTypes
  USE ParallelEnv
  USE ParameterLists
  IMPLICIT NONE
  
  TYPE(ExceptionHandlerType),POINTER :: e
  TYPE(MPI_EnvType) :: mpiTestEnv
  TYPE(ParamType) :: pList, optListLS, optListMat, vecPList
  
#ifdef HAVE_PETSC
#include <finclude/petsc.h>
#undef IS
  PetscErrorCode  :: ierr
  
  CALL PetscInitialize(PETSC_NULL_CHARACTER,ierr)
#endif

  !> set up default parameter list
  CALL optListLS%clear()
  CALL optListLS%add('LinearSolverType->TPLType',NATIVE)   
  CALL optListLS%add('LinearSolverType->solverMethod',1_SNK) ! GE or BICGSTAB
  CALL optListLS%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
  CALL optListLS%add('LinearSolverType->numberOMP',1_SNK)
  CALL optListLS%add('LinearSolverType->timerName','LinearSolver Timer')
  ! set parameters for matrices
  CALL optListLS%add('LinearSolverType->A->MatrixType->n',-1_SNK)
  CALL optListLS%add('LinearSolverType->A->MatrixType->nnz',-1_SNK)
  CALL optListLS%add('LinearSolverType->A->MatrixType->isSym',.FALSE.)
  CALL optListLS%add('LinearSolverType->A->MatrixType->matType',SPARSE)
  ! set parameters for vectors
  CALL optListLS%add('LinearSolverType->x->VectorType->n',-1_SNK)
  CALL optListLS%add('LinearSolverType->b->VectorType->n',-1_SNK)
  
  !Set up optional PL
  CALL optListMat%add('testPL->nnz',-1_SNK)
  CALL optListMat%add('testPL->isSym',.FALSE.)
  CALL optListMat%add('testPL->matType',SPARSE)
  CALL optListMat%add('testPL->MPI_Comm_ID',PE_COMM_SELF)
  
  ! Set up vector parameter list
  CALL vecPList%add('VectorType -> n',2)
  CALL vecPList%add('VectorType -> MPI_Comm_ID',PE_COMM_SELF)
  
  !Configure exception handler for test
  ALLOCATE(e)
  CALL e%setStopOnError(.FALSE.)
  CALL e%setQuietMode(.TRUE.)
  eLinearSolverType => e
  CALL mpiTestEnv%initialize(PE_COMM_SELF)

  WRITE(*,*) '==================================================='
  WRITE(*,*) 'TESTING LINEAR SOLVERS...'
  WRITE(*,*) '==================================================='
  
  WRITE(*,*) 'TESTING NORMS PROCEDURE'
  CALL testNorms()
  
  WRITE(*,*) 'TESTING LINEAR SOLVER TYPES'
  CALL testClear()
  CALL testInit()
  CALL testUpdatedA()
  CALL testDirectSolve()
  CALL testIterativeOthers()
  CALL testIterativeSolve_BICGSTAB()
  CALL testIterativeSolve_CGNR()
  CALL testIterativeSolve_GMRES()

  WRITE(*,*) '==================================================='
  WRITE(*,*) 'TESTING LINEAR SOLVERS PASSED!'
  WRITE(*,*) '==================================================='
  
#ifdef HAVE_PETSC    
  CALL PetscFinalize(ierr)
#else
  CALL mpiTestEnv%finalize()
#endif
  
  DEALLOCATE(e)
  
!
!===============================================================================
CONTAINS
!
!-------------------------------------------------------------------------------
    SUBROUTINE testClear()
      CLASS(LinearSolverType_Base),ALLOCATABLE :: thisLS
    !test Direct
      ALLOCATE(LinearSolverType_Direct :: thisLS)
      
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Direct)
      
        !first build one by hand to test clear
        thisLS%isInit=.TRUE.
        thisLS%solverMethod=BICGSTAB
        thisLS%info=2
        thisLS%isDecomposed=.TRUE.
        CALL thisLS%MPIparallelEnv%initialize(0)
        CALL thisLS%OMPparallelEnv%initialize(1)
        
        ! initialize matrix A
        ALLOCATE(DenseSquareMatrixType :: thisLS%A)
        CALL pList%clear()
        CALL pList%add('MatrixType->n',2_SNK)
        CALL pList%add('MatrixType->isSym',.TRUE.)
        CALL pList%validate(pList,optListMat)
        CALL thisLS%A%init(pList) !2x2, symmetric
        
        ! initialize vector X
        CALL vecPList%set('VectorType -> n',2)
        ALLOCATE(RealVectorType :: thisLS%X)
        CALL thisLS%X%init(vecPList)
        
        ! initialize vector b
        ALLOCATE(RealVectorType :: thisLS%b)
        CALL thisLS%b%init(vecPList)
        
        ! allocate IPIV
        ALLOCATE(thisLS%IPIV(2))
        
        ! initialize matrix M
        ALLOCATE(DenseSquareMatrixType :: thisLS%M)
        CALL pList%clear()
        CALL pList%add('MatrixType->n',10_SNK)
        CALL pList%add('MatrixType->isSym',.TRUE.)
        CALL thisLS%M%init(pList)
        
      ENDSELECT
      
      !call clear
      CALL thisLS%clear()
      
      !check results
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Direct)
        IF(thisLS%isInit .OR. thisLS%solverMethod /= -1                 &
          .OR. thisLS%isDecomposed .OR. ALLOCATED(thisLS%A)             &
          .OR. ALLOCATED(thisLS%X) .OR. thisLS%info /= 0                &
          .OR. ALLOCATED(thisLS%b) .OR. ALLOCATED(thisLS%IPIV)          &
          .OR. ALLOCATED(thisLS%M) .OR. thisLS%MPIparallelEnv%isInit()  & 
          .OR. thisLS%OMPparallelEnv%isInit()) THEN
          WRITE(*,*) 'CALL Direct%clear(...) FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      WRITE(*,*) '  Passed: CALL Direct%clear()'
      DEALLOCATE(thisLS)
      
    !test Iterative
      ALLOCATE(LinearSolverType_Iterative :: thisLS)

      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        !first build one by hand to test
        thisLS%isInit=.TRUE.
        thisLS%solverMethod=1
        thisLS%info=2
        thisLS%normType=2
        thisLS%maxIters=2
        thisLS%iters=2
        thisLS%convTol=2._SRK
        thisLS%residual=2._SRK
        thisLS%isDecomposed=.TRUE.
        CALL thisLS%MPIparallelEnv%initialize(0)
        CALL thisLS%OMPparallelEnv%initialize(1)
        
        ! initialize matrix A
        ALLOCATE(DenseSquareMatrixType :: thisLS%A)
        CALL pList%clear()
        CALL pList%add('MatrixType->n',2_SNK)
        CALL pList%add('MatrixType->isSym',.TRUE.)
        CALL pList%validate(pList,optListMat)
        CALL thisLS%A%init(pList) !2x2, symmetric
        
        ! initialize vector X
        CALL vecPList%set('VectorType -> n',2)
        ALLOCATE(RealVectorType :: thisLS%X)
        CALL thisLS%X%init(vecPList)
        
        ! initialize vector b
        ALLOCATE(RealVectorType :: thisLS%b)
        CALL thisLS%b%init(vecPList)        
        
        ! initialize matrix M
        ALLOCATE(DenseSquareMatrixType :: thisLS%M)
        CALL pList%clear()
        CALL pList%add('MatrixType->n',10_SNK)
        CALL pList%add('MatrixType->isSym',.TRUE.)
        CALL thisLS%M%init(pList)
      ENDSELECT
      
      CALL thisLS%clear()
      
      !check results
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        IF(thisLS%isInit .OR.thisLS%solverMethod == 1                  &
          .OR. ALLOCATED(thisLS%M) .OR. ALLOCATED(thisLS%A)            &
          .OR. ALLOCATED(thisLS%X) .OR. thisLS%info /= 0               &
          .OR. thisLS%MPIparallelEnv%isInit()                          &
          .OR. thisLS%OMPparallelEnv%isInit()                          &
          .OR. thisLS%normType == 2 .OR. thisLS%maxIters == 2          &
          .OR. thisLS%iters == 2 .OR. thisLS%isDecomposed              &
          .OR. thisLS%residual == 2._SRK .OR. thisLS%convTol == 2._SRK) THEN
          WRITE(*,*) 'CALL Iterative%clear() FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      CALL thisLS%clear()
      DEALLOCATE(thisLS)
      WRITE(*,*) '  Passed: CALL Iterative%clear()'
      
    ENDSUBROUTINE testClear
!
!-------------------------------------------------------------------------------
    SUBROUTINE testInit()
      CLASS(LinearSolverType_Base),ALLOCATABLE :: thisLS
      CLASS(ParamType),POINTER :: tmp
   !test Direct
      ALLOCATE(LinearSolverType_Direct :: thisLS)
      
      !Bad input
      CALL pList%clear()
      CALL pList%add('LinearSolverType->matType',SPARSE)
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',-1_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->n',1_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',1_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',1_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',1_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      CALL thisLS%clear()
      
      CALL pList%clear()
      CALL pList%add('LinearSolverType->matType',SPARSE)
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',BICGSTAB)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->n',1_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',1_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',1_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',1_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      CALL thisLS%clear()
      
      CALL pList%clear()
      CALL pList%add('LinearSolverType->matType',SPARSE)
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',BICGSTAB)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->n',1_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',1_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',1_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',1_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      CALL thisLS%clear()
      
      !first test a correct use case, with timer name
      CALL pList%clear()
      CALL pList%add('LinearSolverType->matType',SPARSE)
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',BICGSTAB)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->n',1_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',1_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',1_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',1_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Direct)
        IF(.NOT. (thisLS%isInit .AND. thisLS%solverMethod == 1 &
           .AND. thisLS%MPIparallelEnv%isInit()                &
           .AND. thisLS%OMPparallelEnv%isInit()                &
           .AND. thisLS%SolveTime%getTimerName() == 'testTimer')) THEN
          WRITE(*,*) 'CALL Direct%init(...) FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      CALL thisLS%clear()
      
      !first test a correct use case, with no timer name
      CALL pList%clear()
      CALL pList%add('LinearSolverType->matType',SPARSE)
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',BICGSTAB)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->n',1_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',1_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',1_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',1_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      SELECTTYPE(thisLS); TYPE IS (LinearSolverType_Direct)
        IF(.NOT. (thisLS%isInit .AND. thisLS%solverMethod == 1 &
           .AND. thisLS%MPIparallelEnv%isInit()                &
           .AND. thisLS%OMPparallelEnv%isInit()                &
           .AND. thisLS%SolveTime%getTimerName() == 'LinearSolver Timer')) THEN
          WRITE(*,*) 'CALL Direct%init(...) FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      CALL thisLS%clear()
      DEALLOCATE(thisLS)
      WRITE(*,*) '  Passed: CALL Direct%init(...)'
      
    !test Iterative
      ALLOCATE(LinearSolverType_Iterative :: thisLS)
      
      !Bad input
      CALL pList%clear()
      CALL pList%add('LinearSolverType->matType',SPARSE)
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',-1_SNK)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->n',1_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',1_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',1_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',1_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      CALL thisLS%clear()
      
      CALL pList%clear()
      CALL pList%add('LinearSolverType->matType',SPARSE)
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',BICGSTAB)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->n',1_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',1_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',1_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',1_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      CALL thisLS%clear()
      
      CALL pList%clear()
      CALL pList%add('LinearSolverType->matType',SPARSE)
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',BICGSTAB)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->n',1_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',1_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',1_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',1_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      CALL thisLS%clear()
      
      !first test a correct use case, with timer name
      CALL pList%clear()
      CALL pList%add('LinearSolverType->matType',SPARSE)
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',BICGSTAB)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->n',1_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',1_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',1_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',1_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        IF(.NOT. (thisLS%isInit .AND. thisLS%solverMethod == 1 &
           .AND. thisLS%MPIparallelEnv%isInit() &
           .AND. thisLS%OMPparallelEnv%isInit() &
           .AND. thisLS%SolveTime%getTimerName() == 'testTimer')) THEN
          WRITE(*,*) 'CALL Iterative%init(...) FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      CALL thisLS%clear()
      
      !first test a correct use case, with no timer name
      CALL pList%clear()
      CALL pList%add('LinearSolverType->matType',SPARSE)
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',BICGSTAB)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->n',1_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',1_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',1_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',1_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      SELECTTYPE(thisLS); TYPE IS (LinearSolverType_Iterative)
        IF(.NOT. (thisLS%isInit .AND. thisLS%solverMethod == 1 &
           .AND. thisLS%MPIparallelEnv%isInit() &
           .AND. thisLS%OMPparallelEnv%isInit() &
           .AND. thisLS%SolveTime%getTimerName() == 'LinearSolver Timer')) THEN
          WRITE(*,*) 'CALL Iterative%init(...) FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      CALL thisLS%clear()
      DEALLOCATE(thisLS)
      WRITE(*,*) '  Passed: CALL Iterative%init(...)'

    ENDSUBROUTINE testInit
!
!-------------------------------------------------------------------------------
    SUBROUTINE testUpdatedA()
      CLASS(LinearSolverType_Base),ALLOCATABLE :: thisLS

    !test Direct
      ALLOCATE(LinearSolverType_Direct :: thisLS)
      
      ! initialize matrix M
      ALLOCATE(DenseSquareMatrixType :: thisLS%M)
      CALL pList%clear()
      CALL pList%add('MatrixType->n',10_SNK)
      CALL pList%add('MatrixType->isSym',.TRUE.)
      CALL thisLS%M%init(pList)
      thisLS%isDecomposed=.TRUE.
      
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->matType',SPARSE)
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',LU)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->n',1_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',1_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',1_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',1_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Direct)
        ALLOCATE(thisLS%IPIV(10))
      ENDSELECT
      CALL thisLS%updatedA()
      !Check
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Direct)
        IF(ALLOCATED(thisLS%M) .OR. thisLS%isDecomposed &
          .OR. ALLOCATED(thisLS%IPIV)) THEN
          WRITE(*,*) 'CALL Direct%updatedA() FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      WRITE(*,*) '  Passed: CALL Direct%updatedA()'
      CALL thisLS%clear()
      DEALLOCATE(thisLS)
      
    !test Iterative
      ALLOCATE(LinearSolverType_Iterative :: thisLS)
      
      ! initialize matrix M
      ALLOCATE(DenseSquareMatrixType :: thisLS%M)
      CALL pList%clear()
      CALL pList%add('MatrixType->n',10_SNK)
      CALL pList%add('MatrixType->isSym',.TRUE.)
      CALL thisLS%M%init(pList)
      thisLS%isDecomposed=.TRUE.
      
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->matType',SPARSE)
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',BICGSTAB)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->n',1_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',1_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',1_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',1_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      CALL thisLS%updatedA()
      
      !Check
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        IF(thisLS%isDecomposed .OR. ALLOCATED(thisLS%M)) THEN
          WRITE(*,*) 'CALL Iterative%updatedA() FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      WRITE(*,*) '  Passed: CALL Iterative%updatedA()'
      CALL thisLS%clear()
      DEALLOCATE(thisLS)
      
    ENDSUBROUTINE testUpdatedA
!
!-------------------------------------------------------------------------------
    SUBROUTINE testDirectSolve()
      CLASS(LinearSolverType_Base),ALLOCATABLE :: thisLS
      REAL(SRK),ALLOCATABLE :: dummyvec(:)
      
      ALLOCATE(LinearSolverType_Direct :: thisLS)
      
    !Test GE (Dense-Square)
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',GE)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',DENSESQUARE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',3_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->isSym',.FALSE.)
      CALL pList%add('LinearSolverType->x->VectorType->n',3_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',3_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)

      !A=[1 2 3]  b=[6]   x=[*]
      !  [1 3 2]    [6]     [*]
      !  [1 3 2]    [6]     [*]
      SELECTTYPE(A => thisLS%A); TYPE IS(DenseSquareMatrixType)
        CALL A%set(1,1,1._SRK)
        CALL A%set(1,2,2._SRK)
        CALL A%set(1,3,3._SRK)
        CALL A%set(2,1,1._SRK)
        CALL A%set(2,2,3._SRK)
        CALL A%set(2,3,2._SRK)
        CALL A%set(3,1,1._SRK)
        CALL A%set(3,2,3._SRK)
        CALL A%set(3,3,2._SRK)
      ENDSELECT
      ! set all values to 6
      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(6._SRK)
      ENDSELECT
      CALL thisLS%solve()
      
      !Check the result
      IF(thisLS%info /= -1 ) THEN
        WRITE(*,*) 'CALL Direct%solve() -GE method FAILED!'
        STOP 666
      ENDIF
      CALL thisLS%clear()
      
    ! Test LU (Dense-Square)
      
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',LU)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',DENSESQUARE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',4_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->isSym',.FALSE.)
      CALL pList%add('LinearSolverType->x->VectorType->n',4_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',4_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      !A=[1 2 3 4]  b=[10]   x=[1]
      !  [1 3 2 3]    [9]      [1]
      !  [3 2 3 1]    [9]      [1]
      !  [1 1 1 1]    [4]      [1]
      SELECTTYPE(A => thisLS%A); TYPE IS(DenseSquareMatrixType)
        CALL A%set(1,1,1._SRK)
        CALL A%set(1,2,2._SRK)
        CALL A%set(1,3,3._SRK)
        CALL A%set(1,4,4._SRK)
        CALL A%set(2,1,1._SRK)
        CALL A%set(2,2,3._SRK)
        CALL A%set(2,3,2._SRK)
        CALL A%set(2,4,3._SRK)
        CALL A%set(3,1,3._SRK)
        CALL A%set(3,2,2._SRK)
        CALL A%set(3,3,3._SRK)
        CALL A%set(3,4,1._SRK)
        CALL A%set(4,1,1._SRK)
        CALL A%set(4,2,1._SRK)
        CALL A%set(4,3,1._SRK)
        CALL A%set(4,4,1._SRK)
      ENDSELECT

      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1,10._SRK)
        CALL b%set(2,9._SRK)
        CALL b%set(3,9._SRK)
        CALL b%set(4,4._SRK)
      ENDSELECT
      
      ! solve
      CALL thisLS%solve()
      
      !Check the result
      SELECTTYPE (X=>thisLS%X); TYPE IS(RealVectorType)
        IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
        ALLOCATE(dummyvec(X%n))
        CALL X%get(dummyvec)
        IF(.NOT.((dummyvec(1) .APPROXEQ. 1._SRK) &
           .AND. (dummyvec(2) .APPROXEQ. 1._SRK) &
           .AND. (dummyvec(3) .APPROXEQ. 1._SRK) &
           .AND. (dummyvec(4) .APPROXEQ. 1._SRK) &
           .AND. (thisLS%info == 0) )) THEN
          WRITE(*,*) 'CALL Direct%solve() -LU method FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      CALL thisLS%clear()
     
    ! Test LU (Tridiagonal)

      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',LU)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',TRIDIAG)
      CALL pList%add('LinearSolverType->A->MatrixType->n',3_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->isSym',.TRUE.)
      CALL pList%add('LinearSolverType->x->VectorType->n',3_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',3_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      !A=[ 4 -1  0]
      !  [-1  4 -1]
      !  [ 0 -1  4]
      SELECTTYPE(A => thisLS%A); TYPE IS(TriDiagMatrixType)
        CALL A%set(1,1,4._SRK)
        CALL A%set(1,2,-1._SRK)
        CALL A%set(2,2,4._SRK)
        CALL A%set(2,3,-1._SRK)
        CALL A%set(3,3,4._SRK)
      ENDSELECT

      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1,1._SRK)
        CALL b%set(2,2._SRK)
        CALL b%set(3,3._SRK)
      ENDSELECT

      CALL thisLS%solve()

      !Check M
      !M=[ 4   -1    0]
      !  [-.25 3.75 -1]
      !  [ 0   -.267  3.7333]
      SELECTTYPE(M => thisLS%M); TYPE IS(TriDiagMatrixType)
        IF(.NOT.((M%a(1,1) .APPROXEQ.  0.0_SRK)  &
           .AND. (M%a(1,2) .APPROXEQ. -0.25_SRK) &
           .AND. (M%a(1,3) .APPROXEQ. -0.266666666666666666_SRK) &
           .AND. (M%a(2,1) .APPROXEQ.  0.25_SRK) &
           .AND. (M%a(2,2) .APPROXEQ. (1.0_SRK/3.75_SRK)) &
           .AND. (M%a(2,3) .APPROXEQ. (1.0_SRK/3.7333333333333334_SRK)) &
           .AND. (M%a(3,1) .APPROXEQ. -1.0_SRK)  &
           .AND. (M%a(3,2) .APPROXEQ. -1.0_SRK)  &
           .AND. (M%a(3,3) .APPROXEQ.  0.0_SRK)  &
           .AND. thisLS%isDecomposed)) THEN
          WRITE(*,*) 'CALL Direct%solve() -LU method FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      
      !Check X
      SELECTTYPE(X => thisLS%X); TYPE IS(RealVectorType)
        IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
        ALLOCATE(dummyvec(X%n))
        CALL X%get(dummyvec)
        IF(.NOT. ((dummyvec(1) .APPROXEQ. 0.46428571428571430_SRK) &
           .AND.  (dummyvec(2) .APPROXEQ. 0.85714285714285721_SRK) &
           .AND.  (dummyvec(3) .APPROXEQ. 0.96428571428571430_SRK) &
           .AND.   thisLS%info == 0)) THEN
          WRITE(*,*) 'CALL Direct%solve() -LU method FAILED!'
          STOP 666
        ENDIF
      ENDSELECT

      
      !Reset X, and solve it again
      SELECTTYPE(X => thisLS%X); TYPE IS(RealVectorType)
        CALL X%set(1.0_SRK)
        CALL thisLS%solve()
        IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
        ALLOCATE(dummyvec(X%n))
        CALL X%get(dummyvec)
        IF(.NOT. ((dummyvec(1) .APPROXEQ. 0.46428571428571430_SRK) &
           .AND.  (dummyvec(2) .APPROXEQ. 0.85714285714285721_SRK) &
           .AND.  (dummyvec(3) .APPROXEQ. 0.96428571428571430_SRK) &
           .AND.   thisLS%isDecomposed)) THEN
          WRITE(*,*) 'CALL Direct%solve() -LU method FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      CALL thisLS%A%clear()
      CALL thisLS%clear()
      
      !Sparse matrix, just make sure it could go to CGNR and could
      ! get result, the details will be tested in CGNR
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',GE)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',SPARSE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',9_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',33_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',9_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',9_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      ! initialize matrix A
      SELECTTYPE(A=>thisLS%A); TYPE IS(SparseMatrixType)
        CALL A%setShape(1,1, 4.0_SRK)
        CALL A%setShape(1,2,-1.0_SRK)
        CALL A%setShape(1,4,-1.0_SRK)
        CALL A%setShape(2,1,-1.0_SRK)
        CALL A%setShape(2,2, 4.0_SRK)
        CALL A%setShape(2,3,-1.0_SRK)
        CALL A%setShape(2,5,-1.0_SRK)
        CALL A%setShape(3,2,-1.0_SRK)
        CALL A%setShape(3,3, 4.0_SRK)
        CALL A%setShape(3,6,-1.0_SRK)
        CALL A%setShape(4,1,-1.0_SRK)
        CALL A%setShape(4,4, 4.0_SRK)
        CALL A%setShape(4,5,-1.0_SRK)
        CALL A%setShape(4,7,-1.0_SRK)
        CALL A%setShape(5,2,-1.0_SRK)
        CALL A%setShape(5,4,-1.0_SRK)
        CALL A%setShape(5,5, 4.0_SRK)
        CALL A%setShape(5,6,-1.0_SRK)
        CALL A%setShape(5,8,-1.0_SRK)
        CALL A%setShape(6,3,-1.0_SRK)
        CALL A%setShape(6,5,-1.0_SRK)
        CALL A%setShape(6,6, 4.0_SRK)
        CALL A%setShape(6,9,-1.0_SRK)
        CALL A%setShape(7,4,-1.0_SRK)
        CALL A%setShape(7,7, 4.0_SRK)
        CALL A%setShape(7,8,-1.0_SRK)
        CALL A%setShape(8,5,-1.0_SRK)
        CALL A%setShape(8,7,-1.0_SRK)
        CALL A%setShape(8,8, 4.0_SRK)
        CALL A%setShape(8,9,-1.0_SRK)
        CALL A%setShape(9,6,-1.0_SRK)
        CALL A%setShape(9,8,-1.0_SRK)
        CALL A%setShape(9,9, 4.0_SRK)
      ENDSELECT
      
      ! initialize vector X
      SELECTTYPE(X => thisLS%X); TYPE IS(RealVectorType)
        CALL X%set(1.0_SRK)
      ENDSELECT
      
      ! initialize vector b
      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1.0_SRK)
      ENDSELECT
      
      !solve it
      CALL thisLS%solve()
      
      IF(thisLS%info /= 0) THEN
        WRITE(*,*) thisLS%info
        WRITE(*,*) 'CALL Direct%solve() -GE method FAILED!'
        STOP 666
      ENDIF
      CALL thisLS%clear()
      
      !DenseRect matrix, just make sure that it could go to CGNR.
      !The result will be checked later.
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',GE)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',DENSERECT)
      CALL pList%add('LinearSolverType->A->MatrixType->n',3_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->m',2_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',2_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',3_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      ! A=[1 1]  b=[1]
      !   [1 2]    [2]
      !   [1 3]    [2]
      ! initialize matrix A
      SELECTTYPE(A => thisLS%A); TYPE IS(DenseRectMatrixType)
        CALL A%set(1,1,1._SRK)
        CALL A%set(1,2,1._SRK)
        CALL A%set(2,1,1._SRK)
        CALL A%set(2,2,2._SRK)
        CALL A%set(3,1,1._SRK)
        CALL A%set(3,2,3._SRK)
      ENDSELECT
      
      ! initialize vector b
      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1,1._SRK)
        CALL b%set(2,2._SRK)
        CALL b%set(3,2._SRK)
      ENDSELECT
      
      ! initialize vector X
      SELECTTYPE(X => thisLS%X); TYPE IS(RealVectorType)
      CALL X%set(1.0_SRK)
      ENDSELECT
      
      ! solve it
      CALL thisLS%solve()
      IF(thisLS%info /= 0) THEN
        WRITE(*,*) 'CALL Direct%solve() -GE method FAILED!'
        STOP 666
      ENDIF
      CALL thisLS%clear()

    !Test LU (Dense square matrix and tridiagonal matrix
    !other matrix will raise error message)
      !Dense square
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',LU)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',DENSESQUARE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',3_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->isSym',.FALSE.)
      CALL pList%add('LinearSolverType->x->VectorType->n',3_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',3_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      !Singular case
      !A=[1 0 1 ]  b=[4]    x=[*]
      !  [2 5 -2]    [6]      [*]
      !  [1 0 1 ]    [4]      [*]
      SELECTTYPE(A => thisLS%A); TYPE IS(DenseSquareMatrixType)
        CALL A%set(1,1,1._SRK)
        CALL A%set(1,2,0._SRK)
        CALL A%set(1,3,1._SRK)
        CALL A%set(2,1,2._SRK)
        CALL A%set(2,2,5._SRK)
        CALL A%set(2,3,-2._SRK)
        CALL A%set(3,1,1._SRK)
        CALL A%set(3,2,0._SRK)
        CALL A%set(3,3,1._SRK)
      ENDSELECT
      
      ! initialize vector b
      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1,4._SRK)
        CALL b%set(2,6._SRK)
        CALL b%set(3,4._SRK)
      ENDSELECT
      
      ! solve it
      CALL thisLS%solve()
      IF(thisLS%info /= -1) THEN
        WRITE(*,*) 'CALL Direct%solve() -LU method FAILED!'
        STOP 666
      ENDIF
      CALL thisLS%clear()
      
      ! Normal Case (non-singular)
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',LU)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',DENSESQUARE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',3_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->isSym',.FALSE.)
      CALL pList%add('LinearSolverType->x->VectorType->n',3_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',3_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      !A=[1 0 1 ]  b=[4]    x=[1]
      !  [2 5 -2]    [6]      [2]
      !  [3 6 9 ]    [42]     [3]
      SELECTTYPE(A => thisLS%A); TYPE IS(DenseSquareMatrixType)
        CALL A%set(1,1,1._SRK)
        CALL A%set(1,2,0._SRK)
        CALL A%set(1,3,1._SRK)
        CALL A%set(2,1,2._SRK)
        CALL A%set(2,2,5._SRK)
        CALL A%set(2,3,-2._SRK)
        CALL A%set(3,1,3._SRK)
        CALL A%set(3,2,6._SRK)
        CALL A%set(3,3,9._SRK)
      ENDSELECT

      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1, 4._SRK)
        CALL b%set(2, 6._SRK)
        CALL b%set(3,42._SRK)
      ENDSELECT

      CALL thisLS%solve()
      
      !Check M
      ! M=[3      6   9]
      !   [1/3    2  -2]
      !   [2/3 -1/2  -9]
      SELECTTYPE(M => thisLS%M); TYPE IS(DenseSquareMatrixType)
        IF(.NOT.((M%A(1,1) .APPROXEQ.   3._SRK) &
           .AND. (M%A(1,2) .APPROXEQ.   6._SRK) &
           .AND. (M%A(1,3) .APPROXEQ.   9._SRK) &
           .AND. (M%A(2,1) .APPROXEQ.   1._SRK/3._SRK) &
           .AND. (M%A(2,2) .APPROXEQ.  -2._SRK) &
           .AND. (M%A(2,3) .APPROXEQ.  -2._SRK) &
           .AND. (M%A(3,1) .APPROXEQ.   2._SRK/3._SRK) &
           .AND. (M%A(3,2) .APPROXEQ. -0.5_SRK) &
           .AND. (M%A(3,3) .APPROXEQ.  -9._SRK) )) THEN
          WRITE(*,*) 'CALL Direct%solve() -LU method FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      
      ! check IPIV: IPIV=[3 3 0]
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Direct)
        IF(    thisLS%IPIV(1) /= 3 &
          .OR. thisLS%IPIV(2) /= 3 &
          .OR. thisLS%IPIV(3) /= 0 ) THEN
          WRITE(*,*) 'CALL Direct%solve() -LU method FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      
      ! check X
      SELECTTYPE(X => thisLS%X); TYPE IS(RealVectorType)
        IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
        ALLOCATE(dummyvec(X%n))
        CALL X%get(dummyvec)
        IF(.NOT.((dummyvec(1) .APPROXEQ. 1._SRK) &
           .AND. (dummyvec(2) .APPROXEQ. 2._SRK) &
           .AND. (dummyvec(3) .APPROXEQ. 3._SRK) &
           .AND. (thisLS%info == 0))) THEN
          WRITE(*,*) 'CALL Direct%solve() -LU method FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      
      ! reset right hand side and solve it again
      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1, 4._SRK)
        CALL b%set(2,14._SRK)
        CALL b%set(3,30._SRK)
      ENDSELECT
      thisLS%isDecomposed=.FALSE.
      CALL thisLS%solve()
      
      ! check X
      SELECTTYPE(X => thisLS%X); TYPE IS(RealVectorType)
        IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
        ALLOCATE(dummyvec(X%n))
        CALL X%get(dummyvec)
        IF(.NOT.((dummyvec(1) .APPROXEQ. 3._SRK) &
           .AND. (dummyvec(2) .APPROXEQ. 2._SRK) &
           .AND. (dummyvec(3) .APPROXEQ. 1._SRK) &
           .AND. (thisLS%info == 0))) THEN
          WRITE(*,*) 'CALL Direct%solve() -LU method FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      CALL thisLS%clear()

      !Tridiagonal (singular case)
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',LU)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',TRIDIAG)
      CALL pList%add('LinearSolverType->A->MatrixType->n',3_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->isSym',.TRUE.)
      CALL pList%add('LinearSolverType->x->VectorType->n',3_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',3_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)

      !A=[ 4 -1  0]
      !  [-1  4 -1]
      !  [ 0  0  0]
      SELECTTYPE(A => thisLS%A); TYPE IS(TriDiagMatrixType)
        CALL A%set(1,1,4._SRK)
        CALL A%set(1,2,-1._SRK)
        CALL A%set(2,2,4._SRK)
        CALL A%set(2,3,0._SRK)
        CALL A%set(3,3,0._SRK)
      ENDSELECT
      
      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1,1._SRK)
        CALL b%set(1,2._SRK)
        CALL b%set(1,3._SRK)
      ENDSELECT
      
      CALL thisLS%solve()
      !Check 
      IF(thisLS%info /= -1) THEN
        WRITE(*,*) 'CALL Direct%solve() -LU method FAILED!'
        STOP 666
      ENDIF
      CALL thisLS%clear()
      
      ! Tridiagonal (non-singular case)
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',LU)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',TRIDIAG)
      CALL pList%add('LinearSolverType->A->MatrixType->n',3_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->isSym',.TRUE.)
      CALL pList%add('LinearSolverType->x->VectorType->n',3_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',3_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      !A=[ 4 -1  0]
      !  [-1  4 -1]
      !  [ 0 -1  4]
      SELECTTYPE(A => thisLS%A); TYPE IS(TriDiagMatrixType)
        CALL A%set(1,1,4._SRK)
        CALL A%set(1,2,-1._SRK)
        CALL A%set(2,2,4._SRK)
        CALL A%set(2,3,-1._SRK)
        CALL A%set(3,3,4._SRK)
      ENDSELECT
      
      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1,1._SRK)
        CALL b%set(2,2._SRK)
        CALL b%set(3,3._SRK)
      ENDSELECT
      
      CALL thisLS%solve()
      
      !Check M
      !M=[ 4   -1    0]
      !  [-.25 3.75 -1]
      !  [ 0   -.267  3.7333]
      SELECTTYPE(M => thisLS%M); TYPE IS(TriDiagMatrixType)
        IF(.NOT.((M%a(1,1) .APPROXEQ.  0.0_SRK)  &
           .AND. (M%a(1,2) .APPROXEQ. -.25_SRK)  &
           .AND. (M%a(1,3) .APPROXEQ. -0.266666666666666666_SRK) &
           .AND. (M%a(2,1) .APPROXEQ. 0.25_SRK)  &
           .AND. (M%a(2,2) .APPROXEQ. (1.0_SRK/3.75_SRK)) &
           .AND. (M%a(2,3) .APPROXEQ. (1.0_SRK/3.7333333333333334_SRK)) &
           .AND. (M%a(3,1) .APPROXEQ. -1.0_SRK)  &
           .AND. (M%a(3,2) .APPROXEQ. -1.0_SRK)  &
           .AND. (M%a(3,3) .APPROXEQ.  0.0_SRK)  &
           .AND. thisLS%isDecomposed)) THEN
          WRITE(*,*) 'CALL Direct%solve() -LU method FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      
      !Check X
      SELECTTYPE(X => thisLS%X); TYPE IS(RealVectorType)
        IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
        ALLOCATE(dummyvec(X%n))
        CALL X%get(dummyvec)
        IF(.NOT. ((dummyvec(1) .APPROXEQ. 0.46428571428571430_SRK) &
           .AND.  (dummyvec(2) .APPROXEQ. 0.85714285714285721_SRK) &
           .AND.  (dummyvec(3) .APPROXEQ. 0.96428571428571430_SRK) &
           .AND.   thisLS%info == 0)) THEN
          WRITE(*,*) 'CALL Direct%solve() -LU method FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      
      !Reset X, and solve it again
      SELECTTYPE(X => thisLS%X); TYPE IS(RealVectorType)
        CALL X%set(1.0_SRK)
        CALL thisLS%solve()
        CALL X%get(dummyvec)
        IF(.NOT. ((dummyvec(1) .APPROXEQ. 0.46428571428571430_SRK) &
           .AND.  (dummyvec(2) .APPROXEQ. 0.85714285714285721_SRK) &
           .AND.  (dummyvec(3) .APPROXEQ. 0.96428571428571430_SRK) &
           .AND.   thisLS%isDecomposed)) THEN
          WRITE(*,*) 'CALL Direct%solve() -LU method FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      CALL thisLS%clear()

      !test solve for TriDiag with a non-diagonally dominant A
      !just make sure we get the output statement.
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',LU)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',TRIDIAG)
      CALL pList%add('LinearSolverType->A->MatrixType->n',3_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->isSym',.TRUE.)
      CALL pList%add('LinearSolverType->x->VectorType->n',3_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',3_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      !A=[ O.5  -1    0]
      !  [-1   0.5   -1]
      !  [ 0    -1  0.5]
      SELECTTYPE(A => thisLS%A); TYPE IS(TriDiagMatrixType)
        CALL A%set(1,1,0.5_SRK)
        CALL A%set(1,2,-1._SRK)
        CALL A%set(2,2,0.5_SRK)
        CALL A%set(2,3,-1._SRK)
        CALL A%set(3,3,0.5_SRK)
      ENDSELECT

      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1,1._SRK)
        CALL b%set(2,2._SRK)
        CALL b%set(3,3._SRK)
      ENDSELECT
            
      CALL thisLS%solve()
      CALL thisLS%clear()
      
      !Sparse matrix, just make sure it could go to CGNR and could
      ! get result, the details will be test in CGNR
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',GE)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',SPARSE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',9_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',33_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',9_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',9_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      SELECTTYPE(A => thisLS%A); TYPE IS(SparseMatrixType)
        CALL A%setShape(1,1, 4.0_SRK)
        CALL A%setShape(1,2,-1.0_SRK)
        CALL A%setShape(1,4,-1.0_SRK)
        CALL A%setShape(2,1,-1.0_SRK)
        CALL A%setShape(2,2, 4.0_SRK)
        CALL A%setShape(2,3,-1.0_SRK)
        CALL A%setShape(2,5,-1.0_SRK)
        CALL A%setShape(3,2,-1.0_SRK)
        CALL A%setShape(3,3, 4.0_SRK)
        CALL A%setShape(3,6,-1.0_SRK)
        CALL A%setShape(4,1,-1.0_SRK)
        CALL A%setShape(4,4, 4.0_SRK)
        CALL A%setShape(4,5,-1.0_SRK)
        CALL A%setShape(4,7,-1.0_SRK)
        CALL A%setShape(5,2,-1.0_SRK)
        CALL A%setShape(5,4,-1.0_SRK)
        CALL A%setShape(5,5, 4.0_SRK)
        CALL A%setShape(5,6,-1.0_SRK)
        CALL A%setShape(5,8,-1.0_SRK)
        CALL A%setShape(6,3,-1.0_SRK)
        CALL A%setShape(6,5,-1.0_SRK)
        CALL A%setShape(6,6, 4.0_SRK)
        CALL A%setShape(6,9,-1.0_SRK)
        CALL A%setShape(7,4,-1.0_SRK)
        CALL A%setShape(7,7, 4.0_SRK)
        CALL A%setShape(7,8,-1.0_SRK)
        CALL A%setShape(8,5,-1.0_SRK)
        CALL A%setShape(8,7,-1.0_SRK)
        CALL A%setShape(8,8, 4.0_SRK)
        CALL A%setShape(8,9,-1.0_SRK)
        CALL A%setShape(9,6,-1.0_SRK)
        CALL A%setShape(9,8,-1.0_SRK)
        CALL A%setShape(9,9, 4.0_SRK)
      ENDSELECT
      
      SELECTTYPE(X => thisLS%X); TYPE IS(RealVectorType)
        CALL X%set(1.0_SRK)
      ENDSELECT
      
      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1.0_SRK)
      ENDSELECT
      
      ! solve it
      CALL thisLS%solve()
      
      IF(thisLS%info /= 0) THEN
        WRITE(*,*) 'CALL Direct%solve() -LU method FAILED!'
        STOP 666
      ENDIF
      CALL thisLS%clear()
      
      !DenseRect matrix, just make sure that it could go to CGNR.
      !The result will be checked later.
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',LU)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',DENSERECT)
      CALL pList%add('LinearSolverType->A->MatrixType->n',3_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->m',2_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',2_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',3_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      ! A=[1 1]  b=[1]
      !   [1 2]    [2]
      !   [1 3]    [2]
      SELECTTYPE(A => thisLS%A); TYPE IS(DenseRectMatrixType)
        CALL A%set(1,1,1._SRK)
        CALL A%set(1,2,1._SRK)
        CALL A%set(2,1,1._SRK)
        CALL A%set(2,2,2._SRK)
        CALL A%set(3,1,1._SRK)
        CALL A%set(3,2,3._SRK)
      ENDSELECT
      
      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1,1._SRK)
        CALL b%set(2,2._SRK)
        CALL b%set(3,2._SRK)
      ENDSELECT

      SELECTTYPE(X => thisLS%X); TYPE IS(RealVectorType)
        CALL X%set(1.0_SRK)
      ENDSELECT
      
      !Solve it
      CALL thisLS%solve()
      IF(thisLS%info /= 0) THEN
        WRITE(*,*) 'CALL Direct%solve() -LU method FAILED!'
        STOP 666
      ENDIF
      WRITE(*,*) '  Passed: CALL Direct%solve()'
    !end test of direct solver
      CALL thisLS%clear()
      DEALLOCATE(thisLS)
      
    ENDSUBROUTINE testDirectSolve
!
!-------------------------------------------------------------------------------
    SUBROUTINE testIterativeOthers()
      CLASS(LinearSolverType_Base),ALLOCATABLE :: thisLS
      REAL(SRK),POINTER :: thisX(:),thisX2(:)
      REAL(SRK),ALLOCATABLE :: resid_soln(:),dummyvec(:)
      TYPE(RealVectorType) :: resid
      INTEGER(SIK) :: i
#ifdef HAVE_PETSC
      PetscReal :: rtol,abstol,dtol
      PetscInt  :: maxits,restart
      PetscErrorCode :: ierr
#endif

      ALLOCATE(LinearSolverType_Iterative :: thisLS)
      
    !Test setX0
    
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',BICGSTAB)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',SPARSE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',2_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',1_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',2_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',2_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      SELECTTYPE(X => thisLS%X); TYPE IS(RealVectorType)
        CALL X%set(1,0._SRK)
        CALL X%set(2,0._SRK)
      ENDSELECT
      
      ! initialize X0
      ALLOCATE(thisX2(3))
      thisX2=(/1._SRK,2._SRK,3._SRK/)
      
      !test case that is expected to work, thisX has already been allocated
      SELECTTYPE(thisLS); TYPE IS (LinearSolverType_Iterative)
        CALL thisLS%setX0(thisX2)
        IF(.NOT. (ALLOCATED(thisLS%X) .AND. ASSOCIATED(thisX2) &
           .AND.  thisLS%hasX0)) THEN
          WRITE(*,*) 'CALL Iterative%setX0(...) FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      DEALLOCATE(thisX2)
      CALL thisLS%clear()
      WRITE(*,*) '  Passed: CALL Iterative%setX0(...)'
  
#ifdef HAVE_PETSC
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',PETSC)
      CALL pList%add('LinearSolverType->solverMethod',BICGSTAB)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',SPARSE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',2_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',1_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',2_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',2_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      ! initialize vector X
      SELECTTYPE(X => thisLS%X); TYPE IS(PETScVectorType)
        CALL X%set(1,0._SRK)
        CALL X%set(2,0._SRK)
      ENDSELECT
      
      ! initialize X0
      ALLOCATE(thisX2(3))
      thisX2=(/1._SRK,2._SRK,3._SRK/)
      
      !test case that is expected to work, thisX has already been allocated
      SELECTTYPE(thisLS); TYPE IS (LinearSolverType_Iterative)
        CALL thisLS%setX0(thisX2)
        IF(.NOT. (ALLOCATED(thisLS%X) .AND. ASSOCIATED(thisX2) &
           .AND.  thisLS%hasX0)) THEN
          WRITE(*,*) 'CALL PETScIterative%setX0(...) FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      DEALLOCATE(thisX2)
      CALL thisLS%clear()
      WRITE(*,*) '  Passed: CALL PETScIterative%setX0(...)'
#endif
      
    !Test setConv
      !Bad input
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',BICGSTAB)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',SPARSE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',2_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',1_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',2_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',2_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      SELECTTYPE(thisLS); TYPE IS (LinearSolverType_Iterative)
        CALL thisLS%setConv(-2,-0.1_SRK,-1,-1)
        CALL thisLS%setConv(-2,1.1_SRK,-1,-1)
        !Check if default value is used
        IF(thisLS%maxIters /= 1000_SIK .OR. thisLS%normType /= 2_SIK &
          .OR. thisLS%convTol /= 0.001_SRK .OR. thisLS%nRestart /= 30_SIK) THEN
          WRITE(*,*) thisLS%maxIters, thisLS%normType, thisLS%convTol, thisLS%nRestart
          WRITE(*,*) 'CALL Iterative%setConv(...) FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      
      !Correct input
      SELECTTYPE(thisLS); TYPE IS (LinearSolverType_Iterative)
        CALL thisLS%setConv(1_SIK,0.01_SRK,100_SIK,10_SIK)
        IF(thisLS%maxIters /= 100_SIK .OR. thisLS%normType /= 1_SIK &
          .OR. thisLS%convTol /= 0.01_SRK .OR. thisLS%nRestart /= 10_SIK) THEN
          WRITE(*,*) 'CALL Iterative%setConv(...) FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      CALL thisLS%clear()
      WRITE(*,*) '  Passed: CALL Iterative%setConv(...)'
 
#ifdef HAVE_PETSC
      !Test setConv
      !Bad input
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',PETSC)
      CALL pList%add('LinearSolverType->solverMethod',GMRES) ! GMRES to test nrestart
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',SPARSE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',2_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',1_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',2_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',2_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      SELECTTYPE(thisLS); TYPE IS (LinearSolverType_Iterative)
        CALL thisLS%setConv(-2,-0.1_SRK,-1,-1)
        CALL thisLS%setConv(-2,1.1_SRK,-1,-1)
        !Check if default value is used
        CALL KSPGetTolerances(thisLS%ksp,rtol,abstol,dtol,maxits,ierr)
        CALL KSPGMRESGetRestart(thisLS%ksp,restart,ierr)
        IF(maxits /= 1000_SIK .OR. rtol /= 0.001_SRK &
           .OR. abstol /= 0.001_SRK .OR. restart /= 30_SIK) THEN
          WRITE(*,*) 'CALL PETScIterative%setConv(...) FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      
      !Correct input
      SELECTTYPE(thisLS); TYPE IS (LinearSolverType_Iterative)
        CALL thisLS%setConv(1_SIK,0.01_SRK,100_SIK,10_SIK)
        CALL KSPGetTolerances(thisLS%ksp,rtol,abstol,dtol,maxits,ierr)
        CALL KSPGMRESGetRestart(thisLS%ksp,restart,ierr)
        IF(maxits /= 100_SIK .OR. rtol /= 0.01_SRK &
           .OR. abstol /= 0.01_SRK .OR. restart /= 10_SIK) THEN
           WRITE(*,*) maxits,rtol,abstol,restart
          WRITE(*,*) 'CALL PETScIterative%setConv(...) FAILED2!'
          STOP 666
        ENDIF
      ENDSELECT
      CALL thisLS%clear()
      WRITE(*,*) '  Passed: CALL PETScIterative%setConv(...)'
#endif
      
    !Test getResidual
      !Bad input
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        CALL thisLS%getResidual(resid)
        
        CALL pList%clear()
        CALL pList%add('LinearSolverType->matType',SPARSE)
        CALL pList%add('LinearSolverType->TPLType',NATIVE)
        CALL pList%add('LinearSolverType->solverMethod',BICGSTAB)
        CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
        CALL pList%add('LinearSolverType->numberOMP',1_SNK)
        CALL pList%add('LinearSolverType->timerName','testTimer')
        CALL pList%add('LinearSolverType->A->MatrixType->matType',SPARSE)
        CALL pList%add('LinearSolverType->A->MatrixType->n',5_SNK)
        CALL pList%add('LinearSolverType->A->MatrixType->nnz',1_SNK)
        CALL pList%add('LinearSolverType->x->VectorType->n',5_SNK)
        CALL pList%add('LinearSolverType->b->VectorType->n',5_SNK)
        CALL pList%validate(pList,optListLS)
        CALL thisLS%init(pList)

        CALL thisLS%getResidual(resid)

        CALL vecPList%set('VectorType->n',5_SNK)
        CALL resid%init(vecPList)
        CALL thisLS%getResidual(resid)
        
        CALL thisLS%clear()
        CALL resid%clear()
        
      ENDSELECT
      
      !Correct input
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',BICGSTAB)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',SPARSE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',9_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',33_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',9_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',9_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)

      !A =  4    -1     0    -1     0     0     0     0     0
      !    -1     4    -1     0    -1     0     0     0     0
      !     0    -1     4     0     0    -1     0     0     0
      !    -1     0     0     4    -1     0    -1     0     0
      !     0    -1     0    -1     4    -1     0    -1     0
      !     0     0    -1     0    -1     4     0     0    -1
      !     0     0     0    -1     0     0     4    -1     0
      !     0     0     0     0    -1     0    -1     4    -1
      !     0     0     0     0     0    -1     0    -1     4
      SELECTTYPE(A => thisLS%A); TYPE IS(SparseMatrixType)
        CALL A%setShape(1,1, 4.0_SRK)
        CALL A%setShape(1,2,-1.0_SRK)
        CALL A%setShape(1,4,-1.0_SRK)
        CALL A%setShape(2,1,-1.0_SRK)
        CALL A%setShape(2,2, 4.0_SRK)
        CALL A%setShape(2,3,-1.0_SRK)
        CALL A%setShape(2,5,-1.0_SRK)
        CALL A%setShape(3,2,-1.0_SRK)
        CALL A%setShape(3,3, 4.0_SRK)
        CALL A%setShape(3,6,-1.0_SRK)
        CALL A%setShape(4,1,-1.0_SRK)
        CALL A%setShape(4,4, 4.0_SRK)
        CALL A%setShape(4,5,-1.0_SRK)
        CALL A%setShape(4,7,-1.0_SRK)
        CALL A%setShape(5,2,-1.0_SRK)
        CALL A%setShape(5,4,-1.0_SRK)
        CALL A%setShape(5,5, 4.0_SRK)
        CALL A%setShape(5,6,-1.0_SRK)
        CALL A%setShape(5,8,-1.0_SRK)
        CALL A%setShape(6,3,-1.0_SRK)
        CALL A%setShape(6,5,-1.0_SRK)
        CALL A%setShape(6,6, 4.0_SRK)
        CALL A%setShape(6,9,-1.0_SRK)
        CALL A%setShape(7,4,-1.0_SRK)
        CALL A%setShape(7,7, 4.0_SRK)
        CALL A%setShape(7,8,-1.0_SRK)
        CALL A%setShape(8,5,-1.0_SRK)
        CALL A%setShape(8,7,-1.0_SRK)
        CALL A%setShape(8,8, 4.0_SRK)
        CALL A%setShape(8,9,-1.0_SRK)
        CALL A%setShape(9,6,-1.0_SRK)
        CALL A%setShape(9,8,-1.0_SRK)
        CALL A%setShape(9,9, 4.0_SRK)
      ENDSELECT

      ! initialize vector X
      SELECTTYPE(X => thisLS%X); TYPE IS(RealVectorType)
        CALL X%set(1._SRK)
      ENDSELECT

      ! initialize vector b
      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(10._SRK)
      ENDSELECT
      
      CALL vecPList%set('VectorType->n',9_SNK)
      CALL resid%init(vecPList)
      ALLOCATE(resid_soln(9))
      resid_soln=(/-8._SRK,-9._SRK,-8._SRK,-9._SRK,-10._SRK, &
        -9._SRK,-8._SRK,-9._SRK,-8._SRK/)
      
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        CALL thisLS%getResidual(resid)
      ENDSELECT
      
      DO i=1,resid%n
        IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
        ALLOCATE(dummyvec(resid%n))
        CALL resid%get(dummyvec)
        IF(.NOT.(dummyvec(i) .APPROXEQ. resid_soln(i))) THEN
          WRITE(*,*) 'CALL Iterative%getResidual(...) FAILED!'
          STOP 666
        ENDIF
      ENDDO

      CALL thisLS%clear()
      DEALLOCATE(resid_soln)      
      DEALLOCATE(thisLS)
      WRITE(*,*) '  Passed: CALL Iterative%getResidual(...)'
      
    ENDSUBROUTINE testIterativeOthers
!
!-------------------------------------------------------------------------------
    SUBROUTINE testIterativeSolve_BICGSTAB()
      CLASS(LinearSolverType_Base),ALLOCATABLE :: thisLS
      REAL(SRK),ALLOCATABLE :: thisB(:),dummyvec(:)
      REAL(SRK),POINTER :: thisX(:)
      INTEGER(SIK) :: i
      LOGICAL(SBK) :: match

      ALLOCATE(LinearSolverType_Iterative :: thisLS)

    !With BiCGSTAB
      !The sparse matrix type
      
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',BICGSTAB)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',SPARSE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',9_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',33_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',9_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',9_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      !A =  4    -1     0    -1     0     0     0     0     0
      !    -1     4    -1     0    -1     0     0     0     0
      !     0    -1     4     0     0    -1     0     0     0
      !    -1     0     0     4    -1     0    -1     0     0
      !     0    -1     0    -1     4    -1     0    -1     0
      !     0     0    -1     0    -1     4     0     0    -1
      !     0     0     0    -1     0     0     4    -1     0
      !     0     0     0     0    -1     0    -1     4    -1
      !     0     0     0     0     0    -1     0    -1     4
      SELECTTYPE(A => thisLS%A); TYPE IS(SparseMatrixType)
          CALL A%setShape(1,1, 4.0_SRK)
          CALL A%setShape(1,2,-1.0_SRK)
          CALL A%setShape(1,4,-1.0_SRK)
          CALL A%setShape(2,1,-1.0_SRK)
          CALL A%setShape(2,2, 4.0_SRK)
          CALL A%setShape(2,3,-1.0_SRK)
          CALL A%setShape(2,5,-1.0_SRK)
          CALL A%setShape(3,2,-1.0_SRK)
          CALL A%setShape(3,3, 4.0_SRK)
          CALL A%setShape(3,6,-1.0_SRK)
          CALL A%setShape(4,1,-1.0_SRK)
          CALL A%setShape(4,4, 4.0_SRK)
          CALL A%setShape(4,5,-1.0_SRK)
          CALL A%setShape(4,7,-1.0_SRK)
          CALL A%setShape(5,2,-1.0_SRK)
          CALL A%setShape(5,4,-1.0_SRK)
          CALL A%setShape(5,5, 4.0_SRK)
          CALL A%setShape(5,6,-1.0_SRK)
          CALL A%setShape(5,8,-1.0_SRK)
          CALL A%setShape(6,3,-1.0_SRK)
          CALL A%setShape(6,5,-1.0_SRK)
          CALL A%setShape(6,6, 4.0_SRK)
          CALL A%setShape(6,9,-1.0_SRK)
          CALL A%setShape(7,4,-1.0_SRK)
          CALL A%setShape(7,7, 4.0_SRK)
          CALL A%setShape(7,8,-1.0_SRK)
          CALL A%setShape(8,5,-1.0_SRK)
          CALL A%setShape(8,7,-1.0_SRK)
          CALL A%setShape(8,8, 4.0_SRK)
          CALL A%setShape(8,9,-1.0_SRK)
          CALL A%setShape(9,6,-1.0_SRK)
          CALL A%setShape(9,8,-1.0_SRK)
          CALL A%setShape(9,9, 4.0_SRK)
      ENDSELECT
      
      ! build X0 and set it to 1.0s
      ALLOCATE(thisX(9))
      thisX=1.0_SRK
      
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
          CALL thisLS%setX0(thisX)
      ENDSELECT

      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1.0_SRK)
      ENDSELECT

      !set iterations and convergence information and build/set M
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        CALL thisLS%setConv(2_SIK,1.0E-9_SRK,1000_SIK,30_SIK)
      ENDSELECT
      
      !solve it
      CALL thisLS%solve()
      
      !Store expected solution (from MATLAB) in B
      ALLOCATE(thisB(9))
      thisB(1)=0.6875_SRK
      thisB(2)=0.875_SRK
      thisB(3)=0.6875_SRK
      thisB(4)=0.875_SRK
      thisB(5)=1.125_SRK
      thisB(6)=0.875_SRK
      thisB(7)=0.6875_SRK
      thisB(8)=0.875_SRK
      thisB(9)=0.6875_SRK
      !multiply by 10,000 so we can match first five places.
      thisB=10000.0_SRK*thisB
      match=.TRUE.
      DO i=1,SIZE(thisB)
        SELECTTYPE(X => thisLS%X); TYPE IS(RealVectorType)
          IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
          ALLOCATE(dummyvec(X%n))
          CALL X%get(dummyvec)
          IF(NINT(thisB(i)) /= NINT(10000.0_SRK*dummyvec(i))) THEN
            match=.FALSE.
            EXIT
          ENDIF
        ENDSELECT
      ENDDO
      IF(.NOT. match) THEN
        WRITE(*,*) 'CALL Iterative%solve() -BICGSTAB FAILED!'
        STOP 666
      ENDIF
      
      DEALLOCATE(thisB)
      CALL thisLS%A%clear()
      CALL thisLS%clear()
      
    !test with A being densesquare
    
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',BICGSTAB)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',DENSESQUARE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',9_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->isSym',.TRUE.)
      CALL pList%add('LinearSolverType->x->VectorType->n',9_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',9_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      DO i=1,9
        SELECTTYPE(A => thisLS%A); TYPE IS(DenseSquareMatrixType)
          CALL A%set(i,i,4.0_SRK)
          IF((i < 9).AND.((i /= 3).AND.(i /= 6))) THEN
            CALL A%set(i,i+1,-1.0_SRK)
          ENDIF
          IF(i < 7) THEN
            CALL A%set(i,i+3,-1.0_SRK)
          ENDIF
        ENDSELECT
      ENDDO
      
      !build X0 and set it to 1.0s
      ALLOCATE(thisX(9))
      thisX=1.0_SRK
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        CALL thisLS%setX0(thisX)
      ENDSELECT
      
      !set b
      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1.0_SRK)
      ENDSELECT
      
      !set iterations and convergence information and build/set M
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        CALL thisLS%setConv(2_SIK,1.0E-9_SRK,1000_SIK,30_SIK)
      ENDSELECT
      
      !solve
      CALL thisLS%solve()
      ALLOCATE(thisB(9))
      thisB(1)=0.6875_SRK
      thisB(2)=0.875_SRK
      thisB(3)=0.6875_SRK
      thisB(4)=0.875_SRK
      thisB(5)=1.125_SRK
      thisB(6)=0.875_SRK
      thisB(7)=0.6875_SRK
      thisB(8)=0.875_SRK
      thisB(9)=0.6875_SRK
      thisB=thisB*10000._SRK
      match=.TRUE.
      DO i=1,SIZE(thisB)
        SELECTTYPE(X => thisLS%X); TYPE IS(RealVectorType)
          IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
          ALLOCATE(dummyvec(X%n))
          CALL X%get(dummyvec)
          IF(NINT(thisB(i)) /= NINT(10000.0_SRK*dummyvec(i))) THEN
            match=.FALSE.
            EXIT
          ENDIF
        ENDSELECT
      ENDDO
      IF(.NOT. match) THEN
        WRITE(*,*) 'CALL Iterative%solve() - BiCGSTAB FAILED!'
        STOP 666
      ENDIF
      !test to see how it performs with an already decomposed M
      !reset X to 1.0s
      match=.TRUE.
      DO i=1,SIZE(thisB)
        SELECTTYPE(X => thisLS%X); TYPE IS(RealVectorType)
          CALL X%set(1.0_SRK)
          CALL thisLS%solve()
          IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
          ALLOCATE(dummyvec(X%n))
          CALL X%get(dummyvec)
          IF(NINT(thisB(i)) /= NINT(10000.0_SRK*dummyvec(i))) THEN
            match=.FALSE.
            EXIT
          ENDIF
        ENDSELECT
      ENDDO
      IF(.NOT. match) THEN
        WRITE(*,*) 'CALL Iterative%solve() - BiCGSTAB FAILED!'
        STOP 666
      ENDIF
      CALL thisLS%clear()
      DEALLOCATE(thisB)
      
      ! TriDiagonal matrix, it will go to LU method
      
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',BICGSTAB)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',TRIDIAG)
      CALL pList%add('LinearSolverType->A->MatrixType->n',3_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->isSym',.TRUE.)
      CALL pList%add('LinearSolverType->x->VectorType->n',3_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',3_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      !A=[ 4 -1  0]
      !  [-1  4 -1]
      !  [ 0 -1  4]
      SELECTTYPE(A => thisLS%A); TYPE IS(TriDiagMatrixType)
        CALL A%set(1,1, 4._SRK)
        CALL A%set(1,2,-1._SRK)
        CALL A%set(2,2, 4._SRK)
        CALL A%set(2,3,-1._SRK)
        CALL A%set(3,3, 4._SRK)
      ENDSELECT
      
      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1,1._SRK)
        CALL b%set(2,2._SRK)
        CALL b%set(3,3._SRK)
      ENDSELECT

      CALL thisLS%solve()
      
      IF(thisLS%info /= 0) THEN
        WRITE(*,*) 'CALL Iterative%solve() - BiCGSTAB FAILED!'
        STOP 666
      ENDIF
      CALL thisLS%clear()
      DEALLOCATE(thisX)

      !DenseRect matrix
      
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',BICGSTAB)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',DENSERECT)
      CALL pList%add('LinearSolverType->A->MatrixType->n',3_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->m',2_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',2_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',3_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      ! A=[1 1]  b=[1]
      !   [1 2]    [2]
      !   [1 3]    [2]
      SELECTTYPE(A => thisLS%A); TYPE IS(DenseRectMatrixType)
        CALL A%set(1,1,1._SRK)
        CALL A%set(1,2,1._SRK)
        CALL A%set(2,1,1._SRK)
        CALL A%set(2,2,2._SRK)
        CALL A%set(3,1,1._SRK)
        CALL A%set(3,2,3._SRK)
      ENDSELECT
      
      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1,1._SRK)
        CALL b%set(2,2._SRK)
        CALL b%set(2,2._SRK)
      ENDSELECT
      
      ! initialize X0
      ALLOCATE(thisX(2))
      thisX=1.0_SRK

      !set iterations and convergence information and build/set M
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        CALL thisLS%setX0(thisX)
        CALL thisLS%setConv(2_SIK,1.0E-9_SRK,1000_SIK,30_SIK)
      ENDSELECT
      
      !Solve it
      CALL thisLS%solve()
      
      IF(thisLS%info /= 0) THEN
        WRITE(*,*) 'CALL Iterative%solve() -BiCGSTAB method FAILED!'
        STOP 666
      ENDIF
      
      DEALLOCATE(thisX)
      CALL thisLS%clear()
      WRITE(*,*)'  Passed: CALL Iterative%solver() BICGSTAB'

#ifdef HAVE_PETSC
      !The PETSC sparse matrix type
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',PETSC)
      CALL pList%add('LinearSolverType->solverMethod',BICGSTAB)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',SPARSE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',9_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->isSym',.FALSE.)
      CALL pList%add('LinearSolverType->x->VectorType->n',9_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',9_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      !A =  4    -1     0    -1     0     0     0     0     0
      !    -1     4    -1     0    -1     0     0     0     0
      !     0    -1     4     0     0    -1     0     0     0
      !    -1     0     0     4    -1     0    -1     0     0
      !     0    -1     0    -1     4    -1     0    -1     0
      !     0     0    -1     0    -1     4     0     0    -1
      !     0     0     0    -1     0     0     4    -1     0
      !     0     0     0     0    -1     0    -1     4    -1
      !     0     0     0     0     0    -1     0    -1     4
      SELECTTYPE(A => thisLS%A); TYPE IS(PETScMatrixType)
          CALL A%set(1,1, 4.0_SRK)
          CALL A%set(1,2,-1.0_SRK)
          CALL A%set(1,4,-1.0_SRK)
          CALL A%set(2,1,-1.0_SRK)
          CALL A%set(2,2, 4.0_SRK)
          CALL A%set(2,3,-1.0_SRK)
          CALL A%set(2,5,-1.0_SRK)
          CALL A%set(3,2,-1.0_SRK)
          CALL A%set(3,3, 4.0_SRK)
          CALL A%set(3,6,-1.0_SRK)
          CALL A%set(4,1,-1.0_SRK)
          CALL A%set(4,4, 4.0_SRK)
          CALL A%set(4,5,-1.0_SRK)
          CALL A%set(4,7,-1.0_SRK)
          CALL A%set(5,2,-1.0_SRK)
          CALL A%set(5,4,-1.0_SRK)
          CALL A%set(5,5, 4.0_SRK)
          CALL A%set(5,6,-1.0_SRK)
          CALL A%set(5,8,-1.0_SRK)
          CALL A%set(6,3,-1.0_SRK)
          CALL A%set(6,5,-1.0_SRK)
          CALL A%set(6,6, 4.0_SRK)
          CALL A%set(6,9,-1.0_SRK)
          CALL A%set(7,4,-1.0_SRK)
          CALL A%set(7,7, 4.0_SRK)
          CALL A%set(7,8,-1.0_SRK)
          CALL A%set(8,5,-1.0_SRK)
          CALL A%set(8,7,-1.0_SRK)
          CALL A%set(8,8, 4.0_SRK)
          CALL A%set(8,9,-1.0_SRK)
          CALL A%set(9,6,-1.0_SRK)
          CALL A%set(9,8,-1.0_SRK)
          CALL A%set(9,9, 4.0_SRK)
      ENDSELECT
      
      ! build X0 and set it to 1.0s
      ALLOCATE(thisX(9))
      thisX=1.0_SRK
      
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        CALL thisLS%setX0(thisX)
      ENDSELECT

      SELECTTYPE(b => thisLS%b); TYPE IS(PETScVectorType)
        CALL b%set(1.0_SRK)
      ENDSELECT

      !set iterations and convergence information and build/set M
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        CALL thisLS%setConv(2_SIK,1.0E-9_SRK,1000_SIK,30_SIK)
      ENDSELECT
      
      !solve it
      CALL thisLS%solve()
      
      !Store expected solution (from MATLAB) in B
      ALLOCATE(thisB(9))
      thisB(1)=0.6875_SRK
      thisB(2)=0.875_SRK
      thisB(3)=0.6875_SRK
      thisB(4)=0.875_SRK
      thisB(5)=1.125_SRK
      thisB(6)=0.875_SRK
      thisB(7)=0.6875_SRK
      thisB(8)=0.875_SRK
      thisB(9)=0.6875_SRK
      !multiply by 10,000 so we can match first five places.
      thisB=10000.0_SRK*thisB
      match=.TRUE.
      DO i=1,SIZE(thisB)
        SELECTTYPE(X => thisLS%X); TYPE IS(PETScVectorType)
          IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
          ALLOCATE(dummyvec(X%n))
          CALL X%get(dummyvec)
          IF(NINT(thisB(i)) /= NINT(10000.0_SRK*dummyvec(i))) THEN
            match=.FALSE.
            EXIT
          ENDIF
        ENDSELECT
      ENDDO
      IF(.NOT. match) THEN
        WRITE(*,*) 'CALL PETScIterative%solve() -BICGSTAB FAILED!'
        STOP 666
      ENDIF
      
      DEALLOCATE(thisB)
      CALL thisLS%A%clear()
      CALL thisLS%clear()
      
      !test with A being densesquare
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',PETSC)
      CALL pList%add('LinearSolverType->solverMethod',BICGSTAB)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',DENSESQUARE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',9_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->isSym',.TRUE.)
      CALL pList%add('LinearSolverType->x->VectorType->n',9_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',9_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      DO i=1,9
        SELECTTYPE(A => thisLS%A); TYPE IS(PETScMatrixType)
          CALL A%set(i,i,4.0_SRK)
          IF((i < 9).AND.((i /= 3).AND.(i /= 6))) THEN
            CALL A%set(i,i+1,-1.0_SRK)
          ENDIF
          IF(i < 7) THEN
            CALL A%set(i,i+3,-1.0_SRK)
          ENDIF
        ENDSELECT
      ENDDO
      
      !build X0 and set it to 1.0s
      ALLOCATE(thisX(9))
      thisX=1.0_SRK
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        CALL thisLS%setX0(thisX)
      ENDSELECT
      
      SELECTTYPE(b => thisLS%b); TYPE IS(PETScVectorType)
        CALL b%set(1.0_SRK)
      ENDSELECT
      
      !set iterations and convergence information and build/set M
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        CALL thisLS%setConv(2_SIK,1.0E-9_SRK,1000_SIK,30_SIK)
      ENDSELECT
      
      !solve
      CALL thisLS%solve()
      
      ALLOCATE(thisB(9))
      thisB(1)=0.6875_SRK
      thisB(2)=0.875_SRK
      thisB(3)=0.6875_SRK
      thisB(4)=0.875_SRK
      thisB(5)=1.125_SRK
      thisB(6)=0.875_SRK
      thisB(7)=0.6875_SRK
      thisB(8)=0.875_SRK
      thisB(9)=0.6875_SRK
      thisB=thisB*10000._SRK
      match=.TRUE.
      DO i=1,SIZE(thisB)
        SELECTTYPE(X => thisLS%X); TYPE IS(RealVectorType)
          IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
          ALLOCATE(dummyvec(X%n))
          CALL X%get(dummyvec)
          IF(NINT(thisB(i)) /= NINT(10000.0_SRK*dummyvec(i))) THEN
            match=.FALSE.
            EXIT
          ENDIF
        ENDSELECT
      ENDDO
      IF(.NOT. match) THEN
        WRITE(*,*) 'CALL PETSCIterative%solve() -BiCGSTAB FAILED!'
        STOP 666
      ENDIF
      !test to see how it performs with an already decomposed M
      !reset X to 1.0s
      match=.TRUE.
      DO i=1,SIZE(thisB)
        SELECTTYPE(X => thisLS%X); TYPE IS(RealVectorType)
          CALL X%set(1.0_SRK)
          CALL thisLS%solve()
          IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
          ALLOCATE(dummyvec(X%n))
          CALL X%get(dummyvec)
          IF(NINT(thisB(i)) /= NINT(10000.0_SRK*dummyvec(i))) THEN
            match=.FALSE.
            EXIT
          ENDIF
        ENDSELECT
      ENDDO
      IF(.NOT. match) THEN
        WRITE(*,*) 'CALL PETSCIterative%solve() -BiCGSTAB FAILED!'
        STOP 666
      ENDIF
      WRITE(*,*)'  Passed: CALL PETScIterative%solver() BICGSTAB'
      CALL thisLS%clear()
      DEALLOCATE(thisB)
#endif
 
      DEALLOCATE(thisLS)

    ENDSUBROUTINE testIterativeSolve_BICGSTAB
!
!-------------------------------------------------------------------------------
    SUBROUTINE testIterativeSolve_CGNR()
      CLASS(LinearSolverType_Base),ALLOCATABLE :: thisLS
      REAL(SRK),ALLOCATABLE :: thisB(:),dummyvec(:)
      REAL(SRK),POINTER :: thisX(:)
      INTEGER(SIK) :: i
      LOGICAL(SBK) :: match

      ALLOCATE(LinearSolverType_Iterative :: thisLS)

      !test CGNR
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',CGNR)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',DENSERECT)
      CALL pList%add('LinearSolverType->A->MatrixType->n',2_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->m',3_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',3_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',2_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      !underdetermined matrix: solver%info should return -1.
      ! A=[1 1 1]  b= [1]
      !   [1 2 3]     [3]
      SELECTTYPE(A => thisLS%A); TYPE IS(DenseRectMatrixType)
        CALL A%set(1,1,1._SRK)
        CALL A%set(1,2,1._SRK)
        CALL A%set(1,3,1._SRK)
        CALL A%set(2,1,1._SRK)
        CALL A%set(2,2,2._SRK)
        CALL A%set(2,3,3._SRK)
      ENDSELECT

      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1,1._SRK)
        CALL b%set(2,3._SRK)
      ENDSELECT
      
      ! initialize X0
      ALLOCATE(thisX(3))
      thisX=1.0_SRK
      
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        CALL thisLS%setX0(thisX)
      ENDSELECT
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
          CALL thisLS%setConv(2_SIK,1.0E-9_SRK,1000_SIK,30_SIK)
      ENDSELECT
      
      CALL thisLS%solve()

      IF(thisLS%info /= -1) THEN
        WRITE(*,*)'CALL Iterative%solve() -CGNR FAILED!'
        STOP 666
      ENDIF
      
      DEALLOCATE(thisX)
      CALL thisLS%clear()
      
      !normal or overdetermined matrix should give answers
      
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',CGNR)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',DENSERECT)
      CALL pList%add('LinearSolverType->A->MatrixType->n',3_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->m',2_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',2_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',3_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      ! A=[1 1]  b=[1]
      !   [1 2]    [2]
      !   [1 3]    [2]
      SELECTTYPE(A => thisLS%A); TYPE IS(DenseRectMatrixType)
        CALL A%set(1,1,1._SRK)
        CALL A%set(1,2,1._SRK)
        CALL A%set(2,1,1._SRK)
        CALL A%set(2,2,2._SRK)
        CALL A%set(3,1,1._SRK)
        CALL A%set(3,2,3._SRK)
      ENDSELECT
      
      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1,1._SRK)
        CALL b%set(2,2._SRK)
        CALL b%set(3,2._SRK)
      ENDSELECT
 
      ! initialize X0
      ALLOCATE(thisX(2))
      thisX=1.0_SRK
      
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        CALL thisLS%setX0(thisX)
      ENDSELECT
      !set iterations and convergence information and
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
          CALL thisLS%setConv(2_SIK,1.0E-13_SRK,1000_SIK,30_SIK)
      ENDSELECT
      
      CALL thisLS%solve()
      
      SELECTTYPE(X => thisLS%X); TYPE IS(RealVectorType)
        IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
        ALLOCATE(dummyvec(X%n))
        CALL X%get(dummyvec)
        IF(.NOT.(SOFTEQ(dummyvec(1),2._SRK/3._SRK,1.0E-13_SRK) &
           .AND. SOFTEQ(dummyvec(2),0.5_SRK,1.0E-13_SRK))) THEN
          WRITE(*,*)'CALL Iterative%solve() FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      
      DEALLOCATE(thisX)
      CALL thisLS%clear()
      
      !DenseSquare matrix
      
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',CGNR)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',DENSESQUARE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',3_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->isSym',.TRUE.)
      CALL pList%add('LinearSolverType->x->VectorType->n',3_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',3_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      !A=[ 4 -1  0]
      !  [-1  4 -1]
      !  [ 0 -1  4]
      SELECTTYPE(A => thisLS%A); TYPE IS(DenseSquareMatrixType)
        CALL A%set(1,1,4._SRK)
        CALL A%set(1,2,-1._SRK)
        CALL A%set(2,2,4._SRK)
        CALL A%set(2,3,-1._SRK)
        CALL A%set(3,3,4._SRK)
      ENDSELECT

      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1,1._SRK)
        CALL b%set(2,2._SRK)
        CALL b%set(3,3._SRK)
      ENDSELECT

      !set iterations and convergence information and
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
          CALL thisLS%setConv(2_SIK,1.0E-13_SRK,1000_SIK,30_SIK)
      ENDSELECT
      
      CALL thisLS%solve()
      
      !Check X
      SELECTTYPE(X => thisLS%X); TYPE IS(RealVectorType)
        IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
        ALLOCATE(dummyvec(X%n))
        CALL X%get(dummyvec)
        IF(.NOT. ((dummyvec(1) .APPROXEQ. 0.46428571428571430_SRK) &
           .AND.  (dummyvec(2) .APPROXEQ. 0.85714285714285721_SRK) &
           .AND.  (dummyvec(3) .APPROXEQ. 0.96428571428571430_SRK) &
           .AND.   thisLS%info == 0)) THEN
          WRITE(*,*) 'CALL Iterative%solve() -CGNR method FAILED!'
          STOP 666
        ENDIF
      ENDSELECT

      CALL thisLS%clear()

      !Sparse matrix
      
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',CGNR)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',SPARSE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',3_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',7_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',3_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',3_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)

      !A=[ 4 -1  0]
      !  [-1  4 -1]
      !  [ 0 -1  4]
      SELECTTYPE(A => thisLS%A); TYPE IS(SparseMatrixType)
        CALL A%setShape(1,1, 4._SRK)
        CALL A%setShape(1,2,-1._SRK)
        CALL A%setShape(2,1,-1._SRK)
        CALL A%setShape(2,2, 4._SRK)
        CALL A%setShape(2,3,-1._SRK)
        CALL A%setShape(3,2,-1._SRK)
        CALL A%setShape(3,3, 4._SRK)
      ENDSELECT
     
      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1,1._SRK)
        CALL b%set(2,2._SRK)
        CALL b%set(3,3._SRK)
      ENDSELECT

      ! initialize X0
      ALLOCATE(thisX(3))
      thisX=0._SRK
      
      !set iterations and convergence information and
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        CALL thisLS%setX0(thisX)
        CALL thisLS%setConv(2_SIK,1.0E-13_SRK,1000_SIK,30_SIK)
      ENDSELECT
      
      !Check X
      CALL thisLS%solve()
      SELECTTYPE(X => thisLS%X); TYPE IS(RealVectorType)
        IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
        ALLOCATE(dummyvec(X%n))
        CALL X%get(dummyvec)
        IF(.NOT. ((dummyvec(1) .APPROXEQ. 0.46428571428571430_SRK) &
           .AND.  (dummyvec(2) .APPROXEQ. 0.85714285714285721_SRK) &
           .AND.  (dummyvec(3) .APPROXEQ. 0.96428571428571430_SRK) &
           .AND.   thisLS%info == 0)) THEN
          WRITE(*,*) 'CALL Iterative%solve() -CGNR method FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      
      DEALLOCATE(thisX)
      CALL thisLS%clear()
      
      !TriDiagonal matrix, it will go to LU
      
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',CGNR)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',TRIDIAG)
      CALL pList%add('LinearSolverType->A->MatrixType->n',3_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->isSym',.TRUE.)
      CALL pList%add('LinearSolverType->x->VectorType->n',3_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',3_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)

      !A=[ 4 -1  0]
      !  [-1  4 -1]
      !  [ 0 -1  4]
      SELECTTYPE(A => thisLS%A); TYPE IS(TriDiagMatrixType)
        CALL A%set(1,1,4._SRK)
        CALL A%set(1,2,-1._SRK)
        CALL A%set(2,2,4._SRK)
        CALL A%set(2,3,-1._SRK)
        CALL A%set(3,3,4._SRK)
      ENDSELECT

      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1,1._SRK)
        CALL b%set(2,2._SRK)
        CALL b%set(3,3._SRK)
      ENDSELECT

      CALL thisLS%solve()

      IF(thisLS%info /= 0) THEN
        WRITE(*,*) 'CALL Iterative%solve() - BiCGSTAB FAILED!'
        STOP 666
      ENDIF
      
      CALL thisLS%clear()

      WRITE(*,*)'  Passed: CALL Iterative%solver() CGNR'
      
#ifdef HAVE_PETSC
      ! DenseSquare matrix
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',PETSC)
      CALL pList%add('LinearSolverType->solverMethod',CGNR)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',DENSESQUARE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',3_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->isSym',.TRUE.)
      CALL pList%add('LinearSolverType->x->VectorType->n',3_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',3_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      !A=[ 4 -1  0]
      !  [-1  4 -1]
      !  [ 0 -1  4]
      SELECTTYPE(A => thisLS%A); TYPE IS(PETScMatrixType)
        CALL A%set(1,1,4._SRK)
        CALL A%set(1,2,-1._SRK)
        CALL A%set(2,2,4._SRK)
        CALL A%set(2,3,-1._SRK)
        CALL A%set(3,3,4._SRK)
      ENDSELECT

      SELECTTYPE(b => thisLS%b); TYPE IS(PETScVectorType)
        CALL b%set(1,1._SRK)
        CALL b%set(2,2._SRK)
        CALL b%set(3,3._SRK)
      ENDSELECT

      !set iterations and convergence information and
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
          CALL thisLS%setConv(2_SIK,1.0E-13_SRK,1000_SIK,30_SIK)
      ENDSELECT
      
      CALL thisLS%solve()
      
      !Check X
      SELECTTYPE(X => thisLS%X); TYPE IS(PETScVectorType)
        IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
        ALLOCATE(dummyvec(X%n))
        CALL X%get(dummyvec)
        IF(.NOT. ((dummyvec(1) .APPROXEQ. 0.46428571428571430_SRK) &
           .AND.  (dummyvec(2) .APPROXEQ. 0.85714285714285721_SRK) &
           .AND.  (dummyvec(3) .APPROXEQ. 0.96428571428571430_SRK) &
           .AND.   thisLS%info == 0)) THEN
          WRITE(*,*) 'CALL PETScIterative%solve() -CGNR method FAILED!'
          STOP 666
        ENDIF
      ENDSELECT

      CALL thisLS%clear()

      ! Sparse matrix
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',PETSC)
      CALL pList%add('LinearSolverType->solverMethod',CGNR)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',SPARSE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',3_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->isSym',.TRUE.)
      CALL pList%add('LinearSolverType->x->VectorType->n',3_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',3_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)

      !A=[ 4 -1  0]
      !  [-1  4 -1]
      !  [ 0 -1  4]
      SELECTTYPE(A => thisLS%A); TYPE IS(PETScMatrixType)
        CALL A%set(1,1, 4._SRK)
        CALL A%set(1,2,-1._SRK)
        CALL A%set(2,1,-1._SRK)
        CALL A%set(2,2, 4._SRK)
        CALL A%set(2,3,-1._SRK)
        CALL A%set(3,2,-1._SRK)
        CALL A%set(3,3, 4._SRK)
      ENDSELECT
     
      SELECTTYPE(b => thisLS%b); TYPE IS(PETScVectorType)
        CALL b%set(1,1._SRK)
        CALL b%set(2,2._SRK)
        CALL b%set(3,3._SRK)
      ENDSELECT

      ! initialize X0
      ALLOCATE(thisX(3))
      thisX=0._SRK
      
      !set iterations and convergence information and
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        CALL thisLS%setX0(thisX)
        CALL thisLS%setConv(2_SIK,1.0E-13_SRK,1000_SIK,30_SIK)
      ENDSELECT
      
      !Check X
      CALL thisLS%solve()
      SELECTTYPE(X => thisLS%X); TYPE IS(PETScVectorType)
        IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
        ALLOCATE(dummyvec(X%n))
        CALL X%get(dummyvec)
        IF(.NOT. ((dummyvec(1) .APPROXEQ. 0.46428571428571430_SRK) &
           .AND.  (dummyvec(2) .APPROXEQ. 0.85714285714285721_SRK) &
           .AND.  (dummyvec(3) .APPROXEQ. 0.96428571428571430_SRK) &
           .AND.   thisLS%info == 0)) THEN
          WRITE(*,*) 'CALL PETScIterative%solve() -CGNR method FAILED!'
          STOP 666
        ENDIF
      ENDSELECT
      
      DEALLOCATE(thisX)
      CALL thisLS%clear()

      WRITE(*,*)'  Passed: CALL PETScIterative%solver() CGNR'
#endif

      DEALLOCATE(thisLS)

    ENDSUBROUTINE testIterativeSolve_CGNR
!
!-------------------------------------------------------------------------------
    SUBROUTINE testIterativeSolve_GMRES()
      CLASS(LinearSolverType_Base),ALLOCATABLE :: thisLS
      REAL(SRK),ALLOCATABLE :: thisB(:),dummyvec(:)
      REAL(SRK),POINTER :: thisX(:)
      INTEGER(SIK) :: i
      LOGICAL(SBK) :: match

      ALLOCATE(LinearSolverType_Iterative :: thisLS)

      !With GMRES
      !The sparse matrix type
      
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',GMRES)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',SPARSE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',9_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->nnz',33_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',9_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',9_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      !A =  4    -1     0    -1     0     0     0     0     0
      !    -1     4    -1     0    -1     0     0     0     0
      !     0    -1     4     0     0    -1     0     0     0
      !    -1     0     0     4    -1     0    -1     0     0
      !     0    -1     0    -1     4    -1     0    -1     0
      !     0     0    -1     0    -1     4     0     0    -1
      !     0     0     0    -1     0     0     4    -1     0
      !     0     0     0     0    -1     0    -1     4    -1
      !     0     0     0     0     0    -1     0    -1     4
      SELECTTYPE(A => thisLS%A); TYPE IS(SparseMatrixType)
          CALL A%setShape(1,1, 4.0_SRK)
          CALL A%setShape(1,2,-1.0_SRK)
          CALL A%setShape(1,4,-1.0_SRK)
          CALL A%setShape(2,1,-1.0_SRK)
          CALL A%setShape(2,2, 4.0_SRK)
          CALL A%setShape(2,3,-1.0_SRK)
          CALL A%setShape(2,5,-1.0_SRK)
          CALL A%setShape(3,2,-1.0_SRK)
          CALL A%setShape(3,3, 4.0_SRK)
          CALL A%setShape(3,6,-1.0_SRK)
          CALL A%setShape(4,1,-1.0_SRK)
          CALL A%setShape(4,4, 4.0_SRK)
          CALL A%setShape(4,5,-1.0_SRK)
          CALL A%setShape(4,7,-1.0_SRK)
          CALL A%setShape(5,2,-1.0_SRK)
          CALL A%setShape(5,4,-1.0_SRK)
          CALL A%setShape(5,5, 4.0_SRK)
          CALL A%setShape(5,6,-1.0_SRK)
          CALL A%setShape(5,8,-1.0_SRK)
          CALL A%setShape(6,3,-1.0_SRK)
          CALL A%setShape(6,5,-1.0_SRK)
          CALL A%setShape(6,6, 4.0_SRK)
          CALL A%setShape(6,9,-1.0_SRK)
          CALL A%setShape(7,4,-1.0_SRK)
          CALL A%setShape(7,7, 4.0_SRK)
          CALL A%setShape(7,8,-1.0_SRK)
          CALL A%setShape(8,5,-1.0_SRK)
          CALL A%setShape(8,7,-1.0_SRK)
          CALL A%setShape(8,8, 4.0_SRK)
          CALL A%setShape(8,9,-1.0_SRK)
          CALL A%setShape(9,6,-1.0_SRK)
          CALL A%setShape(9,8,-1.0_SRK)
          CALL A%setShape(9,9, 4.0_SRK)
      ENDSELECT
      
      ! build X0 and set it to 1.0s
      ALLOCATE(thisX(9))
      thisX=1.0_SRK
      
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
          CALL thisLS%setX0(thisX)
      ENDSELECT

      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1.0_SRK)
      ENDSELECT

      !set iterations and convergence information and build/set M
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        CALL thisLS%setConv(2_SIK,1.0E-9_SRK,1000_SIK,30_SIK)
      ENDSELECT
      
      !solve it
      CALL thisLS%solve()
      
      !Store expected solution (from MATLAB) in B
      ALLOCATE(thisB(9))
      thisB(1)=0.6875_SRK
      thisB(2)=0.875_SRK
      thisB(3)=0.6875_SRK
      thisB(4)=0.875_SRK
      thisB(5)=1.125_SRK
      thisB(6)=0.875_SRK
      thisB(7)=0.6875_SRK
      thisB(8)=0.875_SRK
      thisB(9)=0.6875_SRK
      !multiply by 10,000 so we can match first five places.
      thisB=10000.0_SRK*thisB
      match=.TRUE.
      DO i=1,SIZE(thisB)
        SELECTTYPE(X => thisLS%X); TYPE IS(RealVectorType)
          IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
          ALLOCATE(dummyvec(X%n))
          CALL X%get(dummyvec)
          IF(NINT(thisB(i)) /= NINT(10000.0_SRK*dummyvec(i))) THEN
            match=.FALSE.
            EXIT
          ENDIF
        ENDSELECT
      ENDDO
      IF(.NOT. match) THEN
        WRITE(*,*) 'CALL Iterative%solve() -GMRES FAILED!'
        STOP 666
      ENDIF
      
      DEALLOCATE(thisB)
      CALL thisLS%clear()
      
    !test with A being densesquare
    
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',GMRES)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',DENSESQUARE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',9_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->isSym',.TRUE.)
      CALL pList%add('LinearSolverType->x->VectorType->n',9_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',9_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      DO i=1,9
        SELECTTYPE(A => thisLS%A); TYPE IS(DenseSquareMatrixType)
          CALL A%set(i,i,4.0_SRK)
          IF((i < 9).AND.((i /= 3).AND.(i /= 6))) THEN
            CALL A%set(i,i+1,-1.0_SRK)
          ENDIF
          IF(i < 7) THEN
            CALL A%set(i,i+3,-1.0_SRK)
          ENDIF
        ENDSELECT
      ENDDO
      
      !build X0 and set it to 1.0s
      ALLOCATE(thisX(9))
      thisX=1.0_SRK
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        CALL thisLS%setX0(thisX)
      ENDSELECT
      
      !set b
      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1.0_SRK)
      ENDSELECT
      
      !set iterations and convergence information and build/set M
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        CALL thisLS%setConv(2_SIK,1.0E-9_SRK,1000_SIK,30_SIK)
      ENDSELECT
      
      !solve
      CALL thisLS%solve()
      ALLOCATE(thisB(9))
      thisB(1)=0.6875_SRK
      thisB(2)=0.875_SRK
      thisB(3)=0.6875_SRK
      thisB(4)=0.875_SRK
      thisB(5)=1.125_SRK
      thisB(6)=0.875_SRK
      thisB(7)=0.6875_SRK
      thisB(8)=0.875_SRK
      thisB(9)=0.6875_SRK
      thisB=thisB*10000._SRK
      match=.TRUE.
      DO i=1,SIZE(thisB)
        SELECTTYPE(X => thisLS%X); TYPE IS(RealVectorType)
          IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
          ALLOCATE(dummyvec(X%n))
          CALL X%get(dummyvec)
          IF(NINT(thisB(i)) /= NINT(10000.0_SRK*dummyvec(i))) THEN
            match=.FALSE.
            EXIT
          ENDIF
        ENDSELECT
      ENDDO
      IF(.NOT. match) THEN
        WRITE(*,*) 'CALL Iterative%solve() - GMRES FAILED!'
        STOP 666
      ENDIF
      !test to see how it performs with an already decomposed M
      !reset X to 1.0s
      match=.TRUE.
      DO i=1,SIZE(thisB)
        SELECTTYPE(X => thisLS%X); TYPE IS(RealVectorType)
          CALL X%set(1.0_SRK)
          CALL thisLS%solve()
          IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
          ALLOCATE(dummyvec(X%n))
          CALL X%get(dummyvec)
          IF(NINT(thisB(i)) /= NINT(10000.0_SRK*dummyvec(i))) THEN
            match=.FALSE.
            EXIT
          ENDIF
        ENDSELECT
      ENDDO
      IF(.NOT. match) THEN
        WRITE(*,*) 'CALL Iterative%solve() - GMRES FAILED!'
        STOP 666
      ENDIF
      CALL thisLS%A%clear()
      CALL thisLS%clear()
      
      ! TriDiagonal matrix, it will go to LU method
      
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',GMRES)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',TRIDIAG)
      CALL pList%add('LinearSolverType->A->MatrixType->n',3_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->isSym',.TRUE.)
      CALL pList%add('LinearSolverType->x->VectorType->n',3_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',3_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      !A=[ 4 -1  0]
      !  [-1  4 -1]
      !  [ 0 -1  4]
      SELECTTYPE(A => thisLS%A); TYPE IS(TriDiagMatrixType)
        CALL A%set(1,1, 4._SRK)
        CALL A%set(1,2,-1._SRK)
        CALL A%set(2,2, 4._SRK)
        CALL A%set(2,3,-1._SRK)
        CALL A%set(3,3, 4._SRK)
      ENDSELECT
      
      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1,1._SRK)
        CALL b%set(2,2._SRK)
        CALL b%set(3,3._SRK)
      ENDSELECT

      CALL thisLS%solve()
      
      IF(thisLS%info /= 0) THEN
        WRITE(*,*) 'CALL Iterative%solve() - GMRES FAILED!'
        STOP 666
      ENDIF
      CALL thisLS%A%clear()
      CALL thisLS%clear()
      DEALLOCATE(thisX)

    !DenseRect matrix
      
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',NATIVE)
      CALL pList%add('LinearSolverType->solverMethod',GMRES)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',DENSERECT)
      CALL pList%add('LinearSolverType->A->MatrixType->n',3_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->m',2_SNK)
      CALL pList%add('LinearSolverType->x->VectorType->n',2_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',3_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      ! A=[1 1]  b=[1]
      !   [1 2]    [2]
      !   [1 3]    [2]
      SELECTTYPE(A => thisLS%A); TYPE IS(DenseRectMatrixType)
        CALL A%set(1,1,1._SRK)
        CALL A%set(1,2,1._SRK)
        CALL A%set(2,1,1._SRK)
        CALL A%set(2,2,2._SRK)
        CALL A%set(3,1,1._SRK)
        CALL A%set(3,2,3._SRK)
      ENDSELECT
      
      SELECTTYPE(b => thisLS%b); TYPE IS(RealVectorType)
        CALL b%set(1,1._SRK)
        CALL b%set(2,2._SRK)
        CALL b%set(2,2._SRK)
      ENDSELECT
      
      ! initialize X0
      ALLOCATE(thisX(2))
      thisX=1.0_SRK

      !set iterations and convergence information and build/set M
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        CALL thisLS%setX0(thisX)
        CALL thisLS%setConv(2_SIK,1.0E-9_SRK,1000_SIK,30_SIK)
      ENDSELECT
      
      !Solve it
      CALL thisLS%solve()
      
      IF(thisLS%info /= 0) THEN
        WRITE(*,*) 'CALL Iterative%solve() -GMRES method FAILED!'
        STOP 666
      ENDIF
     
      DEALLOCATE(thisB)
      CALL thisLS%clear()
      WRITE(*,*)'  Passed: CALL Iterative%solver() GMRES'
      
#ifdef HAVE_PETSC
      !With GMRES
      !The sparse matrix type
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',PETSC)
      CALL pList%add('LinearSolverType->solverMethod',GMRES)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',SPARSE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',9_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->isSym',.TRUE.)
      CALL pList%add('LinearSolverType->x->VectorType->n',9_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',9_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      !A =  4    -1     0    -1     0     0     0     0     0
      !    -1     4    -1     0    -1     0     0     0     0
      !     0    -1     4     0     0    -1     0     0     0
      !    -1     0     0     4    -1     0    -1     0     0
      !     0    -1     0    -1     4    -1     0    -1     0
      !     0     0    -1     0    -1     4     0     0    -1
      !     0     0     0    -1     0     0     4    -1     0
      !     0     0     0     0    -1     0    -1     4    -1
      !     0     0     0     0     0    -1     0    -1     4
      SELECTTYPE(A => thisLS%A); TYPE IS(PETSCMatrixType)
          CALL A%set(1,1, 4.0_SRK)
          CALL A%set(1,2,-1.0_SRK)
          CALL A%set(1,4,-1.0_SRK)
          CALL A%set(2,1,-1.0_SRK)
          CALL A%set(2,2, 4.0_SRK)
          CALL A%set(2,3,-1.0_SRK)
          CALL A%set(2,5,-1.0_SRK)
          CALL A%set(3,2,-1.0_SRK)
          CALL A%set(3,3, 4.0_SRK)
          CALL A%set(3,6,-1.0_SRK)
          CALL A%set(4,1,-1.0_SRK)
          CALL A%set(4,4, 4.0_SRK)
          CALL A%set(4,5,-1.0_SRK)
          CALL A%set(4,7,-1.0_SRK)
          CALL A%set(5,2,-1.0_SRK)
          CALL A%set(5,4,-1.0_SRK)
          CALL A%set(5,5, 4.0_SRK)
          CALL A%set(5,6,-1.0_SRK)
          CALL A%set(5,8,-1.0_SRK)
          CALL A%set(6,3,-1.0_SRK)
          CALL A%set(6,5,-1.0_SRK)
          CALL A%set(6,6, 4.0_SRK)
          CALL A%set(6,9,-1.0_SRK)
          CALL A%set(7,4,-1.0_SRK)
          CALL A%set(7,7, 4.0_SRK)
          CALL A%set(7,8,-1.0_SRK)
          CALL A%set(8,5,-1.0_SRK)
          CALL A%set(8,7,-1.0_SRK)
          CALL A%set(8,8, 4.0_SRK)
          CALL A%set(8,9,-1.0_SRK)
          CALL A%set(9,6,-1.0_SRK)
          CALL A%set(9,8,-1.0_SRK)
          CALL A%set(9,9, 4.0_SRK)
      ENDSELECT
      
      ! build X0 and set it to 1.0s
      ALLOCATE(thisX(9))
      thisX=1.0_SRK
      
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
          CALL thisLS%setX0(thisX)
      ENDSELECT

      SELECTTYPE(b => thisLS%b); TYPE IS(PETScVectorType)
        CALL b%set(1.0_SRK)
      ENDSELECT

      !set iterations and convergence information and build/set M
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        CALL thisLS%setConv(2_SIK,1.0E-9_SRK,1000_SIK,30_SIK)
      ENDSELECT
      
      !solve it
      CALL thisLS%solve()
      
      !Store expected solution (from MATLAB) in B
      ALLOCATE(thisB(9))
      thisB(1)=0.6875_SRK
      thisB(2)=0.875_SRK
      thisB(3)=0.6875_SRK
      thisB(4)=0.875_SRK
      thisB(5)=1.125_SRK
      thisB(6)=0.875_SRK
      thisB(7)=0.6875_SRK
      thisB(8)=0.875_SRK
      thisB(9)=0.6875_SRK
      !multiply by 10,000 so we can match first five places.
      thisB=10000.0_SRK*thisB
      match=.TRUE.
      DO i=1,SIZE(thisB)
        SELECTTYPE(X => thisLS%X); TYPE IS(PETScVectorType)
          IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
          ALLOCATE(dummyvec(X%n))
          CALL X%get(dummyvec)
          IF(NINT(thisB(i)) /= NINT(10000.0_SRK*dummyvec(i))) THEN
            match=.FALSE.
            EXIT
          ENDIF
        ENDSELECT
      ENDDO
      IF(.NOT. match) THEN
        WRITE(*,*) 'CALL PETScIterative%solve() -BICGSTAB FAILED!'
        STOP 666
      ENDIF
      
      DEALLOCATE(thisB)
      CALL thisLS%A%clear()
      CALL thisLS%clear()
      
    !test with A being densesquare
    
      ! initialize linear system
      CALL pList%clear()
      CALL pList%add('LinearSolverType->TPLType',PETSC)
      CALL pList%add('LinearSolverType->solverMethod',GMRES)
      CALL pList%add('LinearSolverType->MPI_Comm_ID',PE_COMM_SELF)
      CALL pList%add('LinearSolverType->numberOMP',1_SNK)
      CALL pList%add('LinearSolverType->timerName','testTimer')
      CALL pList%add('LinearSolverType->A->MatrixType->matType',DENSESQUARE)
      CALL pList%add('LinearSolverType->A->MatrixType->n',9_SNK)
      CALL pList%add('LinearSolverType->A->MatrixType->isSym',.TRUE.)
      CALL pList%add('LinearSolverType->x->VectorType->n',9_SNK)
      CALL pList%add('LinearSolverType->b->VectorType->n',9_SNK)
      CALL pList%validate(pList,optListLS)
      CALL thisLS%init(pList)
      
      DO i=1,9
        SELECTTYPE(A => thisLS%A); TYPE IS(PETScMatrixType)
          CALL A%set(i,i,4.0_SRK)
          IF((i < 9).AND.((i /= 3).AND.(i /= 6))) THEN
            CALL A%set(i,i+1,-1.0_SRK)
          ENDIF
          IF(i < 7) THEN
            CALL A%set(i,i+3,-1.0_SRK)
          ENDIF
        ENDSELECT
      ENDDO
      
      !build X0 and set it to 1.0s
      ALLOCATE(thisX(9))
      thisX=1.0_SRK
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        CALL thisLS%setX0(thisX)
      ENDSELECT
      
      SELECTTYPE(b => thisLS%b); TYPE IS(PETScVectorType)
        CALL b%set(1.0_SRK)
      ENDSELECT
      
      !set iterations and convergence information and build/set M
      SELECTTYPE(thisLS); TYPE IS(LinearSolverType_Iterative)
        CALL thisLS%setConv(2_SIK,1.0E-9_SRK,1000_SIK,30_SIK)
      ENDSELECT
      
      !solve
      CALL thisLS%solve()
      ALLOCATE(thisB(9))
      thisB(1)=0.6875_SRK
      thisB(2)=0.875_SRK
      thisB(3)=0.6875_SRK
      thisB(4)=0.875_SRK
      thisB(5)=1.125_SRK
      thisB(6)=0.875_SRK
      thisB(7)=0.6875_SRK
      thisB(8)=0.875_SRK
      thisB(9)=0.6875_SRK
      thisB=thisB*10000._SRK
      match=.TRUE.
      DO i=1,SIZE(thisB)
        SELECTTYPE(X => thisLS%X); TYPE IS(PETScVectorType)
          IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
          ALLOCATE(dummyvec(X%n))
          CALL X%get(dummyvec)
          IF(NINT(thisB(i)) /= NINT(10000.0_SRK*dummyvec(i))) THEN
            match=.FALSE.
            EXIT
          ENDIF
        ENDSELECT
      ENDDO
      IF(.NOT. match) THEN
        WRITE(*,*) 'CALL PETScIterative%solve() - GMRES FAILED!'
        STOP 666
      ENDIF
      !test to see how it performs with an already decomposed M
      !reset X to 1.0s
      match=.TRUE.
      DO i=1,SIZE(thisB)
        SELECTTYPE(X => thisLS%X); TYPE IS(PETScVectorType)
          CALL X%set(1.0_SRK)
          CALL thisLS%solve()
          IF(ALLOCATED(dummyvec)) DEALLOCATE(dummyvec)
          ALLOCATE(dummyvec(X%n))
          CALL X%get(dummyvec)
          IF(NINT(thisB(i)) /= NINT(10000.0_SRK*dummyvec(i))) THEN
            match=.FALSE.
            EXIT
          ENDIF
        ENDSELECT
      ENDDO
      IF(.NOT. match) THEN
        WRITE(*,*) 'CALL PETScIterative%solve() - GMRES FAILED!'
        STOP 666
      ENDIF
      CALL thisLS%A%clear()
            
      WRITE(*,*)'  Passed: CALL PETScIterative%solver() GMRES'
#endif
      CALL thisLS%clear()
      DEALLOCATE(thisX)
      
    ENDSUBROUTINE testIterativeSolve_GMRES
!
!-------------------------------------------------------------------------------
    SUBROUTINE testNorms()
      INTEGER(SIK) :: normType
      REAL(SRK),DIMENSION(10) :: x
      REAL(SRK) :: norm
      !set up x
      x=(/0._SRK,-1._SRK,2._SRK,-3._SRK,4._SRK, & 
          -5._SRK,6._SRK,-7._SRK,8._SRK,-9._SRK/)
      !test normType 1
      normType=1 !taxicab norm, just the absolute sum of these.
      !expected answer = 45
      CALL LNorm(x,normType,norm)
      IF(.NOT. (norm .APPROXEQ. 45._SRK)) THEN
        WRITE(*,*) 'CALL LNorm() - 1-NORM FAILED!'
        STOP 666
      ENDIF
      !test normType 2
      normType=2 !Euclidean norm
      !expected answer = 16.88194301613413218312
      CALL LNorm(x,normType,norm)
      IF(.NOT. (norm .APPROXEQ. 16.88194301613413218312_SRK)) THEN
        WRITE(*,*) 'CALL LNorm() - 2-NORM FAILED!'
        STOP 666
      ENDIF
      !test normType -1
      normType=-1 !Infinite norm
      CALL LNorm(x,normType,norm)
      !expected answer = 9.0
      IF(.NOT. (norm .APPROXEQ. 9.0_SRK)) THEN
        WRITE(*,*) 'CALL LNorm() - INFINITE-NORM FAILED!'
        WRITE(*,*) norm
        STOP 666
      ENDIF
      !test normType 3 (just some p-norm)
      normType=3 !L-norm, w/ L=3
      CALL LNorm(x,normType,norm)
      !expected answer = 12.65148997952623864269
      IF(.NOT. (norm .APPROXEQ. 12.65148997952623864269_SRK)) THEN
        WRITE(*,*) 'CALL LNorm() - L-NORM FAILED!'
        WRITE(*,*) norm
        STOP 666
      ENDIF
      !test an invalid norm (<=-2)
      normType=-2
      CALL LNorm(x,normType,norm)
      !expected answer = 0.0
      IF(.NOT. (norm == 0.0_SRK)) THEN
        WRITE(*,*) 'CALL LNorm() - L-NORM FAILED!'
        WRITE(*,*) norm
        STOP 666
      ENDIF
      WRITE(*,*) '  Passed: CALL LNORM(...)'
    ENDSUBROUTINE testNorms

ENDPROGRAM testLinearSolver
