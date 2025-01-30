!
! Copyright (C) 2025 Levin Rug, Willi Schimmel (E-Mail: l.rug@lmu.de)
! See ./SRC/Cminor.f90 for the copyright notice
! See ./LICENSE for license information
! SPDX-License-Identifier: GPL-3.0
!
!================================================================================!
!                                                                                !
!                         This module includes a collection                      !
!                    of sparse matrix calculations for chemical                  !
!                       reaction systems, main format CSR                        !
!                                                                                !
!================================================================================!
!                                                                                
MODULE Sparse_Mod
  ! Contains:
  !   - SYMBOLIC MATRIX * MATRIX
  !   - SYMBOLIC MATRIX + MATRIX
  !   - TRANSPOSE MATRIX
  !   - SPARSE MATRIX * MATRIC
  !   - SPARSE MATRIX +- MATRIC
  !   - SORT ALGORITHM FOR SYMBOLIC MATRIX*MATRIX CALC
  !   - SPARSE JACOBIAN MATRIX CALC
  !   - SPARSE MITER MATRIX CALC
  !   - CONVERT COMPRESSED ROW FORMAT TO ROW-INDEX, COL-INDEX
  !   - PRINT SPARSE MATRIX (compressed Row format)
  !   - PRINT SPARSE MATRIX (ROW-INDEX COLUMN-INDEX)
  !   - SPARSE MATRIX*VECTOR+VECTOR (SPARSE-MATRIX, DENSE VECTORS)
  ! 
  USE UniRnk_Mod,  ONLY: unirnk
  USE Control_Mod, ONLY: nDropletClasses, ZERO, combustion, ONE,     &
                       & nD_spc, nD_Ptr_spc, nD_reac, nD_Ptr_reacs,  &
                       & milli, adiabatic_parcel
  !
  USE Reac_Mod,    ONLY: nDIM, nreac, nspc
  !
  USE Kind_Mod,    ONLY: dp
  !
  IMPLICIT NONE 
  ! 
  INTEGER, PARAMETER, PRIVATE :: inilen=400
  !
  !
  TYPE CSR_Matrix_T            !Compressed Rowindex, standart columnindex
    INTEGER               :: m=0, n=0, nnz=0
    INTEGER, ALLOCATABLE  :: RowPtr(:)
    INTEGER, ALLOCATABLE  :: ColInd(:)
    INTEGER, ALLOCATABLE  :: DiagPtr(:)         ! position of diag entries
    INTEGER, ALLOCATABLE  :: DiagPtr_P(:)       ! permuted pos of diag entries
    INTEGER, ALLOCATABLE  :: DiagPtr_R(:)       ! permuted pos of rate entries
    INTEGER, ALLOCATABLE  :: DiagPtr_C(:)       ! permuted pos of conenctations
    INTEGER, ALLOCATABLE  :: RowVectorPtr(:)    ! for combustion matrix ( -U^T*D_c )
    INTEGER, ALLOCATABLE  :: ColVectorPtr(:)    ! for combustion matrix ( ~K )
    INTEGER               :: XPtr=-42           ! for combustion matrix ( X )
    INTEGER, ALLOCATABLE  :: Permu(:)           ! permutation vector (markowitz)
    INTEGER, ALLOCATABLE  :: InvPer(:)          ! invers permutation (inv markowitz)
    INTEGER, ALLOCATABLE  :: LUperm(:)          ! 
    REAL(dp), ALLOCATABLE :: Val(:)
    INTEGER, ALLOCATABLE  :: ValPtr(:)          ! for vector-valued entries, store start and end points of the entries
    LOGICAL               :: vector_entried=.FALSE. ! logical variable to store if a matrix is supposed to be vector-entried
  END TYPE CSR_Matrix_T
  !
  TYPE SpRowIndColInd_T     !standart Rowindex, standart columnindex
    INTEGER :: m=0,n=0,nnz=0
    INTEGER, ALLOCATABLE :: RowInd(:)
    INTEGER, ALLOCATABLE :: ColInd(:)
    REAL(dp), ALLOCATABLE :: Val(:)
  END TYPE SpRowIndColInd_T
  !
  TYPE SpRowColD_T
    INTEGER :: m=0,n=0
    INTEGER, POINTER :: RowPtr(:,:)=>NULL()
    INTEGER, POINTER :: ColInd(:)=>NULL()
    INTEGER, POINTER :: Permu(:)=>NULL()
    INTEGER, POINTER :: InvPer(:)=>NULL()
    INTEGER, ALLOCATABLE :: Restr(:) ! make pointer ???
    INTEGER :: ep=1
    INTEGER :: last=0
    INTEGER :: len=0
    INTEGER :: nnz=0
  END TYPE SpRowColD_T
  !  
  ! global matrices containing chemical reaktion data (stoech.coefs)
  TYPE(CSR_Matrix_T) ::  A                & ! coef matrix of educts
  &                    , B                & ! coef matrix of products
  &                    , BA               & ! B-A
  &                    , BAT                ! Transpose(B-A)
  
  TYPE(CSR_Matrix_T) :: TB_sparse ! sparse matrix containing thirdbody 
  !
  TYPE(CSR_Matrix_T) :: Jac_CC

  TYPE(CSR_Matrix_T) :: Miter
  TYPE(CSR_Matrix_T) :: LU_Miter                    
  !

  INTERFACE OPERATOR(*)
    MODULE PROCEDURE DAX_sparse
    MODULE PROCEDURE DAX_T_sparse
  END INTERFACE

  !
  CONTAINS
  
  FUNCTION New_CSR(m,n,nnz,ri,ci,val) RESULT(newA)
    INTEGER :: m, n
    INTEGER, OPTIONAL :: nnz
    INTEGER, OPTIONAL :: ri(:), ci(:)
    REAL(dp), OPTIONAL :: val(:)
    TYPE(CSR_Matrix_T) :: newA
    !
    INTEGER :: i, sameRow, cCnt

    newA%m = m
    newA%n = n
    !
    ALLOCATE(newA%RowPtr(m+1))
    newA%RowPtr = 0
    newA%RowPtr(1) = 1
    !
    IF (PRESENT(nnz)) THEN
      ALLOCATE(newA%ColInd(nnz))
      newA%ColInd = -1
      ALLOCATE(newA%Val(nnz))
      newA%Val = ZERO
      newA%nnz = nnz
    END IF
    
    ! if row and column indices are given
    IF (PRESENT(ri).AND.PRESENT(ci)) THEN
      IF (SIZE(ri) /= SIZE(ci)) STOP ' SIZE(ri) /= SIZE(ci) '
      cCnt = 0
      DO i = 1,m
        sameRow = COUNT(ri==i)
        newA%RowPtr(i+1) = newA%RowPtr(i) + sameRow
      END DO
      newA%ColInd = ci
      IF (PRESENT(val)) newA%val = val
    END IF
  END FUNCTION New_CSR
  !
  !
  SUBROUTINE Free_Matrix_CSR(A)
    TYPE(CSR_Matrix_T) :: A
    !
    IF (ALLOCATED(A%RowPtr))    DEALLOCATE(A%RowPtr)
    IF (ALLOCATED(A%ColInd))    DEALLOCATE(A%ColInd)
    IF (ALLOCATED(A%DiagPtr))   DEALLOCATE(A%DiagPtr)
    IF (ALLOCATED(A%DiagPtr_R)) DEALLOCATE(A%DiagPtr_R)
    IF (ALLOCATED(A%DiagPtr_C)) DEALLOCATE(A%DiagPtr_C)
    IF (ALLOCATED(A%Permu))     DEALLOCATE(A%Permu)
    IF (ALLOCATED(A%InvPer))    DEALLOCATE(A%InvPer)  
    IF (ALLOCATED(A%Val))       DEALLOCATE(A%Val)
  END SUBROUTINE Free_Matrix_CSR
  !
  !
  SUBROUTINE Free_SpRowColD(A)
    TYPE(SpRowColD_T) :: A
    !

    A%m=0
    A%n=0
    A%ep=1
    A%last=0
    A%len=0
    A%nnz=0
    IF (ASSOCIATED(A%RowPtr)) NULLIFY(A%RowPtr)
    IF (ASSOCIATED(A%ColInd)) NULLIFY(A%ColInd)
    IF (ASSOCIATED(A%Permu )) NULLIFY(A%Permu)
    IF (ASSOCIATED(A%InvPer)) NULLIFY(A%InvPer)
  END SUBROUTINE Free_SpRowColD
  !
  !
  SUBROUTINE Free_SpRowIndColInd(A)
    TYPE(SpRowIndColInd_T) :: A
    !
    IF (ALLOCATED(A%RowInd)) DEALLOCATE(A%RowInd)
    IF (ALLOCATED(A%ColInd)) DEALLOCATE(A%ColInd)
    IF (ALLOCATED(A%Val))    DEALLOCATE(A%Val)
  END SUBROUTINE Free_SpRowIndColInd   
  !
  !
  FUNCTION SparseID(dim) RESULT(Mat)
    TYPE(CSR_Matrix_T) :: Mat
    INTEGER :: dim
    INTEGER :: i
    !
    Mat = New_CSR(dim,dim,dim)
    DO i=1,dim
      Mat%RowPtr(i+1)=Mat%RowPtr(i)+1
      Mat%ColInd(i)=i
    END DO
    Mat%Val = ONE
  END FUNCTION SparseID
  !

  FUNCTION CSR_to_FULL(CSR) RESULT(Full)
    TYPE(CSR_Matrix_T) :: CSR
    REAL(dp) :: Full(CSR%m,CSR%n)

    INTEGER :: i,j,jj 
    
    Full = ZERO
    DO i=1,CSR%m
      DO jj=CSR%RowPtr(i),CSR%RowPtr(i+1)-1
        j = CSR%ColInd(jj)
        Full(i,j) = CSR%Val(jj)
      END DO
    END DO
  END FUNCTION CSR_to_FULL


  SUBROUTINE CompressIntegerArray(Array)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: Array(:)
    INTEGER, ALLOCATABLE :: tmpArray(:)
    
    INTEGER :: i, N, cnt
    
    N = COUNT(Array/=0)
    ALLOCATE(tmpArray(N))

    cnt = 0
    DO i=1,SIZE(Array)
      IF (Array(i)/=0) THEN
        cnt = cnt + 1
        tmpArray(cnt) = Array(i)
      END IF
    END DO
    Array = [tmpArray]
  END SUBROUTINE CompressIntegerArray



  SUBROUTINE CompressDoubleArray(Array)
    REAL(dp), ALLOCATABLE, INTENT(INOUT) :: Array(:)
    REAL(dp), ALLOCATABLE :: tmpArray(:)
    
    INTEGER :: i, N, cnt
    REAL(dp), PARAMETER :: big = -99999999999999.d0
    
    N = COUNT( Array /= big )
    ALLOCATE(tmpArray(N))

    cnt = 0
    DO i=1,SIZE(Array)
      IF ( Array(i) /= big ) THEN
        cnt = cnt + 1
        tmpArray(cnt) = Array(i)
      END IF
    END DO
    Array = [tmpArray]
  END SUBROUTINE CompressDoubleArray


  FUNCTION FULL_to_CSR(Full) RESULT(CSR)
    REAL(dp) :: Full(:,:)
    TYPE(CSR_Matrix_T) :: CSR

    INTEGER :: i,j,jj 
    INTEGER :: m,n,nnz 

    m = SIZE(Full,1)
    n = SIZE(Full,2)
    nnz = COUNT(Full /= ZERO)

    CSR = New_CSR(m,n,nnz)
    
    jj = 0
    DO i = 1 , m
      CSR%RowPtr(i+1) = CSR%RowPtr(i) + COUNT(Full(i,:)/=ZERO)
      DO j = 1 , n
        IF ( Full(i,j) /= ZERO) THEN
          jj = jj + 1
          CSR%ColInd(jj) = j
          CSR%Val(jj)    = Full(i,j)
        END IF
      END DO
    END DO
    CSR%nnz = jj
  END FUNCTION FULL_to_CSR


  FUNCTION Copy_SpRowIndColInd(orig) RESULT(copy)
    TYPE(SpRowIndColInd_T) :: orig
    TYPE(SpRowIndColInd_T) :: copy

    copy%m = orig%m
    copy%n = orig%n
    copy%nnz = orig%nnz
    copy%RowInd = orig%RowInd
    copy%ColInd = orig%ColInd
    copy%Val    = orig%Val
  END FUNCTION Copy_SpRowIndColInd


  FUNCTION Copy_SpRowColD(orig) RESULT(copy)
    TYPE(SpRowColD_T) :: orig
    TYPE(SpRowColD_T) :: copy

    copy%m = orig%m
    copy%n = orig%n
    copy%nnz = orig%nnz
    copy%RowPtr = orig%RowPtr
    copy%ColInd = orig%ColInd
    copy%Permu  = orig%Permu
    copy%InvPer = orig%InvPer

  END FUNCTION Copy_SpRowColD


  FUNCTION Copy_CSR(orig) RESULT(copy)
    TYPE(CSR_Matrix_T) :: orig
    TYPE(CSR_Matrix_T) :: copy

    copy = New_CSR(orig%m,orig%n,orig%nnz)
    copy%RowPtr = orig%RowPtr
    copy%ColInd = orig%ColInd
    copy%Val    = orig%Val
    Copy%nnz    = orig%nnz
    IF (orig%XPtr/=-42) copy%XPtr = orig%XPtr
    IF (ALLOCATED(orig%DiagPtr))   copy%DiagPtr   = orig%DiagPtr
    IF (ALLOCATED(orig%DiagPtr_P)) copy%DiagPtr_P = orig%DiagPtr_P
    IF (ALLOCATED(orig%DiagPtr_R)) copy%DiagPtr_R = orig%DiagPtr_R
    IF (ALLOCATED(orig%DiagPtr_C)) copy%DiagPtr_C = orig%DiagPtr_C
    IF (ALLOCATED(orig%RowVectorPtr)) copy%RowVectorPtr = orig%RowVectorPtr
    IF (ALLOCATED(orig%ColVectorPtr)) copy%ColVectorPtr = orig%ColVectorPtr
    IF (ALLOCATED(orig%Permu))  copy%Permu  = orig%Permu
    IF (ALLOCATED(orig%InvPer)) copy%InvPer = orig%InvPer
    IF (ALLOCATED(orig%LUperm)) copy%LUperm = orig%LUperm
  END FUNCTION Copy_CSR

  !
  FUNCTION RowColD_to_CSR(SpRow) RESULT(CSR)
    TYPE(CSR_Matrix_T) :: CSR
    TYPE(SpRowColD_T)  :: SpRow
    !
    INTEGER :: i,j,jj,nzrA
    ! 

    CSR = New_CSR(SpRow%n,SpRow%m,SpRow%nnz)

    ALLOCATE(CSR%DiagPtr(nDIM))
    ALLOCATE(CSR%DiagPtr_P(nDIM))
    CSR%DiagPtr   = -11
    CSR%DiagPtr_P = -11

    NzrA=0
    DO i=1,CSR%m
      CSR%RowPtr(i+1)=CSR%RowPtr(i)
      DO jj=SpRow%RowPtr(1,i),SpRow%RowPtr(2,i)
        j=SpRow%ColInd(jj) 
        NzrA=NzrA+1
        CSR%ColInd(NzrA)=j
        IF (i==j) THEN
          CSR%DiagPtr(i)=CSR%RowPtr(i+1)
        END IF
        CSR%RowPtr(i+1)=CSR%RowPtr(i+1)+1
      END DO
    END DO
    
    IF (ASSOCIATED(SpRow%Permu)) THEN
      IF (.NOT.ALLOCATED(CSR%Permu)) THEN
        ALLOCATE(CSR%Permu(CSR%n))
        CSR%Permu(:)=-16
      END IF
      CSR%Permu(:)=SpRow%Permu(:)
    END IF
    IF (ASSOCIATED(SpRow%InvPer)) THEN
      IF (.NOT.ALLOCATED(CSR%InvPer)) THEN
        ALLOCATE(CSR%InvPer(CSR%n))
        CSR%InvPer(:)=-16
      END IF
      CSR%InvPer(:)=SpRow%InvPer(:)
    END IF

    ! permutate pointer to diagonal, rate and conc entries
    CSR%DiagPtr_P = CSR%DiagPtr( CSR%Permu )
    CSR%nnz = SpRow%nnz
  END FUNCTION RowColD_to_CSR

  FUNCTION SpRowIndColInd_to_CSR(SpRiCi) RESULT(CSR)
    ! IN:
    TYPE(SpRowIndColInd_T) :: SpRiCi
    ! OUT:
    TYPE(CSR_Matrix_T)     :: CSR
    ! TEMP:
    INTEGER :: i, ii, jj, n, m, nnz

    m = SpRiCi%m
    n = SpRiCi%n
    nnz = SpRiCi%nnz

    CSR = New_CSR(m,n,nnz)
    
    ii = 0
    DO i = 1 , m
      CSR%RowPtr(i+1) = CSR%RowPtr(i) + COUNT(SpRiCi%RowInd==i)
      DO jj = CSR%RowPtr(i) , CSR%RowPtr(i+1)-1
        ii = ii + 1
        CSR%ColInd(jj) = SpRiCi%ColInd(ii)
        CSR%Val(jj)    = SpRiCi%Val(ii)
      END DO
    END DO
    CSR%nnz = ii
  END FUNCTION SpRowIndColInd_to_CSR
  !
  !
  FUNCTION CSR_to_SpRowColD(CSR) RESULT(SpRowCol)
    TYPE(SpRowColD_T) :: SpRowCol
    TYPE(CSR_Matrix_T) :: CSR
    !
    INTEGER :: nzrCSR
    INTEGER :: AddLen
    INTEGER :: Start,End
    INTEGER :: i,j,jj
    !
    AddLen=10
    SpRowCol%m=CSR%m
    SpRowCol%n=CSR%n
    ALLOCATE(SpRowCol%RowPtr(2,SpRowCol%n))
    nzrCSR=SIZE(CSR%ColInd)
    SpRowCol%len=10*nzrCSR+AddLen*SpRowCol%n
    ALLOCATE(SpRowCol%ColInd(SpRowCol%len))
    SpRowCol%ColInd=0
    ALLOCATE(SpRowCol%Permu(SpRowCol%n))
    SpRowCol%Permu=0
    ALLOCATE(SpRowCol%InvPer(SpRowCol%n))
    SpRowCol%InvPer=0
    !
    Start=1
    DO i=1,CSR%n
      SpRowCol%RowPtr(1,i)=Start
      End=Start+(CSR%RowPtr(i+1)-CSR%RowPtr(i)-1)
      SpRowCol%RowPtr(2,i)=End
      DO jj=CSR%RowPtr(i),CSR%RowPtr(i+1)-1
        j=CSR%ColInd(jj)
        SpRowCol%ColInd(Start)=j
        Start=Start+1
      END DO  
      Start=Start+AddLen-1
    END DO  
    SpRowCol%ep=SpRowCol%RowPtr(2,SpRowCol%n)+1
    SpRowCol%last=SpRowCol%n  
  END FUNCTION CSR_to_SpRowColD
  !
  !
  SUBROUTINE Sort_SpRowColD(A)
    TYPE(SpRowColD_T) :: A
    INTEGER :: i,jj
    DO i=1,A%m
      DO jj=A%RowPtr(1,i),A%RowPtr(2,i)
        A%ColInd(jj)=A%Permu(A%ColInd(jj))
      END DO
      CALL SortVec(A%ColInd(A%RowPtr(1,i):A%RowPtr(2,i)))
    END DO
  END SUBROUTINE Sort_SpRowColD
  !
  !
  SUBROUTINE SymbLU_SpRowColD(A,Permu)
    TYPE(SpRowColD_T) :: A
    INTEGER :: Permu(:)
    !
    INTEGER :: RowPiv(A%n)
    INTEGER :: i,j,l,jj,ip,iPiv
    LOGICAL :: ins
    !
    DO i=1,A%n
      A%InvPer(i)=i
      A%Permu(i)=i
    END DO
    !
    DO i=1,A%n
      ip=A%Permu(Permu(i))
      CALL Swap(A%InvPer(i),A%InvPer(ip))
      CALL Swap(A%RowPtr(1,i),A%RowPtr(1,ip))
      CALL Swap(A%RowPtr(2,i),A%RowPtr(2,ip))
      A%Permu(A%InvPer(i))=i
      A%Permu(A%InvPer(ip))=ip
      IF (A%last==i) THEN
        A%last=ip
      ELSE IF (A%last==ip) THEN
        A%last=i
      END IF
      !
      !   Update
      iPiv=0
      DO jj=A%RowPtr(1,i),A%RowPtr(2,i)
        IF (A%Permu(A%ColInd(jj))>i) THEN
          iPiv=iPiv+1 
          RowPiv(iPiv)=A%ColInd(jj)
        END IF
      END DO
      IF (iPiv>0) THEN
        DO j=i+1,A%n
          DO jj=A%RowPtr(1,j),A%RowPtr(2,j)
            IF (A%Permu(A%ColInd(jj))==i) THEN
              DO l=1,iPiv
                CALL Insert_SpRowColD(A,j,RowPiv(l),ins)
              END DO
              EXIT
            END IF
          END DO
        END DO
      END IF
    END DO
    DO i=1,A%n
      DO jj=A%RowPtr(1,i),A%RowPtr(2,i)
        A%ColInd(jj)=A%Permu(A%ColInd(jj))
      END DO
      CALL SortVec(A%ColInd(A%RowPtr(1,i):A%RowPtr(2,i)))
      A%nnz=A%nnz+SIZE(A%ColInd(A%RowPtr(1,i):A%RowPtr(2,i)))
    END DO
  END SUBROUTINE SymbLU_SpRowColD
  !
  !
  ! this function has to become the new normal Symbolic lu decomposition !
  ! but restr has to be constructed while initializing, so later: delete old function and
  ! rename this one to itself minus _with_Restr
  SUBROUTINE SymbLU_SpRowColD_M(A) 
    TYPE(SpRowColD_T) :: A
    !
    INTEGER :: r(A%n),c(A%n),RowPiv(A%n)
    INTEGER :: i,j,l,jj,ip,ip1(1),iPiv
    REAL(dp) :: md 
    INTEGER ::  rc
    LOGICAL :: ins

    c = 0
    DO i = 1 , A%n
      A%InvPer(i) = i
      A%Permu(i)  = i
      ! Compute initial Markowitz count
      r(i) = A%RowPtr(2,i) - A%RowPtr(1,i) + 1
      DO jj = A%RowPtr(1,i) , A%RowPtr(2,i)
        c(A%ColInd(jj)) = c(A%ColInd(jj)) + 1
      END DO
    END DO
    
    MAIN_LOOP: DO i = 1 , A%n
      ip = 0
      md = 1.d99
      
      DO j = i , A%n
        rc = (r(j)-1) * (c(j)-1)
        IF ( rc <= md .AND. A%Restr(j)>=0) THEN
        !IF ( rc < md .AND. A%Restr(j)>=0) THEN
          md = rc
          ip = j
        END IF
      END DO

      ! only restr species left -> permute restr species among themselves
      IF (ip==0) THEN
        ip1(:)=MINLOC((r(i:A%n)-1)*(c(i:A%n)-1))+(i-1)
        ip=ip1(1)
      END IF
      
      CALL Swap( r(i) , r(ip) )
      CALL Swap( c(i) , c(ip) )
      CALL Swap( A%Restr(i), A%Restr(ip) )
      CALL Swap( A%InvPer(i)   , A%InvPer(ip) )
      CALL Swap( A%RowPtr(1,i) , A%RowPtr(1,ip) )
      CALL Swap( A%RowPtr(2,i) , A%RowPtr(2,ip) )
      
      A%Permu(A%InvPer(i))  = i
      A%Permu(A%InvPer(ip)) = ip

      IF ( A%last == i ) THEN
        A%last = ip
      ELSE IF ( A%last == ip ) THEN
        A%last = i
      END IF
      !
      ! Update
      iPiv = 0
      DO jj = A%RowPtr(1,i) , A%RowPtr(2,i)
        IF ( A%Permu(A%ColInd(jj)) > i ) THEN
          iPiv = iPiv + 1 
          RowPiv(iPiv) = A%ColInd(jj)
          c(A%Permu(A%ColInd(jj))) = c(A%Permu(A%ColInd(jj))) - 1
        END IF
      END DO
      IF ( iPiv > 0 ) THEN
        DO j = i+1 , A%n
          DO jj = A%RowPtr(1,j) , A%RowPtr(2,j)
            IF ( A%Permu(A%ColInd(jj)) == i ) THEN
              r(j) = r(j) - 1 
              c(i) = c(i) - 1
              DO l = 1 , iPiv
                CALL Insert_SpRowColD( A, j, RowPiv(l), ins)
                IF (ins) THEN
                  c(A%Permu(RowPiv(l))) = c(A%Permu(RowPiv(l))) + 1
                  r(j) = r(j) + 1
                END IF
              END DO
              EXIT
            END IF
          END DO
          IF ( c(i) == 1 ) EXIT
        END DO
      END IF
    END DO MAIN_LOOP

    DO i = 1 , A%n
      DO jj = A%RowPtr(1,i) , A%RowPtr(2,i)
        A%ColInd(jj) = A%Permu(A%ColInd(jj))
      END DO
      CALL SortVec( A%ColInd(A%RowPtr(1,i):A%RowPtr(2,i)) )
      A%nnz = A%nnz + SIZE(A%ColInd(A%RowPtr(1,i):A%RowPtr(2,i)))
    END DO
  END SUBROUTINE SymbLU_SpRowColD_M
  !
  !
  SUBROUTINE Swap(i,j)
    INTEGER :: i, j, iTemp
    iTemp=i;    i=j;    j=iTemp
  END SUBROUTINE Swap
  !
  !
  SUBROUTINE SortVec(vec)
    INTEGER :: vec(:)
    !
    INTEGER :: i,itemp,j,n
    n=SIZE(Vec)
    DO i=1,n
      DO j=1,n-i
        IF (vec(j)>vec(j+1)) THEN
          itemp=vec(j)
          vec(j)=vec(j+1)
          vec(j+1)=itemp
        END IF
      END DO
    END DO
  END SUBROUTINE SortVec
  !
  !
  SUBROUTINE SortVecAsc(vec,q)
    INTEGER, INTENT(INOUT) :: vec(:)
    INTEGER, ALLOCATABLE, OPTIONAL :: q(:)
    !
    INTEGER :: i,n,iMin(1)
    INTEGER, ALLOCATABLE :: tmpVec(:)
   
    n = SIZE(Vec)
    ALLOCATE(tmpVec(n));  tmpVec = 0
   
    IF (PRESENT(q).AND..NOT.ALLOCATED(q)) ALLOCATE(q(n))
    
    DO i = 1 , n
      iMin = MINLOC(vec)
      tmpVec(i) = vec(iMin(1))
      vec(iMin(1)) = 99999999
      IF (PRESENT(q)) q(i) = iMin(1)
    END DO
    vec = tmpVec
  END SUBROUTINE SortVecAsc
  !
  SUBROUTINE Insert_SpRowColD(A,iA,jA,ins)
    TYPE(SpRowColD_T) :: A
    INTEGER :: iA,jA
    LOGICAL :: ins
    !
    INTEGER :: itemp,j,l
    INTEGER, ALLOCATABLE :: iWork(:)
    !
    ! Test ob Element (ia,ja) bereits enthalten
    IF (iA==0.OR.jA==0) THEN
      WRITE(*,*) 'iA',iA
      WRITE(*,*) 'jA',jA
      STOP 'STOP Insert_SpRowColD'
    END IF  
    ins=.TRUE.
    DO j=A%RowPtr(1,iA),A%RowPtr(2,iA)
      IF (jA==A%ColInd(j)) THEN
        ins=.FALSE.
      END IF
    END DO
    !
    ! Test auf freien Speicherplatz in der ia-ten
    ! Zeile von a
    ! 
    IF (ins) THEN
      IF (A%ColInd(A%RowPtr(2,iA)+1)/=0) THEN
        ! ja-te Zeile von a wird nach hinten
        itemp=A%ep
        DO l=A%RowPtr(1,iA),A%RowPtr(2,iA)
          A%ColInd(A%ep)=A%ColInd(l)
          A%ColInd(l)=0
          A%ep=A%ep+1
        END DO
        A%RowPtr(2,iA)=A%ep-1
        A%RowPtr(1,iA)=itemp
        ! A%ep=A%ep+1
        A%last=iA
      ENDIF
      A%RowPtr(2,iA)=A%RowPtr(2,iA)+1
      A%ColInd(A%RowPtr(2,iA))=jA
      IF (iA==A%last) A%ep=A%ep+1
    END IF
    !
    IF (A%ep>=A%len-A%m) THEN
      CALL gcmat_SpRowColD(A)
    END IF
    !
    IF (A%ep>=A%len-A%m) THEN
      !   Speicherplatz von A nicht ausreichend
      ALLOCATE(iWork(A%ep))
      iWork(1:A%ep)=A%ColInd(1:A%ep)
      DEALLOCATE(A%ColInd)
      A%len=2*A%len
      ALLOCATE(A%ColInd(A%len))
      A%ColInd(1:A%ep)=iWork(1:A%ep)
      DEALLOCATE(iWork)  
    END IF
  END SUBROUTINE Insert_SpRowColD
  !
  ! 
  SUBROUTINE gcmat_SpRowColD(A)
   !   Externe Variable
   TYPE (SpRowColD_T) :: A
   ! 
   !   gcmat komprimiert eine zeilenorientierte, dynamische
   !   Speicherstruktur einer schwachbesetzten Matrix a
   !   der Ordnung n. Die Spaltenindizes der Nichtnull-
   !   elemente der i-ten Zeile von a sind in A%ColInd(A%RowPtr(i,1)),
   !   A%ColInd(A%RowPtr(1,i)+1)...,A%ColInd(A%RowPtr(2,i)) ent-
   !   halten.
   ! 
   !    Beschreibung der Parameter
   ! 
   ! 
   !    A%n      (i/o) integer
   !                 Dimension of matrix a
   !
   !    A%RowPtr (i/o) integer(2,n)
   !                 Pointerfeld zur Beschreibung von a. A%RowPtr muss
   !                 durch das rufende Programm belegt werden.
   ! 
   !    A%ColInd (i/o) integer(len)
   !   Externe Variable
   !   TYPE (SpRowColD_T) :: A
   ! 
   !   gcmat komprimiert eine zeilenorientierte, dynamische
   !   Speicherstruktur einer schwachbesetzten Matrix a
   !   der Ordnung n. Die Spaltenindizes der Nichtnull-
   !   elemente der i-ten Zeile von a sind in A%ColInd(A%RowPtr(i,1)),
   !   A%ColInd(A%RowPtr(1,i)+1)...,A%ColInd(A%RowPtr(2,i)) ent-
   !   halten.
   ! 
   !    Beschreibung der Parameter
   ! 
   ! 
   !    A%n      (i/o) integer
   !                 Dimension of matrix a
   !
   !    A%RowPtr (i/o) integer(2,n)
   !                 Pointerfeld zur Beschreibung von a. A%RowPtr muss
   !                 durch das rufende Programm belegt werden.
   ! 
   !    A%ColInd (i/o) integer(len)
   !                 Feld zur dynamischen Verwaltung der Spaltenindizes
   !                 der Nichtnullelemente von a. 
   ! 
   !    A%ep     (i/o) integer
   !                 Pointer der auf das erste freie Feld in a verweist,
   !                 d.h A%ColInd(ep),...,A%ColInd(len) sind frei verfuegbar.
   ! 
   !    A%len    (i)   integer
   !                 Gesamtlaenge des Feldes A%ColInd. 
   ! 
   !  Interne Variable
   !
   INTEGER i,iz,j,l,pointr,rowlen,ep,len,m
   !
   m=A%m
   ep=A%ep
   len=A%len
   pointr=1
   i=1
   !
   DO 
     IF (i>=ep) EXIT
       IF (A%ColInd(i).ne.0) THEN
         !
         ! Ermittlung der aktuellen Zeile sowie deren Laenge
         DO l=1,m
           IF (A%RowPtr(1,l).le.i.and.i.le.A%RowPtr(2,l)) THEN
             iz=l
           END IF
         END DO
         rowlen=A%RowPtr(2,iz)-A%RowPtr(1,iz)
         ! 
         ! Setzen der neuen Anfangs- und Endadresse der
         ! aktuellen Zeile
         A%RowPtr(1,iz)=pointr
         A%RowPtr(2,iz)=pointr+rowlen
         DO j=pointr,pointr+rowlen
           A%ColInd(j)=A%ColInd(i)
           i=i+1
         END DO
         i=i-1
         pointr=A%RowPtr(2,iz)+1
       ENDIF
      i=i+1
    END DO
    !
    !   Belegung des freien Teils von A%ColInd mit 0
    !
    ep=pointr
    DO i=1,m
      IF (A%RowPtr(1,i).gt.ep) THEN
        A%RowPtr(1,i)=ep
        A%RowPtr(2,i)=A%RowPtr(1,i)-1
        ep=ep+inilen
      END IF
    END DO
    !        
    DO i=pointr,len
      A%ColInd(i)=0
    END DO
    A%ep=ep
  END SUBROUTINE gcmat_SpRowColD
  !
  ! Transpose Matrix
  SUBROUTINE TransposeSparse(MatAT,MatA)
    TYPE(CSR_Matrix_T), INTENT(OUT) :: MatAT
    TYPE(CSR_Matrix_T), INTENT(IN)  :: MatA
    !
    INTEGER :: i,j                                 ! zählvariabl schleIFen
    INTEGER :: indx
    !
    !
    MatAT = New_CSR(MatA%n,MatA%m)
    !
    !
    DO i=1,MatA%m
      DO j=MatA%RowPtr(i),MatA%RowPtr(i+1)-1
        MatAT%RowPtr(MatA%ColInd(j)+1)=MatAT%RowPtr(MatA%ColInd(j)+1)+1
      END DO
    END DO
    !
    DO i=1,MatAT%m
      !print*, 'RowPtr=',MatAT%RowPtr(i)
      MatAT%RowPtr(i+1)=MatAT%RowPtr(i)+MatAT%RowPtr(i+1)
    END DO
      !print*, 'RowPtr=',MatAT%RowPtr(MatAT%m+1)
    !
    ALLOCATE(MatAT%ColInd(SIZE(MatA%ColInd)))
    ALLOCATE(MatAT%Val(SIZE(MatA%Val)))
    MatAT%ColInd=0
    MatAT%Val=ZERO
    DO i=1,MatA%m
      DO j=MatA%RowPtr(i),MatA%RowPtr(i+1)-1
        indx=MatA%ColInd(j)
        MatAT%ColInd(MatAT%RowPtr(indx))=i
        MatAT%Val(MatAT%RowPtr(indx))=MatA%Val(j)
        MatAT%RowPtr(indx)=MatAT%RowPtr(indx)+1
      END DO
    END DO
    DO i=MatAT%m,1,-1
      MatAT%RowPtr(i+1)=MatAT%RowPtr(i)
    END DO
    MatAT%RowPtr(1)=1
    MatAT%nnz=MatA%RowPtr(MatA%m+1)-1
  END SUBROUTINE TransposeSparse
  !
  !
  ! SYMBOLIC MATRIX * MATRIX
  SUBROUTINE SymbolicMult(A,B,C)
    ! A*B=C
    TYPE(CSR_Matrix_T), INTENT(IN) :: A
    TYPE(CSR_Matrix_T), INTENT(IN) :: B
    TYPE(CSR_Matrix_T), INTENT(OUT) :: C
    !
    INTEGER ::  indx(MAX(A%n,A%m,B%n))
    INTEGER :: i, j, jj, k
    INTEGER :: istart, length, iTemp
    !
      !symbolic matrix multiply c=a*b
    C = New_CSR(A%m,B%n)
    !
    !main loop
    indx=0
    DO i=1,A%m
      iStart=-1
      Length=0
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        j=A%ColInd(jj)
        DO k=B%RowPtr(j),B%RowPtr(j+1)-1
          IF (indx(B%ColInd(k))==0) THEN
            indx(B%ColInd(k))=istart
            istart=B%ColInd(k)
            length=length+1
          END IF
        END DO
      END DO
      C%RowPtr(i+1)=C%RowPtr(i)+Length
      !
      DO j=C%RowPtr(i),C%RowPtr(i+1)-1
        iTemp=iStart
        istart=indx(istart)
        indx(iTemp)=0
      END DO
      indx(i) = 0
    END DO
    !==========================================
    ALLOCATE(C%ColInd(C%RowPtr(C%m+1)-1))
    C%ColInd=0
    !
    DO i=1,A%m
      iStart=-1
      Length=0
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        j=A%ColInd(jj)
        DO k=B%RowPtr(j),B%RowPtr(j+1)-1
          IF (indx(B%ColInd(k))==0) THEN
            indx(B%ColInd(k))=istart
            istart=B%ColInd(k)
            length=length+1
          END IF
        END DO
      END DO
      C%RowPtr(i+1)=C%RowPtr(i)+length
      DO j=C%RowPtr(i),C%RowPtr(i+1)-1
       C%ColInd(j)=istart
       istart=indx(istart)
       indx(C%ColInd(j))=0
      END DO
      indx(i) = 0
    END DO
    !
    DO i=1,C%m
      CALL Sort(C%Colind(C%RowPtr(i):C%RowPtr(i+1)-1))
    END DO
    ALLOCATE(C%Val(C%RowPtr(C%m+1)-1))
    C%Val=ZERO
    C%nnz=C%RowPtr(C%m+1)-1
  END SUBROUTINE SymbolicMult
  !
  ! SYMBOLIC MATRIX + MATRIX  -->  MatC = MatA + MatB
  SUBROUTINE SymbolicAdd(MatC,MatA,MatB)
    TYPE(CSR_Matrix_T), INTENT(IN)  :: MatA
    TYPE(CSR_Matrix_T), INTENT(IN)  :: MatB
    TYPE(CSR_Matrix_T), INTENT(OUT) :: MatC
    !
    INTEGER :: i,ii,j,jj,k,kk
    INTEGER, ALLOCATABLE :: TmpCol(:)
    INTEGER, ALLOCATABLE :: PermVec(:)
    INTEGER, ALLOCATABLE :: Indizes(:)
    INTEGER :: ColLen
    INTEGER :: currentlength
    INTEGER :: sameCnt
    !
    IF (.NOT.((MatA%m==MatB%m).AND.(MatA%n==MatB%n))) THEN
      WRITE(*,*) 'Wrong Matrix Dim'
      WRITE(*,*) 'A: ',MatA%m,MatA%n
      WRITE(*,*) 'B: ',MatB%m,MatB%n
      STOP 'STOP SymbolicAdd'
    END IF
    !
    MatC = New_CSR(MatA%m,MatA%n)
    !
    DO i=1,MatC%m
      currentlength=(MatA%RowPtr(i+1)-MatA%RowPtr(i))+(MatB%RowPtr(i+1)-MatB%RowPtr(i))
      ALLOCATE(Indizes(currentlength))
      Indizes=0
      sameCnt=0
      DO ii=MatA%RowPtr(i),MatA%RowPtr(i+1)-1
        DO j=MatB%RowPtr(i),MatB%RowPtr(i+1)-1
          IF (MatA%ColInd(ii)==MatB%ColInd(j)) THEN
            sameCnt=sameCnt+1
            Indizes(sameCnt)=j
          END IF
        END DO
      END DO
      MatC%RowPtr(i+1)=MatC%RowPtr(i)+currentlength-sameCnt
      DEALLOCATE(Indizes)
    END DO
    !
    ! Allocate ColInd
    ALLOCATE(MatC%ColInd(MatC%RowPtr(MatC%m+1)-1))
    MatC%ColInd=0
    !
    kk=1
    DO i=1,MatC%m
      k=1
      currentlength=MatC%RowPtr(i+1)-MatC%RowPtr(i)
      sameCnt=0
      DO ii=MatA%RowPtr(i),MatA%RowPtr(i+1)-1
        DO j=MatB%RowPtr(i),MatB%RowPtr(i+1)-1
          IF (MatA%ColInd(ii)==MatB%ColInd(j)) THEN
            sameCnt=sameCnt+1
          END IF
        END DO
      END DO
      !
      ALLOCATE(TmpCol(currentlength+sameCnt))
      TmpCol=0
      ALLOCATE(PermVec(currentlength+sameCnt))
      PermVec=0
      !
      DO jj=MatA%RowPtr(i),MatA%RowPtr(i+1)-1
        TmpCol(k)=MatA%ColInd(jj)
        k=k+1
      END DO
      DO jj=MatB%RowPtr(i),MatB%RowPtr(i+1)-1
        TmpCol(k)=MatB%ColInd(jj)
        k=k+1
      END DO
      !
      CALL unirnk(TmpCol,PermVec,ColLen)
      DO j=1,ColLen
        MatC%ColInd(kk)=TmpCol(PermVec(j))
        kk=kk+1
      END DO
      !
      DEALLOCATE(TmpCol)
      DEALLOCATE(PermVec)
    END DO
    !
    ALLOCATE(MatC%Val(MatC%RowPtr(MatC%m+1)-1))
    MatC%Val=ZERO
    MatC%nnz=MatC%RowPtr(MatC%m+1)-1
  END SUBROUTINE SymbolicAdd
  !
  ! SPARSE  MATRIX + MATRIX  -->  MatC = MatA (Sub_Add) MatB 
  SUBROUTINE SparseAdd(MatC,MatA,MatB,Sub)
    TYPE(CSR_Matrix_T), INTENT(IN)    :: MatA, MatB
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: MatC
    CHARACTER, OPTIONAL,  INTENT(IN)    :: Sub
    !
    ! Temp variables
    REAL(dp), ALLOCATABLE :: tmpBVal(:)
    INTEGER :: i, ii, jj ,kk
    !
    ALLOCATE(tmpBVal(MatB%RowPtr(MatB%m+1)-1))
    IF (PRESENT(Sub)) THEN
      IF (Sub=='-') THEN
        tmpBVal = -MatB%Val
      ELSE
        tmpBVal = MatB%Val
      END IF
    END IF
    !
    DO i=1,MatC%m
      DO ii=MatC%RowPtr(i),MatC%RowPtr(i+1)-1
        DO jj=MatA%RowPtr(i),MatA%RowPtr(i+1)-1
          IF (MatA%ColInd(jj)==MatC%ColInd(ii)) THEN
            MatC%Val(ii) = MatC%Val(ii)+MatA%Val(jj)
          END IF
        END DO
        DO kk=MatB%RowPtr(i),MatB%RowPtr(i+1)-1
          IF (MatB%ColInd(kk)==MatC%ColInd(ii)) THEN
            MatC%Val(ii) = MatC%Val(ii)+tmpBVal(kk)
          END IF
        END DO
      END DO
    END DO
    DEALLOCATE(tmpBVal)
  END SUBROUTINE SparseAdd
  !
  ! SORT ALGORITHM FOR SYMBOLIC MATRIX*MATRIX CALC
  SUBROUTINE Sort(v)
    INTEGER :: v(:)
    !
    INTEGER :: i,j,temp
    DO i=1,SIZE(v)
      DO j=i+1,SIZE(v)
        IF (v(i)>v(j)) THEN
          temp=v(i)
          v(i)=v(j)
          v(j)=temp
        END IF
      END DO
    END DO
  END SUBROUTINE Sort
  !
  !
  SUBROUTINE CompressList(ColInd,Val,Type,Name)
    INTEGER,      ALLOCATABLE :: ColInd(:)
    REAL(dp),     ALLOCATABLE, OPTIONAL :: Val(:)
    CHARACTER(*), ALLOCATABLE, OPTIONAL :: Type(:), Name(:)
    !
    INTEGER :: i,j,iList,MemberCol
    INTEGER :: TempListCol(SIZE(ColInd))
    REAL(dp) :: MemberVal
    REAL(dp) :: TempListVal(SIZE(ColInd))
    LOGICAL :: Insert

    CHARACTER(100) :: MemberName
    CHARACTER(10)  :: MemberType
    CHARACTER(100) :: TempListName(SIZE(ColInd))
    CHARACTER(10)  :: TempListType(SIZE(ColInd))
    !
    TempListVal=ZERO
    iList=0
    !
    S1:DO i=1,SIZE(ColInd)
      MemberCol=ColInd(i)
      IF (PRESENT(Val)) MemberVal=Val(i)
      IF (PRESENT(Name)) MemberName=Name(i)
      IF (PRESENT(Type)) MemberType=Type(i)
      Insert=.TRUE.
      S2:DO j=1,iList
        IF (MemberCol==TempListCol(j)) THEN
          Insert=.FALSE.
          IF (PRESENT(Val)) TempListVal(iList)=TempListVal(iList)+MemberVal
          EXIT S2
        END IF
      END DO S2
      IF (Insert) THEN
        iList=iList+1
        TempListCol(iList)=MemberCol
        IF (PRESENT(Val)) TempListVal(iList)=MemberVal
        IF (PRESENT(Name)) TempListName(iList)=MemberName
        IF (PRESENT(Type)) TempListType(iList)=MemberType
      END IF
    END DO S1
    DEALLOCATE(ColInd)
    IF (PRESENT(Val)) DEALLOCATE(Val)
    IF (PRESENT(Name)) DEALLOCATE(Name)
    IF (PRESENT(Type)) DEALLOCATE(Type)
    ALLOCATE(ColInd(1:iList))
    IF (PRESENT(Val)) ALLOCATE(Val(1:iList))
    IF (PRESENT(Name)) ALLOCATE(Name(1:iList))
    IF (PRESENT(Type)) ALLOCATE(Type(1:iList))
    ColInd=TempListCol(1:iList)
    IF (PRESENT(Val)) Val=TempListVal(1:iList)
    IF (PRESENT(Name)) Name=TempListName(1:iList)(:)
    IF (PRESENT(Type)) Type=TempListType(1:iList)(:)
  END SUBROUTINE CompressList
  !
  !
  ! SPARSE JACOBIMATRIX CALC
  PURE SUBROUTINE Jacobian_CC(JacCC,gMat,aMat,rVec,yVec)
    !
    ! jMat = gMat*Dr*aMat*invDy;
    !
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: JacCC
    TYPE(CSR_Matrix_T), INTENT(IN)    :: gMat
    TYPE(CSR_Matrix_T), INTENT(IN)    :: aMat
    REAL(dp),           INTENT(IN)    :: rVec(aMat%m)
    REAL(dp),           INTENT(IN)    :: yVec(aMat%n)
    !
    INTEGER :: i,j,jj,k,kk
    REAL(dp) :: ajj
    REAL(dp) :: temp(MAX(gMat%m,gMat%n,aMat%n))
    !
    temp=ZERO
    !
    JacCC%Val=ZERO
    !
    DO i=1,gMat%m

      DO jj=gMat%RowPtr(i),gMat%RowPtr(i+1)-1
        j   = gMat%ColInd(jj)
        ajj = gMat%Val(jj)*rVec(j)

        DO kk=aMat%RowPtr(j),aMat%RowPtr(j+1)-1
          k = aMat%ColInd(kk)
          temp(k) = temp(k) + ajj*aMat%Val(kk)/yVec(k)
        END DO

      END DO

      DO jj=JacCC%RowPtr(i),JacCC%RowPtr(i+1)-1
        j = JacCC%ColInd(jj)
        JacCC%Val(jj) = temp(j)
        temp(j) = ZERO
      END DO

    END DO
  END SUBROUTINE Jacobian_CC

  ! VAL PTR version for multiple droplet classes
  PURE SUBROUTINE Jacobian_CC_ValPtr(JacCC,gMat,aMat,rVec,yVec)
    !
    ! jMat = gMat*Dr*aMat*invDy;
    !
    USE Reac_Mod, ONLY: nspc2, nreac2
    TYPE(CSR_Matrix_T),    INTENT(INOUT) :: JacCC
    TYPE(CSR_Matrix_T),    INTENT(IN)    :: gMat
    TYPE(CSR_Matrix_T),    INTENT(IN)    :: aMat
    REAL(dp),              INTENT(IN)    :: rVec(nreac2)
    REAL(dp),              INTENT(IN)    :: yVec(nspc2)
    !
    INTEGER :: i,j,jj,k,kk,j1,j2
    REAL(dp) :: ajj(nDropletClasses)
    REAL(dp) :: temp(aMat%n*nDropletClasses)
    !
    temp=ZERO
    !
    JacCC%Val=ZERO
    !
    DO i = 1 , gmat%m
      DO jj = gMat%RowPtr(i), gMat%RowPtr(i+1)-1
        j  = gMat%ColInd(jj)
        j1 = nD_Ptr_reacs(j)
        j2 = nD_Ptr_reacs(j+1)

        ajj(1:j2-j1) = gMat%Val(jj)*rVec(j1:j2-1)

        DO kk = aMat%RowPtr(j), aMat%RowPtr(j+1)-1
          k = aMat%ColInd(kk)
          ! we are at: 
          ! the dependence of species i on species k due to reaction j
          ! NOTE: if species i is nD, then every entry is nD !

          IF ( nD_spc(i) .AND. nD_spc(k) ) THEN
            ! aqueous daq/daq
            temp(1+(k-1)*nDropletClasses:k*nDropletClasses) = temp(1+(k-1)*nDropletClasses:k*nDropletClasses) + &
                                  & ajj(1:nDropletClasses) * aMat%Val(kk) / yVec(nD_Ptr_spc(k):nD_Ptr_spc(k+1)-1)
          ELSE IF ( nD_spc(i) .AND. .NOT. nD_spc(k) ) THEN
            ! henry daq/dgas
            temp(1+(k-1)*nDropletClasses:k*nDropletClasses) = temp(1+(k-1)*nDropletClasses:k*nDropletClasses) + &
                                  & ajj(1:nDropletClasses) * aMat%Val(kk) / yVec(nD_Ptr_spc(k))
          ELSE IF ( (.NOT. nD_spc(i)) .AND. nD_spc(k) ) THEN
            ! henry dgas/daq
            temp(nD_Ptr_spc(k):nD_Ptr_spc(k+1)-1) = temp(nD_Ptr_spc(k):nD_Ptr_spc(k+1)-1) + &
                                                  & ajj(1:nDropletClasses) * aMat%Val(kk) / yVec(nD_Ptr_spc(k):nD_Ptr_spc(k+1)-1)
          ELSE
            ! gaseous dgas/dgas
            ! (if j is a henry reac, ajj is still a vector)
            temp(nD_Ptr_spc(k)) = temp(nD_Ptr_spc(k)) + SUM(ajj(1:j2-j1))*aMat%Val(kk)/yVec(nD_Ptr_spc(k))
          END IF
        END DO
      END DO

      DO jj = JacCC%RowPtr(i), JacCC%RowPtr(i+1)-1
        j  = JacCC%ColInd(jj)
        j1 = JacCC%ValPtr(jj)
        j2 = JacCC%ValPtr(jj+1)
        IF ( nD_spc(i) ) THEN
          JacCC%Val(j1:j2-1) = temp(1+(j-1)*nDropletClasses:j*nDropletClasses)
          temp(1+(j-1)*nDropletClasses:j*nDropletClasses) = ZERO
        ELSE
          JacCC%Val(j1:j2-1) = temp(nD_Ptr_spc(j):nD_Ptr_spc(j+1)-1)
          temp(nD_Ptr_spc(j):nD_Ptr_spc(j+1)-1) = ZERO
        END IF
      END DO
    END DO
  END SUBROUTINE Jacobian_CC_ValPtr

  PURE SUBROUTINE Jacobian_parcel(Jac_parcel, dSeqdmw)
    USE Reac_Mod,    ONLY: nspc2, iAqMassEq2
    REAL(dp), INTENT(IN)  :: dSeqdmw(nDropletClasses)
    REAL(dp), INTENT(OUT) :: Jac_parcel(nDropletClasses+4)

    Jac_parcel = 0.0_dp

    Jac_parcel(iAqMassEq2-nspc2) = dSeqdmw

  END SUBROUTINE Jacobian_parcel
 
  ! JacTC = -1/cv/rho [C_v*dTdt + U^T*JacCC]
  PURE SUBROUTINE Jacobian_TC(JacTC,JacCC,cv,dUdT,dTdt,U,rRho)
    TYPE(CSR_Matrix_T), INTENT(IN)  :: JacCC
    REAL(dp),           INTENT(IN)  :: cv , dTdt, rRho
    REAL(dp),           INTENT(IN)  :: dUdT(:), U(:)
    REAL(dp),           INTENT(OUT) :: JacTC(JacCC%m)
    !
    REAL(dp) :: tmpJacVal(JacCC%m)

    !tmpJacVal = DAX_sparse(JacCC,U)
    tmpJacVal = U * JacCC
    !JacTC = - (dUdT*dTdt + tmpJacVal) / cv * rRho * milli
    JacTC = - (dUdT*dTdt + tmpJacVal) / cv * rRho

  END SUBROUTINE Jacobian_TC
  

  ! JacCT = BAT * Dr * dkdT_over_k
  PURE SUBROUTINE Jacobian_CT(JacCT,gMat,Dr,dkdT_over_k)
    TYPE(CSR_Matrix_T), INTENT(IN)  :: gMat
    REAL(dp),           INTENT(IN)  :: Dr(:), dkdT_over_k(:)
    REAL(dp),           INTENT(OUT) :: JacCT(gMat%m)
    !
    REAL(dp) :: dRatedT(gMat%n)

    ! deriv. of reaction rate with resp. to temperature
    dRatedT = Dr * dkdT_over_k
    JacCT   = DAX_sparse(gMat,dRatedT)

  END SUBROUTINE Jacobian_CT
 

  ! JacTT = -1/cv/rho [dTdT*dcvdT+C_v*dCdt + U^T*JacCC]
  PURE SUBROUTINE Jacobian_TT(JacTT,JacCT,cv,dcvdT,dTdt,dUdT,dcdt,U,rRho)
    REAL(dp), INTENT(IN)    :: JacCT(:)
    REAL(dp), INTENT(IN)    :: cv , dcvdT , dTdt , rRho
    REAL(dp), INTENT(IN)    :: dUdT(:) , dCdt(:) , U(:)
    REAL(dp), INTENT(INOUT) :: JacTT
    !
    !JacTT = - milli * rRho/cv  *                &
    !     &   (  dTdt*dcvdT/cv + SUM(dUdT*dCdt) &
    !     &                    + SUM(U*JacCT)   )
    JacTT = - rRho/cv  *                         &
         &   (  dTdt*dcvdT/rRho + SUM(dUdT*dCdt) &
         &                      + SUM(U*JacCT)   )
    
  END SUBROUTINE Jacobian_TT
 

  ! SPARSE MITER CALCULATION
  PURE SUBROUTINE Miter_Classic(Miter,h,g,J1,J2,J3,J4)
    USE Reac_Mod,    ONLY: nspc, iAqMassEq, iqEq
    REAL(dp),           INTENT(IN)    :: h, g
    TYPE(CSR_Matrix_T), INTENT(IN)    :: J1
    REAL(dp), OPTIONAL, INTENT(IN)    :: J2(:), J3(:), J4
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: Miter
    !
    INTEGER :: i,j,jj,cnt
    REAL(dp) :: hg
   
    Miter%Val = ZERO
    hg  = h*g
    cnt = 0

    DO i = 1,J1%m
      DO jj = Miter%RowPtr(i) , Miter%RowPtr(i+1)-1
        j = Miter%ColInd(jj)
        IF      ( j==i ) THEN
          ! diagonal d(dCdt)/dC
          Miter%Val(jj) = ONE - hg*J1%Val(jj-cnt)
        ELSE IF ( j==Miter%m.AND.combustion) THEN
          cnt = cnt + 1
        ELSE
          ! none diagonal  d(dCdt)/dC
          Miter%Val(jj) = - hg*J1%Val(jj-cnt)
        END IF
      END DO
    END DO

    IF (combustion) THEN
      ! bottom row d(dTdt)/dC
      Miter%Val(Miter%RowVectorPtr) = - hg*J2
      ! right hand column d(dCdt)/dT
      Miter%Val(Miter%ColVectorPtr) = - hg*J3
      ! lower right hand corner d(dTdt)/dT
      Miter%Val(Miter%XPtr) = ONE - hg*J4
    END IF
    IF (adiabatic_parcel) THEN
      ! J2 = dSeqdmw
      
      DO i=1,5
        Miter%Val(Miter%DiagPtr(J1%m+i)) = ONE - hg*J2(i)
      END DO

      ! entry at (iq,iAqMass)
      Miter%Val(Miter%RowPtr(iqEq)) = hg * (milli*J2(iAqMassEq-nspc))
    END IF

  END SUBROUTINE Miter_Classic

  PURE SUBROUTINE Miter_Classic_ValPtr(Miter,h,g,J1,J2)
    USE Reac_Mod,    ONLY: nspc2, nspc, iqEq, iAqMassEq2
    REAL(dp),           INTENT(IN)    :: h, g
    TYPE(CSR_Matrix_T), INTENT(IN)    :: J1
    REAL(dp), OPTIONAL, INTENT(IN)    :: J2(:)
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: Miter
    !
    INTEGER :: i,j,jj,cnt
    REAL(dp) :: hg

    Miter%Val = ZERO
    hg  = h*g
    cnt=0

    DO i = 1,J1%m
      DO jj = Miter%RowPtr(i) , Miter%RowPtr(i+1)-1
        j = Miter%ColInd(jj)
        IF ( j==i ) THEN
          ! diagonal d(dCdt)/dC
          Miter%Val(Miter%ValPtr(jj):Miter%ValPtr(jj+1)-1) = ONE - hg*J1%Val(J1%ValPtr(jj-cnt):J1%ValPtr(jj-cnt+1)-1)
        ELSE
          ! none diagonal  d(dCdt)/dC
          Miter%Val(Miter%ValPtr(jj):Miter%ValPtr(jj+1)-1) = - hg*J1%Val(J1%ValPtr(jj-cnt):J1%ValPtr(jj-cnt+1)-1)
        END IF
      END DO
    END DO

    IF (adiabatic_parcel) THEN
      ! J2 = dSeqdmw, J4 = rho_parcel

      DO i = 1, 5
        Miter%Val(Miter%ValPtr(Miter%DiagPtr(J1%m+i)):Miter%ValPtr(Miter%DiagPtr(J1%m+i)+1)-1) = ONE - hg * J2(nD_Ptr_spc(nspc+i)-nspc2 : nD_Ptr_spc(nspc+i+1)-nspc2-1)
      END DO
      ! entry at (iq,iAqMass)
      Miter%Val(Miter%ValPtr(Miter%RowPtr(iqEq)):Miter%ValPtr(Miter%RowPtr(iqEq)+1)-1) = - hg * (-milli*J2(iAqMassEq2-nspc2))
    END IF

  END SUBROUTINE Miter_Classic_ValPtr

  SUBROUTINE BuildSymbolicClassicMatrix(CL,Jac)
    USE Reac_Mod, ONLY: iqEq, iAqMassEq
    TYPE(CSR_Matrix_T), INTENT(OUT) :: CL
    TYPE(CSR_Matrix_T), INTENT(IN)  :: Jac

    TYPE(CSR_Matrix_T) :: Id
    TYPE(CSR_Matrix_T) :: CL0
    !
    INTEGER :: i, j, jj, ndim, nnz
    !
    !------------------------------------------------------------------------------
    ! --- Set Matrix dimensions and nonzeros 
    !------------------------------------------------------------------------------
    !
    !
    IF ( combustion ) THEN
      ndim = Jac%n + 1                  ! number of rows/coloumns
      nnz  = Jac%RowPtr(Jac%m+1)-1   &  ! nonzeros of Jacobian
      &       + Jac%m + Jac%n + 1       ! Dc,U^T and Dr,~K and X (down right)
    ELSE IF ( adiabatic_parcel ) THEN
      ndim = Jac%n + 5
      nnz  = Jac%RowPtr(Jac%m+1)-1   &  ! nonzeros of Jacobian of concentrations
         & + 5                       &  ! diagonal entries for water masses, T, q, rho, z
         & + 1                          ! entry at (iq,iWaterMasses)
    ELSE
      ndim = Jac%n                      ! number of rows/coloumns
      nnz  = Jac%RowPtr(Jac%m+1)-1      ! nonzeros of Jacobian
    END IF

    CL0 = New_CSR ( ndim , ndim , nnz )
    ID = SparseID( ndim )
    ALLOCATE(CL0%DiagPtr(ndim))
    CL0%DiagPtr  = -12

    IF ( combustion ) THEN
      !
      DO i = 1 , ndim - 1
        CL0%RowPtr(i+1) = CL0%RowPtr(i) + (Jac%RowPtr(i+1)-Jac%RowPtr(i)) + 1

        CL0%ColInd( CL0%RowPtr(i):CL0%RowPtr(i+1)-1 ) =                           & 
        &    [ Jac%ColInd(Jac%RowPtr(i):(Jac%RowPtr(i+1)-1)) , ndim ]
      END DO
      
      ! set last row (full row)
      CL0%RowPtr(ndim+1) = CL0%RowPtr(ndim) + ndim
      CL0%ColInd(CL0%RowPtr(ndim):CL0%RowPtr(ndim+1)-1) = [( i , i = 1 , ndim )]

      CALL SymbolicAdd(CL,Id,CL0)
      ALLOCATE(CL%RowVectorPtr(ndim-1))
      ALLOCATE(CL%ColVectorPtr(ndim-1))
      CL%RowVectorPtr = -1
      CL%ColVectorPtr = -1
      CL%Val          = ONE
    ELSE IF ( adiabatic_parcel ) THEN
      ! copy Jac 
      DO i = 1 , Jac%n
        CL0%RowPtr(i+1) = CL0%RowPtr(i) + (Jac%RowPtr(i+1)-Jac%RowPtr(i))

        CL0%ColInd( CL0%RowPtr(i):CL0%RowPtr(i+1)-1 ) = Jac%ColInd(Jac%RowPtr(i):(Jac%RowPtr(i+1)-1))
      END DO
      ! set diagonal entries for the additional equations only
      DO i = Jac%n+1, ndim
        CL0%RowPtr(i+1) = CL0%RowPtr(i) + 1
        IF (i==iqEq) THEN ! add one more value at iAqMass
          CL0%RowPtr(i+1) = CL0%RowPtr(i+1) + 1
          CL0%ColInd(CL0%RowPtr(i)) = iAqMassEq
        END IF
        CL0%ColInd(CL0%RowPtr(i+1)-1) = i
      END DO
      CALL SymbolicAdd(CL,Id,CL0)
      CL%Val = ONE
    ELSE
      CALL SymbolicAdd(CL,Id,Jac)
    END IF
    !
    ALLOCATE(CL%DiagPtr(ndim))
    ALLOCATE(CL%ValPtr(CL%RowPtr(SIZE(CL%RowPtr))))
    CL%ValPtr(1) = 1
    nnz = 0
    CL%DiagPtr  = -12
    ! get diagonal pointer
    DO i = 1 , ndim
      IF ( combustion .AND. i < ndim ) CL%ColVectorPtr(i)  = CL%RowPtr(i+1) - 1
      DO jj = CL%RowPtr(i) , CL%RowPtr(i+1)-1 
        j = CL%ColInd(jj)
        IF ( i == j ) CL%DiagPtr(i) = jj
        IF (nD_spc(i) .OR. nD_spc(j)) THEN
          nnz = nnz + nDropletClasses
        ELSE
          nnz = nnz + 1
        END IF
        CL%ValPtr(jj+1) = nnz+1
      END DO
    END DO
    ! re-allocate value array with proper nnz (proper=regarding multiple droplet classes)
    DEALLOCATE(CL%Val)
    ALLOCATE(CL%Val(nnz))
    CL%Val = ONE
    IF ( combustion ) THEN
      CL%RowVectorPtr = [( i , i = CL%RowPtr(nDim),CL%RowPtr(nDim+1)-2 )]
      CL%XPtr = CL%DiagPtr(ndim)
    END IF
  END SUBROUTINE BuildSymbolicClassicMatrix
  !
  !
  PURE SUBROUTINE SetLUvaluesCL(LU,A,Permu)
    !
    ! Set values to block matrix
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: LU
    TYPE(CSR_Matrix_T), INTENT(IN)    :: A
    INTEGER,            INTENT(IN)    :: Permu(:)
    !
    LU%Val = ZERO
    LU%Val( Permu ) = A%Val
    !
  END SUBROUTINE SetLUvaluesCL 
  !
  ! Permutes the values in Miter and writes it to LU structur of Miter
  ! Permutation vector is generated in this routine
  SUBROUTINE Get_LU_Permutation(Permutation,LU,A)
    ! INOUT
    TYPE(CSR_Matrix_T) :: LU
    ! IN
    TYPE(CSR_Matrix_T) :: A
    INTEGER            :: nnzA
    ! OUT
    INTEGER, ALLOCATABLE :: Permutation(:)
    ! TEMP
    INTEGER :: i,ip,j,jj,jp,jjp,jp1
 
    nnzA=A%RowPtr(A%m+1)-1
   
    IF (.NOT.ALLOCATED(Permutation)) ALLOCATE(Permutation(nnzA))
    Permutation = -14

    LU%Val = ZERO

    DO i = 1 , A%n  
      ip = LU%Permu(i)
      DO jj = A%RowPtr(i) , A%RowPtr(i+1)-1
        j  = A%ColInd(jj)
        jp = LU%Permu(A%ColInd(jj))
        DO jjp = LU%RowPtr(ip) , LU%RowPtr(ip+1)-1
          jp1 = LU%ColInd(jjP)
          IF ( jp1 == jp ) THEN
            LU%Val(jjP) = A%Val(jj)
            Permutation(jj) = jjP
          END IF  
        END DO  
      END DO  
    END DO
    ALLOCATE(A%LUperm(nnzA),LU%LUperm(nnzA))
    A%LUperm  = Permutation
    LU%LUperm = Permutation

    IF ( combustion ) THEN
      LU%ColVectorPtr = Permutation( A%ColVectorPtr )
      LU%RowVectorPtr = Permutation( A%RowVectorPtr )
      LU%XPtr         = Permutation( A%XPtr )
    END IF

  END SUBROUTINE Get_LU_Permutation
  !
  !
  SUBROUTINE Get_ValPtr_Permutations(ValPtrPermutation, w_InvColInd, bPermu, bInvPermu, &
                                   & bPtr, LU, A, LU_Perm, nD, nDIM2)
    ! ARGUMENTS
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: LU
    TYPE(CSR_Matrix_T), INTENT(IN) :: A
    INTEGER, DIMENSION(:), INTENT(IN) :: LU_Perm
    INTEGER :: nD, nDIM2

    ! OUT:
    INTEGER, ALLOCATABLE, INTENT(OUT) :: ValPtrPermutation(:), w_InvColInd(:,:,:), &
                                       & bPermu(:), bInvPermu(:), bPtr(:)

    ! INTERNAL VARIABLES
    INTEGER :: nnzA, i, j, jj, cnt, jjp, nnzLU, nval
    INTEGER, ALLOCATABLE :: LU_nVals(:)

    nnzA = A%ValPtr(SIZE(A%ValPtr))-1
    nnzLU = LU%RowPtr(SIZE(LU%RowPtr))-1

    ALLOCATE(ValPtrPermutation(nnzA), LU_nVals(nnzLU))
    ValPtrPermutation = -14
    LU_nVals = -14

    DO i = 1 , LU%n
      DO jj = LU%RowPtr(i), LU%RowPtr(i+1)-1
        j = LU%ColInd(jj)
        IF ( nD_spc(LU%InvPer(i)) .OR. nD_spc(LU%InvPer(j)) ) THEN
          LU_nVals(jj) = nD
        ELSE
          LU_nVals(jj) = 1
        END IF
      END DO
    END DO

    ! make LU.Val array proper size (to this point, it has compressed size from RowColD_to_CSR)
    IF (ALLOCATED(LU%Val)) DEALLOCATE(LU%Val)
    ALLOCATE(LU%Val(SUM(LU_nVals)))

    ! write ValPtr array for LU
    IF (ALLOCATED(LU%ValPtr)) DEALLOCATE(LU%ValPtr)
    ALLOCATE(LU%ValPtr(nnzLU+1))
    LU%ValPtr = 0
    LU%ValPtr(1) = 1
    DO i = 2 , SIZE(LU%ValPtr)
      LU%ValPtr(i) = SUM(LU_nVals(1:i-1)) + 1
    END DO

    cnt = 1
    DO jj = 1 , SIZE(LU_Perm)
      jjp = LU_Perm(jj)
      DO i = LU%ValPtr(jjp), LU%ValPtr(jjp+1)-1
        ValPtrPermutation(cnt) = i
        cnt = cnt+1
      END DO
    END DO

    ALLOCATE(w_InvColInd(LU%n, 2, LU%n))
    w_InvColInd = 0

    DO i=1,LU%n
      cnt = 1
      DO jj=LU%RowPtr(i), LU%RowPtr(i+1)-1
        nval = LU%ValPtr(jj+1) - LU%ValPtr(jj)
        w_InvColInd(i, 1, LU%ColInd(jj)) = cnt
        w_InvColInd(i, 2, LU%ColInd(jj)) = nval
        cnt = cnt + nval
      END DO
    END DO

    ALLOCATE(bPermu(nDIM2), bInvPermu(nDIM2), bPtr(LU%n+1))
    cnt = 1
    bPtr(1) = 1
    DO i = 1 , LU%n
      !IF ( LU%InvPer(i)<=nAq ) THEN
      IF ( nD_spc(LU%InvPer(i)) ) THEN
        DO j=1,nD
          bPermu(nD_Ptr_spc(LU%InvPer(i))+j-1) = cnt
          bInvPermu(cnt) = nD_Ptr_spc(LU%InvPer(i))+j-1
          cnt=cnt+1
        END DO
        bPtr(i+1) = cnt
      ELSE
        bPermu(nD_Ptr_spc(LU%InvPer(i))) = cnt
        bInvPermu(cnt) = nD_Ptr_spc(LU%InvPer(i))
        cnt=cnt+1
        bPtr(i+1) = cnt
      END IF
    END DO

  END SUBROUTINE Get_ValPtr_Permutations
  !
  !
  SUBROUTINE WriteSparseMatrix(A,FileName,nr,ns)
    TYPE(CSR_Matrix_T), INTENT(IN) :: A
    CHARACTER(*), OPTIONAL :: FileName
    INTEGER, OPTIONAL :: nr,ns
    !
    INTEGER :: i,j,jj
    !

    ! 
    OPEN(UNIT=99,FILE=ADJUSTL(TRIM(FileName))//'.SparseMat',STATUS='UNKNOWN')

    WRITE(99,*) '###########################################################'
    WRITE(99,*) '##############  Sparse Matrix  Matlab input ###############'
    WRITE(99,*) '###########################################################'
    WRITE(99,*) 
    WRITE(99,*) '###########################################################'
    WRITE(99,*) '#     Name:        ' , ADJUSTL(FileName)
    WRITE(99,*) '#'  
    WRITE(99,*) '#     Dimension:   ' , A%m,' x ',A%n
    WRITE(99,*) '#     Nonzeros:    ' , A%nnz
    WRITE(99,*) '#     nreac, nspc: ' , nr,' , ',ns
    WRITE(99,*) '###########################################################'
    WRITE(99,*)

    WRITE(99,*) '###########################################################'
    WRITE(99,*) '# Sparse Matrix in row/column - index format'
    WRITE(99,*) '###########################################################'
    WRITE(99,*) 'MATRIX'
    DO i=1,A%m
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        j=A%ColInd(jj)
        WRITE(99,'(1X,I12,10X,I12,10X,E23.14E4)') i,j,A%Val(jj)
      END DO
    END DO
    WRITE(99,*)
    WRITE(99,*)

    IF (ALLOCATED(A%DiagPtr)) THEN
      WRITE(99,*) '###########################################################'
      WRITE(99,*) '# Array for diagonal entries'
      WRITE(99,*) '###########################################################'
      WRITE(99,*) 'DIAG_PTR'
      DO i=1,SIZE(A%DiagPtr)
        WRITE(99,'(1X,I6,10X,I6)') i,A%DiagPtr(i)
      END DO
      WRITE(99,*)
      WRITE(99,*)
    END IF

    IF (ALLOCATED(A%Permu)) THEN
      WRITE(99,*) '##############################################################'
      WRITE(99,*) '# Array for LU permutation and inverse permutation (symmetric)'
      WRITE(99,*) '##############################################################'
      WRITE(99,*) 'DIAG_PERMUTATION'
      DO i=1,SIZE(A%Permu)
        WRITE(99,'(1X,I6,10X,I6,10X,I6)') i,A%Permu(i), A%InvPer(i)
      END DO
      WRITE(99,*)
      WRITE(99,*)
    END IF

    IF (ALLOCATED(A%LUperm)) THEN
      WRITE(99,*) '###########################################################'
      WRITE(99,*) '# Array for LU permutation of all nonzeros '
      WRITE(99,*) '###########################################################'
      WRITE(99,*) 'LU_PERMUTATION'
      DO i=1,SIZE(A%LUperm)
        WRITE(99,'(1X,I6,10X,I6)') i,A%LUperm(i)
      END DO
      WRITE(99,*)
      WRITE(99,*)
    END IF

    IF (combustion) THEN
      IF (ALLOCATED(A%RowVectorPtr)) THEN
        WRITE(99,*) '###########################################################'
        WRITE(99,*) '# Array for indices pointing to the row vector             '
        WRITE(99,*) '#                        (only with combustion) '
        WRITE(99,*) '###########################################################'
        WRITE(99,*) 'ROW_PTR'
        DO i=1,SIZE(A%RowVectorPtr)
          WRITE(99,'(1X,I12,10X,I12)') i,A%RowVectorPtr(i)
        END DO
        WRITE(99,*)
        WRITE(99,*)
      END IF
      IF (ALLOCATED(A%ColVectorPtr)) THEN
        WRITE(99,*) '###########################################################'
        WRITE(99,*) '# Array for indices pointing to the column vector          '
        WRITE(99,*) '#                        (only with combustion) '
        WRITE(99,*) '###########################################################'
        WRITE(99,*) 'COL_PTR'
        DO i=1,SIZE(A%ColVectorPtr)
          WRITE(99,'(1X,I12,10X,I12)') i,A%ColVectorPtr(i)
        END DO
        WRITE(99,*)
        WRITE(99,*)
      END IF
      IF (A%XPtr/=-42) THEN
        WRITE(99,*) '###########################################################'
        WRITE(99,*) '# Index pointing to the X value           '
        WRITE(99,*) '#                        (only with combustion) '
        WRITE(99,*) '###########################################################'
        WRITE(99,*) 'X_PTR'
        WRITE(99,'(1X,I12,10X,I12)') i,A%XPtr
        WRITE(99,*)
        WRITE(99,*)
      ELSE
        WRITE(99,*) 'X_PTR is not definded'
      END IF
    END IF

    CLOSE(99)
    WRITE(*,*) '  Writing matrices to file: ',TRIM(FileName)//'.SparseMat'
  END SUBROUTINE WriteSparseMatrix
  !
  ! sparse vector * matrix
  PURE FUNCTION DAX_T_sparse(X,A) RESULT(Rhs)
    TYPE(CSR_Matrix_T), INTENT(IN) :: A
    REAL(dp),           INTENT(IN) :: X(A%m)
    REAL(dp)                       :: Rhs(A%n)
    !TEMP
    INTEGER :: i, j, jj  ! RowPtr(i), RowPtr(i+1)-1

    Rhs = ZERO
    DO i=1,A%m
      DO jj=A%RowPtr(i), A%RowPtr(i+1)-1
        j = A%ColInd(jj)
        Rhs(j) = Rhs(j) + X(i) * A%Val(jj)
      END DO
    END DO

  END FUNCTION DAX_T_sparse

  ! sparse matrix * vector
  PURE FUNCTION DAX_sparse(A,X) RESULT(Rhs)
    TYPE(CSR_Matrix_T), INTENT(IN) :: A
    REAL(dp),           INTENT(IN) :: X(:)
    REAL(dp)                       :: Rhs(A%m)
    !TEMP
    INTEGER :: i, rp_i, rp_i1, jj  ! RowPtr(i), RowPtr(i+1)-1

    Rhs = ZERO
    DO i=1,A%m
      rp_i   = A%RowPtr(i)  ;  rp_i1  = A%RowPtr(i+1)-1
      !Rhs(i) = SUM(A%Val(rp_i:rp_i1)*X(A%ColInd(rp_i:rp_i1)))
      DO jj=rp_i, rp_i1
        Rhs(i) = Rhs(i) + A%Val(jj) * X(A%ColInd(jj))
      END DO
    END DO

  END FUNCTION DAX_sparse

  ! normal sparse matrix (nSpc x nReacs) * (nReacs_nD) vector-valued reaction vector
  PURE FUNCTION MULT_BAT_Rate_ValPtr(BAT, Rate) RESULT(Rhs)

    USE Control_Mod, ONLY: nD_spc, nD_reac, nD_Ptr_Spc, nD_Ptr_reacs
    USE Reac_Mod,    ONLY: nspc2

    TYPE(CSR_Matrix_T), INTENT(IN) :: BAT
    REAL(dp),           INTENT(IN) :: Rate(:)
    REAL(dp)                       :: Rhs(nspc2)
    !TEMP
    INTEGER :: i, j, jj
   
    Rhs = ZERO
    DO i=1,BAT%m
      DO jj = BAT%RowPtr(i) , BAT%RowPtr(i+1)-1
        j = BAT%ColInd(jj)
        IF ( nD_reac(j) .EQV. nD_spc(i) ) THEN
          ! both nD or both not nD 
          ! (aqua -> aqua or gas -> gas or gas -> aqua ) 
          Rhs(nD_Ptr_spc(i):nD_Ptr_spc(i+1)-1) = Rhs(nD_Ptr_spc(i):nD_Ptr_spc(i+1)-1) &
                                             & + BAT%Val(jj) * Rate(nD_Ptr_reacs(j):nD_Ptr_reacs(j+1)-1)
        ELSE
          ! henry reac aqua -> gas
          Rhs(nD_Ptr_spc(i)) = Rhs(nD_Ptr_spc(i)) &
                           & + BAT%Val(jj) * SUM(Rate(nD_Ptr_reacs(j):nD_Ptr_reacs(j+1)-1))
        END IF
      END DO
    END DO
  END FUNCTION MULT_BAT_Rate_ValPtr
 
  ! multiplies ValPtr matrix with fitting vector-valued-vector (matrix quadratic)
  PURE FUNCTION MATMUL_nD(A, X) RESULT(Rhs)
    USE Reac_Mod, ONLY: nspc2
    TYPE(CSR_Matrix_T), INTENT(IN) :: A
    REAL(dp),           INTENT(IN) :: X(:)
    REAL(dp)                       :: Rhs(nspc2)
    !TEMP
    INTEGER :: i, j, jj
    
    Rhs = ZERO

    DO i=1,A%m
      DO jj = A%RowPtr(i) , A%RowPtr(i+1)-1
        j = A%ColInd(jj)
        IF (nD_spc(i) .AND. nD_spc(j)) THEN
          Rhs(nD_Ptr_spc(i):nD_Ptr_spc(i+1)-1) = Rhs(nD_Ptr_spc(i):nD_Ptr_spc(i+1)-1) &
                                             & + A%Val(A%ValPtr(jj):A%ValPtr(jj+1)-1) * X(nD_Ptr_spc(j):nD_Ptr_spc(j+1)-1)
        ELSE IF (nD_spc(i) .AND. .NOT. nD_spc(j)) THEN
          Rhs(nD_Ptr_spc(i):nD_Ptr_spc(i+1)-1) = Rhs(nD_Ptr_spc(i):nD_Ptr_spc(i+1)-1) &
                                             & + A%Val(A%ValPtr(jj):A%ValPtr(jj+1)-1) * X(nD_Ptr_spc(j))
        ELSE
          Rhs(nD_Ptr_spc(i)) = Rhs(nD_Ptr_spc(i)) &
                           & + SUM(A%Val(A%ValPtr(jj):A%ValPtr(jj+1)-1) * X(nD_Ptr_spc(j):nD_Ptr_spc(j+1)-1))
        END IF
      END DO
    END DO
  END FUNCTION MATMUL_nD
  !
  !
  PURE SUBROUTINE SparseLU(A)
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: A
    !
    REAL(dp) :: w(A%n)
    REAL(dp) :: alpha
    INTEGER :: i,j,jj,kk
    !INTEGER :: OPCOUNT

    !OPCOUNT=0
    !
    DO i=1,A%n
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        w(A%ColInd(jj))=A%Val(jj)
      END DO
      DO jj=A%RowPtr(i),A%DiagPtr(i)-1
        j=A%ColInd(jj)
        alpha=w(j)/A%Val(A%DiagPtr(j))
        !OPCOUNT=OPCOUNT+1
        w(j)=alpha
        DO kk=A%DiagPtr(j)+1,A%RowPtr(j+1)-1
          w(A%ColInd(kk))=w(A%ColInd(kk))-alpha*A%Val(kk)
          !OPCOUNT=OPCount+1
        END DO
      END DO
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        A%Val(jj)=w(A%ColInd(jj))
      END DO
    END DO
    !print*, 'opcount=', opcount
    !stop
  END SUBROUTINE SparseLU
  !
  !
  PURE SUBROUTINE SparseLU_ValPtr(A, w_guide)
  TYPE(CSR_Matrix_T), INTENT(INOUT) :: A
  INTEGER, DIMENSION(A%n,2,A%n), INTENT(IN) :: w_guide
  !
  ! column may have only vector entries
  REAL(dp) :: w(nDropletClasses*A%n)
  ! alpha may have nD entries or one, we don't distinguish
  REAL(dp) :: alpha(nDropletClasses)
  INTEGER :: i,j,jj,kk,cnt,k,j1,j2,k1,k2,nval,n_alpha

  DO i = 1 , A%n
    cnt=0
    DO jj = A%RowPtr(i), A%RowPtr(i+1)-1
      cnt  = w_guide(i, 1, A%ColInd(jj))
      nval = w_guide(i, 2, A%ColInd(jj))
      w(cnt:cnt+nval-1) = A%Val(A%ValPtr(jj):A%ValPtr(jj+1)-1)
    END DO

    DO jj = A%RowPtr(i), A%DiagPtr(i)-1
      j    = A%ColInd(jj)
      j1   = w_guide(i, 1, j)
      nval = w_guide(i, 2, j)
      j2   = j1 + nval-1

      alpha(1:nval) = w(j1:j2) / A%Val(A%ValPtr(A%DiagPtr(j)):A%ValPtr(A%DiagPtr(j)+1)-1)
      w(j1:j2) = alpha(1:nval)
      n_alpha = nval
      DO kk = A%DiagPtr(j)+1 , A%RowPtr(j+1)-1
        k    = A%ColInd(kk)
        k1   = w_guide(i, 1, k)
        nval = w_guide(i, 2, k)
        k2   = k1 + nval-1

        CALL SubMult(w(k1:k2), alpha(1:n_alpha), A%Val(A%ValPtr(kk):A%ValPtr(kk+1)-1), nval, n_alpha, A%ValPtr(kk+1)-A%ValPtr(kk))
      END DO
    END DO

    DO jj = A%RowPtr(i) , A%RowPtr(i+1)-1
      j1   = w_guide(i, 1, A%ColInd(jj))
      nval = w_guide(i, 2, A%ColInd(jj))
      j2   = j1 + nval-1

      A%Val(A%ValPtr(jj):A%ValPtr(jj+1)-1) = w(j1:j2)
    END DO
  END DO
  END SUBROUTINE SparseLU_ValPtr
  !
  !
  PURE SUBROUTINE SolveSparse(LU,rhs)
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: LU
    REAL(dp),           INTENT(INOUT) :: Rhs(:)
    !
    INTEGER :: i,jj
    REAL(dp) :: b(LU%n)

   
    DO i = 1, LU%n
      b(LU%Permu(i)) = Rhs(i)
    END DO

    !--  L-solve
    DO i=2,LU%n
      DO jj=LU%RowPtr(i),LU%DiagPtr(i)-1
        b(i)=b(i)-LU%Val(jj)*b(LU%ColInd(jj))
      END DO
    END DO

    !--  U-solve
    DO i=LU%n,1,-1
      DO jj=LU%DiagPtr(i)+1,LU%RowPtr(i+1)-1
        b(i)=b(i)-LU%Val(jj)*b(LU%ColInd(jj))
      END DO

      b(i)=b(i)/LU%Val(LU%DiagPtr(i))
    END DO

    !--  Back-Permutation of solution
    DO i = 1, LU%n
      Rhs(LU%InvPer(i)) = b(i)
    END DO
  END SUBROUTINE SolveSparse
  !
  !
  PURE SUBROUTINE SolveSparse_ValPtr(LU, rhs, bPermu, bInvPermu, bPtr)
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: LU
    REAL(dp),           INTENT(INOUT) :: rhs(:)
    INTEGER,            INTENT(IN)    :: bPermu(:), bInvPermu(:), bPtr(:)
    !INTEGER :: nD, nAq

    INTEGER :: i, jj
    REAL(dp) :: b(SIZE(bPermu))
    b(bPermu) = rhs

    !-- L-solve
    DO i = 1 , LU%n
      DO jj = LU%RowPtr(i), LU%DiagPtr(i)-1
        CALL SubMult(b(bPtr(i):bPtr(i+1)-1), & 
                   & LU%Val(LU%ValPtr(jj):LU%ValPtr(jj+1)-1), &
                   & b(bPtr(LU%ColInd(jj)):bPtr(LU%ColInd(jj)+1)-1), &
                   & bPtr(i+1)-bPtr(i), &
                   & LU%ValPtr(jj+1)-LU%ValPtr(jj), &
                   & bPtr(LU%ColInd(jj)+1)-bPtr(LU%ColInd(jj)))
      END DO
    END DO

    !--  U-solve
    DO i = LU%n , 1 , -1
      DO jj = LU%DiagPtr(i)+1 , LU%RowPtr(i+1)-1
        CALL SubMult(b(bPtr(i):bPtr(i+1)-1), & 
                   & LU%Val(LU%ValPtr(jj):LU%ValPtr(jj+1)-1), &
                   & b(bPtr(LU%ColInd(jj)):bPtr(LU%ColInd(jj)+1)-1), &
                   & bPtr(i+1)-bPtr(i), &
                   & LU%ValPtr(jj+1)-LU%ValPtr(jj), &
                   & bPtr(LU%ColInd(jj)+1)-bPtr(LU%ColInd(jj)))
      END DO
      b(bPtr(i):bPtr(i+1)-1) = b(bPtr(i):bPtr(i+1)-1) / LU%Val(LU%ValPtr(LU%DiagPtr(i)):LU%ValPtr(LU%DiagPtr(i)+1)-1)
    END DO

    rhs( bInvPermu ) = b
  
  END SUBROUTINE SolveSparse_ValPtr
  !
  !
  PURE SUBROUTINE SubMult(a, b, c, n_a, n_b, n_c)
    REAL(dp), DIMENSION(n_a), INTENT(INOUT) :: a
    REAL(dp), DIMENSION(n_b), INTENT(IN) :: b
    REAL(dp), DIMENSION(n_c), INTENT(IN) :: c
    INTEGER, INTENT(IN) :: n_a, n_b, n_c

    INTEGER :: i
    
    ! possible cases: n_a  n_b  n_c
    !         case 1:   1    1    1
    !         case 2:  nD   nD   nD
    !         case 3:  nD   nD    1
    !         case 4:   1   nD   nD
    !
    ! NOT possible (will not happen): nD 1 nD 

    IF      (n_a==n_c) THEN ! case 1 or case 2
      a = a - b*c
    ELSE IF (n_a==n_b) THEN ! case 3
      a = a - b*c(1)
    ELSE                    ! case 4
      DO i=1,n_b
        a = a - b(i)*c(i)
      END DO
    END IF
  END SUBROUTINE SubMult

  SUBROUTINE CSR_2_Empty_ValPtr( CSR, vec_entry, n )
    ! this routine changes a csr matrix, maintaining the sparsity pattern 
    ! but activating vector entries

    ! CSR matrix to become ValPtr (vector-entried)
    TYPE(CSR_Matrix_T), INTENT(INOUT) :: CSR
    ! vector with logical values determining if a row/col should become vector-entried
    LOGICAL, DIMENSION(:), INTENT(IN) :: vec_entry
    ! number of entries for a vector entry
    INTEGER :: n

    INTEGER :: i, jj, j


    ALLOCATE( CSR%ValPtr( CSR%RowPtr(SIZE(CSR%RowPtr)) ) )
    CSR%ValPtr = 0
    CSR%ValPtr(1) = 1

    DO i = 1 , CSR%m
      DO jj = CSR%RowPtr(i), CSR%RowPtr(i+1)-1
        j = CSR%ColInd(jj)
        IF ( vec_entry(i) .OR. vec_entry(j) ) THEN
          CSR%ValPtr(jj+1) = CSR%ValPtr(jj) + n
        ELSE
          CSR%ValPtr(jj+1) = CSR%ValPtr(jj) + 1
        END IF
      END DO
    END DO

    DEALLOCATE( CSR%Val )
    ALLOCATE( CSR%Val( CSR%ValPtr(SIZE(CSR%ValPtr))-1 ) )

    CSR%vector_entried = .TRUE.

  END SUBROUTINE CSR_2_Empty_ValPtr

END MODULE Sparse_Mod
