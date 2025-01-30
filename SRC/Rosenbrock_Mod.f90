!
! Copyright (C) 2025 Levin Rug, Willi Schimmel (E-Mail: l.rug@lmu.de)
! See ./SRC/Cminor.f90 for the copyright notice
! See ./LICENSE for license information
! SPDX-License-Identifier: GPL-3.0
!
!=========================================================================!
!                                                                         !         
!                                                                         !
!          Rosenbrock Modul for integrating one time step                 ! 
!                                                                         ! 
!                                                                         !
!=========================================================================!

MODULE Rosenbrock_Mod
  !
  USE Kind_Mod,         ONLY: dp
  USE Sparse_Mod,       ONLY: SolveSparse, Miter, LU_Miter, Jac_CC, BAT, MATMUL_nD, OPERATOR(*), Miter_Classic, &
                            & Jacobian_CC_ValPtr, Jacobian_CT, Jacobian_TC, Jacobian_TT, Miter_Classic_ValPtr,  &
                            & SetLUvaluesCL, SolveSparse_ValPtr, SparseLU_ValPtr, CSR_Matrix_T, Copy_CSR,       &
                            & A, MULT_BAT_Rate_ValPtr, SparseLU, Jacobian_CC, Jacobian_parcel, SparseAdd, SparseID
  !
  USE Chemsys_Mod,      ONLY: nsr
  !
  USE Meteo_Mod,        ONLY: pressure_from_height, qsatw
  !
  USE Rates_Mod,        ONLY: ReactionRates, UpdateTempArray, InternalEnergy, DiffInternalEnergy,               &
                            & Diff2InternalEnergy, UpdateEmission, rhs_T_cond_and_parcel, rhs_rho, rhs_z,       &
                            & rhs_condensation, rhs_q_parcel
  !
  USE Control_Mod,      ONLY: Out, Error_Est, ONE, aTolAll, ConcDataPrint, eps, Start_Timer, End_Timer, Tspan,  &
                            & FluxDataPrint, iStpConc, iStpFlux, mONE, Tspan_tot, TimeSetValues, m_parcel,      &
                            & RTolRow, StpConc, StpFlux, combustion, TimeConcWrite, TimeErrCalc, TimeFac,       &
                            & TimeJac, TimeFluxWrite, TimeJacobianA, TimeRhsCalc, TimeSolve, ZERO, rho_parcel0, &
                            & AtolAqua, AtolGas, AtolTemp, maxStp, minStp, rTWO, maxErrorCounter, pressure, RH, &
                            & nDropletClasses, adiabatic_parcel
  !
  USE Reac_Mod,         ONLY: y_name, nspc, SCperm, nDIM, nDIM2, iGs, iGs2, iAs2, hasGasSpc, hasAquaSpc, iqEq,  &
                            & nspc2, rRho, rho, nreac2, nDIM, nreac, rNspc, iAqMassEq2, iTeq2, iqEq2, iRhoEq2,  &
                            & iZeq2, DropletClasses, MolMass, iAqMassEq
  !
  USE ChemKinInput_Mod, ONLY: MassAveMixSpecHeat

  IMPLICIT NONE
  !
  !
  ! Rosenbrock-Parameter
  TYPE RosenbrockMethod_T
    INTEGER :: Order                                ! Classical approximation order of the method
    INTEGER :: nStage                               ! Number of stages
    INTEGER :: pinterp                              ! Interpolation order
    REAL(dp) :: ga                            ! Diagonalentry gamma
    REAL(dp) :: pow                           ! needed for sitepsize control pow=1/nstage
    REAL(dp), ALLOCATABLE :: Asum(:)          ! Row sum of A
    REAL(dp), ALLOCATABLE :: Alpha(:,:)       ! Propagation table, strictly lower triangular
    REAL(dp), ALLOCATABLE :: a(:,:)           ! Propagation table, strictly lower triangular (converted Alpha)
    REAL(dp), ALLOCATABLE :: Gamma(:,:)       ! Stage table, lower triangular with nonzero diagonal
    REAL(dp), ALLOCATABLE :: iGamma(:,:)      ! inverse Stage table
    REAL(dp), ALLOCATABLE :: C(:,:)           ! Stage table, lower triangular with nonzero diagonal (converted Gamma)
    REAL(dp), ALLOCATABLE :: B(:)             ! Step completion table
    REAL(dp), ALLOCATABLE :: m(:)             ! Step completion table(converted B)
    REAL(dp), ALLOCATABLE :: Be(:)            ! Step completion table for embedded method of order one less
    REAL(dp), ALLOCATABLE :: me(:)            ! Step completion table for embedded method of order one less (converted Be)
  END TYPE RosenbrockMethod_T


  TYPE(RosenbrockMethod_T) :: ROS

  REAL(dp), PRIVATE :: timerStart
  REAL(dp), ALLOCATABLE :: rRate0(:)

  INTEGER, ALLOCATABLE :: LU_Perm(:), LU_Perm_ValPtr(:), w_InvColInd(:,:,:), &
                        & bPermu(:), bInvPermu(:), bPtr(:)

  CONTAINS
  !

  SUBROUTINE Rosenbrock(YNew, err, ierr, Y0, t, h)
    !--------------------------------------------------------
    ! Input:
    !   - Y0............. concentrations at Time = t
    !   - t.............. time
    !   - h.............. step size
    !
    REAL(dp),          INTENT(IN) :: Y0(nDIM2)
    REAL(dp),          INTENT(IN) :: t, h
    !--------------------------------------------------------
    ! Output:
    !   - Ynew........... new concentratinos 
    !   - err............ error calc with embedded formula.
    !
    REAL(dp), INTENT(OUT) :: YNew(nDIM2)
    REAL(dp), INTENT(OUT) :: err
    INTEGER , INTENT(OUT) :: ierr(1,1)
    !-------------------------------------------------------
    ! TemporarY variables:
    !
    REAL(dp), DIMENSION(nDIM2)   :: Yhat
    REAL(dp) :: k( nDIM2 , ROS%nStage )
    REAL(dp), ALLOCATABLE :: rhs(:)

    !
    INTEGER :: iStg

    ALLOCATE(rhs(nDIM2))

    ! Initial settings
    rhs = ZERO
    k   = ZERO

    LOOP_n_STAGES:  DO iStg = 1 , ROS%nStage

      IF ( iStg==1 ) THEN ! calculate Miter matrix, rhs and LU decomposition

        CALL Assemble_Miter_and_First_Rhs_Classic(rhs, Y0, t, h)
        ! --- LU - Decomposition ---
        CALL Start_Timer(timerStart)
        IF (nDropletClasses>1) THEN
          CALL SparseLU_ValPtr( LU_Miter, w_InvColInd )
        ELSE
          CALL SparseLU( LU_Miter )
        END IF
        CALL End_Timer(TimeFac, timerStart)

      ELSE ! IF (iStg>1): no Miter re-calculation, only rhs update

        CALL Assemble_Interstitial_Rhs(rhs, Y0, t, h, k, iStg)

      END IF

      CALL Start_Timer(timerStart)
      IF (nDropletClasses>1) THEN
        CALL SolveSparse_ValPtr( LU_Miter, rhs, bPermu, bInvPermu, bPtr)
      ELSE
        CALL SolveSparse( LU_Miter , rhs )
      END IF
      CALL End_Timer(TimeSolve, timerStart)

      ! write to k(:,iStg)
      CALL Assign_Interstitial_State(k, rhs, iStg)

    END DO  LOOP_n_STAGES

    !--- Update Concentrations (+Temperatur)
    YNew = Y0
    DO iStg=1,ROS%nStage; YNew = YNew + ROS%m(iStg)*k(:,iStg); END DO

    CALL Start_Timer(timerStart)
    ! embedded formula for err calc ord-1
    YHat = Y0
    DO iStg=1,ROS%nStage; YHat = YHat + ROS%me(iStg)*k(:,iStg); END DO

    CALL ERROR( err , ierr , YNew , YHat , ATolAll , RTolROW )
    CALL End_Timer(TimeErrCalc, timerStart)

    maxErrorCounter(ierr(1,1)) = maxErrorCounter(ierr(1,1)) + 1

  END SUBROUTINE Rosenbrock

  SUBROUTINE Assemble_Interstitial_Rhs(rhs, Y0, t, h, k, iStg)
    ! IN:
    REAL(dp) :: t, h
    REAL(dp) :: Y0(ndim2)
    REAL(dp) :: k( nDIM2 , ROS%nStage )
    INTEGER :: iStg
    ! OUT:
    REAL(dp) :: rhs(ndim2)

    REAL(dp) :: tt, Y(ndim2), Rate(nreac2), Emiss(nspc2)
    INTEGER :: jStg

    tt  = t + ROS%Asum(iStg) * h
    Y   = Y0

    DO jStg = 1 , iStg-1
      Y = Y + ROS%a(iStg,jStg) * k(:,jStg)
    END DO

    ! Update Rates at  (t + SumA*h) , and  (Y + A*)k
    IF (combustion) THEN
      CALL ReactionRates( Y , Rate )
    ELSE
      CALL ReactionRates( tt , Y , Rate )
    END IF

    CALL UpdateEmission(Emiss, Y)

    CALL Assemble_Rhs_Classic(rhs, Rate, Emiss, Y, h, k=k, iStg=iStg)

  END SUBROUTINE Assemble_Interstitial_Rhs

  SUBROUTINE Assemble_Miter_and_First_Rhs_Classic(rhs, Y0, t, h)
    ! IN:
    REAL(dp), INTENT(IN) :: Y0(nDIM2)
    REAL(dp) :: t, h
    ! OUT:
    REAL(dp) :: rhs(ndim2)

    REAL(dp), DIMENSION(nspc2)   :: U, dUdT, dCdt, Emiss, Jac_CT, Jac_TC
    REAL(dp), DIMENSION(nDropletClasses+4) :: Jac_parcel
    REAL(dp), DIMENSION(nDropletClasses) :: dSeqdmw
    REAL(dp) :: dTdt, cv, dcvdT, Y(ndim2), Jac_TT
    REAL(dp) :: Rate(nreac2), dkdT_over_k(nreac2)


    Rate    = ZERO
    Y       = Y0
    Emiss   = ZERO

    IF ( combustion ) THEN
      CALL ReactionRates( Y0, Rate, dkdT_over_k )
      dkdT_over_k = ROS%ga*dkdT_over_k
    ELSE
      Y   = MAX( ABS(Y0)  , eps ) * SIGN( ONE , Y0 )  ! concentrations =/= 0
      CALL ReactionRates( t, Y, Rate )
    END IF

    CALL UpdateEmission(Emiss,Y)

    ! calculate right hand side
    IF ( combustion ) THEN
      CALL Assemble_Rhs_Classic(rhs, Rate, Emiss, Y, h, dCdt_out=dCdt, dTdt_out=dTdt)
    ELSE IF ( adiabatic_parcel ) THEN
      CALL Assemble_Rhs_Classic(rhs, Rate, Emiss, Y, h, dSeqdmw_out=dSeqdmw)
    ELSE 
      CALL Assemble_Rhs_Classic(rhs, Rate, Emiss, Y, h)
    END IF

    ! --- Update matrix procedure
    CALL Start_Timer(TimeJacobianA)

    ! d(dcdt)/dc
    IF (nDropletClasses>1) THEN
      CALL Jacobian_CC_ValPtr( Jac_CC , BAT  , A , Rate , Y )
    ELSE
      CALL Jacobian_CC( Jac_CC , BAT  , A , Rate , Y )
    END IF

    IF ( combustion ) THEN
    !
      CALL Calculate_Heat_Sources(U, dUdT, cv, Y0(nDIM2), Y0(1:nspc2), dcvdT=dcvdT)
      CALL Jacobian_CT( Jac_CT , BAT , Rate , dkdT_over_k )
      CALL Jacobian_TC( Jac_TC , Jac_CC , cv , dUdT , dTdt , U , rRho)
      CALL Jacobian_TT( Jac_TT , Jac_CT , cv , dcvdT , dTdt , dUdT , dCdt , U , rRho)
      !
      CALL Miter_Classic( Miter , h , ROS%ga , Jac_CC , Jac_TC , Jac_CT , Jac_TT )
    !
    ELSE IF ( adiabatic_parcel ) THEN
    !
      CALL Jacobian_parcel(Jac_parcel, dSeqdmw)
      IF (nDropletClasses>1) THEN
        CALL Miter_Classic_ValPtr( Miter , h , ROS%ga , Jac_CC, J2=Jac_parcel )
      ELSE
        CALL Miter_Classic( Miter , h , ROS%ga , Jac_CC, J2=Jac_parcel )
      END IF
    !
    ELSE ! only chemistry box-model
    !
      IF (nDropletClasses>1) THEN
        CALL Miter_Classic_ValPtr( Miter , h , ROS%ga , Jac_CC )
      ELSE
        CALL Miter_Classic( Miter , h , ROS%ga , Jac_CC )
      END IF
    !
    END IF
    Out%npds = Out%npds + 1

    CALL End_Timer(TimeJac, TimeJacobianA)

    ! finish by permuting for better LU decomposition
    CALL Start_Timer(timerStart)
    CALL SetLUvaluesCL( LU_Miter , Miter , LU_Perm_ValPtr )
    CALL End_Timer(TimeSetValues, timerStart)

    ! write data to files if needed
    CALL StreamWriteData(t, h, Rate, Y0)

  END SUBROUTINE Assemble_Miter_and_First_Rhs_Classic

  SUBROUTINE Assemble_Rhs_Classic(rhs, Rate, Emiss, Y, h, k, iStg, dCdt_out, dTdt_out, dSeqdmw_out)
    ! IN:
    REAL(dp) :: Rate(nreac2), Emiss(nspc2), Y(ndim2), h
    REAL(dp), OPTIONAL :: k( nDIM2 , ROS%nStage )
    INTEGER, OPTIONAL :: iStg
    ! OUT:
    REAL(dp) :: rhs(ndim2)
    REAL(dp), OPTIONAL :: dCdt_out(nspc2), dTdt_out, dSeqdmw_out(nDropletClasses)

    REAL(dp) :: cv, dCdt(nspc2), dTdt, dmdt(nDropletClasses), drhodt, dqdt, dzdt
    REAL(dp), DIMENSION(nspc2) :: U, dUdT
    INTEGER :: jStg

    CALL Start_Timer(timerStart)
    IF (nDropletClasses>1) THEN
      dCdt = MULT_BAT_Rate_ValPtr(BAT , Rate) + Emiss
    ELSE
      dCdt = BAT * Rate + Emiss
    END IF

    rhs(1:nspc2) = dCdt
    IF (combustion) THEN
      CALL Calculate_Heat_Sources(U, dUdT, cv, Y(nDIM2), Y(1:nspc2))
      dTdt = - SUM(U*dCdt) * rRho / cv
      rhs( nDIM2 ) = dTdt
    END IF

    IF (adiabatic_parcel) THEN

      pressure = pressure_from_height(Y(izEq2))           ! pressure of parcel instantaneously adjusts to ambient conditions
      RH       = Y(iqEq2) / qsatw(Y(iTeq2), pressure)
      dzdt     = rhs_z()
      IF ( PRESENT(dSeqdmw_out) ) THEN ! calculate and output dSeqdm (for Jacobian in first step of every time step)
        CALL rhs_condensation(dmdt, DropletClasses, Y, RH, dSeqdmw_out=dSeqdmw_out)
      ELSE
        CALL rhs_condensation(dmdt, DropletClasses, Y, RH)
      END IF
      dTdt     = rhs_T_cond_and_parcel(dmdt, dzdt)
      drhodt   = rhs_rho(Y(iRhoEq2), Y(iTeq2), dzdt, pressure, dTdt)
      dqdt     = rhs_q_parcel(dmdt)

      rhs(1:nspc2)    = rhs(1:nspc2) + Y(1:nspc2)*drhodt/Y(iRhoEq2)
      rhs(iAqMassEq2) = dmdt
      rhs(iTeq2)      = dTdt
      rhs(iRhoEq2)    = drhodt
      rhs(iqEq2)      = dqdt
      rhs(iZeq2)      = dzdt

    END IF

    rhs = h * rhs
    IF ( PRESENT(iStg) ) THEN
      DO jStg=1,iStg-1
        rhs = rhs + ROS%C(iStg,jStg)*k(:,jStg)
      END DO
    END IF
    IF (PRESENT(dCdt_out)) dCdt_out = dCdt
    IF (PRESENT(dTdt_out)) dTdt_out = dTdt

    CALL End_Timer(TimeRhsCalc, timerStart)
  END SUBROUTINE Assemble_Rhs_Classic

  SUBROUTINE Calculate_Heat_Sources(U, dUdT, cv, T, MoleConc, dcvdT)
    ! IN:
    REAL(dp) :: T, cv, MoleConc(nspc2)
    REAL(dp), OPTIONAL :: dcvdT
    ! OUT:
    REAL(dp), DIMENSION(nspc2) :: U, dUdT

    REAL(dp) :: Tarr(10), d2UdT2(nspc2)

    Tarr = UpdateTempArray ( T )       
    CALL InternalEnergy    ( U       , Tarr)  
    CALL DiffInternalEnergy( dUdT    , Tarr)              
    CALL MassAveMixSpecHeat( cv      , dUdT    , MoleConc=MoleConc)
    IF ( PRESENT(dcvdT) ) THEN
      CALL Diff2InternalEnergy ( d2UdT2  , Tarr)
      CALL MassAveMixSpecHeat  ( dcvdT   , d2UdT2  , MoleConc=MoleConc)
    END IF
  END SUBROUTINE Calculate_Heat_Sources


  SUBROUTINE ERROR(err,ierr,ynew,yhat,ATol,RTol)
    !
    REAL(dp) :: err
    REAL(dp), DIMENSION(nDIM2) :: ynew, yhat, ATol
    REAL(dp) :: RTol
    !
    REAL(dp) :: scalTol(nDIM2), e_n(nDIM2), ymax(nDIM2)
    INTEGER :: ierr(1,1)
    !
    ymax      = MAX(ABS(yhat),ABS(ynew))
    scalTol   = ONE / ( ATol + ymax*RTol )  ! scaling strategy
    e_n       = ABS( ynew - yhat ) * scalTol      ! local error est.
    ierr(1,1) = MAXLOC( e_n , 1 )           ! max error component
    !
    IF ( Error_Est == 2 ) THEN
      !err = SUM( e_n*e_n ) 
      !err = SQRT(err) / rNspc

      err = SUM( e_n*e_n ) * rNspc   ! euclid norm
      !err = SQRT(SUM( e_n*e_n ) * rNspc)  ! euclid norm
    ELSE
      err = MAXVAL( e_n )     ! maximum norm
    END IF

    !err = MAX(err,1.0e-10_dp)
    !
  END SUBROUTINE ERROR


  SUBROUTINE Assign_Interstitial_State(k, rhs, iStg)
    REAL(dp) :: k( nDIM2 , ROS%nStage ), rhs(:)

    INTEGER :: iStg

    k( 1:nDIM2 , iStg ) = rhs

  END SUBROUTINE Assign_Interstitial_State

  SUBROUTINE StreamWriteData(t, h, Rate, Conc)
    USE IO_Mod,   ONLY: StreamWriteFluxes, StreamWriteConcentrations

    REAL(dp) :: t, h, Rate(nreac2), Conc(ndim2)

    IF ( FluxDataPrint .AND. t - Tspan_tot(1) + (Tspan_tot(2)-Tspan_tot(1)) >= StpFlux*REAL(iStpFlux,dp) ) THEN
      CALL Start_Timer(timerStart)
      CALL StreamWriteFluxes(Rate,t,h)
      !WRITE(*,*) 'Captured fluxes at time ',t,'. StpFlux = ',StpFlux,' iStpFlux = ',iStpFlux
      CALL End_Timer(TimeFluxWrite, timerStart)
    END IF

    IF ( ConcDataPrint .AND. t - Tspan_tot(1) + (Tspan_tot(2)-Tspan_tot(1)) >= StpConc*REAL(iStpConc,dp) ) THEN
      CALL Start_Timer(timerStart)
      CALL StreamWriteConcentrations(Conc)
      !WRITE(*,*) 'Captured concentrations at time ',t,'. StpConc = ',StpConc,' iStpConc = ',iStpConc
      CALL End_Timer(TimeConcWrite, timerStart)
    END IF

  END SUBROUTINE StreamWriteData

  !
  !=======================================================
  !       chose one of the methods in ~/METHODS/
  !=======================================================
  SUBROUTINE SetRosenbrockMethod(RCo,method)
    !-------------------------------------------------------------
    ! Input: method ... string with Rosenbrock method path
    TYPE(RosenbrockMethod_T), INTENT(OUT) :: RCo
    CHARACTER(*) :: method
    !-------------------------------------------------------------
    ! Output:
    !   - RosenbrockMethod_T ... coefficients of the chosen one 
    !-------------------------------------------------------------
    ! Temporary variables: 
    REAL(dp), ALLOCATABLE :: ID(:,:)
    INTEGER, ALLOCATABLE :: IPIV(:)
    INTEGER :: INFO

    CHARACTER(20) :: tmethod
    !
    INTEGER :: i
    !
    tmethod = ADJUSTL(method(INDEX(method,'/')+1:INDEX(method,'.')-1))
    SELECT CASE (tmethod)
      CASE ('Ros2AMF')       
        INCLUDE 'METHODS/Ros2AMF.ros'
      CASE ('Ros3w')         
        INCLUDE 'METHODS/Ros3w.ros'
      CASE ('Ros3Dw')        
        INCLUDE 'METHODS/Ros3Dw.ros'
      CASE ('Ros3Pw')        
        INCLUDE 'METHODS/Ros3Pw.ros'
      CASE ('Ros34PW1a')     
        INCLUDE 'METHODS/Ros34PW1a.ros'
      CASE ('Ros34PW2')      
        INCLUDE 'METHODS/Ros34PW2.ros'
      CASE ('Ros34PW3')      
        INCLUDE 'METHODS/Ros34PW3.ros'
      CASE ('Rodas3')  
        INCLUDE 'METHODS/Rodas3.ros'
      CASE ('TSRosW2P')      
        INCLUDE 'METHODS/TSRosW2P.ros'
      CASE ('TSRosW2M')      
        INCLUDE 'METHODS/TSRosW2M.ros'
      CASE ('TSRosWRA34PW2') 
        INCLUDE 'METHODS/TSRosWRA34PW2.ros'
      CASE ('TSRosWSandu3')  
        INCLUDE 'METHODS/TSRosWSandu3.ros'
      CASE DEFAULT
        WRITE(*,*) '    Unknown Method:  ',method
        WRITE(*,*) '    Use Rodas3 instead.'
        INCLUDE 'METHODS/Rodas3.ros'
    END SELECT
    !INCLUDE RosenbrockMethod
    !
    ! converting the butcher tableau 
    ! automatic transformation to avoid mat*vec in ROW methode
    ALLOCATE(ID(ROS%nStage,ROS%nStage) , ROS%iGamma(ROS%nStage,ROS%nStage))
    ROS%pow    = ONE / (ROS%Order+ONE)
    ROS%iGamma = ZERO
    ID         = ZERO
    DO i=1,ROS%nStage
      ROS%iGamma(i,i) = ONE
      ID(i,i) = ONE
    END DO

    ALLOCATE(RCo%Asum(RCo%nStage))
    DO i=1,RCo%nStage
      RCo%Asum(i) = SUM(RCo%Alpha(i,:))
    END DO
    
    ! calculate the inverse matrix of gamma
    !
    ! CAUTION ROS%Gamma (IN) =/= ROS%Gamma (OUT) !!!
    !
    ALLOCATE(IPIV(ROS%nStage))
    CALL dgesv(  ROS%nStage,     &        ! # linear eqations
    &            ROS%nStage,     &        ! # RHS (coloums)
    &            ROS%Gamma,      &        ! Matrix A of A*A^(-1)=ID
    &            ROS%nStage,     &        ! leading dimension of A (nStage)
    &            IPIV,           &        ! pivot indices of dimension nStage
    &            ROS%iGamma,     &        ! Matrix ID of A*A^(-1)=ID
    &            ROS%nStage,     &        ! leading dimension of RHS
    &            INFO)                    ! INFO (integer) if INFO=0 succsessfull	
    !
    IF ( INFO/= 0 ) WRITE(*,*) 'Error while calc row-method parameter'
    !       
    ALLOCATE(ROS%a(ROS%nStage,ROS%nStage))
    ROS%a = ZERO
    ROS%a = ROS%ga*MATMUL(ROS%Alpha, ROS%iGamma)
    !  
    ALLOCATE(ROS%C(ROS%nStage,ROS%nStage))
    ROS%C = ZERO
    ROS%C = ID - ROS%ga * ROS%iGamma
    FORALL (i=1:ROS%nStage) ROS%C(i,i)=ZERO
    !  
    ALLOCATE(ROS%m(ROS%nStage))
    ROS%B = MATMUL(ROS%B, ROS%iGamma)
    ROS%m = ROS%ga * ROS%B
    !
    IF (.NOT.ROS%nStage==1) THEN
      ALLOCATE(ROS%me(ROS%nStage))
      ROS%Be = MATMUL(ROS%Be,ROS%iGamma)
      ROS%me = ROS%ga * ROS%Be(:)
    END IF
    !
    DEALLOCATE(ID)
    DEALLOCATE(IPIV)
  END SUBROUTINE SetRosenbrockMethod
  !
  !
  !==========================================================
  ! Calculates an initial stepsize based on 2. deriv.
  ! This function originates from Curtis 1980, 
  ! "The FACSIMILE numerical integrator for stiff initial value problems"
  !==========================================================
  FUNCTION InitialStepSize(Jac,Rate,t,Y_in,pow) RESULT(h)
    !------------------------------------------------- 
    ! Input:
    !        - public variables
    !        - Tspan 
    !        - Y0  ( initial vector )
    !        - Jacobian matrix
    TYPE(CSR_Matrix_T), INTENT(IN) :: Jac
    REAL(dp), INTENT(IN) :: t, pow
    REAL(dp), INTENT(IN) :: Y_in(:)
    REAL(dP), INTENT(INOUT) :: Rate(:)
    !-------------------------------------------------
    ! Output:
    !        - initial step size
    REAL(dp) :: h 
    !-------------------------------------------------
    !
    ! Temp vars:
    REAL(dp) :: tdel, rh, absh
    REAL(dp), DIMENSION(nDIM2) :: Y
    REAL(dp), DIMENSION(nspc2) :: wt, DfDt, Tmp, f0, f1, zeros, Emiss
    REAL(dp), DIMENSION(nreac2)   :: dkdT_over_k
    REAL(dp), ALLOCATABLE     :: thresh(:)

    zeros = ZERO
    Y = Y_in 

    !---- Compute an initial step size h using Yp=Y'(t) 
    CALL UpdateEmission(Emiss,Y)
    !f0 = BAT * Rate + Emiss
    f0 = MULT_BAT_Rate_ValPtr(BAT , Rate) + Emiss

    ALLOCATE( thresh(nDIM2) )
    IF ( combustion ) THEN
      thresh(iGs)  = AtolGas / RTolROW
      thresh(nDIM) = AtolTemp / RTolROW
    ELSE
      IF ( hasGasSpc  ) thresh(iGs2) = AtolGas / RTolROW
      IF ( hasAquaSpc ) thresh(iAs2) = AtolAqua / RTolROW
    END IF

    wt   = MAX( ABS(Y(1:nspc2)) , thresh(1:nspc2) )
    rh   = ( 1.25_dp * MAXVAL(ABS(f0/wt)) )/(RTolRow**pow)
    absh = MIN( maxStp , Tspan(2)-Tspan(1) )
    IF ( absh * rh > ONE )  absh = ONE / rh

    !---- Compute Y''(t) and a better initial step size
    h = absh
    tdel  = t + MIN( SQRT(eps) * MAX( ABS(t) , ABS(t+h) ) , absh )

    IF (combustion) THEN
      CALL ReactionRates( Y , Rate , dkdT_over_k )
    ELSE
      CALL ReactionRates( t+tdel , Y , Rate)
    END IF
    Out%nRateEvals = Out%nRateEvals + 1

    CALL UpdateEmission(Emiss,Y)

    !f1   = BAT * Rate + Emiss
    f1   = MULT_BAT_Rate_ValPtr(BAT , Rate) + Emiss

    ! DfDt : total derivative of f w.r.t. t
    ! = df/dt (following line) + df/dc * dc/dt 
    !                            = Jac * dc/dt = Jac * f
    DfDt = ( f1 - f0 ) / tdel

    !tmp  = Jac * f0 
    tmp  = MATMUL_nD(Jac, f0) 
    DfDt = DfDt + Tmp

    rh   = 1.25_dp * SQRT( rTWO * MAXVAL( ABS(DfDt/wt) ) ) / RTolRow**pow
  
    absh = MIN( maxStp , Tspan(2)-Tspan(1) )
    IF ( absh*rh > ONE )  absh = ONE/rh
    h = MAX( absh , minStp )

  END FUNCTION InitialStepSize

END MODULE Rosenbrock_Mod
