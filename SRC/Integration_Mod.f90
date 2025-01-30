!
! Copyright (C) 2025 Levin Rug, Willi Schimmel (E-Mail: l.rug@lmu.de)
! See ./SRC/Cminor.f90 for the copyright notice
! See ./LICENSE for license information
! SPDX-License-Identifier: GPL-3.0
!
!=========================================================================!
!                                                                         !         
!                   Module for time integration with diagonal             !
!                      implicite Rosenbrock-Wanner-Methods                !         
!                                                                         !         
!=========================================================================!
!
MODULE Integration_Mod
  !
  USE Kind_Mod,         ONLY: dp
  USE Control_Mod,      ONLY: ZERO, Out, combustion, eps, minStp, maxStp, NetcdfPrint,     &
                            & ONE, rTEN, StpNetCdf, tBegin, tEnd, combustion, TimeNetCdf,  &
                            & TimeIntegration, TimeNetCDFA, TWO, WaitBar, pressure,        &
                            & q_parcel, T_parcel, rho_parcel, updraft_velocity, z_parcel,  &
                            & adiabatic_parcel, nDropletClasses, RH, Pi43, m_parcel,       &
                            & milli, rTHREE, Pi34, kilo, Start_Timer, End_Timer
  !
  USE Reac_Mod,         ONLY: nDIM, nDIM2, Diag_Index, nspc, nspc2, ns_gas, iAqMassEq2,    &
                            & iRhoEq2, iZeq2, ns_AQUA, DropletClasses, iTeq2, iqEq2
  !
  USE ChemKinInput_Mod, ONLY: Density
  USE Sparse_Mod,       ONLY: Jac_CC, BAT, A
  USE Rosenbrock_Mod,   ONLY: Rosenbrock, Ros
  USE NetCDF_Mod,       ONLY: NetCDF, SetOutputNcdf, StepNetCDF
  USE Meteo_Mod,        ONLY: cv, pressure_from_height, qsatw, rho_H2O, get_wet_radii,     &
                            & LWC_array
  !
  USE Rates_Mod,        ONLY: rhs_condensation, rhs_T_cond_and_parcel, rhs_q_parcel,       &
                            & rhs_rho, UpdateTempArray

  IMPLICIT NONE
  !
  INTEGER :: iBar=-1              ! waitbar increment
  !
  CONTAINS
  !
  !======================================================================= 
  !===================    Time Integration Routine  ======================
  !======================================================================= 
  SUBROUTINE Integrate(y_iconc, Temperature, h0, Tspan)
    !--------------------------------------------------------------------
    ! Input:
    !   - y0 ............. Initial vector
    !   - R0 ............. reaction rates at time=Tspan(1)
    !   - Tspan .......... (/ SimulationTimeStart , SimulationTimeEnd /)
    !   - RtolRow ........ rel. tolerance for Rosenbrock-Wanner-Method
    !   - method.......... Rosenbrock-Wanner-Method
    !   - PrintSpc ....... print spc PrintSpc(1:3)
    REAL(dp), INTENT(INOUT) :: y_iconc(nspc2)
    REAL(dp), INTENT(INOUT) :: h0, Temperature
    REAL(dp) :: Tspan(2)

    ! Temporary variables:
    !
    REAL(dp) :: Y0(nDIM2), Y(nDIM2)       ! old, current y vector
    !
    REAL(dp) :: t             ! current time
    REAL(dp) :: time_int = 0.0d0

    REAL(dp) :: h, tnew, tmp
    REAL(dp) :: error
    REAL(dp) :: tmp_tb
    INTEGER  :: ierr(1,1)
    ! 

    LOGICAL :: done=.FALSE.
    LOGICAL :: failed
    !

    done = .FALSE.

    Y0(1:nspc2) = y_iconc(:)
    Y(1:nspc2)  = y_iconc(:)

    IF ( adiabatic_parcel ) THEN
      Y0(iAqMassEq2) = DropletClasses%waterMass
      Y0(iTeq2)      = T_parcel
      Y0(iqEq2)      = q_parcel
      Y0(iRhoEq2)    = rho_parcel
      Y0(iZeq2)      = z_parcel
    END IF

    IF ( combustion ) THEN
      !--- initial temperature
      Y0(nDIM) = Temperature
      Y(nDIM)  = Temperature

      Y0 = MAX( ABS(Y0), 1.0e-100_dp ) * SIGN( ONE , Y0 )    ! |y| >= eps
      Y  = MAX( ABS(Y) , 1.0e-100_dp ) * SIGN( ONE , Y  )    ! |y| >= eps
    END IF

    t = Tspan(1)
    tnew = Tspan(1)
    h = h0

    ! this is for the waitbar
    tmp_tb = (tEnd-tBegin) * 0.01_dp

    CALL Start_Timer(time_int)
            
    MAIN_LOOP_ROSENBROCK_METHOD: DO 
          
      h = MIN( maxStp, MAX( minStp , h ) )

        !-- Stretch the step if within 5% of tfinal-t.
        IF ( 1.05_dp * h >= Tspan(2) - t ) THEN
          h = ABS(Tspan(2) - t)
          done = .TRUE.
        END IF

        DO    ! Evaluate the formula.
          !-- LOOP FOR ADVANCING ONE STEP.
          failed = .FALSE.      ! no failed attempts

          ! Rosenbrock Timestep 
          CALL Rosenbrock(  Y             & ! new concentration
          &               , error         & ! error value
          &               , ierr          & ! max error component
          &               , Y0            & ! current concentration 
          &               , t             & ! current time
          &               , h             ) ! stepsize

          tnew  = t + h
          IF (done) THEN
            tnew  = Tspan(2)    ! Hit end point exactly.
            h     = tnew-t      ! Purify h.
          END IF
          Out%ndecomps   = Out%ndecomps   + 1
          Out%nRateEvals = Out%nRateEvals + ROS%nStage
          Out%nSolves    = Out%nSolves    + ROS%nStage

          failed = (error > ONE)

          IF (failed) THEN      ! failed step
            ! Accept the solution only if the weighted error is no more than the
            ! tolerance one.  Estimate an h that will yield an error of rtol on
            ! the next step or the next try at taking this step, as the case may be,
            ! and use 0.8 of this value to avoid failures.
            Out%nfailed  = Out%nfailed+1
            IF ( h <= minStp ) THEN
              STOP '....Integration_Mod '
            END IF

            h    = MAX( minstp , h * MAX( rTEN, 0.8_dp * error**(-ROS%pow) ) )
            done = .FALSE.
          ELSE ! successful step
            EXIT
          END IF
      END DO
      Out%nsteps = Out%nsteps + 1

      IF ( adiabatic_parcel ) THEN
        rho_parcel               = Y(iRhoEq2)
        T_parcel                 = Y(iTeq2)
        DropletClasses%waterMass = Y(iAqMassEq2)
        DropletClasses%wetRadius = get_wet_radii()
        q_parcel                 = Y(iqEq2)
        z_parcel                 = Y(izEq2)
        pressure                 = pressure_from_height(z_parcel) ! pressure of parcel instantaneously adjusts to ambient conditions
        RH                       = q_parcel / qsatw(T_parcel, pressure)
      ELSE IF ( ns_AQUA>0 ) THEN
        DropletClasses%waterMass = kilo * LWC_array(t) / rho_parcel
        DropletClasses%wetRadius = get_wet_radii()
      END IF

      ! --- save to NetCDF file
      IF ( NetCdfPrint ) THEN
        IF ( done .OR. StpNetCDF < ZERO ) THEN 
          CALL Start_Timer(TimeNetCDFA)
          IF (combustion) Temperature = Y(nDIM)
          CALL SetOutputNCDF( NetCDF, tnew , h , Y , Temperature )
          CALL StepNetCDF( NetCDF )
          CALL End_Timer(TimeNetCDF, TimeNetCDFA)
        END IF
      END IF 

      !-- Call progress bar.
      IF ( WaitBar .AND. t-tBegin > iBar*tmp_tb ) CALL Progress(iBar)

      !-- If there were no failures compute a new h.
      tmp = 1.25_dp * error**ROS%pow
      IF ( TWO * tmp > ONE ) THEN
        h = h / tmp
      ELSE
        h = h * TWO
      END IF
          
      !-- Advance the integration one step.
      t  = tnew
      Y0 = Y

      !-- Termination condition for the main loop.
      IF ( done ) EXIT  MAIN_LOOP_ROSENBROCK_METHOD

    END DO MAIN_LOOP_ROSENBROCK_METHOD  ! MAIN LOOP
    
    ! save values for next initial vector
    IF ( combustion ) Temperature = Y(nDIM)
    y_iconc = Y(1:nspc2)
    h0 = h

    CALL End_Timer(TimeIntegration, time_int)

  END SUBROUTINE Integrate

  ! The Progressbar
  SUBROUTINE Progress(j)
    !
    INTEGER(4)    :: j,k
    CHARACTER(69) :: bar="          Start Integration...........    ???% |                    |"

    iBar = iBar + 1 ! iBar runs from 0 to 100
    !
    WRITE(bar(43:45),'(I3)') j
    !
    DO k=1,j/5
      bar(48+k:48+k)="*"
    END DO
    ! print the progress bar.
    WRITE(*,'(A1,A69,$)') char(13), bar
  END SUBROUTINE Progress
  !

END MODULE Integration_Mod
