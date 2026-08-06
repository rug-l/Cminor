!
!  ---------------------------------------------------------------------------
! |                                                                           |
! |  The Chemical Mechanism Integrator (Cminor) - a Fortran software package  | 
! |  for flexible and detailed studies of chemical kinetic systems.           |
! |                                                                           |
! |  Copyright (C) 2025 Levin Rug, Willi Schimmel (E-Mail: l.rug@lmu.de)      |
! |                                                                           |
! |  This program is free software: you can redistribute it and/or modify     |
! |  it under the terms of the GNU General Public License as published by     |
! |  the Free Software Foundation, either version 3 of the License, or        |
! |  (at your option) any later version.                                      |
! |                                                                           |
! |  This program is distributed in the hope that it will be useful,          |
! |  but WITHOUT ANY WARRANTY; without even the implied warranty of           |
! |  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the             |
! |  GNU General Public License for more details.                             |
! |                                                                           |
! |  You should have received a copy of the GNU General Public License        |
! |  along with this program (File ''Cminor/LICENSE'').                       |
! |  If not, see https://www.gnu.org/licenses/gpl-3.0.txt.                    |
! |  SPDX-License-Identifier: GPL-3.0                                         |
! |                                                                           |
!  ---------------------------------------------------------------------------
!
!
! Cminor_Driver_Mod — box-model orchestration (init / run / reduce / finalize).
! Option 1: Integrate() still owns the Rosenbrock loop + NetCDF writes.
! Later (option 2): extract one Rosenbrock step and move the time loop here.
!
MODULE Cminor_Driver_Mod

  USE Kind_Mod,         ONLY: dp

  USE InitRoutines_Mod, ONLY: InitRun

  USE Rates_Mod,        ONLY: ReactionRates, RateCnt

  USE Sparse_Mod,       ONLY: SymbLU_SpRowColD_M, CSR_Matrix_T, SpRowColD_T, A, B, SparseAdd,   &
                            & BA, BAT, SparseID, RowColD_to_CSR, Miter, LU_Miter, Jac_CC,       &
                            & CSR_to_SpRowColD, BuildSymbolicClassicMatrix,                    &
                            & CSR_2_Empty_ValPtr, Free_Matrix_CSR, Free_SpRowColD, SymbolicAdd, &
                            & Get_LU_Permutation, Get_ValPtr_Permutations, Jacobian_CC_ValPtr,  &
                            & SymbLU_SpRowColD, SymbolicMult, TransposeSparse,                  &
                            & WriteSparseMatrix

  USE Chemsys_Mod,      ONLY: Print_ChemFile, ReactionSystem, nsr, nsr2,                         &
                            & InputChemicalData, Read_Diag, make_ChemSys_1_to_nD_arrays,        &
                            & Read_DiagEmiss, Read_INI_file, ReadSystem, Setup_SpeciesOrder

  USE Integration_Mod,  ONLY: Integrate

  USE Rosenbrock_Mod,   ONLY: SetRosenbrockMethod, ROS, LU_Perm, LU_Perm_ValPtr, w_InvColInd,   &
                            & bPermu, bInvPermu, bPtr, InitialStepSize

  USE Control_Mod,      ONLY: nD_spc, Tspan, ZERO, clock_rate, clock_maxcount, SysUnit,         &
                            & TimeSymbolic, TimeNetCDF, Timer_Finish, Start_Timer, End_Timer,   &
                            & TimeJac, combustion, tBegin, tEnd, T_parcel, q_parcel,            &
                            & Simulation, StpNetCDF, Temperature0, SysFile, RunFile, AtolGas,   &
                            & Pressure0, Pa_to_dyncm2, Ordering, ONE, Out, adiabatic_parcel,    &
                            & ODEsolver, NetCDFprint, nDropletClasses, MWFile, MWUnit, eps,     &
                            & maxErrorCounter, MatrixPrint, ATolAll, Atolq, AtolRho, Atolz,     &
                            & iStpFlux, InitFile, FluxUnit, FluxMetaUnit, AtolWaterMass, LABEL, &
                            & FluxFile, FluxMetaFile, FluxDataPrint, z_parcel, rho_parcel,      &
                            & DataUnit, DataFile, ConcUnit, ConcFile, TimerNetCDF, Timer_Read,  &
                            & ConcMetaUnit, ConcMetaFile, ConcDataPrint, ChemUnit, AtolTemp,    &
                            & AtolAqua, ChemFile, Time_Finish, Time_Read, TimerSymbolic, RH,    &
                            & Reduction, ReduceOnly
#ifdef ISSA
  USE ISSA_Reduce_Mod, ONLY: RunISSAReduction
#endif

  USE Reac_Mod,         ONLY: InitValAct, y_name, bGs, y_emi, y_depos, iAs2, iGs2, nDIM2, rho,  &
                            & bHr, rnspc, bAs2, SwitchTemp, rRho, MW, rMW, ddelGFEdT,          &
                            & nspc2, nspc, ns_GAS, nr_special, GFE, dGFEdT, delGFE,            &
                            & nreac, nreac2, iRhoeq2, iAqMassEq2, DropletClasses,              &
                            & nDIM, InitValKAT, iGs, ns_KAT, iTeq2, iqEq2, izeq2, hasAquaSpc,   &
                            & hasGasSpc, iR

  USE IO_Mod,           ONLY: Logo, OpenFile_wseq, OpenFile_wstream, Output_Statistics,         &
                            & Print_Run_Param, ShowMaxErrorCounter, Matrix_Statistics

  USE ChemKinInput_Mod, ONLY: Read_Elements, MoleFr_To_MassFr, MoleFr_To_MoleConc, Density,     &
                            & GetSpeciesNames, Read_MolecularWeights,                          &
                            & Read_Reaction, Read_Species, Read_ThermoData,                     &
                            & Setup_ThirdBodyIndex, Setup_ReactionIndex

  USE NetCDF_Mod,       ONLY: InitNetCDF, SetOutputNCDF, StepNetCDF, NetCDF

  USE fparser,          ONLY: initf, parsef

  USE Meteo_Mod,        ONLY: Set_pseudoLWCbounds, LWCb, LWC_array

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cminor_initialize, cminor_run, cminor_reduce, cminor_finalize

  INTEGER(8) :: StartTimer
  REAL(dp), ALLOCATABLE :: Y(:), Y_initial(:)
  REAL(dp), ALLOCATABLE :: Rate(:)
  REAL(dp), ALLOCATABLE :: MoleFrac(:), MassFrac(:), MoleConc(:)
  REAL(dp) :: Press_in_dyncm2
  TYPE(CSR_Matrix_T) :: Id, tmpJacCC
  TYPE(SpRowColD_T) :: temp_LU_Dec
  INTEGER, ALLOCATABLE :: PivOrder(:)
  CHARACTER(14) :: fmt0 = '  [molec/cm3] '
  CHARACTER(1) :: inpt = ''
  LOGICAL :: done = .FALSE.
  REAL(dp) :: h0
  REAL(dp) :: dummyNumber(1), dummyLWC(1)
  REAL(dp), ALLOCATABLE :: nDropletsNCDF(:), LWCsNCDF(:)

CONTAINS

  !================================================================
  ! Public driver API
  !================================================================

  SUBROUTINE cminor_initialize(run_path)
    CHARACTER(*), INTENT(IN) :: run_path

    CALL SYSTEM_CLOCK(count_rate = clock_rate)
    CALL SYSTEM_CLOCK(count_max  = clock_maxcount)
    CALL Logo
    RunFile = TRIM(ADJUSTL(run_path))

    CALL InitRun
    WRITE(*,777) 'Initialize run-file .......... done'

    LWCb = Set_pseudoLWCbounds()

    CALL Start_Timer(Timer_Read)
    CALL Start_Timer(Timer_Finish)

    !--- mechanism + INI (Mode 1 host reuses) --------------------------
    WRITE(*,777,ADVANCE='NO') 'Reading sys-file .............'
    IF ( combustion ) THEN
      CALL Init_Combustion
    ELSE
      CALL Init_Atmosphere
    END IF

    WRITE(*,777) 'Reading ini-file ............. done'
    WRITE(*,777) 'Printing chem-file ........... done'
    WRITE(*,*)
    IF ( nDropletClasses <= 1 ) THEN
      WRITE(*,'(10X,A,I6)') '    Number of Reactions = ', nreac
      WRITE(*,'(10X,A,I6)') '    Number of Species   = ', nspc
    ELSE
      WRITE(*,'(10X,A,I6,A4,I7,A1)') '    Number of Reactions (with droplet classes) = ', nreac, '   (', nreac2, ')'
      WRITE(*,'(10X,A,I6,A4,I7,A1)') '    Number of Species   (with droplet classes) = ', nspc , '   (', nspc2 , ')'
    END IF

    CALL Read_Diag( NetCDF%spc_Pos , NetCDF%spc_Phase , InitFile )
    CALL Read_DiagEmiss( NetCDF%Emiss_Pos, InitFile )
    CALL End_Timer(Timer_Read, Time_Read)

    !--- species order, ODE dims, Atol, Y0 (Mode 1 host reuses) --------
    CALL Init_State

    !--- box I/O only (skip for Mode 1 host) ---------------------------
    CALL Init_NetCDF_File

    WRITE(*,777,ADVANCE='NO') 'Symbolic-phase................'
    CALL Start_Timer(TimerSymbolic)

    !--- sparse BA/BAT + Jac_CC pattern (Mode 1 host reuses) -----------
    CALL SymbolicAdd( BA , B , A )
    CALL SparseAdd  ( BA , B , A, '-' )
    CALL TransposeSparse( BAT , BA )
    CALL SymbolicMult( BAT , A , tmpJacCC )
    Id = SparseID( nspc )
    CALL SymbolicAdd( Jac_CC , Id , tmpJacCC )
    CALL Free_Matrix_CSR( Id )
    ALLOCATE(maxErrorCounter(nDIM2))
    maxErrorCounter = 0

    !--- box Rosenbrock LU only (Mode 1 keeps host GS/LU) --------------
    CALL SetRosenbrockMethod( ROS , ODEsolver )
    CALL BuildSymbolicClassicMatrix( Miter, Jac_CC )
    CALL Init_Symbolic_LU
    WRITE(*,*) 'done'
    WRITE(*,*)

    CALL Write_Sparse_Matrices
    WRITE(*,777) 'Matrix Statistics: '
    WRITE(*,*)
    CALL Matrix_Statistics(A,B,BAT,tmpJacCC,Miter,LU_Miter)
    CALL Free_Matrix_CSR( tmpJacCC )
    CALL Free_SpRowColD( temp_LU_Dec )
    CALL End_Timer(TimerSymbolic, TimeSymbolic)

    !--- initialize flux and conc files for reduction (Mode 1 host reuses) ---
    IF ( Simulation .AND. .NOT. ReduceOnly ) THEN
      IF ( FluxDataPrint ) THEN
        iStpFlux = 0
        FluxFile     = 'OUTPUT/flux_'//TRIM(LABEL)//'.dat'
        FluxMetaFile = 'OUTPUT/fluxmeta_'//TRIM(LABEL)//'.dat'
        CALL OpenFile_wStream(FluxUnit,FluxFile);       CLOSE(FluxUnit)
        CALL OpenFile_wSeq(FluxMetaUnit,FluxMetaFile);  CLOSE(FluxMetaUnit)
      END IF

      IF ( ConcDataPrint ) THEN
        ConcFile     = 'OUTPUT/conc_'//TRIM(LABEL)//'.dat'
        ConcMetaFile = 'OUTPUT/concmeta_'//TRIM(LABEL)//'.dat'
        CALL OpenFile_wStream(ConcUnit,ConcFile);       CLOSE(ConcUnit)
        CALL OpenFile_wSeq(ConcMetaUnit,ConcMetaFile);  CLOSE(ConcMetaUnit)
      END IF
    END IF


777 FORMAT(10X,A)
  END SUBROUTINE cminor_initialize

  SUBROUTINE cminor_run()
    ! t0 Rate + Jac_CC ValPtr fill (box h0 / stats); host will call rates per cell instead
    
    ! --- compute reaction rates and jacobian at t0 ---
    ! Note: rhs of the ODE system is not computed at this point, however,
    ! for icon coupling: f_chem = BAT * Rate, J_chem = Jac_CC
    RateCnt = 0
    IF (combustion) THEN
      CALL ReactionRates( Y_initial, Rate )
    ELSE
      CALL ReactionRates( Tspan(1), Y_initial, Rate )
    END IF
    Y = MAX( ABS(InitValAct) , eps ) * SIGN( ONE , InitValAct )

    CALL Start_Timer(StartTimer)
    CALL CSR_2_Empty_ValPtr( Jac_CC, nD_spc, nDropletClasses )
    CALL Jacobian_CC_ValPtr(Jac_CC , BAT , A , Rate , Y)
    Out%npds = Out%npds + 1
    CALL End_Timer(StartTimer, TimeJac)

    ! Main loop (box)
    done = .FALSE.
    IF ( Simulation .AND. .NOT. ReduceOnly ) THEN
      IF ( FluxDataPrint ) THEN
        iStpFlux = 0
        FluxFile     = 'OUTPUT/flux_'//TRIM(LABEL)//'.dat'
        FluxMetaFile = 'OUTPUT/fluxmeta_'//TRIM(LABEL)//'.dat'
        CALL OpenFile_wStream(FluxUnit,FluxFile);       CLOSE(FluxUnit)
        CALL OpenFile_wSeq(FluxMetaUnit,FluxMetaFile);  CLOSE(FluxMetaUnit)
      END IF

      IF ( ConcDataPrint ) THEN
        ConcFile     = 'OUTPUT/conc_'//TRIM(LABEL)//'.dat'
        ConcMetaFile = 'OUTPUT/concmeta_'//TRIM(LABEL)//'.dat'
        CALL OpenFile_wStream(ConcUnit,ConcFile);       CLOSE(ConcUnit)
        CALL OpenFile_wSeq(ConcMetaUnit,ConcMetaFile);  CLOSE(ConcMetaUnit)
      END IF

      Tspan = [tBegin, tEnd]
      h0 = InitialStepSize( Jac_CC, Rate, Tspan(1), Y_initial, ROS%pow )

      DO
        CALL Integrate ( InitValAct   &
        &              , Temperature0 &
        &              , h0           &
        &              , Tspan        )

        IF (Tspan(2) == tEnd .OR. done) EXIT
        Tspan = [Tspan(2), Tspan(2)+StpNetCDF]

        IF ( Tspan(2) >= tEnd ) THEN
          TSpan(2) = tEnd
          done = .TRUE.
        END IF
      END DO

      CALL End_Timer(Timer_Finish, Time_Finish)
      CALL Output_Statistics
    END IF
  END SUBROUTINE cminor_run

  SUBROUTINE cminor_reduce()
#ifdef ISSA
    IF ( Reduction ) THEN
      IF ( ReduceOnly ) THEN
        FluxFile     = 'OUTPUT/flux_'//TRIM(LABEL)//'.dat'
        FluxMetaFile = 'OUTPUT/fluxmeta_'//TRIM(LABEL)//'.dat'
      END IF
      CALL RunISSAReduction()
    END IF
#endif
  END SUBROUTINE cminor_reduce

  SUBROUTINE cminor_finalize()
    WRITE(*,*); WRITE(*,*)
    CALL ShowMaxErrorCounter()
  END SUBROUTINE cminor_finalize

  !================================================================
  ! Phase helpers
  !================================================================

  SUBROUTINE Init_Combustion
    INTEGER :: i

    CALL Read_Elements( SysFile , SysUnit )
    CALL Read_Species ( SysFile , SysUnit )
    CALL Read_Reaction( SysFile , SysUnit )

    WRITE(*,*) 'done   ---->  Solve Gas Energy Equation '

    CALL Print_ChemFile   ( ReactionSystem , ChemFile , ChemUnit , .TRUE. )

    CALL GetSpeciesNames( ChemFile , y_name )
    CALL Read_ThermoData( SwitchTemp , DataFile , DataUnit , nspc )

    bGs(1) = 1
    bGs(2) = ns_GAS
    iGs = [(i, i=bGs(1),bGs(2))]
    iGs2 = iGs
    hasGasSpc = .TRUE.

    CALL make_ChemSys_1_to_nD_arrays()
    CALL Setup_ThirdBodyIndex
    CALL Setup_ReactionIndex

    ALLOCATE( InitValAct(ns_GAS+1) , InitValKat(ns_KAT) , y_emi(ns_GAS) , y_depos(ns_GAS))
    y_emi   = ZERO
    y_depos = ZERO

    ALLOCATE( GFE(nspc)    , DGFEdT(nspc)   &
    &       , DelGFE(nreac), DDelGFEdT(nreac) )
    GFE      = ZERO; DGFEdT    = ZERO
    DelGFE   = ZERO; DDelGFEdT = ZERO

    IF ( MWFile /= '' ) THEN
      CALL Read_MolecularWeights(MW,MWFile,MWUnit,nspc)

      rMW = [ONE / MW]

      ALLOCATE( MoleFrac(ns_GAS) , MassFrac(ns_GAS) )
      MoleFrac = ZERO; MassFrac  = ZERO

      MoleFrac = 1.0e-20_dp
      CALL Read_INI_file( InitFile , MoleFrac, InitValKat , 'GAS' , 'INITIAL' )

      Press_in_dyncm2 = Pressure0 * Pa_to_dyncm2

      MassFrac = MoleFr_To_MassFr  ( MoleFrac )
      MoleConc = MoleFr_To_MoleConc( MoleFrac                &
      &                            , Press = Press_in_dyncm2 &
      &                            , Temp  = Temperature0    )
    ELSE
      WRITE(*,*)
      WRITE(*,777) '    No molecular weights are given.  '
      WRITE(*,777) '    Make sure the initial values are given in [mol/cm3] !'
      WRITE(*,*)
      ALLOCATE( MoleConc(ns_GAS) )
      MoleConc = 1.0e-20_dp
      CALL Read_INI_file( InitFile , MoleConc, InitValKat , 'GAS' , 'INITIAL' )
    END IF

    rho  = Density( MoleConc )
    rRho = ONE/rho
    InitValAct(1:ns_GAS) = MoleConc
    InitValAct(ns_GAS+1) = Temperature0

777 FORMAT(10X,A)
  END SUBROUTINE Init_Combustion

  SUBROUTINE Init_Atmosphere
    INTEGER :: i

    CALL ReadSystem( SysFile )
    WRITE(*,*) 'done'

    CALL Print_ChemFile( ReactionSystem, ChemFile, ChemUnit, .FALSE. )

    IF ( nr_special > 0 ) THEN
      CALL initf( nr_special )
      DO i = 1,nr_special
        CALL parsef ( i, ReactionSystem(iR%iSPECIAL(i))%Special%Formula    &
        &              , ReactionSystem(iR%iSPECIAL(i))%Special%cVariables )
      END DO
    END IF

    CALL InputChemicalData( InitFile , DataFile )
  END SUBROUTINE Init_Atmosphere

  SUBROUTINE Init_State
    ! FO/SO/HO maps, ODE dims, print, Atol, Y0
    CALL Setup_SpeciesOrder(A)

    nsr  = nspc  + nreac
    nsr2 = nspc2 + nreac2

    IF ( combustion ) THEN
      nDIM  = nspc + 1
      nDIM2 = nspc2 + 1
    ELSE IF ( adiabatic_parcel ) THEN
      nDIM  = nspc + 5
      nDIM2 = nspc2 + nDropletClasses + 4
    ELSE
      nDIM  = nspc
      nDIM2 = nspc2
    END IF
    rNspc = ONE / REAL(nspc2, KIND=dp)

    CALL Print_Run_Param
    WRITE(*,777) 'Initial values:   '
    WRITE(*,*)
    IF ( combustion ) fmt0 = '  [mol/cm3]'
    IF ( hasGasSpc )  WRITE(*,798) 'gaseous', SUM(InitValAct( bGs(1):bGs(2) )) , fmt0
    IF ( hasAquaSpc ) WRITE(*,798) 'aqueous', SUM(InitValAct( bAs2(1):bAs2(2) )) , fmt0
    WRITE(*,800) SUM(y_emi) , '  [molec/cm3/s]'
    IF ( combustion ) THEN
      WRITE(*,801) Temperature0
      WRITE(*,802) Pressure0
      WRITE(*,803) rho
    END IF
    WRITE(*,*)

    ALLOCATE( Y(nspc2), Rate(nreac2), ATolAll(nDIM2), Y_initial(nDIM2) )
    ATolAll = 1.0_dp
    IF ( hasGasSpc )  ATolAll(iGs2) = ATolGas
    IF ( hasAquaSpc ) ATolAll(iAs2) = ATolAqua
    IF ( combustion ) ATolAll(nDIM2) = ATolTemp
    IF ( adiabatic_parcel ) THEN
      ATolAll(iAqMassEq2) = ATolWaterMass
      ATolAll(iqEq2)      = ATolq
      ATolAll(izEq2)      = ATolz
      ATolAll(iTeq2)      = ATolTemp
      ATolAll(iRhoEq2)    = ATolRho
    END IF

    IF ( combustion ) THEN
      Y_initial(1:nspc2) = InitValAct(1:nspc2)
      Y_initial(nDIM)    = Temperature0
    ELSE IF ( adiabatic_parcel ) THEN
      Y_initial(1:nspc2)    = InitValAct
      Y_initial(iAqMassEq2) = DropletClasses%waterMass
      Y_initial(iTeq2)      = T_parcel
      Y_initial(iqEq2)      = q_parcel
      Y_initial(iRhoEq2)    = rho_parcel
      Y_initial(iZeq2)      = z_parcel
    ELSE
      Y_initial = InitValAct
    END IF

777 FORMAT(10X,A)
798 FORMAT(10X,'    Sum Initval (',A7,')      =  ', Es8.2, A)
800 FORMAT(10X,'    Sum Emissions (gaseous)    =  ', Es8.2, A)
801 FORMAT(10X,'    Temperature                =  ', Es8.2,'  [K]')
802 FORMAT(10X,'    Pressure                   =  ', Es8.2,'  [Pa]')
803 FORMAT(10X,'    Reactor density            =  ', Es8.2,'  [kg/cm3]')
  END SUBROUTINE Init_State

  SUBROUTINE Init_NetCDF_File
    IF ( NetCdfPrint ) THEN
      CALL Start_Timer(TimerNetCDF)
      CALL InitNetcdf
      IF (ALLOCATED(DropletClasses%Number)) THEN
        nDropletsNCDF = DropletClasses%Number
        LWCsNCDF      = LWC_array(Tspan(1))
      ELSE
        nDropletsNCDF = dummyNumber
        LWCsNCDF      = dummyLWC
      END IF
      CALL SetOutputNCDF( NetCDF,             &
                        & Tspan(1),           &
                        & ZERO,               &
                        & InitValAct,         &
                        & Temperature0,       &
                        & rho_parcel,         &
                        & q_parcel,           &
                        & pressure0,          &
                        & z_parcel,           &
                        & RH,                 &
                        & nDropletsNCDF,      &
                        & LWCsNCDF            )
      CALL StepNetCDF( NetCDF )
      CALL End_Timer(TimerNetCDF, TimeNetCDF)
    END IF
  END SUBROUTINE Init_NetCDF_File

  SUBROUTINE Init_Symbolic_LU
    INTEGER :: i, j, jj

    temp_LU_Dec = CSR_to_SpRowColD(Miter)

    IF ( Ordering ) THEN
      ALLOCATE(temp_LU_Dec%Restr(Miter%m))
      temp_LU_Dec%Restr = 0
      IF (hasAquaSpc .AND. nDropletClasses>1) THEN
        DO i=bHr(1), bHr(2)
          DO jj = A%RowPtr(i), A%RowPtr(i+1)-1
            j = A%ColInd(jj)
            IF (.NOT. nD_spc(j)) THEN
              temp_LU_Dec%Restr(j) = -1
            END IF
          END DO
        END DO
      END IF
      CALL SymbLU_SpRowColD_M( temp_LU_Dec )
    ELSE
      ALLOCATE(PivOrder(temp_LU_Dec%n))
      PivOrder = -90
      PivOrder(     1 : nDIM     ) = [(i , i = 1     , nDim )]
      CALL SymbLU_SpRowColD(temp_LU_Dec , PivOrder)
    END IF

    LU_Miter = RowColD_to_CSR( temp_LU_Dec )
    CALL Get_LU_Permutation( LU_Perm, LU_Miter, Miter )
    CALL Get_ValPtr_Permutations(  LU_Perm_ValPtr, w_InvColInd, bPermu, bInvPermu, bPtr, &
                                  &  LU_Miter, Miter, LU_Perm, nDropletClasses, nDIM2)
  END SUBROUTINE Init_Symbolic_LU

  SUBROUTINE Write_Sparse_Matrices
    IF (MatrixPrint) THEN
      CALL WriteSparseMatrix(A,'MATRICES/alpha_'//TRIM(LABEL), nreac, nspc)
      CALL WriteSparseMatrix(B,'MATRICES/beta_'//TRIM(LABEL), nreac, nspc)
      CALL WriteSparseMatrix(tmpJacCC,'MATRICES/JAC_'//TRIM(LABEL), nreac, nspc)
      CALL WriteSparseMatrix(BA,'MATRICES/BA_'//TRIM(LABEL), nreac, nspc)
      CALL WriteSparseMatrix(Miter,'MATRICES/Miter_'//TRIM(LABEL), nreac, nspc)
      CALL WriteSparseMatrix(LU_Miter,'MATRICES/LU_Miter_'//TRIM(LABEL), nreac, nspc)
      WRITE(*,777,ADVANCE='NO') '  Continue? [y/n]';  READ(*,*) inpt
      IF (inpt/='y') STOP
    END IF
777 FORMAT(10X,A)
  END SUBROUTINE Write_Sparse_Matrices

  ! SUBROUTINE Init_ReactionRates_t0
  !   RateCnt = 0
  !   IF (combustion) THEN
  !     CALL ReactionRates( Y_initial, Rate )
  !   ELSE
  !     CALL ReactionRates( Tspan(1), Y_initial, Rate )
  !   END IF
  !   Y = MAX( ABS(InitValAct) , eps ) * SIGN( ONE , InitValAct )

  !   CALL Start_Timer(StartTimer)
  !   CALL CSR_2_Empty_ValPtr( Jac_CC, nD_spc, nDropletClasses )
  !   CALL Jacobian_CC_ValPtr(Jac_CC , BAT , A , Rate , Y)
  !   Out%npds = Out%npds + 1
  !   CALL End_Timer(StartTimer, TimeJac)
  ! END SUBROUTINE Init_ReactionRates_t0

END MODULE Cminor_Driver_Mod
