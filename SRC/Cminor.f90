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
PROGRAM Cminor
  !
  USE Kind_Mod,         ONLY: dp

  USE InitRoutines_Mod, ONLY: InitRun

  USE Rates_Mod,        ONLY: ReactionRates, RateCnt

  USE Sparse_Mod,       ONLY: SymbLU_SpRowColD_M, CSR_Matrix_T, SpRowColD_T, A, B, SparseAdd,   &
                            & BA, BAT, SparseID, RowColD_to_CSR, Miter, LU_Miter, Jac_CC,       &
                            & CSR_to_SpRowColD, Copy_CSR, BuildSymbolicClassicMatrix,           &
                            & CSR_2_Empty_ValPtr, Free_Matrix_CSR, Free_SpRowColD, SymbolicAdd, &
                            & Get_LU_Permutation, Get_ValPtr_Permutations, Jacobian_CC_ValPtr,  &
                            & SymbLU_SpRowColD, SymbolicMult, TransposeSparse,                  &
                            & WriteSparseMatrix

  USE Chemsys_Mod,      ONLY: Print_ChemFile, ReactionSystem, nsr, nsr2, PositionSpeciesAll,    &
                            & InputChemicalData, Read_Diag, make_ChemSys_1_to_nD_arrays,        &
                            & Read_DiagEmiss, Read_INI_file, ReadSystem, Setup_SpeciesOrder

  USE Integration_Mod,  ONLY: Integrate

  USE Rosenbrock_Mod,   ONLY: SetRosenbrockMethod, ROS, LU_Perm, LU_Perm_ValPtr, w_InvColInd,   &
                            & bPermu, bInvPermu, bPtr, InitialStepSize

  USE Control_Mod,      ONLY: nD_spc, List, Tspan, ZERO, clock_rate, clock_maxcount, SysUnit,   &
                            & TimeSymbolic, TimeNetCDF, timer_finish, Start_Timer, End_Timer,   &
                            & TimeJac, time_read, combustion, tBegin, tEnd, T_parcel, q_parcel, &
                            & Simulation, StpNetCDF, Temperature0, SysFile, RunFile, AtolGas,   &
                            & Pressure0, Pa_to_dyncm2, Ordering, ONE, Out, adiabatic_parcel,    &
                            & ODEsolver, NetCDFprint, nDropletClasses, MWFile, MWUnit, eps,     &
                            & maxErrorCounter, MatrixPrint, ATolAll, Atolq, AtolRho, Atolz,     &
                            & iStpFlux, InitFile, FluxUnit, FluxMetaUnit, AtolWaterMass, LABEL, &
                            & FluxFile, FluxMetaFile, FluxDataPrint, z_parcel, rho_parcel,      &
                            & DataUnit, DataFile, ConcUnit, ConcFile, clock_maxcount_int,       &
                            & ConcMetaUnit, ConcMetaFile, ConcDataPrint, ChemUnit, AtolTemp,    &
                            & AtolAqua, ChemFile, clock_rate_int

  USE Reac_Mod,         ONLY: InitValAct, y_name, bGs, y_emi, y_depos, iAs2, iGs2, nDIM2, rho,  &
                            & bHr, rnspc, bAs2, SwitchTemp, rRho, MW, rMW, PHOTO, ddelGFEdT,    &
                            & nspc2, nspc, ns_GAS, nr_special, nr_G_photo, GFE, dGFEdT, delGFE, &
                            & nr_A_photo, nreac, nreac2, iRhoeq2, iAqMassEq2, DropletClasses,   &
                            & nDIM, InitValKAT, iGs, ns_KAT, iTeq2, iqEq2, izeq2, hasAquaSpc,   &
                            & hasGasSpc, iR

  USE IO_Mod,           ONLY: Logo, OpenFile_wseq, OpenFile_wstream, Output_Statistics,         &
                            & Print_Run_Param, ShowMaxErrorCounter, ConvertTime,                &
                            & Matrix_Statistics

  USE ChemKinInput_Mod, ONLY: Read_Elements, MoleFr_To_MassFr, MoleFr_To_MoleConc, Density,     &
                            & GetSpeciesNames, GetSpeciesNames, Read_MolecularWeights,          &
                            & Read_Reaction, Read_Species, Read_ThermoData,                     &
                            & Setup_ThirdBodyIndex, Setup_ReactionIndex

  USE NetCDF_Mod,       ONLY: InitNetCDF, SetOutputNCDF, StepNetCDF, NetCDF

  USE fparser,          ONLY: initf, parsef

  USE Meteo_Mod,        ONLY: Set_pseudoLWCbounds, LWCb


  IMPLICIT NONE
  !
  CHARACTER(80)   :: Filename0 = ''        ! *.run file
  !
  INTEGER  :: i,j,jj

  ! NetCDF stuff
  REAL(dp) :: StartTimer
  REAL(dp), ALLOCATABLE :: Y(:), Y_initial(:)
  ! reaction rate array + part. derv. rate over temperatur vector
  REAL(dp), ALLOCATABLE :: Rate(:)
  !

  ! converting mole to mass to conc
  REAL(dp), ALLOCATABLE   :: MoleFrac(:), MassFrac(:), MoleConc(:)
  REAL(dp)                :: Press_in_dyncm2

  ! matrices for symbolic phase
  TYPE(CSR_Matrix_T)     :: Id , tmpJacCC  ! compressed row
  TYPE(SpRowColD_T)      :: temp_LU_Dec    ! sparse-LU matrix format
  ! permutation vector/ pivot order for LU-decomp
  INTEGER, ALLOCATABLE :: PivOrder(:)

  ! format string
  CHARACTER(14) :: fmt0 = '  [molec/cm3] '

  CHARACTER(1) :: inpt=''
  LOGICAL :: done=.FALSE.

  ! step size
  REAL(dp) :: h0

  !
  !================================================================
  !                        MAIN Programm
  !================================================================
  !
  ! initialize clock for timing
  CALL SYSTEM_CLOCK(count_rate = clock_rate_int    )
  CALL SYSTEM_CLOCK(count_max  = clock_maxcount_int)
  clock_rate     = REAL(clock_rate_int    , dp)
  clock_maxcount = REAL(clock_maxcount_int, dp)

  !----------------------------------------------------------------
  ! --- Print the program logo  
  CALL Logo

  !----------------------------------------------------------------
  ! --- Read run control parameters (which runfile)
  CALL getarg( 1 , FileName0 )             
  IF ( FileName0 == '' ) THEN
    WRITE(*,777,ADVANCE='NO') 'Input RUNFilename: '; READ(*,*)   FileName0
  END IF
  RunFile = TRIM(ADJUSTL(FileName0))

  !================================================================
  !===                     Initialization
  !================================================================
  !
  !----------------------------------------------------------------
  ! --- Read run-file
  CALL InitRun
  WRITE(*,777) 'Initialize run-file .......... done'
 
  !----------------------------------------------------------------
  ! --- set cloud intervall
  LWCb = Set_pseudoLWCbounds()

  !----------------------------------------------------------------
  !  --- read the .sys data, save coefs in sparse matrix
  CALL Start_Timer(Time_Read)
  CALL Start_Timer(Timer_Finish)

  WRITE(*,777,ADVANCE='NO') 'Reading sys-file .............'

  IF ( combustion ) THEN

    CALL Read_Elements( SysFile , SysUnit )
    CALL Read_Species ( SysFile , SysUnit )
    CALL Read_Reaction( SysFile , SysUnit )

    WRITE(*,*) 'done   ---->  Solve Gas Energy Equation '

   !-----------------------------------------------------------------------
   ! --- print reactions and build A, B and (B-A) structure
   !-----------------------------------------------------------------------
    CALL Print_ChemFile   ( ReactionSystem , ChemFile , ChemUnit , .TRUE. )

    CALL GetSpeciesNames( ChemFile , y_name )
    CALL Read_ThermoData( SwitchTemp , DataFile , DataUnit , nspc )

    ! set boundaries for gaseous species
    bGs(1) = 1
    bGs(2) = ns_GAS
    iGs = [(i, i=bGs(1),bGs(2))]
    iGs2 = iGs
    hasGasSpc = .TRUE.
    
    !-----------------------------------------------------------------------
    ! --- create guide arrays for droplet classes and aqueous species and reacs
    ! --- (these arrays are trivial but (unfortunately) needed if no aqueous phase is present,
    ! --- somewhen they should make it possible to say nDropletClasses=0 to cut aqueous
    ! --- processes out of a mechanism)
    CALL make_ChemSys_1_to_nD_arrays()

    !--- gather correct indices cause ChemKin mechanisms read in unsorted
    CALL Setup_ThirdBodyIndex
    CALL Setup_ReactionIndex
   
    !--- Read initial values
    ALLOCATE( InitValAct(ns_GAS+1) , InitValKat(ns_KAT) , y_emi(ns_GAS) , y_depos(ns_GAS))
    y_emi   = ZERO
    y_depos = ZERO 

    !--- malloc gibbs energy, derivates
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

    ! Initialising reactor density
    rho  = Density( MoleConc )
    rRho = ONE/rho       ! in [cm3/kg]
    InitValAct(1:ns_GAS) = MoleConc
    InitValAct(ns_GAS+1) = Temperature0

  ELSE ! atmospheric mechanism

    CALL ReadSystem( SysFile )
    WRITE(*,*) 'done  ---->  Fix Temperature'

    !-----------------------------------------------------------------------
    ! --- print reactions and build A, B and (B-A) structure
    !-----------------------------------------------------------------------
    CALL Print_ChemFile( ReactionSystem, ChemFile, ChemUnit, .FALSE. )

    !-----------------------------------------------------------------------
    ! --- initialize fpraser for reactions with special rate formula
    !-----------------------------------------------------------------------
    IF ( nr_special > 0 ) THEN
     
      ! Initialize function parser for n special functions
      CALL initf( nr_special ) 
      
      ! Parse and bytecompile ith function string 
      DO i = 1,nr_special
        CALL parsef ( i, ReactionSystem(iR%iSPECIAL(i))%Special%Formula    &
        &              , ReactionSystem(iR%iSPECIAL(i))%Special%cVariables )
      END DO
    END IF

    !-----------------------------------------------------------------------
    ! --- Input of initial data and thermodynamic properties
    CALL InputChemicalData( InitFile , DataFile )

  END IF

  WRITE(*,777) 'Reading ini-file ............. done'
  WRITE(*,777) 'Printing chem-file ........... done'
  WRITE(*,*)
  IF ( nDropletClasses<=1 ) THEN
    WRITE(*,'(10X,A,I6)') '    Number of Reactions = ', nreac
    WRITE(*,'(10X,A,I6)') '    Number of Species   = ', nspc
  ELSE
    WRITE(*,'(10X,A,I6,A4,I7,A1)') '    Number of Reactions (with droplet classes) = ', nreac, '   (', nreac2, ')'
    WRITE(*,'(10X,A,I6,A4,I7,A1)') '    Number of Species   (with droplet classes) = ', nspc , '   (', nspc2 , ')'
  END IF

  !-----------------------------------------------------------------------
  ! --- Read species for diagnose (print species concs to NetCDF file)
  CALL Read_Diag( NetCDF%spc_Pos , NetCDF%spc_Phase , InitFile )
  !-----------------------------------------------------------------------

  !--------------------------------------------------------------------
  ! --- Read Emissions for diagnose (print emission values to NetCDF file)
  CALL Read_DiagEmiss( NetCDF%Emiss_Pos, InitFile )
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! --- Timer
  CALL End_Timer(Time_Read)
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! --- this is for the new mass action product routine 
  CALL Setup_SpeciesOrder(A)
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! --- Dimension initialization for the unknowns and matrices
  !-----------------------------------------------------------------------
  !
  nsr = nspc + nreac
  nsr2 = nspc2 + nreac2

  IF ( combustion ) THEN
    nDIM  = nspc+1
    nDIM2 = nspc2+1
  ELSE IF ( adiabatic_parcel ) THEN
    ! add equations for condensation of every droplet class and T, q, rho, z
    nDIM  = nspc + 5
    nDIM2 = nspc2 + nDropletClasses + 4
  ELSE
    nDIM  = nspc
    nDIM2 = nspc2
  END IF

  rNspc = ONE/REAL(nspc2,KIND=dp)  ! rNspc for error calculation

  ! true if mechanism contains photolytic reactions
  PHOTO = (nr_G_photo+nr_A_photo) > 0

  !-----------------------------------------------------------------------
  !--- print input parameter, method, tols, concs, temp, pressure,....
  !-----------------------------------------------------------------------
  CALL Print_Run_Param
  
  WRITE(*,777) 'Initial values:   '
  WRITE(*,*)
    
  IF ( combustion ) fmt0 = '  [mol/cm3]'
  IF ( hasGasSpc )   WRITE(*,798) 'gaseous', SUM(InitValAct( bGs(1):bGs(2) )) , fmt0
  IF ( hasAquaSpc )  WRITE(*,798) 'aqueous', SUM(InitValAct( bAs2(1):bAs2(2) )) , fmt0
  WRITE(*,800) SUM(y_emi) , '  [molec/cm3/s]'

  IF ( combustion ) THEN
    WRITE(*,801) Temperature0
    WRITE(*,802) Pressure0
    WRITE(*,803) rho
  END IF
  WRITE(*,*)

  !-----------------------------------------------------------------------
  ! ---  Allocate arrays for some varables
  !-----------------------------------------------------------------------
  ALLOCATE( Y(nspc2), Rate(nreac2) )

  ! All species get the same abs tolerance
  ALLOCATE( ATolAll(nDIM2) )
  ATolAll = 1.0_dp ! default

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

  !-----------------------------------------------------------------------
  ! --- Initialize NetCDF output file
  !-----------------------------------------------------------------------
  IF ( NetCdfPrint ) THEN
    CALL Start_Timer(TimeNetCDF)
    CALL InitNetcdf
    CALL SetOutputNCDF( NetCDF, Tspan(1), ZERO, InitValAct, Temperature0 )
    CALL StepNetCDF( NetCDF )
    CALL End_Timer(TimeNetCDF)
  END IF

  !-----------------------------------------------------------------------
  ! --- Symbolic phase 
  !    - generating sparse matrices of stoechiometic coeffs (beta-alpha)^T 
  !    - generating sparse matrices of Jacobian matrix
  !    - generating sparse matrices of LU decomposition of the Jacobian
  !-----------------------------------------------------------------------
  WRITE(*,777,ADVANCE='NO') 'Symbolic-phase................'
  CALL Start_Timer(TimeSymbolic)

  CALL SymbolicAdd( BA , B , A )      ! symbolic addition:    BA = B + A
  CALL SparseAdd  ( BA , B , A, '-' ) ! numeric subtraction:  BA = B - A
  CALL TransposeSparse( BAT , BA )    ! transpose BA:        BAT = Transpose(BA) 

  !---- Get Rosenbrock Parameters
  CALL SetRosenbrockMethod( ROS , ODEsolver )

  ! create symbolic Jacobian
  CALL SymbolicMult( BAT , A , tmpJacCC )  
  Id = SparseID( nspc )                    
  CALL SymbolicAdd( Jac_CC , Id , tmpJacCC )
  CALL Free_Matrix_CSR( Id )

  ALLOCATE(maxErrorCounter(nDIM2))
  maxErrorCounter = 0

  ! Set symbolic structure of iteration matrix for Row-Method
  ! also assign constant matrix parts like (beta-alpha)^T, alpha
  CALL BuildSymbolicClassicMatrix( Miter, Jac_CC )

  ! Permutation given by Markowitz Ordering strategie
  temp_LU_Dec = CSR_to_SpRowColD(Miter) 

  IF ( Ordering ) THEN

    ! test for better permutation 
    ALLOCATE(temp_LU_Dec%Restr(Miter%m))
    temp_LU_Dec%Restr = 0
    ! go through henry reacs, mark non-vector (gaseous) counterparts as restricted
    ! this is needed for vector lu decomposition
    IF (hasAquaSpc .AND. nDropletClasses>1) THEN
      DO i=bHr(1), bHr(2) ! go through henry reactions
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

    ! Pivots bis nreac in reihenfolge 1,2,...,nreac
    PivOrder(     1 : nDIM     ) = [(i , i = 1     , nDim )]
    CALL SymbLU_SpRowColD(temp_LU_Dec , PivOrder)
  END IF

  ! converting back to csr format
  LU_Miter = RowColD_to_CSR( temp_LU_Dec )

  ! Get the permutation vector LU_Perm and map values of Miter
  ! to the permuted LU matrix
  CALL Get_LU_Permutation( LU_Perm, LU_Miter, Miter )

  ! Get ValPtr Permutations for vector-entried sparse matrix (aqueous chemistry)
  CALL Get_ValPtr_Permutations(  LU_Perm_ValPtr, w_InvColInd, bPermu, bInvPermu, bPtr, &
                                &  LU_Miter, Miter, LU_Perm, nDropletClasses, nDIM2)

  WRITE(*,*) 'done' ! Symbolic phase
  WRITE(*,*) 

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

  WRITE(*,777) 'Matrix Statistics: '
  WRITE(*,*) 
  CALL Matrix_Statistics(A,B,BAT,tmpJacCC,Miter,LU_Miter)

  CALL Free_Matrix_CSR( tmpJacCC )
  CALL Free_SpRowColD( temp_LU_Dec )

  CALL End_Timer(TimeSymbolic)

  ALLOCATE(Y_initial(nDIM2))
  IF ( combustion ) THEN
    Y_initial(1:nspc2)    = InitValAct(1:nspc2)
    Y_initial(nDIM)       = Temperature0
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

  ! ---- Calculate first reaction rates
  RateCnt = 0
  IF (combustion) THEN
    CALL ReactionRates( Y_initial, Rate )
  ELSE
    CALL ReactionRates( Tspan(1), Y_initial, Rate )
  END IF
  Y = MAX( ABS(InitValAct) , eps ) * SIGN( ONE , InitValAct )    ! |y| >= eps

  ! ---- Calculate values of Jacobian
  CALL Start_Timer(StartTimer)
  CALL CSR_2_Empty_ValPtr( Jac_CC, nD_spc, nDropletClasses )
  CALL Jacobian_CC_ValPtr(Jac_CC , BAT , A , Rate , Y)
  Out%npds = Out%npds + 1
  CALL End_Timer(TimeJac, StartTimer)

  IF ( Simulation ) THEN
    IF ( FluxDataPrint ) THEN
      FluxFile     = 'OUTPUT/flux_'//TRIM(LABEL)//'.dat'
      FluxMetaFile = 'OUTPUT/fluxmeta_'//TRIM(LABEL)//'.dat'
    END IF

    ! open file to save the fluxes 
    IF ( FluxDataPrint ) THEN
      iStpFlux = 0
      CALL OpenFile_wStream(FluxUnit,FluxFile);       CLOSE(FluxUnit)
      CALL OpenFile_wSeq(FluxMetaUnit,FluxMetaFile);  CLOSE(FluxMetaUnit)
    END IF
    
    IF ( ConcDataPrint ) THEN
      ConcFile     = 'OUTPUT/conc_'//TRIM(LABEL)//'.dat'
      ConcMetaFile = 'OUTPUT/concmeta_'//TRIM(LABEL)//'.dat'
    END IF

    ! open file to save the fluxes 
    IF ( ConcDataPrint ) THEN
      iStpFlux = 0
      CALL OpenFile_wStream(ConcUnit,ConcFile);       CLOSE(ConcUnit)
      CALL OpenFile_wSeq(ConcMetaUnit,ConcMetaFile);  CLOSE(ConcMetaUnit)
    END IF

    !-----------------------------------------------------------------------
    ! --- Start the integration routine 
    !-----------------------------------------------------------------------
    IF ( StpNetCDF < ZERO ) THEN
      Tspan = [tBegin, tEnd]
    ELSE 
      Tspan = [tBegin, tBegin+StpNetCDF]
    END IF

    !---- Calculate a first stepsize based on 2nd deriv.
    h0 = InitialStepSize( Jac_CC, Rate, Tspan(1), Y_initial, ROS%pow )

    DO
      CALL Integrate ( InitValAct   &  ! initial concentrations activ species
      &              , Temperature0 &  ! initial temperature
      &              , h0           &  ! reaction rates at time=t0
      &              , Tspan        )  ! integration invervall

      IF (Tspan(2) == tEnd .OR. done) EXIT
      Tspan = [Tspan(2), Tspan(2)+StpNetCDF]
      
      ! --- Hit end point exactly.
      IF ( Tspan(2) >= tEnd ) THEN
        TSpan(2) = tEnd
        done = .TRUE.
      END IF

    END DO

    ! --- stop timer and print output statistics
    CALL End_Timer(Timer_Finish)

    CALL Output_Statistics

  END IF
  
  WRITE(*,*); WRITE(*,*)

  CALL ShowMaxErrorCounter()

  !================================================================
  !==  FORMAT Statements
  !================================================================
  !
  777  FORMAT(10X,A)
  798  FORMAT(10X,'    Sum Initval (',A7,')      =  ', Es8.2, A)
  800  FORMAT(10X,'    Sum Emissions (gaseous)    =  ', Es8.2, A)
  801  FORMAT(10X,'    Temperature                =  ', Es8.2,'  [K]') 
  802  FORMAT(10X,'    Pressure                   =  ', Es8.2,'  [Pa]')
  803  FORMAT(10X,'    Reactor density            =  ', Es8.2,'  [kg/cm3]')
END PROGRAM Cminor
