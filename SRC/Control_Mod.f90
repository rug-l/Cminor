!
! Copyright (C) 2025 Levin Rug, Willi Schimmel (E-Mail: l.rug@lmu.de)
! See ./SRC/Cminor.f90 for the copyright notice
! See ./LICENSE for license information
! SPDX-License-Identifier: GPL-3.0
!
!-------------------------------------------------------------
!--  Modul for the Saving Control Values
!-------------------------------------------------------------
 MODULE Control_Mod
   USE Kind_Mod, ONLY: dp
   USE iso_fortran_env
!
!-----------------------------------------------------------------
!---  Scenario
!-----------------------------------------------------------------
!
!--- Identifier for scenario
      CHARACTER(80) :: LABEL        = ''      ! Identifier for scenario

!--- system clock rate and max count
      INTEGER(8)  :: clock_rate, clock_maxcount

!--- Files
      CHARACTER(80) :: RunFile      = ''          & ! Simulation data file
&                    , SysFile      = ''          & ! Chemical mechanism
&                    , ChemFile     = ''          & ! Chemical mechanism in .chem format
&                    , InitFile     = ''          & ! Initial concentrations
&                    , DataFile     = ''          & ! Gas and Aqueous data
&                    , MWFile       = ''          & ! molecular weights of species 
&                    , NetcdfFile   = ''          & ! NetCDF output file
&                    , ODEsolver    = ''            ! Method for Rosenbrock Integration

      CHARACTER(80) :: FluxMetaFile = ''         ! meta data for unformatted flux data
      CHARACTER(80) :: FluxFile     = ''         ! flux data (unformatted)
      CHARACTER(80) :: ConcMetaFile = ''         ! meta data for unformatted concentration data
      CHARACTER(80) :: ConcFile     = ''         ! concentration data (unformatted)

!
!--- Unit Numbers
      INTEGER, PARAMETER :: RunUnit      = 101  & 
&                         , MetUnit      = 102  & 
&                         , SysUnit      = 103  & 
&                         , ChemUnit     = 104  & 
&                         , MWUnit       = 105  & 
&                         , InitUnit     = 106  & 
&                         , DataUnit     = 107  & 
&                         , FluxMetaUnit = 109  &
&                         , FluxUnit     = 110  & 
&                         , ConcMetaUnit = 212  &
&                         , ConcUnit     = 213
!
!-- Set Levels and Parameters for Processes
      INTEGER  :: nDropletClasses  ! number of droplet classes
      ! arrays to mark and guide through vector-valued species and reactions arrays 
      ! in case of multiple droplet classes (when nDropletClasses>1)
      ! see: SUBROUTINE make_ChemSys_1_to_nD_arrays
      LOGICAL, ALLOCATABLE :: nD_spc(:), nD_KAT(:), nD_reac(:)
      INTEGER, ALLOCATABLE :: nD_Ptr_spc(:), nD_Ptr_KAT(:), nD_Ptr_reacs(:)


      REAL(dp) :: LWCLevelmin       & ! Lower level for LWC
&               , LWCLevelmax       & ! Upper level for LWC
&               , activation_radius & ! radius for droplets to be activated (below, eq radius is used always)
&               , updraft_velocity  & ! assumed updraft velocity (for adiabatic_parcel processes)
&               , RH0               & ! initial relative humidity (for adiabatic_parcel processes)
&               , RH                & ! current relative humidity (for adiabatic_parcel processes)
&               , q_parcel          & ! amount of water vapor in box (for adiabatic_parcel processes)
&               , rho_parcel        & ! density of parcel in kg/m3
&               , rho_parcel0       & ! initial density of parcel in kg/m3
&               , m_parcel          & ! mass of parcel in kg (for adiabatic_parcel processes)
&               , V_parcel0=1.0_dp  & ! volume of parcel at beginning in m3
&               , z_parcel          & ! current height of parcel (for adiabatic_parcel processes)
&               , Temperature0      & ! Initial temperature
&               , T_parcel          & ! temperature of parcel (for adiabatic_parcel processes)
&               , Pressure0         & ! Initial pressure
&               , Pressure            ! current pressure


!-- print matrices of the system (alpha, beta, miter, lu_miter, permutvector,..)
      LOGICAL :: MatrixPrint        & ! print certain matrices to file
&              , NetCdfPrint        & ! print out concs and other stuff to netcdf file
&              , DropletClassPrint  & ! print out concs and other stuff to netcdf file
&              , constLWC           & ! true if lwc value is fixed
&              , combustion         & ! if combustion=.TRUE. simulate combustion mechanism
&              , pHSet              & ! Initial pH by charge balance (1=on, 0=off)
&              , adiabatic_parcel   & ! trigger for adiabatic_parcel processes (requires updraft velocity)
&              , WaitBar            & ! ladebalken im terminal bei simulation (=1, default=0)
&              , FluxDataPrint      & ! writing flux data and analyse after simulaiton -> print new reaction file
&              , ConcDataPrint      & ! writing concentration data and analyse after simulaiton -> print new reaction file
&              , Simulation           ! calculation of species concentration 

      INTEGER :: Error_Est         ! error estimation 1 = inf norm  , 2 = euklid norm
    
!--- Output of reaction fluxes/concentrations/life times
      REAL(dp)                    :: StpFlux  & ! time interval to capture flux
                                   , StpConc    ! time interval to capture concentrations

      INTEGER                     :: iStpFlux=0 & ! number of written flux steps
&                                  , iStpConc=0   ! number of written concentration steps

!-----------------------------------------------------------------
!---  Times
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!--- Times:  < 0 in seconds, > 0 in hours, = 0  meteorology
      REAL(dp) :: tBegin              & ! Model start time
&               , tEnd                & ! Model end time
&               , StpNetcdf             ! Time step for Netcdf output

!--- NetCDF globals      
      INTEGER, ALLOCATABLE :: iNcdfGas(:), iNcdfAqua(:), iNcdfEmiss(:)
      INTEGER, ALLOCATABLE :: iNcdfGas2(:), iNcdfAqua2(:), iNcdfEmiss2(:)
      INTEGER, ALLOCATABLE :: iNCout_G(:), iNCout_A_l(:), iNCout_A_m3(:), iNCout_E(:)
      INTEGER, ALLOCATABLE :: iNCout_A_l_Ptr(:), iNCout_A_m3_Ptr(:)
      INTEGER  :: nNcdfGas=0, nNcdfAqua=0, nNcdfEmiss=0
      INTEGER  :: nNcdfGas2=0, nNcdfAqua2=0, nNcdfEmiss2=0
!

!---  Photolysis
      INTEGER :: iDate       ! Current date
      REAL(dp) :: rlat, rlon ! Latitude, Longitude

!---  Liquid water content parameter
      REAL(dp) :: LWCconst=2.0d-8       ! [l/m3] no cloud 

!--   dust factor for damping of photolysis rates, measured JNO2
      REAL(dp) :: Dust = 1.0d0

!--   minimal and maximal numerical step size
      REAL(dp) :: minStp = 1.0d-25
      REAL(dp) :: maxStp = 100.0d0

!-- initialize timers and time variables
      INTEGER(8) :: Timer_Start,   &
                  & Timer_Finish,  &
                  & TimerNetCDF,   &
                  & TimerRates,    &
                  & TimerJacobian, &
                  & TimerSymbolic, &
                  & Timer_Read,    &
                  & TimerNumeric

      REAL(dp) :: Tspan(2)
      REAL(dp) :: Tspan_tot(2)

      REAL(dp) :: Time_Finish=0.0d0
      REAL(dp) :: Time_Read=0.0d0
      REAL(dp) :: TimeRates=0.0d0
      REAL(dp) :: TimeSymbolic=0.0d0
      REAL(dp) :: TimeFac=0.0d0
      REAL(dp) :: TimeSolve=0.0d0
      REAL(dp) :: TimeJac=0.0d0
      REAL(dp) :: TimeNetCDF=0.0d0
      REAL(dp) :: TimeErrCalc=0.0d0
      REAL(dp) :: TimeFluxWrite=0.0d0
      REAL(dp) :: TimeFluxRead=0.0d0
      REAL(dp) :: TimeRhsCalc=0.0d0
      REAL(dp) :: TimeConcWrite=0.0d0
      REAL(dp) :: TimeConcRead=0.0d0
      REAL(dp) :: TimeSetValues=0.0d0
      REAL(dp) :: TimeIntegration=0.0d0
      REAL(dp) :: TimeNumeric=0.0d0
!

!--- type for some statistics
      TYPE Output_T
        REAL(dp), ALLOCATABLE :: y(:)  ! y-vector at Tend
        INTEGER :: nsteps     = 0      ! # succ. steps
        INTEGER :: nfailed    = 0      ! # failed steps
        INTEGER :: nRateEvals = 0      ! # Rate evaluation
        INTEGER :: npds       = 0      ! # Jacobian evaluation
        INTEGER :: ndecomps   = 0      ! # LU factorisation
        INTEGER :: nsolves    = 0      ! # solved lin algebra
        REAL(dp) :: Ttimestep = 0.0d0  ! mean Time for one ROW step
      END TYPE Output_T

      TYPE(Output_T) :: Out
!-----------------------------------------------------------------
!---  Numerics
!-----------------------------------------------------------------
!
!
!--- Tolerances for ROW Scheme
      REAL(dp) :: RtolROW        & ! Relative tolerance for ROW method
&               , AtolGas        & ! Absolute tolerance for gas phase
&               , AtolAqua       & ! Absolute tolerance for liquid phase
&               , AtolTemp       & ! Absolute tolerance for Temperatur
&               , AtolWaterMass  & ! Absolute tolerance for water masses
&               , Atolq          & ! Absolute tolerance for mixing ratio
&               , Atolz          & ! Absolute tolerance for height z
&               , AtolRho          ! Absolute tolerance for parcel density

      LOGICAL  :: Ordering     ! default = .TRUE. (Markowitz)

      REAL(dp), ALLOCATABLE :: ATolAll(:)

!
!-----------------------------------------------------------------
!---  Set Constants and Unit Conversion
!-----------------------------------------------------------------
!
!---  constants
    REAL(dp), PARAMETER :: HOUR      = 3600.0              &
&                        , hourday   = 86400.0             &
&                        , secday    = 4.32d04             &
&                        , mONE      = -1.d0               &
&                        , Pi        = 4.0d0*ATAN(1.0d0)   &
&                        , DR        = Pi / 180.d0         &  
&                        , PiHalf    = 2.0d0*ATAN(1.0d0)   & 
&                        , Pi34      = 3.0d0/4.0d0/Pi      & 
&                        , Pi43      = 4.0d0*Pi/3.0d0      & 
&                        , eps       = EPSILON(1.0d0)      &  ! such that 1+eps>1 with working precision
&                        , epsL      = 1000*EPSILON(1.0d0) &  ! such that 1+eps>1 with working precision
&                        , small     = TINY(1.0d0)         &  ! smallest pos. real value
&                        , big       = HUGE(1.0d0)            ! largest pos. real value

!
!--- Real number constants
  REAL(dp), PARAMETER ::  ZERO    =     0.0d0   & 
&                      ,  ONE     =     1.0d0   & 
&                      ,  TWO     =     2.0d0  ,   rTWO     =   ONE/tWO     &
&                      ,  THREE   =     3.0d0  ,   rTHREE   =   ONE/THREE   &
&                      ,  FOUR    =     4.0d0  ,   rFOUR    =   ONE/FOUR    &
&                      ,  FIVE    =     5.0d0  ,   rFIVE    =   ONE/FIVE    &
&                      ,  SIX     =     6.0d0  ,   rSIX     =   ONE/SIX     &
&                      ,  SEVEN   =     7.0d0  ,   rSEVEN   =   ONE/SEVEN   &
&                      ,  EIGHT   =     8.0d0  ,   rEIGHT   =   ONE/EIGHT   &
&                      ,  NINE    =     9.0d0  ,   rNINE    =   ONE/NINE    &
&                      ,  TEN     =    10.0d0  ,   rTEN     =   ONE/TEN     &
&                      ,  ELEVN   =    11.0d0  ,   rELEVN   =   ONE/ELEVN   &
&                      ,  TWELV   =    12.0d0  ,   rTWELV   =   ONE/TWELV   &
&                      , TWENTY   =    20.0d0  ,   rTWENTY  =   ONE/TWENTY  &
&                      , mTHIRTY  =   -30.0d0  &
&                      , rm300    = mONE/300.d0,   r300     =   ONE/300.d0
!
!--- Orders of magnitude
  REAL(dp), PARAMETER :: nano     =     1.0d-09    &
&                      , micro    =     1.0d-06    &
&                      , milli    =     1.0d-03    &
&                      , kilo     =     1.0d+03    &
&                      , mega     =     1.0d+06    &
&                      , giga     =     1.0d+09
!
!--- Natural logarithms
  REAL(dp), PARAMETER ::   ln10   =     LOG(TEN)    &
&                      ,  rln10   = ONE/LOG(TEN)
!
!--- minimum values if there is no sun
  REAL(dp), PARAMETER ::    EyChiZmin  =  9.357d-14
!
!
!--- Unit Conversion constants
!
!     Pressure
  REAL(dp), PARAMETER :: bar_to_dyncm2 = 1.0d06
  REAL(dp), PARAMETER :: dyncm2_to_Pa  = 1.0d-01
  REAL(dp), PARAMETER :: Pa_to_dyncm2  = 1.0d+01
  REAL(dp), PARAMETER :: bar_to_Pa     = 1.0d05
  REAL(dp), PARAMETER :: atm_to_Pa     = 101325.0d0
!
!     Energy
  REAL(dp), PARAMETER :: cal_to_joule  = 4.184d0
  REAL(dp), PARAMETER :: joule_to_cal  = ONE / cal_to_joule
  REAL(dp), PARAMETER :: joule_to_kcal = milli * joule_to_cal
  REAL(dp), PARAMETER :: kcal_to_joule = kilo * cal_to_joule
  REAL(dp), PARAMETER :: joule_to_erg  = TEN * mega
  REAL(dp), PARAMETER :: erg_to_joule  = rTEN * micro
!

!--- Unit Conversion
    INTEGER :: GasUnit, AquaUnit, GasRateUnit

    INTEGER :: UnitGas=0
    INTEGER :: UnitAqua=0


!-----------------------------------------------------------------
!---  Output control 
!-----------------------------------------------------------------

!
!-- lengths of input for chemical reaction mechnism
    CHARACTER(3), PARAMETER :: LenLineChar = '400'
    INTEGER, PARAMETER :: LenLine=400
    INTEGER, PARAMETER :: LenName=100
    INTEGER, PARAMETER :: LenType=20

!-- Type declaration 
    TYPE List
      INTEGER,   ALLOCATABLE  :: List(:)
      Real(dp),  ALLOCATABLE  :: ListE(:)
      Real(dp),  ALLOCATABLE  :: ListP(:)
      INTEGER                :: len
    END TYPE List

    TYPE Chain
      CHARACTER(100), ALLOCATABLE :: sName(:,:)
      INTEGER,        ALLOCATABLE :: sIdx(:,:)
      Type(List),     ALLOCATABLE :: rIdx(:)
    END TYPE Chain

    INTEGER, ALLOCATABLE :: maxErrorCounter(:)

    CONTAINS

      SUBROUTINE Start_Timer(StartTime)
        INTEGER(8) :: StartTime

        ! saves current time point
        CALL SYSTEM_CLOCK(COUNT=StartTime)

      END SUBROUTINE Start_Timer

      SUBROUTINE End_Timer(StartTime, MeasuredTime)
        REAL(dp) :: MeasuredTime
        INTEGER(8) :: StartTime

        INTEGER(8)  :: EndTime

        CALL SYSTEM_CLOCK(COUNT=EndTime)

        MeasuredTime = MeasuredTime + REAL(EndTime - StartTime, dp) / REAL(clock_rate, dp)
      END SUBROUTINE End_Timer

      SUBROUTINE ask_for_continuing()
        CHARACTER(LenName) :: inpt
        WRITE(*,'(A)',ADVANCE='NO') '  Continue? [y = continue, else stop] ';  READ(*,*) inpt
        IF (inpt/='y') STOP
        WRITE(*,*) ''
      END SUBROUTINE

 END MODULE Control_Mod
