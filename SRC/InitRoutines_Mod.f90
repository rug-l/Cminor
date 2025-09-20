!
! Copyright (C) 2025 Levin Rug, Willi Schimmel (E-Mail: l.rug@lmu.de)
! See ./SRC/Cminor.f90 for the copyright notice
! See ./LICENSE for license information
! SPDX-License-Identifier: GPL-3.0
!
MODULE InitRoutines_Mod   
   
   USE Kind_Mod, ONLY: dp

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE InitRun()
!==================================================
!===  Reading and Setting of Run Control Parameters
!==================================================
      USE Control_Mod, ONLY: RunFile, RunUnit, SysFile, ATolGas, ATolAqua, ATolTemp,    &
                           & LABEL, ChemFile, ConcDataPrint, ConstLWC, T_parcel, milli, &
                           & DataFile, Dust, Atolq, AtolRho, Atolz, AtolWaterMass,      &
                           & error_est, FluxDataPrint, iDate, rho_parcel, rho_parcel0,  &
                           & InitFile, LWCConst, updraft_velocity, q_parcel, V_parcel0, &
                           & LWCLevelMin, LWCLevelMax, MatrixPrint, MaxStp, WaitBar,    &
                           & minStp, MWFile, nDropletClasses, NetCdfFile, NetCdfPrint,  &
                           & ODEsolver, Ordering, RH0, RH, m_parcel, epsL, Pressure,    &
                           & phSet, pressure0, rlat, rlon, adiabatic_parcel, ZERO,      &
                           & rTolROW, Simulation, StpConc, StpFlux, StpNetCdf, Tspan,   &
                           & tBegin, tEnd, Temperature0, combustion, Tspan_tot,         &
                           & activation_radius, DropletClassPrint

      USE Meteo_Mod,   ONLY: R, molw_air, qsatw, H2O, RefRH, N2, O2, RefM_dry, RefO2,   &
                           & H2, RefH2, RefH2O, RefM, RefTemp, RefPressure, esatw,      &
                           & mol2part, RefN2, SI_Gas, alpha_H2O, beta_H2O

      INTEGER        :: io_stat
      CHARACTER(400) :: io_msg = ''

!-----------------------------------------------------------------
!--- NAMELISTS
      NAMELIST /SCENARIO/  LABEL ,     &
      &                    WaitBar , &
      &                    combustion , &
      &                    Simulation

      NAMELIST /FILES/  SysFile ,    &
      &                 DataFile ,   &
      &                 InitFile ,   &
      &                 MWFile

      NAMELIST /TIMES/  tBegin, tEnd

      NAMELIST /METEO/  LWCLevelmin ,      &
      &                 LWCLevelmax ,      &
      &                 nDropletClasses ,  &
      &                 activation_radius ,&
      &                 adiabatic_parcel , &
      &                 RH0 ,              &
      &                 updraft_velocity , &
      &                 pHSet ,            &
      &                 iDate ,            &
      &                 rlat ,             &
      &                 rlon ,             &
      &                 dust ,             &
      &                 Temperature0 ,     &
      &                 alpha_H2O,         &
      &                 beta_H2O,          &
      &                 Pressure0

      NAMELIST /NUMERICS/  RtolROW ,      &
      &                    AtolGas ,      &
      &                    AtolAqua ,     &
      &                    AtolTemp,      &
      &                    AtolWaterMass, &
      &                    Atolq,         &
      &                    Atolz,         &
      &                    AtolRho,       &
      &                    minStp ,       &
      &                    maxStp ,       &
      &                    ODEsolver ,    &
      &                    Error_Est ,    &
      &                    Ordering

      NAMELIST /OUTPUT/  NetCdfFile ,       &
      &                  StpNetcdf ,        &
      &                  StpFlux ,          &
      &                  StpConc ,          &
      &                  MatrixPrint,       &
      &                  DropletClassPrint, &
      &                  FluxDataPrint,     &
      &                  ConcDataPrint

!
!===================================================================
!===  Set and Read Simulation Values
!===================================================================
!
!--- Open run control file
      OPEN(UNIT=RunUnit,FILE=RunFile,IOSTAT=io_stat,IOMSG=io_msg)
      CALL ErrorCheck(io_stat,io_msg,'opening run-file')

!-----------------------------------------------------------------
!---  Scenario
!-----------------------------------------------------------------
!
!--- Set Default Values for SCENARIO Namelist

      LABEL = ''
      WaitBar  = .TRUE.
      combustion = .FALSE.
      Simulation = .TRUE.

!--- Read SCENARIO namelist
      READ(RunUnit,SCENARIO,IOSTAT=io_stat,IOMSG=io_msg)
      CALL ErrorCheck(io_stat,io_msg,'reading SCENARIO list')

!-----------------------------------------------------------------
!---  Files
!-----------------------------------------------------------------
      READ(RunUnit,FILES,IOSTAT=io_stat,IOMSG=io_msg)
      CALL ErrorCheck(io_stat,io_msg,'reading FILES list')
      
!--- Adjust Filenames
      CALL FileNameCheck(SysFile,'SysFile')
      CALL FileNameCheck(DataFile,'DataFile')
      CALL FileNameCheck(InitFile,'InitFile')
      ChemFile   = ADJUSTL(SysFile(:INDEX(SysFile,'.sys')-1)//'.chem')
      MWFile     = ADJUSTL(MWFile)

      IF ( TRIM(LABEL) == '' ) THEN
        LABEL = ADJUSTL(SysFile(INDEX(SysFile,'/')+1:INDEX(SysFile,'.sys')-1))
      ELSE
        LABEL = ADJUSTL(LABEL)
      END IF

!-----------------------------------------------------------------
!---  Times
!-----------------------------------------------------------------
!
!--- Set Default Values for TIMES Namelist

      tBegin = 0.0_dp
      tEnd   = 0.0_dp

!--- Read TIMES namelist
      READ(RunUnit,TIMES,IOSTAT=io_stat,IOMSG=io_msg)
      CALL ErrorCheck(io_stat,io_msg,'reading TIMES list')

      IF ( tBegin >= tEnd ) THEN
        WRITE(*,*) '  tBegin >= tEnd  '
        STOP 
      ELSE
        Tspan = [tBegin , tEnd]
        Tspan_tot = [tBegin , tEnd]
      END IF
!
!-----------------------------------------------------------------
!---  Meteorology
!-----------------------------------------------------------------
!
!--- Set Default Values for METEO Namelist
      pHSet             = .TRUE.
      adiabatic_parcel  = .FALSE.
      RH0               = RefRH
      nDropletClasses   = 1
      activation_radius = 1.0e-12_dp
      updraft_velocity  = 1.0_dp
      LwcLevelmin       = 2.0e-08_dp 
      LwcLevelmax       = 3.0e-04_dp 
      constLWC          = .FALSE. 
      beta_H2O          = 1.0
      alpha_H2O         = 0.0415

      idate  = 011027
      rlat   = 50.65_dp
      rlon   = 10.77_dp
      Dust   = 1.0_dp

      Temperature0 = 280.0_dp
      IF (combustion) THEN
        Pressure0 = 200000.0_dp
      ELSE
        Pressure0 = 85000.0_dp
      END IF

      REWIND(RunUnit)
      READ(RunUnit,METEO,IOSTAT=io_stat,IOMSG=io_msg)
      CALL ErrorCheck(io_stat,io_msg,'reading METEO list')

      IF (LWCLevelmin<ZERO .OR. LWCLevelmax<ZERO) STOP 'LWC is given negative in the *.run-file, this is not allowed!'

      ! real comparison: LWCLevelmin == LWCLevelmax ?
      IF ( ABS(LWCLevelmin-LWCLevelmax)<epsL*MAX(LWCLevelmin, LWCLevelmax) ) THEN
        constLWC = .TRUE.
        LWCconst = LWCLevelmin
      END IF


      ! -- Set parcel/box-model conditions
      !
      ! set variant T and p variables to initial value
      T_parcel = Temperature0
      Pressure = Pressure0
      !
      ! initial mixing ratio of parcel
      RH = RH0
      q_parcel = RH * qsatw(T_parcel, Pressure0)
      !
      ! calculate density of air in parcel according to ideal gas law
      rho_parcel  = milli * Pressure0 * molw_air / R / Temperature0 ! kg/m3
      ! moist air
      !rho_parcel  = milli * Pressure0 * molw_air / R / (Temperature0*(1+0.61*q_parcel)) ! kg/m3
      rho_parcel0 = rho_parcel
      ! calculate mass of parcel in kg
      m_parcel    = rho_parcel * V_parcel0

      H2O      = mol2part * RH*esatw(Temperature0)  / SI_Gas / Temperature0   ! passive species H2O [molec/cm3]

      ! assign reference values
      RefH2O   = H2O * Temperature0 / RefTemp * RefPressure / Pressure0 ! H2O conc if parcel had reference temperature and pressure
      RefM     = mol2part * RefPressure / SI_Gas / RefTemp
      RefM_dry = RefM - RefH2O
      RefN2    = 0.7808     * RefM_dry    ! passive species N2 [molec/cm3]
      RefO2    = 0.2095     * RefM_dry    ! passive species O2 [molec/cm3]
      RefH2    = 0.00000055 * RefM_dry    ! passive species H2 [molec/cm3]

      ! assign running values
      N2 = RefN2 * ( RefTemp / Temperature0 * Pressure0 / RefPressure )
      O2 = RefO2 * ( RefTemp / Temperature0 * Pressure0 / RefPressure )
      H2 = RefH2 * ( RefTemp / Temperature0 * Pressure0 / RefPressure )

!
!-----------------------------------------------------------------
!---  Numerics
!-----------------------------------------------------------------
!
!--- Set Default Values for NUMERICS Namelist
      RtolROW     = 1.0e-5_dp   ! Relative tolerance For ROW
      AtolGas     = 1.0e-7_dp   ! Absolute tolerance for gas phase
      AtolAqua    = 1.0e-7_dp   ! Absolute tolerance for liquid phase
      AtolTemp    = 1.0e-3_dp   ! Absolute tolerance for temperature
      AtolWaterMass=1.0e-8_dp   ! Absolute tolerance for water masses
      Atolq       = 1.0e-6_dp   ! Absolute tolerance for mixing ratio
      Atolz       = 1.0e-1_dp   ! Absolute tolerance for height z
      AtolRho     = 1.0e-4_dp   ! Absolute tolerance for parcel density
      Error_Est   = 2           ! error estimation default 2-norm
      minStp      = 1.0e-20_dp  ! minimum timestep of ROW method in [sec]
      maxStp      = 250.0_dp    ! maximum timestep of ROW method in [sec]
      ODEsolver   = 'Rodas3'    ! ROW scheme
      Ordering    = .TRUE.      ! sparse LU, no numerical pivoting

!--- Read NUMERICS namelist
      READ(RunUnit,NUMERICS,IOSTAT=io_stat,IOMSG=io_msg)
      CALL ErrorCheck(io_stat,io_msg,'reading NUMERICS list')

      ! Test that Rosenbrock tolerance > 0
      IF ( RtolROW <= ZERO ) THEN
        WRITE(*,*) '  RtolROW must be a positive scalar!'
        STOP
      END IF

      ! Test that absolute tolerance for gas and aqua species is > 0
      IF ( AtolGas <= ZERO .OR. AtolAqua<=ZERO .OR. AtoLTemp<=ZERO ) THEN
        WRITE(*,*) '  ATols must be positive!'
        STOP
      END IF

      ! Test if maximum stepsize is not to small/big
      IF ( maxStp <= ZERO ) THEN
        maxStp = tEnd - tBegin
      ELSE IF ( maxStp > tEnd-tBegin ) THEN
        WRITE(*,*) '  Maximum stepsize = ',maxStp, ' to high!'
        STOP
      END IF
      IF ( minStp <= 1.e-50_dp ) THEN
        WRITE(*,*) '  Minimums stepsize = ',minStp, ' to low!'
        STOP
      ELSE IF ( minStp > maxStp ) THEN
        WRITE(*,*) '  Minimums stepsize = ', minStp, ' > ', maxStp
        STOP
      END IF


!-----------------------------------------------------------------
!---  Output of Data
!-----------------------------------------------------------------
!
!--- Set Default Values for OUTPUT Namelist
      StpNetcdf     = -1.0_dp      ! Time step for Netcdf output      [in sec]
      StpFlux       = -1.0_dp
      StpConc       = -1.0_dp
      MatrixPrint   = .FALSE.
      DropletClassPrint = .FALSE.
      NetCdfPrint   = .TRUE.
      FluxDataPrint = .FALSE.
      ConcDataPrint = .FALSE.
!
!--- Read OUTPUT namelist
      READ(RunUnit,OUTPUT,IOSTAT=io_stat,IOMSG=io_msg)
      CALL ErrorCheck(io_stat,io_msg,'reading OUTPUT list')

      NetCDFFile = ADJUSTL(NetCDFFile)
      IF ( TRIM(NetCdfFile) == '' ) NetCdfPrint = .FALSE.   ! no output if no filename is declared

   END SUBROUTINE InitRun

  SUBROUTINE ErrorCheck(io_stat,io_msg,cause)
    INTEGER      :: io_stat
    CHARACTER(*) :: io_msg, cause
    IF ( io_stat>0 ) THEN
      WRITE(*,*) '   ERROR while '//cause//'  ::  ',io_stat,'  '//TRIM(io_msg)
      STOP
    END IF
  END SUBROUTINE ErrorCheck

  SUBROUTINE FileNameCheck(Name,miss)
    CHARACTER(*) :: Name
    CHARACTER(*) :: miss
    LOGICAL      :: ex

    INQUIRE(FILE=TRIM(Name), EXIST=ex)
    
    IF ( TRIM(Name) == '' .OR. .NOT.ex ) THEN
      WRITE(*,*); WRITE(*,*)
      WRITE(*,'(10X,A)') 'ERROR    Missing:  '//TRIM(miss)
      WRITE(*,'(10X,A)') '         FileName: '//TRIM(Name)
      WRITE(*,*); WRITE(*,*)
      STOP
    ELSE
      Name = ADJUSTL(Name)
    END IF
  END SUBROUTINE FileNameCheck

END MODULE InitRoutines_Mod
