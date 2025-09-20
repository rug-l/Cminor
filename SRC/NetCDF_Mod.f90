!
! Copyright (C) 2025 Levin Rug, Willi Schimmel (E-Mail: l.rug@lmu.de)
! See ./SRC/Cminor.f90 for the copyright notice
! See ./LICENSE for license information
! SPDX-License-Identifier: GPL-3.0
!
MODULE NetCDF_Mod
!--- Netcdf Output
  USE Kind_Mod,    ONLY: dp
  USE Meteo_Mod,   ONLY: mol2part, Zenith, LWC_array, get_wet_radii
  USE netcdf
  USE Reac_Mod,    ONLY: ns_GAS, Diag_Index, Diag_Name, Diag_LongName, DropletClasses,         &
                       & hasAquaSpc, hasPartiSpc, hasSolidSpc, hasPhotoReac, nDIM, nDIM2,      &
                       & nspc, nFrac, y_emi, y_name, MolMass, iAs, Hp_ind, OHm_ind, ns_AQUA,   &
                       & hasGasSpc
  !
  USE Control_Mod, ONLY: LenName, iNCout_G, iNCout_A_l, iNCout_A_m3, iNCout_A_m3_Ptr,          &
                       & Pi, tBegin, tEnd, combustion, TWO, HOUR, nDropletClasses, NetCdfFile, &
                       & nNcdfAqua, nNcdfEmiss, rlat, rlon, iNcdfAqua, iNcdfAqua2, iNcdfEmiss, &
                       & iNcdfGas, iNcdfGas2, Pi34, rTHREE, UnitGas,       &
                       & nNcdfGas, nD_Ptr_spc, mega, iNCout_A_l_Ptr, &
                       & adiabatic_parcel
  !
  IMPLICIT NONE

  TYPE NetCDF_T
    INTEGER                   :: n_Out         ! number of diagnosis species
    INTEGER                   :: iTime         ! Time step counter
    REAL(dp)                  :: Time          ! Time value
    INTEGER,      ALLOCATABLE :: Spc_Pos(:)    ! Species Positions (indices)
    INTEGER,      ALLOCATABLE :: Emiss_Pos(:)  ! Emission Positions (indices)
    CHARACTER(1), ALLOCATABLE :: Spc_Phase(:)  ! Phase of species g=gas , a=aqua, s=solid, p=parti
    REAL(dp),     ALLOCATABLE :: Spc_Conc(:)   ! Species concentrations at iTime
    REAL(dp),     ALLOCATABLE :: Emiss_Conc(:) ! Species concentrations at iTime
    REAL(dp),     ALLOCATABLE :: AquaSumSpc(:) ! Sum of an aqueous species over all droplet classes at iTime
    REAL(dp),     ALLOCATABLE :: AquaSumDroplet(:) ! Sum of all aqueous species in a droplet classes at iTime
    REAL(dp),     ALLOCATABLE :: DissolvedMass(:)  ! Mass of all dissolved species
    REAL(dp),     ALLOCATABLE :: DissolvedMolecules(:)  ! Sum of all dissolved molecules per droplet
    REAL(dp)                  :: Temperature   ! Temperature value at iTime
    REAL(dp)                  :: rho_air       ! density of air at iTime
    REAL(dp)                  :: q             ! water vapor mixing ratio at iTime
    REAL(dp)                  :: pressure      ! pressure at iTime
    REAL(dp)                  :: z             ! height at iTime
    REAL(dp)                  :: RH            ! relative humidity at iTime
    REAL(dp)                  :: LWC           ! Liquid water content value at iTime
    REAL(dp)                  :: StepSize      ! Step size value at iTime
    REAL(dp)                  :: GasConc       ! Sum of all gaseous species at iTime
    REAL(dp)                  :: SulfConc      ! Sum of all sulphuric species at iTime
    REAL(dp)                  :: Zenith        ! Zenith value at iTime
    REAL(dp),     ALLOCATABLE :: WetRadius(:)  ! Size of cloud dropletts
    INTEGER                   :: MaxErrorSpc   ! Index of species which has max error
    REAL(dp)                  :: ROWerror      ! Error value ROS step at iTime
  END TYPE NetCDF_T

  TYPE(NetCDF_T) :: NetCDF

  INTEGER  :: ncid
  !
  !

  INTEGER, ALLOCATABLE :: Diag_varid(:), WetRadius_varid(:), pH_varid(:), EmissDiag_varid(:)  &
  &                     , AquaSumSpc_varid(:), DissolvedMolecules_varid(:), DissolvedMass_varid(:), AquaSumDroplet_varid(:)
  INTEGER :: rho_parcel_varid, z_parcel_varid, pressure_varid, q_parcel_varid, RH_varid
  INTEGER :: x_varid, y_varid, rec_varid, LWC_varid
  INTEGER :: StepSize_varid
  INTEGER :: Temperature_varid, zenith_varid

  INTEGER, ALLOCATABLE :: pH_ind(:)


  CONTAINS

  !==================================================================
  !===  Initialization of Netcdf Output File
  !==================================================================
  SUBROUTINE InitNetCDF

    ! -- Temporary variables --
    INTEGER :: j, iDr
    !
    !    ============================================================
    !    Variable Attribute Names
    !    ============================================================
    CHARACTER (LEN = *), PARAMETER :: REC_NAME     = "time"
    CHARACTER (LEN = *), PARAMETER :: LON_NAME     = "lon"
    CHARACTER (LEN = *), PARAMETER :: LAT_NAME     = "lat"
    CHARACTER (LEN = *), PARAMETER :: LON_LONGNAME = "longitude"
    CHARACTER (LEN = *), PARAMETER :: LAT_LONGNAME = "latitude"
  !    ============================================================
  !    Variable Attribute Units
  !    ============================================================
    CHARACTER(60),  ALLOCATABLE    :: DIAG_UNITS(:)
    CHARACTER(60)                  :: Emiss_DIAG_UNIT
    CHARACTER(200), ALLOCATABLE    :: Diag_LongName(:)
    CHARACTER (LEN = *), PARAMETER :: NC_UNITS     = "units"
    CHARACTER (LEN = *), PARAMETER :: LON_UNITS    = "degrees_east"
    CHARACTER (LEN = *), PARAMETER :: LAT_UNITS    = "degrees_north"
    !    ============================================================
    !    Dimension ID
    !    ============================================================
    INTEGER, PARAMETER :: NDIMS  = 1
    INTEGER            :: dimIDs(NDIMS)
    INTEGER            :: rec_dimid
    !    ============================================================
    CHARACTER(LenName) :: tmpName
    INTEGER            :: iDiagSpc
    !
    NetCDF%iTime = 0
    NetCDF%n_Out = nNcdfGas + 2*nFrac*nNcdfAqua
    NetCDF%n_Out = nNcdfGas + 2*nDropletClasses*nNcdfAqua

    ! output array containing diagnosis species (and temperature if combustion mechanism)
    ALLOCATE(NetCDF%Spc_Conc(NetCDF%n_Out), NetCDF%WetRadius(nDropletClasses)) 

    ! -- Allocate Netcdf Names etc.
    !
    ALLOCATE( Diag_LongName(NetCDF%n_Out) , Diag_Name(NetCDF%n_Out)   , Diag_Index(NetCDF%n_Out) &
    &       , Diag_UNITS(NetCDF%n_Out)    , Diag_varid(NetCDF%n_Out)    &
    &       , WetRadius_varid(nDropletClasses) , pH_varid(nDropletClasses) , pH_ind(nDropletClasses)  &
    &       , AquaSumDroplet_varid(nDropletClasses), EmissDiag_varid(nNcdfEmiss)                      &
    &       , DissolvedMolecules_varid(nDropletClasses), DissolvedMass_varid(nDropletClasses) )

    Diag_varID      = 0
    WetRadius_varid = 0
    pH_ind          = 0

    j = 0
    ! gaseous species
    ALLOCATE(iNCout_G(0))
    DO iDiagSpc = 1 , nNcdfGas
      j = j + 1
      tmpName = ADJUSTL(y_name(iNcdfGas(iDiagSpc)))
      Diag_Index(j)    = iNcdfGas2(iDiagSpc)
      Diag_Name(j)     = TRIM(tmpName)
      Diag_LongName(j) = TRIM(tmpName)
      IF ( combustion ) THEN
        DIAG_UNITS(j)    = "mol/mol"
      ELSE
        DIAG_UNITS(j)    = "molec/cm3"
      END IF
  
      CALL check_name_slash(Diag_Name(j))
      CALL check_name_bracket(Diag_Name(j))
      iNCout_G = [iNCout_G, j]
    END DO
  
    ! aqueous species
    ALLOCATE(iNCout_A_l(0),iNCout_A_m3(0), iNCout_A_m3_Ptr(nNcdfAqua), iNCout_A_l_Ptr(nNcdfAqua))

    DO iDiagSpc = 1 , nNcdfAqua
      iNCout_A_l_Ptr(iDiagSpc)  = 1 + (iDiagSpc-1)*nDropletClasses
      iNCout_A_m3_Ptr(iDiagSpc) = 1 + (iDiagSpc-1)*nDropletClasses
      DO iDr = 1 , nDropletClasses
        j = j + 1
        WRITE(tmpName,'(A,I0)') TRIM(y_name(iNcdfAqua(iDiagSpc)))//'_',iDr
        ! diag index needs to be the true index (non-compressed)
        Diag_Index(j)    = iNcdfAqua2((iDiagSpc-1)*nDropletClasses+iDr)
        Diag_Name(j)     = TRIM(tmpName)//'_l'
        Diag_LongName(j) = TRIM(tmpName)//'AQUA'
        DIAG_UNITS(j)    = "mol/l"
        iNCout_A_l  = [iNCout_A_l, j]

        ! this has to be here to capture the mol/l index (not mol/m3) 
        IF (INDEX(tmpName,'Hp_')>0) THEN
          pH_ind(iDr) = j
        END IF
  
        j = j + 1
        Diag_Index(j)    = iNcdfAqua2((iDiagSpc-1)*nDropletClasses+iDr)
        Diag_Name(j)     = TRIM(tmpName)//'_m3'
        Diag_LongName(j) = TRIM(tmpName)//'AIR'
        DIAG_UNITS(j)    = "mol/m3"
        iNCout_A_m3 = [iNCout_A_m3, j]
  
        CALL check_name_slash(Diag_Name(j))
        CALL check_name_bracket(Diag_Name(j))
      END DO
    END DO

    Emiss_DIAG_UNIT = "molec/cm3"

    ! ============================================================
    ! --  create Netcdf file (for each size bin / fraction)
    ! ============================================================
    ! Create new NetCDF File
      
    CALL check(  NF90_CREATE(  TRIM(NetCDFFile) &   ! NetCDF output file name
    &                        , NF90_CLOBBER     &   ! Overwrite existing file with the same name
    &                        , ncid     )    )      ! Returned netCDF ID
    !
    ! ============================================================
    ! --  define dimension variables
    ! ============================================================
    ! Define the dimensions. The record dimension is defined to have
    ! unlimited length - it can grow as needed. In this example it is
    ! the time dimension.
    CALL check(  NF90_DEF_DIM(  ncid           & ! NetCDF ID, from prev call NF90_CREATE 
    &                         , REC_NAME          & ! Time is the record dimension
    &                         , NF90_UNLIMITED    & ! space for "unlimited" time steps
    &                         , rec_dimID     )   ) ! returned dimension ID
    !
    ! ============================================================
    ! --  define coordinate variables
    ! ============================================================
    ! Define the coordinate variables. We will only define coordinate
    ! variables for lat and lon.  Ordinarily we would need to provide
    ! an array of dimension IDs for each variable's dimensions, but
    ! since coordinate variables only have one dimension, we can
    ! simply provide the address of that dimension ID (lat_dimid) and
    ! similarly for (lon_dimid).
    !
    dimIDs = (/ rec_dimID /)
    CALL check(  NF90_DEF_VAR(  ncid, REC_NAME,       NF90_DOUBLE, dimIDs, rec_varid  ) )
    CALL check(  NF90_DEF_VAR(  ncid, TRIM(LON_NAME), NF90_DOUBLE, dimIDs, x_varid    ) )
    CALL check(  NF90_DEF_VAR(  ncid, TRIM(LAT_NAME), NF90_DOUBLE, dimIDs, y_varid    ) )
    !
    ! ============================================================
    ! --  define data variables
    ! ============================================================
    DO j = 1 , NetCDF%n_Out
      CALL check( NF90_DEF_VAR( ncid, TRIM(Diag_Name(j)), NF90_DOUBLE, dimIDs, Diag_varid(j)) )
    END DO

    ! lwc
    CALL check ( NF90_DEF_VAR( ncid, 'Step_Size' , NF90_DOUBLE, dimIDS, StepSize_varid ) )
    IF ( hasPhotoReac ) CALL check ( NF90_DEF_VAR( ncid, 'Zenith' , NF90_DOUBLE, dimIDS, zenith_varid ) )
      
    IF ( ns_AQUA>0 ) THEN
      CALL check ( NF90_DEF_VAR( ncid, 'LWC_Level' , NF90_DOUBLE, dimIDS, LWC_varid ) )
      ALLOCATE(NetCDF%DissolvedMolecules(nDropletClasses), NetCDF%DissolvedMass(nDropletClasses), NetCDF%AquaSumDroplet(nDropletClasses))

      DO iDr=1,nDropletClasses
        WRITE(tmpName,'(A,I0)') 'wetRadius_',iDr
        CALL check ( NF90_DEF_VAR( ncid, TRIM(tmpName), NF90_DOUBLE, dimIDS, wetRadius_varid(iDr) ) )
        WRITE(tmpName,'(A,I0)') 'AllAquaSpcSum_',iDr
        CALL check ( NF90_DEF_VAR( ncid, TRIM(tmpName), NF90_DOUBLE, dimIDS, AquaSumDroplet_varid(iDr) ) )
        WRITE(tmpName,'(A,I0)') 'DissolvedMass_',iDr
        CALL check ( NF90_DEF_VAR( ncid, TRIM(tmpName), NF90_DOUBLE, dimIDS, DissolvedMass_varid(iDr) ) )
        WRITE(tmpName,'(A,I0)') 'DissolvedMolecules_',iDr
        CALL check ( NF90_DEF_VAR( ncid, TRIM(tmpName), NF90_DOUBLE, dimIDS, DissolvedMolecules_varid(iDr) ) )
          
        IF (pH_ind(iDr)/=0) THEN
          WRITE(tmpName,'(A,I0)') 'pH_Value_',iDr
          CALL check ( NF90_DEF_VAR( ncid, TRIM(tmpName), NF90_DOUBLE, dimIDS, pH_varid(iDr) ) )
        END IF
      END DO
      IF (nDropletClasses>1) THEN
        ALLOCATE(AquaSumSpc_varid(nNcdfAqua))
        ALLOCATE(NetCDF%AquaSumSpc(nNcdfAqua))
        DO iDiagSpc = 1,nNcdfAqua
          WRITE(tmpName,'(A,I0)') TRIM(y_name(iNcdfAqua(iDiagSpc)))//'_Sum'
          CALL check ( NF90_DEF_VAR( ncid, TRIM(tmpName), NF90_DOUBLE, dimIDS, AquaSumSpc_varid(iDiagSpc) ) )
        END DO
      END IF
    END IF

    CALL check ( NF90_DEF_VAR( ncid, 'Temperature' , NF90_DOUBLE, dimIDS, Temperature_varid ) )

    DO j = 1 , nNcdfEmiss
      CALL check( NF90_DEF_VAR( ncid, 'Emiss_'//TRIM(y_name(iNcdfEmiss(j))), NF90_DOUBLE, dimIDs, EmissDiag_varid(j)) )
    END DO

    IF ( adiabatic_parcel ) THEN
      CALL check( NF90_DEF_VAR( ncid, 'rho_parcel', NF90_DOUBLE, dimIDs, rho_parcel_varid) )
      CALL check( NF90_DEF_VAR( ncid, 'q_parcel',   NF90_DOUBLE, dimIDs, q_parcel_varid) )
      CALL check( NF90_DEF_VAR( ncid, 'pressure',   NF90_DOUBLE, dimIDs, pressure_varid) )
      CALL check( NF90_DEF_VAR( ncid, 'z_parcel',   NF90_DOUBLE, dimIDs, z_parcel_varid) )
      CALL check( NF90_DEF_VAR( ncid, 'RH_parcel',  NF90_DOUBLE, dimIDs, RH_varid) )
    END IF

    !
    ! ============================================================
    ! --  define attributes (name, unit...) of the variables
    ! ============================================================
    ! x = longitude
    CALL check( NF90_PUT_ATT( ncid, x_varid, "standard_name", TRIM(LON_LONGNAME) ) )
    CALL check( NF90_PUT_ATT( ncid, x_varid, "long_name",     TRIM(LON_LONGNAME) ) )
    CALL check( NF90_PUT_ATT( ncid, x_varid, NC_UNITS,        TRIM(LON_UNITS)    ) )
    CALL check( NF90_PUT_ATT( ncid, x_varid, "axis",          "x"                ) )
    CALL check( NF90_PUT_ATT( ncid, x_varid, "_CoordinateAxisType", "lon"        ) )
    ! y = latitude
    CALL check( NF90_PUT_ATT( ncid, y_varid, "standard_name", TRIM(LAT_LONGNAME) ) )
    CALL check( NF90_PUT_ATT( ncid, y_varid, "long_name",     TRIM(LAT_LONGNAME) ) )
    CALL check( NF90_PUT_ATT( ncid, y_varid, NC_UNITS,        TRIM(LAT_UNITS)    ) )
    CALL check( NF90_PUT_ATT( ncid, y_varid, "axis",          "y"                ) )
    CALL check( NF90_PUT_ATT( ncid, y_varid, "_CoordinateAxisType", "lat"        ) )
    ! rec = time
    CALL check( NF90_PUT_ATT( ncid, rec_varid, "standard_name",       "time"    ) )
    CALL check( NF90_PUT_ATT( ncid, rec_varid, "long_name",           REC_NAME  ) )
    IF (combustion) THEN
      CALL check( NF90_PUT_ATT( ncid, rec_varid, NC_UNITS, "seconds" ) )
    ELSE
      CALL check( NF90_PUT_ATT( ncid, rec_varid, NC_UNITS, "hours" ) )
    END IF
    CALL check( NF90_PUT_ATT( ncid, rec_varid, "_CoordinateAxisType", "time"    ) )
      
    ! Diagnose species
    DO j=1,NetCDF%n_Out
      CALL check( NF90_PUT_ATT(ncid, Diag_varid(j), NC_UNITS, Diag_UNITS(j)) )
      CALL check( NF90_PUT_ATT(ncid, Diag_varid(j), "long_name", TRIM(Diag_LongName(j))) )
      CALL check( NF90_PUT_ATT(ncid, Diag_varid(j), "_CoordinateAxes", "time") )
    END DO

    ! Diagnose Emissions
    DO j=1,nNcdfEmiss
      CALL check( NF90_PUT_ATT(ncid, EmissDiag_varid(j), NC_UNITS, Emiss_Diag_UNIT) )
      CALL check( NF90_PUT_ATT(ncid, EmissDiag_varid(j), "long_name", 'Emiss_'//TRIM(y_name(iNcdfEmiss(j)))) )
      CALL check( NF90_PUT_ATT(ncid, EmissDiag_varid(j), "_CoordinateAxes", "time") )
    END DO

    IF ( adiabatic_parcel ) THEN
      CALL check( NF90_PUT_ATT(ncid, rho_parcel_varid, NC_UNITS, '[kg/m3]' ) )  
      CALL check( NF90_PUT_ATT(ncid, rho_parcel_varid, "long_name", 'density of adiabatic parcel') ) 
      CALL check( NF90_PUT_ATT(ncid, rho_parcel_varid, "_CoordinateAxes", "time") )
      CALL check( NF90_PUT_ATT(ncid, q_parcel_varid  , NC_UNITS, '[kg/kg]' ) )  
      CALL check( NF90_PUT_ATT(ncid, q_parcel_varid  , "long_name", 'mixing ratio of adiabatic parcel') ) 
      CALL check( NF90_PUT_ATT(ncid, q_parcel_varid  , "_CoordinateAxes", "time") )
      CALL check( NF90_PUT_ATT(ncid, pressure_varid  , NC_UNITS, '[Pa]' ) )  
      CALL check( NF90_PUT_ATT(ncid, pressure_varid  , "long_name", 'pressure of adiabatic parcel') ) 
      CALL check( NF90_PUT_ATT(ncid, pressure_varid  , "_CoordinateAxes", "time") )
      CALL check( NF90_PUT_ATT(ncid, z_parcel_varid  , NC_UNITS, '[m]' ) )  
      CALL check( NF90_PUT_ATT(ncid, z_parcel_varid  , "long_name", 'height of adiabatic parcel') ) 
      CALL check( NF90_PUT_ATT(ncid, z_parcel_varid  , "_CoordinateAxes", "time") )
      CALL check( NF90_PUT_ATT(ncid, RH_varid        , NC_UNITS, '[none]' ) )  
      CALL check( NF90_PUT_ATT(ncid, RH_varid        , "long_name", 'relative humidity in adiabatic parcel') ) 
      CALL check( NF90_PUT_ATT(ncid, RH_varid        , "_CoordinateAxes", "time") )
    END IF

    ! stepsize
    CALL check( NF90_PUT_ATT(ncid, StepSize_varid, NC_UNITS, '[sec]' ) )  
    CALL check( NF90_PUT_ATT(ncid, StepSize_varid, "long_name", 'step size in [sec]') ) 
    CALL check( NF90_PUT_ATT(ncid, StepSize_varid, "_CoordinateAxes", "time") )
    ! pH value
    IF ( ns_AQUA>0 ) THEN
      ! lwc
      CALL check( NF90_PUT_ATT(ncid, LWC_varid, NC_UNITS, '[l/m3]' ) )  
      CALL check( NF90_PUT_ATT(ncid, LWC_varid, "long_name", '[liter/m3]') ) 
      CALL check( NF90_PUT_ATT(ncid, LWC_varid, "_CoordinateAxes", "time") )
      DO iDr=1,nDropletClasses
        CALL check( NF90_PUT_ATT(ncid, wetRadius_varid(iDr), NC_UNITS, '[m]' ) )  
        CALL check( NF90_PUT_ATT(ncid, wetRadius_varid(iDr), "long_name", 'wet droplet radius [m]') ) 
        CALL check( NF90_PUT_ATT(ncid, wetRadius_varid(iDr), "_CoordinateAxes", "time") )
        ! pH value
        IF (pH_ind(iDr)/=0) THEN
          CALL check( NF90_PUT_ATT(ncid, pH_varid(iDr), NC_UNITS, '[-]' ) )  
          CALL check( NF90_PUT_ATT(ncid, pH_varid(iDr), "long_name", 'pH-Value ( = -log10[Hp] )') ) 
          CALL check( NF90_PUT_ATT(ncid, pH_varid(iDr), "_CoordinateAxes", "time") )
        END IF
        ! sum of all species in the droplet 
        CALL check( NF90_PUT_ATT(ncid, AquaSumDroplet_varid(iDr), NC_UNITS, '[mol/l]' ) )  
        CALL check( NF90_PUT_ATT(ncid, AquaSumDroplet_varid(iDr), "Long_name", 'sum aqueus conc [mol/l]') ) 
        CALL check( NF90_PUT_ATT(ncid, AquaSumDroplet_varid(iDr), "_CoordinateAxes", "time") )
        ! mass of all species in the droplet 
        CALL check( NF90_PUT_ATT(ncid, DissolvedMass_varid(iDr), NC_UNITS, '[g]' ) )  
        CALL check( NF90_PUT_ATT(ncid, DissolvedMass_varid(iDr), "Long_name", 'mass of all dissolved species [g]') ) 
        CALL check( NF90_PUT_ATT(ncid, DissolvedMass_varid(iDr), "_CoordinateAxes", "time") )
        ! sum of molecules of all species in the droplet 
        CALL check( NF90_PUT_ATT(ncid, DissolvedMolecules_varid(iDr), NC_UNITS, '[molec/droplet]' ) )  
        CALL check( NF90_PUT_ATT(ncid, DissolvedMolecules_varid(iDr), "Long_name", 'sum of all dissolved molecules [molecules/droplet]') ) 
        CALL check( NF90_PUT_ATT(ncid, DissolvedMolecules_varid(iDr), "_CoordinateAxes", "time") )
      END DO
      IF (nDropletClasses>1) THEN
        DO iDiagSpc = 1, nNcdfAqua
          CALL check( NF90_PUT_ATT(ncid, AquaSumSpc_varid(iDiagSpc), NC_UNITS, "mol/m3") )
          CALL check( NF90_PUT_ATT(ncid, AquaSumSpc_varid(iDiagSpc), "long_name", "sum of the species in all droplet classes") )
          CALL check( NF90_PUT_ATT(ncid, AquaSumSpc_varid(iDiagSpc), "_CoordinateAxes", "time") )
        END DO
      END IF
    END IF

    ! save temperature
    CALL check( NF90_PUT_ATT(ncid, Temperature_varid, NC_UNITS, '[K]' ) )  
    CALL check( NF90_PUT_ATT(ncid, Temperature_varid, "long_name", 'Temperature in Kelvin') ) 
    CALL check( NF90_PUT_ATT(ncid, Temperature_varid, "_CoordinateAxes", "temp") )

    ! ============================================================
    ! --  End define mode
    ! ============================================================
    !
    CALL check( NF90_ENDDEF( ncid ) )
    !
    ! ============================================================
    ! --  Close netcdf file
    ! ============================================================
    !
    CALL check( NF90_CLOSE( ncid ) )

    WRITE(*,'(10X,A)') 'Init NetCDF................... done'
    !      
    !
    ! ============================================================
    CONTAINS
    !  
    ! Check if species name contains slashes, if -> replace all '/' with '_'
    SUBROUTINE check_name_slash(Name)
      CHARACTER(100) :: Name
      INTEGER        :: pos

      DO
        pos = INDEX(TRIM(ADJUSTL(Name)),'/')
        IF ( pos > 0 ) THEN
          Name = TRIM(Name(1:pos-1)//'_'//Name(pos+1:))
        ELSE
          EXIT
        END IF
      END DO
    END SUBROUTINE check_name_slash

    ! Check if the first character is a bracket '[', if -> replace '[' with '_['
    SUBROUTINE check_name_bracket(Name)
      CHARACTER(100) :: Name
      IF ( (INDEX(TRIM(ADJUSTL(Name)),'['))==1) Name = '_'//TRIM(ADJUSTL(Name))
    END SUBROUTINE check_name_bracket
    
    
    SUBROUTINE check(STATUS)
      INTEGER, INTENT ( in) :: STATUS
      ! 
      IF(STATUS /= nf90_noerr) THEN
        WRITE(*,*) '  Error with NetCDF data ', trim(nf90_strerror(STATUS)), STATUS
        STOP 
      END IF
    END SUBROUTINE check
  
END SUBROUTINE InitNetCDF
!
!
SUBROUTINE SetOutputNCDF(NCDF,Time,StpSize,Conc,Temp,rho_air,q,pressure,z,RH,nDroplets,LWCs)
  ! OUT:
  TYPE(NetCDF_T) :: NCDF

  REAL(dp), INTENT(IN) :: Conc(:)
  REAL(dp), INTENT(IN) :: Time, StpSize, Temp, rho_air, q, pressure, RH, z
  REAL(dp), DIMENSION(nDropletClasses) :: LWCs, nDroplets

  !-- internal variable
  INTEGER :: j
  REAL(dp), PARAMETER    :: y_Min = 1.e-40_dp    
  REAL(dp), ALLOCATABLE :: tConc(:), tG(:), tA(:)

  !==================================================================
  !===  Saving Output
  !==================================================================

  DO j=1,SIZE(Conc)
    IF ( ISNAN(Conc(j)) ) THEN
      WRITE(*,*); WRITE(*,*); WRITE(*,*)
      WRITE(*,'(10X,A)')        '  ERROR:  Species concentration is NaN ! '
      WRITE(*,'(10X,A)')        '  -------------------------------------- '
      WRITE(*,'(10X,A,I0)')     '  NCDF idx      =  ', j
      WRITE(*,'(10X,A,A)')      '  Species name  =  ', y_name(j)
      WRITE(*,'(10X,A,Es12.4)') '  Species val   =  ', Conc(j)
      WRITE(*,*); WRITE(*,*); WRITE(*,*)
      STOP
    END IF
  END DO

  tConc = [ MAX(Conc,y_Min) ]             ! minimum for logarithmic plot
  IF ( combustion ) THEN
    tConc(1:nspc) = MoleConc_to_MoleFr(tConc(1:nspc))
    tG = [ tConc(iNcdfGas) ]     ! in mole fractions [-]
  ELSE
    tG = [ tConc(iNcdfGas2)   ] !  molec/cm3
    tA = [ tConc(iNcdfAqua2) / mol2part ]
  END IF

  NCDF%Temperature = Temp

  IF (adiabatic_parcel) THEN
    NCDF%rho_air     = rho_air
    NCDF%q           = q
    NCDF%pressure    = pressure
    NCDF%z           = z
    NCDF%RH          = RH
  END IF

  NCDF%Time        = Time
  NCDF%StepSize    = StpSize

  ! NCDF%Spc_Conc = 
  ! ------------------------------------------------------ !
  ! |   GAS   |   AQ_1  |   AQ_2   |   ....   |   AQ_n   | !   
  ! ------------------------------------------------------ !
   
  IF (hasGasSpc) NCDF%Spc_Conc(iNCout_G) = tG

  IF (hasAquaSpc) THEN
    NCDF%LWC = SUM(LWCs)

    NCDF%WetRadius = get_wet_radii(LWCs=LWCs, given_rho=rho_air)

    ! for mol/l unit, divide all species by dropletclass-specific LWC
    DO j = 1 , nDropletClasses
      NCDF%Spc_Conc(iNCout_A_l(iNCout_A_l_Ptr+j-1))  = tA(iNCout_A_l_Ptr+j-1)/LWCs(j)
      NCDF%AquaSumDroplet(j) = SUM(NCDF%Spc_Conc(iNCout_A_l(iNCout_A_l_Ptr+j-1))) - (tConc(nD_Ptr_spc(Hp_ind)+j-1)+tConc(nD_Ptr_spc(OHm_ind)+j-1))/mol2part/LWCs(j)
      NCDF%DissolvedMass(j)  = ( SUM(MolMass(iAs)*tConc(nD_Ptr_spc(iAs) +j-1) / mol2part)   &
                             &   - 1.0  * tConc(nD_Ptr_spc(Hp_ind) +j-1) / mol2part    &
                             &   - 17.0 * tConc(nD_Ptr_spc(OHm_ind)+j-1) / mol2part )  &
                             &  / (nDroplets(j)*rho_air)
      NCDF%DissolvedMolecules(j) = (SUM(tConc(nD_Ptr_spc(iAs)+j-1)) - (tConc(nD_Ptr_spc(Hp_ind)+j-1)+tConc(nD_Ptr_spc(OHm_ind)+j-1))) &
                               & * mega / (nDroplets(j)*rho_air)
    END DO
    NCDF%Spc_Conc(iNCout_A_m3) = tA
    IF (nDropletClasses>1) THEN
      DO j = 1 , nNcdfAqua
        NCDF%AquaSumSpc(j) = SUM(tA(1+(j-1)*nDropletClasses:j*nDropletClasses))
      END DO
    END IF
  END IF

  IF ( hasPhotoReac ) NCDF%Zenith = Zenith(Time)

  ! this is, unfortunately, not linearly interpolated, since y_emi is updated
  ! while integrating, and not before and after, which would be neccessary to interpolate
  IF (nNcdfEmiss>0) NCDF%Emiss_Conc = y_emi(iNcdfEmiss)

  CONTAINS
  
  FUNCTION MoleConc_to_MoleFr(MoleConc) RESULT(MoleFr)
    REAL(dp), ALLOCATABLE :: MoleFr(:)   ! Mole fraction [mol/mol]
    REAL(dp), INTENT(IN)  :: MoleConc(:) ! Mole concentration  [mol/cm3]

    MoleFr = MoleConc / SUM(MoleConc)

  END FUNCTION MoleConc_to_MoleFr

END SUBROUTINE SetOutputNCDF
  
  
  
  SUBROUTINE StepNetcdf( NCDF )
  !
  !==================================================================
  !===  Initialization of ASCII Output File
  !==================================================================
  !  
  ! IN:
  TYPE(NetCDF_T) :: NCDF

  REAL(dp) :: pH(nDropletClasses)
  
  REAL(dp) :: timeLoc
  !
  !-- internal variable
  INTEGER :: j, iDr

  NetCDF%iTime = NetCDF%iTime + 1
  !
  ! ====================================================================================
  ! == Open Netcdf-File in write mode to modify data ===================================
  ! ====================================================================================

  CALL check( NF90_OPEN( TRIM(NetcdfFile), NF90_WRITE, ncid ) ) 

  ! ====================================================================================
  ! == Write new data to netcdf-file ===================================================
  ! ====================================================================================
  IF ( combustion ) THEN
    timeLoc = NCDF%Time       ! time in seconds for combustion mechanism
  ELSE
    timeLoc = NCDF%Time/HOUR  ! time output in hours
  END IF

  CALL check( NF90_PUT_VAR( ncid, rec_varid, timeLoc, start=(/NCDF%iTime/) ) )
  CALL check( NF90_PUT_VAR( ncid, x_varid,   rlon,    start=(/NCDF%iTime/) ) )
  CALL check( NF90_PUT_VAR( ncid, y_varid,   rlat,    start=(/NCDF%iTime/) ) )
    
  !----------------------------------------------------------------
  ! Write concentration
  DO j=1,NCDF%n_Out
    CALL check( NF90_PUT_VAR( ncid, Diag_varid(j), NCDF%Spc_Conc(j), start = (/NCDF%iTime/) ) )
  END DO
    
  !----------------------------------------------------------------
  ! Write emissions
  DO j=1,nNcdfEmiss
    CALL check( NF90_PUT_VAR( ncid, EmissDiag_varid(j), NCDF%Emiss_Conc(j), start = (/NCDF%iTime/) ) )
  END DO

  ! Write stepsize
  CALL check( NF90_PUT_VAR( ncid, StepSize_varid, NCDF%StepSize, start = (/NCDF%iTime/) ) )
    
  ! write zenith, error value stepssize control, index of max. error species
  IF ( hasPhotoReac ) CALL check( NF90_PUT_VAR( ncid, zenith_varid, NCDF%Zenith, start = (/NCDF%iTime/) ) )

  ! write pH-values, droplet radii, liquid water content
  IF ( ns_AQUA>0 ) THEN
    DO iDr=1,nDropletClasses
      IF (pH_ind(iDr)>=LBOUND(NCDF%Spc_Conc,1) .AND. pH_ind(iDr)<=UBOUND(NCDF%Spc_Conc,1)) THEN
        pH(iDr) = -LOG10( NCDF%Spc_Conc(pH_ind(iDr)) )
      CALL check( NF90_PUT_VAR( ncid, pH_varid(iDr), pH(iDr), start = (/NCDF%iTime/) ) )
      END IF
      CALL check( NF90_PUT_VAR( ncid, wetRadius_varid(iDr), NCDF%WetRadius(iDr), start = (/NCDF%iTime/) ) )
      CALL check( NF90_PUT_VAR( ncid, AquaSumDroplet_varid(iDr), NCDF%AquaSumDroplet(iDr), start = (/NCDF%iTime/) ) )
      CALL check( NF90_PUT_VAR( ncid, DissolvedMass_varid(iDr), NCDF%DissolvedMass(iDr), start = (/NCDF%iTime/) ) )
      CALL check( NF90_PUT_VAR( ncid, DissolvedMolecules_varid(iDr), NCDF%DissolvedMolecules(iDr), start = (/NCDF%iTime/) ) )
    END DO
    CALL check( NF90_PUT_VAR( ncid, LWC_varid, NCDF%LWC, start = (/NCDF%iTime/) ) )
    IF (nDropletClasses>1) THEN
      DO j=1,nNcdfAqua
        CALL check( NF90_PUT_VAR( ncid, AquaSumSpc_varid(j), NCDF%AquaSumSpc(j), start = (/NCDF%iTime/) ) )
      END DO
    END IF
  END IF

  IF ( adiabatic_parcel ) THEN
    CALL check( NF90_PUT_VAR( ncid, rho_parcel_varid, NCDF%rho_air , start = (/NCDF%iTime/) ) )
    CALL check( NF90_PUT_VAR( ncid, q_parcel_varid  , NCDF%q       , start = (/NCDF%iTime/) ) )
    CALL check( NF90_PUT_VAR( ncid, pressure_varid  , NCDF%pressure, start = (/NCDF%iTime/) ) )
    CALL check( NF90_PUT_VAR( ncid, z_parcel_varid  , NCDF%z       , start = (/NCDF%iTime/) ) )
    CALL check( NF90_PUT_VAR( ncid, RH_varid        , NCDF%RH      , start = (/NCDF%iTime/) ) )
  END IF

  ! Write temperatue value
  CALL check( NF90_PUT_VAR( ncid, Temperature_varid, NCDF%Temperature, start = (/NCDF%iTime/) ) )

  ! ====================================================================================
  ! == Close netcdf-file ===============================================================
  ! ====================================================================================
  !
  CALL check( nf90_sync(ncid) )
  CALL check( nf90_close(ncid) )

  !
  CONTAINS
    SUBROUTINE check(STATUS)
      INTEGER, INTENT ( in) :: STATUS
      !
      IF(STATUS /= nf90_noerr) THEN
        PRINT *, '  Error StepNetCDF :: ' , trim(nf90_strerror(STATUS))
        STOP "Stopped at StepNetCDF"
      END IF
    END SUBROUTINE check
    !
  END SUBROUTINE StepNetcdf

END MODULE NetCDF_Mod

