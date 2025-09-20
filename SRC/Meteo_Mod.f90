!
! Copyright (C) 2025 Levin Rug, Willi Schimmel (E-Mail: l.rug@lmu.de)
! See ./SRC/Cminor.f90 for the copyright notice
! See ./LICENSE for license information
! SPDX-License-Identifier: GPL-3.0
!
MODULE Meteo_Mod

  USE Kind_Mod,    ONLY: dp
  USE Control_Mod, ONLY: hour, hourday, TWO, rTWO, THREE, rlat, rlon, FOUR, rFOUR, iDate,   &
                       & FIVE, dr, ZERO, Pi, ONE, LWCConst, micro, milli, adiabatic_parcel, &
                       & nDropletClasses, Pressure0, Temperature0, joule_to_cal, kilo,      &
                       & joule_to_erg, atm_to_Pa, rho_parcel, rho_parcel0, rTHREE, Pi34,    &
                       & mega, Pi43, q_parcel, atm_to_Pa, V_parcel0

  IMPLICIT NONE

  REAL(dp), PARAMETER :: Cp   = 1004.0D0 ! heat capacity of air at constant pressure J/kg/K
  REAL(dp), PARAMETER :: Cp_vapor = 1897.0D0 ! heat capacity of water vapor at constant pressure J/kg/K
  REAL(dp), PARAMETER :: Cv   = 717.0D0
  REAL(dp), PARAMETER :: Rd   = Cp-Cv
  REAL(dp), PARAMETER :: g    = 9.81 ! gravitational acceleration m/s2
  REAL(dp), PARAMETER :: L_v   = 2.5104e+6 ! latent heat of condensation J/kg

  REAL(dp)            :: N2                               &   ! passive species N2
  &                    , O2                               &   ! passive species O2
  &                    , H2O                              &   ! passive species H2O
  &                    , H2                               &   ! passive species H2
  &                    , N2O2                                 ! passive species N2 + O2

  REAL(dp), PARAMETER :: GasConst_R = 0.0820574d0           &   ! [in l*atm/mol/K]
  &                    , RefTemp    = 298.15D0              &
  &                    , RefRH      = 0.6                   &
  &                    , InvRefTemp = 1.0D0/RefTemp         &
  &                    , aH2OmolperL= 5.55D01               &
  &                    , rho_H2O    = 1000000               & ! density of water in g/m3
  &                    , molw_H2O   = 18.015                & ! molecular weight of water in g/mol
  &                    , molw_air   = 28.9652d0             & ! molar mass of dry air in g/mol
  &                    , h2o_push   = 6.023d20/molw_h2o     &
  &                    , Rv         = 461.5                 & ! gas constant for water in J/kg/K
  &                    , R_air      = 287.1                 & ! gas constant for dry air in J/kg/K
  &                    , RefPressure= 101325.0              & ! reference atmospheric pressure [Pa]
  &                    , SI_am      = 1.66053892173d-27     & ! Atomic mass unit  [kg]  
  &                    , SI_na      = 6.0221412927d+23      & ! Avogadro's number [1/mol]
  &                    , SI_kB      = 1.380648813d-23       & ! Bolzmann constant [J/K]
  &                    , SI_Gas     = SI_na * SI_kB         & ! Gas constant      [J/mol/K]
  &                    , mol2part   = SI_na*micro           &
  &                    , aH2O       = aH2OmolperL*mol2part

  ! passive species reference values are initialized in InitRoutines_Mod according to user specified values
  REAL(dp)            :: RefM_dry  & ! reference number of molecules per cm3 of dry air
  &                    , RefM      & ! reference number of molecules per cm3 of moist air
  &                    , RefH2O    & ! passive species H2O [molec/cm3] at reference pressure and temperature
  &                    , RefN2     & ! passive species N2  [molec/cm3] at reference pressure and temperature
  &                    , RefO2     & ! passive species O2  [molec/cm3] at reference pressure and temperature
  &                    , RefH2     & ! passive species H2  [molec/cm3] at reference pressure and temperature
  &                    , alpha_H2O & ! accommodation coefficient for water condensation
  &                    , beta_H2O    ! condensation coefficient for water condensation

!     Universal gas constant [J / mol / K]
  REAL(dp), PARAMETER :: R         = SI_Gas
  REAL(dp), PARAMETER :: rR        = ONE/R
!
!     Universal gas constant, calorie units [cal / mol / K]
  REAL(dp), PARAMETER :: Rcal      = R * joule_to_cal
  REAL(dp), PARAMETER :: rRcal     = ONE/Rcal
!
!     Universal gas constant, CGS units [erg / mol / K]
  REAL(dp), PARAMETER :: Rerg      = R * joule_to_erg
  REAL(dp), PARAMETER :: rRerg      = ONE/Rerg
!
!     Standard pressure [Pa]
  REAL(dp), PARAMETER :: Patm      = atm_to_Pa
  REAL(dp), PARAMETER :: rPatm     = ONE/Patm

  ! threshold for suppressing reactions in droplets with sum of concentrations higher than that
  REAL(dp), PARAMETER :: eps_ionic_strength = 0.02 ! [mol/l]

!-- more LWC stuff for pseudo function 
  REAL(dp) :: LWCb(6) ! boundaries for linear pseudo lwc function

CONTAINS 

  !=========================================================================!
  !                  calculate sun 
  !=========================================================================!
  FUNCTION Zenith(Time) RESULT(Chi)
    !-----------------------------------------------------------------------!
    ! Input:
    !   - Time
    REAL(dp) :: Time
    !-----------------------------------------------------------------------!
    ! Output:
    !   - sun angle chi
    REAL(dp) :: Chi
    !-----------------------------------------------------------------------!
    ! Temporary variables:
    !INTEGER :: IDAT
    REAL(dp) :: LBGMT, LZGMT
    REAL(dp) :: ML
    ! 
    REAL(dp) :: GMT
    REAL(dp) :: RLT, RPHI
    !    
    INTEGER  :: IIYEAR, IYEAR, IMTH, IDAY, IIY, NYEARS, LEAP, NOLEAP
    REAL(dp) :: YREF,YR
    !   
    INTEGER  :: I, IJ, JD, IJD, IN
    REAL(dp) :: D, RML, W, WR, EC, EPSI, PEPSI, YT, CW, SW, SSW  & 
    &         , EYT, FEQT1, FEQT2, FEQT3, FEQT4, FEQT5, FEQT6 &
    &         , FEQT7, FEQT, EQT
    !         
    REAL(dp) :: REQT, RA, RRA, TAB, RDECL, DECL, ZPT, CSZ, ZR    &
    &         , CAZ, RAZ, AZIMUTH
    !           
    INTEGER :: IMN(12)
    DATA IMN/31,28,31,30,31,30,31,31,30,31,30,31/
    !
    !----------------------------------------------------------------------!
    !
    ! set GMT
    GMT = Time / HOUR
    !
    !  convert to radians
    RLT = rlat*DR
    RPHI = rlon*DR
    !
    !  parse date
    IIYEAR = iDate/10000
    IYEAR = 19*100 + IIYEAR
    IF (IIYEAR <= 50) IYEAR = IYEAR + 100 
    IMTH = (iDate - IIYEAR*10000)/100
    IDAY = iDate - IIYEAR*10000 - IMTH*100
    !
    !  identify and correct leap years
    IIY = (IIYEAR/4)*4
    IF(IIY.EQ.IIYEAR) IMN(2) = 29
    !
    !  count days from Dec.31,1973 to Jan 1, YEAR, then add to 2,442,047.5
    YREF =  2442047.5_dp
    NYEARS = IYEAR - 1974
    LEAP = (NYEARS+1)/4
    IF(NYEARS.LE.-1) LEAP = (NYEARS-2)/4
    NOLEAP = NYEARS - LEAP
    YR = YREF + 365.0_dp*NOLEAP + 366.0_dp*LEAP
    !
    IJD = 0
    IN = IMTH - 1
    IF(IN.EQ.0) GO TO 40
    DO 30 I=1,IN
    IJD = IJD + IMN(I)
  30   CONTINUE
    IJD = IJD + IDAY
    GO TO 50
  40   IJD = IDAY
  50   IJ = IYEAR - 1973
    !
    !      print julian days current "ijd"
    JD = IJD + INT(YR - YREF)
    D = JD + GMT/24.0_dp
    !
    !      calc geom mean longitude
    ML = 279.2801988_dp + .9856473354_dp*D + 2.267e-13_dp*D*D
    RML = ML*DR
    !
    !      calc equation of time in sec
    !      w = mean long of perigee
    !      e = eccentricity
    !      epsi = mean obliquity of ecliptic
    W = 282.4932328_dp + 4.70684e-5_dp*D + 3.39e-13_dp*D*D
    WR = W*DR
    EC = 1.6720041e-2_dp - 1.1444e-9_dp*D - 9.4e-17_dp*D*D
    EPSI = 23.44266511_dp - 3.5626e-7_dp*D - 1.23e-15_dp*D*D
    PEPSI = EPSI*DR
    YT = (TAN(PEPSI*rTWO))**2
    CW = COS(WR)
    SW = SIN(WR)
    SSW = SIN(TWO*WR)
    EYT = TWO*EC*YT
    FEQT1 = SIN(RML)*(-EYT*CW - TWO*EC*CW)
    FEQT2 = COS(RML)*(TWO*EC*SW - EYT*SW)
    FEQT3 = SIN(TWO*RML)*(YT - (FIVE*EC*EC*rFOUR)*(CW*CW-SW*SW))
    FEQT4 = COS(TWO*RML)*(FIVE*EC**2*SSW*rFOUR)
    FEQT5 = SIN(THREE*RML)*(EYT*CW)
    FEQT6 = COS(THREE*RML)*(-EYT*SW)
    FEQT7 = -SIN(FOUR*RML)*(rTWO*YT*YT)
    FEQT = FEQT1 + FEQT2 + FEQT3 + FEQT4 + FEQT5 + FEQT6 + FEQT7
    EQT = FEQT*13751.0_dp
    !
    !   convert eq of time from sec to deg
    REQT = EQT/240.0_dp
    !
    !   calc right ascension in rads
    RA = ML - REQT
    RRA = RA*DR
    !
    !   calc declination in rads, deg
    TAB = 0.43360_dp*SIN(RRA)
    RDECL = ATAN(TAB)
    DECL = RDECL/DR
    !
    !   calc local hour angle
    LBGMT = 12.0_dp - EQT/HOUR + rlon*24.0_dp/360.0_dp
    LZGMT = 15.0_dp*(GMT - LBGMT)
    ZPT = LZGMT*DR
    CSZ = SIN(RLT)*SIN(RDECL) + COS(RLT)*COS(RDECL)*COS(ZPT)
    ZR = ACOS(CSZ)
    ! 
    !   calc local solar azimuth
    CAZ = (SIN(RDECL) - SIN(RLT)*COS(ZR))/(COS(RLT)*SIN(ZR))
    RAZ = ACOS(CAZ)
    AZIMUTH = RAZ/DR
    !
    !--- set Zenith Angle
    Chi =  1.745329252e-02_dp * ZR/DR
  END FUNCTION Zenith

  FUNCTION UpdateSun(Time) RESULT(Sun)
  !--------------------------------------------------------------------
    ! Input:
    !   - Time
    REAL(dp) :: Time
    !--------------------------------------------------------------------!
    ! Output:
    !   - Sun
    REAL(dp)  :: Sun
    !--------------------------------------------------------------------!
    ! Temporary variables:
    REAL(dp), PARAMETER :: SunRise=4.50_dp, SunSet=19.50_dp
    REAL(dp) :: Thour, Tlocal, Ttmp
    !
    Thour  = MODULO(Time / HOUR,24.0_dp)
    Tlocal = Thour - FLOOR(Thour/hourday)*hourday
    !
    IF( (Tlocal>=SunRise) .AND. (Tlocal<=SunSet) ) THEN
      Ttmp = (TWO*Tlocal-SunRise-SunSet) / (SunSet-SunRise);
      IF ( Ttmp >ZERO ) THEN
        Ttmp =  Ttmp * Ttmp
      ELSE
        Ttmp = -Ttmp * Ttmp
      END IF
      Sun = (ONE+COS(Pi*Ttmp)) * rTWO
    ELSE
      Sun = ZERO
    END IF
  END FUNCTION UpdateSun
    
  FUNCTION Set_pseudoLWCbounds() RESULT(bounds)
    USE Control_Mod, ONLY: tBegin, HOUR
    REAL(dp) :: bounds(6)
    ! --- set cloud intervall
    bounds(1) = tBegin * HOUR 
    bounds(2) = bounds(1) + 1.00_dp * HOUR
    bounds(3) = bounds(2) + 0.25_dp * HOUR
    bounds(4) = bounds(3) + 9.50_dp * HOUR
    bounds(5) = bounds(4) + 0.25_dp * HOUR
    bounds(6) = bounds(5) + 1.00_dp * HOUR
    ! alternative bounds
    !bounds(1) = tBegin * HOUR 
    !bounds(2) = bounds(1) + 2.00_dp * HOUR
    !bounds(3) = bounds(2) + 2.00_dp * HOUR
    !bounds(4) = bounds(3) + 4.00_dp * HOUR
    !bounds(5) = bounds(4) + 2.00_dp * HOUR
    !bounds(6) = bounds(5) + 2.00_dp * HOUR
  END FUNCTION Set_pseudoLWCbounds

  FUNCTION pseudoLWC(RealTime)  RESULT(LWC)
    USE Kind_Mod,     ONLY: dp
    USE Control_Mod,  ONLY: constLWC, LWCLevelmin, LWCLevelmax

    REAL(dp) :: LWC
    REAL(dp) :: RealTime
    REAL(dp) :: Time
    !
    !         LWC
    !           /|\       fake LWC function (periodically)
    !            |
    !            |
    ! maxlwclvl _|_ _ _ ____             ____
    !            |     /|  |\           /    \
    !            |    /      \         /      \
    !            |   /  |  |  \       /        \
    !            |  /          \     /          \
    ! minlwclvl _|_/    |  |    \___/            \__...... 
    !------------+----------------------------------------------->
    !            | |    |  |    | |                               Time
    !     LWCb   1 2    3  4    5 6
    !
    IF ( constLWC ) THEN
      ! Constant LWC level 
      LWC = LWCconst
    ELSE
      ! calculate LWC level
      Time = MODULO(RealTime,LWCb(6))
      !
      LWC = LwcLevelmin

      ! --- Constant minimum LWC level
      IF     (LWCb(1)<=Time.AND.Time<LWCb(2)) THEN
        LWC = LwcLevelmin

      ! --- linear increase of LWC level till maximum
      ELSEIF (LWCb(2)<=Time.AND.Time<LWCb(3)) THEN
        LWC = (Time-LWCb(2))/(LWCb(3)-LWCb(2))*LwcLevelmax  +  &
        &     (LWCb(3)-Time)/(LWCb(3)-LWCb(2))*LwcLevelmin

      ! --- Constant maximum LWC level
      ELSEIF (LWCb(3)<=Time.AND.Time<LWCb(4)) THEN
        LWC = LwcLevelmax

      ! --- linear decrease of LWC level till minimum
      ELSEIF (LWCb(4)<=Time.AND.Time<LWCb(5)) THEN
        LWC = (Time-LWCb(4))/(LWCb(5)-LWCb(4))*LwcLevelmin  +  &
        &     (LWCb(5)-Time)/(LWCb(5)-LWCb(4))*LwcLevelmax

      ! --- Constant minium LWC level
      ELSEIF (LWCb(5)<=Time.AND.Time<=LWCb(6)) THEN
        LWC = LwcLevelmin
      END IF
    END IF
  END FUNCTION

  ! turn total lwc into an array of lwcs for each droplet class
  FUNCTION LWC_array(RealTime, waterMasses) RESULT(LWCs)
    USE Reac_Mod,    ONLY: DropletClasses
    USE Control_Mod, ONLY: nDropletClasses

    REAL(dp) :: RealTime, LWCs(nDropletClasses)
    REAL(dp), OPTIONAL :: waterMasses(nDropletClasses)

    REAL(dp) :: LWC_total

    IF (adiabatic_parcel) THEN ! copy class specific water mass
      ! * rho_parcel to convert water masses from kg/kg to l/m3
      IF (PRESENT(waterMasses)) THEN
        LWCs = waterMasses * milli * rho_parcel
      ELSE 
        LWCs = DropletClasses%waterMass * milli * rho_parcel
      END IF
    ELSE ! get water masses by pre-defined LWC function
      LWC_total = pseudoLWC(RealTime)

      LWCs(DropletClasses%active)   = LWC_total * DropletClasses%LWC_portion(DropletClasses%active)

      LWCs(DropletClasses%inactive) = DropletClasses%inactive_LWC
    END IF

  END FUNCTION LWC_array


  !**************************************************************************!
  !***
  !***  Computation of Hp concentration (pH-Value)
  !***
  !**************************************************************************!
  !
  FUNCTION HpValue_from_Electroneutrality(ConcAqua, LWC) RESULT(HpConc)
    USE Reac_Mod,    ONLY: ns_AQUA, Hp_ind, Charge
    USE Control_Mod, ONLY: ZERO
    USE Kind_Mod,    ONLY: dp

    REAL(dp) :: HpConc, IonCharge, LWC
    REAL(dp) :: ConcAqua(ns_AQUA)
    INTEGER :: j

    IonCharge = ZERO
    DO j=1,ns_AQUA
      IF (j/=Hp_ind .AND. ConcAqua(j)>1E-15_dp) THEN
        IonCharge = IonCharge + ConcAqua(j) * Charge(j)
      END IF
    END DO

    ! Old (does not equilibriate Water diss, i.e. adds only Hp or OHm, not both)
    !HpConc = - IonCharge

    ! New (to have zero net charge and Hp*OHm=10^-14)
    HpConc    = ZERO

    ! quadratic solution for 
    !   IonCharge + Hp + OHm = 0 
    ! and 
    !   Hp/(LWC*mol2part) * OHm/(LWC*mol2part) = 10**-14
    HpConc = - IonCharge/2 + SQRT((IonCharge/2)**2 + 1E-14_dp * LWC**2 * mol2part**2)

  END FUNCTION HpValue_from_Electroneutrality

  FUNCTION pressure_from_height(z) RESULT(pressure)
    ! IN:
    REAL(dp) :: z ! height
    ! OUT:
    REAL(dp) :: pressure ! corresponding pressure

    ! barometric height formula
    pressure = Pressure0 * EXP(-g*z / ( R_air * Temperature0 ))

  END FUNCTION pressure_from_height

!-- saturation water vapor mixing ratio (kg/kg)
    FUNCTION qsatw(T,p) RESULT(qs)
      REAL(dp) :: qs

      REAL(dp), INTENT(IN) :: T, p

      qs = R_air / Rv * esatw(T) / ( p - esatw(T) )

    END FUNCTION qsatw
!
!-- saturation water vapor mixing ratio (kg/kg)
    FUNCTION dqsatwdT(T,p) RESULT(dqsdT)
      REAL(dp) :: dqsdT

      REAL(dp), INTENT(IN) :: T, p

       dqsdT = R_air / Rv * desatwdT(T)*p / ( p - esatw(T) )**2

    END FUNCTION dqsatwdT
!
!-- saturation water vapor pressure (Pa) (Flatau et.al, 1992, JAM)
    FUNCTION esatw(T) RESULT(es)
      REAL(dp) :: es

      REAL(dp), INTENT(IN) :: T

      REAL(dp) :: dT
      REAL(dp) :: a(0:8) = (/ 6.11239921, 0.443987641, 0.142986287e-1,            &
                           0.264847430e-3, 0.302950461e-5, 0.206739458e-7,     &
                           0.640689451e-10, -0.952447341e-13, -0.976195544e-15/)

      dT = T - 273.15_dp
      es = a(0)+dT*(a(1)+dT*(a(2)+dT*(a(3)+dT*(a(4)+dT*(a(5)+dT*(a(6)+dT*(a(7)+a(8)*dT)))))))
      es = es * 100.0_dp

    END FUNCTION esatw

    FUNCTION desatwdT(T) RESULT(desdT)
      REAL(dp) :: desdT
      REAL(dp), INTENT(IN) :: T ! temperature (K)

      REAL(dp) :: a(0:8) = (/ 0.443956472    ,  0.285976452e-1 ,  0.794747212e-3,  &
      &                       0.121167162e-4 ,  0.103167413e-6 ,  0.385208005e-9,  &
      &                      -0.604119582e-12, -0.792933209e-14, -0.599634321e-17 /)
      REAL(dp) :: dT

      dT = T - 273.15_dp

      desdT = a(0)+dT*(a(1)+dT*(a(2)+dT*(a(3)+dT*(a(4)+dT*(a(5)+dT*(a(6)+dT*(a(7)+a(8)*dT))))))) 
    END FUNCTION desatwdT

    ! calculate S_eq, i.e. value of Köhler curve 
    SUBROUTINE calculate_all_S_eq(S_eq, m_w, r, n_s, T, dSeqdm)
      REAL(dp), DIMENSION(nDropletClasses), INTENT(IN) :: m_w, r, n_s
      REAL(dp)                            , INTENT(IN) :: T

      REAL(dp), DIMENSION(nDropletClasses), INTENT(OUT) :: S_eq
      REAL(dp), DIMENSION(nDropletClasses), INTENT(OUT), OPTIONAL :: dSeqdm

      REAL(dp), DIMENSION(nDropletClasses) :: afactor, n_w, solute_effect, curvature_effect &
                                         &  , dsolute_effectdm, dcurvature_effectdm

      !-- Kelvin term (curvature)
      afactor = 2.0 * sigma_air_liq(T) / ( milli * rho_H2O * Rv * T )
      curvature_effect = EXP( afactor / r )

      !-- Raoult term (solution)
      ! number of water molecules
      n_w = m_w / molw_H2O * SI_na
      solute_effect = n_w / ( n_w + n_s )

      ! value of Köhler curve
      S_eq = solute_effect * curvature_effect

      IF ( PRESENT(dSeqdm) ) THEN
        dcurvature_effectdm = afactor / (4*Pi*rho_H2O*r**4) * curvature_effect
        dsolute_effectdm = n_s / ( n_w + n_s )**2 / molw_H2O * SI_na
      
        dSeqdm = dsolute_effectdm * curvature_effect + solute_effect * dcurvature_effectdm
      END IF

    END SUBROUTINE calculate_all_S_eq

    FUNCTION calculate_one_S_eq(m_w, r, n_s, T) RESULT (S_eq)
      REAL(dp) :: m_w, r, n_s, T

      REAL(dp) :: S_eq

      REAL(dp) :: afactor, n_w, solute_effect, curvature_effect

      !-- Kelvin term (curvature)
      afactor = 2.0 * sigma_air_liq(T)  / ( milli * rho_H2O * Rv * T )
      curvature_effect = EXP( afactor / r )

      !-- Raoult term (solution)
      ! number of water molecules
      n_w = m_w / molw_H2O * SI_na
      solute_effect = n_w / ( n_w + n_s )

      ! value of Köhler curve
      S_eq = solute_effect * curvature_effect

    END FUNCTION calculate_one_S_eq

  SUBROUTINE calculate_saturated_radius(r_solution, V_solution_l, V_H2O_l, droplet_conc, V_solute, T)
    REAL(dp) :: r_solution, V_solution_l, V_H2O_l, droplet_conc(:), V_solute, T

    REAL(dp) :: a, b, V_H2O, V_solution    &
    &         , n_solute, r_H2O

    n_solute = SUM( droplet_conc * mega ) ! * mega for molec/cm3 -> molec/m3

    b = n_solute / ( SI_na * (FOUR/THREE) * Pi * rho_H2O / molw_H2O )

    a = 2 * sigma_air_liq(T) / (rho_H2O * Rv*0.001 * T)

    r_H2O = SQRT(b/a)
    V_H2O = (FOUR/THREE) * Pi * r_H2O**3

    ! cheeky correction, which happens to be a very good approximation
    ! just calculate radius of solution as radius of sphere of volume
    ! equal to sum of aerosol volumes and water volume
    ! (this is "cheeky", because sqrt(b/a) is derived under the assumption, that
    !  r_solution=r_H2O)
    V_solution   = V_H2O + V_solute ! ( = ideal solution assumption)
    r_solution   = (Pi34 * V_solution)**rTHREE

    V_H2O_l      = V_H2O*kilo       ! m3 -> liter unit for lwc
    V_solution_l = V_solution*kilo  ! m3 -> liter unit for lsc

  END SUBROUTINE calculate_saturated_radius

  SUBROUTINE calculate_RH_radius(r_eq, m_w, RH, droplet_conc, T_parcel, v_solute)
    REAL(dp) :: RH, T_parcel, v_solute
    REAL(dp), DIMENSION(:) :: droplet_conc

    REAL(dp) :: r_eq, m_w

    REAL(dp) :: rel_threshold, r, S_eq, n_s, m_pre, m_max, m_min
    INTEGER :: iter, max_steps = 500

    IF (RH>1) THEN ! solution not unique, abort
      STOP 'STOP in calc_RH_radius: RH greater than 1, köhler curve not invertible.'
    END IF

    ! exit criterion if borders of bisection have smaller relative difference than the following threshold
    rel_threshold = 1E-10

    ! set initial borders of bisection
    m_min = 1E-30
    m_max = 1E-3

    n_s = SUM(droplet_conc)

    ! do bisection
    DO iter = 1, max_steps
      ! try logarithmic mid point
      m_w = EXP( 0.5 * ( LOG(m_min) + LOG(m_max) ) )
      
      ! calculate corresponding solution radius (assuming ideal solution, i.e. volumes add up)
      !r = (Pi34 * (m_w / rho_h2o))**rTHREE
      r = (Pi34 * (m_w / rho_h2o + v_solute))**rTHREE

      S_eq = calculate_one_S_eq(m_w, r, n_s, T_parcel)

      IF ( S_eq - RH > ZERO ) THEN
        m_pre = m_max
        m_max = m_w
      ELSE
        m_pre = m_min
        m_min = m_w
      END IF

      IF ( ABS(1-LOG(m_w)/LOG(m_pre)) < rel_threshold ) EXIT

    END DO

    r_eq = r

  END SUBROUTINE calculate_RH_radius

  ! calculate critical RH and r, i.e. where dKöhler/dr = 0
  SUBROUTINE calculate_critical_RH_and_r(critical_S, critical_r, droplet_conc, T_parcel, v_solute)
    REAL(dp) :: T_parcel, v_solute
    REAL(dp), DIMENSION(:) :: droplet_conc

    REAL(dp) :: critical_S, critical_r

    REAL(dp) :: rel_threshold, r_min, r_max, r, S_eq_plus, S_eq_minus, n_s, r_pre, dKoehlerdr, h, r_minus, r_plus
    REAL(dp) :: m_w_minus, m_w_plus
    INTEGER :: iter, max_steps = 500

    ! exit criterion if borders of bisection have smaller relative difference than the following threshold
    rel_threshold = 1E-15

    ! relative step size for approximate derivative of köhler curve
    h = 1E-15

    ! set initial borders of bisection
    r_min = 1.01_dp*((Pi34*v_solute)**rTHREE)
    r_max = 1E-2

    n_s = SUM(droplet_conc)

    ! do bisection
    DO iter = 1, max_steps
      ! try logarithmic mid point

      r       = EXP( 0.5 * ( LOG(r_min) + LOG(r_max) ) )
      r_plus  = r + h/TWO
      r_minus = r - h/TWO

      m_w_minus = (Pi43 * r_minus**THREE - v_solute) * rho_h2o
      m_w_plus  = (Pi43 * r_plus**THREE  - v_solute) * rho_h2o

      S_eq_minus = calculate_one_S_eq(m_w_minus, r_minus, n_s, T_parcel)
      S_eq_plus  = calculate_one_S_eq(m_w_plus, r_plus, n_s, T_parcel)

      dKoehlerdr = ( S_eq_plus - S_eq_minus ) / (LOG(r_plus)-LOG(r_minus))

      IF ( dKoehlerdr > ZERO ) THEN
        r_pre = r_min
        r_min = r
      ELSE
        r_pre = r_max
        r_max = r
      END IF

      IF ( ABS(1-LOG(r)/LOG(r_pre)) < rel_threshold ) EXIT
    END DO
    
    critical_S = ( S_eq_plus + S_eq_minus ) / TWO
    critical_r = r
  END SUBROUTINE calculate_critical_RH_and_r
!
!-- Surface tension between liquid water and air (in J/m2)
    FUNCTION sigma_air_liq(tabs) RESULT(sigma)
       REAL(dp) :: sigma

       REAL(dp), INTENT(IN) :: tabs

       REAL(dp)             :: tabs_c

       tabs_c = tabs - 273.15_dp
!
!--    Pruppacher and Klett (1997), Eq. 5-12 
!      (this works only in a limited temperature range)
       !sigma = 75.93 + 0.115 * tabs_c    + 6.818e-2 * tabs_c**2 +      &
       !                     6.511e-3 * tabs_c**3 + 2.933e-4 * tabs_c**4 +      &
       !                     6.283e-6 * tabs_c**5 + 5.285e-8 * tabs_c**6
       !sigma = sigma_air_liq * 1.0E-3

!--    Straka (2009)
       !sigma = 0.0761 - 0.000155 * tabs_c

       ! Vargaftik et al., international tables of surface tension of water
       sigma = 0.07564_dp - 0.00014606_dp * tabs_c

    END FUNCTION sigma_air_liq

    FUNCTION thermal_conductivity(T) RESULT(thermalconductivity)
      REAL(dp) :: thermalconductivity
      REAL(dp) :: T

      ! thermal conductivity of air (Rogers and Yau, Table 7.1)
      thermalconductivity = 7.94048e-5_dp * T + 0.00227011_dp

    END FUNCTION thermal_conductivity

    FUNCTION diff_coeff(T, pressure) RESULT(diffcoeff)
      REAL(dp) :: diffcoeff
      REAL(dp) :: T, pressure

      ! Molecular diffusivity of water vapor in air in m2/s (Hall und Pruppacher, 1976)
      diffcoeff = 0.211E-4_dp * ( T / 273.15_dp )**1.94_dp * ( 101325.0_dp / pressure )

    END FUNCTION diff_coeff

    FUNCTION get_wet_radii(waterMass, LWCs) RESULT(radii)
      USE Reac_Mod,    ONLY: DropletClasses

      REAL(dp), OPTIONAL :: waterMass(nDropletClasses), LWCs(nDropletClasses)

      REAL(dp) :: radii(nDropletClasses)

      IF ( PRESENT(LWCs) ) THEN
        ! LWCs is l/m3
        !            ->g/m3      ->g/kg     ->m3/kg ->m3                    ->m
        radii = (Pi34*kilo*LWCs/(rho_parcel*rho_H2O*DropletClasses%Number))**(rTHREE)
      ELSE IF ( PRESENT(waterMass) ) THEN
        ! waterMass is g/kg
        !                       ->m3/kg ->m3                    ->m
        radii = (Pi34*waterMass/(rho_H2O*DropletClasses%Number))**(rTHREE)
      ELSE
        ! DropletClasses%waterMass is in g/kg
        !                                      ->m3/kg ->m3                    ->m
        radii = (Pi34*DropletClasses%waterMass/(rho_H2O*DropletClasses%Number))**(rTHREE)
      END IF

    END FUNCTION get_wet_radii

END MODULE Meteo_Mod
