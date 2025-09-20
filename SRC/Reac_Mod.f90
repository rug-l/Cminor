!
! Copyright (C) 2025 Levin Rug, Willi Schimmel (E-Mail: l.rug@lmu.de)
! See ./SRC/Cminor.f90 for the copyright notice
! See ./LICENSE for license information
! SPDX-License-Identifier: GPL-3.0
!
!==============================================================
!===  MODULE for the Description of Chemistry, Deposition
!===  and Emissions
!==============================================================
!
  MODULE Reac_Mod
    USE Kind_Mod,    ONLY: dp
    USE Control_Mod, ONLY: LenName

    CHARACTER(9) :: measure_gas_ph(2)=(/"molec/cm3", "mol/m3   "/)
    CHARACTER(9) :: measure_aqua_ph(1)=(/"mol/l"/)
    CHARACTER(12):: units(2)

!--------------------------------------------------------------
!---  Reaction Mechanism
!--------------------------------------------------------------
!
!---  reaction types
    TYPE def_para
      CHARACTER(9) :: name_type   ! name of the reaction type
      INTEGER      :: n_par       ! number of parameters for corresponding reaction type
      INTEGER      :: n_reac      ! number of corresponding reactions
      LOGICAL      :: act         ! true if n_reac>0
    END TYPE def_para

    TYPE val_para
      INTEGER,  ALLOCATABLE :: iR(:,:)
      REAL(dp), ALLOCATABLE :: vR(:,:)
    END TYPE val_para


    INTEGER, PARAMETER ::           &
    &            iPHOTABC   = 1  &
    &         ,  iPHOTMCM   = 2  &
    &         ,  iPHOTAB    = 3  &
    &         ,  iCONST     = 4  &
    &         ,  iTEMP      = 5  &
    &         ,  iTEMP1     = 6  &
    &         ,  iTEMP2     = 7  &
    &         ,  iTEMP3     = 8  &
    &         ,  iTEMP4     = 9  &
    &         ,  iASPEC1    = 10  &
    &         ,  iASPEC2    = 11  &
    &         ,  iASPEC3    = 12  &
    &         ,  iDCONST    = 13  &
    &         ,  iDTEMP     = 14  &
    &         ,  iDTEMP2    = 15  &
    &         ,  iDTEMP3    = 16  &
    &         ,  iDTEMP4    = 17  &
    &         ,  iDTEMP5    = 18  &
    &         ,  iT1H2O     = 19  &
    &         ,  iS4H2O     = 20  &
    &         ,  iTROE      = 21  &
    &         ,  iTROEQ     = 22  &
    &         ,  iTROEF     = 23  &
    &         ,  iTROEQF    = 24  &
    &         ,  iTROEXP    = 25  &
    &         ,  iTROEMCM   = 26  &
    &         ,  iSPEC1     = 27  &
    &         ,  iSPEC2     = 28  &
    &         ,  iSPEC3     = 29  &
    &         ,  iSPEC4     = 30  &
    &         ,  iSPEC1MCM  = 31  &
    &         ,  iSPEC2MCM  = 32  &
    &         ,  iSPEC3MCM  = 33  &
    &         ,  iSPEC4MCM  = 34  &
    &         ,  iSPEC5MCM  = 35  &
    &         ,  iSPEC6MCM  = 36  &
    &         ,  iSPEC7MCM  = 37  &
    &         ,  iSPEC8MCM  = 38  &
    &         ,  iSPEC9MCM  = 39  &
    &         ,  iHOM1      = 40  &
    &         ,  iPHOTO     = 41  &
    &         ,  iPHOTO2    = 42  &
    &         ,  iPHOTO3    = 43  &
    &         ,  iSPECIAL   = 44

    INTEGER, PARAMETER :: nReacTypes = 44

    TYPE(val_para), DIMENSION(nReacTypes) :: reac_val

    TYPE(def_para), DIMENSION(nReacTypes) :: reac_par = &
                   (/def_para("PHOTABC",   3,  0, .FALSE.), &
                     def_para("PHOTMCM",   3,  0, .FALSE.), &
                     def_para("PHOTAB",    2,  0, .FALSE.), &
                     def_para("CONST",     1,  0, .FALSE.), &
                     def_para("TEMP",      3,  0, .FALSE.), &
                     def_para("TEMP1",     2,  0, .FALSE.), &
                     def_para("TEMP2",     2,  0, .FALSE.), &
                     def_para("TEMP3",     2,  0, .FALSE.), &
                     def_para("TEMP4",     2,  0, .FALSE.), &
                     def_para("ASPEC1",    2,  0, .FALSE.), &
                     def_para("ASPEC2",    3,  0, .FALSE.), &
                     def_para("ASPEC3",    2,  0, .FALSE.), &
                     def_para("DCONST",    2,  0, .FALSE.), &
                     def_para("DTEMP",     3,  0, .FALSE.), &
                     def_para("DTEMP2",    4,  0, .FALSE.), &
                     def_para("DTEMP3",    4,  0, .FALSE.), &
                     def_para("DTEMP4",    3,  0, .FALSE.), &
                     def_para("DTEMP5",    3,  0, .FALSE.), &
                     def_para("T1H2O",     2,  0, .FALSE.), &
                     def_para("S4H2O",     4,  0, .FALSE.), &
                     def_para("TROE",      4,  0, .FALSE.), &
                     def_para("TROEQ",     6,  0, .FALSE.), &
                     def_para("TROEF",     5,  0, .FALSE.), &
                     def_para("TROEQF",    7,  0, .FALSE.), &
                     def_para("TROEXP",    5,  0, .FALSE.), &
                     def_para("TROEMCM",  10,  0, .FALSE.), &
                     def_para("SPEC1",     2,  0, .FALSE.), &
                     def_para("SPEC2",     2,  0, .FALSE.), &
                     def_para("SPEC3",     6,  0, .FALSE.), &
                     def_para("SPEC4",     4,  0, .FALSE.), &
                     def_para("SPEC1MCM",  3,  0, .FALSE.), &
                     def_para("SPEC2MCM",  3,  0, .FALSE.), &
                     def_para("SPEC4MCM",  4,  0, .FALSE.), &
                     def_para("SPEC3MCM",  2,  0, .FALSE.), &
                     def_para("SPEC5MCM",  4,  0, .FALSE.), &
                     def_para("SPEC6MCM",  4,  0, .FALSE.), &
                     def_para("SPEC7MCM",  6,  0, .FALSE.), &
                     def_para("SPEC8MCM",  4,  0, .FALSE.), &
                     def_para("SPEC9MCM",  4,  0, .FALSE.), &
                     def_para("HOM1",      4,  0, .FALSE.), &
                     def_para("PHOTO",     4,  0, .FALSE.), &
                     def_para("PHOTO2",    4,  0, .FALSE.), &
                     def_para("PHOTO3",    4,  0, .FALSE.), &
                     def_para("SPECIAL",   0,  0, .FALSE.)  /) 


!---  reaction structures
    TYPE reactant
       INTEGER :: i_spc
       REAL(dp) :: d_koef
    END TYPE reactant

    TYPE reaction
       CHARACTER(8):: str_class
       CHARACTER(12):: str_type
       CHARACTER(20), POINTER :: factor

       INTEGER :: type
       INTEGER :: anz_p
       INTEGER :: n_so
       INTEGER :: n_si
       INTEGER :: n_siso
       INTEGER :: n_so_a
       INTEGER :: n_si_a
       INTEGER :: n_block
       INTEGER :: diag_flag
       INTEGER :: nstr

       INTEGER, POINTER :: so(:)
       INTEGER, POINTER :: si(:)
       INTEGER, POINTER :: struct(:,:)
       INTEGER, POINTER :: reac_si_nr(:)
       INTEGER, POINTER :: reac_so_nr(:)

       REAL(dp), POINTER  :: last_rate(:,:)
       REAL(dp), POINTER  :: back_rate(:,:)
       REAL(dp), POINTER  :: v(:,:)
       REAL(dp), POINTER  :: w(:,:)
       REAL(dp), POINTER :: dparam(:)  

       LOGICAL :: odd = .FALSE.

       REAL(dp) :: fac_exp
       REAL(dp) :: fac_A             

       TYPE (reactant), POINTER :: reactant(:)

       TYPE (reaction), POINTER :: next
       TYPE (reaction), POINTER :: next_all
    END TYPE reaction

    TYPE reactype
      CHARACTER(80) :: name=''
      INTEGER       :: n=0
      LOGICAL       :: exist=.FALSE.
    END TYPE reactype

!--------------------------------------------------------------
!--   dimensions
    INTEGER :: nt=0           ! unused?
    INTEGER :: nFrac=0        ! number of aquatic droplett classes
    INTEGER :: nspc=0         ! Number of all species excluding katalytic/passive species
    INTEGER :: nspc2=0        ! same as nspc, but with n droplet classes
    INTEGER :: nreac=0        ! number of all reactions, HENRY/DISS counts as two reactions
    INTEGER :: nreac2=0       ! same as nr, but with n droplet classes
    INTEGER :: neq=0          ! number of all reactions, HENRY/DISS counts as two reactions
    INTEGER :: neq2=0         ! same as neq, but with n droplet classes
    INTEGER :: nPhases=0      ! number of involed species phases
    INTEGER :: nDIM=0         ! Dimension of ODE system
    INTEGER :: nDIM2=0        ! same as nDIM, but with n droplet classes
    REAL(dp) :: rNspc         ! real value 1/nspc for error calcualtion

    ! number of species in each phase/state
    INTEGER :: ns_KAT = 0, ns_KAT2 = 0
    INTEGER :: ns_GAS = 0, ns_AQUA = 0

    ! number of each reaction class
    INTEGER :: nr_gas = 0, nr_aqua = 0, nr_henry = 0, nr_diss = 0, nr_special = 0, nr_liquid = 0 ! = nr_diss+nr_aqua

    LOGICAL :: hasGasSpc=.FALSE., hasAquaSpc=.FALSE., hasSolidSpc=.FALSE., hasPartiSpc=.FALSE.
    LOGICAL :: hasGasReac  =.FALSE., hasAquaReac=.FALSE. &
           & , hasHenryReac=.FALSE., hasDissReac=.FALSE., hasLiquidReac=.FALSE.

    LOGICAL :: hasPhotoReac=.FALSE.

    ! number of each reaction type
    INTEGER :: nr_G_photo = 0, nr_G_const = 0, nr_G_temp = 0, nr_G_troe =0, nr_G_spec = 0, nr_G_lind = 0
    INTEGER :: nr_A_photo = 0, nr_A_const = 0, nr_A_temp = 0, nr_A_spec =0
    INTEGER :: nr_S_temp  = 0, nr_S_equi  = 0, nr_S_spec = 0
    INTEGER :: nr_SimpTB  = 0, nr_press   = 0
    INTEGER :: nr_G_special = 0, nr_A_special = 0,nr_S_special = 0,nr_P_special = 0,nr_M_special = 0
    INTEGER :: nr_H_special = 0, nr_D_special = 0, nr_D_Temp

    INTEGER :: nr_PHOTabc = 0, nr_PHOTab = 0, nr_PHOTmcm = 0, nr_CONST = 0
    INTEGER :: nr_TEMP = 0, nr_TEMP1 = 0, nr_TEMP2 = 0, nr_TEMP3 = 0, nr_TEMP4 = 0
    INTEGER :: nr_TROE  = 0, nr_TROEf = 0, nr_TROEq = 0, nr_TROEqf = 0, nr_TROExp = 0, nr_TROEmcm = 0
    INTEGER :: nr_SPEC1 = 0, nr_SPEC2 = 0, nr_SPEC3 = 0, nr_SPEC4 = 0
    INTEGER :: nr_SPEC1mcm = 0, nr_SPEC2mcm = 0, nr_SPEC3mcm = 0, nr_SPEC4mcm = 0
    INTEGER :: nr_SPEC5mcm = 0, nr_SPEC6mcm = 0, nr_SPEC7mcm = 0, nr_SPEC8mcm = 0, nr_SPEC9mcm = 0
    INTEGER :: nr_S4H2O  = 0, nr_T1H2O   = 0
    INTEGER :: nr_ASPEC1 = 0, nr_ASPEC2  = 0, nr_ASPEC3   = 0
    INTEGER :: nr_DTEMP  = 0, nr_DTEMP2  = 0, nr_DTEMP3   = 0, nr_DTEMP4 = 0, nr_DTEMP5 = 0, nr_DCONST = 0
    INTEGER :: nr_HENRYga = 0, nr_HENRYag  = 0
    INTEGER :: nr_FACTOR = 0, nr_FAC_H2  = 0, nr_FAC_O2N2 = 0, nr_FAC_M  = 0, nr_FAC_O2 = 0, nr_FAC_N2 = 0
    INTEGER :: nr_FAC_H2O = 0, nr_FAC_RO2 = 0, nr_FAC_O2O2 = 0, nr_FAC_aH2O = 0, nr_FAC_RO2aq = 0
    INTEGER :: nr_HOaqua = 0, nr_SOaqua = 0, nr_TOaqua = 0      ! higher order aqueous reactions
    INTEGER :: nr_PHOTOkpp = 0, nr_PHOTO2kpp = 0, nr_PHOTO3kpp = 0
    INTEGER :: nr_HOM1 = 0

    LOGICAL :: PHOTO=.FALSE.

!    INTEGER :: nreakstemp,nreaksequi,nreaksspec
!--    define indices of special species
    INTEGER :: Hp_ind,      &   ! Index Hp (aqua)
             & OHm_ind,     &   ! Index OHm  (aqua)
             & H2O_ind,     &   ! Index H2O (gas)
             & H2_ind,      &   ! Index H2O (gas)
             & N2_ind,      &   ! Index H2O (gas)
             & O2_ind,      &   ! Index H2O (gas)
             & aH2O_ind,    &   ! Index aH2O (aqua)
             & Temp_ind         ! Index Temperatur

    INTEGER, ALLOCATABLE :: gaseous_passive_ind(:)   ! indices of gaseous passive species

!--    mass transfer coefficient
    REAL(dp) :: dkmt = 1.0d3      ! Standard Mass Trasfer Coefficient
!
!--------------------------------------------------------------
!--    species names 
    CHARACTER(LenName), ALLOCATABLE :: y_name(:)
!
!--------------------------------------------------------------
!-- Netcdf output and diagnostic variables 
    CHARACTER(100), ALLOCATABLE    :: Diag_Name(:)
    INTEGER, ALLOCATABLE           :: Diag_Index(:)
    CHARACTER(60), ALLOCATABLE  :: Diag_Name_Netcdf(:) &
&                                , Diag_LongName(:)     

!--------------------------------------------------------------
!--    Peroxyradicals
    LOGICAL              :: hasRO2
    INTEGER, ALLOCATABLE :: RO2(:)          ! Species-Index of peroxyradical
    INTEGER              :: nRO2            ! Number of Gasphase-Peroxyradicals in mechanism
!--    aqueous phase
    LOGICAL              :: hasRO2aq
    INTEGER, ALLOCATABLE :: RO2aq(:)        ! Species-Index of peroxyradical
    INTEGER              :: nRO2aq          ! Number of aqueos phase peroxy radicals in mechanism
!
!--------------------------------------------------------------
!--    Aerosol species properties
    REAL(dp), ALLOCATABLE :: Charge(:)       ! charge of ions
    !REAL(dp), ALLOCATABLE :: SolubInd(:)     ! solubility index
    REAL(dp), ALLOCATABLE :: MolMass(:)      ! molar mass of species
    !REAL(dp), ALLOCATABLE :: OrgIndex(:)     ! carbon atoms
    !CHARACTER(2), ALLOCATABLE :: CC(:)      ! compound class

!--------------------------------------------------------------
!--    input arrays
    REAL(dp), ALLOCATABLE :: InitValAct(:) ! active species' concentrations
    REAL(dp), ALLOCATABLE :: InitValKat(:), InitValKat_Ref(:) ! passive species' concentrations

    REAL(dp), ALLOCATABLE :: henry_diff(:)
    REAL(dp), ALLOCATABLE :: henry_accom(:)

!--------------------------------------------------------------
!--    deposition and emissions
    REAL(dp), ALLOCATABLE :: y_emi(:), y_emi_1nD(:)
    REAL(dp), ALLOCATABLE :: y_depos(:), y_depos_1nD(:)

!--------------------------------------------------------------

    CHARACTER(2),  ALLOCATABLE :: ThAtoms(:)           ! (:,4)
    REAL(dp),      ALLOCATABLE :: lowA(:),lowB(:),lowC(:),lowD(:),lowE(:),lowF(:),lowG(:)
    REAL(dp),      ALLOCATABLE :: highA(:),highB(:),highC(:),highD(:),highE(:),highF(:),highG(:)
    REAL(dp),      ALLOCATABLE :: SwitchTemp(:)
    !
    REAL(dp), ALLOCATABLE :: GFE(:), DGFEdT(:)
    REAL(dp), ALLOCATABLE :: DelGFE(:), DDelGFEdT(:)
    !
    REAL(dp), ALLOCATABLE :: MW(:)    ! molecular weight (combustion)
    REAL(dp), ALLOCATABLE :: rMW(:)    !1/ molecular weight (combustion)
    !
    ! more speedchem stuff
    REAL(dp) :: rho, rRho   ! rho = density, rRho=kilo/rho
    INTEGER, ALLOCATABLE :: SCperm(:)

    !
    ! index arrays for different reaction types
    ! Types for vectorised version
    TYPE ReacTypeIndex_CK
      INTEGER, ALLOCATABLE :: iArr(:),  iLind(:), iTroe(:)
      INTEGER, ALLOCATABLE :: iEqui(:), iXrev(:),  iTBody(:)
      INTEGER, ALLOCATABLE :: iTBodyExtra(:), iHigh(:), iLow(:)
      INTEGER :: nArr, nLind, nTroe
      INTEGER :: nEqui,nXrev, nTBody
      INTEGER :: nTBodyExtra, nHigh, nLow
    END TYPE ReacTypeIndex_CK

    TYPE(ReacTypeIndex_CK) :: RTind

    !
    ! index arrays for different reaction types
    ! Types for vectorised version
    TYPE ReacTypeIndex_TR
      INTEGER, ALLOCATABLE :: iPHOTabc(:),  iPHOTab(:),   iPHOTmcm(:), iCONST(:)  &
      &                     , iTEMP(:),    iTEMP1(:),    iTEMP2(:),    iTEMP3(:), iTEMP4(:)  &
      &                     , iTROE(:),     iTROEf(:),    iTROEq(:),   iTROEqf(:), iTROExp(:) &
      &                     , iTROEmcm(:),  iSPEC1(:),    iSPEC2(:),   iSPEC3(:),  iSPEC4(:)  &
      &                     , iSPEC1mcm(:), iSPEC2mcm(:), iSPEC3mcm(:), iSPEC4mcm(:)  &
      &                     , iSPEC5mcm(:), iSPEC6mcm(:), iSPEC7mcm(:), iSPEC8mcm(:), iSPEC9mcm(:)  &
      &                     , iS4H2O(:),    iT1H2O(:),    iASPEC1(:),   iASPEC2(:)    &
      &                     , iASPEC3(:),   iDCONST(:,:),   iDTEMP(:,:) &
      &                     , iDTEMP2(:,:), iDTEMP3(:,:), iDTEMP4(:,:), iDTEMP5(:,:) &
      &                     , iHENRY(:,:) , iHOM1(:) &
      &                     , iFAC_H2(:),  iFAC_O2N2(:), iFAC_M(:),    iFAC_O2(:),   iFAC_N2(:) &
      &                     , iFAC_H2O(:), iFAC_RO2(:),  iFAC_O2O2(:), iFAC_aH2O(:), iFAC_RO2aq(:) &
      &                     , HOaqua(:)      &
      &                     , iHOaqua(:) &
      &                     , iSOaqua(:) &
      &                     , iTOaqua(:) &
      &                     , iPHOTOkpp(:), iPHOTO2kpp(:), iPHOTO3kpp(:) &
      &                     , iPHOTOStWe(:) &
      &                     , iSPECIAL(:)
      REAL(dp), ALLOCATABLE :: PHOTabc(:,:)   &
      &                     ,  PHOTab(:,:)    &
      &                     ,  PHOTmcm(:,:)   &
      &                     ,  CONST(:)       &
      &                     ,  TEMP(:,:), TEMP1(:,:), TEMP2(:,:), TEMP3(:,:), TEMP4(:,:) &
      &                     ,  TROE(:,:)      &
      &                     ,  TROEf(:,:)     &
      &                     ,  TROEq(:,:)     &
      &                     ,  TROEqf(:,:)    &
      &                     ,  TROExp(:,:)    &
      &                     ,  TROEmcm(:,:)   &
      &                     ,  SPEC1(:,:), SPEC2(:,:)  &
      &                     ,  SPEC3(:,:)     &
      &                     ,  SPEC4(:,:)     &
      &                     ,  SPEC1mcm(:,:), SPEC2mcm(:,:)  &
      &                     ,  SPEC3mcm(:,:)   &
      &                     ,  SPEC4mcm(:,:), SPEC5mcm(:,:), SPEC6mcm(:,:), SPEC8mcm(:,:) &
      &                     ,  SPEC7mcm(:,:)  &
      &                     ,  SPEC9mcm(:,:)  &
      &                     ,  S4H2O(:,:)     &
      &                     ,  T1H2O(:,:)     &
      &                     ,  ASPEC1(:,:), ASPEC3(:,:) &
      &                     ,  ASPEC2(:,:)&
      &                     ,  DTEMP(:,:), DTEMP5(:,:)  &
      &                     ,  DTEMP2(:,:), DTEMP3(:,:), DTEMP4(:,:) &
      &                     ,  DCONST(:,:)    &
      &                     ,  HENRY(:,:)     &
      &                     ,  PHOTOkpp(:), PHOTO2kpp(:), PHOTO3kpp(:) &
      &                     ,  PHOTOStWe(:)   &
      &                     ,  HOM1(:,:)
    END TYPE ReacTypeIndex_TR

    TYPE(ReacTypeIndex_TR) :: iR


    TYPE ReacTypeParameter_CK
      ! different kinds of arrhenius parameter
      REAL(dp), ALLOCATABLE :: A(:),    b(:),    E(:)
      REAL(dp), ALLOCATABLE :: A0(:),   b0(:),   E0(:)
      REAL(dp), ALLOCATABLE :: AX(:),   bX(:),   EX(:)
      REAL(dp), ALLOCATABLE :: Ainf(:), binf(:), Einf(:)
      ! Troe parameter
      REAL(dp), ALLOCATABLE :: T1(:), T2(:), T3(:), T4(:)
    END TYPE ReacTypeParameter_CK

    TYPE(ReacTypeParameter_CK) :: RTpar

    INTEGER, ALLOCATABLE  :: iFO_kat(:,:) ! reaction number where stoech coef == ONE
    INTEGER, ALLOCATABLE  :: iFO(:,:) ! reaction number where stoech coef == ONE
    INTEGER, ALLOCATABLE  :: iSO(:,:) ! reactions where species have second order 
    INTEGER, ALLOCATABLE  :: iHO(:,:) ! higher order reactions
    REAL(dp), ALLOCATABLE :: aHO(:)   ! higher order reactions contains also fractions and noninteger values
    INTEGER :: nFirst_order=0, nSecond_order=0, nHigher_order=0, nFirst_orderKAT=0

    INTEGER, SAVE :: bGs(2),bAs(2)        ! phase boundaries species
    INTEGER, SAVE :: bGr(2),bHr(2),bAr(2) ! phase boundaries reactions

    INTEGER, ALLOCATABLE :: iGs(:),iAs(:)        ! indices phases species
    INTEGER, ALLOCATABLE :: iGr(:),iHr(:),iAr(:) ! indices phases reactions

    ! second versions are the same but considering n droplet classes (non-compressed versions)
    INTEGER, SAVE :: bGs2(2), bAs2(2)          ! phase boundaries species
    INTEGER, SAVE :: bGr2(2), bHr2(2), bAr2(2) ! phase boundaries reactions

    INTEGER, ALLOCATABLE :: iGs2(:),iAs2(:)         ! indices phases species
    INTEGER, ALLOCATABLE :: iGr2(:),iHr2(:),iAr2(:) ! indices phases reactions

    ! indices of (aqueous) reactions to suppress if solute conc in droplet exceeds threshold (e.g. 0.1 mol/l)
    INTEGER, ALLOCATABLE :: iCutOffReacs(:)

    ! indices of parcel equations
    INTEGER :: iAqMassEq, iTeq, iqEq, iRhoEq, iZeq, iTeq2, iqEq2, iRhoEq2, iZeq2
    INTEGER, ALLOCATABLE :: iAqMassEq2(:)

  TYPE AFrac_T
    CHARACTER(LenName), ALLOCATABLE :: Species(:)
    REAL(dp),           ALLOCATABLE :: Frac1(:)     ! [g/g]
  END TYPE AFrac_T

  TYPE(AFRAC_T), ALLOCATABLE :: AFrac(:)

  TYPE Modes_T
    REAL(dp), ALLOCATABLE :: Radius(:)         ! [m] radius particle
    REAL(dp), ALLOCATABLE :: Number(:)         ! [#/m3]
    REAL(dp), ALLOCATABLE :: Std_Deviation(:)  ! [-]
    REAL(dp), ALLOCATABLE :: Density(:)        ! [g/m3]
  END TYPE Modes_T

  TYPE(Modes_T) :: Mode

  TYPE DropletClasses_T
    REAL(dp), ALLOCATABLE :: dryRadius(:)   ! [m]         radius aerosol
    REAL(dp), ALLOCATABLE :: dryVolume(:)   ! [m3]        volume of dry aerosol
    REAL(dp), ALLOCATABLE :: waterMass(:)   ! [g/kg]      total mass of water of superdroplet (prognostic variable)
    REAL(dp), ALLOCATABLE :: wetRadius(:)   ! [m]         radius of solution droplet
    REAL(dp), ALLOCATABLE :: saturated_air_radius(:)  ! [m]         radius haze droplet at RH=100% (not-activated aerosol)
    REAL(dp), ALLOCATABLE :: saturated_air_LWC(:) ! [l/m3]      amount of water in haze droplet
    REAL(dp), ALLOCATABLE :: saturated_air_LSC(:)     ! [l/m3]      liquid solution content
    REAL(dp), ALLOCATABLE :: inactive_LWC(:)   ! [l/m3]      irrelevant amount of water in inactive droplet
    REAL(dp), ALLOCATABLE :: inactive_r     ! [m]         irrelevant amount of water in inactive droplet
    REAL(dp), ALLOCATABLE :: Number(:)      ! [#/kg]      number of droplets per class
    REAL(dp), ALLOCATABLE :: Number0(:)     ! [#/kg]      number of droplets per class at beginning
    REAL(dp), ALLOCATABLE :: Conc(:,:)      ! [molec/cm3] concentrations in this droplet (at beginning)
    REAL(dp), ALLOCATABLE :: LWC_portion(:) !             Factor of LWC each droplet class gets
    REAL(dp), ALLOCATABLE :: critical_RH(:) !             critical relative humidity (activation RH     dKöhler/dr = 0)
    REAL(dp), ALLOCATABLE :: critical_r(:)  ! [m]         critical radius            (activation radius dKöhler/dr = 0)
    REAL(dp)              :: totalNumber    ! [#/kg]      total number of all droplets
    INTEGER,  ALLOCATABLE :: active(:)      !             indices of droplet classes that get activated (r>r_min)
    INTEGER,  ALLOCATABLE :: inactive(:)    !             indices of non-activated droplet classes (will remain at haze state)
  END TYPE DropletClasses_T

  TYPE(DropletClasses_T) :: DropletClasses
 
END MODULE Reac_Mod
