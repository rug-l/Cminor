!
! Copyright (C) 2025 Levin Rug, Willi Schimmel (E-Mail: l.rug@lmu.de)
! See ./SRC/Cminor.f90 for the copyright notice
! See ./LICENSE for license information
! SPDX-License-Identifier: GPL-3.0
!
!===================================================================================!
!                                                                                   !
!                     Module for reading the chemical system                        !
!                 file .sys and building the coeficient matrices                    !
!                                                                                   !
!===================================================================================!
MODULE Chemsys_Mod
  !
  USE Kind_Mod,          ONLY: dp
  
  USE Meteo_Mod,         ONLY: mol2part, HpValue_from_Electroneutrality, pseudoLWC, RefM, &
                             & O2, N2, H2O, N2O2, H2, RefTemp, RefH2, RefO2, LWC_array,   &
                             & SI_na, rho_H2O, molw_H2O, Rv, SI_gas, calculate_RH_radius, &
                             & calculate_saturated_radius, calculate_critical_RH_and_r,   &
                             & RefN2, RefH2O, RefH2, RefPressure, get_wet_radii
   
  USE Sparse_Mod,        ONLY: CSR_Matrix_T, A, B, BA, CompressList, CompressDoubleArray, &
                             & CompressIntegerArray, SortVecAsc
   
  USE LexicalStringSort, ONLY: StringSort

  USE hashtbl
  
  USE UniRnk_Mod,        ONLY: unirnk
   
  USE Control_Mod,       ONLY: ZERO, LenType, LenLine, iNcdfGas, iNcdfAqua, nD_reac, Pi,  &
                             & iNcdfGas2, iNcdfAqua2, nNcdfGas, UnitGas, tBegin, RH0,     &
                             & nNcdfAqua, nNcdfGas2, nNcdfAqua2, Pressure0, milli, Pi34,  &
                             & ChemUnit, ChemFile, iNcdfEmiss, Temperature0, q_parcel,    &
                             & THREE, FOUR, kilo, ask_for_continuing, z_parcel, ONE,      &
                             & nDropletClasses, rTWO, nD_spc, nD_KAT, T_parcel, Pi43,     &
                             & combustion, mONE, SysFile, LenLineChar, UnitAqua, mega,    &
                             & tEnd, HOUR, hourday, iDate, rlat, rlon, rho_parcel, TWO,   &
                             & phSet, nD_Ptr_spc, nD_Ptr_KAT, nD_Ptr_reacs, nNcdfEmiss,   &
                             & activation_radius, rTHREE, LWCLevelmin, LWCLevelmax,       &
                             & adiabatic_parcel, DropletClassPrint
   
  USE Reac_Mod           ! ONLY: nearly all of reac_mod.. 
   
  USE InputTool_Mod,     ONLY: OpenIniFile, LineFile, InputUnit_Initials, ClearIniFile,   &
                             & RewindFile, CloseIniFile
  !
  !
  IMPLICIT NONE
  !
  !
  INTEGER, PARAMETER ::  maxLENinActDuct=9
  ! 
  TYPE Duct_T
    CHARACTER(LenName) :: Species=''
    CHARACTER(LenType) :: Type
    REAL(dp)           :: Koeff
    INTEGER            :: iSpecies=0
  END TYPE Duct_T

  TYPE Special_T
    INTEGER                         :: nVariables = 0
    INTEGER,            ALLOCATABLE :: iVariables(:)
    CHARACTER(LenName), ALLOCATABLE :: cVariables(:)
    CHARACTER(LenLine)              :: Formula = ''
    LOGICAL                         :: Temp = .FALSE.
  END TYPE Special_T
  !
  ! LIST FORM
  TYPE Reaction_T
    CHARACTER(LenType)        :: Type, TypeConstant
    CHARACTER(LenName)        :: Comment
    CHARACTER(LenLine)        :: Line1, Line2, Line3, Line4
    CHARACTER(LenName)        :: Factor
    TYPE(Duct_T),     POINTER :: Educt(:)=>NULL(), Product(:)=>NULL()
    REAL(dp),     ALLOCATABLE :: Constants(:)
    TYPE(Duct_T),     POINTER :: InActEduct(:)=>NULL(), InActProduct(:)=>NULL()
    TYPE(Special_T)           :: Special
    INTEGER                   :: nInActEd=0, nInActPro=0
    TYPE(Reaction_T), POINTER :: Next=>NULL()
  END TYPE Reaction_T
  !
  ! ARRAY FORM
  TYPE ReactionStruct_T
    CHARACTER(LenType)              :: Type,  TypeConstant
    CHARACTER(LenLine)              :: Line1='' , Line2='' , Line3='', Line4=''
    LOGICAL                         :: bR = .FALSE. , brX = .FALSE. 
    CHARACTER(LenName)              :: Factor = ''
    CHARACTER(LenName)              :: Comment = ''
    CHARACTER(2)                    :: direction = ''
    REAL(dp)                        :: SumAqCoef     
    TYPE(Special_T)                 :: Special
    TYPE(Duct_T)  ,     ALLOCATABLE :: Educt(:), Product(:)
    REAL(dp),           ALLOCATABLE :: Constants(:)
    REAL(dp),           ALLOCATABLE :: LowConst(:), HighConst(:), TroeConst(:) ! combustion press dep reactions
    REAL(dp),           ALLOCATABLE :: InActEduct(:), InActProduct(:)
    INTEGER                         :: nInActEd = 0, nInActPro = 0, nActEd = 0, nActPro = 0
    INTEGER                         :: nConst = 0
    INTEGER                         :: HenrySpc = 0
    LOGICAL                         :: TB = .FALSE. , TBextra=.FALSE.
    INTEGER,            ALLOCATABLE :: TBidx(:)
    CHARACTER(LenName), ALLOCATABLE :: TBspc(:)
    REAL(dp),           ALLOCATABLE :: TBalpha(:)
    CHARACTER(LenName), ALLOCATABLE :: InActEductSpc(:), InActProductSpc(:)
  END TYPE ReactionStruct_T
  !
  !
  TYPE ListReaction_T
    TYPE(Reaction_T), POINTER :: Start=>NULL()
    TYPE(Reaction_T), POINTER :: End=>NULL()
    INTEGER :: LenList=0
  END TYPE ListReaction_T
  !
  TYPE Species_T
    CHARACTER(LenName) :: Species=''
    LOGICAL            :: isHenry=.FALSE.
    REAL(dp)     :: Hf=0.0d0, Gf=0.0d0, Cp=0.0d0
  END TYPE Species_T

  

  
  TYPE(Reaction_T), POINTER   :: System
  TYPE(ListReaction_T), SAVE  :: ListRGas, ListRHenry, ListRAqua,        &
  &                              ListRDiss
  !
  TYPE(hash_tbl_sll)          :: ListAqua, ListGas,           &
  &                              ListNonReac, ListAtoms
  TYPE(hash_tbl_sll)          :: ListFamilies
  !
  TYPE(Species_T), ALLOCATABLE, TARGET :: ListAqua2(:), ListGas2(:),     &
  &                               ListNonReac2(:)
  INTEGER :: InputUnit=10
  !
  CHARACTER(37), PARAMETER :: SetSpecies='ABCDEFGHIJKLMNOPQRSTUVWXYZapsc[]()=+*'

  TYPE Element_T
    CHARACTER(5) :: Element=''
  END TYPE Element_T

  TYPE(Element_T) :: Elements(11)=(/       &
  &                    Element_T('(')      &
  &                    ,Element_T(')')     &
  &                    ,Element_T('exp')   &
  &                    ,Element_T('+')     &
  &                    ,Element_T('-')     &
  &                    ,Element_T('*')     &
  &                    ,Element_T('/')     &
  &                    ,Element_T('**')    &
  &                    ,Element_T('abs')   &
  &                    ,Element_T('sqrt')  &
  &                    ,Element_T('log')   /)
 
  INTEGER :: nsr                        ! # activ species + all Reactions
  INTEGER :: nsr2                       ! # activ species + all Reactions
 
  !
  CHARACTER(20) :: Filename
  CHARACTER(20) :: IniName
  !
  TYPE(Reaction_T), POINTER :: Current
  TYPE(ReactionStruct_T), ALLOCATABLE :: ReactionSystem(:)
  TYPE(ListReaction_T), ALLOCATABLE :: CompleteReactionList(:)
  !
  !
  REAL(dp), ALLOCATABLE :: Emis(:)          & ! emission values
  &                      , InitValInAct(:)    ! initial values inactiv spc
  !
  !
  REAL(dp), ALLOCATABLE :: sumBAT(:)          ! sum_j=1,n_s (b_ij-a_ij),  i el. N_R

  INTEGER :: fNumber = 0
  !
  CONTAINS
  ! ------------------------------------
  ! -----------SUBROUTINES--------------
  ! ------------------------------------
  !
  SUBROUTINE SortReactionList(ReacStructOut,ReacStructIn)
    TYPE(ReactionStruct_T), ALLOCATABLE :: ReacStructIn(:)
    TYPE(ReactionStruct_T), ALLOCATABLE :: ReacStructOut(:)
    !
    INTEGER :: i
    CHARACTER(20), ALLOCATABLE :: ReacTypeSorted(:)
    INTEGER, ALLOCATABLE :: iReacTypeSorted(:)
    !
    ! sort the reaction list --> TypeConstant
    ALLOCATE(ReacTypeSorted(nreac))
    ALLOCATE(iReacTypeSorted(nreac))  
    DO i=1,nreac
      ReacTypeSorted(i)=ReacStructIn(i)%Type
    END DO
    CALL StringSort(ReacTypeSorted,iReacTypeSorted)
    ALLOCATE(ReacStructOut(nreac))
    !
    DO i=1,SIZE(ReacStructIN)
      ReacStructOut(i)=ReacStructIn(iReacTypeSorted(i))
    END DO
    DEALLOCATE(ReacStructIn)
    DEALLOCATE(iReacTypeSorted)
  END SUBROUTINE SortReactionList
  !
  !
  SUBROUTINE ReadReaction(Out)
    LOGICAL :: Out
    !
    INTEGER :: iLine,PosColon,Pos,is
    CHARACTER(LenLine) :: LocLine
    CHARACTER(LenLine) :: Line(1:4)
    CHARACTER(20) :: CLASS
    CHARACTER(40) :: TypeR
    INTEGER :: idxFAC
    
    iLine = 0
    Line = ''
 
    DO

      READ( InputUnit , '(A'//TRIM(ADJUSTL(LenLineChar))//')' , IOSTAT=is ) LocLine

      IF ( ABS(is) > 0 ) EXIT

      idxFAC = INDEX(LocLine,'$')
      IF ( idxFAC > 0 ) THEN
        SELECT CASE (TRIM(LocLine(idxFAC:)))
          CASE ('$H2','$O2N2','$M','$O2','$N2','$H2O','$O2O2','$aH2O','$+M','$(+M)','$RO2','$RO2aq')
            IF ( TRIM(LocLine(idxFAC:)) == '$H2'    ) nr_FAC_H2    = nr_FAC_H2    + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$O2N2'  ) nr_FAC_O2N2  = nr_FAC_O2N2  + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$M'     ) nr_FAC_M     = nr_FAC_M     + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$O2'    ) nr_FAC_O2    = nr_FAC_O2    + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$N2'    ) nr_FAC_N2    = nr_FAC_N2    + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$H2O'   ) nr_FAC_H2O   = nr_FAC_H2O   + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$O2O2'  ) nr_FAC_O2O2  = nr_FAC_O2O2  + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$RO2'   ) nr_FAC_RO2   = nr_FAC_RO2   + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$RO2aq' ) nr_FAC_RO2aq = nr_FAC_RO2aq + 1
            IF ( TRIM(LocLine(idxFAC:)) == '$aH2O'  ) nr_FAC_aH2O  = nr_FAC_aH2O  + 1
            nr_FACTOR = nr_FACTOR + 1
          CASE DEFAULT
            WRITE(*,*) '  Unknown FACTOR:  ', TRIM(LocLine(idxFAC:)), '   at Line:  ???'
        END SELECT
      END IF

      ! if no comment or blank line then
      IF ( ADJUSTL(LocLine(1:1)) /= '#'       .AND.  &
      &    ADJUSTL(LocLine(1:7)) /= 'COMMENT' .AND.  &
      &    LEN(TRIM(LocLine)) > 0 ) THEN
        iLine = iLine + 1
        Line(iLine) = LocLine
        IF ( iLine == 4 ) EXIT
      END IF

    END DO

    IF ( iLine >= 3 ) THEN
      Pos = SCAN(Line(1),'#')
      IF ( Pos > 0 ) Line(1) = Line(1)(:Pos-1)
      
      ! if new reaction line starts go one line back
      IF ( INDEX(Line(4),'CLASS') > 0 ) BACKSPACE(InputUnit)
      
      ! read the reaction TYPE
      PosColon = Index(Line(1),':')
      CLASS    = ADJUSTL(Line(1)(PosColon+1:))
      
      ! count the number of each reaction type
      SELECT CASE (CLASS)
        
        CASE ('GAS') ! gaseous phase reactions

          nr_gas = nr_gas + 1
          CALL InsertReaction( ListRGas , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('PHOTO','PHOTO2','PHOTO3','PHOTAB','PHOTABC','PHOTMCM')
              nr_G_photo = nr_G_photo + 1
              IF ( TypeR == 'PHOTAB'  ) nr_PHOTab  = nr_PHOTab  + 1
              IF ( TypeR == 'PHOTABC' ) nr_PHOTabc = nr_PHOTabc + 1
              IF ( TypeR == 'PHOTMCM' ) nr_PHOTmcm = nr_PHOTmcm + 1
              IF ( TypeR == 'PHOTO'   ) nr_PHOTOkpp  = nr_PHOTOkpp  + 1
              IF ( TypeR == 'PHOTO2'  ) nr_PHOTO2kpp = nr_PHOTO2kpp + 1
              IF ( TypeR == 'PHOTO3'  ) nr_PHOTO3kpp = nr_PHOTO3kpp + 1
            CASE ('CONST')
              nr_G_const = nr_G_const + 1
              nr_CONST   = nr_CONST   + 1
            CASE ('TEMP','TEMP1','TEMP2','TEMP3','TEMP4')
              nr_G_temp  = nr_G_temp + 1
              IF ( TypeR == 'TEMP' )  nr_TEMP  = nr_TEMP  + 1
              IF ( TypeR == 'TEMP1' ) nr_TEMP1 = nr_TEMP1 + 1
              IF ( TypeR == 'TEMP2' ) nr_TEMP2 = nr_TEMP2 + 1
              IF ( TypeR == 'TEMP3' ) nr_TEMP3 = nr_TEMP3 + 1
              IF ( TypeR == 'TEMP4' ) nr_TEMP4 = nr_TEMP4 + 1
            CASE ('TROE','TROEF','TROEQ','TROEQF','TROEXP','TROEMCM')
              nr_G_troe = nr_G_troe + 1
              IF ( TypeR == 'TROE'    ) nr_TROE    = nr_TROE    + 1
              IF ( TypeR == 'TROEF'   ) nr_TROEf   = nr_TROEf   + 1
              IF ( TypeR == 'TROEQ'   ) nr_TROEq   = nr_TROEq   + 1
              IF ( TypeR == 'TROEQF'  ) nr_TROEqf  = nr_TROEqf  + 1
              IF ( TypeR == 'TROEXP'  ) nr_TROExp  = nr_TROExp  + 1
              IF ( TypeR == 'TROEMCM' ) nr_TROEmcm = nr_TROEmcm + 1
            CASE ('SPEC1','SPEC2','SPEC3','SPEC4','SPEC1MCM',  &
            &     'SPEC2MCM','SPEC3MCM','SPEC4MCM','SPEC5MCM', &
            &     'SPEC6MCM','SPEC7MCM','SPEC8MCM','SPEC9MCM'  )
              nr_G_spec = nr_G_spec + 1
              IF ( TypeR == 'SPEC1' ) nr_SPEC1 = nr_SPEC1 + 1
              IF ( TypeR == 'SPEC2' ) nr_SPEC2 = nr_SPEC2 + 1
              IF ( TypeR == 'SPEC3' ) nr_SPEC3 = nr_SPEC3 + 1
              IF ( TypeR == 'SPEC4' ) nr_SPEC4 = nr_SPEC4 + 1
              IF ( TypeR == 'SPEC1MCM' ) nr_SPEC1mcm = nr_SPEC1mcm + 1
              IF ( TypeR == 'SPEC2MCM' ) nr_SPEC2mcm = nr_SPEC2mcm + 1
              IF ( TypeR == 'SPEC3MCM' ) nr_SPEC3mcm = nr_SPEC3mcm + 1
              IF ( TypeR == 'SPEC4MCM' ) nr_SPEC4mcm = nr_SPEC4mcm + 1
              IF ( TypeR == 'SPEC5MCM' ) nr_SPEC5mcm = nr_SPEC5mcm + 1
              IF ( TypeR == 'SPEC6MCM' ) nr_SPEC6mcm = nr_SPEC6mcm + 1
              IF ( TypeR == 'SPEC7MCM' ) nr_SPEC7mcm = nr_SPEC7mcm + 1
              IF ( TypeR == 'SPEC8MCM' ) nr_SPEC8mcm = nr_SPEC8mcm + 1
              IF ( TypeR == 'SPEC9MCM' ) nr_SPEC9mcm = nr_SPEC9mcm + 1
            CASE ('S4H2O')
              nr_S4H2O = nr_S4H2O + 1
            CASE ('T1H2O')
              nr_T1H2O = nr_T1H2O + 1
            CASE ('SPECIAL')
              nr_SPECIAL = nr_SPECIAL + 1
              nr_G_special = nr_G_special + 1
            CASE ('HOM1')
              nr_HOM1 = nr_HOM1 + 1
            CASE DEFAULT
              WRITE(*,*) '  Unknown gaseous reaction: ', TypeR
          END SELECT

        CASE ('HENRY')        ! phase transfer pseudo-reactions 

          nr_henry = nr_henry + 1
          CALL InsertReaction( ListRHenry , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('TEMP','TEMP1','TEMP2','TEMP3','TEMP4')
              IF ( TypeR == 'TEMP' )  nr_TEMP  = nr_TEMP  + 1
              IF ( TypeR == 'TEMP1' ) nr_TEMP1 = nr_TEMP1 + 1
              IF ( TypeR == 'TEMP2' ) nr_TEMP2 = nr_TEMP2 + 1
              IF ( TypeR == 'TEMP3' ) nr_TEMP3 = nr_TEMP3 + 1
              IF ( TypeR == 'TEMP4' ) nr_TEMP4 = nr_TEMP4 + 1
            CASE ('CONST')
              nr_CONST = nr_CONST + 1
            CASE ('SPECIAL')
              nr_SPECIAL = nr_SPECIAL + 1
              nr_H_special = nr_H_special + 1
            CASE DEFAULT
              WRITE(*,*) '  Unknown phase transfer reaction: ', TypeR
          END SELECT

        CASE ('AQUA')         ! aquatic phase reactions 

          nr_aqua = nr_aqua + 1
          CALL InsertReaction( ListRAqua , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('PHOTO','PHOTO2','PHOTO3','PHOTAB','PHOTABC','PHOTMCM')
              nr_A_photo = nr_A_photo + 1
              IF ( TypeR == 'PHOTAB'  ) nr_PHOTab  = nr_PHOTab  + 1
              IF ( TypeR == 'PHOTABC' ) nr_PHOTabc = nr_PHOTabc + 1
              IF ( TypeR == 'PHOTMCM' ) nr_PHOTmcm = nr_PHOTmcm + 1
              IF ( TypeR == 'PHOTO'   ) nr_PHOTOkpp  = nr_PHOTOkpp  + 1
              IF ( TypeR == 'PHOTO2'  ) nr_PHOTO2kpp = nr_PHOTO2kpp + 1
              IF ( TypeR == 'PHOTO3'  ) nr_PHOTO3kpp = nr_PHOTO3kpp + 1
            CASE ('CONST')
              nr_A_const = nr_A_const + 1
              nr_CONST   = nr_CONST   + 1
            CASE ('TEMP','Temp1''TEMP2','TEMP3','TEMP4')
              nr_A_temp  = nr_A_temp + 1
              IF ( TypeR == 'TEMP' )  nr_TEMP  = nr_TEMP  + 1
              IF ( TypeR == 'TEMP1' ) nr_TEMP1 = nr_TEMP1 + 1
              IF ( TypeR == 'TEMP2' ) nr_TEMP2 = nr_TEMP2 + 1
              IF ( TypeR == 'TEMP3' ) nr_TEMP3 = nr_TEMP3 + 1
              IF ( TypeR == 'TEMP4' ) nr_TEMP4 = nr_TEMP4 + 1
            CASE ('ASPEC1','ASPEC2','ASPEC3')
              nr_A_spec  = nr_A_spec + 1
              IF ( TypeR == 'ASPEC1' ) nr_ASPEC1 = nr_ASPEC1 + 1
              IF ( TypeR == 'ASPEC2' ) nr_ASPEC2 = nr_ASPEC2 + 1
              IF ( TypeR == 'ASPEC3' ) nr_ASPEC3 = nr_ASPEC3 + 1
            CASE ('SPECIAL')
              nr_SPECIAL = nr_SPECIAL + 1
              nr_A_special = nr_A_special + 1
            CASE ('HOM1')
              nr_HOM1 = nr_HOM1 + 1
            CASE DEFAULT
              WRITE(*,*) '  Unknown aqueous reaction: ', TypeR
          END SELECT

        CASE ('DISS')        ! fast aquatic phase equil. reactions 

          nr_diss = nr_diss + 1
          CALL InsertReaction( ListRDiss , Line , TypeR )

          SELECT CASE (TypeR)
            CASE ('DCONST','DTEMP','DTEMP2','DTEMP3','DTEMP4','DTEMP5')
              IF ( TypeR == 'DCONST'    ) nr_DCONST = nr_DCONST + 1
              IF ( TypeR == 'DTEMP'     ) nr_DTEMP  = nr_DTEMP  + 1
              IF ( TypeR == 'DTEMP2'    ) nr_DTEMP2 = nr_DTEMP2 + 1
              IF ( TypeR == 'DTEMP3'    ) nr_DTEMP3 = nr_DTEMP3 + 1
              IF ( TypeR == 'DTEMP4'    ) nr_DTEMP4 = nr_DTEMP4 + 1
              IF ( TypeR == 'DTEMP5'    ) nr_DTEMP5 = nr_DTEMP5 + 1
            CASE ('SPECIAL')
              nr_SPECIAL = nr_SPECIAL + 1
              nr_D_special = nr_D_special + 1
            CASE DEFAULT
              WRITE(*,*) '  Unknown dissociation reaction: ', TypeR
          END SELECT

        CASE DEFAULT
          WRITE(*,*) '  Unknown reaction CLASS: ', CLASS
          STOP
      END SELECT

      Out = .FALSE.
    ELSE
      Out = .TRUE.
    END IF

  END SUBROUTINE ReadReaction
  !
  !
  SUBROUTINE CompressParty(Ducts,Perm,Len)
    INTEGER, ALLOCATABLE :: Ducts(:)
    INTEGER, ALLOCATABLE :: Perm(:)
    INTEGER :: Len
    !
    ! sort ColInd and Val for acc column indx
    Perm=0
    Len=0
    CALL unirnk(Ducts,Perm,Len)
    Ducts=Ducts(Perm)
    CALL CompressList(Ducts)
  END SUBROUTINE CompressParty
  !
  !
  SUBROUTINE PrintSpecies(ListName,Unit)
    TYPE(Species_T) :: ListName(:)
    INTEGER :: Unit
    !
    INTEGER :: i

    DO i=1,SIZE(ListName)
      WRITE(Unit,*) "'"//TRIM(ListName(i)%Species)//"'"
    END DO
  END SUBROUTINE PrintSpecies
  !
  !
  SUBROUTINE SpcIdx(ListName,idx)
    TYPE(Species_T) :: ListName(:)
    INTEGER :: idx
    !
    INTEGER :: i
    !
    DO i=1,SIZE(ListName)
      IF (i==idx) THEN
        WRITE(*,*) "'"//TRIM(ListName(idx)%Species)//"'"
      END IF
    END DO
  END SUBROUTINE SpcIdx
  !
  !
  SUBROUTINE PrintHeadSpecies(Filename,Unit)
    INTEGER :: Unit
    CHARACTER(*) :: Filename
    !
    CHARACTER(8) :: Date
    CHARACTER(10) :: Time
    INTEGER(8) :: Value(8)
    !
    CALL DATE_AND_TIME(Date,Time,VALUES=Value)
    !
    WRITE(Unit,*) ' ==========================================================='
    WRITE(Unit,*) ' ========  0-dim Simulation of chemical mechanisms  ========'
    WRITE(Unit,*) ' ========     Output -  Chemical Reaction Data      ========'
    WRITE(Unit,*) ' ==========================================================='
    WRITE(Unit,*) ''
    WRITE(Unit,*) ' Created:             ', Date(7:8),'.',Date(5:6),'.',Date(1:4)
    WRITE(Unit,*) ' Chemical Mechanism:  ', TRIM(ADJUSTL(FileName))
    WRITE(Unit,*) ''
    WRITE(Unit,*) ' =================     Units         ======================='
    WRITE(Unit,*) ''
    IF (UnitGas==0) THEN
      WRITE(Unit,*) ' Gas Phase Units:     molec/cm3'
    ELSE
      WRITE(Unit,*) ' Gas Phase Units:     mol/m3'
    END IF
    IF (UnitAqua==0) THEN
      WRITE(Unit,*) ' Aqueous Phase Units: mol/l'
    END IF
    WRITE(Unit,*) ''
    WRITE(Unit,*) ' =================    Numbers        ======================='
    WRITE(Unit,*) ''
    WRITE(Unit,*) ns_GAS &
                 +ns_AQUA &
                 +ns_KAT,    '     Number of Species'
    WRITE(Unit,*) ns_GAS,    '     No. of gaseous species'
    WRITE(Unit,*) ns_AQUA,   '     No. of aqueous species'
    WRITE(Unit,*) ns_KAT,'     Number of Non-reactive Species '
    WRITE(Unit,*) ''
    WRITE(Unit,*) ' =================   Species Names   ======================='
    WRITE(Unit,*) ''
  END SUBROUTINE PrintHeadSpecies
  !
  !
  SUBROUTINE PrintFinalReactions(Unit)
    INTEGER :: Unit
    !
    WRITE(Unit,*) ''
    WRITE(Unit,*) ''
    WRITE(Unit,*) '========================================================='
    WRITE(Unit,*) '========             End of .chem File           ========'
    WRITE(Unit,*) '========================================================='
  END SUBROUTINE PrintFinalReactions
  !
  !
  SUBROUTINE PrintHeadReactions(Unit)
    INTEGER :: Unit
  
    nreac = nr_gas  + 2*nr_henry + nr_aqua   &
    &     + 2*nr_diss

    nr_D_Temp = nr_DTEMP  + nr_DTEMP2 &
    &         + nr_DTEMP3 + nr_DTEMP4 + nr_DTEMP5

    WRITE(Unit,*) ''
    WRITE(Unit,*) ' ================   Description of Reactions   =============='
    WRITE(Unit,*) ''
    WRITE(Unit,*) nreac,       '        NREAK  : Number of Reactions'
    WRITE(Unit,*) nr_gas,      '        NGAS   : Gas phase reactions'
    WRITE(Unit,*) nr_G_photo,  '           Gaseous PHOTO - type reactions'
    WRITE(Unit,*) nr_G_const,  '           Gaseous CONST - type reactions'
    WRITE(Unit,*) nr_G_temp,   '           Gaseous TEMP - type reactions'
    WRITE(Unit,*) nr_SimpTB,   '           Gaseous Simple three-body - type reactions'
    WRITE(Unit,*) nr_G_lind,   '           Gaseous Lindemann - type reactions'
    WRITE(Unit,*) nr_G_troe,   '           Gaseous TROE - type reactions'
    WRITE(Unit,*) nr_G_spec,   '           Gaseous SPEC - type reactions'
    WRITE(Unit,*) nr_G_special,'           Gaseous SPECIAL formula - type reactions'
    WRITE(Unit,*) nr_henry,    '        NHENRY : Henry Equilib. reactions'
    WRITE(Unit,*) nr_diss,     '        NDISS  : Dissociation reactions'
    WRITE(Unit,*) nr_DCONST,   '           Aqueous DCONST - type reactions'
    WRITE(Unit,*) nr_D_TEMP,   '           Aqueous DTEMP - type reactions'
    WRITE(Unit,*) nr_D_special,'           Aqueous SPECIAL formula reactions'
    WRITE(Unit,*) nr_aqua,     '        NAQUA  : Aquatic Equilib. reactions'
    WRITE(Unit,*) nr_A_photo,  '           Aqueous PHOTO - type reactions'
    WRITE(Unit,*) nr_A_const,  '           Aqueous CONST - type reactions'
    WRITE(Unit,*) nr_A_temp,   '           Aqueous TEMP - type reactions'
    WRITE(Unit,*) nr_A_spec,   '           Aqueous SPEC - type reactions'
    WRITE(Unit,*) nr_A_special,'           Aqueous SPECIAL formula - type reactions'
    WRITE(Unit,*)
    WRITE(Unit,*) ' ======================  Reactions   ========================'
    WRITE(Unit,*) ''
  END SUBROUTINE PrintHeadReactions
  !
  SUBROUTINE Print_SysFile(RS,IndexSet,NewName)
    TYPE(ReactionStruct_T), INTENT(IN) :: RS(:)
    INTEGER,      OPTIONAL, INTENT(IN) :: IndexSet(:)
    CHARACTER(*), OPTIONAL             :: NewName
    
    INTEGER :: iR, j
    INTEGER :: nR
    
    OPEN(UNIT=989,FILE=ADJUSTL(TRIM(NewName)),STATUS='UNKNOWN')
  
    WRITE(989,'(A)') '# ================= '//TRIM(NewName)//' ================='
    WRITE(989,'(A)') '# = Please copy the data into your sys-file for ='
    WRITE(989,'(A)') '# =============== chemical input. ==============='
    WRITE(989,'(A)') '#'
    WRITE(989,'(A)') '#  ===================   Unit options   ======================'
    WRITE(989,'(A)') ''
    WRITE(989,'(A)') 'UNIT GAS    0   #    Gas phase units     (0 = molec/cm3, 1 = mol/m3)'
    WRITE(989,'(A)') 'UNIT AQUA   0   #    Aqueous phase units (0 = mol/l)'
    WRITE(989,'(A)') ''
    WRITE(989,'(A)') '#'

 
    IF (PRESENT(IndexSet)) THEN
      nR = SIZE(IndexSet)
      j = 0

      ! Print IndexSet of reactions
      DO 
        j = j + 1
        iR = IndexSet(j)

        IF ( .NOT. RS(iR)%bR ) THEN
          ! if direction is --> then
          CALL PrintReac(iR)
        ELSE
          ! for equilibrium reactions we need to check if the correct one is printed
          ! by correct one i mean the reaction which occurs in the original system
          ! not the reversed reaction string, because the evaluation of the rate 
          ! constant would be wrong if the reversed reaction string is printed
          IF ( iR-1 /= IndexSet(j-1) ) THEN
            CALL PrintReac(iR-1)
          END IF
        END IF

        IF ( j >= nR ) EXIT
      END DO
    ELSE
      nR = SIZE(RS)

      ! Print all reaction in RS list
      iR = 0
      DO 
        iR = iR + 1
        CALL PrintReac(iR)
        IF ( iR >= nR ) EXIT
        IF ( MAXVAL(INDEX(RS(iR)%Type,['DISS','HENR'])) > 0 ) iR = iR + 1
      END DO
    
    END IF
    CLOSE(989)

      CONTAINS
        
      SUBROUTINE PrintReac(iR)
        INTEGER :: iR
        WRITE(989,'(A)') 'CLASS: '//TRIM(RS(iR)%Type)
        WRITE(989,'(A)') TRIM(RS(iR)%Line1)
        IF ( TRIM(RS(iR)%Line3) /= '' ) WRITE(989,'(A)') TRIM(RS(iR)%Line3)
        !IF ( TRIM(RS(i)%Special%Formula) /= '' ) WRITE(989,'(A,I0)') 'SPECIAL: '//TRIM(RS(i)%Special%Formula)//';  ',RS(i)%Special%nVariables
        IF ( TRIM(RS(iR)%Factor)  /= '' ) WRITE(989,'(A)') 'FACTOR: '//TRIM(RS(iR)%Factor)
        WRITE(989,'(A)') 
      END SUBROUTINE PrintReac
  END SUBROUTINE Print_SysFile

  FUNCTION RemoveSpaces(String) RESULT(StringOut)
    ! replaces multiple spaces in string by one space
    CHARACTER(*), INTENT(IN) :: String
    CHARACTER(LEN(String))   :: StringOUT

    INTEGER :: i

    StringOut = TRIM(String)

    DO
      i = INDEX(TRIM(StringOut),'  ')
      IF ( i == 0 ) EXIT
      StringOut = TRIM(ADJUSTL(StringOut(:i-1))) &
      &           //' '//                        &
      &           TRIM(ADJUSTL(StringOut(i+1:)))
    END DO

  END FUNCTION RemoveSpaces

  !
  SUBROUTINE Print_ChemFile(RS,File,Unit,CK)
    ! IN:
    TYPE(ReactionStruct_T), ALLOCATABLE :: RS(:)
    CHARACTER(*) :: File
    INTEGER      :: Unit
    LOGICAL      :: CK
    ! TEMP:
    INTEGER      :: io_stat
    INTEGER      :: i,j,m,iR
    INTEGER      :: nEduct,nProd
    TYPE(Duct_T), ALLOCATABLE :: ActiveEduct(:)
    TYPE(Duct_T), ALLOCATABLE :: ActiveProduct(:)
    !
    INTEGER      :: nnzA, nnzB
    !
    INTEGER,  ALLOCATABLE :: tmpIdx(:)
    REAL(dp), ALLOCATABLE :: tmpVal(:)
    INTEGER,  ALLOCATABLE :: permutation(:)
    INTEGER :: newLen


     
    !-----------------------------------------------------------------------
    ! --- Build the reaction system
    !-----------------------------------------------------------------------
    IF ( .NOT.CK ) THEN
      CALL AllListsToArray( RS            &
      &                   , ListRGas    , ListRHenry  &
      &                   , ListRAqua   , ListRDiss   )
    END IF

    OPEN(unit=Unit, file=File, status='replace', action='write', access='sequential', iostat=io_stat)
    IF ( io_stat /= 0 ) WRITE(*,*) '  ERROR creating chem-file :: ',io_stat
    REWIND(ChemUnit)
    !----------------------------------------------------------------
    ! ---  build the coeficient matrices and write .chem
    CALL PrintHeadSpecies ( File , Unit )
    
    IF ( ns_GAS   > 0 ) CALL PrintSpecies( ListGas2     , Unit )
    IF ( ns_AQUA  > 0 ) CALL PrintSpecies( ListAqua2    , Unit )
    IF ( ns_KAT   > 0 ) CALL PrintSpecies( ListNonReac2 , Unit )
     
    CALL PrintHeadReactions( Unit )


    !-----------------------------------------------------------------------
    ! --- print reactions and build A, B and (B-A) structure
    !-----------------------------------------------------------------------
   
    ! set matrix dimensions
    A%m  = nreac;   A%n  = nspc
    B%m  = nreac;   B%n  = nspc
    BA%m = nreac;   BA%n = nspc
    
    ! Standart alloc
    ALLOCATE(A%RowPtr(A%m+1),B%RowPtr(B%m+1),BA%RowPtr(BA%m+1))
    A%RowPtr  = 0;    A%RowPtr(1)  = 1
    B%RowPtr  = 0;    B%RowPtr(1)  = 1
    BA%RowPtr = 0;    BA%RowPtr(1) = 1

    ! find maximal duct number and allocate
    nEduct=0
    nProd=0
    DO iR=1,nreac
      IF (SIZE(RS(iR)%Educt)>nEduct)  nEduct = SIZE(RS(iR)%Educt)
      IF (SIZE(RS(iR)%Product)>nProd) nProd  = SIZE(RS(iR)%Product)
    END DO
    ALLOCATE(ActiveEduct(nEduct),ActiveProduct(nProd))
    
    DO iR=1,nreac
      ! count activ educts in reaction iR
      nEduct = 0
      DO i=1,SIZE(RS(iR)%Educt)
        SELECT CASE(RS(iR)%Educt(i)%Type)
          CASE ('Gas','Aqua','GAS')
            nEduct = nEduct + 1
            ActiveEduct(nEduct) = RS(iR)%Educt(i)
          CASE ('Inert', 'Dummy')
            ! inactive
          CASE DEFAULT
            WRITE(*,*) 'Unknown Duct Type in Print_ChemFile! Info: '
            WRITE(*,*) '  Type ', TRIM(RS(iR)%Educt(i)%Type)
            WRITE(*,*) '  Species ', TRIM(RS(iR)%Educt(i)%Species)
            WRITE(*,*) '  Reaction Number ', iR
            WRITE(*,*) '  Reaction Lines '
            WRITE(*,*) TRIM(RS(iR)%Line1)
            WRITE(*,*) TRIM(RS(iR)%Line2)
            WRITE(*,*) TRIM(RS(iR)%Line3)
            WRITE(*,*) TRIM(RS(iR)%Line4)
            STOP
        END SELECT
      END DO
      ! count activ products in reaction iR
      nProd = 0
      DO i=1,SIZE(RS(iR)%Product)
        SELECT CASE(RS(iR)%Product(i)%Type)
          CASE ('Gas','Aqua','GAS')
            nProd = nProd + 1
            ActiveProduct(nProd) = RS(iR)%Product(i)
          CASE ('Inert', 'Dummy')
            ! inactive
          CASE DEFAULT
            WRITE(*,*) 'Unknown Duct Type in Print_ChemFile! Info: '
            WRITE(*,*) '  Type ', TRIM(RS(iR)%Product(i)%Type)
            WRITE(*,*) '  Species ', TRIM(RS(iR)%Product(i)%Species)
            WRITE(*,*) '  Reaction Number ', iR
            WRITE(*,*) '  Reaction Lines '
            WRITE(*,*) TRIM(RS(iR)%Line1)
            WRITE(*,*) TRIM(RS(iR)%Line2)
            WRITE(*,*) TRIM(RS(iR)%Line3)
            WRITE(*,*) TRIM(RS(iR)%Line4)
            STOP
        END SELECT
      END DO
   
      !iR = iR + 1
      WRITE(Unit,*)
      WRITE(Unit,'(A,I6,A)') '#-----------', iR ,'. Reaction ----------- '
    
      WRITE(Unit,*) TRIM(RS(iR)%Type)//'   '//TRIM(RS(iR)%TypeConstant)
     
      WRITE(Unit,*) SIZE(RS(iR)%Educt), SIZE(RS(iR)%Product), nEduct, nProd
      
      ! Educt Matrix A
      IF ( nEduct > 1 ) THEN
        ALLOCATE( tmpIdx(nEduct), tmpVal(nEduct) )
        tmpIdx = 0;  tmpVal = ZERO
        DO j=1,nEduct
          tmpIdx(j) = PositionSpeciesAll(ActiveEduct(j)%Species)
          tmpVal(j) = ActiveEduct(j)%Koeff
        END DO
        CALL CompressList(tmpIdx,tmpVal)
        A%RowPtr(iR+1) = A%RowPtr(iR) + SIZE(tmpIdx)
        DEALLOCATE( tmpIdx , tmpVal )
      ELSE
        A%RowPtr(iR+1) = A%RowPtr(iR) + nEduct
      END IF

      ! Product Matrix B
      IF (nProd>1) THEN
        ALLOCATE( tmpIdx(nProd), tmpVal(nProd) )
        tmpIdx = 0;  tmpVal = ZERO
        DO j=1,nProd
          tmpIdx(j) = PositionSpeciesAll(ActiveProduct(j)%Species)
          tmpVal(j) = ActiveProduct(j)%Koeff
        END DO
        CALL CompressList(tmpIdx,tmpVal)
        B%RowPtr(iR+1) = B%RowPtr(iR) + SIZE(tmpIdx)
        DEALLOCATE( tmpIdx , tmpVal )
      ELSE
        B%RowPtr(iR+1) = B%RowPtr(iR) + nProd
      END IF

      ! ----------------------------------------------------
      ! SpeciesIndx Educt=> 1:#Educt of Reaction iR
      ! SpeciesIndx Product=> #Educt+1:#Educt+#Product of Reaction iR
      ! #active educts of Reaction
      WRITE(Unit,*) (PositionSpeciesAll(RS(iR)%Educt(i)%Species),  &
      &             i=1,SIZE(RS(iR)%Educt)),                       &
      &             (PositionSpeciesAll(RS(iR)%Product(i)%Species),&
      &             i=1,SIZE(RS(iR)%Product)),                     &
      &             nEduct+nProd
      !
      !----------------------------------------------------
      ! Tuple: (SpeciesIndex,-Coeffs) for all active Educts (left)
      ! Tuple: (SpeciesIndex,+Coeffs) for all active Products (right)
      WRITE(Unit,'(*(7X,I5,3X,F6.3))', ADVANCE='NO')                     &
      &                   ( PositionSpeciesAll(ActiveEduct(i)%Species)   &
      &                  ,  -ActiveEduct(i)%Koeff,i=1,nEduct )           &
      &                  ,( PositionSpeciesAll(ActiveProduct(i)%Species) &
      &                  ,   ActiveProduct(i)%Koeff,i=1,nProd )
      WRITE(Unit,*)
      !
      IF (RS(iR)%TypeConstant=='SPECIAL') THEN
        WRITE(Unit,*) TRIM(RS(iR)%Line3)
      ELSE
        WRITE(Unit,*) SIZE(RS(iR)%Constants), RS(iR)%Constants
      END IF

      IF (RS(iR)%Factor /= '') WRITE(Unit,*) 'FACTOR:  ',RS(iR)%Factor

      SELECT CASE (RS(iR)%Factor)
        CASE ('$RO2');   hasRO2   = .TRUE.
        CASE ('$RO2aq'); hasRO2aq = .TRUE.
      END SELECT
      
      IF (CK) WRITE(Unit,*) 'EXTRA1:  ',ADJUSTL(TRIM(RS(iR)%Line2))
      IF (CK) WRITE(Unit,*) 'EXTRA2:  ',ADJUSTL(TRIM(RS(iR)%Line3))
    END DO
   
    ! loop again to set ColInd and Val on A and B
    nnzA = 0
    nnzB = 0
  
    ALLOCATE( A%ColInd(A%RowPtr(A%m+1)-1) , A%Val(A%RowPtr(A%m+1)-1) &
    &       , B%ColInd(B%RowPtr(B%m+1)-1) , B%Val(B%RowPtr(B%m+1)-1) )
    A%ColInd = 0; A%Val = ZERO
    B%ColInd = 0; B%Val = ZERO
    
    ALLOCATE(sumBAT(A%m)); sumBAT = ZERO
    !
    DO iR = 1,nreac
      nEduct = 0
      DO i=1,SIZE(RS(iR)%Educt)
        SELECT CASE(RS(iR)%Educt(i)%Type)
          CASE ('Gas','Aqua','GAS')
            nEduct = nEduct + 1
            ActiveEduct(nEduct) = RS(iR)%Educt(i)
        END SELECT
      END DO
      nProd = 0
      DO i=1,SIZE(RS(iR)%Product)
        SELECT CASE(RS(iR)%Product(i)%Type)
          CASE ('Gas','Aqua','GAS')
            nProd = nProd + 1
            ActiveProduct(nProd) = RS(iR)%Product(i)
        END SELECT
      END DO
      
      ! set ColInd and Val on A and B
      IF (nEduct>1) THEN
        ALLOCATE(tmpIdx(nEduct),tmpVal(nEduct),permutation(nEduct))
        tmpIdx = 0;  tmpVal = ZERO;   permutation = 0
        DO j=1,nEduct
          tmpIdx(j) = PositionSpeciesAll(ActiveEduct(j)%Species)
          tmpVal(j) = ActiveEduct(j)%Koeff
        END DO
        
        ! sort ColInd and Val for acc column indx
        CALL unirnk(tmpIdx,permutation,newLen)
        tmpIdx = tmpIdx(permutation); tmpVal = tmpVal(permutation)
        CALL CompressList(tmpIdx,tmpVal)
        !

        DO m=1,newLen
          nnzA = nnzA + 1
          A%ColInd(nnzA) = tmpIdx(m)
          A%Val(nnzA) = tmpVal(m)
          sumBAT(iR)  = sumBAT(iR) - tmpVal(m)
        END DO
        DEALLOCATE(tmpIdx,tmpVal,permutation)
      ELSE
        ! reactions with only one educt
        DO m=1,nEduct
          nnzA = nnzA + 1
          A%ColInd(nnzA) = PositionSpeciesAll(ActiveEduct(m)%Species)
          A%Val(nnzA) = ActiveEduct(m)%Koeff
          sumBAT(iR)  = sumBAT(iR) - ActiveEduct(m)%Koeff
        END DO
      END IF
      !
      IF (nProd>1) THEN
        ALLOCATE(tmpIdx(nProd),tmpVal(nProd),permutation(nProd))
        tmpIdx = 0;  tmpVal = ZERO;   permutation = 0
        DO j=1,nProd
          tmpIdx(j) = PositionSpeciesAll(ActiveProduct(j)%Species)
          tmpVal(j) = ActiveProduct(j)%Koeff
        END DO
        CALL unirnk(tmpIdx,permutation,newLen)
        tmpIdx = tmpIdx(permutation); tmpVal = tmpVal(permutation)
        CALL CompressList(tmpIdx,tmpVal)
        !
        DO m=1,newLen
          nnzB = nnzB+1
          B%ColInd(nnzB) = tmpIdx(m)
          B%Val(nnzB) = tmpVal(m)
          sumBAT(iR)  = sumBAT(iR) + tmpVal(m)
        END DO
        DEALLOCATE(tmpIdx,tmpVal,permutation)
      ELSE
        DO m=1,nProd
          nnzB = nnzB + 1
          B%ColInd(nnzB) = PositionSpeciesAll(ActiveProduct(m)%Species)
          B%Val(nnzB) = ActiveProduct(m)%Koeff
          sumBAT(iR)  = sumBAT(iR) + ActiveProduct(m)%Koeff
        END DO
      END IF
    END DO

    A%nnz = nnzA
    B%nnz = nnzB

    CALL PrintFinalReactions( Unit )
    CLOSE(ChemUnit)
  END SUBROUTINE Print_ChemFile

  SUBROUTINE Setup_SpeciesOrder(A)
    
    TYPE(CSR_Matrix_T), INTENT(IN) :: A
    INTEGER :: iR, jj
    INTEGER :: nnz, cnt
    INTEGER, ALLOCATABLE :: tmpFO1(:), tmpFO2(:)
    INTEGER, ALLOCATABLE :: tmpSO1(:), tmpSO2(:)
    INTEGER, ALLOCATABLE :: tmpHO1(:), tmpHO2(:)
    REAL(dp), ALLOCATABLE :: atmpHO(:)
    REAL(dp), PARAMETER :: big = -99999999999999._dp

    nnz = A%nnz

    ALLOCATE( tmpFO1(nnz), tmpFO2(nnz), &
            & tmpSO1(nnz), tmpSO2(nnz), &
            & tmpHO1(nnz), tmpHO2(nnz), &
            & atmpHO(nnz)               )

    tmpFO1  = 0;      tmpFO2  = 0
    tmpSO1  = 0;      tmpSO2  = 0
    tmpHO1  = 0;      tmpHO2  = 0
    atmpHO  = big

    cnt = 0
    DO iR = 1 , nreac
      DO jj = A%RowPtr(iR) , A%RowPtr(iR+1)-1
        cnt = cnt + 1
        IF      (A%Val(jj) == ONE) THEN
          tmpFO1(cnt) = iR
          tmpFO2(cnt) = A%ColInd(jj)
        ELSE IF (A%Val(jj) == TWO) THEN
          tmpSO1(cnt) = iR
          tmpSO2(cnt) = A%ColInd(jj)
        ELSE
          tmpHO1(cnt) = iR
          tmpHO2(cnt) = A%ColInd(jj)
          atmpHO(cnt) = A%Val(jj)
        END IF
      END DO
    END DO

    CALL CompressIntegerArray(tmpFO1); CALL CompressIntegerArray(tmpFO2)
    CALL CompressIntegerArray(tmpSO1); CALL CompressIntegerArray(tmpSO2)
    CALL CompressIntegerArray(tmpHO1); CALL CompressIntegerArray(tmpHO2)
    CALL CompressDoubleArray(atmpHO)

    CALL Save_SpeciesOrder(iFO, tmpFO1, tmpFO2)
    CALL Save_SpeciesOrder(iSO, tmpSO1, tmpSO2)
    CALL Save_SpeciesOrder(iHO, tmpHO1, tmpHO2, atmpHO, aHO)

    nFirst_order  = SIZE(iFO, DIM=1)
    nSecond_order = SIZE(iSO, DIM=1)
    nHigher_order = SIZE(iHO, DIM=1)

    CONTAINS
     SUBROUTINE Save_SpeciesOrder(arr, tmp1, tmp2, a_tmp, a)
       INTEGER, ALLOCATABLE :: arr(:,:)
       INTEGER :: tmp1(:), tmp2(:)
       REAL(dp), OPTIONAL :: a_tmp(:)
       REAL(dp), ALLOCATABLE, OPTIONAL :: a(:)
       
       INTEGER :: cnt, i, j

       ! counting number of reactions taking into account vector-valued (aqueous) reactions
       cnt = 0
       DO i = 1 , SIZE(tmp1)
         IF ( nD_reac(tmp1(i)) ) THEN
           cnt = cnt + nDropletClasses
         ELSE
           cnt = cnt + 1
         END IF
       END DO
       ALLOCATE(arr(cnt, 2))
       IF (PRESENT(a_tmp)) ALLOCATE(a(cnt))

       cnt = 1
       DO i = 1 , SIZE(tmp1)
         DO j = 1 , nD_Ptr_reacs(tmp1(i)+1) - nD_Ptr_reacs(tmp1(i))
           arr(cnt, 1) = nD_Ptr_reacs(tmp1(i)) + j-1
           IF (nD_Spc(tmp2(i))) THEN
             arr(cnt, 2) = nD_Ptr_spc(tmp2(i)) + j-1
           ELSE
             arr(cnt, 2) = nD_Ptr_spc(tmp2(i))
           END IF
           IF (PRESENT(a_tmp)) a(cnt) = a_tmp(i)
           cnt = cnt+1
         END DO
       END DO
     
   END SUBROUTINE Save_SpeciesOrder

  END SUBROUTINE Setup_SpeciesOrder


  SUBROUTINE InputChemicalData(InitFileName,DataFileName)
    CHARACTER(*) :: InitFileName, DataFileName

    INTEGER :: iSpc, i, lb, ub, j
    INTEGER :: io_stat
    REAL(dp) :: LWCs(nDropletClasses)

    REAL(dp), ALLOCATABLE :: InitValKat_1nD(:), y_emi_1nD(:), y_depos_1nD(:), dummy_passiv(:)
    !
    ! for pH Start
    REAL(dp) :: kappa

    !
    ALLOCATE( InitValAct(nspc2), dummy_passiv(ns_KAT), InitValKat_1nD(ns_KAT), InitValKat_Ref(ns_KAT2), InitValKat(ns_KAT2),  &
            & y_emi_1nD(nspc), y_depos_1nD(nspc),  y_emi(nspc2) , y_depos(nspc2), gaseous_passive_ind(0) )
    ! this is for mass transfer (accom , diffus term)
    ALLOCATE( henry_diff(nspc+ns_KAT), henry_accom(nspc+ns_KAT) )

    InitValAct     = 1.0e-20_dp
    InitValKat_1nD = 1.0e-20_dp
    y_emi_1nD      = ZERO; y_emi   = ZERO
    y_depos_1nD    = ZERO; y_depos = ZERO
    henry_diff     = ZERO; henry_diff(1:ns_GAS)  = 5.0e-6_dp
    henry_accom    = ZERO; henry_accom(1:ns_GAS) = 5.0e-5_dp

    !=== Set  Chemical DATA  
    !===  (Molar Mass, Charges)
    ALLOCATE( Charge(ns_AQUA)   &      ! charge of ions
    &       , MolMass(ns_GAS+ns_AQUA)) ! molar mass of species

    !---  Set default values
    Charge   = ZERO
    MolMass  = ZERO

    !
    !--- Read thermodynamic data,....
    CALL Read_SpeciesData( henry_diff, henry_accom, DataFileName )
    !
    !
    !=========================================================================
    !===  Read and save species names in lists
    !=========================================================================
    !
    !-- Open .chem-file and skip the first 22 lines (head)
    OPEN(unit=ChemUnit, file=ChemFile, status='old', action='read', access='sequential', iostat=io_stat)
    IF ( io_stat /= 0 ) WRITE(*,*) '  ERROR opening chem-file :: ',io_stat
    REWIND(ChemUnit)

    DO i=1,22;  READ(ChemUnit,*);  END DO ! skip the first 22 lines

    !---  set indices for pH and water dissociation 
    Hp_ind    = 0
    OHm_ind   = 0
    aH2O_ind  = 0
    H2O_ind   = 0
    H2_ind    = 0
    O2_ind    = 0
    N2_ind    = 0
    Temp_ind  = nspc + 1
    
    !=========================================================================
    !===  Read and Split Species Names and Initital Values
    !=========================================================================
    
    ALLOCATE(y_name(nspc+ns_KAT))

    ! read gaseous phase species
    IF (ns_GAS>0) THEN
      bGs(1) = 1
      bGs(2) = ns_GAS
      iGs = [(i, i=bGs(1),bGs(2))]
      bGs2(1) = 1
      bGs2(2) = ns_GAS
      iGs2 = [(i, i=bGs2(1),bGs2(2))]

      DO iSpc = bGs(1),bGs(2)
        READ(ChemUnit,*,IOSTAT=io_stat)  y_name(iSpc)
      END DO

      InitValAct(iGs) = 1.e-20_dp
      CALL Read_INI_file( InitFileName , InitValAct, InitValKat_1nD , 'GAS' , 'INITIAL' )

      N2O2    = N2 + O2  ! passive species N2+O2

      !---  Read gas phase emission values
      CALL Read_INI_file( InitFileName, y_emi_1nD  , dummy_passiv, 'GAS', 'EMISS' )
      CALL Read_INI_file( InitFileName, y_depos_1nD, dummy_passiv, 'GAS', 'DEPOS' )

    END IF

    ! read aqueous phase species
    IF (ns_AQUA>0) THEN
      CALL Read_Modes( Mode , InitFileName )
      CALL Read_AFrac( AFrac, InitFileName )

      bAs(1) = ns_GAS + 1
      bAs(2) = ns_GAS + ns_AQUA
      iAs = [(i, i=bAs(1),bAs(2))]

      bAs2(1) = ns_GAS + 1
      bAs2(2) = ns_GAS + ns_AQUA*nDropletClasses
      iAs2 = [(i, i=bAs2(1),bAs2(2))]

      DO iSpc = bAs(1),bAs(2)
        READ(ChemUnit,*,IOSTAT=io_stat)  y_name(iSpc)
        IF ( y_name(iSpc) == 'Hp'    ) Hp_ind    = iSpc
        IF ( y_name(iSpc) == 'OHm'   ) OHm_ind   = iSpc
      END DO

      InitValAct(iAs) = 1.e-16_dp

      CALL Read_INI_file( InitFileName , InitValAct, InitValKat_1nD , 'AQUA' , 'INITIAL' )
      !
      CALL make_droplet_classes_by_modes(DropletClasses, AFrac, Mode, InitValAct(iAs))
      !
      ! calculate concentrations by dissolved aerosol
      DO i = 1 , nDropletClasses
        DO j = 1 , ns_AQUA
          ! this makes vector entries (A_1, A_2, .., B_1, B_2, .., C_1, C_2, ..)
          InitValAct(bAs(1) + (i-1) + (j-1)*nDropletClasses) = DropletClasses%Conc(i,j)
        END DO
      END DO
    END IF

    ! read catalytic phase species
    IF (ns_KAT>0) THEN

      lb = ns_GAS + ns_AQUA + 1
      ub = ns_GAS + ns_AQUA + ns_KAT
    
      DO iSpc = lb,ub
        READ(ChemUnit,*,IOSTAT=io_stat)  y_name(iSpc)
        henry_diff(iSpc)  = 5.0e-6_dp
        henry_accom(iSpc) = 5.0e-5_dp
        IF ( y_name(iSpc) == '[aH2O]' ) aH2O_ind = iSpc-(ns_AQUA+ns_GAS)
        IF ( y_name(iSpc) == '[H2O]'  )  H2O_ind = iSpc-(ns_AQUA+ns_GAS)
        IF ( y_name(iSpc) == '[O2]'   )   O2_ind = iSpc-(ns_AQUA+ns_GAS)
        IF ( y_name(iSpc) == '[N2]'   )   N2_ind = iSpc-(ns_AQUA+ns_GAS)
        IF ( y_name(iSpc) == '[H2]'   )   H2_ind = iSpc-(ns_AQUA+ns_GAS)
        IF (y_name(iSpc)(1:2) /= '[a') gaseous_passive_ind = [gaseous_passive_ind, iSpc-lb+1]
      END DO

      IF (H2O_ind>0) InitValKat_1nD(H2O_ind) = H2O
      IF (H2_ind>0)  InitValKat_1nD(H2_ind)  = H2
      IF (O2_ind>0)  InitValKat_1nD(O2_ind)  = O2
      IF (N2_ind>0)  InitValKat_1nD(N2_ind)  = N2

    END IF
  
    IF ( ns_AQUA>0 ) THEN
      IF (Hp_ind==0)   WRITE(*,*) '   ReadChem...Warning: Hp not in mechanism!' 
      IF (OHm_ind==0)  WRITE(*,*) '   ReadChem...Warning: OHm not in mechanism!' 
      IF (aH2O_ind==0) WRITE(*,*) '   ReadChem...Warning: [aH2O] not in mechanism!' 
    END IF

    REWIND(ChemUnit);  CLOSE(ChemUnit)

    !
    !--- Boundaries for reaction phases
    bGr = [ 1                  , nr_gas ]
    bHr = [ nr_gas+1           , nr_gas+2*nr_henry ]
    bAr = [ nr_gas+2*nr_henry+1, nr_gas+2*nr_henry+nr_liquid ]
 
    iGr = [(i, i=bGr(1),bGr(2))]
    iHr = [(i, i=bHr(1),bHr(2))]
    iAr = [(i, i=bAr(1),bAr(2))]

    bGr2 = [ 1                                  , nr_gas ]
    bHr2 = [ nr_gas+1                           , nr_gas+2*nr_henry*nDropletClasses ]
    bAr2 = [ nr_gas+2*nr_henry*nDropletClasses+1, nr_gas+2*nr_henry*nDropletClasses+nr_liquid*nDropletClasses ]

    iGr2 = [(i, i=bGr2(1),bGr2(2))]
    iHr2 = [(i, i=bHr2(1),bHr2(2))]
    iAr2 = [(i, i=bAr2(1),bAr2(2))]
    
    IF ( adiabatic_parcel ) THEN 
      iAqMassEq  = nspc+1; iTeq = nspc+2; iqEq = nspc+3; iRhoEq = nspc+4; iZeq = nspc+5
      iAqMassEq2 = (/ (nspc2+i, i=1, nDropletClasses) /)
      iTeq2 = nspc2+nDropletClasses+1; iqEq2 = nspc2+nDropletClasses+2; iRhoEq2 = nspc2+nDropletClasses+3; iZeq2 = nspc2+nDropletClasses+4
    END IF

    !-----------------------------------------------------------------------
    ! --- create guide arrays for droplet classes and aqueous species and reacs
    ! --- (these arrays are trivial but (unfortunately) needed if no aqueous phase is present,
    ! ---  somewhen they should make it possible to say nDropletClasses=0 to cut aqueous
    ! ---  processes out of a mechanism)
    CALL make_ChemSys_1_to_nD_arrays()

    CALL repeat_values_nD_vec(InitValKat, InitValKat_1nD, nD_Ptr_KAT)
    ! create reference array to enable adjusting to changing environmental conditions
    InitValKAT_Ref = InitValKAT
    InitValKAT_Ref(nD_Ptr_KAT(gaseous_passive_ind)) = InitValKAT_Ref(nD_Ptr_KAT(gaseous_passive_ind)) / ( RefTemp / Temperature0 * Pressure0 / RefPressure )
    CALL repeat_values_nD_vec(y_emi  , y_emi_1nD  , nD_Ptr_Spc)
    CALL repeat_values_nD_vec(y_depos, y_depos_1nD, nD_Ptr_Spc)

    !---  Initial pH by charge balance 
    IF ( ns_AQUA>0 ) THEN
      LWCs  = LWC_array(tBegin)
      DO i = 1 , ns_KAT
        IF ( nD_KAT(i) ) THEN ! aqueous passive species
          InitValKat(nD_Ptr_KAT(i):nD_Ptr_KAT(i+1)-1) = InitValKat(nD_Ptr_KAT(i):nD_Ptr_KAT(i+1)-1) * LWCs * mol2part
        END IF
      END DO

      IF ( pHSet .AND. (Hp_ind*OHm_ind)>0 ) THEN
        DO i = 1 , nDropletClasses
          Kappa = HpValue_from_Electroneutrality( InitValAct(nD_Ptr_spc(iAs)+i-1), LWCs(i))
          IF (Kappa/=ZERO) THEN
            IF ( Kappa > ZERO )  THEN
              InitValAct(nD_Ptr_spc(Hp_ind)+i-1)  = Kappa
              InitValAct(nD_Ptr_spc(OHm_ind)+i-1) = 1E-14_dp * LWCs(i)**2 * mol2part**2 / Kappa
            ELSE 
              InitValAct(nD_Ptr_spc(OHm_ind)+i-1) = InitValAct(nD_Ptr_spc(OHm_ind)+i-1) + InitValAct(nD_Ptr_spc(Hp_ind)+i-1) - Kappa
            END IF
          END IF
        END DO
      ELSE
        ! set pH=7
        InitValAct(nD_Ptr_spc(Hp_ind):nD_Ptr_spc(Hp_ind+1)-1)   = 1E-7_dp * LWCs * mol2part
        InitValAct(nD_Ptr_spc(OHm_ind):nD_Ptr_spc(OHm_ind+1)-1) = 1E-7_dp * LWCs * mol2part
      END IF
    END IF

  END SUBROUTINE InputChemicalData


  SUBROUTINE Read_SpeciesData(y_diff,y_acc,FileName)
    REAL(dp) :: y_acc(:) , y_diff(:)
    CHARACTER(*) :: FileName
    !
    !
    CHARACTER(LenName) :: SpeciesName
    INTEGER :: iPos, i
    LOGICAL :: Back=.FALSE.
    REAL(dp) :: mm, alpha, dg, c1
    REAL(dp) :: nue
    INTEGER :: slash
   
    CALL OpenIniFile(FileName)
    !
    i=0
    !-----------------------------------------------------------
    ! --- Gas Phase thermodynamic data
    !
    GAS: DO 
      CALL LineFile( Back                         &
      &            , Start1 = 'BEGIN_DATAGAS'     &
      &            , End    = 'END_DATAGAS'       &
      &            , Name1  = SpeciesName         &
      &            , R1 = mm, R2 = alpha, R3 = dg )
      IF (Back)   EXIT
      !
      iPos = PositionSpeciesAll(SpeciesName)
      IF ( iPos > 0) THEN
        IF (iPos<=ns_GAS+ns_AQUA ) MolMass(iPos) = mm
        IF (alpha==ZERO .AND. dg==ZERO) CYCLE GAS 
        !
        y_diff(iPos)  = 1.0_dp/(3.0_dp*dg)
        ! the E+03 is to switch the molar mass mm from g/mol to kg/mol (it is mm/1000)
        nue = SQRT(8.0e+03_dp*SI_Gas/Pi/mm)
        y_acc(iPos) = 4.0_dp/(3.0_dp*alpha*nue)
        !
      END IF
    END DO GAS
    CALL RewindFile
    CALL ClearIniFile

    !
    !-----------------------------------------------------------
    ! --- Aqua Phase thermodynamic data
    !
    IF ( ns_AQUA>0 ) THEN
      AQUA: DO 
        i = i + 1
        CALL LineFile( Back                     &
        &            , Start1 = 'BEGIN_DATAQUA' &
        &            , End    = 'END_DATAQUA'   &
        &            , Name1  = SpeciesName     &
        &            , R1 = mm, R2 = alpha      &
        &            , rLen = 2                 )
        IF (Back)  EXIT
        !
        iPos = PositionSpeciesAll(SpeciesName)
        IF ( iPos>0 ) THEN 
          IF (iPos<=ns_GAS+ns_AQUA ) THEN
            iPos = iPos - ns_GAS
            MolMass(iPos+ns_GAS) = mm
            Charge(iPos)   = alpha
          ELSE
            ! do not save aqueous passive species infos
          END IF
        END IF
      END DO AQUA
      CALL RewindFile
      CALL ClearIniFile
    END IF
    !
    !-----------------------------------------------------------
    ! --- Gas Phase RO2
    !
    IF ( hasRO2 ) THEN
      CALL LineFile( Back, Start1='BEGIN_DATARO2',    &
      &              End='END_DATARO2',               &
      &              Name1=SpeciesName,               &
      &              R1=c1)

      ALLOCATE(RO2(INT(c1)));  RO2 = 0 
      i = 0
      DO
        CALL LineFile( Back, Start1='BEGIN_DATARO2',  &
        &              End='END_DATARO2',             &
        &              Name1=SpeciesName)
     
        IF ( Back ) EXIT
        slash = INDEX(SpeciesName,'_')
        IF ( slash>0 ) SpeciesName(slash:slash)='/'
        
        IF (PositionSpeciesAll(SpeciesName) > 0) THEN
          i = i + 1
          RO2(i) = PositionSpeciesAll(SpeciesName)
        END IF
      END DO
      CALL RewindFile
      CALL ClearIniFile
      CALL CompressIntegerArray(RO2); nRO2 = SIZE(RO2)

    END IF
    !
    !-----------------------------------------------------------
    ! --- Aqua Phase RO2
    !
    IF ( hasRO2aq ) THEN
      CALL LineFile( Back, Start1='BEGIN_DATARO2aq',  &
      &              End='END_DATARO2aq',             &
      &              Name1=SpeciesName,               &
      &              R1=c1)
      
      ALLOCATE(RO2aq(INT(c1)));      RO2aq = 0
      
      i=0
      DO
        CALL LineFile( Back, Start1='BEGIN_DATARO2aq',  &
        &              End='END_DATARO2aq',             &
        &              Name1=SpeciesName)
        IF (Back) EXIT
        IF (PositionSpeciesALL(SpeciesName)>0) THEN
          i = i + 1
          RO2aq(i)=PositionSpeciesAll(SpeciesName)
        END IF
      END DO
      CALL CompressIntegerArray(RO2aq); nRO2aq = SIZE(RO2aq)
    END IF

    CALL CloseIniFile
      
  END SUBROUTINE Read_SpeciesData
  !


  SUBROUTINE Read_INI_file(FileName,Activ,Passiv,Phase,section)
    CHARACTER(*) :: FileName
    !
    REAL(dp)      :: Activ(:)
    REAL(dp)      :: Passiv(:)
    CHARACTER(*)  :: Phase, section
    !
    INTEGER :: iPos
    REAL(dp) :: c1
    CHARACTER(LenName) :: SpeciesName
    LOGICAL :: Back=.FALSE.

    ! read initial values
    CALL OpenIniFile(FileName)
    DO
      CALL LineFile( Back,                       &
      &              Start1 = 'BEGIN_'//Phase,   &
      &              Start2 = 'BEGIN_'//section, &
      &              End    = 'END_'//section,   &
      &              Name1  = SpeciesName,       &
      &              R1     = c1                 )
      IF (Back) EXIT

      IF (TRIM(SpeciesName)=='[H2O]') THEN
        WRITE(*,*) ''
        WRITE(*,*) '  NOTE :: Initial value for [H2O] is being read. This value is already assigned by'
        WRITE(*,*) '  the relative humidity specified in the *.run-file (if not specified, default value is used).'
        WRITE(*,*) '  The values, which followed the ideal gas law, are now overwritten by the value given here.'
        WRITE(*,*) '  If adiabatic_parcel is activated, this value is afterwards overwritten by q_v and has no'
        WRITE(*,*) '  influence at all.'
        WRITE(*,*) ''
        H2O = c1
        RefH2O = H2O / ( RefTemp / Temperature0 * Pressure0 / RefPressure )
        RefM   = RefN2 + RefO2 + RefH2O
      END IF
      IF (TRIM(SpeciesName)=='[N2]') THEN
        WRITE(*,*) ''
        WRITE(*,*) '  NOTE :: Initial value for [N2] is being read. This value is already assigned by'
        WRITE(*,*) '  pressure and temperature specified in the *.run-file (if not specified, default value is'
        WRITE(*,*) '  used) according to the ideal gas law. The values are now overwritten by the value given here.'
        WRITE(*,*) ''
        N2 = c1
        RefN2 = N2 / ( RefTemp / Temperature0 * Pressure0 / RefPressure )
        RefM  = RefN2 + RefO2 + RefH2O
      END IF
      IF (TRIM(SpeciesName)=='[O2]') THEN
        WRITE(*,*) ''
        WRITE(*,*) '  NOTE :: Initial value for [O2] is being read. This value is already assigned by'
        WRITE(*,*) '  pressure and temperature specified in the *.run-file (if not specified, default value is'
        WRITE(*,*) '  used) according to the ideal gas law. The values are now overwritten by the value given here.'
        WRITE(*,*) ''
        O2 = c1
        RefO2 = O2 / ( RefTemp / Temperature0 * Pressure0 / RefPressure )
        RefM  = RefN2 + RefO2 + RefH2O
      END IF
      IF (TRIM(SpeciesName)=='[H2]') THEN
        WRITE(*,*) ''
        WRITE(*,*) '  NOTE :: Initial value for [H2] is being read. This value is already assigned by'
        WRITE(*,*) '  pressure and temperature specified in the *.run-file (if not specified, default value is'
        WRITE(*,*) '  used) according to the ideal gas law. The values are now overwritten by the value given here.'
        WRITE(*,*) ''
        H2 = c1
        RefH2 = H2 / ( RefTemp / Temperature0 * Pressure0 / RefPressure )
      END IF

      iPos = PositionSpeciesAll(SpeciesName)
      IF (iPos>0) THEN
        SpeciesName = ADJUSTL(SpeciesName)
        IF (SpeciesName(1:1)=='['.AND. SpeciesName(LEN_TRIM(SpeciesName):LEN_TRIM(SpeciesName))==']' .AND. LEN_TRIM(SpeciesName)<maxLENinActDuct) THEN
          Passiv(iPos-nspc) = c1
        ELSE
          ! for the aqueous phase, the initial values should be in mol/l and they are
          ! converted to molec/cm3 while creating the droplet classes
          Activ(iPos) = Activ(iPos) + c1
        END IF
      END IF
    END DO
    CALL CloseIniFile
  END SUBROUTINE Read_INI_file

  !
  SUBROUTINE Read_Diag(DiagSpc_Pos,DiagSpc_Phase,FileName)
    INTEGER, ALLOCATABLE :: DiagSpc_Pos(:)
    CHARACTER(1), ALLOCATABLE :: DiagSpc_Phase(:)
    CHARACTER(*) :: FileName
    !
    INTEGER :: iPos, cnt, i
    CHARACTER(LenName) :: SpeciesName
    LOGICAL :: Back
    !
    ! Read initial values of aqua spc
    CALL OpenIniFile(FileName)
    cnt=0
    DO
      CALL LineFile( Back,                &
      &              Start1='BEGIN_DIAG', &
      &              End   ='END_DIAG',   &
      &              Name1 =SpeciesName   )
      IF (Back)   EXIT
      !
      IF ( ADJUSTL(SpeciesName(1:1)) /= '#' .AND. &
        & PositionSpeciesAll(SpeciesName) > 0    ) cnt = cnt + 1
    END DO
    CALL CloseIniFile
    !
    ALLOCATE(DiagSpc_Pos(cnt),DiagSpc_Phase(cnt))
    ALLOCATE(iNcdfGas(0),iNcdfAqua(0))
    DiagSpc_Pos   = 0
    DiagSpc_Phase = '-'
    !
    CALL OpenIniFile(FileName)
    cnt=0
    DO
      CALL LineFile( Back,                &
      &              Start1='BEGIN_DIAG', &
      &              End   ='END_DIAG',   &
      &              Name1 =SpeciesName   )
      IF (Back)   EXIT
      !
      IF (ADJUSTL(SpeciesName(1:1))/='#') THEN
        iPos = PositionSpeciesAll(SpeciesName)
        IF ( iPos > 0 ) THEN
          cnt = cnt + 1
          DiagSpc_Pos(cnt) = iPos
          IF      ( bGs(1) <= iPos .AND. iPos <= bGs(2) ) THEN
            DiagSpc_Phase(cnt)='g'
            iNcdfGas = [iNcdfGas,iPos]
            nNcdfGas = nNcdfGas + 1
          ELSE IF ( bAs(1) <= iPos .AND. iPos <= bAs(2) ) THEN
            DiagSpc_Phase(cnt)='a'
            iNcdfAqua = [iNcdfAqua,iPos]
            nNcdfAqua = nNcdfAqua + 1
          ELSE
            DiagSpc_Phase(cnt)='k'
          END IF
        END IF
      END IF
    END DO

    ! fill non-compressed arrays, i.e. regarding n droplet classes
    ALLOCATE(iNcdfGas2(nNcdfGas), iNcdfAqua2(nNcdfAqua*nDropletClasses))
    nNcdfGas2 = nNcdfGas
    DO i=1,nNcdfGas
      iNcdfGas2(i)    = nD_Ptr_spc(iNcdfGas(i))
    END DO
    nNcdfAqua2 = nNcdfAqua*nDropletClasses
    DO i=1,nNcdfAqua
      DO cnt=1,nDropletClasses
        iNcdfAqua2((i-1)*nDropletClasses+cnt) = nD_Ptr_spc(iNcdfAqua(i)) + cnt - 1
      END DO
    END DO

    CALL CloseIniFile
  END SUBROUTINE Read_Diag

  SUBROUTINE Read_DiagEmiss(DiagEmiss_Pos, FileName)
    INTEGER, ALLOCATABLE :: DiagEmiss_Pos(:)
    CHARACTER(*) :: FileName
    !
    INTEGER :: iPos, cnt
    CHARACTER(50) :: SpeciesName
    LOGICAL :: Back
    !
    CALL OpenIniFile(FileName)
    cnt=0
    DO
      CALL LineFile( Back,                &
      &              Start1='BEGIN_SAVEEMISS', &
      &              End   ='END_SAVEEMISS',   &
      &              Name1 =SpeciesName   )
      IF (Back)   EXIT
      !
      iPos = PositionSpeciesAll(SpeciesName)
      IF ( ADJUSTL(SpeciesName(1:1)) /= '#' .AND. &
        & iPos > 0 .AND. bGs(1) <= iPos .AND. iPos <= bGs(2) ) THEN 
          cnt = cnt + 1
      END IF
    END DO
    CALL CloseIniFile
    ! 
    ALLOCATE(DiagEmiss_Pos(cnt))
    DiagEmiss_Pos = 0
    !
    ALLOCATE(iNcdfEmiss(0))
    !
    CALL OpenIniFile(FileName)
    cnt = 0
    DO
      CALL LineFile( Back,                &
      &              Start1='BEGIN_SAVEEMISS', &
      &              End   ='END_SAVEEMISS',   &
      &              Name1 =SpeciesName   )
      IF (Back)   EXIT
      !
      IF (ADJUSTL(SpeciesName(1:1))/='#') THEN
        iPos = PositionSpeciesAll(SpeciesName)
        IF ( iPos > 0 ) THEN
          IF      ( bGs(1) <= iPos .AND. iPos <= bGs(2) ) THEN
            cnt = cnt + 1
            DiagEmiss_Pos(cnt) = iPos
            iNcdfEmiss = [iNcdfEmiss,iPos]
            nNcdfEmiss = nNcdfEmiss + 1
          ELSE
            WRITE(*,*) '  Please use gaseous emissions only. Omitting emission: '//TRIM(SpeciesName)//'.'
          END IF
        END IF
      END IF
    END DO
    CALL CloseIniFile

  END SUBROUTINE Read_DiagEmiss

  SUBROUTINE Read_AFrac(AFrac,FileName) 
    CHARACTER(*)  :: FileName
    TYPE(AFrac_T), ALLOCATABLE :: AFrac(:)

    INTEGER       :: cnt, iFrac, cnt2
    REAL(dp)      :: c1
    CHARACTER(LenName) :: SpeciesName, iFrac_str
    CHARACTER(LenName), ALLOCATABLE :: missing_names(:)
    LOGICAL       :: Back=.FALSE.

    ! Read AFRAC values
    ALLOCATE(AFrac(nFrac))

    DO iFrac = 1, nFrac
      IF (iFrac<10) THEN
        WRITE(iFrac_str,'(I1)') iFrac
      ELSE
        WRITE(iFrac_str,'(I2)') iFrac
      END IF
      CALL OpenIniFile(FileName)
      cnt=0
      cnt2=0
      DO
        CALL LineFile( Back,                                            &
        &              Start1 ='BEGIN_AQUA',                            &
        &              Start2 ='BEGIN_AFRAC'//TRIM(ADJUSTL(iFrac_str)), &
        &              End    ='END_AFRAC'//TRIM(ADJUSTL(iFrac_str)),   &
        &              Name1  =SpeciesName,                             &
        &              R1=c1 )
        IF (Back) EXIT
        cnt2=cnt2+1
        IF (PositionSpeciesAll(TRIM(SpeciesName))>0) cnt=cnt+1
      END DO
      CALL CloseIniFile
      ALLOCATE(AFrac(iFrac)%Species(cnt), AFrac(iFrac)%Frac1(cnt), missing_names(cnt2-cnt))
      !
      CALL OpenIniFile(FileName)
      cnt=0
      cnt2=0
      DO
        CALL LineFile( Back,                                            &
        &              Start1 ='BEGIN_AQUA',                            &
        &              Start2 ='BEGIN_AFRAC'//TRIM(ADJUSTL(iFrac_str)), &
        &              End    ='END_AFRAC'//TRIM(ADJUSTL(iFrac_str)),   &
        &              Name1  =SpeciesName,                             &
        &              R1=c1 )
        IF (Back) EXIT
        !
        IF (PositionSpeciesAll(TRIM(SpeciesName))>0) THEN
          cnt=cnt+1
          cnt2=cnt2+1
          AFrac(iFrac)%Species(cnt) = TRIM(SpeciesName)
          AFrac(iFrac)%Frac1(cnt)   = c1
        ELSE
          cnt2=cnt2+1
          missing_names(cnt2-cnt) = SpeciesName
        END IF
      END DO
      CALL CloseIniFile
      IF (ABS(SUM(AFrac(iFrac)%Frac1))-ONE>.01) THEN
        WRITE(*,'(A, I2, A, F7.2, A)') 'Read_AFrac :: Note: AFrac ',iFrac, ' has ', ABS(SUM(AFrac(iFrac)%Frac1)), ' g/g mass (different from 1.0).'
        CALL ask_for_continuing
      END IF
    END DO
    IF ( SIZE(missing_names)>0 ) THEN
      WRITE(*,*) ''
      WRITE(*,*) ' NOTE :: Aerosol species '
      DO cnt = 1, SIZE(missing_names)
        WRITE(*,*) '       '//TRIM(missing_names(cnt))
      END DO
      WRITE(*,*) '     is not part of the chemical mechanism and can therefore not be tracked.'
      WRITE(*,*) '     This may cause errors in charge balance, pH value and droplet activation behaviour.'
      WRITE(*,*) '     I recommend placing dummy reactions containing this species '
      WRITE(*,*) '     in the mechanism to enable tracking, as in the following example. '
      WRITE(*,*) '       CLASS: AQUA'
      WRITE(*,*) '       '//TRIM(missing_names(1))//' = '//TRIM(missing_names(1))
      WRITE(*,*) '       CONST: A: 1.0'
      WRITE(*,*) "     Also don't forget to add the corresponding values to the *.dat file."
      WRITE(*,*) '     '
      WRITE(*,*) ''
      !CALL ask_for_continuing
    END IF
  END SUBROUTINE Read_AFrac

  SUBROUTINE Read_Modes(Mode,FileName)
    CHARACTER(*) :: FileName
    TYPE(Modes_T) :: Mode
    !
    INTEGER :: cnt
    REAL(dp) :: c1,c2,c3,c4
    LOGICAL :: Back
    REAL(dp) :: LWC
    !
    ! Read SPEK values
    CALL OpenIniFile(FileName)
    CALL LineFile( Back,                 &
    &              Start1 ='BEGIN_AQUA', &
    &              Start2 ='BEGIN_SPEK', &
    &              End    ='END_SPEK',   &
    &              R1=c1                 )
    !
    nFrac = INT(c1)
    ALLOCATE( Mode%Radius(nFrac) , Mode%Number(nFrac)    &
    &       , Mode%Density(nFrac), Mode%Std_Deviation(nFrac) ) 
    !
    LWC = pseudoLWC(tBegin)
    cnt = 0
    REWIND(InputUnit_Initials)
    CALL ClearIniFile()
    DO 
      CALL LineFile( Back,                 &
      &              Start1 ='BEGIN_AQUA', &
      &              Start2 ='BEGIN_SPEK', &
      &              End    ='END_SPEK',   &
      &              R1=c1, R2=c2, R3=c3, R4=c4   )
      IF (Back)   EXIT
      !
      IF (c1 < ONE) THEN
        cnt = cnt + 1
        IF (cnt>nFrac) EXIT ! if more mode specifications than requested are there, stop reading
        
        IF (c1<=0.0 .OR. c2<=0.0 .OR. c3<=0.0 .OR. c4<=1.0 ) THEN
          WRITE(*,*) 'Improper Mode given (Number', cnt, '):', c1, c2, c3, c4
          STOP
        END IF

        Mode%Radius(cnt)        = REAL(c1,KIND=dp)
        Mode%Number(cnt)        = REAL(c2,KIND=dp)
        Mode%Density(cnt)       = REAL(c3,KIND=dp)
        Mode%Std_Deviation(cnt) = REAL(c4,KIND=dp)
        
        WRITE(*,'(A23,I1,A3,E6.2,E6.2,E6.2,F6.2)') '  Aerosol Mode Number (', cnt, '):'
        WRITE(*,'(A24,Es10.4,A4)') '    Mean Radius       : ', c1, ' [m]'
        WRITE(*,'(A24,Es10.4,A7)') '    Number            : ', c2, ' [m^-3]'
        WRITE(*,'(A24,Es10.4,A8)') '    Density           : ', c3, ' [g/m^3]'
        WRITE(*,'(A24,F5.2,A4)')   '    Standard deviation: ', c4, ' [-]'
        WRITE(*,*) ''

      END IF
    END DO
    CALL CloseIniFile
  END SUBROUTINE Read_Modes

  !
  SUBROUTINE InsertReaction(List,Line,TypeR)
    TYPE(ListReaction_T) :: List
    CHARACTER(*) :: Line(1:4)
    CHARACTER(*) :: TypeR
    !
    INTEGER :: PosColon,PosEqual,PosFactor
    CHARACTER(LenLine) :: Left,Right
    TYPE(Reaction_T), POINTER :: Reaction
    CHARACTER(LenName), ALLOCATABLE :: Ducts(:)

    !
    IF (ASSOCIATED(List%Start)) THEN
      ALLOCATE(List%End%Next)
      List%End => List%End%Next
    ELSE
      ALLOCATE(List%End)
      List%Start => List%End
    END IF
    List%LenList = List%LenList + 1
    Reaction => List%End
    PosColon = Index(Line(1),':')
    Reaction%Type  = ADJUSTL(Line(1)(PosColon+1:))
    Reaction%Line1 = TRIM(Line(2))
    Reaction%Line3 = TRIM(Line(3))
    PosEqual = Index(Reaction%Line1,' = ')
    IF (PosEqual==0) THEN
      WRITE(*,*); WRITE(*,*)
      WRITE(*,'(10X,A,I0,A)') 'ERROR: Missing separator " = " in reaction ',fNumber,':  '//TRIM(Line(2)) 
      WRITE(*,*); WRITE(*,*)
      STOP 
    ELSE
      fNumber = fNumber + 1
    END IF
   
    ! extract educts
    Left = Reaction%Line1(:PosEqual-1)
    CALL ExtractSpecies( Left, Reaction%Educt,     &
    &                    Reaction%InActEduct,      &
    &                    Reaction%nInActEd,        &
    &                    .TRUE.                    )! geändert

    ! extract products
    Right = Reaction%Line1(PosEqual+3:)
    CALL ExtractSpecies( Right, Reaction%Product,  &
    &                    Reaction%InActProduct,    &
    &                    Reaction%nInActPro,       &
    &                    .FALSE.                   )! geändert

    IF ( INDEX(Line(3),'SPECIAL') > 0 ) THEN
      Ducts = [Reaction%Educt(:)%Species , Reaction%Product(:)%Species]
      CALL ExtractConstants(Line(3),Ducts,Reaction%Constants,Reaction%TypeConstant,Reaction%Special)
    ELSE
      ! extract constants
      CALL ExtractConstants(Line(3),Ducts,Reaction%Constants,Reaction%TypeConstant)
    END IF
    Reaction%Line2  = Line(3)
    PosFactor = INDEX(Line(4),'FACTOR: ')

    IF (PosFactor > 0) THEN
      Reaction%Line4  = TRIM(Line(4))
      Reaction%Factor = TRIM(Line(4)(PosFactor+8:)) 
    ELSE
      Reaction%Line4  = ''
      Reaction%Factor = ''
    END IF
    !
    TypeR = ADJUSTL(Reaction%TypeConstant)
  END SUBROUTINE InsertReaction
  !
  !
  SUBROUTINE ReadUnits
    !
    INTEGER :: Pos
    CHARACTER(LenLine) :: LocLine
    INTEGER :: ios
    CHARACTER(400) :: iom
    LOGICAL :: comments
    !
    REWIND(InputUnit)
    DO 
      READ(InputUnit,'(A'//LenLineChar//')',IOSTAT=ios,IOMSG=iom) LocLine
      IF ( ios /= 0 ) THEN
        WRITE(*,*) ' Error Reading Units'
        WRITE(*,*) ' Error Message  ::  ',TRIM(iom)
      ELSE
        comments = (ADJUSTL(LocLine(1:1)) == '#'        .OR. &
        &           ADJUSTL(LocLine(1:7)) == 'COMMENT'  .OR. &
        &           LEN_TRIM(LocLine)     == 0 )
        IF ( comments ) CYCLE

        Pos = INDEX(LocLine,'CLASS')
        IF ( Pos > 0 ) THEN
          BACKSPACE(InputUnit)
          EXIT
        END IF
        Pos = INDEX(LocLine,'GAS')
        IF ( Pos > 0 ) READ(LocLine(Pos+3:),*) UnitGas 
        Pos = INDEX(LocLine,'AQUA')
        IF ( Pos > 0 ) THEN
          READ(LocLine(Pos+4:),*) UnitAQUA 
        END IF
      END IF
    END DO

  END SUBROUTINE ReadUnits
  !
  !
  SUBROUTINE ExtractConstants(String,Ducts,Constants,Type,SpecialForm)
    CHARACTER(*) :: String
    REAL(dp), ALLOCATABLE :: Constants(:)
    CHARACTER(*) :: Type
    !
    INTEGER :: PosColon,PosName,PosComment,PosSemiColon
    INTEGER :: i
    CHARACTER(10) :: DummyString
    CHARACTER(LEN(String)) :: LocString
    CHARACTER(LEN(String)) :: LocString1
    REAL(dp) :: Dummy
    INTEGER :: is

    ! this is for the new special formula input
    CHARACTER(LenLine)   :: SString
    TYPE(Special_T), OPTIONAL :: SpecialForm
    CHARACTER(LenName), ALLOCATABLE :: Ducts(:)
    INTEGER, ALLOCATABLE :: iSortedDucts(:)
    INTEGER :: nvs, cnt, idxDuct, lenDuct
    !
    LocString  = String
    String     = ''
    PosColon   = INDEX(LocString,':')
    Type       = LocString(:PosColon-1)
    LocString  = ADJUSTL(LocString(PosColon+1:))
    PosComment = INDEX(LocString,'#')

    IF ( PosComment > 0 ) LocString = LocString(:PosComment-1)
    
    LocString1 = LocString

    IF ( Type /= 'SPECIAL' ) THEN
      ALLOCATE(Constants(0))
      DO
        PosColon = INDEX(LocString1,':')
        IF ( PosColon > 0 ) THEN
          LocString1 = ADJUSTL(LocString1(PosColon+1:))
          READ( LocString1 , * , IOSTAT=is ) Dummy, DummyString
          PosName   = INDEX(LocString1,TRIM(DummyString))
          Constants = [Constants , Dummy]
          IF ( PosName > 0 ) LocString1 = LocString1(PosName:)
        ELSE
          EXIT
        END IF
      END DO

    ELSE
      ! special rate constant formula
      IF ( PosColon==0 ) THEN
        WRITE(*,*) '  Missing separator ":" for TypeConstant SPECIAL'
        STOP
      END IF

      PosSemiColon = Index(LocString,';')
      IF ( PosSemiColon==0 ) THEN
        WRITE(*,*) '  Missing separator ";" for TypeConstant SPECIAL'
        STOP
      END IF

      ! save formula
      SString = ADJUSTL(LocString1(:PosSemiColon-1))

      ! read number of variables in formula
      READ( LocString1(PosSemiColon+1:) , * , IOSTAT=is ) nvs

      SpecialForm%nVariables = nvs      
      SpecialForm%Formula  = ADJUSTL(SString)
      IF (INDEX(SString,'TEMP')>0) SpecialForm%Temp = .TRUE.

      ALLOCATE(SpecialForm%cVariables(nvs))
      ALLOCATE(iSortedDucts(SIZE(Ducts)) )

      CALL StringSort(Ducts,iSortedDucts)

      cnt = 0
      DO i = SIZE(Ducts), 1 , -1
        idxDuct = INDEX(SString,TRIM(Ducts(i)))

        ! if a species is catalytic ( on both sides of the reaction equation )
        IF ( i > 1 ) THEN
          IF ( Ducts(i)==Ducts(i-1) ) CYCLE
        END IF

        IF ( idxDuct > 0 ) THEN
          cnt = cnt + 1
          SpecialForm%cVariables(cnt) = ADJUSTL(Ducts(i))

          DO
            idxDuct = INDEX(SString,TRIM(Ducts(i)))
            IF ( idxDuct == 0 ) EXIT
            lenDuct = LEN_TRIM(Ducts(i))
            SString(idxDuct:idxDuct+lenDuct) = REPEAT('_',lenDuct)
          END DO

        END IF
      END DO

      IF ( SpecialForm%Temp ) THEN
        SpecialForm%cVariables(cnt+1) = 'TEMP'
      ELSE
        WRITE(*,*) '      Warning: Missing temperature "TEMP" in SPECIAL Formula :: '&
        &           ,TRIM(SpecialForm%Formula)
      END IF

    END IF
  END SUBROUTINE ExtractConstants
  !
  !
  SUBROUTINE ExtractSpecies(String,Duct,InActDuct,NumInActDucts,isEduct)
    ! IN:
    CHARACTER(*) :: String
    LOGICAL :: isEduct
    ! OUT:
    TYPE(Duct_T), POINTER :: Duct(:)
    TYPE(Duct_T), POINTER :: InActDuct(:)
    INTEGER :: NumInActDucts
    !
    INTEGER :: PosMinus,PosPlus,NumSpec,PosSpecies
    REAL(dp) :: PreFac
    CHARACTER(LenName) :: Species
    CHARACTER(LEN(String)) :: LocString
    INTEGER :: sbL, sbR
    INTEGER :: ios
    LOGICAL :: dummy=.FALSE.

    !
    LocString = String
    NumSpec   = 1

    ! count species
    DO
     PosPlus  = INDEX(LocString,' + ')
     PosMinus = INDEX(LocString,' - ')

     IF ( PosPlus>0 .AND. (PosMinus==0 .OR. PosMinus>PosPlus) ) THEN
       LocString = LocString(PosPlus+3:)
       NumSpec   = NumSpec + 1
     END IF

     IF ( PosMinus>0 .AND. (PosPlus==0 .OR. PosPlus>PosMinus) ) THEN
       LocString = LocString(PosMinus+3:)
       NumSpec   = NumSpec + 1
     END IF

     IF ( PosPlus==0 .AND. PosMinus==0 ) EXIT
    END DO

    ALLOCATE( Duct(NumSpec) , InActDuct(NumSpec) )
    LocString = String
    NumSpec   = 0
    DO
      PosPlus  = INDEX(LocString,' + ')
      PosMinus = INDEX(LocString,' - ')
      IF ( PosMinus>0 .AND. PosMinus<PosPlus ) THEN
        PreFac = mONE
      ELSE
        PreFac = ONE
      END IF
      IF      ( PosPlus>0  .AND. (PosPlus<PosMinus.OR.PosMinus==0) ) THEN
        Species   = ADJUSTL(LocString(:PosPlus-1))
        LocString = LocString(PosPlus+3:)
      ELSE IF ( PosMinus>0 .AND. (PosMinus<PosPlus.OR.PosPlus==0) ) THEN
        Species   = ADJUSTL(LocString(:PosMinus-1))
        LocString = LocString(PosMinus+3:)
      ELSE
        Species = ADJUSTL(LocString)
      END IF

      PosSpecies = SCAN(Species,SetSpecies)
      NumSpec = NumSpec + 1

      ! check if there is a dummy species e.g. 0.000 with no species
      ! or species: (dummy)
      dummy = SCAN(TRIM(ADJUSTL(Species)) , SetSpecies) == 0 .OR. &
            & INDEX(TRIM(ADJUSTL(Species)) , '(dummy)') /= 0 .OR. &
            & INDEX(TRIM(Species) , ')') - INDEX(TRIM(Species) , '(')+1 == LEN_TRIM(Species)

      IF (PosSpecies==1) THEN           
        sbL = INDEX(TRIM(Species),'[');  sbR = INDEX(TRIM(Species),']')

        ! check if species if passive
        IF ( sbL==1 .AND. LEN_TRIM(Species)==sbR-sbL+1 ) THEN
          ! works if there's just one InActEduct
          InActDuct(1)%Koeff   = PreFac
          InActDuct(1)%Species = Species
          NumInActDucts        = NumInActDucts + 1
          IF (NumInActDucts>1 .AND. isEduct) THEN
                  WRITE(*,*) "ERROR :: More than one inactive educt. Can't handle this!!"
                  STOP 'too many passive educts'
          END IF
        END IF
        Duct(NumSpec)%Koeff   = PreFac
        Duct(NumSpec)%Species = Species
      ELSE

        IF (.NOT.dummy) THEN
          READ( Species(1:PosSpecies-1) ,*, IOSTAT=ios) Duct(NumSpec)%Koeff
          IF ( ios == 0 ) THEN
            Duct(NumSpec)%Koeff   = PreFac * Duct(NumSpec)%Koeff
            sbL = INDEX(TRIM(Species(PosSpecies:)),'[');  sbR = INDEX(TRIM(Species(PosSpecies:)),']')
            IF ( sbL==1 .AND. LEN_TRIM(Species(PosSpecies:))==sbR-sbL+1 .AND. isEduct ) THEN
              WRITE(*,*) 'WARNING :: Higher order reaction involving passive species: ', TRIM(Species)
              WRITE(*,*) '           Passiveness of species is not recognized in that case. ABORT!'
            END IF
            Duct(NumSpec)%Species = Species(PosSpecies:)
          ELSE 
            WRITE(*,*) ' IOSTAT = ', ios
            WRITE(*,*) ' Error reading species and coefficients :: ', Species
            STOP
          END IF
        ELSE
          Duct(NumSpec)%Koeff   = ZERO
          Duct(NumSpec)%Species = '(dummy)'
        END IF
      END IF

      ! Syntax check for missing separators, e.g. "+CO2" instead of "+ CO2" or "CO2+" instead of "CO2 +"
      IF  ( INDEX(TRIM(ADJUSTL(Duct(NumSpec)%Species)),' ',.TRUE.) > 0 ) THEN
        WRITE(*,*); WRITE(*,*)
        WRITE(*,777) 'ERROR: Missing white space in reaction substring: '//TRIM(String)
        WRITE(*,777) '       Species = '//TRIM(Duct(NumSpec)%Species)
        WRITE(*,777) '       Check syntax in '//TRIM(SysFile)//'.sys!'
        WRITE(*,*); WRITE(*,*)
        STOP 
      END IF

      ! if no dummy species was found then add new species to hash table
      CALL InsertSpecies( Duct(NumSpec)%Species, Duct(NumSpec)%Type )

      ! if there are no further species exit the loop
      IF ( PosPlus==0 .AND. PosMinus==0 ) EXIT

    END DO
    777 FORMAT(10X,A)
  END SUBROUTINE ExtractSpecies
  
  !
  SUBROUTINE InsertSpecies(Species,Type)
    CHARACTER(*) :: Species
    CHARACTER(*) :: Type

    INTEGER :: len

    len = LEN_TRIM(Species)
    Species = TRIM(ADJUSTL(Species))

    IF (Species(1:1)=='a' .OR. Species(len:len)=='p' .OR. Species(len:len)=='m') THEN
      CALL InsertHash( ListAqua   , TRIM(ADJUSTL(Species)) , ns_AQUA)
      Type = 'Aqua'
    ELSE IF ( Species(1:1)=='[' .AND. Species(len:len)==']' &
            &                   .AND. len<maxLENinActDuct   ) THEN
      CALL InsertHash( ListNonReac, TRIM(ADJUSTL(Species)) , ns_KAT)
      Type = 'Inert'
    ELSE IF (SCAN(TRIM(ADJUSTL(Species)) , SetSpecies) == 0 .OR. &
            & INDEX(TRIM(ADJUSTL(Species)) , '(dummy)') /= 0 .OR. &
            & (INDEX(TRIM(Species) , ')') - INDEX(TRIM(Species) , '(')+1 == LEN_TRIM(Species) .AND. &
            &   LEN_TRIM(Species)>1) ) THEN
      Type = 'Dummy'
    ELSE
      CALL InsertHash( ListGas ,TRIM(ADJUSTL(Species)), ns_GAS)
      Type = 'Gas'
    END IF
  END SUBROUTINE InsertSpecies
  !
  !
  FUNCTION PositionListSpecies(Species)
    TYPE(Species_T), POINTER :: PositionListSpecies
    CHARACTER(*) :: Species
    !
    INTEGER :: PositionSpecies
    !
    PositionSpecies=0
    PositionListSpecies=>NULL()
    IF (Species(1:1)=='a'.OR.SCAN(Species,'pm')>0) THEN
      PositionSpecies=GetHash(ListAqua,TRIM(ADJUSTL(Species)))
      IF (PositionSpecies>0) PositionListSpecies=>ListAqua2(PositionSpecies)
    ELSE IF (Species(1:1)=='['.AND.Species(LEN(TRIM(Species)):LEN(TRIM(Species)))==']') THEN
      PositionSpecies=GetHash(ListNonReac,TRIM(ADJUSTL(Species)))
      IF (PositionSpecies>0) PositionListSpecies=>ListNonReac2(PositionSpecies)
    ELSE
      PositionSpecies=GetHash(ListGas,TRIM(ADJUSTL(Species)))
      IF (PositionSpecies>0) PositionListSpecies=>ListGas2(PositionSpecies)
    END IF
  END FUNCTION PositionListSpecies
  !
  !
  FUNCTION PositionSpeciesGas(Species)
    CHARACTER(*) :: Species
    !
    INTEGER :: PositionSpeciesGas
    ! 
    PositionSpeciesGas=0
    PositionSpeciesGas=GetHash(ListGas,TRIM(ADJUSTL(Species)))
  END FUNCTION
  !
  !
  FUNCTION PositionSpecies(Species)
    CHARACTER(*) :: Species
    !
    INTEGER :: PositionSpecies
    ! 
    PositionSpecies=0
    IF (Species(1:1)=='a'.OR.SCAN(Species,'pm')>0) THEN
      PositionSpecies=GetHash(ListAqua,TRIM(ADJUSTL(Species)))       &
      &               + ns_GAS
    ELSE IF (Species(1:1)=='['.AND.                                  &
    &        Species(LEN(TRIM(Species)):LEN(TRIM(Species)))==']') THEN
      PositionSpecies=GetHash(ListNonReac,TRIM(ADJUSTL(Species)))  
    ELSE
      PositionSpecies=GetHash(ListGas,TRIM(ADJUSTL(Species)))
    END IF
  END FUNCTION PositionSpecies

  !
  FUNCTION PositionSpeciesAll(Species) RESULT(Pos)
    CHARACTER(*) :: Species
    !
    INTEGER :: Pos
    INTEGER :: len

    Pos = 0
    len = LEN_TRIM(Species)

    ! Combustion system
    IF ( combustion ) THEN
      Pos = -1
      Pos = GetHash(ListGas,TRIM(ADJUSTL(Species)))
    ELSE ! tropospheric system
      ! AQUA 
      IF ( Species(1:1)=='a'.OR.SCAN(Species,'pm')>0 ) THEN
        Pos = GetHash(ListAqua,TRIM(ADJUSTL(Species)))
        IF (Pos>0) Pos = Pos + ns_GAS
 
      ! NonReac
      ELSE IF ( Species(1:1) == '[' .AND. Species(len:len) == ']' .AND. &
      &        len < maxLENinActDuct) THEN
        Pos = GetHash(ListNonReac,TRIM(ADJUSTL(Species)))    
        IF (Pos>0) Pos = Pos + ns_GAS + ns_AQUA
        
      ! GAS
      ELSE
        Pos = GetHash(ListGas,TRIM(ADJUSTL(Species)))
      END IF
    END IF
  END FUNCTION PositionSpeciesAll
  !
  FUNCTION PositionAtom(Atom) RESULT(Pos)
    CHARACTER(*) :: Atom
    !
    INTEGER :: Pos
    ! 
    Pos = 0
    Pos = GetHash(ListAtoms,TRIM(ADJUSTL(Atom)))
  END FUNCTION PositionAtom
  !
  !
  SUBROUTINE OpenFile(FileName,Type)
    CHARACTER(*) :: Filename
    CHARACTER(*) :: Type
    !
    LOGICAL :: ExistFile
    !
    INQUIRE(FILE=TRIM(Filename)//'.'//TRIM(Type),EXIST=ExistFile)
    IF (ExistFile) THEN
      OPEN(UNIT=InputUnit,FILE=TRIM(Filename)//'.'//TRIM(Type),STATUS='UNKNOWN')
    END IF
  END SUBROUTINE OpenFile
  !
  !
  SUBROUTINE CloseFile(FileName,Type)
    CHARACTER(*) :: Filename
    CHARACTER(*) :: Type
    !
    LOGICAL :: ExistFile
    !
    INQUIRE(FILE=TRIM(Filename)//'.'//TRIM(Type),EXIST=ExistFile)
    IF (ExistFile) CLOSE(UNIT=InputUnit)
  END SUBROUTINE CloseFile
  !
  !
  SUBROUTINE ReadSystem(FileName)
    CHARACTER(*) :: Filename
    !
    LOGICAL :: Out

    FileName = FileName(:INDEX(FileName,'.')-1)
    !
    CALL InitHashTable(ListAqua,100)
    CALL InitHashTable(ListGas,100)
    CALL InitHashTable(ListNonReac,100)
    CALL OpenFile(FileName,'sys')
    CALL ReadUnits
    DO
      CALL ReadReaction(Out); IF (Out) EXIT
    END DO
    CALL CloseFile(FileName,'sys')
    ALLOCATE(ListGas2(ns_GAS))
    CALL HashTableToList(ListGas,ListGas2)
    CALL SortList(ListGas2)
    CALL ListToHashTable(ListGas2,ListGas)
    ALLOCATE(ListAqua2(ns_AQUA))
    CALL HashTableToList(ListAqua,ListAqua2)
    CALL SortList(ListAqua2)
    CALL ListToHashTable(ListAqua2,ListAqua)
    ALLOCATE(ListNonReac2(ns_KAT))
    CALL HashTableToList(ListNonReac,ListNonReac2)
    CALL SortList(ListNonReac2)
    CALL ListToHashTable(ListNonReac2,ListNonReac)
    !
  END SUBROUTINE ReadSystem
  !
  !
  SUBROUTINE SortList(List)
    TYPE(Species_T) :: List(:)
    !
    TYPE(Species_T) :: Temp
    INTEGER :: i,j
    !
    DO i=1,SIZE(List)
      DO j=1,SIZE(List)-i
        IF (List(j+1)%Species<List(j)%Species) THEN
          Temp=List(j+1)
          List(j+1)=List(j)
          List(j)=Temp
        END IF
      END DO
    END DO
    DO i=1,SIZE(List)
      IF (List(i)%Species=='OHm') THEN
        Temp=List(i)
        IF (i==SIZE(List)) EXIT
        DO j=i,SIZE(List)-1
          List(j)=List(j+1)
        END DO
        List(SIZE(List))=Temp
        EXIT
      END IF
    END DO
    DO i=1,SIZE(List)
      IF (List(i)%Species=='Hp') THEN
        Temp=List(i)
        IF (i==SIZE(List)) EXIT
        DO j=i,SIZE(List)-1
          List(j)=List(j+1)
        END DO
        List(SIZE(List))=Temp
        EXIT
      END IF
    END DO
  END SUBROUTINE SortList
  !
  !
  SUBROUTINE ListToHashTable(List,Table)
    TYPE(Species_T) :: List(:)
    TYPE(hash_tbl_sll) :: Table
    !
    INTEGER :: i
    DO i=1,SIZE(List)
      CALL Table%put(TRIM(ADJUSTL(List(i)%Species)),i)
    END DO
  END SUBROUTINE ListToHashTable
  !
  !
  SUBROUTINE HashTableToList(Table,List)
    TYPE(hash_tbl_sll) :: Table
    TYPE(Species_T) :: List(:)
    !
    INTEGER :: i,j
    TYPE(sllist), POINTER :: child => NULL()
    !
    DO i=LBOUND(table%vec,dim=1), UBOUND(table%vec,dim=1)
      IF (ALLOCATED(table%vec(i)%key)) THEN
        DO j=1,SIZE(table%vec(i)%key)
          List(table%vec(i)%Val)%Species(j:j)=table%vec(i)%key(j)
        END DO
      END IF
      Child=>table%vec(i)%Child
      DO
        IF (ASSOCIATED(Child)) THEN
          DO j=1,SIZE(Child%key)
            List(Child%Val)%Species(j:j)=Child%key(j)
          END DO
          Child=>Child%Child
        ELSE
          EXIT
        END IF
      END DO
    END DO
  END SUBROUTINE HashTableToList
  !
  !
  SUBROUTINE ReactionListToArray(ReacList,ReacArr)
    TYPE(ListReaction_T) :: ReacList
    TYPE(ReactionStruct_T), ALLOCATABLE, INTENT(OUT) :: ReacArr(:)
    INTEGER :: i, j, TmpArraySize
    INTEGER :: ListLen=0
    !
    Current=>ReacList%Start
    DO 
      IF (.NOT.ASSOCIATED(Current)) EXIT
      ListLen=ListLen+1
      Current=>Current%Next
    END DO
    !
    ALLOCATE(ReacArr(ListLen))
    !
    Current=>ReacList%Start
    i=1
    !
    DO
      IF (ASSOCIATED(Current)) THEN
        ReacArr(i)%Type=Current%Type
        ReacArr(i)%TypeConstant=Current%TypeConstant
        ReacArr(i)%Line1=Current%Line1
        ReacArr(i)%Line2=Current%Line2
        ReacArr(i)%Factor=Current%Factor
        !
        TmpArraySize=SIZE(Current%Educt)
        ALLOCATE(ReacArr(i)%Educt(TmpArraySize))
        DO j=1,TmpArraySize
          ReacArr(i)%Educt(j)%Species=Current%Educt(j)%Species
          ReacArr(i)%Educt(j)%Type=Current%Educt(j)%Type
          ReacArr(i)%Educt(j)%Koeff=Current%Educt(j)%Koeff
        END DO
        !
        TmpArraySize=SIZE(Current%Product)
        ALLOCATE(ReacArr(i)%Product(TmpArraySize))
        DO j=1,TmpArraySize
          ReacArr(i)%Product(j)%Species=Current%Product(j)%Species
          ReacArr(i)%Product(j)%Type=Current%Product(j)%Type
          ReacArr(i)%Product(j)%Koeff=Current%Product(j)%Koeff
        END DO
        !
        ALLOCATE(ReacArr(i)%Constants(SIZE(Current%Constants)))
        ReacArr(i)%Constants=Current%Constants
        !
        TmpArraySize=SIZE(Current%Educt)
        ALLOCATE(ReacArr(i)%InActEduct(TmpArraySize))
        ALLOCATE(ReacArr(i)%InActEductSpc(TmpArraySize))
        DO j=1,TmpArraySize
          ReacArr(i)%InActEduct(j)=Current%InActEduct(j)%Koeff
          ReacArr(i)%InActEductSpc(j)=Current%InActEduct(j)%Species
        END DO
        !
        TmpArraySize=SIZE(Current%InActProduct)
        ALLOCATE(ReacArr(i)%InActProduct(TmpArraySize))
        DO j=1,TmpArraySize
          ReacArr(i)%InActProduct(j)=Current%InActProduct(j)%Koeff    
        END DO
        !    
        ReacArr(i)%nInActEd=Current%nInActEd
        ReacArr(i)%nInActPro=Current%nInActPro
      ELSE
        EXIT
      END IF
      Current=>Current%Next
      i=i+1
    END DO
  END SUBROUTINE ReactionListToArray
  !
  !
  !
  SUBROUTINE AllListsToArray(ReacStruct,LGas,LHenry,LAqua,LDiss)
    ! OUT:
    TYPE(ReactionStruct_T), ALLOCATABLE :: ReacStruct(:)
    ! IN:
    TYPE(ListReaction_T) :: LGas, LHenry, LAqua, LDiss
    
    INTEGER :: i, j, iList, iEq
    INTEGER :: nList

    INTEGER :: icnt(48), icntFAC(10), iHen

    ! #Reactions
    nreac  = nr_gas + 2*nr_henry + nr_aqua + 2*nr_diss
    nreac2 = nr_gas + (2*nr_henry + nr_aqua + 2*nr_diss)*nDropletClasses

    ! #Species
    nspc  = ns_GAS + ns_AQUA
    nspc2 = ns_GAS + ns_AQUA*nDropletClasses

    ns_KAT2 = 0
    DO i = 1 , ns_KAT
      IF (TRIM(ListNonReac2(i)%Species(1:2))=='[a') THEN
        ns_KAT2 = ns_KAT2 + nDropletClasses
      ELSE
        ns_KAT2 = ns_KAT2 + 1
      END IF
    END DO
    
    nr_liquid = nr_aqua + 2*nr_diss

    !-----------------------------------------------------------------------
    ! --- Set logicals
    !-----------------------------------------------------------------------
    IF ( ns_GAS   > 0 ) hasGasSpc   = .TRUE.; nPhases = nPhases + 1
    IF ( ns_AQUA  > 0 ) hasAquaSpc  = .TRUE.; nPhases = nPhases + 1

    IF ( nr_gas    > 0 ) hasGasReac    = .TRUE.
    IF ( nr_aqua   > 0 ) hasAquaReac   = .TRUE.
    IF ( nr_henry  > 0 ) hasHenryReac  = .TRUE.
    IF ( nr_diss   > 0 ) hasDissReac   = .TRUE.
    IF ( nr_liquid > 0 ) hasLiquidReac = .TRUE.

    hasPhotoReac  = (nr_G_photo+nr_A_photo) > 0
    
    nList = 0
    IF (ASSOCIATED(LGas%Start))    nList=nList+1
    IF (ASSOCIATED(LHenry%Start))  nList=nList+1
    IF (ASSOCIATED(LAqua%Start))   nList=nList+1
    IF (ASSOCIATED(LDiss%Start))   nList=nList+1
    ALLOCATE(CompleteReactionList(nList))
    nList = 0
    IF (ASSOCIATED(LGas%Start))   THEN
      nList=nList+1; CompleteReactionList(nList)=LGas
    END IF
    IF (ASSOCIATED(LHenry%Start)) THEN
      nList=nList+1; CompleteReactionList(nList)=LHenry
    END IF
    IF (ASSOCIATED(LAqua%Start))  THEN
      nList=nList+1; CompleteReactionList(nList)=LAqua
    END IF
    IF (ASSOCIATED(LDiss%Start))  THEN
      nList=nList+1; CompleteReactionList(nList)=LDiss
    END IF

    ALLOCATE( ReacStruct(nreac) )
  
    CALL AllocateRTarrays
 
    i=1
    iHen = 0
    icnt = 0
    icntFAC = 0

    DO iList=1,nList
      Current => CompleteReactionList(iList)%Start
      DO WHILE (ASSOCIATED(Current)) 
        ReacStruct(i)%Type   = Current%Type
        ReacStruct(i)%Line1  = Current%Line1
        ReacStruct(i)%Line2  = Current%Line2
        ReacStruct(i)%Line3  = Current%Line3

        ReacStruct(i)%Line4  = Current%Line4
        ReacStruct(i)%Factor = Current%Factor
        ReacStruct(i)%TypeConstant = Current%TypeConstant
         
        CALL Setup_iFACTOR( i, icntFAC, Current%Factor)

        ! forward direaction
        CALL Setup_ReacParameter( i , icnt ,           &
                                & Current%TypeConstant,&
                                & Current%Constants,   &
                                & TRIM(Current%Line1)  )
        !
        ReacStruct(i)%nActEd  = SIZE(Current%Educt)
        ReacStruct(i)%nActPro = SIZE(Current%Product)
        ALLOCATE( ReacStruct(i)%Educt(ReacStruct(i)%nActEd),   &
                & ReacStruct(i)%Product(ReacStruct(i)%nActPro) )
            
        ! aus A + B -> .5 C + .5 D + .5 C + .5 D ===> A + B -> C + D

        DO j = 1 , ReacStruct(i)%nActEd
          ReacStruct(i)%Educt(j)%Species  = Current%Educt(j)%Species
          ReacStruct(i)%Educt(j)%Type     = Current%Educt(j)%Type
          ReacStruct(i)%Educt(j)%Koeff    = Current%Educt(j)%Koeff
          ReacStruct(i)%Educt(j)%iSpecies = PositionSpeciesAll(Current%Educt(j)%Species)
        END DO
   
        ReacStruct(i)%Product = CleanUpDucts(Current%Product)
        ReacStruct(i)%nActPro = SIZE(ReacStruct(i)%Product)

        IF ( ReacStruct(i)%TypeConstant == 'SPECIAL' ) THEN
          j = SIZE(Current%Special%cVariables)
          IF (Current%Special%Temp) THEN
            ALLOCATE(ReacStruct(i)%Special%cVariables(j),ReacStruct(i)%Special%iVariables(j-1))
          ELSE
            ALLOCATE(ReacStruct(i)%Special%cVariables(j),ReacStruct(i)%Special%iVariables(j))
          END IF

          ReacStruct(i)%Special%nVariables = j
          ReacStruct(i)%Special%Formula  = Current%Special%Formula
          ReacStruct(i)%Special%Temp     = Current%Special%Temp
          DO j = 1,ReacStruct(i)%Special%nVariables
            ReacStruct(i)%Special%cVariables(j) = Current%Special%cVariables(j)
            IF ( Current%Special%cVariables(j) /= 'TEMP' ) THEN 
              ReacStruct(i)%Special%iVariables(j) = PositionSpeciesAll(Current%Special%cVariables(j))
            END IF
          END DO
        ELSE
          ReacStruct(i)%nConst  = SIZE(Current%Constants)
          ALLOCATE(ReacStruct(i)%Constants(ReacStruct(i)%nConst))
          ReacStruct(i)%Constants = Current%Constants
        END IF
        !
        IF ( Current%Type == 'HENRY' ) THEN
          ReacStruct(i)%direction = 'GA'
          ReacStruct(i)%HenrySpc  = PositionSpeciesAll(ReacStruct(i)%Educt(1)%Species)

          icnt(29) = icnt(29) + 1
          iR%iHENRY(icnt(29),1) = i
          iR%iHENRY(icnt(29),2) = PositionSpeciesAll(ReacStruct(i)%Educt(1)%Species)
          iR%iHENRY(icnt(29),3) = i + 1
          iR%iHENRY(icnt(29),4) = PositionSpeciesAll(ReacStruct(i)%Product(1)%Species)
        END IF
        !
        !
        ReacStruct(i)%nInActEd  = Current%nInActEd
        ReacStruct(i)%nInActPro = Current%nInActPro
        ALLOCATE( ReacStruct(i)%InActEduct(Current%nInActEd),      &
                & ReacStruct(i)%InActEductSpc(Current%nInActEd),   & 
                & ReacStruct(i)%InActProduct(Current%nInActPro),   &
                & ReacStruct(i)%InActProductSpc(Current%nInActPro) )

        DO j = 1 , Current%nInActEd
          ReacStruct(i)%InActEduct(j)    = Current%InActEduct(j)%Koeff
          ReacStruct(i)%InActEductSpc(j) = Current%InActEduct(j)%Species
          nFirst_orderKAT = nFirst_orderKAT + 1
        END DO
        DO j = 1 , Current%nInActPro
          ReacStruct(i)%InActProduct(j)    = Current%InActProduct(j)%Koeff    
          ReacStruct(i)%InActProductSpc(j) = Current%InActProduct(j)%Species
        END DO

        ! ----------
        ! sum educt coefficients to detect order of reaction for unit correction if needed
        ReacStruct(i)%SumAqCoef = SUM(Current%Educt%Koeff)
        ! if factor is RO2aq we need an additional unit correction to make the RO2aq concentration to mol/l
        IF (TRIM(ReacStruct(i)%Factor)=="$RO2aq") THEN
          ReacStruct(i)%SumAqCoef = ReacStruct(i)%SumAqCoef + ONE
        END IF
        ! (the only other aqueous factor is aH2O, given directly in mol/l in rate calculation)
        ! ----------

        IF ( ReacStruct(i)%Type=='AQUA'.OR. ReacStruct(i)%Type=='DISS' ) THEN
          !IF ( ReacStruct(i)%SumAqCoef /= ONE ) nr_HOaqua = nr_HOaqua + 1
          IF ( ReacStruct(i)%SumAqCoef == TWO   ) THEN
            nr_SOaqua = nr_SOaqua + 1
          ELSE IF ( ReacStruct(i)%SumAqCoef == THREE ) THEN
            nr_TOaqua = nr_TOaqua + 1
          ELSE IF ( ReacStruct(i)%SumAqCoef /= ONE   ) THEN
            nr_HOaqua = nr_HOaqua + 1
          END IF
        END IF
        !
        ! for equilibrium reactions save <-- direction
        SELECT CASE (Current%Type)
          CASE ('DISS','HENRY')
            i=i+1
            iEq = INDEX(Current%Line1,' = ')
            ReacStruct(i)%Type   = Current%Type
            ReacStruct(i)%Line1  = TRIM(Current%Line1(iEq+3:))//' = '//TRIM(Current%Line1(:iEq))
            ReacStruct(i)%Line2  = 'reverse reaction'
            ReacStruct(i)%bR     = .TRUE.
            ReacStruct(i)%Line3  = Current%Line3
            ReacStruct(i)%TypeConstant = Current%TypeConstant
            ReacStruct(i)%Factor = ''

            ReacStruct(i)%nActEd  = SIZE(Current%Product)
            ReacStruct(i)%nActPro = SIZE(Current%Educt)
            ALLOCATE( ReacStruct(i)%Educt(ReacStruct(i)%nActEd),   &
                    & ReacStruct(i)%Product(ReacStruct(i)%nActPro) )

            DO j=1,ReacStruct(i)%nActEd
              ReacStruct(i)%Educt(j)%Species  = Current%Product(j)%Species
              ReacStruct(i)%Educt(j)%Type     = Current%Product(j)%Type
              ReacStruct(i)%Educt(j)%Koeff    = Current%Product(j)%Koeff
              ReacStruct(i)%Educt(j)%iSpecies = PositionSpeciesAll(Current%Product(j)%Species)
            END DO
            DO j=1,ReacStruct(i)%nActPro
              ReacStruct(i)%Product(j)%Species  = Current%Educt(j)%Species
              ReacStruct(i)%Product(j)%Type     = Current%Educt(j)%Type
              ReacStruct(i)%Product(j)%Koeff    = Current%Educt(j)%Koeff
              ReacStruct(i)%Product(j)%iSpecies = PositionSpeciesAll(Current%Educt(j)%Species)
            END DO
            
            IF ( ReacStruct(i)%TypeConstant == 'SPECIAL' ) THEN
              ReacStruct(i)%Special%nVariables = ReacStruct(i-1)%Special%nVariables
              ReacStruct(i)%Special%Formula    = ReacStruct(i-1)%Special%Formula
              ReacStruct(i)%Special%Temp       = ReacStruct(i-1)%Special%Temp
              ReacStruct(i)%Special%cVariables = ReacStruct(i-1)%Special%cVariables
              ReacStruct(i)%Special%iVariables = ReacStruct(i-1)%Special%iVariables
            ELSE
              ReacStruct(i)%nConst  = SIZE(Current%Constants)
              ALLOCATE(ReacStruct(i)%Constants(ReacStruct(i)%nConst))
              ReacStruct(i)%Constants = Current%Constants
            END IF
            
            IF ( Current%Type == 'HENRY' ) THEN
              ReacStruct(i)%direction = 'AG'
              ReacStruct(i)%HenrySpc  = PositionSpeciesAll(ReacStruct(i)%Product(1)%Species)
            END IF
            !
            ReacStruct(i)%nInActEd  = Current%nInActPro
            ReacStruct(i)%nInActPro = Current%nInActEd
            ALLOCATE( ReacStruct(i)%InActEduct(Current%nInActPro),    &
                    & ReacStruct(i)%InActEductSpc(Current%nInActPro), & 
                    & ReacStruct(i)%InActProduct(Current%nInActEd),   &
                    & ReacStruct(i)%InActProductSpc(Current%nInActEd) )

            DO j = 1 , Current%nInActPro
              ReacStruct(i)%InActEduct(j)    = Current%InActProduct(j)%Koeff
              ReacStruct(i)%InActEductSpc(j) = Current%InActProduct(j)%Species
              nFirst_orderKAT = nFirst_orderKAT + 1
            END DO
            DO j = 1 , Current%nInActEd
              ReacStruct(i)%InActProduct(j)    = Current%InActEduct(j)%Koeff    
              ReacStruct(i)%InActProductSpc(j) = Current%InActEduct(j)%Species
            END DO
            
            ReacStruct(i)%SumAqCoef = SUM(Current%Product%Koeff)

            IF ( ReacStruct(i)%Type=='AQUA'.OR. ReacStruct(i)%Type=='DISS' ) THEN
              !IF ( ReacStruct(i)%SumAqCoef /= ONE ) nr_HOaqua = nr_HOaqua + 1
              IF ( ReacStruct(i)%SumAqCoef == TWO   ) THEN
                nr_SOaqua = nr_SOaqua + 1
              ELSE IF ( ReacStruct(i)%SumAqCoef == THREE ) THEN
                nr_TOaqua = nr_TOaqua + 1
              ELSE IF ( ReacStruct(i)%SumAqCoef /= ONE   ) THEN
                nr_HOaqua = nr_HOaqua + 1
              END IF
            END IF
            !
        END SELECT
        !
        Current=>Current%Next
        i=i+1
      END DO
    END DO
    
    ! build the array for mass action products of inactive species
    ALLOCATE( iFO_kat(nfirst_orderKAT,2) )
    nFirst_orderKAT = 0

    ! counting the aquatic reactions with more than one educt and
    ! specifying reactions to be cut off eventually, but this isn't currently done
    ! for numerical life-keeping reasons
    !ALLOCATE( iR%iHOaqua(nr_HOaqua), iR%HOaqua(nr_HOaqua), iCutOffReacs(0) )
    ALLOCATE( iR%iHOaqua(nr_HOaqua), iR%HOaqua(nr_HOaqua),  &
            & iR%iSOaqua(nr_SOaqua), iR%iTOaqua(nr_TOaqua), &
            & iCutOffReacs(0)                               )
    nr_SOaqua = 0
    nr_TOaqua = 0
    nr_HOaqua = 0

    DO i=1,nreac

      IF ( ReacStruct(i)%nInActEd > 0 ) THEN
        nFirst_orderKAT = nFirst_orderKAT + 1
        iFO_kat(nfirst_orderKAT,1) = i
        iFO_kat(nfirst_orderKAT,2) = PositionSpeciesAll(ReacStruct(i)%InActEductSpc(1)) - nspc
      END IF

      IF ( ReacStruct(i)%Type=='AQUA'.OR. ReacStruct(i)%Type=='DISS' ) THEN
        !IF ( ReacStruct(i)%SumAqCoef /= ONE ) THEN
        !  nr_HOaqua = nr_HOaqua + 1
        !  iR%iHOaqua(nr_HOaqua) = i
        !  ! iR%HOaqua saves the exponent for unit correction of the rate constant in Rates_Mod
        !  ! this has to be one less than the sum of educt coefficients, see the comment in Rates_Mod
        !  iR%HOaqua(nr_HOaqua) = ReacStruct(i)%SumAqCoef - ONE
        !END IF
        IF ( ReacStruct(i)%SumAqCoef == TWO ) THEN
          nr_SOaqua = nr_SOaqua + 1
          iR%iSOaqua(nr_SOaqua) = i
        ELSE IF ( ReacStruct(i)%SumAqCoef == THREE ) THEN
          nr_TOaqua = nr_TOaqua + 1
          iR%iTOaqua(nr_TOaqua) = i
        ELSE IF ( ReacStruct(i)%SumAqCoef /= ONE ) THEN
          nr_HOaqua = nr_HOaqua + 1
          iR%iHOaqua(nr_HOaqua) = i
          ! iR%HOaqua saves the exponent for unit correction of the rate constant in Rates_Mod
          ! this has to be one less than the sum of educt coefficients, see the comment in Rates_Mod
          IF (ABS(INT(ReacStruct(i)%SumAqCoef - ONE) - (ReacStruct(i)%SumAqCoef - ONE)) > 1.0e-12_dp) THEN
            WRITE(*,*) 'Broken order stoichiometric coefficients for educts are suppressed right now,&
                      & because integer exponentiation is faster. To enable, turn iR%HOaqua into REAL.'
            STOP
          END IF
          !iR%HOaqua(nr_HOaqua) = ReacStruct(i)%SumAqCoef - ONE
          iR%HOaqua(nr_HOaqua) = NINT(ReacStruct(i)%SumAqCoef - ONE)
        END IF
      END IF

      !IF (ReacStruct(i)%Type=='AQUA' .OR. ReacStruct(i)%Type=='DISS' .OR. ReacStruct(i)%Type=='HENRY') THEN
      IF (ReacStruct(i)%Type=='AQUA') THEN
        iCutOffReacs = [iCutOffReacs, i]
      END IF

    END DO

    CONTAINS 

      FUNCTION CleanUpDucts(DuctsIN) RESULT(Ducts)
        TYPE(Duct_T), INTENT(IN)  :: DuctsIN(:)
        TYPE(Duct_T), ALLOCATABLE :: Ducts(:)

        INTEGER :: i, n, Len, nDupes
        INTEGER,        ALLOCATABLE :: tmp_iSpc(:), tmp_iSpc_sort(:) 
        CHARACTER(LenType), ALLOCATABLE :: tmp_Type(:)
        CHARACTER(LenName), ALLOCATABLE :: tmp_cSpc(:)
        REAL(dp),       ALLOCATABLE :: tmp_Koeff(:)
        INTEGER,        ALLOCATABLE :: Perm(:)

        
        n = SIZE(DuctsIN)
        
        IF ( n > 1 ) THEN
        
          ALLOCATE(tmp_cSpc(n), tmp_Type(n), tmp_Koeff(n), tmp_iSpc(n), tmp_iSpc_sort(n))
          DO i = 1 , n 
            tmp_cSpc(i)  = DuctsIN(i)%Species
            tmp_Type(i)  = DuctsIN(i)%Type
            tmp_Koeff(i) = DuctsIN(i)%Koeff
            tmp_iSpc(i)  = PositionSpeciesAll(DuctsIN(i)%Species)
          END DO 

          tmp_iSpc_sort = tmp_iSpc
          CALL SortVecAsc(tmp_iSpc_sort,Perm)

          DO i = 1 , n-1
            IF ( tmp_iSpc_sort(i) == tmp_iSpc_sort(i+1) ) THEN
              tmp_Koeff(Perm(i+1))  = tmp_Koeff(Perm(i+1)) + tmp_Koeff(Perm(i))
              tmp_iSpc(Perm(i)) = -1              
            END IF
          END DO

          nDupes = COUNT(tmp_iSpc==-1)

          IF ( nDupes == 0 ) THEN

            ALLOCATE(Ducts(n))
            DO i = 1 , n
              Ducts(i)%iSpecies = tmp_iSpc(i)
              Ducts(i)%Species  = tmp_cSpc(i)(:)
              Ducts(i)%Type     = tmp_Type(i)(:)
              Ducts(i)%Koeff    = tmp_Koeff(i)
            END DO 

          ELSE 

            ALLOCATE(Ducts(n-nDupes))
            Len = 0
            
            DO i = 1 , n
              IF ( tmp_iSpc(i) > 0 ) THEN
                Len = Len + 1
                Ducts(Len)%iSpecies = tmp_iSpc(i)
                Ducts(Len)%Species  = tmp_cSpc(i)(:)
                Ducts(Len)%Type     = tmp_Type(i)(:)
                Ducts(Len)%Koeff    = tmp_Koeff(i)
              END IF
            END DO
              
          END IF

        ELSE
          ALLOCATE(Ducts(1))
          Ducts%iSpecies = PositionSpeciesAll(DuctsIN(1)%Species)
          Ducts%Species  = DuctsIN%Species
          Ducts%Type     = DuctsIN%Type
          Ducts%Koeff    = DuctsIN%Koeff
        END IF
      END FUNCTION CleanUpDucts
  END SUBROUTINE AllListsToArray
  !
  SUBROUTINE Setup_iFACTOR(iReac,icntFAC,Factor)
    CHARACTER(*), INTENT(IN) :: Factor
    INTEGER,      INTENT(IN) :: iReac
    INTEGER,      INTENT(INOUT) :: icntFAC(:)

    SELECT CASE (Factor)
      CASE ('$H2')   ; icntFAC(1)=icntFAC(1)+1  ; iR%iFAC_H2(icntFAC(1))=iReac
      CASE ('$O2N2') ; icntFAC(2)=icntFAC(2)+1  ; iR%iFAC_O2N2(icntFAC(2))=iReac
      CASE ('$M')    ; icntFAC(3)=icntFAC(3)+1  ; iR%iFAC_M(icntFAC(3))=iReac
      CASE ('$O2')   ; icntFAC(4)=icntFAC(4)+1  ; iR%iFAC_O2(icntFAC(4))=iReac
      CASE ('$N2')   ; icntFAC(5)=icntFAC(5)+1  ; iR%iFAC_N2(icntFAC(5))=iReac
      CASE ('$H2O')  ; icntFAC(6)=icntFAC(6)+1  ; iR%iFAC_H2O(icntFAC(6))=iReac
      CASE ('$RO2')  ; icntFAC(7)=icntFAC(7)+1  ; iR%iFAC_RO2(icntFAC(7))=iReac   ; hasRO2=.TRUE.
      CASE ('$O2O2') ; icntFAC(8)=icntFAC(8)+1  ; iR%iFAC_O2O2(icntFAC(8))=iReac
      CASE ('$aH2O') ; icntFAC(9)=icntFAC(9)+1  ; iR%iFAC_aH2O(icntFAC(9))=iReac
      CASE ('$RO2aq'); icntFAC(10)=icntFAC(10)+1; iR%iFAC_RO2aq(icntFAC(10))=iReac; hasRO2aq=.TRUE.
    END SELECT
  END SUBROUTINE Setup_iFACTOR


  SUBROUTINE Setup_ReacParameter(iReac,icnt,TypeR,C,Line1)
    REAL(dp), INTENT(IN) :: C(:)
    CHARACTER(*),   INTENT(IN) :: TypeR
    CHARACTER(*),   INTENT(IN) :: Line1
    INTEGER,        INTENT(IN) :: iReac
    INTEGER,        INTENT(INOUT) :: icnt(48)

    SELECT CASE ( TRIM(TypeR) )
      CASE ('PHOTABC')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
        icnt(3)=icnt(3)+1; iR%iPHOTabc(icnt(3))=iReac; iR%PHOTabc(icnt(3),:)=C 
      CASE ('PHOTMCM')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
        icnt(4)=icnt(4)+1; iR%iPHOTmcm(icnt(4))=iReac; iR%PHOTmcm(icnt(4),:)=C 
      CASE ('PHOTAB')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(2)=icnt(2)+1; iR%iPHOTab(icnt(2))=iReac;  iR%PHOTab(icnt(2),:)=C
      CASE ('CONST')
        IF ( SIZE(C)<1 ) CALL ErrorMSG(iReac,Line1)
        icnt(1)=icnt(1)+1; iR%iCONST(icnt(1))=iReac;   iR%CONST(icnt(1))=C(1)
      CASE ('TEMP')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
        icnt(46)=icnt(46)+1; iR%iTEMP(icnt(46))=iReac; iR%TEMP(icnt(46),:)=C 
      CASE ('TEMP1')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(5)=icnt(5)+1; iR%iTEMP1(icnt(5))=iReac; iR%TEMP1(icnt(5),:)=C 
      CASE ('TEMP2')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(6)=icnt(6)+1; iR%iTEMP2(icnt(6))=iReac; iR%TEMP2(icnt(6),:)=C 
      CASE ('TEMP3')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(7)=icnt(7)+1; iR%iTEMP3(icnt(7))=iReac; iR%TEMP3(icnt(7),:)=C 
      CASE ('TEMP4')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(8)=icnt(8)+1; iR%iTEMP4(icnt(8))=iReac; iR%TEMP4(icnt(8),:)=C 
      CASE ('TROE')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
        icnt(9)=icnt(9)+1; iR%iTROE(icnt(9))=iReac;  iR%TROE(icnt(9),:)=C 
      CASE ('TROEF')
        IF ( SIZE(C)<5 ) CALL ErrorMSG(iReac,Line1)
        icnt(10)=icnt(10)+1; iR%iTROEf(icnt(10))=iReac;   iR%TROEf(icnt(10),:)=C  
      CASE ('TROEQ')
        IF ( SIZE(C)<6 ) CALL ErrorMSG(iReac,Line1)
        icnt(11)=icnt(11)+1; iR%iTROEq(icnt(11))=iReac;   iR%TROEq(icnt(11),:)=C 
      CASE ('TROEQF')
        IF ( SIZE(C)<7 ) CALL ErrorMSG(iReac,Line1)
        icnt(12)=icnt(12)+1; iR%iTROEqf(icnt(12))=iReac;  iR%TROEqf(icnt(12),:)=C 
      CASE ('TROEXP')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
        icnt(13)=icnt(13)+1; iR%iTROExp(icnt(13))=iReac;  iR%TROExp(icnt(13),:)=C 
      CASE ('TROEMCM')
        IF ( SIZE(C)<10 ) CALL ErrorMSG(iReac,Line1)
        icnt(14)=icnt(14)+1; iR%iTROEmcm(icnt(14))=iReac; iR%TROEmcm(icnt(14),:)=C 
      CASE ('SPEC1')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(15)=icnt(15)+1; iR%iSPEC1(icnt(15))=iReac;   iR%SPEC1(icnt(15),:)=C 
      CASE ('SPEC2')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(16)=icnt(16)+1; iR%iSPEC2(icnt(16))=iReac;   iR%SPEC2(icnt(16),:)=C 
      CASE ('SPEC3')
        IF ( SIZE(C)<6 ) CALL ErrorMSG(iReac,Line1)
        icnt(17)=icnt(17)+1; iR%iSPEC3(icnt(17))=iReac;   iR%SPEC3(icnt(17),:)=C 
      CASE ('SPEC4')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
        icnt(18)=icnt(18)+1; iR%iSPEC4(icnt(18))=iReac;   iR%SPEC4(icnt(18),:)=C 
      CASE ('SPEC1MCM')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
        icnt(19)=icnt(19)+1; iR%iSPEC1mcm(icnt(19))=iReac; iR%SPEC1mcm(icnt(19),:)=C 
      CASE ('SPEC2MCM')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
        icnt(20)=icnt(20)+1; iR%iSPEC2mcm(icnt(20))=iReac; iR%SPEC2mcm(icnt(20),:)=C 
      CASE ('SPEC3MCM')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(21)=icnt(21)+1; iR%iSPEC3mcm(icnt(21))=iReac; iR%SPEC3mcm(icnt(21),:)=C 
      CASE ('SPEC4MCM')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
        icnt(22)=icnt(22)+1; iR%iSPEC4mcm(icnt(22))=iReac; iR%SPEC4mcm(icnt(22),:)=C 
      CASE ('SPEC5MCM')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
        icnt(23)=icnt(23)+1; iR%iSPEC5mcm(icnt(23))=iReac; iR%SPEC5mcm(icnt(23),:)=C 
      CASE ('SPEC6MCM')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
        icnt(24)=icnt(24)+1; iR%iSPEC6mcm(icnt(24))=iReac; iR%SPEC6mcm(icnt(24),:)=C 
      CASE ('SPEC7MCM')
        IF ( SIZE(C)<6 ) CALL ErrorMSG(iReac,Line1)
        icnt(25)=icnt(25)+1; iR%iSPEC7mcm(icnt(25))=iReac; iR%SPEC7mcm(icnt(25),:)=C 
      CASE ('SPEC8MCM')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
        icnt(26)=icnt(26)+1; iR%iSPEC8mcm(icnt(26))=iReac; iR%SPEC8mcm(icnt(26),:)=C 
      CASE ('SPEC9MCM')
        IF ( SIZE(C)<10 ) CALL ErrorMSG(iReac,Line1)
        icnt(47)=icnt(47)+1; iR%iSPEC9mcm(icnt(47))=iReac; iR%SPEC9mcm(icnt(47),:)=C 
      CASE ('S4H2O')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
        icnt(27)=icnt(27)+1; iR%iS4H2O(icnt(27))=iReac;  iR%S4H2O(icnt(27),:)=C 
      CASE ('T1H2O')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(28)=icnt(28)+1; iR%iT1H2O(icnt(28))=iReac;  iR%T1H2O(icnt(28),:)=C 
      CASE ('ASPEC1')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(30)=icnt(30)+1; iR%iASPEC1(icnt(30))=iReac; iR%ASPEC1(icnt(30),:)=C 
      CASE ('ASPEC2')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
        icnt(31)=icnt(31)+1; iR%iASPEC2(icnt(31))=iReac; iR%ASPEC2(icnt(31),:)=C 
      CASE ('ASPEC3')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(32)=icnt(32)+1; iR%iASPEC3(icnt(32))=iReac; iR%ASPEC3(icnt(32),:)=C 
      CASE ('DCONST')
        IF ( SIZE(C)<2 ) CALL ErrorMSG(iReac,Line1)
        icnt(34)=icnt(34)+1
        iR%iDCONST(icnt(34),1)=iReac
        iR%iDCONST(icnt(34),2)=iReac+1
        iR%DCONST(icnt(34),:)=C 
      CASE ('DTEMP')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
        icnt(35)=icnt(35)+1; iR%iDTEMP(icnt(35),1)=iReac
        iR%iDTEMP(icnt(35),2)=iReac+1; iR%DTEMP(icnt(35),:)=C 
      CASE ('DTEMP2')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
        icnt(36)=icnt(36)+1; iR%iDTEMP2(icnt(36),1)=iReac
        iR%iDTEMP2(icnt(36),2)=iReac+1 ; iR%DTEMP2(icnt(36),:)=C
      CASE ('DTEMP3')
        IF ( SIZE(C)<4 ) CALL ErrorMSG(iReac,Line1)
        icnt(37)=icnt(37)+1; iR%iDTEMP3(icnt(37),1)=iReac
        iR%iDTEMP3(icnt(37),2)=iReac+1; iR%DTEMP3(icnt(37),:)=C
      CASE ('DTEMP4')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
        icnt(38)=icnt(38)+1; iR%iDTEMP4(icnt(38),1)=iReac
        iR%iDTEMP4(icnt(38),2)=iReac+1; iR%DTEMP4(icnt(38),:)=C
      CASE ('DTEMP5')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
        icnt(39)=icnt(39)+1; iR%iDTEMP5(icnt(39),1)=iReac
        iR%iDTEMP5(icnt(39),2)=iReac+1; iR%DTEMP5(icnt(39),:)=C
      CASE ('PHOTO')
        IF ( SIZE(C)<1 ) CALL ErrorMSG(iReac,Line1)
         icnt(41)=icnt(41)+1; iR%iPHOTOkpp(icnt(41))=iReac; iR%PHOTOkpp(icnt(41))=C(1)
      CASE ('PHOTO2')
        IF ( SIZE(C)<1 ) CALL ErrorMSG(iReac,Line1)
        icnt(42)=icnt(42)+1; iR%iPHOTO2kpp(icnt(42))=iReac; iR%PHOTO2kpp(icnt(42))=C(1)
      CASE ('PHOTO3')
        IF ( SIZE(C)<1 ) CALL ErrorMSG(iReac,Line1)
        icnt(43)=icnt(43)+1; iR%iPHOTO3kpp(icnt(43))=iReac; iR%PHOTO3kpp(icnt(43))=C(1)
      CASE ('SPECIAL')
        icnt(44)=icnt(44)+1; iR%iSPECIAL(icnt(44))=iReac
      CASE ('HOM1')
        IF ( SIZE(C)<3 ) CALL ErrorMSG(iReac,Line1)
        icnt(45)=icnt(45)+1; iR%iHOM1(icnt(45))=iReac; iR%HOM1(icnt(45),:)=C 
      CASE DEFAULT
        WRITE(*,*) ''
        WRITE(*,*) ' Reaction Type unknown:  ',TRIM(TypeR),'  --> check input file'
        WRITE(*,*) ''
    END SELECT

    CONTAINS

      SUBROUTINE ErrorMSG(iR,Line)
        INTEGER      :: iR
        CHARACTER(*) :: Line
        WRITE(*,*);  WRITE(*,*)
        WRITE(*,*) '  ERROR -- > check parameter in reaction: ', iR ,'  ::  '//Line
        WRITE(*,*);  WRITE(*,*)
        STOP
      END SUBROUTINE ErrorMSG

  END SUBROUTINE Setup_ReacParameter
  
  
  SUBROUTINE CheckConstants(RS)
    TYPE(ReactionStruct_T) :: RS(:)     ! reaction system
    CHARACTER(LenType) :: Const_T       ! constant type for reaction i
    !
    INTEGER :: i,j
    !

    DO i = 1,SIZE(RS)
      Const_T = RS(i)%TypeConstant
      DO j=1,nReacTypes
        IF ( TRIM(reac_par(j)%name_type) == ADJUSTL(TRIM(Const_T)) ) EXIT
      END DO
      IF ( SIZE(RS(i)%Constants) /= reac_par(j)%n_par ) THEN
        WRITE(*,*) 'ERROR: Wrong number of constants:'
        WRITE(*,*) '----->  reaction:     ',i, '   ', TRIM(RS(i)%Line1)
        WRITE(*,*) '----->  desired #consts: ', reac_par(j)%n_par, j
        WRITE(*,*) '----->  actual  #consts: ', SIZE(RS(i)%Constants)
        WRITE(*,*) '       Check sys-file for syntax errors!'
        STOP
      END IF

    END DO
  END SUBROUTINE CheckConstants


  SUBROUTINE AllocateRTarrays()

    ! allocate index arrays and parameter arrays for the new vectorized version
    ALLOCATE( iR%iCONST(nr_CONST), iR%CONST(nr_CONST))

    ALLOCATE( iR%iPHOTabc(nr_PHOTabc), iR%iPHOTab(nr_PHOTab), iR%iPHOTmcm(nr_PHOTmcm))
    ALLOCATE( iR%PHOTabc(nr_PHOTabc,3), iR%PHOTab(nr_PHOTab,2), iR%PHOTmcm(nr_PHOTmcm,3))

    ALLOCATE( iR%iTEMP(nr_TEMP),  iR%iTEMP1(nr_TEMP1),  iR%iTEMP2(nr_TEMP2),  iR%iTEMP3(nr_TEMP3),  iR%iTEMP4(nr_TEMP4))
    ALLOCATE( iR%TEMP(nr_TEMP,3), iR%TEMP1(nr_TEMP1,2), iR%TEMP2(nr_TEMP2,2), iR%TEMP3(nr_TEMP3,2), iR%TEMP4(nr_TEMP4,2))

    ALLOCATE( iR%iTROE(nr_TROE),     iR%iTROEf(nr_TROEf),   iR%iTROEq(nr_TROEq),   &
            & iR%iTROEqf(nr_TROEqf), iR%iTROExp(nr_TROExp), iR%iTROEmcm(nr_TROEmcm))
    ALLOCATE( iR%TROE(nr_TROE,4),     iR%TROEf(nr_TROEf,5),   iR%TROEq(nr_TROEq,6),   &
            & iR%TROEqf(nr_TROEqf,7), iR%TROExp(nr_TROExp,5), iR%TROEmcm(nr_TROEmcm,10))

    ALLOCATE( iR%iSPEC1(nr_SPEC1), iR%iSPEC2(nr_SPEC2), &
            & iR%iSPEC3(nr_SPEC3), iR%iSPEC4(nr_SPEC4))
    ALLOCATE( iR%SPEC1(nr_SPEC1,2), iR%SPEC2(nr_SPEC2,2), iR%SPEC3(nr_SPEC3,6), iR%SPEC4(nr_SPEC4,4))

    ALLOCATE( iR%iSPEC1mcm(nr_SPEC1mcm), iR%iSPEC2mcm(nr_SPEC2mcm),&
            & iR%iSPEC3mcm(nr_SPEC3mcm), iR%iSPEC4mcm(nr_SPEC4mcm),&
            & iR%iSPEC5mcm(nr_SPEC5mcm), iR%iSPEC6mcm(nr_SPEC6mcm),&
            & iR%iSPEC7mcm(nr_SPEC7mcm), iR%iSPEC8mcm(nr_SPEC8mcm),&
            & iR%iSPEC9mcm(nr_SPEC9mcm) )
    ALLOCATE( iR%SPEC1mcm(nr_SPEC1mcm,3), iR%SPEC2mcm(nr_SPEC2mcm,3),&
            & iR%SPEC3mcm(nr_SPEC3mcm,2), iR%SPEC4mcm(nr_SPEC4mcm,4),&
            & iR%SPEC5mcm(nr_SPEC5mcm,4), iR%SPEC6mcm(nr_SPEC6mcm,4),&
            & iR%SPEC7mcm(nr_SPEC7mcm,6), iR%SPEC8mcm(nr_SPEC8mcm,4),&
            & iR%SPEC9mcm(nr_SPEC9mcm,10) )

    ALLOCATE( iR%iS4H2O(nr_S4H2O), iR%iT1H2O(nr_T1H2O))
    ALLOCATE( iR%S4H2O(nr_S4H2O,4), iR%T1H2O(nr_T1H2O,2))

    ALLOCATE( iR%iASPEC1(nr_ASPEC1), iR%iASPEC2(nr_ASPEC2),& 
            & iR%iASPEC3(nr_ASPEC3))
    ALLOCATE( iR%ASPEC1(nr_ASPEC1,2), iR%ASPEC2(nr_ASPEC2,3),&
            & iR%ASPEC3(nr_ASPEC3,2) ) 

    ALLOCATE( iR%iDTEMP(nr_DTEMP,2),   iR%iDTEMP2(nr_DTEMP2,2),                           &
            & iR%iDTEMP3(nr_DTEMP3,2), iR%iDTEMP4(nr_DTEMP4,2), iR%iDTEMP5(nr_DTEMP5,2) )
    ALLOCATE( iR%DTEMP(nr_DTEMP,3),    iR%DTEMP2(nr_DTEMP2,4),                           &
            & iR%DTEMP3(nr_DTEMP3,4),  iR%DTEMP4(nr_DTEMP4,4),  iR%DTEMP5(nr_DTEMP5,3) )

    ALLOCATE( iR%iDCONST(nr_DCONST,2) )
    ALLOCATE( iR%DCONST(nr_DCONST,2) )

    ALLOCATE( iR%iFAC_H2(nr_FAC_H2), iR%iFAC_O2N2(nr_FAC_O2N2), iR%iFAC_M(nr_FAC_M),        &
            & iR%iFAC_O2(nr_FAC_O2), iR%iFAC_N2(nr_FAC_N2), iR%iFAC_H2O(nr_FAC_H2O),        &
            & iR%iFAC_RO2(nr_FAC_RO2), iR%iFAC_O2O2(nr_FAC_O2O2), iR%iFAC_aH2O(nr_FAC_aH2O),&
            & iR%iFAC_RO2aq(nr_FAC_RO2aq) )

    ALLOCATE( iR%iHENRY(nr_henry,4) )

    ALLOCATE( iR%iPHOTOkpp(nr_PHOTOkpp), iR%iPHOTO2kpp(nr_PHOTO2kpp), &
            & iR%iPHOTO3kpp(nr_PHOTO3kpp),  iR%PHOTOkpp(nr_PHOTOkpp), &
            & iR%PHOTO2kpp(nr_PHOTO2kpp),  iR%PHOTO3kpp(nr_PHOTO3kpp) )

    ALLOCATE( iR%iSpecial(nr_special) ) 

    ALLOCATE( iR%iHOM1(nr_HOM1) , iR%HOM1(nr_HOM1,3) ) 

  END SUBROUTINE AllocateRTarrays

  SUBROUTINE PrintReaction(iR,Unit)
    INTEGER :: iR
    INTEGER :: Unit

    WRITE(Unit,*) ''
    WRITE(Unit,*) ' ********************************************************************************************'
    WRITE(Unit,*) '  Reaction Number   :: ', iR
    WRITE(Unit,*) '  Reaction Class    :: ', TRIM(ReactionSystem(iR)%Type)
    WRITE(Unit,*) '  Constant Type     :: ', TRIM(ReactionSystem(iR)%TypeConstant)
    WRITE(Unit,*) '  Reaction          :: ', TRIM(ReactionSystem(iR)%Line1)
    WRITE(Unit,*) '  Order of Reaction :: ', INT(SUM(ReactionSystem(iR)%Educt%Koeff))
    WRITE(Unit,*) '  Factor            :: ', TRIM(ReactionSystem(iR)%Factor)
    WRITE(Unit,*) '  Constants         :: ', ReactionSystem(iR)%Constants
    WRITE(Unit,*) ' ********************************************************************************************'
    WRITE(Unit,*) ''
  END SUBROUTINE PrintReaction

  SUBROUTINE make_ChemSys_1_to_nD_arrays()
    INTEGER :: i, cnt
    INTEGER :: nDIM     ! this instance of nDIM unfortunately has to be created here redundantly,
                        ! because the global nDIM is created after we call this subroutine, alternatively,
                        ! the global nDIM could be assigned here, which would be inconsistent with the
                        ! assigning of the other nDIM-like variables, which are assigned afterwards in Cminor.f90

    nDIM = nspc
    IF (combustion) nDIM = nDIM+1
    IF (adiabatic_parcel) nDIM = nDIM+5

    ALLOCATE( nD_spc(nDIM), nD_KAT(ns_KAT), nD_reac(nreac),            &
            & nD_Ptr_spc(nDIM+1), nD_Ptr_KAT(ns_KAT+1), nD_Ptr_reacs(nreac+1) )

    nD_spc   = .FALSE.
    nD_KAT   = .FALSE.
    nD_reac  = .FALSE.

    ! SPECIES ARRAYS
    !
    ! set all aqueous species to vector-valued species
    IF ( ALLOCATED(iAs) ) nD_spc(iAs)     = .TRUE.
    IF ( adiabatic_parcel )   nD_spc(iAqMassEq) = .TRUE.

    ! create Start-/End-Pointer Array to know where all instances of a species are in a
    ! full value array (i.e. in an array where each vector-valued species has a vector-entry)
    cnt=1
    nD_Ptr_spc(1) = 1
    DO i = 1 , nDIM
      IF ( i>=bAs(1) .AND. i<=bAs(2) .AND. ALLOCATED(iAs) ) THEN
        cnt = cnt + nDropletClasses
      ELSE IF ( adiabatic_parcel .AND. i==iAqMassEq ) THEN
        cnt = cnt + nDropletClasses
      ELSE 
        cnt = cnt + 1
      END IF
      nD_Ptr_spc(i+1) = cnt
    END DO

    ! KAT ARRAYS
    !
    cnt=1
    nD_Ptr_KAT(1) = 1
    DO i = 1, ns_KAT
      IF (INDEX(y_name(nspc+i), '[a')>0) THEN
        cnt = cnt + nDropletClasses
        nD_KAT(i) = .TRUE.
      ELSE
        cnt = cnt + 1
      END IF
      nD_Ptr_KAT(i+1) = cnt
    END DO

    ! REAC ARRAYS
    ! 
    cnt = 1
    nD_Ptr_reacs(1) = 1
    DO i = 1 , nreac
      IF ( ReactionSystem(i)%Type == 'HENRY' .OR. &
         & ReactionSystem(i)%Type == 'DISS'  .OR. &
         & ReactionSystem(i)%Type == 'AQUA'       ) THEN
        
        nD_reac(i) = .TRUE.
        cnt = cnt + nDropletClasses
      ELSE
        cnt = cnt + 1
      END IF
      nD_Ptr_reacs(i+1) = cnt
    END DO

  END SUBROUTINE make_ChemSys_1_to_nD_arrays

  SUBROUTINE make_droplet_classes_by_modes(DropletClasses, AFrac, Mode, Ini_AqSpc)
    TYPE(DropletClasses_T), INTENT(OUT) :: DropletClasses
  
    TYPE(AFRAC_T), DIMENSION(:), INTENT(IN) :: AFrac
    TYPE(Modes_T), INTENT(IN) :: Mode
    REAL(dp) :: Ini_AqSpc(ns_AQUA)
        
    REAL(dp) :: qexp, min_r, max_r, tmp1, tmp2, r_step
    REAL(dp), ALLOCATABLE :: radius_bounds(:,:), Number0_m3(:)
    INTEGER :: i, j, k, iSpc, jSpc, n_same_spc, n_present_spc, n_distinct_modes, iMode, iPos, tmp_int
    INTEGER, ALLOCATABLE :: compatible_modes(:), nClassesPerMode(:)

    ! determine compatible modes 
    ! ( compatible modes are modes with same species )
    ! compatible_modes(i)=j means mode i is compatible with new mode j
    ALLOCATE(compatible_modes(nFrac))
    compatible_modes = 0
    ! n_distinct_modes serves as id during the loop, afterwards it is the total number of distinct modes
    n_distinct_modes = 0
    DO i = 1 , nFrac
      IF ( compatible_modes(i)==0 ) THEN         ! mode i is not yet compatible with another one
        n_distinct_modes = n_distinct_modes + 1  ! a new distinct mode has to be created
        compatible_modes(i) = n_distinct_modes   ! mode i is compatible with itself
        DO j = i+1 , nFrac                       ! search for modes compatible with i
          IF (compatible_modes(j)==0 .AND. SIZE(AFrac(i)%Species)==SIZE(AFrac(j)%Species)) THEN
            ! mode i and j may be compatible, compare species
            n_same_spc = 0
            n_present_spc = 0
            DO iSpc = 1 , SIZE(AFrac(i)%Species)
              IF ( AFrac(i)%Frac1(iSpc)>ZERO ) THEN
                n_present_spc = n_present_spc + 1
                DO jSpc = 1 , SIZE(AFrac(j)%Species)
                  ! check if aerosol j contains the same species
                  IF ( TRIM(AFrac(i)%Species(iSpc))==TRIM(AFrac(j)%Species(jSpc)) .AND. &
                     & AFrac(j)%Frac1(jSpc) == AFrac(i)%Frac1(iSpc)             ) THEN
                    n_same_spc = n_same_spc + 1
                    EXIT
                  END IF
                END DO
                ! if species i was not found in mode j, mode i and j are not compatible
                IF (n_same_spc<n_present_spc) EXIT
              END IF
            END DO
            ! if the following holds, we found that mode i and j contain the same species, so they are compatible
            IF (n_same_spc==n_present_spc) compatible_modes(j) = n_distinct_modes
          END IF
        END DO
      END IF
    END DO

    ! determine lower and upper bound of finally considered dry radii of the distinct modes
    ALLOCATE(radius_bounds(2, n_distinct_modes))
    ! exponent for computing p-quantiles (p=0.995)
    ! qexp = inverf(2*p-1) = -inverf(2*(1-p)-1)
    qexp  = 1.821386367718449
    ! p=0.9995
    !qexp = 2.32675

    ! (in the following way we get at least (p-(1-p))*100% of each mode (99% for p=0.995))
    DO iMode = 1 , n_distinct_modes
      min_r = 1.d99
      max_r = -1.d99
      DO i = 1 , nFrac
        IF ( compatible_modes(i) == iMode ) THEN
          ! compute lower and upper p-quantile of mode
          tmp1 = EXP( LOG(Mode%Radius(i)) - SQRT(2*LOG(Mode%Std_deviation(i))**2)*qexp )
          !tmp2 = EXP( LOG(Mode%Radius(i)) + SQRT(2*LOG(Mode%Std_deviation(i))**2)*qexp )
          ! use mass quantile for upper boundary (mass quantile from Feingold and Levin, 1986)
          tmp2 = EXP( LOG(Mode%Radius(i)) + THREE*LOG(Mode%Std_deviation(i))**2 + SQRT(2*LOG(Mode%Std_deviation(i))**2)*qexp )
          ! track the lowest / highest of quantiles
          IF ( tmp1 < min_r ) min_r = tmp1
          IF ( tmp2 > max_r ) max_r = tmp2
        END IF
      END DO
      radius_bounds(1, iMode) = min_r
      radius_bounds(2, iMode) = max_r
    END DO

    ! now distribute the budget of droplet classes
    ! according to: a new mode gets more droplet classes if it has a higher number of old modes
    !               to ensure higher resolution for e.g. bimodal distributions -
    !               if nDropletClasses is not divisible by n_Frac,
    !               give the remaining budget to the modes with lower number of old modes
    ALLOCATE(nClassesPerMode(n_distinct_modes))
    ! determine droplet classes per old mode
    tmp_int = nDropletClasses / nFrac
    nClassesPerMode = 0
    DO i = 1 , nFrac
      ! if old mode i is in new mode, add tmp_int
      nClassesPerMode(compatible_modes(i)) = nClassesPerMode(compatible_modes(i)) + tmp_int
    END DO
    IF ( MODULO(nDropletClasses, nFrac)>0 ) THEN
      DO i = 1, MODULO(nDropletClasses, nFrac)
        nClassesPerMode(i) = nClassesPerMode(i) + 1
      END DO
    END IF

    ALLOCATE( DropletClasses%dryRadius(nDropletClasses)     &
    &       , DropletClasses%dryVolume(nDropletClasses)     &
    &       , DropletClasses%waterMass(nDropletClasses)     &
    &       , DropletClasses%wetRadius(nDropletClasses)     &
    &       , DropletClasses%saturated_air_radius(nDropletClasses)    &
    &       , DropletClasses%critical_RH(nDropletClasses)   &
    &       , DropletClasses%critical_r(nDropletClasses)    &
    &       , DropletClasses%saturated_air_LWC(nDropletClasses) &
    &       , DropletClasses%saturated_air_LSC(nDropletClasses)       &
    &       , DropletClasses%LWC_portion(nDropletClasses)   &
    &       , DropletClasses%active(0)                      &
    &       , DropletClasses%inactive(0)                    &
    &       , DropletClasses%inactive_LWC(nDropletClasses)  &
    &       , DropletClasses%Number(nDropletClasses)        &
    &       , Number0_m3(nDropletClasses)                   &
    &       , DropletClasses%Conc(nDropletClasses, ns_AQUA) )

    j = 1
    DO iMode = 1, n_distinct_modes
      min_r = LOG(radius_bounds(1, iMode))
      max_r = LOG(radius_bounds(2, iMode))

      ! distribute classes in logarithmically equidistant steps
      r_step = (max_r - min_r) / nClassesPerMode(iMode)
        
      DO i = 1, nClassesPerMode(iMode)
        ! radius = middle of logarithmic interval
        DropletClasses%dryRadius(j) = EXP( min_r + rTWO*r_step + (i-1)*r_step )
        DropletClasses%dryVolume(j) = (FOUR/THREE) * Pi * DropletClasses%dryRadius(j)**THREE
        DropletClasses%Number(j)    = 0
        DropletClasses%Conc(j, :) = 1.e-16_dp
        DO k = 1 , nFrac
          IF (compatible_modes(k) == iMode) THEN ! = old mode k is in new mode iMode
            ! the number of particles is CDF(min_r+i*r_step)-CDF(min_r+(i-1)*r_step)
            tmp1 = 0.0
            tmp1 = tmp1 + lognormal_cdf(EXP(min_r+i*r_step)    , Mode%Number(k), Mode%Radius(k), Mode%Std_Deviation(k))   &
                    &         - lognormal_cdf(EXP(min_r+(i-1)*r_step), Mode%Number(k), Mode%Radius(k), Mode%Std_Deviation(k))
            DropletClasses%Number(j) = DropletClasses%Number(j) + tmp1

            ! calculate concentrations by dissolved aerosol
            DO iSpc = 1 , SIZE(AFrac(k)%Species)
              iPos = PositionSpeciesAll(AFrac(k)%Species(iSpc))
              IF ( iPos > 0 ) THEN
                DropletClasses%Conc(j, iPos-bAs(1)+1) = DropletClasses%Conc(j, iPos-bAs(1)+1)            &
                &                                     + 1                                                &  ! [#/m3]
                &                                     * DropletClasses%dryVolume(j)                      &  ! [m3]
                &                                     * Mode%Density(k)                                  &  ! [g/m3]
                &                                     * AFrac(k)%Frac1(iSpc)                             &  ! [g/g]
                &                                     / MolMass(iPos)                                    &  ! 1/[g/mol]
                &                                     * mol2Part                                            ! convert from mol/m3 to molec/cm3
              END IF
            END DO
          END IF
        END DO

        j = j+1
      END DO
    END DO

    Number0_m3 = DropletClasses%Number
    ! convert class numbers to per kg air
    DropletClasses%Number      = DropletClasses%Number / rho_parcel
    DropletClasses%totalNumber = SUM(DropletClasses%Number)

    ! determine classes with sufficient dry radius to be activated (dry radius larger than a given threshold)
    DO i = 1 , nDropletClasses
      ! all droplets smaller than activation_radius will stay on equilibrium radius
      ! based on simplified köhler theory (r_w_eq = sqrt(b/a))
      IF (DropletClasses%dryRadius(i) > activation_radius .OR. adiabatic_parcel) THEN
        DropletClasses%active = [DropletClasses%active, i]
      ELSE
        DropletClasses%inactive = [DropletClasses%inactive, i]
      END IF
    END DO

    ! the following needs to be done before multiplying with the number of droplets
    ! as we need concentrations for one droplet here to determine its eq radius
    DO i = 1 , nDropletClasses
      ! determine equilibrium radius for S=1
      ! according to cheeky method
      CALL calculate_saturated_radius( DropletClasses%saturated_air_radius(i), &
                                     & DropletClasses%saturated_air_LSC(i),    &
                                     & DropletClasses%saturated_air_LWC(i),    &
                                     & DropletClasses%Conc(i,:),     &
                                     & DropletClasses%dryVolume(i),  &
                                     & T_parcel)
    END DO

    ! set inactive radius to .1µm (doesnt matter, all reactions are suppressed)
    DropletClasses%inactive_r   = 0.1e-6_dp
    DropletClasses%inactive_LWC = milli * rho_h2o * Pi43 * DropletClasses%inactive_r**THREE * Number0_m3
    DropletClasses%wetRadius(DropletClasses%inactive) = DropletClasses%inactive_r
    DropletClasses%waterMass(DropletClasses%inactive) = DropletClasses%inactive_LWC(DropletClasses%inactive)

    ! correct lwc for number of droplets
    DropletClasses%saturated_air_LWC = DropletClasses%saturated_air_LWC*Number0_m3
    DropletClasses%saturated_air_LSC = DropletClasses%saturated_air_LSC*Number0_m3

    ! calculate critical S and r (just for info)
    DO i = 1 , nDropletClasses
      CALL calculate_critical_RH_and_r( DropletClasses%critical_RH(i),  &
                                      & DropletClasses%critical_r(i),   &
                                      & DropletClasses%Conc(i, :)*mega, &  ! * mega to come from molec/cm3 to one aerosol sphere!
                                      & T_parcel,                       &
                                      & DropletClasses%dryVolume(i)     )
    END DO

    ! initialization of wet radii
    IF ( adiabatic_parcel ) THEN
      DO i = 1, nDropletClasses
        CALL calculate_RH_radius( DropletClasses%wetRadius(i),    &
                                & DropletClasses%waterMass(i),    &
                                & RH0,                            &
                                & DropletClasses%Conc(i, :)*mega, &  ! * mega to come from molec/cm3 to one aerosol sphere!
                                & T_parcel,                       &
                                & DropletClasses%dryVolume(i)     )
      END DO
      DropletClasses%waterMass = DropletClasses%waterMass * DropletClasses%Number
    ELSE
      ! set LWC portions
      DropletClasses%LWC_portion(DropletClasses%inactive) = ZERO
      DO i = 1 , SIZE(DropletClasses%active)
        DropletClasses%LWC_portion(DropletClasses%active(i)) = &
                  & Number0_m3(DropletClasses%active(i))/SUM(Number0_m3(DropletClasses%active))
      END DO

      DropletClasses%wetRadius = get_wet_radii(LWCs=LWC_array(tBegin))
    END IF

    ! add given initial values (in mol/l)
    DO i = 1, ns_AQUA
      IF ( Ini_AqSpc(i)>2.0e-16_dp ) THEN
        DropletClasses%Conc(:, i) = DropletClasses%Conc(:, i) + Ini_AqSpc(i)*LWC_array(tBegin)*mol2part
      END IF
    END DO

    ! multiply with number of droplets in droplet class
    DO i = 1 , nDropletClasses
      DO j = 1, bAs(2)-bAs(1)+1
        IF (DropletClasses%Conc(i, j)>1.E-15_dp) DropletClasses%Conc(i, j) = Number0_m3(i) * DropletClasses%Conc(i, j)
      END DO
    END DO

    DropletClasses%Number0 = DropletClasses%Number

    IF (DropletClassPrint) THEN
      WRITE(*,*) ''
      WRITE(*,*) '  Droplet class information general:'
      WRITE(*,*) ''
      WRITE(*,*) '    Total Number per m^3: ', DropletClasses%totalNumber*rho_parcel
      WRITE(*,*) '    Total Number per kg : ', DropletClasses%totalNumber
      IF (.NOT. adiabatic_parcel) THEN
        WRITE(*,*) '    Active Classes:   ', DropletClasses%active
        WRITE(*,*) '    Inactive Classes: ', DropletClasses%inactive
      END IF
      DO i = 1 , nDropletClasses
        WRITE(*,*) ''
        WRITE(*,*) '  Droplet class information specific:'
        WRITE(*,*) '    Class ', i
        WRITE(*,*) ''
        WRITE(*,*) '      Dry Radius:        ', DropletClasses%dryRadius(i), ' [m]'
        WRITE(*,*) '      Critical RH:       ', DropletClasses%critical_RH(i), ' [-]'
        WRITE(*,*) '      Critical Radius:   ', DropletClasses%critical_r(i), ' [m]'
        WRITE(*,*) '      Radius at 100% RH: ', DropletClasses%saturated_air_radius(i), ' [m]'
        WRITE(*,*) '      Wet Radius:        ', DropletClasses%wetRadius(i), ' [m]'
        WRITE(*,*) '      Water Mass:        ', DropletClasses%waterMass(i), ' [g/kg]'
        WRITE(*,*) '      Number:            ', DropletClasses%Number(i), ' [#/kg]'
        WRITE(*,*) ''
        WRITE(*,*) '      Concentrations in molec/cm^3:'
        DO j = 1 , ns_AQUA
          IF (DropletClasses%Conc(i,j)>1.e-16_dp) WRITE(*,'(A8,A40,A2,Es10.4)') '        ', TRIM(y_name(iAs(j))), '  ', DropletClasses%Conc(i,j)
        END DO
        WRITE(*,*) ''
      END DO
    END IF

   CONTAINS 

     FUNCTION lognormal_cdf(x, n, m, d) RESULT(F_x)
       ! returns CDF of lognormal distribution
       !
       ! argument
       REAL(dp) :: x
       ! total number, mean, standard deviation
       REAL(dp) :: n, m, d
 
       ! result
       REAL(dp) :: F_x

       F_x = n * rTWO * (ONE + ERF( (LOG(x) - LOG(m)) / (LOG(d)*SQRT(TWO)) ))
     END FUNCTION lognormal_cdf
  END SUBROUTINE make_droplet_classes_by_modes

  PURE SUBROUTINE repeat_values_nD_vec(vec_nD, vec, nD_entries)
    REAL(dp), INTENT(IN)  :: vec(:)
    REAL(dp), INTENT(OUT) :: vec_nD(:)
    INTEGER, INTENT(IN)   :: nD_entries(:)

    INTEGER :: i
           
    DO i = 1, SIZE(vec)
      vec_nD(nD_entries(i):nD_entries(i+1)-1) = vec(i)
    END DO
  END SUBROUTINE repeat_values_nD_vec

END MODULE Chemsys_Mod
