!
! Copyright (C) 2025 Levin Rug, Willi Schimmel (E-Mail: l.rug@lmu.de)
! See ./SRC/Cminor.f90 for the copyright notice
! See ./LICENSE for license information
! SPDX-License-Identifier: GPL-3.0
!
  !=========================================================================!
  !                                                                         !         
  !                                                                         !
  !               Calculating the Rates of the chemical system              ! 
  !                                                                         ! 
  !                                                                         !
  !=========================================================================!
  !
  MODULE Rates_Mod
    !
    USE Kind_Mod,         ONLY: dp

    USE Reac_Mod,         ONLY: nspc, RTind, iR, RTpar, DropletClasses, y_emi, nspc,    &
                              & nspc2, nr_henry, nDIM, nreac, nDIM2, nreac2, iSO, iFO,  &
                              & SwitchTemp, y_depos,     &
                              & nr_PHOTabc, nr_PHOTmcm, nr_TROEmcm, Hp_ind, nr_TROE,    &
                              & nr_TROEq, nr_TROEf, nr_TROEqf, nr_SPEC3, nr_TROExp,     &
                              & nr_SPEC5mcm, nr_SPEC6mcm, nr_SPEC7mcm, nr_SPEC8mcm,     &
                              & nr_SPEC9mcm, nr_ASPEC1, nr_ASPEC2, nr_ASPEC3, nr_CONST, &
                              & nr_DCONST, nr_DTEMP, nr_DTEMP2, nr_DTEMP3, nr_DTEMP4,   &
                              & nr_DTEMP5, nr_HOM1, nr_PHOTAB, H2O_ind, DelGFE, iTeq2,  &
                              & nr_PHOTO2kpp, nr_PHOTO3kpp, nr_PHOTOkpp, iRhoEq2, aHO,  &
                              & nr_S4H2O, nr_SPEC1, nr_SPEC1MCM, nr_SPEC2, nr_SPEC2MCM, &
                              & nr_SPEC4, nr_SPEC4MCM, nr_SPEC3MCM, nr_T1H2O, nr_TEMP,  &
                              & nr_TEMP, nr_TEMP1, nr_TEMP2, nr_TEMP3, nr_TEMP4, iHO,   &
                              & nr_FAC_H2, nr_FAC_H2O, nr_FAC_M, nr_FAC_N2, nr_FAC_O2,  &
                              & nr_FAC_O2N2, nr_FAC_O2O2, nr_FAC_RO2, nr_FAC_RO2aq,     &
                              & RO2, RO2aq, dkmt, nFirst_Order, nFirst_OrderKAT, iAs,   &
                              & nSecond_Order, nHigher_Order, nr_Factor,         &
                              & nr_Special, ns_Aqua, PHOTO, DGFEdT, GFE,      &
                              & dDelGFEdT, henry_diff, gaseous_passive_ind, iFO_Kat,    &
                              & henry_accom, InitValKAT_Ref, InitValKAT, aH2O_ind,      &
                              & iAqMassEq2, nr_FAC_aH2O, iCutOffReacs, OHm_ind,         &
                              & DropletClasses_T, nr_HOaqua, nr_SOaqua, nr_TOaqua
     
    USE Control_Mod,      ONLY: nD_Ptr_spc, adiabatic_parcel, rho_parcel, m_parcel, Pi, &
                              & nD_Ptr_reacs, nD_Ptr_KAT, nD_reac, ZERO, TimeRates,     &
                              & TimerRates, ONE, mega, Temperature0, rTEN, rln10, ln10,  &
                              & rTHREE, Pi34, nDropletClasses, milli, r300, PiHalf,     &
                              & mTHIRTY, EyChiZmin, Dust, nD_KAT, Start_Timer, rFIVE,   &
                              & rTWO, rTWENTY, rTWELV, rSIX, TWO, THREE, End_Timer,     &
                              & rFOUR, FOUR, TEN, kilo, pressure, q_parcel, T_parcel,   &
                              & updraft_velocity, Out
     
    USE Sparse_Mod,       ONLY: BA, TB_sparse, DAX_sparse

    USE Chemsys_Mod,      ONLY: ReactionSystem, sumBAT, ReactionStruct_T,               &
                              & repeat_values_nD_vec
     
    USE Meteo_Mod,        ONLY: mol2part, UpdateSun, Zenith, refTemp, RefM, molw_H2O,   &
                              & aH2O, gasConst_R, InvRefTemp, pseudoLWC, RefPressure,   &
                              & N2, N2O2, H2O, O2, LWC_array, rRcal, rPatm, R, L_v, H2, &
                              & calculate_all_S_eq, rho_H2O, Rv, thermal_conductivity,  &
                              & Cp, g, R_air, esatw, qsatw, Cv, RefH2O, RefN2, RefH2,   &
                              & aH2OmolperL, diff_coeff, molw_air, RefO2, SI_Gas,       &
                              & eps_ionic_strength, dqsatwdT, get_wet_radii, alpha_H2O, &
                              & beta_H2O, dpdz

    USE ChemKinInput_Mod, ONLY: lowA, lowB, lowC, lowD, lowE, lowF, lowG, highA, highB, &
                              & highC, highD, highE, highF, highG
    !
    IMPLICIT NONE
    !
    !
    REAL(dp) :: rFacEq
    REAL(dp), PARAMETER :: dTroe =  0.14_dp

    INTEGER :: RateCnt

    INTERFACE ReactionRates
      MODULE PROCEDURE ReactionRates_Combustion
      MODULE PROCEDURE ReactionRates_Atmosphere
    END INTERFACE ReactionRates
    !
    CONTAINS
    !
    !
    !======================================================================!
    !      Calculate the Rates for current concentraions and time
    !======================================================================!
    SUBROUTINE ReactionRates_Combustion(Y_in,Rate,dKdT_over_k, k_out)
      !--------------------------------------------------------------------!
      ! Input: 
      REAL(dp), INTENT(IN) :: Y_in(nDIM)
      !--------------------------------------------------------------------!
      ! Output:
      REAL(dp), INTENT(OUT) :: Rate(nreac)
      REAL(dp), INTENT(OUT), OPTIONAL :: dKdT_over_k(nreac)
      REAL(dp), INTENT(OUT), OPTIONAL :: k_out(nreac)
      !--------------------------------------------------------------------!
      ! Temporary variables:
      REAL(dp) :: Conc(nspc)
      REAL(dp) :: T(10) ! array for storing Temperature in various ways (clean, reciprocal, roots, ..)
      REAL(dp) :: K(nreac)
      REAL(dp) :: Prod(nreac) 
      REAL(dp) :: Meff(nreac)
      !
      REAL(dp) :: tmpK(nreac), dtmpK(nreac)

      ! temp arrays for chemkin input (temp depended)
      REAL(dp) :: MeffX(RTind%nTBodyExtra)
      REAL(dp) :: F_PD(nreac), DF_PDdT(nreac)
      REAL(dp) :: Dk0dT(nreac), DkinfdT(nreac), F_PDLind(nreac), dF_PDLind_dT(nreac)
      REAL(dp) :: k0(RTind%nLow), kinf(RTind%nHigh), k0M(RTind%nLow)
      REAL(dp) :: rkinfpk0M(RTind%nLow)
      REAL(dp) :: rKeq(RTind%nEqui), DeRdT(RTind%nEqui)
      REAL(dp) :: rRcT
      ! arrays to apply difference quotient for troe factor
      REAL(dp) :: k02(RTind%nLow), kinf2(RTind%nHigh), k0M2(RTind%nLow), F_PDLind2(nreac), T2(10), hT=0.1_dp

      !==================================================================!
      !===============calc rates for ReactionSystem======================!
      !==================================================================!
 
      CALL Start_Timer(TimerRates)

      Rate = ZERO
      tmpK  = ZERO

      IF (PRESENT(dKdT_over_k)) THEN
        dkdT_over_k = ZERO
        dtmpK = ZERO
      END IF

      !--- Concentration
      Conc = Y_in(1:nspc)
      T = UpdateTempArray ( Y_in(nDIM) )

      CALL GibbsFreeEnergie( GFE , T )
      CALL CalcDeltaGibbs  ( DelGFE )
      rFacEq = mega * R * T(1) * rPatm  ! in [cm3/mol]
      !
      IF (PRESENT(dKdT_over_k)) THEN
        CALL CalcDiffGibbsFreeEnergie( DGFEdT , T )
        CALL CalcDiffDeltaGibbs( DDelGFEdT )
      END IF

      !*************************************************************
      ! --- compute rate of all reactions ---
      !*************************************************************
      !
      ! initializing vector components
      Meff = ONE; DF_PDdT = Zero; 

      ! calculate effective molecularity of the reactions
      IF ( RTind%nTBody>0 ) THEN
        Meff( RTind%iTBody ) = SUM( Conc )

        IF ( RTind%nTBodyExtra>0) THEN
          MeffX = DAX_sparse( TB_sparse, Conc )
          Meff(RTind%iTBodyExtra) = Meff(RTind%iTBodyExtra) - MeffX
        END IF

      END IF

      ! ====== Compute the rate constant for specific reaction type
      !
      rRcT = rRcal*T(6)

      ! normal Arrhenius
      IF ( RTind%nArr>0 )  THEN
        tmpK(RTind%iArr)  = RTpar%A*EXP(RTpar%b*T(8) - RTpar%E*rRcT)
        IF (PRESENT(dKdT_over_k)) dtmpK(RTind%iArr) = (RTpar%b + RTpar%E*rRcT)*T(6)
      END IF

      ! Backward reaction with explicitly given Arrhenius parameters
      IF ( RTind%nXrev>0 ) THEN
        tmpK(RTind%iXrev) = RTpar%AX*EXP(RTpar%bX*T(8) - RTpar%EX*rRcT)
        IF (PRESENT(dKdT_over_k)) dtmpK(RTind%iXrev) = (RTpar%bX + RTpar%EX*rRcT)*T(6)
      END IF

      ! reduced pressure value for falloff reactions
      F_PD = ONE;        dF_PDLind_dT = ONE
      Dk0dT = ZERO;      DkinfdT  = ZERO

      IF ( RTind%nLow>0 ) THEN
        ! Low pressure Arrhenius
        k0    = RTpar%A0*EXP(RTpar%b0*T(8) - RTpar%E0*rRcT)
        ! High pressure Arrhenius
        kinf    = RTpar%Ainf*EXP(RTpar%binf*T(8) - RTpar%Einf*rRcT)

        k0M = k0 * Meff(RTind%iLow)
        F_PDLind = ONE
        F_PDLind(RTind%iLow) = k0M / (kinf+k0M)

        tmpK(RTind%iLow)   = kinf

        IF (PRESENT(dKdT_over_k)) THEN
          rkinfpk0M = ONE / (kinf + k0M)
          Dk0dT(RTind%iLow)    = (RTpar%b0   + RTpar%E0  *rRcT) * T(6)
          DkinfdT(RTind%iHigh) = (RTpar%binf + RTpar%Einf*rRcT) * T(6)

          dF_PDLind_dT(RTind%iLow) = kinf*rkinfpk0M * (Dk0dT(RTind%iLow) - DkinfdT(RTind%iHigh))

          dtmpK(RTind%iHigh) = DkinfdT(RTind%iHigh)
        END IF
      END IF

      ! reverse reaction, calculating the equi constant 
      IF (RTind%nEqui>0) THEN
        ! the '-1's are to use the values of the forward reaction (see Perini et al. 2012)
        rKeq = EXP(DelGFE(RTind%iEqui-1)) * rFacEq**(sumBAT(RTind%iEqui-1))
        tmpK(RTind%iEqui)  = tmpK(RTind%iEqui) * rKeq

        IF (PRESENT(dKdT_over_k)) THEN
          DeRdT = sumBAT(RTind%iEqui-1)*T(6) + DDelGFEdT(RTind%iEqui-1) 
          dtmpK(RTind%iEqui) = dtmpK(RTind%iEqui) + DeRdT
        END IF
      END IF

      ! rate constants according to Lindemann's form
      IF ( RTind%nLind>0 ) THEN
        F_PD(RTind%iLind)    = F_PDLind(RTind%iLind)
        IF (PRESENT(dKdT_over_k)) DF_PDdT(RTind%iLind) = dF_PDLind_dT(RTind%iLind)
        Meff(RTind%iLind)   = ONE
      END IF

      ! rate constants according to Troe's form
      IF ( RTind%nTroe>0 ) THEN
        F_PD(RTind%iTroe)    = TroeFactorVec(T, F_PDLind)*F_PDLind(RTind%iTroe)
        IF (PRESENT(dKdT_over_k)) THEN
          T2 = UpdateTempArray ( Y_in(nDIM) + hT )
          ! Low pressure Arrhenius
          k02 = RTpar%A0*EXP(RTpar%b0*T2(8) - RTpar%E0*rRcal*T2(6))
          ! High pressure Arrhenius
          kinf2 = RTpar%Ainf*EXP(RTpar%binf*T2(8) - RTpar%Einf*rRcal*T2(6))
          k0M2 = k02 * k0M / k0
          F_PDLind2 = ONE
          F_PDLind2(RTind%iLow) = k0M2 / (kinf2+k0M2)

          DF_PDdT(RTind%iTroe) = dF_PDLind_dT(RTind%iTroe) + (TroeFactorVec(T2, F_PDLind2)-TroeFactorVec(T, F_PDLind))/hT
        END IF
        Meff(RTind%iTroe)   = ONE
      END IF

      k = F_PD * tmpK
      IF (PRESENT(dKdT_over_k)) dkdT_over_k = DF_PDdT + dtmpK

      ! ==== Law of mass action productories
      Prod = MassActionProducts( Conc )

      Rate  = Meff * k * Prod

      IF (PRESENT(k_out)) k_out = k

      Out%nRateEvals = Out%nRateEvals + 1
      CALL End_Timer(TimerRates, TimeRates)
    END SUBROUTINE ReactionRates_Combustion

    !======================================================================!
    !      Calculate the Rates for current concentraions and time
    !======================================================================!
    SUBROUTINE ReactionRates_Atmosphere(Time,Y_in,Rate)
      USE fparser, ONLY: evalf
      !--------------------------------------------------------------------!
      ! Input: 
      REAL(dp), INTENT(IN) :: Time
      REAL(dp), INTENT(IN) :: Y_in(nDIM2)
      !--------------------------------------------------------------------!
      ! Output:
      REAL(dp), INTENT(OUT) :: Rate(nreac2)
      !--------------------------------------------------------------------!
      ! Temporary variables:
      REAL(dp) :: Conc(nspc2)
      REAL(dp) :: chi(3), LWCs(nDropletClasses)
      REAL(dp) :: T(10), Temp
      REAL(dp) :: Meff(nreac2), k(nreac2) , Prod(nreac2)
      REAL(dp) :: mAir, passive_conversion,                &
                &   aqua_unit_conversion(nDropletClasses), &
                & SOaqua_unit_conversion(nDropletClasses), &
                & TOaqua_unit_conversion(nDropletClasses)

      REAL(dp) :: tHenry(nr_HENRY,2,nDropletClasses)
      INTEGER  :: j,j_nD,i
      !==================================================================!
      !===============calc rates for ReactionSystem======================!
      !==================================================================!

      CALL Start_Timer(TimerRates)

      Rate = ZERO

      ! --- Compute zenith for photo reactions
      IF ( PHOTO ) THEN
        ! tropos syntax calculation of zenith
        chi(1) = Zenith(Time)
        !WRITE(*,*) ' chi = ', chi
        chi(2) = ONE/COS(chi(1))
        ! kkp syntax calculation of sun for photo reactions
        chi(3) = UpdateSun(Time)
      END IF

      Conc = Y_in(1:nspc2)
      Temp = T_parcel
      T = UpdateTempArray( Temp )
      
      !*************************************************************
      ! --- compute rate of all reactions (gas,henry,aqua,diss) ---
      !*************************************************************

      ! compute aqueous parameters if needed
      IF ( ns_AQUA>0 .OR. nr_FAC_aH2O>0) THEN 
        LWCs = LWC_array(Time)
        DropletClasses%wetRadius = get_wet_radii(LWCs=LWCs)
      END IF

      ! ====== Update passive species
      ! adjust according to ideal gas law and current pressure/temperature
      passive_conversion = RefTemp * T(6) * Pressure / RefPressure
      InitValKat(nD_Ptr_KAT(gaseous_passive_ind)) = InitValKat_Ref(nD_Ptr_KAT(gaseous_passive_ind)) * passive_conversion
      N2   = RefN2 * passive_conversion
      O2   = RefO2 * passive_conversion
      H2   = RefH2 * passive_conversion
      mAir = RefM  * passive_conversion
      IF ( adiabatic_parcel ) THEN
        H2O = q_parcel / molw_H2O / ( 1 / molw_air + q_parcel / molw_H2O ) &
        &   * rho_parcel*kilo/molw_air * mol2part
        IF (H2O_ind>0) InitValKat(H2O_ind) = H2O
        !H2O = q_parcel * rho_parcel*kilo / molw_H2O
      ELSE
        H2O = RefH2O * passive_conversion
      END IF
      N2O2 = N2 + O2

      ! ====== Computing effective molecularity 
      Meff = ONE
      IF ( nr_FACTOR > 0 ) Meff = EffectiveMolecularity( Conc , mAir )

      ! ====== Compute the rate constant for specific reaction type
      k = ComputeRateConstant( T, chi, mAir, Conc, LWCs )

      ! ====== Compute mass transfer coefficient 
      IF ( nr_HENRY > 0 ) THEN
        thenry = MassTransfer( k(nD_Ptr_reacs(iR%iHENRY(:,1))), T, LWCs )
        DO i = 1 , nDropletClasses
          k(nD_Ptr_reacs(iR%iHENRY(:,1))+i-1) = thenry(:,1,i)
          k(nD_Ptr_reacs(iR%iHENRY(:,3))+i-1) = thenry(:,2,i)
        END DO
      END IF

      ! ====== Special Reactions
      IF ( nr_special > 0 ) THEN 
        DO i = 1,nr_SPECIAL
          ! this looks weird, but assuming a special reaction to have either only scalar or nD species
          ! it should work, because as (if) j_nD grows, the indices in Conc have to grow equally
          j = iR%iSPECIAL(i)
          DO j_nD = nD_Ptr_reacs(j), nD_Ptr_reacs(j+1)-1
            IF (ReactionSystem(j)%Special%Temp) THEN
              k(j_nD) = evalf(i, [ Conc(nD_Ptr_spc(ReactionSystem(j)%Special%iVariables)+j_nD-nD_Ptr_reacs(j)),T(1)])
            ELSE
              k(j_nD) = evalf(i, Conc(nD_Ptr_spc(ReactionSystem(j)%Special%iVariables)+j_nD-nD_Ptr_reacs(j)))
            END IF
          END DO
        END DO
      END IF

      ! ====== Correct unit of concentrations for higher order aqueous reactions
      IF ( ns_AQUA > 0 ) THEN

        InitValKat(nD_Ptr_KAT(aH2O_ind):nD_Ptr_KAT(aH2O_ind+1)-1) = aH2O*LWCs
        
        IF (nr_SOaqua>0 .OR. nr_TOaqua>0 .OR. nr_HOaqua>0) aqua_unit_conversion   = 1/(LWCs*mol2part)
        IF (nr_SOaqua>0)                                   SOaqua_unit_conversion = aqua_unit_conversion
        IF (nr_TOaqua>0)                                   TOaqua_unit_conversion = aqua_unit_conversion*aqua_unit_conversion

        ! iR%HOaqua is the order of the reaction minus one, which results in dividing all educts (which are like LWC*c_aq*mol2part
        ! with c_aq in mol/l) by LWC*mol2part except one, this results in having left one c_aq*mol2part at the rhs,
        ! i.e. also at the lhs, leaving us with the derivative of not c_aq but c_aq*LWC*mol2part, which is what we want,
        ! since we store this, not c_aq.
        ! meanwhile, by dividing the rhs we respect the units of the reaction constants which want mol/l concentrations
        ! NOTE: 
        ! this can be done in two ways, ordered by droplet classes (Version 1) or by reactions (Version 2)
        ! for sulfur oxidation, Version 2 is about 20 times faster
        !
        ! Version 1:
        !DO i = 1 , nDropletClasses
        !  !k(nD_Ptr_reacs(iR%iHOaqua)+i-1) = k(nD_Ptr_reacs(iR%iHOaqua)+i-1) / (LWCs(i)*mol2part)**iR%HOaqua
        ! IF (nr_SOaqua>0) k(nD_Ptr_reacs(iR%iSOaqua)+i-1) = k(nD_Ptr_reacs(iR%iSOaqua)+i-1) * SOaqua_unit_conversion(i)
        ! IF (nr_TOaqua>0) k(nD_Ptr_reacs(iR%iTOaqua)+i-1) = k(nD_Ptr_reacs(iR%iTOaqua)+i-1) * TOaqua_unit_conversion(i)
        ! IF (nr_HOaqua>0) k(nD_Ptr_reacs(iR%iHOaqua)+i-1) = k(nD_Ptr_reacs(iR%iHOaqua)+i-1) *   aqua_unit_conversion(i)**iR%HOaqua
        !
        !  ! cut off reactions if summed up solute conc (ionic strength) exceeds threshold
        !  !IF ( (SUM(Conc(nD_Ptr_Spc(iAs)+i-1)) - Conc(nD_Ptr_spc(Hp_ind)+i-1) - Conc(nD_Ptr_spc(OHm_ind)+i-1))/mol2part/LWCs(i) > eps_ionic_strength ) THEN
        !  !  k(nD_Ptr_Reacs(iCutOffReacs)+i-1) = ZERO
        !  !  WRITE(*,*) 'WARNING :: cutting off reactions in haze droplets can cause numeric fatality!'
        !  !END IF
        !
        !END DO
        ! Version 1 end
        !
        ! Version 2:
        DO i=1, nr_SOaqua
          k(nD_Ptr_reacs(iR%iSOaqua(i)):nD_Ptr_reacs(iR%iSOaqua(i)+1)-1) = k(nD_Ptr_reacs(iR%iSOaqua(i)):nD_Ptr_reacs(iR%iSOaqua(i)+1)-1) * SOaqua_unit_conversion
        END DO
        DO i=1, nr_TOaqua
          k(nD_Ptr_reacs(iR%iTOaqua(i)):nD_Ptr_reacs(iR%iTOaqua(i)+1)-1) = k(nD_Ptr_reacs(iR%iTOaqua(i)):nD_Ptr_reacs(iR%iTOaqua(i)+1)-1) * TOaqua_unit_conversion
        END DO
        DO i=1, nr_HOaqua
          k(nD_Ptr_reacs(iR%iHOaqua(i)):nD_Ptr_reacs(iR%iHOaqua(i)+1)-1) = k(nD_Ptr_reacs(iR%iHOaqua(i)):nD_Ptr_reacs(iR%iHOaqua(i)+1)-1) *   aqua_unit_conversion**iR%HOaqua(i)
        END DO
        ! Version 2 end


        ! set everything to zero in inactive droplet classes
        IF (.NOT. adiabatic_parcel) THEN
          DO i=1,SIZE(DropletClasses%inactive)
            k(nD_Ptr_Reacs(iCutOffReacs)+DropletClasses%inactive(i)-1) = ZERO
          END DO
        END IF
      END IF

      ! ==== Law of mass action productories
      Prod = MassActionProducts( Conc )

      Rate = Meff * k * Prod
      RateCnt = RateCnt + 1

      !CALL Debug_Rates(ReactionSystem,Time,Meff,k,Prod,Rate)

      Out%nRateEvals = Out%nRateEvals + 1
      CALL End_Timer(TimerRates, TimeRates)

    END SUBROUTINE ReactionRates_Atmosphere

    PURE FUNCTION MassActionProducts(Conc) RESULT(Prod)
      REAL(dp)             :: Prod(nreac2)
      REAL(dp), INTENT(IN) :: Conc(nspc2)
      INTEGER :: i

      Prod = ONE
      !
      ! stoechometric coefficients equal 1
      DO i=1,nFirst_order
        Prod(iFO(i,1)) = Prod(iFO(i,1)) * Conc(iFO(i,2))
      END DO
      !
      ! stoechometric coefficients equal 2
      DO i=1,nSecond_order
        Prod(iSO(i,1)) = Prod(iSO(i,1)) * Conc(iSO(i,2)) * ABS(Conc(iSO(i,2)))
      END DO
      !
      ! stoechometric coefficients not equal 1 or 2
      DO i=1,nHigher_order
        Prod(iHO(i,1)) = Prod(iHO(i,1)) * Conc(iHO(i,2)) ** aHO(i)
      END DO
      !
      ! if there are passive (katalytic) species e.g. [N2], [O2] or [aH2O]
      DO i=1,nfirst_orderKAT
        IF (nD_reac(iFO_KAT(i,1))) THEN
          IF (nD_KAT(iFO_KAT(i,2))) THEN
            Prod(nD_Ptr_reacs(iFO_kat(i,1)):nD_Ptr_reacs(iFO_kat(i,1)+1)-1) = Prod(nD_Ptr_reacs(iFO_kat(i,1)):nD_Ptr_reacs(iFO_kat(i,1)+1)-1) &
                                                                          & * InitValKat(nD_Ptr_KAT(iFO_kat(i,2)):nD_Ptr_KAT(iFO_kat(i,2)+1)-1)
          ELSE
            ! e.g. Henry [O2] = aO2
            Prod(nD_Ptr_reacs(iFO_kat(i,1)):nD_Ptr_reacs(iFO_kat(i,1)+1)-1) = Prod(nD_Ptr_reacs(iFO_kat(i,1)):nD_Ptr_reacs(iFO_kat(i,1)+1)-1) &
                                                                          & * InitValKat(nD_Ptr_KAT(iFO_kat(i,2)))
          END IF
        ELSE
          Prod(nD_Ptr_reacs(iFO_kat(i,1))) = Prod(nD_Ptr_reacs(iFO_kat(i,1))) * InitValKat(nD_Ptr_KAT(iFO_kat(i,2)))
        END IF
      END DO

    END FUNCTION MassActionProducts

    PURE FUNCTION MassTransfer(kin,T,LWCs) RESULT(k)
      !REAL(dp)             :: k(nr_HENRY,2,nFrac)
      REAL(dp)             :: k(nr_HENRY,2,nDropletClasses)
      REAL(dp), INTENT(IN) :: kin(nr_HENRY)
      REAL(dp), INTENT(IN) :: T(:), LWCs(nDropletClasses)
      ! TEMP
      !REAL(dp) :: kmt(nr_HENRY,nFrac)
      REAL(dp) :: kmt(nr_HENRY,nDropletClasses)
      REAL(dp) :: term_diff(nr_HENRY), term_accom(nr_HENRY)
      !REAL(dp) :: r_wet(nFrac)
      !REAL(dp) :: r_wet2(nFrac)
      REAL(dp) :: r_wet(nDropletClasses)
      REAL(dp) :: r_wet2(nDropletClasses)
      INTEGER :: i
      !
      !---------------------------------------------------------------------------
      term_diff  = henry_diff(  iR%iHENRY(:,2) )          ! diffusion term
      term_accom = henry_accom( iR%iHENRY(:,2) ) * T(10)  ! accom term
      !--------------------------------------------------------------------------!

      r_wet  = DropletClasses%wetRadius
      r_wet2 = r_wet * r_wet

      !--  mass transfer coefficient
      kmt = dkmt  ! set minimal transfer coefficient
      DO  i = 1,nr_HENRY 
        IF (term_diff(i) /= ZERO) THEN
          kmt(i,:) = term_diff(i)*r_wet2 + term_accom(i)*r_wet
        END IF
      END DO 
      kmt = ONE / kmt


      ! units: multiply eq. (3) from Jaruga 2018 (which has units mol/l) with LWC*mol2part for change in molec/cm3
      ! direction aq->gas:
      !   needs no unit conversion
      ! direction gas->aq:
      !   is then kmt*LWC*mol2part*c_inf_l where c_inf_l is in mol/l, therefore convert
      !                     molec/cm3 -> mol/m3  -> mol/dm3 = mol/l
      !           c_inf_l = c_inf_ours / mol2part * milli
      !   we end up with kmt * LWC * mol2part * c_inf_l = kmt * LWC * c_inf_ours * milli

      ! direction GasSpecies-->AquaSpecies
      !k(:,1,:) = milli * kmt(:,:) * LWC
      DO i = 1 , nDropletClasses
        k(:,1,i) = milli * kmt(:,i) * LWCs(i)
      END DO

      ! direction AquaSpecies-->GasSpecies
      !DO i=1,nFrac
      DO i=1,nDropletClasses
        k(:,2,i) = kmt(:,i) / (kin(:) * GasConst_R * T(1))  ! (...) = HenryConst*GasConstant[l*atm/mol/K]*Temperatur
      END DO

    END FUNCTION MassTransfer

    PURE FUNCTION EffectiveMolecularity(Conc,mAir) RESULT(M)
      !OUT
      REAL(dp) :: M(nreac2)
      !IN
      REAL(dp), INTENT(IN) :: Conc(:)
      REAL(dp), INTENT(IN) :: mAir
      !
      INTEGER :: iDr

      M = ONE

      IF(nr_FAC_M>0)     M(nD_Ptr_reacs(iR%iFAC_M))    = mAir
      IF(nr_FAC_RO2>0)   M(nD_Ptr_reacs(iR%iFAC_RO2))  = SUM(Conc(nD_Ptr_spc(RO2)))

      IF(nr_FAC_O2>0)    M(nD_Ptr_reacs(iR%iFAC_O2))   = O2
      IF(nr_FAC_O2O2>0)  M(nD_Ptr_reacs(iR%iFAC_O2O2)) = O2*O2
      IF(nr_FAC_N2>0)    M(nD_Ptr_reacs(iR%iFAC_N2))   = N2
      IF(nr_FAC_O2N2>0)  M(nD_Ptr_reacs(iR%iFAC_O2N2)) = O2*N2
      IF(nr_FAC_H2O>0)   M(nD_Ptr_reacs(iR%iFAC_H2O))  = H2O
      IF(nr_FAC_H2>0)    M(nD_Ptr_reacs(iR%iFAC_H2))   = H2


      ! ---------------------
      ! -- aqueous factors --
      ! aqueous reactions need either mol/l input or unit correction (done above in ReactionRates_Atmosphere)
      ! RO2aq are marked as higher order while reading the mechanisma and experience unit correction
      IF (nr_FAC_RO2aq>0) THEN
        DO iDr = 1 , nDropletClasses
          M(nD_Ptr_reacs(iR%iFAC_RO2aq)+iDr-1) = SUM(Conc(nD_Ptr_spc(RO2aq)+iDr-1))
        END DO
      END IF
      ! aH2O has constant concentration and is directly given as mol/l value
      IF (nr_FAC_aH2O>0) THEN
        DO iDr = 1 , nDropletClasses
          M(nD_Ptr_reacs(iR%iFAC_aH2O)+iDr-1) = aH2OmolperL
        END DO
      END IF
      ! ---------------------

    END FUNCTION EffectiveMolecularity

    PURE FUNCTION ComputeRateConstant(T,chi,mAir,Conc,LWCs) RESULT(k_nD)

      REAL(dp)                :: k(nreac)
      REAL(dp)                :: k_nD(nreac2)
      REAL(dp), INTENT(IN)    :: mAir, chi(:)
      REAL(dp), INTENT(IN)    :: T(:)
      REAL(dp), INTENT(IN)    :: Conc(nspc2)
      REAL(dp), INTENT(IN)    :: LWCs(nDropletClasses)

      ! Photoabc tempo parameter
      REAL(dp), DIMENSION(nr_PHOTabc) :: ChiZabc, yChiZabc, EyChiZabc
      REAL(dp), DIMENSION(nr_PHOTmcm) :: ChiZmcm, yChiZmcm
      INTEGER :: i
      ! Aspec tempo parameter
      REAL(dp), PARAMETER   :: x = 13.0d0
      ! Troe tempo parameters
      REAL(dp), DIMENSION(nr_TROE)    :: k1, k2, log10_k1k2
      REAL(dp), DIMENSION(nr_TROEq)   :: k1q, k2q, k3q, log10_k1k2q
      REAL(dp), DIMENSION(nr_TROEf)   :: k1f, k2f, log10_k1k2f
      REAL(dp), DIMENSION(nr_TROEqf)  :: k1qf, k2qf, k3qf, log10_k1k2qf
      REAL(dp), DIMENSION(nr_TROExp)  :: k1xp, k2xp, log10_k1k2xp
      REAL(dp), DIMENSION(nr_TROEmcm) :: k1mcm, k2mcm, Fc, tmpTROE
      REAL(dp), PARAMETER   :: n = 0.75_dp, d = 1.27_dp
      REAL(dp) :: F, Tr300, Hp_molperl(nDropletClasses)
      ! Spec tempo parameters
      REAL(dp), DIMENSION(nr_SPEC3)    :: k1s, k2s ,k3s
      REAL(dp), DIMENSION(nr_SPEC5mcm) :: k1smcm5, k2smcm5
      REAL(dp), DIMENSION(nr_SPEC6mcm) :: k1smcm6, k2smcm6
      REAL(dp), DIMENSION(nr_SPEC7mcm) :: k1smcm7, k2smcm7
      REAL(dp), DIMENSION(nr_SPEC8mcm) :: k1smcm8, k2smcm8
      REAL(dp), DIMENSION(nr_SPEC9mcm) :: k1smcm9, k2smcm9, k3smcm9

      k    = ZERO
      k_nD = ZERO
      
      !--------------------------------------------------------------------------!
      !---  Photolysis
      !--------------------------------------------------------------------------!
      !
      !---  Roeth
      IF (nr_PHOTAB>0) THEN
        IF ( chi(1) < PiHalf ) THEN
          k(iR%iPHOTab) = Dust * iR%PHOTab(:,1)*EXP(-iR%PHOTab(:,2)*chi(2))
        END IF
      END IF

      !---  Derwent and Hough (1988)
      IF (nr_PHOTabc>0) THEN
        IF ( chi(1) < PiHalf ) THEN
          ChiZabc = chi(1) * iR%PHOTabc(:,3) 
          DO i = 1,nr_PHOTabc
            IF (ChiZabc(i) < PiHalf) THEN
              yChiZabc(i) = iR%PHOTabc(i,2) * (One - One/COS(ChiZabc(i)))
              IF ( yChiZabc(i) > mTHIRTY ) THEN
                EyChiZabc(i) = EXP(yChiZabc(i))
              ELSE
                EyChiZabc(i) = EyChiZmin   ! = 9.357d-14  
              END IF
            ELSE
              EyChiZabc(i) = EyChiZmin   ! = 9.357d-14 
            END IF
          END DO
          k(iR%iPHOTabc) = Dust * iR%PHOTabc(:,1) * EyChizabc
        END IF
      END IF

      !---  MCM version
      IF (nr_PHOTMCM>0) THEN
        IF ( chi(1) < PiHalf ) THEN
          ChiZmcm  = EXP( -iR%PHOTmcm(:,3) * chi(2) )
          yChiZmcm = chi(1) ** iR%PHOTmcm(:,2)
          k(iR%iPHOTmcm) = Dust * iR%PHOTmcm(:,1) * yChiZmcm * ChiZmcm
        END IF
      END IF

      !--------------------------------------------------------------------------!
      !---  Constant
      !--------------------------------------------------------------------------!
      !
      IF (nr_CONST>0) k(iR%iCONST) = iR%CONST

      !--------------------------------------------------------------------------!
      !---  Temperature-Dependent  (Arrhenius)
      !--------------------------------------------------------------------------!
      !
      IF (nr_TEMP>0) THEN
        k(iR%iTEMP)  = iR%TEMP(:,1) * T(1)**iR%TEMP(:,2) * EXP(-iR%TEMP(:,3)*T(6))
      END IF
      IF (nr_TEMP1>0) THEN
        k(iR%iTEMP1) = iR%TEMP1(:,1)*EXP(-iR%TEMP1(:,2)*T(6))
      END IF
      IF (nr_TEMP2>0) THEN
        k(iR%iTEMP2) = iR%TEMP2(:,1)*T(2)*EXP(-iR%TEMP2(:,2)*T(6))
      END IF
      IF (nr_TEMP3>0) THEN
        k(iR%iTEMP3) = iR%TEMP3(:,1)*EXP(iR%TEMP3(:,2)*(T(6)-InvRefTemp))
      END IF
      IF (nr_TEMP4>0) THEN
        k(iR%iTEMP4) = iR%TEMP4(:,1)*T(1)*EXP(-iR%TEMP4(:,2)*T(6))
      END IF

      !--------------------------------------------------------------------------!
      !---  Dissociation reactions (equilibrium reactions)
      !--------------------------------------------------------------------------!
      !
      !---  constant
      IF (nr_DCONST>0) THEN
        k(iR%iDCONST(:,2)) = iR%DCONST(:,2) ! backward reactions
        k(iR%iDCONST(:,1)) = iR%DCONST(:,1) * k(iR%iDCONST(:,2)) ! forward reactions
      END IF

      !---  temperatur-dependent
      IF (nr_DTEMP>0)  THEN 
        k(iR%iDTEMP(:,2)) = iR%DTEMP(:,3) ! backward reactions
        k(iR%iDTEMP(:,1)) = iR%DTEMP(:,1)*EXP(iR%DTEMP(:,2)*(T(6)-InvRefTemp)) * k(iR%iDTEMP(:,2)) ! forward reactions
      END IF

      !---  temperatur-dependent forward and backward reaction
      IF (nr_DTEMP2>0) THEN 
        k(iR%iDTEMP2(:,2)) = iR%DTEMP2(:,3)*EXP(iR%DTEMP2(:,4)*(T(6)-InvRefTemp)) ! backward reactions
        k(iR%iDTEMP2(:,1)) = iR%DTEMP2(:,1)*EXP(iR%DTEMP2(:,2)*(T(6)-InvRefTemp)) * k(iR%iDTEMP2(:,2))! forward reactions
      END IF
      
      !--- Jacobson(1999, Table B7, p. 603-605)
      IF (nr_DTEMP3>0) THEN 
        k(iR%iDTEMP3(:,2)) = iR%DTEMP3(:,4) ! backward reactions
        k(iR%iDTEMP3(:,1)) = iR%DTEMP3(:,1) * EXP(iR%DTEMP3(:,2)*(RefTemp*T(6)-ONE)         &
        &      + iR%DTEMP3(:,3)*(One-RefTemp*T(6) + LOG10(RefTemp*T(6)))) * k(iR%iDTEMP3(:,2)) ! forward reactions
      END IF

      !---  temperatur-dependent forward and backward reaction
      IF (nr_DTEMP4>0) THEN 
        k(iR%iDTEMP4(:,2)) = ONE ! backward reactions
        k(iR%iDTEMP4(:,1)) = iR%DTEMP4(:,1)*EXP(iR%DTEMP4(:,2)        &
        &      * (T(1)*InvRefTemp-One) + iR%DTEMP4(:,3)            &
        &      * (ONE+LOG10(T(1)*InvRefTemp)-T(1)*InvRefTemp)) * k(iR%iDTEMP4(:,2)) ! forward reactions
      END IF

      !---  temperatur-dependent forward and backward reaction
      IF (nr_DTEMP5>0) THEN 
        k(iR%iDTEMP5(:,2)) = ONE ! backward reactions
        k(iR%iDTEMP5(:,1)) = iR%DTEMP5(:,1)*(T*InvRefTemp)**iR%DTEMP4(:,2)  &
        &      * EXP(iR%DTEMP5(:,3)*(T(6)-InvRefTemp)) * k(iR%iDTEMP5(:,2))! forward reactions
      END IF

      !--------------------------------------------------------------------------!
      !---  Humidity-Dependent
      !--------------------------------------------------------------------------!
      !
      !---  for reaction O1D = HO+HO 
      IF (nr_T1H2O>0) k(iR%iT1H2O) = iR%T1H2O(:,1)*H2O*EXP(-iR%T1H2O(:,2)*T(6))
      
      !---  for reaction HO2+HO2 = H2O2
      IF (nr_S4H2O>0) THEN 
        k(iR%iS4H2O) = iR%S4H2O(:,1)*EXP(iR%S4H2O(:,2)*T(6))      & 
        &            + iR%S4H2O(:,3)*EXP(iR%S4H2O(:,4)*T(6))*mAir
      END IF
     
      !--------------------------------------------------------------------------!
      !---  Special reactions for Gas Phase Chemistry 
      !--------------------------------------------------------------------------!
      !
      Tr300 = T(1)*r300

      !---  Troe
      IF (nr_TROE>0) THEN   
        F = 0.6_dp
        k1 = iR%TROE(:,1) * (Tr300)**(-iR%TROE(:,2))*mAir
        k2 = iR%TROE(:,3) * (Tr300)**(-iR%TROE(:,4))
        log10_k1k2  = LOG10(k1/k2)
        k(iR%iTROE) = k1/(One+k1/k2)*F**(One/(One+log10_k1k2*log10_k1k2))
      END IF

      !---  Troe equilibrium
      IF (nr_TROEQ>0) THEN  
        F = 0.6_dp
        k1q = iR%TROEq(:,1)*(Tr300)**(-iR%TROEq(:,2))*mAir
        k2q = iR%TROEq(:,3)*(Tr300)**(-iR%TROEq(:,4))
        log10_k1k2q = LOG10(k1q/k2q)
        k3q = k1q/(One+k1q/k2q)*F**(One/(One+log10_k1k2q*log10_k1k2q))
        k(iR%iTROEq) = k3q/(iR%TROEq(:,5)*EXP(iR%TROEq(:,6)*T(6)))
      END IF
      
      !---  Troe with variable F factor
      IF (nr_TROEF>0) THEN  
        k1f = iR%TROEf(:,1)*(Tr300)**(-iR%TROEf(:,2))*mAir
        k2f = iR%TROEf(:,3)*(Tr300)**(-iR%TROEf(:,4))
        log10_k1k2f = LOG10(k1f/k2f)
        k(iR%iTROEf) = k1f/(One+k1f/k2f)*iR%TROEf(:,5)**(One/(One+log10_k1k2f*log10_k1k2f))
      END IF
      
      !---  Troe equilibrium with variable F factor
      IF (nr_TROEQF>0) THEN 
        k1qf = iR%TROEqf(:,1) * (Tr300)**(-iR%TROEqf(:,2))*mAir
        k2qf = iR%TROEqf(:,3) * (Tr300)**(-iR%TROEqf(:,4))
        log10_k1k2qf = LOG10(k1qf/k2qf)
        k3qf = k1qf/(One+k1qf/k2qf)*iR%TROEqf(:,7)**(One/(One+log10_k1k2qf*log10_k1k2qf))
        k(iR%iTROEqf) = k3qf/(iR%TROEqf(:,5)*EXP(iR%TROEqf(:,6)*T(6)))          
      END IF
      
      !---  modified Troe with variable F factor
      IF (nr_TROEXP>0) THEN 
        k1xp = iR%TROExp(:,1)*EXP(-iR%TROExp(:,2)*T(6))*mAir
        k2xp = iR%TROExp(:,3)*EXP(-iR%TROExp(:,4)*T(6))
        log10_k1k2xp  = LOG10(k1xp/k2xp)
        k(iR%iTROExp) = k1xp/(One+k1xp/k2xp)*iR%TROExp(:,5)**(One/(One+log10_k1k2xp*log10_k1k2xp)) 
      END IF
      
      !---  MCM-reaction TROE
      IF (nr_TROEMCM>0) THEN 
        k1mcm = iR%TROEmcm(:,1)*(T(1)*InvRefTemp)**iR%TROEmcm(:,2)*EXP(iR%TROEmcm(:,3)*T(6))*mAir
        k2mcm = iR%TROEmcm(:,4)*(T(1)*InvRefTemp)**iR%TROEmcm(:,5)*EXP(iR%TROEmcm(:,6)*T(6))
        Fc    = iR%TROEmcm(:,7)*EXP(iR%TROEmcm(:,8)*T(6))+iR%TROEmcm(:,9)*EXP(T(1)/iR%TROEmcm(:,10))
        tmpTROE = LOG10(k1mcm/k2mcm)/(n-d*LOG10(Fc))
        k(iR%iTROEmcm) = k1mcm/(One+k1mcm/k2mcm)*Fc**(One/(One+tmpTROE*tmpTROE))
      END IF

      !--------------------------------------------------------------------------!
      !---  Special Types  (Gas Phase: Density-Dependent in RADM/RACM)
      !--------------------------------------------------------------------------!
      !
      !---  (CO+HO = HO2+CO2)
      IF (nr_SPEC1>0) k(iR%iSPEC1) = iR%SPEC1(:,1)*(ONE+mAir*iR%SPEC1(:,2))

      !---  (O3PX + [O2] = O3)
      IF (nr_SPEC2>0) k(iR%iSPEC2) = iR%SPEC2(:,1)*(Tr300)**iR%SPEC2(:,2)*mAir
      
      !---  (HNO3+HO = NO3)
      IF (nr_SPEC3>0) THEN 
        k1s = iR%SPEC3(:,1)*EXP(iR%SPEC3(:,2)*T(6))
        k2s = iR%SPEC3(:,3)*EXP(iR%SPEC3(:,4)*T(6))
        k3s = iR%SPEC3(:,5)*EXP(iR%SPEC3(:,6)*T(6))*mAir
        k(iR%iSPEC3) = k1s+k3s/(One+k3s/k2s)
      END IF

      !---  (HO2+HO2 = H2O2)
      IF (nr_SPEC4>0) THEN 
        k(iR%iSPEC4) = iR%SPEC4(:,1)*EXP(iR%SPEC4(:,2)*T(6))+iR%SPEC4(:,3)*EXP(iR%SPEC4(:,4)*T(6))*mAir 
      END IF

      !---  SPEC1 for MCM
      IF (nr_SPEC1MCM>0) THEN 
        k(iR%iSPEC1mcm) = iR%SPEC1mcm(:,1)*(ONE+mAir*iR%SPEC1mcm(:,2)*Tr300/iR%SPEC1mcm(:,3))
      END IF

      !---  SPEC2 for MCM (HO2 + O3 = OH)
      IF (nr_SPEC2MCM>0) THEN 
        k(iR%iSPEC2mcm) = iR%SPEC2mcm(:,1)*(Tr300)**iR%SPEC2mcm(:,2)*EXP(iR%SPEC2mcm(:,3)*T(6))
      END IF

      !---  OH + CO = HO2
      IF (nr_SPEC3MCM>0) THEN 
        k(iR%iSPEC3mcm) = iR%SPEC3mcm(:,1)*(ONE+mAir/iR%SPEC3mcm(:,2))
      END IF

      !---  SPEC4 for MCM
      IF (nr_SPEC4MCM>0) THEN 
        k(iR%iSPEC4mcm) = iR%SPEC4mcm(:,1)*(ONE+iR%SPEC4mcm(:,2) &
        &               * EXP(iR%SPEC4mcm(:,3)*T(6))*H2O)*EXP(iR%SPEC4mcm(:,4)*T(6))
      END IF

      !---  (SPEC5 for MCM)
      IF (nr_SPEC5MCM>0) THEN 
        F = 0.21_dp
        k1smcm5 = iR%SPEC5mcm(:,1)*mAir*F*EXP(iR%SPEC5mcm(:,2)*T(6))
        k2smcm5 = iR%SPEC5mcm(:,3)*mAir*F*EXP(iR%SPEC5mcm(:,4)*T(6))
        k(iR%iSPEC5mcm) = k1smcm5*(One-k2smcm5)
      END IF

      !---  (SPEC6 for MCM)
      IF (nr_SPEC6MCM>0) THEN 
        k1smcm6 = iR%SPEC6mcm(:,1)*EXP(iR%SPEC6mcm(:,2)*T(6))
        k2smcm6 = iR%SPEC6mcm(:,3)*EXP(iR%SPEC6mcm(:,4)*T(6))
        k(iR%iSPEC6mcm) = k1smcm6*(ONE-k2smcm6)
      END IF

      !---  (SPEC7 for MCM)
      IF (nr_SPEC7MCM>0) THEN 
        k1smcm7 = iR%SPEC7mcm(:,1)*EXP(iR%SPEC7mcm(:,2)*T(6))
        k2smcm7 = iR%SPEC7mcm(:,3)*EXP(iR%SPEC7mcm(:,4)*T(6))
        k(iR%iSPEC7mcm) = k1smcm7*(iR%SPEC7mcm(:,5)-iR%SPEC7mcm(:,6)/(One+k2smcm7))
      END IF

      !---  (SPEC8 for MCM)
      IF (nr_SPEC8MCM>0) THEN 
        F = 0.21_dp
        k1smcm8 = iR%SPEC8mcm(:,1)*mAir*F*EXP(iR%SPEC8mcm(:,2)*T(6))
        k2smcm8 = iR%SPEC8mcm(:,3)*mAir*F*EXP(iR%SPEC8mcm(:,4)*T(6))
        k(iR%iSPEC8mcm) = k1smcm8/(One+k2smcm8)*T(6)
      END IF

      !---  (SPEC9 for MCM)
      IF (nr_SPEC9MCM>0) THEN 
        F = 0.21_dp
        k1smcm9 = iR%SPEC9mcm(:,1)*mAir*F*EXP(iR%SPEC8mcm(:,2)*T(6)) &
        &       / (ONE  + iR%SPEC9mcm(:,3)*mAir*F*EXP(iR%SPEC8mcm(:,4)*T(6)))
        k2smcm9 = iR%SPEC9mcm(:,5)*mAir*F*EXP(iR%SPEC9mcm(:,6)*T(6)) &
        &       / ((ONE + iR%SPEC9mcm(:,7)*mAir*F*EXP(iR%SPEC9mcm(:,8)*T(6)))*T(1))
        k3smcm9 = iR%SPEC9mcm(:,9)*EXP(iR%SPEC9mcm(:,10)*T(6))
        k(iR%iSPEC9mcm) =  (k1smcm9 * k3smcm9) / (k2smcm9 + k3smcm9)
      END IF

      IF (nr_HOM1>0) THEN
        k(iR%iHOM1) = iR%HOM1(:,1) * EXP(iR%HOM1(:,2)/T(1)) * EXP(iR%HOM1(:,3)/(T(3)))
      END IF

      ! ************************************************************************
      ! *** KPP photolytic reactions
      ! ************************************************************************
      !
      IF (nr_PHOTOkpp>0)  THEN 
        k(iR%iPHOTOkpp)  = iR%PHOTOkpp(:) * Chi(3)
      END IF
      IF (nr_PHOTO2kpp>0) THEN 
        k(iR%iPHOTO2kpp) = iR%PHOTO2kpp(:) * Chi(3)*Chi(3)
      END IF
      IF (nr_PHOTO3kpp>0) THEN 
        k(iR%iPHOTO3kpp) = iR%PHOTO3kpp(:) * Chi(3)*Chi(3)*Chi(3)
      END IF 
      

      ! ====== repeat values of aqueous rate constants since they are equal in all droplet classes
      CALL repeat_values_nD_vec(k_nD, k, nD_Ptr_reacs)


      ! ====== fill droplet class specific rate constants

      !--------------------------------------------------------------------------!
      !---  Special reactions for aqueous phase chemistry
      !--------------------------------------------------------------------------!
      !
      !---  pH dependent
      IF (nr_ASPEC1>0 .OR. nr_ASPEC2>0 .OR. nr_ASPEC3>0) THEN
        Hp_molperl = Conc(nD_Ptr_spc(Hp_ind):nD_Ptr_spc(Hp_ind+1)-1) / mol2part / LWCs
      END IF

      IF (nr_ASPEC1>0) THEN
        DO i = 1, nDropletClasses
          k_nD(nD_Ptr_reacs(iR%iASPEC1)+i-1) = Hp_molperl(i)*(iR%ASPEC1(:,1)*EXP(iR%ASPEC1(:,2)       & 
          &                                    * (T(6)-InvRefTemp)))/(ONE + x*Hp_molperl(i))
        END DO
      END IF
      IF (nr_ASPEC2>0) THEN 
        DO i = 1, nDropletClasses
          k_nD(nD_Ptr_reacs(iR%iASPEC2)+i-1) = Hp_molperl(i)**iR%ASPEC2(:,2)                          &
          &                                    * (iR%ASPEC2(:,1)*EXP(iR%ASPEC2(:,3)*(T(6)-InvRefTemp)))
        END DO
      END IF
      IF (nr_ASPEC3>0) THEN 
        DO i = 1, nDropletClasses
          k_nD(nD_Ptr_reacs(iR%iASPEC3)+i-1) = iR%ASPEC3(:,1)*EXP(iR%ASPEC3(:,2)*(-LOG10(Hp_molperl(i))))
        END DO
      END IF

    END FUNCTION ComputeRateConstant



    SUBROUTINE UpdateEmission(Emissions, Conc)
      REAL(dp), INTENT(OUT) :: Emissions(nspc2), Conc(nspc2)

      ! emissions and deposition combined
      Emissions = y_emi - y_depos*Conc

    END SUBROUTINE UpdateEmission

    !   ***************************************************************
    !   Fitting polynomials for combustion mechanisms
    !   ***************************************************************
    PURE SUBROUTINE GibbsFreeEnergie(Gibbs,T)
    !   ** Species nondimensional gibbs potentials                   **
      REAL(dp), INTENT(INOUT) :: Gibbs(:)
      REAL(dp), INTENT(IN)    :: T(:)
      !
      Gibbs(:)=ZERO
      ! WHERE wird bald abgeschaft, -> vektorisieren
      WHERE (SwitchTemp>T(1))
        Gibbs = lowA*(ONE-T(8)) - rTWO*lowB*T(1) - rSIX*lowC*T(2)        &
        &        - rTWELV*lowD*T(3) - rTWENTY*lowE*T(4) + lowF*T(6) - lowG 
      ELSEWHERE
        Gibbs = highA*(ONE-T(8)) - rTWO*highB*T(1) - rSIX*highC*T(2)         &
        &        - rTWELV*highD*T(3) - rTWENTY*highE*T(4) + highF*T(6) - highG 
      END WHERE
    END SUBROUTINE GibbsFreeEnergie
    !
    !
    PURE SUBROUTINE CalcDiffGibbsFreeEnergie(DGibbsdT,T)
      REAL(dp), INTENT(INOUT) :: DGibbsdT(:)
      REAL(dp), INTENT(IN)    :: T(:)
      !
      DGibbsdT(:)=ZERO
      ! WHERE wird bald abgeschaft, -> vektorisieren
      WHERE (SwitchTemp>T(1))
        DGibbsdT=-(lowA*T(6) + 0.5d0*lowB + lowC*T(1)/3.0d0 +            &
        &              0.25d0*lowD*T(2) + 0.2d0*lowE*T(3) + lowF*T(7) )
      ELSEWHERE
        DGibbsdT=-(highA*T(6) + 0.5d0*highB + highC*T(1)/3.0d0 +         &
        &              0.25d0*highD*T(2) + 0.2d0*highE*T(3) + highF*T(7) )
      END WHERE  
    END SUBROUTINE CalcDiffGibbsFreeEnergie
    !
    !
    PURE SUBROUTINE CalcDeltaGibbs(DelGibbs)
      REAL(dp), INTENT(INOUT) :: DelGibbs(:)
      !
      INTEGER :: iR
      INTEGER :: from, to
      !
      DelGibbs = ZERO         
      !
      DO iR=1,nreac
        from = BA%RowPtr(iR);   to = BA%RowPtr(iR+1)-1
        DelGibbs(iR) = DelGibbs(iR)   &
        &            + SUM( BA%Val(from:to) * GFE(BA%ColInd(from:to)) )
      END DO

    END SUBROUTINE CalcDeltaGibbs
    

    PURE SUBROUTINE CalcDiffDeltaGibbs(DiffDelGibbs)
      REAL(dp), INTENT(INOUT) :: DiffDelGibbs(:)
      !
      INTEGER :: iR, jS, jj
      !
      DiffDelGibbs(:)=ZERO         
      DO iR=1,BA%m
        DO jj=BA%RowPtr(iR),BA%RowPtr(iR+1)-1
          jS=BA%ColInd(jj)
          DiffDelGibbs(iR)=DiffDelGibbs(iR)+BA%Val(jj)*DGFEdT(jS)      
        END DO
      END DO
    END SUBROUTINE CalcDiffDeltaGibbs
    
    
    SUBROUTINE scTHERMO(C,H,S,T)
      !
      ! IN
      REAL(dp) :: T(8)
      !
      ! OUT
      REAL(dp) :: C(nspc)       ! molar heat capacities at constant pressure
      REAL(dp) :: H(nspc)       ! the standardstate molar enthalpy
      REAL(dp) :: S(nspc)       ! standard-state entropy at 298 K
      !
      REAL(dp) :: dHdT(nspc)    ! Enthaply derivative in dT [J/mol/K^2]
      REAL(dp) :: dGdT(nspc)    ! Gibbs potential derivative in dT [J/mol/K^2]
      REAL(dp) :: dCvdT(nspc)   ! Constant volume specific heat derivative in dT [J/mol/K]
      
      WHERE (SwitchTemp>T(1))
        C = lowA + lowB*T(1) + lowC*T(2) + lowD*T(3) + lowE*T(4)
        H = lowA + rTWO*lowB*T(1) + rTHREE*lowC*T(2) + rFOUR*lowD*T(3)    &
        &        + rFIVE*lowE*T(4) + lowF*T(6)
        S = lowA*T(8) + lowB*T(1) + rTWO*lowC*T(2) + rTHREE*lowD*T(3)     &
        &        + rFIVE*lowE*T(4) + lowG
        !
        dCvdT = R * (highB + TWO*highC*T(1) + THREE*highD*T(2) + FOUR*highE*T(3))
        dHdT = R * (highA + highB*T(1) + highC*T(2) + highD*T(3) + highE*T(4))
        dGdT = - R * (highG + highA*T(8) + highB*T(1) + rTWO*highC*T(2)      &
        &                   + rTHREE*highD*T(3) +rFOUR*highE*T(4))
      ELSEWHERE
        C = highA + highB*T(1) + highC*T(2) + highD*T(3) + highE*T(4)
        H = highA + rTWO*highB*T(1) + rTHREE*highC*T(2) + rFOUR*highD*T(3)    &
        &         + rFIVE*highE*T(4) + highF*T(6)
        S = highA*T(8) + highB*T(1) + rTWO*highC*T(2) + rTHREE*highD*T(3)     &
        &         + rFIVE*highE*T(4) + highG
        !
        dCvdT = R * (highB + TWO*highC*T(1) + THREE*highD*T(2) + FOUR*highE*T(3))
        dHdT = R * (highA + highB*T(1) + highC*T(2) + highD*T(3) + highE*T(4))
        dGdT = - R * (highG + highA*T(8) + highB*T(1) + rTWO*highC*T(2)      &
        &                   + rTHREE*highD*T(3) +rFOUR*highE*T(4))
      END WHERE  
    END SUBROUTINE scTHERMO

    FUNCTION TroeFactorVec(T, F_PDLind)
      !OUT
      REAL(dp) :: TroeFactorVec(RTind%nTroe)
      !IN
      REAL(dp) :: T(:), F_PDLind(nreac)

      !TEMP
      REAL(dp) :: Fcent(RTind%nTroe),       &
                & logF_Troe(RTind%nTroe),   &
                & FACtroe(RTind%nTroe),     &
                & log10_Pr(RTind%nTroe),    &
                & log10_Fcent(RTind%nTroe), &
                & cTroe(RTind%nTroe),       &
                & n1Troe(RTind%nTroe)

      Fcent = (ONE - RTpar%T1) * EXP(-T(1)/RTpar%T2) &
            &      + RTpar%T1  * EXP(-T(1)/RTpar%T3) &
            &      +             EXP(-T(6)*RTpar%T4)

      log10_Fcent = LOG10(Fcent)
      log10_Pr    = LOG10(F_PDLind(RTind%iTroe)/(ONE-F_PDLind(RTind%iTroe)))

      cTroe  = -0.40e0_dp - 0.67e0_dp*log10_Fcent
      n1Troe =  0.75e0_dp - 1.27e0_dp*log10_Fcent
      
      FACtroe = (log10_Pr+cTroe) / (n1Troe-dTroe*(log10_Pr+cTroe))

      logF_Troe  = log10_Fcent / (ONE+FACtroe*FACtroe)
      TroeFactorVec = TEN**logF_Troe

    END FUNCTION TroeFactorVec

    PURE FUNCTION UpdateTempArray(Temperature) RESULT(TempArr)
      REAL(dp), INTENT(IN) :: Temperature 
      REAL(dp) :: TempArr(10)
      !
      INTEGER :: i
      !
      TempArr(1) = Temperature                 ! T
      DO i=2,5
        TempArr(i) = TempArr(i-1)*Temperature  ! T^2 ... T^5
      END DO
      TempArr(6)  = ONE / Temperature          ! 1/T
      TempArr(7)  = ONE / TempArr(2)           ! 1/T^2
      TempArr(8)  = LOG(Temperature)           ! ln(T)
      TempArr(9)  = SQRT(Temperature)          ! sqrt(T)
      TempArr(10) = ONE / TempArr(9)           ! 1/sqrt(T)
    END FUNCTION UpdateTempArray
    
    !-------------------------------------------------------------------------
    !---  Species internal energies in moles [J/mol]  
    !-------------------------------------------------------------------------
    PURE SUBROUTINE InternalEnergy(U,T)
      REAL(dp), INTENT(INOUT) :: U(nspc)   
      REAL(dp), INTENT(IN)    :: T(:)
      
      WHERE (SwitchTemp>T(1))
        U  = (  (lowA-ONE)*T(1) + rTWO*lowB*T(2)  + rTHREE*lowC*T(3)    & 
        &     + rFOUR*lowD*T(4) + rFIVE*lowE*T(5) + lowF )
      ELSEWHERE
        U  = (  (highA-ONE)*T(1) + rTWO*highB*T(2)  + rTHREE*highC*T(3) & 
        &     + rFOUR*highD*T(4) + rFIVE*highE*T(5) + highF )
      END WHERE
      U    = U * R
    END SUBROUTINE InternalEnergy
    !
    !-------------------------------------------------------------------------------
    !--- Nondimensionl Derivatives of specific heats at constant volume in moles [-]
    !-------------------------------------------------------------------------------
    PURE SUBROUTINE DiffInternalEnergy(dUdT,T)
      REAL(dp), INTENT(INOUT) :: dUdT(nspc)   
      REAL(dp), INTENT(IN)    :: T(:)                      
      !
      WHERE (SwitchTemp>T(1))
        dUdT  = (lowA + lowB*T(1) + lowC*T(2)     & 
        &             + lowD*T(3) + lowE*T(4)) 
      ELSEWHERE
        dUdT  = (highA + highB*T(1) + highC*T(2)  & 
        &              + highD*T(3) + highE*T(4))  
      END WHERE
      dUdT = (dUdT - ONE)        ! speedchem SCmodule.f ~2618
    END SUBROUTINE DiffInternalEnergy
    
    
    PURE SUBROUTINE Diff2InternalEnergy(d2UdT2,T)
      REAL(dp), INTENT(INOUT) :: d2UdT2(nspc)     !Constant volume specific heat’s derivative [J/mol/K2] 
      REAL(dp), INTENT(IN)    :: T(:)                      
      !
      WHERE (SwitchTemp>T(1))
        d2UdT2 = lowB + TWO*lowC*T(1) + THREE*lowD*T(2) + FOUR*lowE*T(3) 
      ELSEWHERE    
        d2UdT2 = highB + TWO*highC*T(1) + THREE*highD*T(2) + FOUR*highE*T(3)
      END WHERE
      d2UdT2 = d2UdT2 * R
    END SUBROUTINE Diff2InternalEnergy


  ! calculates change of density due to ideal gas law and instantaneous adjustment to ambient pressure
  !FUNCTION rhs_rho(rho, T, dzdt, p, dTdt, dqdt, q, z) RESULT(drhodt)
  !  REAL(dp) :: rho, T, dzdt, p, dTdt, dqdt, q, z
  FUNCTION rhs_rho(T, dzdt, p, dTdt, z) RESULT(drhodt)
    REAL(dp) :: T, dzdt, p, dTdt, z

    REAL(dp) :: drhodt !, Tv_fac

    !Tv_fac = 1+0.61*q

    ! moist air
    !drhodt = dpdz(z) * dzdt / R_air / (T*Tv_fac) &
    !       & - dTdt * p / R_air / (T**2 * Tv_fac) - 0.61 * p / (R_air*T*Tv_fac**2) * dqdt
    ! dry air
    !drhodt = - rho * g * dzdt / R_air / T - dTdt * p / R_air / T**2
    drhodt = dpdz(z) * dzdt / R_air / T - dTdt * p / R_air / T**2
    !drhodt = ZERO

  END FUNCTION rhs_rho

  ! calculates rate of temperature change due to adiabatically rising parcel (updraft velocity w) and condensation dmdt
  !FUNCTION rhs_T_cond_and_parcel(dmdt, dzdt, q) RESULT (dTdt)
  !  REAL(dp) :: dmdt(nDropletClasses), dzdt, q
  FUNCTION rhs_T_cond_and_parcel(dmdt, dzdt) RESULT (dTdt)
    REAL(dp) :: dmdt(nDropletClasses), dzdt

    REAL(dp) :: dTdt !, Cp_m

    !Cp_m = Cp_moist(q)

    dTdt = - g * dzdt / Cp + L_v * milli * SUM(dmdt) / Cp

  END FUNCTION rhs_T_cond_and_parcel

  FUNCTION rhs_q_parcel(dmdt) RESULT(dqdt)
    REAL(dp) :: dmdt(nDropletClasses)

    REAL(dp) :: dqdt

    ! dmdt is in g/kg/s -> convert to kg/kg/s
    dqdt = - milli * SUM(dmdt)

  END FUNCTION rhs_q_parcel

  SUBROUTINE rhs_condensation(dmdt, DropletClasses, Y, RH, dSeqdmw_out)
    TYPE(DropletClasses_T), INTENT(IN)  :: DropletClasses
    REAL(dp)              , INTENT(IN)  :: Y(nDIM2), RH

    REAL(dp), DIMENSION(nDropletClasses), INTENT(OUT) :: dmdt
    REAL(dp), DIMENSION(nDropletClasses), OPTIONAL, INTENT(OUT) :: dSeqdmw_out

    REAL(dp), DIMENSION(nDropletClasses) :: m_w, dSeqdm

    REAL(dp) :: F_k, F_d, e_s, r_alpha, r_beta

    REAL(dp), DIMENSION(nDropletClasses) :: G_pre, S_eq, r, n_s
    INTEGER :: i

    ! saturation vapor pressure
    e_s = esatw(T_parcel)

    ! water mass in g
    m_w = DropletClasses%waterMass / DropletClasses%Number

    ! radius of droplets
    r = DropletClasses%wetRadius

    DO i = 1, nDropletClasses
      ! sum concentrations scaled to one droplet (mega for cm3->m3)
      n_s(i) = SUM(Y(nD_Ptr_spc(iAs)+i-1)) * mega / DropletClasses%Number(i)
    END DO

    IF ( PRESENT(dSeqdmw_out) ) THEN
      CALL calculate_all_S_eq(S_eq, m_w, r, n_s, T_parcel, dSeqdm)
    ELSE
      CALL calculate_all_S_eq(S_eq, m_w, r, n_s, T_parcel)
    END IF

    ! relaxation radius Bohrer
    r_alpha = SQRT(2*pi / Rv / T_parcel) * diff_coeff(T_parcel, pressure) / alpha_H2O
    r_beta  = thermal_conductivity(T_parcel) / (beta_H2O * Pressure) * SQRT(2*Pi*R_air*T_parcel) / (Cv+R_air/2)
    ! factors for heat conduction and vapor diffusion
    F_k = L_v / ( thermal_conductivity(T_parcel) * T_parcel ) * &
    &     ( L_v / ( Rv * T_parcel ) - 1.0 )

    F_d = Rv * T_parcel / ( e_s * diff_coeff(T_parcel, pressure) )

    G_pre = 1.0 / ( F_k*S_eq * (r+r_beta) + F_d * (r+r_alpha) )

    ! Köhler derivative term in mass form
    ! (kilo for unit g/s/kg, /kg comes from DropletClasses%Number, which is in #/kg)
    dmdt = DropletClasses%Number * kilo * 4*Pi*r**2 * G_pre * (RH - S_eq)

    IF ( PRESENT(dSeqdmw_out) ) THEN
      dSeqdmw_out = kilo * 4*Pi*r**2 * G_pre * (-dSeqdm)
    END IF

  END SUBROUTINE rhs_condensation

  FUNCTION rhs_z() RESULT(dzdt)
    REAL(dp) :: dzdt

    dzdt = updraft_velocity
  END FUNCTION rhs_z

 END MODULE Rates_Mod
  
