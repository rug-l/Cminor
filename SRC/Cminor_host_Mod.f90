!
! Copyright (C) 2025 Levin Rug, Willi Schimmel (E-Mail: l.rug@lmu.de)
! See ./SRC/Cminor.f90 for the copyright notice
! See ./LICENSE for license information
! SPDX-License-Identifier: GPL-3.0
!
! Mode 1 ICON/MUSCAT host API: chem RHS f and sparse Jacobian J values.
! Prereq: mechanism + symbolic A/BAT/Jac_CC already built (e.g. cminor_initialize
! through Jac_CC). Do not call Integrate / Miter / LU here.
! f = BAT * Rate (no Emiss). J = Jacobian_CC* (pure df/dy).
!
MODULE Cminor_host_Mod

  USE Kind_Mod, ONLY: dp

  USE Control_Mod, ONLY: T_parcel, Pressure, nDropletClasses, nD_spc, eps, ONE, &
                       & combustion

  USE Reac_Mod, ONLY: nspc2, nreac2

  USE Rates_Mod, ONLY: ReactionRates

  USE Sparse_Mod, ONLY: A, BAT, Jac_CC, OPERATOR(*), MULT_BAT_Rate_ValPtr, &
                      & Jacobian_CC, Jacobian_CC_ValPtr, CSR_2_Empty_ValPtr

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cminor_host_set_state, cminor_host_f, cminor_host_jac, &
          & cminor_host_f_and_jac, cminor_host_jac_nnz

  REAL(dp), ALLOCATABLE :: Rate_work(:), Y_work(:)

CONTAINS

  SUBROUTINE ensure_work()
    IF (.NOT. ALLOCATED(Rate_work)) THEN
      IF (nreac2 <= 0 .OR. nspc2 <= 0) STOP 'cminor_host: call init before f/J'
      ALLOCATE(Rate_work(nreac2), Y_work(nspc2))
    END IF
  END SUBROUTINE ensure_work

  SUBROUTINE ensure_valptr()
    IF (nDropletClasses > 1 .AND. .NOT. Jac_CC%vector_entried) THEN
      CALL CSR_2_Empty_ValPtr(Jac_CC, nD_spc, nDropletClasses)
    END IF
  END SUBROUTINE ensure_valptr

  SUBROUTINE rates_at(t, y)
    REAL(dp), INTENT(IN) :: t
    REAL(dp), INTENT(IN) :: y(:)

    CALL ensure_work()
    IF (combustion) THEN
      CALL ReactionRates(y, Rate_work)
    ELSE
      CALL ReactionRates(t, y, Rate_work)
    END IF
  END SUBROUTINE rates_at

  SUBROUTINE rhs_from_rate(f)
    REAL(dp), INTENT(OUT) :: f(:)

    IF (SIZE(f) < nspc2) STOP 'cminor_host_f: f too small'
    IF (nDropletClasses > 1) THEN
      f(1:nspc2) = MULT_BAT_Rate_ValPtr(BAT, Rate_work)
    ELSE
      f(1:nspc2) = BAT * Rate_work
    END IF
  END SUBROUTINE rhs_from_rate

  SUBROUTINE jac_from_rate(y)
    REAL(dp), INTENT(IN) :: y(:)

    CALL ensure_valptr()
    Y_work = MAX(ABS(y(1:nspc2)), eps) * SIGN(ONE, y(1:nspc2))
    IF (nDropletClasses > 1) THEN
      CALL Jacobian_CC_ValPtr(Jac_CC, BAT, A, Rate_work, Y_work)
    ELSE
      CALL Jacobian_CC(Jac_CC, BAT, A, Rate_work, Y_work)
    END IF
  END SUBROUTINE jac_from_rate

  !--- set meteo globals read by ReactionRates_Atmosphere -----------------
  SUBROUTINE cminor_host_set_state(Temp, Press)
    REAL(dp), INTENT(IN) :: Temp
    REAL(dp), INTENT(IN) :: Press

    T_parcel = Temp
    Pressure = Press
  END SUBROUTINE cminor_host_set_state

  INTEGER FUNCTION cminor_host_jac_nnz()
    CALL ensure_valptr()
    IF (.NOT. ALLOCATED(Jac_CC%Val)) STOP 'cminor_host_jac_nnz: Jac_CC not built'
    cminor_host_jac_nnz = SIZE(Jac_CC%Val)
  END FUNCTION cminor_host_jac_nnz

  !--- f = BAT * R(y); no emissions --------------------------------------
  SUBROUTINE cminor_host_f(t, y, f)
    REAL(dp), INTENT(IN)  :: t
    REAL(dp), INTENT(IN)  :: y(:)
    REAL(dp), INTENT(OUT) :: f(:)

    CALL rates_at(t, y)
    CALL rhs_from_rate(f)
  END SUBROUTINE cminor_host_f

  !--- sparse Jac values into jval (CSR order of Jac_CC%Val) --------------
  SUBROUTINE cminor_host_jac(t, y, jval)
    REAL(dp), INTENT(IN)  :: t
    REAL(dp), INTENT(IN)  :: y(:)
    REAL(dp), INTENT(OUT) :: jval(:)

    CALL rates_at(t, y)
    CALL jac_from_rate(y)
    IF (SIZE(jval) < SIZE(Jac_CC%Val)) STOP 'cminor_host_jac: jval too small'
    jval(1:SIZE(Jac_CC%Val)) = Jac_CC%Val
  END SUBROUTINE cminor_host_jac

  !--- one Rate eval → f and J -------------------------------------------
  SUBROUTINE cminor_host_f_and_jac(t, y, f, jval)
    REAL(dp), INTENT(IN)  :: t
    REAL(dp), INTENT(IN)  :: y(:)
    REAL(dp), INTENT(OUT) :: f(:)
    REAL(dp), INTENT(OUT) :: jval(:)

    CALL rates_at(t, y)
    CALL rhs_from_rate(f)
    CALL jac_from_rate(y)
    IF (SIZE(jval) < SIZE(Jac_CC%Val)) STOP 'cminor_host_f_and_jac: jval too small'
    jval(1:SIZE(Jac_CC%Val)) = Jac_CC%Val
  END SUBROUTINE cminor_host_f_and_jac

END MODULE Cminor_host_Mod
