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
! Thin CLI driver (KPP-style). Orchestration lives in Cminor_Driver_Mod.
!
PROGRAM Cminor

  USE Control_Mod,       ONLY: F90_PATH_MAX
  USE Cminor_Driver_Mod, ONLY: cminor_initialize, cminor_run, cminor_reduce, cminor_finalize

  IMPLICIT NONE

  CHARACTER(F90_PATH_MAX) :: runfile = ''
  INTEGER :: ios

  IF ( COMMAND_ARGUMENT_COUNT() >= 1 ) THEN
    CALL GET_COMMAND_ARGUMENT(1, runfile)
  END IF
  IF ( LEN_TRIM(runfile) == 0 ) THEN
    WRITE(*,'(A)', ADVANCE='NO') 'Input RUNFilename: '
    READ(*,'(A)', IOSTAT=ios) runfile
    IF ( ios /= 0 .OR. LEN_TRIM(runfile) == 0 ) STOP 'No run file given'
  END IF

  CALL cminor_initialize(runfile)
  CALL cminor_run()
  CALL cminor_reduce()
  CALL cminor_finalize()

END PROGRAM Cminor
