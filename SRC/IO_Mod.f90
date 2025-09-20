!
! Copyright (C) 2025 Levin Rug, Willi Schimmel (E-Mail: l.rug@lmu.de)
! See ./SRC/Cminor.f90 for the copyright notice
! See ./LICENSE for license information
! SPDX-License-Identifier: GPL-3.0
!
MODULE IO_Mod
  IMPLICIT NONE
  !
  CONTAINS
  SUBROUTINE Logo()
    WRITE(*,*) 
    WRITE(*,777) ! in Georgi16
    WRITE(*,777) "***************************************************************"
    WRITE(*,777) "*     ____                                                    *"
    WRITE(*,777) "*    6MMMMb/                 68b                              *"
    WRITE(*,777) "*   8P    YM                 Y89                              *"
    WRITE(*,777) "*  6M      Y ___  __    __   ___ ___  __     _____   ___  __  *"
    WRITE(*,777) "*  MM        `MM 6MMb  6MMb  `MM `MM 6MMb   6MMMMMb  `MM 6MM  *"
    WRITE(*,777) "*  MM         MM69 `MM69 `Mb  MM  MMM9 `Mb 6M'   `Mb  MM69 "//'"'//"  *"
    WRITE(*,777) "*  MM         MM'   MM'   MM  MM  MM'   MM MM     MM  MM'     *"
    WRITE(*,777) "*  MM         MM    MM    MM  MM  MM    MM MM     MM  MM      *"
    WRITE(*,777) "*  YM      6  MM    MM    MM  MM  MM    MM MM     MM  MM      *"
    WRITE(*,777) "*   8b    d9  MM    MM    MM  MM  MM    MM YM.   ,M9  MM      *"
    WRITE(*,777) "*    YMMMM9  _MM_  _MM_  _MM__MM__MM_  _MM_ YMMMMM9  _MM_     *"
    WRITE(*,777) "*                                                             *"
    WRITE(*,777) "***************************************************************"
    WRITE(*,*) 
    WRITE(*,777) "               The Chemical Mechanism Integrator               " 
    WRITE(*,777) "         Copyright (C) 2025 Levin Rug, Willi Schimmel          "
    WRITE(*,*) 
    WRITE(*,777) "        This program comes with ABSOLUTELY NO WARRANTY,        "
    WRITE(*,777) "  it is free software, and you are welcome to redistribute it  "
    WRITE(*,777) "     under certain conditions, for details see ./LICENSE,      "
    WRITE(*,777) "       and the copyright statement in ./SRC/Cminor.f90.        "
    WRITE(*,*) 
    WRITE(*,*) 
    WRITE(*,*) ; WRITE(*,*)
    777 FORMAT(10X,A)
  END SUBROUTINE Logo

  !
  SUBROUTINE Print_Run_Param()
    USE Control_Mod, ONLY: Simulation, SysFile, NetcdfFile, InitFile, ODEsolver, RtolROW &
    &                    , AtolGas, AtolAqua, AtolTemp, combustion, Error_Est   &
    &                    , Ordering, ATolWaterMass, ATolq, ATolz, ATolRho, adiabatic_parcel
    !
    USE Reac_Mod,    ONLY: ns_AQUA

    IF ( INDEX(SysFile,'.sys')==0)  SysFile = TRIM(SysFile)//'.sys'

    IF (Simulation) THEN
      WRITE(*,*)
      WRITE(*,777)   'Run - Parameter:'
      WRITE(*,*)
      WRITE(*,777)   '    Mechanism:             '//TRIM(SysFile)
      IF (NetCdfFile /= '') THEN
        WRITE(*,777)   '    NetCDF-File:           '//TRIM(NetCdfFile)
      ELSE
        WRITE(*,777)   '    NetCDF-File:           *** no NetCDF output ***'
      END IF
      WRITE(*,777)   '    Initial value file:    '//TRIM(InitFile)
      WRITE(*,777)   '    ODE solver:            '//TRIM(ODEsolver)
      IF (Error_Est==2) THEN
        WRITE(*,777)   '    Error Estimation:      Euklid Norm'
      ELSE
        WRITE(*,777)   '    Error Estimation:      Maximum Norm'
      END IF
      IF (Ordering .EQV. .FALSE.) THEN
        WRITE(*,777)   '    Solve Linear Systems:  No sparsity-preserving ordering, decomposing the original matrices.'
      END IF
      WRITE(*,777)
      WRITE(*,777)   'Tolerance:   '
      WRITE(*,777)
      WRITE(*,'(10X,A,2X,Es8.2)')   '    Relative Rosenbrock        = ',RtolROW
      WRITE(*,'(10X,A,2X,Es8.2)')   '    Absolute (gaseous species) = ',AtolGas
      IF (ns_AQUA>0) WRITE(*,'(10X,A,2X,Es8.2)')   '    Absolute (aqueous species) = ',AtolAqua
      IF ( combustion ) THEN
        WRITE(*,'(10X,A,2X,Es8.2)')   '    Absolute Temperature       = ',AtolTemp
      END IF
      IF (adiabatic_parcel) THEN
        WRITE(*,'(10X,A,2X,Es8.2)')   '    Absolute Temperature       = ',AtolTemp
        WRITE(*,'(10X,A,2X,Es8.2)')   '    Absolute Water Mass        = ',AtolWaterMass
        WRITE(*,'(10X,A,2X,Es8.2)')   '    Absolute mixing ratio      = ',Atolq
        WRITE(*,'(10X,A,2X,Es8.2)')   '    Absolute height            = ',Atolz
        WRITE(*,'(10X,A,2X,Es8.2)')   '    Absolute air density       = ',AtolRho
      END IF
      WRITE(*,*) 
    END IF
    777 FORMAT(10X,A)
  END SUBROUTINE Print_Run_Param
  !
  !
  SUBROUTINE Output_Statistics
    USE Control_Mod, ONLY: Time_Read, TimeRates, TimeSymbolic, TimeIntegration, TimeFac, &
                         & TimeSolve, TimeJac, Time_Finish, TimeNetCDF, TimeErrCalc,    &
                         & TimeRhsCalc, TimeFluxWrite, TimeConcWrite, TimeSetValues,     &
                         & Out, ConcDataPrint, FluxDataPrint
    !
    CHARACTER(8) :: unit(14)
    !
    TimeIntegration = TimeIntegration - TimeNetCDF - TimeFluxWrite - TimeConcWrite

    CALL ConvertTime(Time_Read,unit(1)(:))
    CALL ConvertTime(TimeSymbolic,unit(2)(:))
    CALL ConvertTime(TimeNetCDF,unit(3)(:))
    CALL ConvertTime(TimeFluxWrite,unit(4)(:))
    CALL ConvertTime(TimeConcWrite,unit(5)(:))
    CALL ConvertTime(TimeFac,unit(6)(:))
    CALL ConvertTime(TimeRhsCalc,unit(7)(:))
    CALL ConvertTime(TimeSolve,unit(8)(:))
    CALL ConvertTime(TimeRates,unit(9)(:))
    CALL ConvertTime(TimeJac,unit(10)(:))
    CALL ConvertTime(TimeErrCalc,unit(11)(:))
    CALL ConvertTime(TimeIntegration,unit(12)(:))
    CALL ConvertTime(Time_Finish,unit(13)(:))
    CALL ConvertTime(TimeSetValues,unit(14)(:))
    !
    299 format(10X,A,3X,F9.5,A)
    298 format(10X,A,3X,I10)
    777 FORMAT(10X,A)
    WRITE(*,*);  WRITE(*,*);    WRITE(*,*)
    WRITE(*,777) 'Statistics (Numbers):'; 
    WRITE(*,*)
    WRITE(*,298) '    successful time steps   =', Out%nsteps
    WRITE(*,298) '    failed time steps       =', Out%nfailed
    WRITE(*,298) '    rate evaluations        =', Out%nRateEvals
    WRITE(*,298) '    Jacobian calculations   =', Out%npds
    WRITE(*,298) '    LU factorizations       =', Out%ndecomps
    WRITE(*,298) '    solved linear systems   =', Out%nsolves
    WRITE(*,*);  WRITE(*,*)
    WRITE(*,777)   'Statistics (Time):'
    WRITE(*,*)
    WRITE(*,299) '    reading mechanism       =', Time_Read,unit(1)
    WRITE(*,299) '    symbolic phase          =', TimeSymbolic,unit(2)
    WRITE(*,299) '    writing NetCDF-File     =', TimeNetCDF,unit(3)
    IF (FluxDataPrint) WRITE(*,299) '    writing flux-dataset    =', TimeFluxWrite,unit(4)
    IF (ConcDataPrint) WRITE(*,299) '    writing conc-dataset    =', TimeConcWrite,unit(5) ; WRITE(*,*)
    WRITE(*,299) '            factorization   =', TimeFac  ,unit(6)
    WRITE(*,299) '          + right-hand side =', TimeRhsCalc  ,unit(7)
    WRITE(*,299) '          + set LU values   =', TimeSetValues  ,unit(14)
    WRITE(*,299) '          + linear systems  =', TimeSolve,unit(8)
    WRITE(*,299) '          + reaction rates  =', TimeRates,unit(9)
    WRITE(*,299) '          + Jacobian        =', TimeJac  ,unit(10)
    WRITE(*,299) '          + error calc      =', TimeErrCalc  ,unit(11)
    WRITE(*,777) '    ------------------------=----------------------'
    WRITE(*,299) '    integration/step        =', TimeIntegration/(Out%nsteps+Out%nfailed),unit(12)
    WRITE(*,299) '    integration             =', TimeIntegration,unit(12); WRITE(*,*)
    WRITE(*,299) '    total runtime           =', Time_Finish,unit(13)
    WRITE(*,*);  WRITE(*,*);  WRITE(*,*)
  END SUBROUTINE
  !

  SUBROUTINE  Matrix_Statistics(A,B,BAT,Jac,M,LUM)
    
    USE Sparse_Mod, ONLY: CSR_Matrix_T
    TYPE(CSR_Matrix_T) :: A,B,BAT,Jac,M,LUM
    297 format(10X,A)
    298 format(10X,A18,3(I12,A2))

    WRITE(*,297) '                 |     rows    |    colums   |      nnz    |'
    WRITE(*,297) ' ----------------+-------------+-------------+-------------+-'
    WRITE(*,298) '           alpha |', A%m,   ' |',A%n,     ' |',A%nnz,   ' |'
    WRITE(*,298) '            beta |', B%m,   ' |',B%n,     ' |',B%nnz,   ' |'
    WRITE(*,298) '  (beta-alpha)^T |', BAT%m, ' |',BAT%n,   ' |',BAT%nnz, ' |'
    WRITE(*,298) '  Jacobian (= J) |', Jac%m, ' |',Jac%n,   ' |',Jac%nnz, ' |'
    WRITE(*,298) '       I - h*g*J |', M%m,   ' |',M%n,     ' |',M%nnz,   ' |'
    WRITE(*,298) '   LU(I - h*g*J) |', LUM%m, ' |',LUM%n,   ' |',LUM%nnz, ' |'
    WRITE(*,*)
    WRITE(*,*)

  END SUBROUTINE Matrix_Statistics

  SUBROUTINE ShowMaxErrorCounter()
    USE Control_Mod, ONLY: maxErrorCounter, LABEL, RtolROW, AtolGas, AtolAqua, &
    &                      AtolTemp, nDropletClasses, nD_Ptr_spc, combustion, adiabatic_parcel, combustion

    USE Reac_Mod,    ONLY: nspc, y_name, iAqMassEq2, iTeq2, iRhoEq2, iqEq2, iZeq2, nDIM2


    INTEGER :: i, ii

    OPEN(UNIT=99,FILE='OUTPUT/LocalErrorMaxima.log',STATUS='UNKNOWN')
    WRITE(99,*) '         Mechanism: ', TRIM(LABEL)
    WRITE(99,*) '  Tolerance   rel.: ', RtolROW
    WRITE(99,*) '      (gas)   abs.: ', AtolGas
    WRITE(99,*) '      (aqua)  abs.: ', AtolAqua
    WRITE(99,*) '      (temp)  abs.: ', AtolTemp
    WRITE(99,*)

    IF (nDropletClasses>1) THEN
      DO i=1,nspc
        DO ii = nD_Ptr_spc(i),nD_Ptr_spc(i+1)-1
          IF ( maxErrorCounter(ii) > 0 ) THEN
            WRITE(99,'(A,2X,I8,5X,A,2X,A,1X,I0)') ' number of local error maximum = ' , maxErrorCounter(ii), TRIM(y_name(i)), "DropletClass", ii-nD_Ptr_spc(i)+1
          END IF
        END DO
      END DO
    ELSE
      DO i=1,nspc
          IF ( maxErrorCounter(i) > 0 ) THEN
            WRITE(99,'(A,2X,I8,5X,A)') ' number of local error maximum = ' , maxErrorCounter(i), TRIM(y_name(i))
          END IF
      END DO
    END IF
    
    IF (combustion) WRITE(99,'(A,2X,I8,5X,A)') ' number of local error maximum = ' , maxErrorCounter(nDIM2), "Temperature"
    
    IF (adiabatic_parcel) THEN
      DO i=1,nDropletClasses
        IF ( maxErrorCounter(iAqMassEq2(i)) > 0 ) THEN
          WRITE(99,'(A,2X,I8,5X,A,1X,I0)') ' number of local error maximum = ' , maxErrorCounter(iAqMassEq2(i)), "Water mass Droplet Class", i
        END IF
      END DO
      IF ( maxErrorCounter(iTeq2) > 0 ) THEN
        WRITE(99,'(A,2X,I8,5X,A)') ' number of local error maximum = ' , maxErrorCounter(iTeq2), "Temperature of parcel"
      END IF
      IF ( maxErrorCounter(iRhoEq2) > 0 ) THEN
        WRITE(99,'(A,2X,I8,5X,A)') ' number of local error maximum = ' , maxErrorCounter(iRhoEq2), "Density of parcel"
      END IF
      IF ( maxErrorCounter(iqEq2) > 0 ) THEN
        WRITE(99,'(A,2X,I8,5X,A)') ' number of local error maximum = ' , maxErrorCounter(iqEq2), "Mixing ratio (q)"
      END IF
      IF ( maxErrorCounter(iZeq2) > 0 ) THEN
        WRITE(99,'(A,2X,I8,5X,A)') ' number of local error maximum = ' , maxErrorCounter(iZeq2), "Height (z)"
      END IF
    END IF

    CLOSE(99)
  END SUBROUTINE ShowMaxErrorCounter

  SUBROUTINE OpenFile_wStream(UnitNr,FileName)
    INTEGER,      INTENT(IN) :: UnitNr
    CHARACTER(*), INTENT(IN) :: FileName
    INTEGER :: io_stat
    OPEN(unit=UnitNr, file=FileName, status='replace', action='write', access='stream', iostat=io_stat)
    CALL file_err(FileName,io_stat)
  END SUBROUTINE OpenFile_wStream

  SUBROUTINE OpenFile_wSeq(UnitNr,FileName)
    INTEGER,      INTENT(IN) :: UnitNr
    CHARACTER(*), INTENT(IN) :: FileName
    INTEGER :: io_stat
    OPEN(unit=UnitNr, file=FileName, status='replace', action='write', access='sequential', iostat=io_stat)
    CALL file_err(FileName,io_stat)
  END SUBROUTINE OpenFile_wSeq

  SUBROUTINE OpenFile_rSeq(UnitNr,FileName)
    INTEGER,      INTENT(IN) :: UnitNr
    CHARACTER(*), INTENT(IN) :: FileName
    INTEGER :: io_stat
    OPEN(unit=UnitNr, file=FileName, status='old', action='read', access='sequential', iostat=io_stat)
    CALL file_err(FileName,io_stat)
  END SUBROUTINE OpenFile_rSeq

  SUBROUTINE OpenFile_rStream(UnitNr,FileName)
    INTEGER,      INTENT(IN) :: UnitNr
    CHARACTER(*), INTENT(IN) :: FileName
    INTEGER :: io_stat
    OPEN(unit=UnitNr, file=FileName, status='old', action='read', access='stream', iostat=io_stat)
    CALL file_err(FileName,io_stat)
  END SUBROUTINE OpenFile_rStream

  SUBROUTINE StreamWriteConcentrations(Conc)
    USE Kind_Mod,    ONLY: dp
    USE Control_Mod, ONLY: ConcUnit, ConcFile, ConcMetaUnit, ConcMetaFile, iStpConc
    REAL(dp), DIMENSION(:), INTENT(IN) :: Conc

    INTEGER :: io_stat, io_pos
    CHARACTER(100) :: io_msg

    OPEN(unit=ConcUnit,      file=ConcFile,  status='old',   action='write', &
    &    position='append', access='stream', iostat=io_stat, iomsg=io_msg    )
    CALL file_err(ConcFile,io_stat,io_msg)
    INQUIRE(ConcUnit, POS=io_pos)
    WRITE(ConcUnit) Conc
    CLOSE(ConcUnit)

    iStpConc   = iStpConc + 1
    OPEN(unit=ConcMetaUnit, file=ConcMetaFile, status='old', action='write', position='append')
    WRITE(ConcMetaUnit,*) iStpConc, io_pos
    CLOSE(ConcMetaUnit)

  END SUBROUTINE StreamWriteConcentrations

  SUBROUTINE file_err(filename,io_stat,io_msg)
    CHARACTER(Len=*), INTENT(in) :: filename
    INTEGER         , INTENT(in) :: io_stat
    CHARACTER(Len=*), INTENT(in), OPTIONAL :: io_msg
    IF (io_stat /= 0) THEN
      WRITE(*,"(79('!'))")
      WRITE(*,'(A,I0)')    'ERROR operating on file:  '//TRIM(filename)//'  with io status:  ',io_stat 
      IF (PRESENT(io_msg)) WRITE(*,'(A)')       'Message:  '//TRIM(io_msg)
      WRITE(*,"(79('!'))")
      WRITE(*,*)'Exit ...'
      STOP
    END IF
  END SUBROUTINE file_err

  ! writing reaction rates , time and stepsize h to file via stream access
  SUBROUTINE StreamWriteFluxes(Rate,t,h)
    USE Kind_Mod,    ONLY: dp
    USE Control_Mod, ONLY: FluxUnit, FluxFile, FluxMetaUnit, FluxMetaFile, iStpFlux
    REAL(dp) :: Rate(:)
    REAL(dp) :: t , h

    INTEGER :: io_stat, io_pos
    CHARACTER(100) :: io_msg

    OPEN(unit=FluxUnit,      file=FluxFile,  status='old',   action='write', &
    &    position='append', access='stream', iostat=io_stat, iomsg=io_msg    )
    CALL file_err(FluxFile,io_stat,io_msg)
    INQUIRE(FluxUnit, POS=io_pos)
    WRITE(FluxUnit) Rate
    CLOSE(FluxUnit)

    iStpFlux   = iStpFlux + 1
    OPEN(unit=FluxMetaUnit, file=FluxMetaFile, status='old', action='write', position='append')
    WRITE(FluxMetaUnit,*) iStpFlux, io_pos ,t , h
    CLOSE(FluxMetaUnit)

  END SUBROUTINE StreamWriteFluxes

  SUBROUTINE ConvertTime(t,fmt)
    USE Kind_Mod, ONLY: dp
    REAL(dp),   INTENT(INOUT) :: t    ! in seconds
    CHARACTER(8), INTENT(OUT) :: fmt

    IF ( t > 60.0_dp) THEN
      t   = t/60.0_dp
      fmt = ' [min]'
      IF ( t > 60.0_dp) THEN
        t   = t/60.0_dp
        fmt = ' [hours]'
        IF ( t > 24.0_dp) THEN
          t   = t/24.0_dp
          fmt = ' [days]'
        END IF
      END IF
    ELSE
      fmt = ' [sec]'
    END IF
    
  END SUBROUTINE ConvertTime
END MODULE IO_Mod

