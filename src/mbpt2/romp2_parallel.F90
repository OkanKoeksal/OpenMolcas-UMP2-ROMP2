!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2026, Okan Koeksal                                     *
!***********************************************************************

module ROMP2_Parallel
use Definitions, only: iwp,u6
use Para_Info, only: nProcs,MyRank,Is_Real_Par
implicit none
private
public :: ROMP2_IsParallel,ROMP2_IsRoot,ROMP2_Any,ROMP2_SyncError
contains
logical function ROMP2_IsParallel()
  ROMP2_IsParallel=(nProcs>1).and.Is_Real_Par()
end function
logical function ROMP2_IsRoot()
  ROMP2_IsRoot=(.not.ROMP2_IsParallel()).or.(MyRank==0)
end function
logical function ROMP2_Any(Local)
  ! All ranks must call this function at the same control-flow point.
  logical,intent(in) :: Local
  integer(kind=iwp) :: Flag
  Flag=0
  if (Local) Flag=1
  if (ROMP2_IsParallel()) call gaIgOP_SCAL(Flag,'max')
  ROMP2_Any=Flag/=0
end function
subroutine ROMP2_SyncError(Code,Message)
  integer(kind=iwp),intent(inout) :: Code
  character(len=*),intent(inout) :: Message
  if (ROMP2_Any(Code/=0)) then
    if (Code/=0) write(u6,'(A,I6,2A)') ' ROMP2 rank ',MyRank,': ',trim(Message)
    Code=1
    Message='ROMP2 failed on at least one rank; see preceding diagnostics.'
  end if
end subroutine
end module ROMP2_Parallel
