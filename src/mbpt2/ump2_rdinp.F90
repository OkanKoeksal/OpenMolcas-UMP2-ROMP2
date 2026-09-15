!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine UMP2_RdInp(ierr,Message)

use Definitions, only: iwp
use spool, only: SpoolInp, Close_LuSpool
use UMP2_Global, only: nSym, nOccA, nOccB, nFro, &
                       ExplicitFrozen, iPL, EnergyReady

implicit none

integer(kind=iwp), intent(out) :: ierr
character(len=*), intent(out) :: Message

integer(kind=iwp) :: LuSpool, EqPos, FrozenInput
logical(kind=iwp) :: HaveFrozen, Parsed

character(len=4) :: Command
character(len=180) :: Line, RawLine, ValueLine

integer(kind=iwp), external :: iPrintLevel
logical(kind=iwp), external :: Reduce_Prt
character(len=180), external :: Get_Ln

ierr = 0
Message = ''

! Invalidate old options and results before reading new input.
! Publish the frozen-core selection only after parsing succeeds.

nFro = 0
ExplicitFrozen = .false.
EnergyReady = .false.

HaveFrozen = .false.
FrozenInput = 0

! UMP2_Read_Reference must have completed before this routine.

if (nSym /= 1) then
  ierr = 1
  Message = 'UMP2 input processing requires a C1 reference.'
  return
end if

if ((nOccA < 0) .or. (nOccB < 0)) then
  ierr = 1
  Message = 'Invalid UHF occupation counts before input processing.'
  return
end if

! Follow the existing MBPT2 print-level convention.

iPL = iPrintLevel(-1)
if (Reduce_Prt() .and. (iPL < 3)) iPL = 0

! Use the same input-spooling convention as RdInp.
! Only one of the RHF and UHF readers should run per invocation.

LuSpool = 17
call SpoolInp(LuSpool)
rewind(LuSpool)
call RdNLst(LuSpool,'MBPT2')

ReadCommands: do

  RawLine = Get_Ln(LuSpool)
  Line = RawLine
  call StdFmt(Line,Command)

  select case (Command)

    case ('FROZ')

      if (HaveFrozen) then
        ierr = 1
        Message = 'Frozen was specified more than once.'
        exit ReadCommands
      end if

      ! Normally OpenMolcas supplies the value on the next
      ! processed line, as in the existing RdInp routine.
      ! Also accept an equals sign if it survives preprocessing.

      EqPos = index(RawLine,'=')

      if (EqPos > 0) then
        ValueLine = RawLine(EqPos+1:)
        if (len_trim(ValueLine) == 0) then
          ValueLine = Get_Ln(LuSpool)
        end if
      else
        ValueLine = Get_Ln(LuSpool)
      end if

      call ParseFrozen(ValueLine,FrozenInput,Parsed)

      if (.not. Parsed) then
        ierr = 1
        Message = 'Frozen requires exactly one nonnegative integer.'
        exit ReadCommands
      end if

      if (FrozenInput > min(nOccA,nOccB)) then
        ierr = 1
        Message = 'Frozen exceeds the occupied count of one spin set.'
        exit ReadCommands
      end if

      HaveFrozen = .true.

    case ('END ')

      exit ReadCommands

    case default

      ierr = 1
      Message = 'Unsupported UMP2 keyword: '//trim(Command)// &
                '. This initial version accepts Frozen only.'
      exit ReadCommands

  end select

end do ReadCommands

! Close the spool on both successful parsing and detected errors.

call Close_LuSpool(LuSpool)

if (ierr /= 0) return

if (.not. HaveFrozen) then
  ierr = 1
  Message = 'Specify Frozen explicitly; use Frozen = 0 for all electrons.'
  return
end if

nFro = FrozenInput
ExplicitFrozen = .true.

contains

subroutine ParseFrozen(Text,Value,Success)

  character(len=*), intent(in) :: Text
  integer(kind=iwp), intent(out) :: Value
  logical(kind=iwp), intent(out) :: Success

  character(len=len(Text)) :: Token
  integer(kind=iwp) :: FirstDigit, LastDigit, k, IOStatus

  Value = 0
  Success = .false.

  Token = trim(adjustl(Text))
  LastDigit = len_trim(Token)

  if (LastDigit == 0) return

  FirstDigit = 1
  if (Token(1:1) == '+') FirstDigit = 2

  if (FirstDigit > LastDigit) return

  ! Reject multiple values, negative numbers, decimal values,
  ! repetition syntax, and null list-directed values.
  !
  ! An ordinary list-directed integer read alone could silently
  ! accept extra fields or a slash without assigning the value.

  do k=FirstDigit,LastDigit
    if (index('0123456789',Token(k:k)) == 0) return
  end do

  read(Token,*,iostat=IOStatus) Value

  if (IOStatus /= 0) then
    Value = 0
    return
  end if

  if (Value < 0) then
    Value = 0
    return
  end if

  Success = .true.

end subroutine ParseFrozen

end subroutine UMP2_RdInp