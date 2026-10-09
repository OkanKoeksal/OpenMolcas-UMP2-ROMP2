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

subroutine UMP2_RdInp(ierr,Message,SymFrozen,SymOcc)

use Definitions, only: iwp
use spool, only: SpoolInp, Close_LuSpool
use UMP2_Global, only: nSym, nOccA, nOccB, nFro, &
                       ExplicitFrozen, iPL, EnergyReady

implicit none

integer(kind=iwp), intent(out) :: ierr
! Optional irrep-resolved input for the symmetry-blocked open-shell path.
integer(kind=iwp), intent(out), optional :: SymFrozen(:)
integer(kind=iwp), intent(in), optional :: SymOcc(:,:)
character(len=*), intent(out) :: Message

integer(kind=iwp) :: LuSpool, EqPos, FrozenInput, nInput, s
integer(kind=iwp) :: FrozenCounts(8)
logical(kind=iwp) :: HaveFrozen, Parsed, SymmetryInput

character(len=4) :: Command
character(len=180) :: Line, RawLine, ValueLine

integer(kind=iwp), external :: iPrintLevel
logical(kind=iwp), external :: Reduce_Prt
character(len=180), external :: Get_Ln

ierr = 0
Message = ''

! Invalidate old options and results before reading new input.
! Publish the frozen-core selection only after parsing succeeds.

SymmetryInput=present(SymFrozen).and.present(SymOcc)
if (present(SymFrozen).neqv.present(SymOcc)) then
  ierr=1
  Message='Incomplete symmetry frozen-orbital input arguments.'
  return
end if
if (SymmetryInput) then
  nInput=size(SymFrozen)
  if ((nInput<1).or.(nInput>8).or.(size(SymOcc,1)/=nInput).or.(size(SymOcc,2)/=2)) then
    ierr=1
    Message='Invalid symmetry frozen-orbital input dimensions.'
    return
  end if
  SymFrozen=0
else
  nInput=1
  nFro = 0
  ExplicitFrozen = .false.
  EnergyReady = .false.
end if

HaveFrozen = .false.
FrozenInput = 0
FrozenCounts=0

! UMP2_Read_Reference must have completed before this routine.

if ((.not.SymmetryInput).and.(nSym /= 1)) then
  ierr = 1
  Message = 'UMP2 input processing requires a C1 reference.'
  return
end if

if ((.not.SymmetryInput).and.((nOccA < 0) .or. (nOccB < 0))) then
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

      call ParseFrozen(ValueLine,FrozenCounts,nInput,Parsed)

      if (.not. Parsed) then
        ierr = 1
        Message = 'Frozen requires one nonnegative integer per irrep.'
        exit ReadCommands
      end if

      if (SymmetryInput) then
        do s=1,nInput
          if (FrozenCounts(s)>minval(SymOcc(s,:))) then
            ierr=1
            Message='Frozen exceeds the occupied count in an irrep.'
            exit
          end if
        end do
        if (ierr/=0) exit ReadCommands
      else
        FrozenInput=FrozenCounts(1)
        if (FrozenInput > min(nOccA,nOccB)) then
          ierr = 1
          Message = 'Frozen exceeds the occupied count of one spin set.'
          exit ReadCommands
        end if
      end if

      HaveFrozen = .true.

    case ('END ')

      exit ReadCommands

    case default

      ierr = 1
      Message = 'Unsupported UMP2 keyword: '//trim(Command)// &
                '. Supported keyword: Frozen.'
      exit ReadCommands

  end select

end do ReadCommands

! Close the spool on both successful parsing and detected errors.

call Close_LuSpool(LuSpool)

if (ierr /= 0) return

if (.not. HaveFrozen) then
  ierr = 1
  Message = 'Specify Frozen explicitly, with one count per irrep (zero to freeze none).'
  return
end if

if (SymmetryInput) then
  SymFrozen=FrozenCounts(1:nInput)
else
  nFro = FrozenInput
  ExplicitFrozen = .true.
end if

contains

subroutine ParseFrozen(Text,Values,nExpected,Success)

  character(len=*), intent(in) :: Text
  integer(kind=iwp), intent(out) :: Values(:)
  integer(kind=iwp), intent(in) :: nExpected
  logical(kind=iwp), intent(out) :: Success

  character(len=len(Text)) :: Token
  integer(kind=iwp) :: k, first, last, count, ios

  Values=0
  Success=.false.
  last=len_trim(Text)
  k=1
  count=0
  do while (k<=last)
    if (index(' ,'//achar(9),Text(k:k))/=0) then
      k=k+1
      cycle
    end if
    first=k
    if (Text(k:k)=='+') k=k+1
    if (k>last) return
    if (index('0123456789',Text(k:k))==0) return
    do while (k<=last)
      if (index('0123456789',Text(k:k))==0) exit
      k=k+1
    end do
    if (k<=last) then
      if (index(' ,'//achar(9),Text(k:k))==0) return
    end if
    count=count+1
    if (count>nExpected) return
    Token=Text(first:k-1)
    read(Token,*,iostat=ios) Values(count)
    if (ios/=0) return
    if (Values(count)<0) return
  end do
  Success=(count==nExpected)

end subroutine ParseFrozen

end subroutine UMP2_RdInp