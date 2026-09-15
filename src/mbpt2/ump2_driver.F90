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


subroutine UMP2_Driver(ireturn)

use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use Definitions, only: wp, iwp, u6
use Para_Info, only: nProcs
use UnixInfo, only: SuperName
use Data_Structures, only: SBA_Type, Allocate_DT, Deallocate_DT
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Cholesky, only: NumCho, nDimRS, ChoNSym => nSym, ChoNBas => nBas
use UMP2_Global, only: UMP2_Clean, nSym, nBas, nOccA, nOccB, nFro, &
                       nOA, nOB, nVA, nVB, EOccA, EOccB, EVirA, EVirB, &
                       DoCholesky, DenTol, ESCF, EAA, EBB, EAB, ECorr, &
                       ETotal, EnergyReady
use UMP2_Cholesky_Transform, only: UMP2_Transform_AO_Batch

implicit none
#include "warnings.h"

integer(kind=iwp), intent(out) :: ireturn
integer(kind=iwp) :: ierr, FinalRC, nCho, nRead, iVec, iRedC
integer(kind=iwp) :: Available, Skip(8), Dims(1), iTol
integer(kind=iwp), external :: IsStructure, Cho_X_GetTol
logical(kind=iwp) :: IsDF, ChoReady, HaveAO
real(kind=wp) :: NAElements, NBElements, Required
real(kind=wp), allocatable :: RedVec(:), La(:,:,:), Lb(:,:,:)
type(SBA_Type), target :: AO
character(len=256) :: Message

! Development implementation: serial C1 canonical UHF, Cholesky,
! explicit frozen core, energy only. All transformed vectors remain
! in memory; AO vectors are expanded one at a time.
!
! This routine must be entered before the original RdMBPT/RdInp.
! It exclusively owns Cholesky initialization for this MBPT2 invocation.

ireturn = _RC_GENERAL_ERROR_
ChoReady = .false.
HaveAO = .false.
Message = ''
call UMP2_Clean()

if (nProcs /= 1) then
  call Fail(_RC_NOT_AVAILABLE_,'Initial UMP2 implementation requires MOLCAS_NPROCS=1.')
  return
end if

if ((IsStructure() == 1) .or. &
    (index(SuperName,'numerical_gradient') == 1)) then
  call Fail(_RC_NOT_AVAILABLE_,'Initial UMP2 implementation supports standalone energies only.')
  return
end if

call UMP2_Read_Reference(ierr,Message)
if (ierr /= 0) then
  call Fail(_RC_INPUT_ERROR_,Message)
  return
end if

call UMP2_RdInp(ierr,Message)
if (ierr /= 0) then
  call Fail(_RC_INPUT_ERROR_,Message)
  return
end if

call UMP2_Setup_Spaces(ierr,Message)
if (ierr /= 0) then
  call Fail(_RC_INPUT_ERROR_,Message)
  return
end if

call DecideOnCholesky(DoCholesky)
call DecideOnDF(IsDF)
if ((.not. DoCholesky) .or. IsDF) then
  call Fail(_RC_NOT_AVAILABLE_,'Initial UMP2 requires SEWARD Cholesky; conventional and RI runs are unsupported.')
  return
end if

write(u6,'(/,A)') ' Canonical UMP2: serial C1 Cholesky energy'
write(u6,'(A,I8)') ' Alpha occupied orbitals, including core: ',nOccA
write(u6,'(A,I8)') ' Beta occupied orbitals, including core:  ',nOccB
write(u6,'(A,I8)') ' Frozen occupied orbitals per spin:       ',nFro
write(u6,'(A,2I8)') ' Correlated occupied orbitals (A,B):      ',nOA,nOB
write(u6,'(A,2I8)') ' Virtual orbitals (A,B):                  ',nVA,nVB

! No vector cache: reserve memory for transformed spin sets.
call Cho_X_Init(ierr,0.0_wp)
if (ierr /= 0) then
  ! Cho_X_Init can leave partially initialized state on failure.
  ! Follow the existing Cholesky drivers' fatal-initialization policy;
  ! do not call Cho_X_Final on an incomplete initialization.
  write(u6,'(A,I8)') ' UMP2: Cho_X_Init failed, return code ',ierr
  call UMP2_Clean()
  call Quit(_RC_CHO_INI_)
  return
end if
ChoReady = .true.

if ((ChoNSym /= 1) .or. (ChoNBas(1) /= nBas)) then
  call Fail(_RC_CHO_INI_,'Cholesky and UHF reference dimensions disagree.')
  return
end if
nCho = NumCho(1)
if (nCho <= 0) then
  call Fail(_RC_CHO_INI_,'No Cholesky vectors are available.')
  return
end if
if (.not. allocated(nDimRS)) then
  call Fail(_RC_CHO_INI_,'Cholesky reduced-set dimensions are missing.')
  return
end if
if ((size(nDimRS,1) < 1) .or. (size(nDimRS,2) < 1)) then
  call Fail(_RC_CHO_INI_,'Invalid Cholesky reduced-set dimensions.')
  return
end if
nRead = maxval(nDimRS(1,:))
if (nRead <= 0) then
  call Fail(_RC_CHO_INI_,'Empty Cholesky reduced sets.')
  return
end if

! Check products in real arithmetic before allocating integer-indexed
! arrays and the bit-count arithmetic in stdalloc. This estimate includes
! headroom for matrix transformations and array-section temporaries.
NAElements = real(nVA,wp)*real(nOA,wp)*real(nCho,wp)
NBElements = real(nVB,wp)*real(nOB,wp)*real(nCho,wp)
if (max(NAElements,NBElements,real(nBas,wp)**2) > real(huge(nRead)-1,wp)/real(storage_size(0.0_wp),wp)) then
  call Fail(_RC_MEMORY_ERROR_,'UMP2 array dimensions exceed the supported integer indexing range.')
  return
end if
Required = NAElements+NBElements+real(nRead,wp) &
           +6.0_wp*real(nBas,wp)**2 &
           +3.0_wp*real(nBas,wp)*real(max(1,nOA,nOB),wp)
call mma_maxDBLE(Available)
if (Required > 0.75_wp*real(Available,wp)) then
  write(u6,'(A,ES16.6)') ' Estimated UMP2 working storage (real words): ',Required
  write(u6,'(A,I16)') ' Available real words:                      ',Available
  call Fail(_RC_MEMORY_ERROR_,'Insufficient memory for the in-memory UMP2 implementation.')
  return
end if

call mma_allocate(La,nVA,nOA,nCho,label='UMP2 La')
call mma_allocate(Lb,nVB,nOB,nCho,label='UMP2 Lb')
call mma_allocate(RedVec,nRead,label='UMP2 reduced vector')

! Case 0 gives A2(nBas*nBas,1) and A3(nBas,nBas,1) aliases.
! Allocate_SBA does NOT initialize ipOff: set it explicitly.
Dims(1) = nBas
call Allocate_DT(AO,Dims,Dims,1,1,1,0,Label='UMP2 full AO vector')
HaveAO = .true.
AO%ipOff = 1
Skip = 0
Skip(1) = 1
iRedC = -1

write(u6,'(A,I8)') ' Cholesky vectors:                       ',nCho

do iVec=1,nCho
  ! iSwap=2 produces full square matrices in the actual Cho_Reordr
  ! implementation. One vector per call avoids layout ambiguities
  ! between the differently documented batch orderings.
  call Cho_X_getVfull(ierr,RedVec,nRead,iVec,1,1,2,iRedC,AO,Skip,.true.)
  if (ierr /= 0) then
    write(u6,'(A,2I10)') ' UMP2 vector/read error: ',iVec,ierr
    call Fail(_RC_IO_ERROR_READ_,'Could not read and expand an AO Cholesky vector.')
    return
  end if

  call UMP2_Transform_AO_Batch(AO%SB(1)%A3, &
                              La(:,:,iVec:iVec),Lb(:,:,iVec:iVec),ierr,Message)
  if (ierr /= 0) then
    call Fail(_RC_GENERAL_ERROR_,Message)
    return
  end if
end do

! All Cholesky indices are now present. Never sum energies computed
! independently from individual Cholesky-vector batches.
call UMP2_Energy(La,Lb,EOccA,EVirA,EOccB,EVirB, &
                 nOA,nVA,nOB,nVB,nCho,DenTol,EAA,EBB,EAB,ierr)
if (ierr /= 0) then
  write(u6,'(A,I8)') ' UMP2 energy-kernel return code: ',ierr
  call Fail(_RC_GENERAL_ERROR_,'UMP2 energy evaluation failed; check orbital-energy denominators.')
  return
end if
ECorr = EAA+EBB+EAB
ETotal = ESCF+ECorr
if (.not. all(ieee_is_finite([ESCF,EAA,EBB,EAB,ECorr,ETotal]))) then
  call Fail(_RC_GENERAL_ERROR_,'UMP2 produced a nonfinite energy.')
  return
end if

! Close files and release vector storage before publishing a result.
call ReleaseVectors(FinalRC)
if (FinalRC /= 0) then
  call Fail(_RC_CHO_RUN_,'Cholesky finalization failed.')
  return
end if
EnergyReady = .true.

write(u6,'(/,A,F24.12)') ' UHF reference energy:     ',ESCF
write(u6,'(A,F24.12)')   ' UMP2 AA correlation:     ',EAA
write(u6,'(A,F24.12)')   ' UMP2 BB correlation:     ',EBB
write(u6,'(A,F24.12)')   ' UMP2 AB correlation:     ',EAB
write(u6,'(A,F24.12)')   ' UMP2 correlation energy: ',ECorr
write(u6,'(A,F24.12)')   ' UMP2 total energy:       ',ETotal
call PrintResult(u6,'(6X,A,T50,F19.10)','Total MBPT2 energy',0,'',[ETotal],1)
call Store_Energies(1,[ETotal],1)
! Distinguish this energy-only route from the RHF property machinery.
call Put_cArray('Relax Method','UMP2    ',8)
call Put_iScalar('mp2prpt',0)
iTol = Cho_X_GetTol(8)
call Add_Info('E_MP2',[ETotal],1,iTol)
call Add_Info('E_UMP2_AA',[EAA],1,iTol)
call Add_Info('E_UMP2_BB',[EBB],1,iTol)
call Add_Info('E_UMP2_AB',[EAB],1,iTol)
call xFlush(u6)
call UMP2_Clean()
ireturn = _RC_ALL_IS_WELL_

contains

subroutine ReleaseVectors(Code)
  integer(kind=iwp), intent(out) :: Code
  Code = 0
  if (HaveAO) then
    call Deallocate_DT(AO)
    HaveAO = .false.
  end if
  call mma_deallocate(RedVec,safe='*')
  call mma_deallocate(La,safe='*')
  call mma_deallocate(Lb,safe='*')
  if (ChoReady) then
    call Cho_X_Final(Code)
    ChoReady = .false.
  end if
end subroutine ReleaseVectors

subroutine Fail(Code,Text)
  integer(kind=iwp), intent(in) :: Code
  character(len=*), intent(in) :: Text
  integer(kind=iwp) :: CleanupRC
  write(u6,'(/,A)') ' UMP2 error: '//trim(Text)
  call ReleaseVectors(CleanupRC)
  if (CleanupRC /= 0) write(u6,'(A,I8)') ' Cholesky cleanup return code: ',CleanupRC
  call UMP2_Clean()
  ireturn = Code
end subroutine Fail

end subroutine UMP2_Driver
