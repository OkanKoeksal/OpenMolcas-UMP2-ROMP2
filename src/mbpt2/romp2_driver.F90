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

subroutine ROMP2_Driver(ireturn)
use Definitions, only: wp,iwp,u6
use UMP2_Global
use ROMP2_Reference, only: ROMP2_Read,ROMP2_Fock
use ROMP2_Cholesky, only: ROMP2_Cho_Open,ROMP2_Cho_Close,ROMP2_Cho_Energy
use ROMP2_Semicanonical, only: ROMP2_Semi
use UMP2_Conventional, only: UMP2_Conventional_Energy
use stdalloc, only: mma_allocate,mma_deallocate
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
#include "warnings.h"
integer(kind=iwp),intent(out) :: ireturn
integer(kind=iwp) :: ierr,i,cleanupRC
real(kind=wp) :: E1A,E1B
real(kind=wp),allocatable :: FA(:,:),FB(:,:)
character(len=256) :: Message
ireturn=_RC_INPUT_ERROR_
call UMP2_Clean()
call ROMP2_Read(ierr,Message)
if (ierr/=0) goto 900
call UMP2_RdInp(ierr,Message)
if (ierr/=0) goto 900
if ((nFro<0).or.(nFro>min(nOccA,nOccB))) then
  Message='ROMP2 Frozen must be between zero and the number of doubly occupied reference orbitals.'
  goto 900
end if
if (DoCholesky) then
  call ROMP2_Cho_Open(ierr,Message)
  if (ierr/=0) goto 900
end if
call mma_allocate(FA,nBas,nBas,label='ROMP2 Fock alpha')
call mma_allocate(FB,nBas,nBas,label='ROMP2 Fock beta')
call ROMP2_Fock(FA,FB,ierr,Message)
if (ierr/=0) goto 900
call ROMP2_Semi(CAlpha,FA,nBas,nOccA,nFro,EOrbA,E1A,ierr)
if (ierr/=0) then
  Message='ROMP2 alpha semicanonicalization or singles evaluation failed.'
  goto 900
end if
call ROMP2_Semi(CBeta,FB,nBas,nOccB,nFro,EOrbB,E1B,ierr)
if (ierr/=0) then
  Message='ROMP2 beta semicanonicalization or singles evaluation failed.'
  goto 900
end if
call mma_deallocate(FA)
call mma_deallocate(FB)
call UMP2_Setup_Spaces(ierr,Message)
if (ierr/=0) goto 900
if (DoCholesky) then
  call ROMP2_Cho_Energy(ierr,Message)
else
  call UMP2_Conventional_Energy(ierr,Message)
end if
if (ierr/=0) goto 900
ECorr=E1A+E1B+EAA+EBB+EAB
ETotal=ESCF+ECorr
if (.not.ieee_is_finite(ETotal)) then
  Message='Nonfinite ROMP2 total energy.'
  goto 900
end if
! Finalize before publishing any energy as successful.
call ROMP2_Cho_Close(ierr)
if (ierr/=0) then
  Message='ROMP2 Cholesky finalization failed.'
  goto 900
end if
if (DoCholesky) then
  write(u6,'(/,A)') ' Semicanonical ROHF-MBPT(2): serial C1 Cholesky'
else
  write(u6,'(/,A)') ' Semicanonical ROHF-MBPT(2): serial C1 conventional'
end if
write(u6,'(A,I8)') ' Frozen occupied orbitals per spin: ',nFro
write(u6,'(A,2I8)') ' Total occupied orbitals (A,B):     ',nOccA,nOccB
write(u6,'(A,2I8)') ' Correlated occupied orbitals (A,B):',nOA,nOB
write(u6,'(A,2I8)') ' Virtual orbitals (A,B):            ',nVA,nVB
if (nFro>0) then
  write(u6,'(A)') ' Core convention: lowest occupied orbitals AFTER full spin semicanonicalization.'
  write(u6,'(A)') ' Excluded index       alpha energy / Eh        beta energy / Eh'
  do i=1,nFro
    write(u6,'(I15,2F25.12)') i,EOrbA(i),EOrbB(i)
  end do
else
  write(u6,'(A)') ' All-electron correlation (no frozen orbitals).'
end if
write(u6,'(A,F24.12)') ' ROHF reference energy: ',ESCF
write(u6,'(A,F24.12)') ' ROMP2 alpha singles:  ',E1A
write(u6,'(A,F24.12)') ' ROMP2 beta singles:   ',E1B
write(u6,'(A,F24.12)') ' ROMP2 AA doubles:     ',EAA
write(u6,'(A,F24.12)') ' ROMP2 BB doubles:     ',EBB
write(u6,'(A,F24.12)') ' ROMP2 AB doubles:     ',EAB
write(u6,'(A,F24.12)') ' ROMP2 correlation:    ',ECorr
write(u6,'(A,F24.12)') ' ROMP2 total energy:   ',ETotal
call Store_Energies(1,[ETotal],1)
call Put_cArray('Relax Method','ROMP2   ',8)
call Put_iScalar('mp2prpt',0)
call Add_Info('E_ROMP2',[ETotal],1,8)
call Add_Info('E_ROMP2_SINGLES',[E1A+E1B],1,8)
call Add_Info('E_ROMP2_AA',[EAA],1,8)
call Add_Info('E_ROMP2_BB',[EBB],1,8)
call Add_Info('E_ROMP2_AB',[EAB],1,8)
ireturn=_RC_ALL_IS_WELL_
900 continue
call ROMP2_Cho_Close(cleanupRC)
if (cleanupRC/=0) write(u6,'(A,I8)') ' ROMP2 Cholesky cleanup return code: ',cleanupRC
if (ireturn/=_RC_ALL_IS_WELL_) write(u6,'(A)') ' ROMP2 error: '//trim(Message)
call mma_deallocate(FA,safe='*')
call mma_deallocate(FB,safe='*')
call UMP2_Clean()
end subroutine ROMP2_Driver
