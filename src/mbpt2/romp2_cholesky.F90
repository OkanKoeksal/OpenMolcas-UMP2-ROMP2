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

module ROMP2_Cholesky
! Serial/MPI C1 Cholesky lifecycle shared by the reference-Fock and doubles
! stages. AO vectors are streamed twice; no four-index AO tensor is built.
use Definitions, only: wp,iwp,u6
use Para_Info, only: nProcs,MyRank
use ROMP2_Parallel
use Data_Structures, only: SBA_Type,Allocate_DT,Deallocate_DT
use Cholesky, only: NumCho,nDimRS,ChoNSym=>nSym,ChoNBas=>nBas
use UMP2_Global, only: nBas,nOA,nOB,nVA,nVB,EOccA,EOccB,EVirA,EVirB,DenTol,EAA,EBB,EAB
use UMP2_Cholesky_Transform, only: UMP2_Transform_AO_Batch
use stdalloc, only: mma_allocate,mma_deallocate,mma_maxDBLE
use linalg_mod, only: mult
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
private
public :: ROMP2_Cho_Open,ROMP2_Cho_Close,ROMP2_Cho_Fock,ROMP2_Cho_Energy
logical,save :: ChoReady=.false.,HaveAO=.false.
integer(kind=iwp),save :: nVec=0,nGlobal=0,Offset=0,nRead=0,iRedC=-1
integer(kind=iwp),save :: Skip(8)=0
real(kind=wp),allocatable,save :: RedVec(:)
type(SBA_Type),target,save :: AO
contains

subroutine ROMP2_Cho_Open(ierr,Message)
  integer(kind=iwp),intent(out) :: ierr
  character(len=*),intent(out) :: Message
  integer(kind=iwp) :: Dims(1),Available,nWorkers
  integer(kind=iwp),allocatable :: Counts(:)
  real(kind=wp) :: Required
#include "warnings.h"
  ierr=1
  Message='ROMP2 Cholesky initialization failed.'
  if (ROMP2_Any(ChoReady.or.HaveAO)) then
    Message='ROMP2 Cholesky interface is already initialized.'
    return
  end if
  call Cho_X_Init(ierr,0.0_wp)
  if (ROMP2_Any(ierr/=0)) then
    ! A partial initialization must not finalize on only a subset of ranks.
    Message='ROMP2 Cho_X_Init failed on at least one rank.'
    call Quit(_RC_CHO_INI_)
    return
  end if
  ChoReady=.true.
  ierr=1
  if (ROMP2_Any((ChoNSym/=1).or.(ChoNBas(1)/=nBas))) then
    Message='Cholesky and ROHF reference dimensions disagree.'
    return
  end if
  nVec=NumCho(1) ! LOCAL count, possibly zero.
  if (ROMP2_Any(nVec<0)) then
    Message='Invalid local Cholesky count.'
    return
  end if
  nWorkers=1
  if (ROMP2_IsParallel()) nWorkers=nProcs
  call mma_allocate(Counts,nWorkers,label='ROMP2 vector counts')
  Counts=0
  if (ROMP2_IsParallel()) then
    Counts(MyRank+1)=nVec
    call gaIgOP(Counts(1),nWorkers,'+')
  else
    Counts(1)=nVec
  end if
  if (sum(real(Counts,wp))>real(huge(nGlobal)-1,wp)) then
    call mma_deallocate(Counts)
    Message='Global Cholesky count exceeds integer range.'
    return
  end if
  nGlobal=sum(Counts)
  Offset=0
  if (ROMP2_IsParallel().and.(MyRank>0)) Offset=sum(Counts(1:MyRank))
  call mma_deallocate(Counts)
  if (nGlobal<=0) then
    Message='The global Cholesky vector count is zero.'
    return
  end if
  if (ROMP2_Any(.not.allocated(nDimRS))) then
    Message='Missing Cholesky reduced-set dimensions.'
    return
  end if
  if (ROMP2_Any((size(nDimRS,1)<1).or.(size(nDimRS,2)<1))) then
    Message='Invalid Cholesky reduced-set dimensions.'
    return
  end if
  nRead=maxval(nDimRS(1,:))
  if (ROMP2_Any((nVec>0).and.(nRead<=0))) then
    Message='Empty reduced sets on a vector-owning rank.'
    return
  end if
  nRead=max(1,nRead)
  Required=real(nRead,wp)+real(nBas,wp)**2
  call mma_maxDBLE(Available)
  if (ROMP2_Any((Required>0.75_wp*real(Available,wp)).or. &
      (Required>real(huge(Available)-1,wp)/real(storage_size(0.0_wp),wp)))) then
    Message='Insufficient memory or integer range for ROMP2 AO Cholesky reader.'
    return
  end if
  call mma_allocate(RedVec,nRead,label='ROMP2 reduced vector')
  Dims(1)=nBas
  call Allocate_DT(AO,Dims,Dims,1,1,1,0,Label='ROMP2 AO Cholesky vector')
  HaveAO=.true.
  AO%ipOff=1
  Skip=0
  Skip(1)=1
  iRedC=-1
  if (ROMP2_IsRoot()) write(u6,'(A,I10)') ' ROMP2 Cholesky vectors (global): ',nGlobal
  ierr=0
  Message=''
end subroutine ROMP2_Cho_Open

subroutine ReadVector(iVec,ierr,Message)
  integer(kind=iwp),intent(in) :: iVec
  integer(kind=iwp),intent(out) :: ierr
  character(len=*),intent(out) :: Message
  real(kind=wp) :: Scale
  ierr=1
  Message='ROMP2 Cholesky reader is not initialized.'
  if ((.not.ChoReady).or.(.not.HaveAO)) return
  ! Same full symmetric AO contract and reduced-set handling as UMP2.
  call Cho_X_getVfull(ierr,RedVec,nRead,iVec,1,1,2,iRedC,AO,Skip,.true.)
  if (ierr/=0) then
    Message='ROMP2 Cholesky AO-vector read failed.'
    return
  end if
  ierr=1
  Message='Nonfinite ROMP2 AO Cholesky vector.'
  if (.not.all(ieee_is_finite(AO%SB(1)%A3))) return
  Scale=max(1.0_wp,maxval(abs(AO%SB(1)%A3(:,:,1))))
  Message='ROMP2 AO Cholesky vector is not a full symmetric matrix.'
  if (maxval(abs(AO%SB(1)%A3(:,:,1)-transpose(AO%SB(1)%A3(:,:,1))))>1.0e-12_wp*Scale) return
  ierr=0
  Message=''
end subroutine ReadVector

subroutine ROMP2_Cho_Fock(DA,DB,FA,FB,ierr,Message)
  ! FA/FB enter as H. Both densities include ALL occupied orbitals,
  ! including the core: freezing changes correlation, not the reference.
  ! J_P = L_P * sum_uv [(DA+DB)_uv L_P,uv]
  ! K_sigma,P = L_P D_sigma L_P^T; L_P is checked symmetric.
  real(kind=wp),intent(in) :: DA(nBas,nBas),DB(nBas,nBas)
  real(kind=wp),intent(inout) :: FA(nBas,nBas),FB(nBas,nBas)
  integer(kind=iwp),intent(out) :: ierr
  character(len=*),intent(out) :: Message
  real(kind=wp),allocatable :: Work(:,:),Exchange(:,:)
  real(kind=wp) :: Charge
  integer(kind=iwp) :: iVec,Available
  ierr=1
  Message='ROMP2 Cholesky interface is not initialized for Fock construction.'
  if (ROMP2_Any((.not.ChoReady).or.(.not.HaveAO))) return
  call mma_maxDBLE(Available)
  if (ROMP2_Any(2.0_wp*real(nBas,wp)**2>0.75_wp*real(Available,wp))) then
    Message='Insufficient memory for ROMP2 Cholesky Fock work arrays.'
    return
  end if
  call mma_allocate(Work,nBas,nBas,label='ROMP2 Cholesky density work')
  call mma_allocate(Exchange,nBas,nBas,label='ROMP2 Cholesky exchange')
  ! H is present only on root, so the reduction includes it exactly once.
  if (.not.ROMP2_IsRoot()) then
    FA=0.0_wp
    FB=0.0_wp
  end if
  ierr=0
  iRedC=-1
  do iVec=1,nVec
    call ReadVector(iVec,ierr,Message)
    if (ierr/=0) exit
    Charge=sum(DA*AO%SB(1)%A3(:,:,1))+sum(DB*AO%SB(1)%A3(:,:,1))
    FA=FA+Charge*AO%SB(1)%A3(:,:,1)
    FB=FB+Charge*AO%SB(1)%A3(:,:,1)
    call mult(AO%SB(1)%A3(:,:,1),DA,Work)
    call mult(Work,AO%SB(1)%A3(:,:,1),Exchange)
    FA=FA-Exchange
    call mult(AO%SB(1)%A3(:,:,1),DB,Work)
    call mult(Work,AO%SB(1)%A3(:,:,1),Exchange)
    FB=FB-Exchange
  end do
  call mma_deallocate(Work)
  call mma_deallocate(Exchange)
  ! Local vector loops have unequal lengths: synchronize only after the loop.
  call ROMP2_SyncError(ierr,Message)
  if (ierr/=0) return
  if (ROMP2_IsParallel()) then
    call GADGOp(FA(1,1),size(FA),'+')
    call GADGOp(FB(1,1),size(FB),'+')
  end if
  if (ROMP2_Any(logical( &
      (.not.all(ieee_is_finite(FA))).or.(.not.all(ieee_is_finite(FB))), &
      kind=kind(.true.)))) then
    ierr=1
    Message='Nonfinite ROMP2 Cholesky Fock matrix.'
    return
  end if
  Message=''
end subroutine ROMP2_Cho_Fock

subroutine ROMP2_Cho_Energy(ierr,Message)
  ! CAlpha/CBeta have ALREADY been independently semicanonicalized.
  ! The shared UMP2 transformation excludes nFro occupied orbitals.
  integer(kind=iwp),intent(out) :: ierr
  character(len=*),intent(out) :: Message
  integer(kind=iwp) :: Available,iVec,GlobalVec
  real(kind=wp) :: NAElements,NBElements,Required
  real(kind=wp),allocatable :: La(:,:,:),Lb(:,:,:)
  ierr=1
  Message='ROMP2 Cholesky interface is not initialized for doubles.'
  if (ROMP2_Any((.not.ChoReady).or.(.not.HaveAO))) return
  NAElements=real(nVA,wp)*real(nOA,wp)*real(nGlobal,wp)
  NBElements=real(nVB,wp)*real(nOB,wp)*real(nGlobal,wp)
  if (ROMP2_Any(max(NAElements,NBElements)> &
      real(huge(Available)-1,wp)/real(storage_size(0.0_wp),wp))) then
    Message='ROMP2 transformed vectors exceed supported integer indexing range.'
    return
  end if
  Required=NAElements+NBElements+3.0_wp*real(nBas,wp)*real(max(1,nOA,nOB),wp)
  call mma_maxDBLE(Available)
  if (ROMP2_Any(Required>0.75_wp*real(Available,wp))) then
    Message='Insufficient memory for in-core ROMP2 occupied-virtual Cholesky vectors.'
    return
  end if
  call mma_allocate(La,nVA,nOA,nGlobal,label='ROMP2 alpha OV Cholesky')
  call mma_allocate(Lb,nVB,nOB,nGlobal,label='ROMP2 beta OV Cholesky')
  La=0.0_wp
  Lb=0.0_wp
  ierr=0
  iRedC=-1
  do iVec=1,nVec
    call ReadVector(iVec,ierr,Message)
    if (ierr/=0) exit
    GlobalVec=Offset+iVec
    call UMP2_Transform_AO_Batch(AO%SB(1)%A3, &
                               La(:,:,GlobalVec:GlobalVec),Lb(:,:,GlobalVec:GlobalVec),ierr,Message)
    if (ierr/=0) exit
  end do
  call ROMP2_SyncError(ierr,Message)
  if (ierr/=0) goto 900
  if (ROMP2_IsParallel()) then
    if (size(La)>0) call GADGOp(La(1,1,1),size(La),'+')
    if (size(Lb)>0) call GADGOp(Lb(1,1,1),size(Lb),'+')
  end if
  if (ROMP2_Any(logical( &
      (.not.all(ieee_is_finite(La))).or.(.not.all(ieee_is_finite(Lb))), &
      kind=kind(.true.)))) then
    ierr=1
    Message='Nonfinite assembled Cholesky vectors.'
    goto 900
  end if
  call UMP2_Energy(La,Lb,EOccA,EVirA,EOccB,EVirB,nOA,nVA,nOB,nVB,nGlobal,DenTol,EAA,EBB,EAB,ierr)
  if (ierr/=0) then
    Message='ROMP2 Cholesky doubles failed; check orbital-energy denominators.'
  else if (.not.all(ieee_is_finite([EAA,EBB,EAB]))) then
    ierr=1
    Message='Nonfinite ROMP2 Cholesky doubles energy.'
  else
    Message=''
  end if
900 continue
  call mma_deallocate(La)
  call mma_deallocate(Lb)
end subroutine ROMP2_Cho_Energy

subroutine ROMP2_Cho_Close(ierr)
  integer(kind=iwp),intent(out) :: ierr
  ierr=0
  if (HaveAO) then
    call Deallocate_DT(AO)
    HaveAO=.false.
  end if
  call mma_deallocate(RedVec,safe='*')
  if (ChoReady) then
    call Cho_X_Final(ierr)
    ChoReady=.false.
  end if
  nVec=0
  nGlobal=0
  Offset=0
  nRead=0
  iRedC=-1
end subroutine ROMP2_Cho_Close
end module ROMP2_Cholesky
