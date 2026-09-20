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

module Open_Shell_Symmetry
! Abelian symmetry energy path. C1 continues to use the existing drivers.
! AO pair metrics and occupied-virtual pairs are grouped by product irrep.
! Conventional metrics and transformed vectors are replicated; computation
! is distributed. This is an in-core implementation, with memory checks.
use Definitions, only: wp,iwp,u6
use Symmetry_Info, only: Mul
use Para_Info, only: nProcs,MyRank
use ROMP2_Parallel, only: ROMP2_Any,ROMP2_IsRoot,ROMP2_IsParallel,ROMP2_SyncError
use ROMP2_Semicanonical, only: ROMP2_Semi
use stdalloc, only: mma_allocate,mma_deallocate,mma_maxDBLE
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
private
public :: Open_Shell_Symmetry_Driver
integer(kind=iwp) :: ns,nb,bs(8),off(8),occ(8,2),frozen(8),workers,rankid
integer(kind=iwp),allocatable :: symorb(:),pairindex(:,:),qindex(:,:,:)
real(kind=wp),allocatable :: coeff(:,:,:),eps(:,:),hcore(:,:),overlap(:,:),density(:,:,:),fock(:,:,:)
real(kind=wp) :: eref,singles(2)
logical :: restricted
logical(kind=iwp) :: cho
character(len=256) :: message
! Pair ordering matches ORDINT: high-irrep AO outer, low-irrep AO inner.
type Pair_Block
  integer(kind=iwp) :: np=0,nq=0,nlocal=0,nglobal=0,voff=0
  integer(kind=iwp),allocatable :: mu(:),nu(:),qi(:),qa(:),spin(:)
  real(kind=wp),allocatable :: metric(:,:),vectors(:,:),projected(:,:),kernel(:,:)
end type Pair_Block
type(Pair_Block) :: blocks(8)
contains

subroutine Open_Shell_Symmetry_Driver(method,ireturn)
  character(len=*),intent(in) :: method
  integer(kind=iwp),intent(out) :: ireturn
  integer(kind=iwp) :: ierr,g,spin,s,n,first,last,rc
  real(kind=wp) :: e(3),total,corr
  real(kind=wp),allocatable :: c(:,:),f(:,:),en(:)
  logical :: cho_open
#include "warnings.h"
  ireturn=_RC_INPUT_ERROR_
  message='Symmetry reference initialization failed.'
  cho_open=.false.
  restricted=(method=='CASSCF  ').or.(method=='RASSCF  ')
  workers=1
  rankid=0
  if (ROMP2_IsParallel()) then
    workers=nProcs
    rankid=MyRank
  end if
  call ReadReference(ierr)
  call ROMP2_SyncError(ierr,message)
  if (ierr/=0) goto 900
  ! RUNFILE records are logically global, but make rank zero authoritative
  ! for all real reference data used by distributed contractions.  This
  ! avoids rank-local RUNFILE state from changing the MPI result.
  call SynchronizeReference()
  call ReadFrozen(ierr)
  call ROMP2_SyncError(ierr,message)
  if (ierr/=0) goto 900
  ! Likewise use the root input parse as the single frozen-space definition.
  call SynchronizeFrozen()
  call MakePairs(ierr)
  if (ierr/=0) goto 900
  if (cho) then
    call ReadCholesky(ierr,cho_open)
  else
    call ReadConventional(ierr)
  end if
  if (ierr/=0) goto 900
  singles=0.0_wp
  if (restricted) then
    call ReferenceFock(ierr)
    if (ierr/=0) goto 900
    ! Each irrep and spin is diagonalized independently. Occupied and
    ! virtual subspaces never mix, including in degenerate spectra.
    ierr=0
    if (ROMP2_IsRoot()) then
      do spin=1,2
        do s=1,ns
          n=bs(s)
          if (n==0) cycle
          first=off(s)+1
          last=off(s)+n
          call mma_allocate(c,n,n,label='symmetry semi coefficients')
          call mma_allocate(f,n,n,label='symmetry semi Fock')
          call mma_allocate(en,n,label='symmetry semi energies')
          c=coeff(first:last,first:last,spin)
          f=fock(first:last,first:last,spin)
          call ROMP2_Semi(c,f,n,occ(s,spin),frozen(s),en,corr,rc)
          if (rc/=0) ierr=rc
          coeff(first:last,first:last,spin)=c
          eps(first:last,spin)=en
          singles(spin)=singles(spin)+corr
          call mma_deallocate(c)
          call mma_deallocate(f)
          call mma_deallocate(en)
        end do
      end do
    end if
    message='Symmetry-block semicanonicalization or singles failed.'
    call ROMP2_SyncError(ierr,message)
    if (ierr/=0) goto 900
    if (ROMP2_IsParallel()) then
      if (.not.ROMP2_IsRoot()) then
        coeff=0.0_wp
        eps=0.0_wp
        singles=0.0_wp
      end if
      call GADGOp(coeff(1,1,1),size(coeff),'+')
      call GADGOp(eps(1,1),size(eps),'+')
      call GADGOp(singles(1),2,'+')
    end if
  end if
  call mma_deallocate(hcore,safe='*')
  call mma_deallocate(overlap,safe='*')
  call mma_deallocate(density,safe='*')
  call mma_deallocate(fock,safe='*')
  call TransformPairs(ierr)
  if (ierr/=0) goto 900
  call Doubles(e,ierr)
  if (ierr/=0) goto 900
  corr=sum(e)+sum(singles)
  total=eref+corr
  if (ROMP2_Any(logical(.not.ieee_is_finite(total),kind=kind(.true.)))) then
    message='Nonfinite symmetry MP2 total energy.'
    goto 900
  end if
  if (cho_open) then
    call Cho_X_Final(ierr)
    cho_open=.false.
    message='Symmetry Cholesky finalization failed.'
    call ROMP2_SyncError(ierr,message)
    if (ierr/=0) goto 900
  end if
  if (ROMP2_IsRoot()) then
    write(u6,'(/,A,I3,A,I5)') ' Open-shell MP2 symmetry blocks: ',ns,'; MPI processes: ',workers
    if (cho) then
      write(u6,'(A)') ' Integral route: symmetry-blocked Cholesky'
    else
      write(u6,'(A)') ' Integral route: symmetry-blocked conventional'
    end if
    write(u6,'(A,8I6)') ' Basis functions per irrep: ',bs(1:ns)
    write(u6,'(A,8I6)') ' Alpha occupied per irrep: ',occ(1:ns,1)
    write(u6,'(A,8I6)') ' Beta occupied per irrep:  ',occ(1:ns,2)
    write(u6,'(A,8I6)') ' Frozen orbitals per irrep: ',frozen(1:ns)
    write(u6,'(A,8I8)') ' AO pairs per product irrep: ',blocks(1:ns)%np
    write(u6,'(A,8I8)') ' OV pairs per product irrep: ',blocks(1:ns)%nq
    write(u6,'(A)') ' Frozen convention: lowest occupied energies within each irrep and spin.'
    if (restricted) then
      write(u6,'(A,F24.12)') ' ROHF reference energy: ',eref
      write(u6,'(A,F24.12)') ' ROMP2 alpha singles:  ',singles(1)
      write(u6,'(A,F24.12)') ' ROMP2 beta singles:   ',singles(2)
      write(u6,'(A,F24.12)') ' ROMP2 AA doubles:     ',e(1)
      write(u6,'(A,F24.12)') ' ROMP2 BB doubles:     ',e(2)
      write(u6,'(A,F24.12)') ' ROMP2 AB doubles:     ',e(3)
      write(u6,'(A,F24.12)') ' ROMP2 correlation:    ',corr
      write(u6,'(A,F24.12)') ' ROMP2 total energy:   ',total
    else
      write(u6,'(A,F24.12)') ' UHF reference energy:     ',eref
      write(u6,'(A,F24.12)') ' UMP2 AA correlation:     ',e(1)
      write(u6,'(A,F24.12)') ' UMP2 BB correlation:     ',e(2)
      write(u6,'(A,F24.12)') ' UMP2 AB correlation:     ',e(3)
      write(u6,'(A,F24.12)') ' UMP2 correlation energy: ',corr
      write(u6,'(A,F24.12)') ' UMP2 total energy:       ',total
    end if
  end if
  call Store_Energies(1,[total],1)
  call Put_iScalar('mp2prpt',0)
  if (restricted) then
    call Put_cArray('Relax Method','ROMP2   ',8)
    call Add_Info('E_ROMP2',[total],1,8)
  else
    call Put_cArray('Relax Method','UMP2    ',8)
    call Add_Info('E_MP2',[total],1,8)
  end if
  ireturn=_RC_ALL_IS_WELL_
900 continue
  if (cho_open) call Cho_X_Final(rc)
  if ((ireturn/=_RC_ALL_IS_WELL_).and.ROMP2_IsRoot()) write(u6,'(A)') ' Symmetry MP2 error: '//trim(message)
  do g=1,8
    call mma_deallocate(blocks(g)%mu,safe='*')
    call mma_deallocate(blocks(g)%nu,safe='*')
    call mma_deallocate(blocks(g)%qi,safe='*')
    call mma_deallocate(blocks(g)%qa,safe='*')
    call mma_deallocate(blocks(g)%spin,safe='*')
    call mma_deallocate(blocks(g)%metric,safe='*')
    call mma_deallocate(blocks(g)%vectors,safe='*')
    call mma_deallocate(blocks(g)%projected,safe='*')
    call mma_deallocate(blocks(g)%kernel,safe='*')
    blocks(g)%np=0
    blocks(g)%nq=0
    blocks(g)%nlocal=0
    blocks(g)%nglobal=0
    blocks(g)%voff=0
  end do
  call mma_deallocate(coeff,safe='*')
  call mma_deallocate(eps,safe='*')
  call mma_deallocate(symorb,safe='*')
  call mma_deallocate(pairindex,safe='*')
  call mma_deallocate(qindex,safe='*')
  call mma_deallocate(hcore,safe='*')
  call mma_deallocate(overlap,safe='*')
  call mma_deallocate(density,safe='*')
  call mma_deallocate(fock,safe='*')
end subroutine Open_Shell_Symmetry_Driver

logical function MemoryOK(elements)
  real(kind=wp),intent(in) :: elements
  integer(kind=iwp) :: available
  call mma_maxDBLE(available)
  MemoryOK=.not.ROMP2_Any((elements>0.65_wp*real(available,wp)).or. &
                       (elements>real(huge(available),wp)/real(storage_size(0.0_wp),wp)))
  if (.not.MemoryOK) message='Insufficient memory or integer range for in-core symmetry MP2.'
end function MemoryOK

subroutine ReadCounts(label,values,ierr)
  character(len=*),intent(in) :: label
  integer(kind=iwp),intent(out) :: values(8),ierr
  integer(kind=iwp) :: n
  logical(kind=iwp) :: found
  ierr=1
  values=0
  call Qpg_iArray(label,found,n)
  if (.not.found) return
  if ((n<ns).or.(n>8)) return
  call Get_iArray(label,values,n)
  if (any(values(1:ns)<0)) return
  ierr=0
end subroutine ReadCounts

subroutine ReadReference(ierr)
  integer(kind=iwp),intent(out) :: ierr
  integer(kind=iwp) :: norb(8),nf(8),ndel(8),active(8),mult,nel,s,spin,i,j,k,l,n
  logical(kind=iwp) :: found,df
  real(kind=wp),allocatable :: buffer(:)
  character(len=20) :: label
  ierr=1
  call Get_iScalar('nSym',ns)
  if ((ns/=2).and.(ns/=4).and.(ns/=8)) return
  call DecideOnCholesky(cho)
  call DecideOnDF(df)
  message='Density fitting is not implemented in symmetry open-shell MP2.'
  if (df) return
  message='Missing or invalid symmetry reference counts.'
  call ReadCounts('nBas',bs,ierr)
  if (ierr/=0) return
  call ReadCounts('nIsh',occ(:,1),ierr)
  if (ierr/=0) return
  call ReadCounts('nFro',nf,ierr)
  if (ierr/=0) return
  ierr=1
  message='Reference-frozen orbitals are unsupported; freeze only in MBPT2.'
  if (any(nf(1:ns)/=0)) return
  if (restricted) then
    call ReadCounts('nAsh',active,ierr)
    if (ierr/=0) return
    call ReadCounts('nDel',ndel,ierr)
    if (ierr/=0) return
    ierr=1
    message='ROMP2 symmetry requires full orbital spaces and a maximum-spin CAS(n,n) determinant.'
    if (any(ndel(1:ns)/=0)) return
    call Get_iScalar('Multiplicity',mult)
    call Get_iScalar('nActel',nel)
    if ((nel/=sum(active(1:ns))).or.(mult/=nel+1)) return
    occ(:,2)=occ(:,1)
    occ(:,1)=occ(:,1)+active
    call Get_dScalar('Last energy',eref)
  else
    call ReadCounts('nIsh_ab',occ(:,2),ierr)
    if (ierr/=0) return
    call ReadCounts('nOrb',norb,ierr)
    if (ierr/=0) return
    ierr=1
    message='Symmetry UMP2 currently requires no deleted or linearly dependent basis directions.'
    if (any(norb(1:ns)/=bs(1:ns))) return
    call Get_dScalar('SCF energy',eref)
  end if
  ierr=1
  if (.not.ieee_is_finite(eref)) return
  if (any(occ(1:ns,1)>bs(1:ns)).or.any(occ(1:ns,2)>bs(1:ns))) return
  nb=sum(bs(1:ns))
  if (nb<=0) return
  if (.not.MemoryOK(5.0_wp*real(nb,wp)**2+3.0_wp*real(nb,wp))) return
  off=0
  do s=2,ns
    off(s)=off(s-1)+bs(s-1)
  end do
  call mma_allocate(coeff,nb,nb,2,label='symmetry MO coefficients')
  call mma_allocate(eps,nb,2,label='symmetry MO energies')
  call mma_allocate(symorb,nb,label='orbital irreps')
  coeff=0.0_wp
  eps=0.0_wp
  do s=1,ns
    symorb(off(s)+1:off(s)+bs(s))=s
  end do
  do spin=1,2
    label='SCF orbitals'
    if (spin==2) label='SCF orbitals_ab'
    if (restricted) label='RASSCF orbitals'
    call Qpg_dArray(trim(label),found,n)
    message='Missing or inconsistent symmetry MO record.'
    if (.not.found) return
    if (n/=sum(bs(1:ns)**2)) return
    call mma_allocate(buffer,n,label='symmetry coefficient input')
    call Get_dArray(trim(label),buffer,n)
    if (.not.all(ieee_is_finite(buffer))) then
      call mma_deallocate(buffer)
      return
    end if
    k=0
    do s=1,ns
      do j=1,bs(s)
        do i=1,bs(s)
          k=k+1
          coeff(off(s)+i,off(s)+j,spin)=buffer(k)
        end do
      end do
    end do
    call mma_deallocate(buffer)
    if (restricted) cycle
    label='OrbE'
    if (spin==2) label='OrbE_ab'
    call Qpg_dArray(trim(label),found,l)
    if (.not.found) return
    if (l/=nb) return
    call Get_dArray(trim(label),eps(:,spin),nb)
    if (.not.all(ieee_is_finite(eps(:,spin)))) return
  end do
  ierr=0
end subroutine ReadReference

subroutine SynchronizeReference()
  ! Explicitly replicate root reference data before any MPI-distributed
  ! Fock, transformation, or doubles work.  The C1 drivers are untouched;
  ! this routine is used only by the non-C1 symmetry path.
  real(kind=wp) :: scalar(1)
  if (.not.ROMP2_IsParallel()) return
  scalar=0.0_wp
  if (ROMP2_IsRoot()) then
    scalar(1)=eref
  else
    coeff=0.0_wp
    eps=0.0_wp
  end if
  call GADGOp(coeff(1,1,1),size(coeff),'+')
  call GADGOp(eps(1,1),size(eps),'+')
  call GADGOp(scalar(1),1,'+')
  eref=scalar(1)
end subroutine SynchronizeReference

subroutine SynchronizeFrozen()
  ! Parse MBPT2 input on every rank for normal error handling, then make
  ! rank zero's accepted counts authoritative for all subsequent maps.
  if (.not.ROMP2_IsParallel()) return
  if (.not.ROMP2_IsRoot()) frozen(1:ns)=0
  call gaIgOP(frozen(1),ns,'+')
end subroutine SynchronizeFrozen

subroutine ReadFrozen(ierr)
  use spool, only: SpoolInp,Close_LuSpool
  integer(kind=iwp),intent(out) :: ierr
  integer(kind=iwp) :: lu,eq,k,start,n,ios,value,s
  character(len=180) :: raw,line,values,token
  character(len=4) :: command
  character(len=180),external :: Get_Ln
  logical :: seen
  ierr=1
  frozen=0
  seen=.false.
  message='Specify Frozen with one nonnegative count per irrep in RUNFILE order.'
  lu=17
  call SpoolInp(lu)
  rewind(lu)
  call RdNLst(lu,'MBPT2')
  do
    raw=Get_Ln(lu)
    line=raw
    call StdFmt(line,command)
    if (command=='END ') exit
    if ((command/='FROZ').or.seen) goto 900
    eq=index(raw,'=')
    values=''
    if (eq>0) values=raw(eq+1:)
    if (len_trim(values)==0) values=Get_Ln(lu)
    ! Strict parsing: reject missing, repeated, negative, or extra counts.
    n=0
    k=1
    do while (k<=len_trim(values))
      if ((values(k:k)==' ').or.(values(k:k)==achar(9)).or.(values(k:k)==',')) then
        k=k+1
        cycle
      end if
      start=k
      do while (k<=len_trim(values))
        if (index('0123456789',values(k:k))==0) exit
        k=k+1
      end do
      if (k==start) goto 900
      n=n+1
      if (n>ns) goto 900
      token=values(start:k-1)
      value=0
      read(token,*,iostat=ios) value
      if ((ios/=0).or.(value<0)) goto 900
      frozen(n)=value
    end do
    if (n/=ns) goto 900
    do s=1,ns
      if (frozen(s)>minval(occ(s,:))) goto 900
    end do
    seen=.true.
  end do
  if (seen) ierr=0
900 continue
  call Close_LuSpool(lu)
end subroutine ReadFrozen

subroutine MakePairs(ierr)
  integer(kind=iwp),intent(out) :: ierr
  integer(kind=iwp) :: g,s,t,i,j,p,a,spin,mu,nu
  ierr=1
  if (.not.MemoryOK(5.0_wp*real(nb,wp)**2)) return
  call mma_allocate(pairindex,nb,nb,label='AO pair map')
  call mma_allocate(qindex,nb,nb,2,label='OV pair map')
  pairindex=0
  qindex=0
  do g=1,ns
    p=0
    do s=1,ns
      t=Mul(g,s)
      if (s<t) cycle
      if (s==t) then
        p=p+bs(s)*(bs(s)+1)/2
      else
        p=p+bs(s)*bs(t)
      end if
    end do
    blocks(g)%np=p
    call mma_allocate(blocks(g)%mu,p,label='pair first AO')
    call mma_allocate(blocks(g)%nu,p,label='pair second AO')
    p=0
    do s=1,ns
      t=Mul(g,s)
      if (s<t) cycle
      do i=1,bs(s)
        do j=1,bs(t)
          if ((s==t).and.(j>i)) cycle
          p=p+1
          mu=off(s)+i
          nu=off(t)+j
          blocks(g)%mu(p)=mu
          blocks(g)%nu(p)=nu
          pairindex(mu,nu)=p
          pairindex(nu,mu)=p
        end do
      end do
    end do
    p=0
    do spin=1,2
      do s=1,ns
        t=Mul(g,s)
        p=p+(occ(s,spin)-frozen(s))*(bs(t)-occ(t,spin))
      end do
    end do
    blocks(g)%nq=p
    call mma_allocate(blocks(g)%qi,p,label='pair occupied MO')
    call mma_allocate(blocks(g)%qa,p,label='pair virtual MO')
    call mma_allocate(blocks(g)%spin,p,label='pair spin')
    p=0
    do spin=1,2
      do s=1,ns
        t=Mul(g,s)
        do i=frozen(s)+1,occ(s,spin)
          do a=occ(t,spin)+1,bs(t)
            p=p+1
            blocks(g)%qi(p)=off(s)+i
            blocks(g)%qa(p)=off(t)+a
            blocks(g)%spin(p)=spin
            qindex(off(s)+i,off(t)+a,spin)=p
          end do
        end do
      end do
    end do
  end do
  ierr=0
end subroutine MakePairs

subroutine ReadConventional(ierr)
  use TwoDat, only: AuxTwo
  integer(kind=iwp),intent(out) :: ierr
  integer(kind=iwp) :: rc,rc2,lu,fs,fb(8),skip(8),g,s,t,u,v,np,nq,p,opt,nmat,p0,q0
  logical(kind=iwp) :: square,exists
  logical :: opened
  real(kind=wp),allocatable :: buffer(:)
  real(kind=wp) :: needed
  ierr=1
  opened=.false.
  needed=0.0_wp
  do g=1,ns
    needed=needed+real(blocks(g)%np,wp)**2
  end do
  if (.not.MemoryOK(needed+real(nb,wp)**2+4.0_wp)) return
  do g=1,ns
    call mma_allocate(blocks(g)%metric,blocks(g)%np,blocks(g)%np,label='symmetry AO pair metric')
    blocks(g)%metric=0.0_wp
  end do
  call mma_allocate(buffer,nb*nb+1,label='ORDINT symmetry row')

  ! In a real MPI calculation ORDINT is distributed over the ranks.
  ! Therefore every rank must open and read its local ORDINT contribution.
  ! Each packed AO-pair row is reconstructed by a global sum before it is
  ! copied into the replicated symmetry-block metric.
  rc=0
  call f_Inquire('ORDINT',exists)
  if ((.not.exists).or.AuxTwo%Opn) rc=1
  if (ROMP2_Any(rc/=0)) goto 710

  lu=43
  call OpnOrd(rc,0,'ORDINT',lu)
  if (rc==0) opened=.true.
  if (ROMP2_Any(rc/=0)) goto 710

  call GetOrd(rc,square,fs,fb,skip)
  if (rc==0) then
    if ((fs/=ns).or.any(fb(1:ns)/=bs(1:ns)).or.any(skip(1:ns)/=0)) rc=1
  end if
  if (ROMP2_Any(rc/=0)) goto 710

  do g=1,ns
    do s=1,ns
      t=Mul(g,s)
      if ((s<t).or.(bs(s)*bs(t)==0)) cycle
      np=bs(s)*bs(t)
      if (s==t) np=bs(s)*(bs(s)+1)/2
      p0=pairindex(off(s)+1,off(t)+1)-1
      do u=1,ns
        v=Mul(g,u)
        if ((u<v).or.(bs(u)*bs(v)==0)) cycle
        ! Packed ORDINT stores only the lower symmetry-pair triangle.
        if ((.not.square).and.(s*(s-1)/2+t<u*(u-1)/2+v)) cycle
        nq=bs(u)*bs(v)
        if (u==v) nq=bs(u)*(bs(u)+1)/2
        q0=pairindex(off(u)+1,off(v)+1)-1
        opt=1
        do p=1,np
          buffer=0.0_wp
          rc=0
          call RdOrd_(rc,opt,s,t,u,v,buffer,nq+1,nmat)
          opt=2
          if (rc==0) then
            if (nmat/=1) rc=1
            if (.not.all(ieee_is_finite(buffer(1:nq)))) rc=1
          end if
          if (ROMP2_Any(rc/=0)) goto 710

          if (ROMP2_IsParallel()) call GADGOp(buffer(1),nq,'+')

          blocks(g)%metric(p0+p,q0+1:q0+nq)=buffer(1:nq)
          if (.not.square) blocks(g)%metric(q0+1:q0+nq,p0+p)=buffer(1:nq)
        end do
      end do
    end do
  end do

710 continue
  rc2=0
  if (opened) then
    call ClsOrd(rc2)
    opened=.false.
  end if
  if (rc==0 .and. rc2/=0) rc=rc2

  call mma_deallocate(buffer)
  message='Symmetry ORDINT read failed or dimensions disagree with the reference.'
  call ROMP2_SyncError(rc,message)
  if (rc/=0) return

  ierr=0
end subroutine ReadConventional

subroutine ReadCholesky(ierr,initialized)
  use Cholesky, only: NumCho,nDimRS,ChoNSym=>nSym,ChoNBas=>nBas
  use Data_Structures, only: SBA_Type,Allocate_DT,Deallocate_DT
  integer(kind=iwp),intent(out) :: ierr
  logical,intent(out) :: initialized
  integer(kind=iwp) :: g,s,t,p,j,rc,nread,redset,skip(8),pos,mu,nu,iswap
  integer(kind=iwp),allocatable :: counts(:,:)
  real(kind=wp),allocatable :: red(:)
  real(kind=wp) :: needed
  type(SBA_Type),target :: ao
#include "warnings.h"
  initialized=.false.
  call Cho_X_Init(rc,0.0_wp)
  if (ROMP2_Any(rc/=0)) then
    ! Initialization can fail partially: do not finalize on a subset.
    call Quit(_RC_CHO_INI_)
    ierr=1
    return
  end if
  initialized=.true.
  ierr=1
  message='Cholesky symmetry dimensions do not match the reference.'
  if (ROMP2_Any(ChoNSym/=ns)) return
  if (ROMP2_Any(any(ChoNBas(1:ns)/=bs(1:ns)))) return
  call mma_allocate(counts,ns,workers,label='symmetry Cholesky ownership')
  counts=0
  counts(:,rankid+1)=NumCho(1:ns)
  if (ROMP2_IsParallel()) call gaIgOP(counts(1,1),size(counts),'+')
  if (any(counts<0).or.(sum(real(counts,wp))>real(huge(pos)-1,wp))) then
    call mma_deallocate(counts)
    return
  end if
  needed=0.0_wp
  do g=1,ns
    blocks(g)%nlocal=counts(g,rankid+1)
    blocks(g)%nglobal=sum(counts(g,:))
    blocks(g)%voff=sum(counts(g,1:rankid))
    needed=needed+real(blocks(g)%np,wp)*real(blocks(g)%nlocal,wp)
  end do
  call mma_deallocate(counts)
  message='No Cholesky vectors found in any symmetry sector.'
  if (sum(blocks(1:ns)%nglobal)<=0) return
  if (ROMP2_Any(.not.allocated(nDimRS))) return
  nread=max(1,maxval(nDimRS))
  if (.not.MemoryOK(needed+real(nread,wp)+real(nb,wp)**2)) return
  call mma_allocate(red,nread,label='symmetry reduced vector')
  skip=1
  do g=1,ns
    call mma_allocate(blocks(g)%vectors,blocks(g)%np,blocks(g)%nlocal,label='local symmetry AO vectors')
    call Allocate_DT(ao,bs(1:ns),bs(1:ns),1,g,ns,0,Label='symmetry full vector')
    ! Allocate_DT with iCase=0 lays out every symmetry block SB(s)
    ! consecutively in A0, including both orientations for g > 1.
    ! Cho_X_getVfull(iSwap=0) writes only the unique s >= Mul(g,s)
    ! blocks, but its ipOff values must still point to the actual starts
    ! of those blocks in the full iCase=0 allocation.  Therefore advance
    ! over every allocated SB(s), not only the unique Cholesky blocks.
    pos=1
    do s=1,ns
      t=Mul(g,s)
      ao%ipOff(s)=pos
      pos=pos+bs(s)*bs(t)
    end do
    redset=-1
    rc=0
    iswap=0
    if (g==1) iswap=2
    do j=1,blocks(g)%nlocal
      call Cho_X_getVfull(rc,red,nread,j,1,g,iswap,redset,ao,skip,.true.)
      if (rc/=0) exit
      do p=1,blocks(g)%np
        mu=blocks(g)%mu(p)
        nu=blocks(g)%nu(p)
        s=symorb(mu)
        t=symorb(nu)
        blocks(g)%vectors(p,j)=ao%SB(s)%A3(mu-off(s),nu-off(t),1)
      end do
      if (.not.all(ieee_is_finite(blocks(g)%vectors(:,j)))) then
        rc=1
        exit
      end if
    end do
    call Deallocate_DT(ao)
    message='Symmetry Cholesky vector read failed.'
    ! Rank-local vector loops have unequal lengths; synchronize afterwards.
    call ROMP2_SyncError(rc,message)
    if (rc/=0) exit
  end do
  call mma_deallocate(red)
  if (rc/=0) return
  ierr=0
end subroutine ReadCholesky

subroutine ReadOne(ierr)
  use OneDat, only: AuxOne,sNoNuc,sNoOri
  integer(kind=iwp),intent(out) :: ierr
  integer(kind=iwp) :: rc,rc2,opt,comp,sy,s,i,j,k,pass,np
  logical :: opened
  real(kind=wp),allocatable :: buffer(:)
  character(len=8) :: label
  np=sum(bs(1:ns)*(bs(1:ns)+1)/2)
  call mma_allocate(buffer,np+4,label='symmetry ONEINT buffer')
  hcore=0.0_wp
  overlap=0.0_wp
  rc=0
  opened=.false.
  if (ROMP2_IsRoot()) then
    if (.not.AuxOne%Opn) then
      call OpnOne(rc,0,'ONEINT',41)
      opened=.true.
      if (rc/=0) goto 710
    end if
    do pass=1,2
      opt=ibset(ibset(0,sNoOri),sNoNuc)
      comp=1
      sy=1
      label='OneHam  '
      if (pass==2) label='Mltpl  0'
      call RdOne(rc,opt,label,comp,buffer,sy)
      if (rc/=0) exit
      if (.not.all(ieee_is_finite(buffer(1:np)))) then
        rc=1
        exit
      end if
      k=0
      do s=1,ns
        do i=1,bs(s)
          do j=1,i
            k=k+1
            if (pass==1) then
              hcore(off(s)+i,off(s)+j)=buffer(k)
              hcore(off(s)+j,off(s)+i)=buffer(k)
            else
              overlap(off(s)+i,off(s)+j)=buffer(k)
              overlap(off(s)+j,off(s)+i)=buffer(k)
            end if
          end do
        end do
      end do
    end do
710 continue
    if (opened) then
      call ClsOne(rc2,0)
      if (rc2/=0) rc=rc2
    end if
  end if
  call mma_deallocate(buffer)
  message='Symmetry ONEINT read failed.'
  call ROMP2_SyncError(rc,message)
  ierr=rc
  if (rc/=0) return
  if (ROMP2_IsParallel()) then
    call GADGOp(hcore(1,1),size(hcore),'+')
    call GADGOp(overlap(1,1),size(overlap),'+')
  end if
end subroutine ReadOne

real(kind=wp) function AOIntegral(i,j,k,l) result(value)
  integer(kind=iwp),intent(in) :: i,j,k,l
  integer(kind=iwp) :: g
  value=0.0_wp
  g=Mul(symorb(i),symorb(j))
  if (g/=Mul(symorb(k),symorb(l))) return
  value=blocks(g)%metric(pairindex(i,j),pairindex(k,l))
end function AOIntegral

subroutine ReferenceFock(ierr)
  integer(kind=iwp),intent(out) :: ierr
  integer(kind=iwp) :: s,t,g,spin,n,m,first,last,other,i,j,k,l,p,v,mu,nu
  real(kind=wp),allocatable :: work(:,:),mat(:,:),vec(:,:),exchange(:,:)
  real(kind=wp) :: charge,enuc,check,grad,oi,oj
  ierr=1
  if (.not.MemoryOK(11.0_wp*real(nb,wp)**2+4.0_wp)) return
  call mma_allocate(hcore,nb,nb,label='symmetry core Hamiltonian')
  call mma_allocate(overlap,nb,nb,label='symmetry overlap')
  call mma_allocate(density,nb,nb,2,label='symmetry spin densities')
  call mma_allocate(fock,nb,nb,2,label='symmetry spin Fock')
  call mma_allocate(work,nb,nb,label='symmetry Fock work')
  call mma_allocate(mat,nb,nb,label='symmetry MO Fock')
  call mma_allocate(vec,nb,nb,label='symmetry AO vector')
  call mma_allocate(exchange,nb,nb,label='symmetry exchange')
  call ReadOne(ierr)
  if (ierr/=0) goto 900
  ierr=1
  density=0.0_wp
  ! Overlap validation and densities are independently symmetry-blocked.
  grad=0.0_wp
  do s=1,ns
    n=bs(s)
    if (n==0) cycle
    first=off(s)+1
    last=off(s)+n
    call DGEMM_('N','N',n,n,n,1.0_wp,overlap(first,first),nb,coeff(first,first,1),nb,0.0_wp,work,nb)
    call DGEMM_('T','N',n,n,n,1.0_wp,coeff(first,first,1),nb,work,nb,0.0_wp,mat,nb)
    do i=1,n
      mat(i,i)=mat(i,i)-1.0_wp
    end do
    grad=max(grad,maxval(abs(mat(1:n,1:n))))
    do spin=1,2
      call DGEMM_('N','T',n,n,occ(s,spin),1.0_wp,coeff(first,first,spin),nb, &
                  coeff(first,first,spin),nb,0.0_wp,density(first,first,spin),nb)
    end do
  end do
  message='Symmetry ROHF orbitals are not orthonormal.'
  if (ROMP2_Any(logical((.not.ieee_is_finite(grad)).or.(grad>1.0e-7_wp),kind=kind(.true.)))) goto 900
  fock=0.0_wp
  if (ROMP2_IsRoot()) then
    fock(:,:,1)=hcore
    fock(:,:,2)=hcore
  end if
  if (cho) then
    do g=1,ns
      do v=1,blocks(g)%nlocal
        vec=0.0_wp
        do p=1,blocks(g)%np
          mu=blocks(g)%mu(p)
          nu=blocks(g)%nu(p)
          vec(mu,nu)=blocks(g)%vectors(p,v)
          vec(nu,mu)=blocks(g)%vectors(p,v)
        end do
        if (g==1) then
          charge=sum((density(:,:,1)+density(:,:,2))*vec)
          fock(:,:,1)=fock(:,:,1)+charge*vec
          fock(:,:,2)=fock(:,:,2)+charge*vec
        end if
        do s=1,ns
          t=Mul(g,s)
          n=bs(s)
          m=bs(t)
          if (n*m==0) cycle
          first=off(s)+1
          last=off(s)+n
          other=off(t)+1
          do spin=1,2
            call DGEMM_('N','N',n,m,m,1.0_wp,vec(first,other),nb,density(other,other,spin),nb,0.0_wp,work,nb)
            call DGEMM_('N','T',n,n,m,1.0_wp,work,nb,vec(first,other),nb,0.0_wp,exchange,nb)
            fock(first:last,first:last,spin)=fock(first:last,first:last,spin)-exchange(1:n,1:n)
          end do
        end do
      end do
    end do
  else
    ! The conventional AO metric is replicated.  Build the small reference
    ! Fock deterministically on rank zero, exactly as in the serial path,
    ! then distribute it below.  MPI work-sharing remains in TransformPairs
    ! and Doubles, where it materially reduces the MP2 cost.
    if (ROMP2_IsRoot()) then
      do i=1,nb
        s=symorb(i)
        do j=off(s)+1,off(s)+bs(s)
          do k=1,nb
            t=symorb(k)
            do l=off(t)+1,off(t)+bs(t)
              charge=(density(k,l,1)+density(k,l,2))*AOIntegral(i,j,k,l)
              fock(i,j,1)=fock(i,j,1)+charge-density(k,l,1)*AOIntegral(i,k,j,l)
              fock(i,j,2)=fock(i,j,2)+charge-density(k,l,2)*AOIntegral(i,k,j,l)
            end do
          end do
        end do
      end do
    end if
  end if
  if (ROMP2_IsParallel()) call GADGOp(fock(1,1,1),size(fock),'+')
  message='Nonfinite symmetry reference Fock matrix.'
  if (ROMP2_Any(logical(.not.all(ieee_is_finite(fock)),kind=kind(.true.)))) goto 900
  call Get_dScalar('PotNuc',enuc)
  check=enuc+0.5_wp*(sum(density(:,:,1)*(hcore+fock(:,:,1)))+sum(density(:,:,2)*(hcore+fock(:,:,2))))
  message='Reconstructed symmetry ROHF energy differs from RASSCF.'
  if (ROMP2_Any(logical((.not.ieee_is_finite(check)).or.(abs(check-eref)>1.0e-7_wp), &
                       kind=kind(.true.)))) goto 900
  ! Stationarity test in the shared spatial ROHF basis, before rotation.
  grad=0.0_wp
  vec=0.0_wp
  do s=1,ns
    n=bs(s)
    if (n==0) cycle
    first=off(s)+1
    do spin=1,2
      call DGEMM_('N','N',n,n,n,1.0_wp,fock(first,first,spin),nb,coeff(first,first,1),nb,0.0_wp,work,nb)
      call DGEMM_('T','N',n,n,n,1.0_wp,coeff(first,first,1),nb,work,nb,0.0_wp,mat,nb)
      do i=1,n
        oi=0.0_wp
        if (i<=occ(s,spin)) oi=1.0_wp
        do j=1,n
          oj=0.0_wp
          if (j<=occ(s,spin)) oj=1.0_wp
          vec(off(s)+i,off(s)+j)=vec(off(s)+i,off(s)+j)+(oi-oj)*mat(i,j)
        end do
      end do
    end do
  end do
  grad=maxval(abs(vec))
  message='ROHF symmetry stationarity check failed; tighten RASSCF convergence.'
  if (ROMP2_Any(logical((.not.ieee_is_finite(grad)).or.(grad>1.0e-5_wp),kind=kind(.true.)))) goto 900
  ierr=0
900 continue
  call mma_deallocate(work)
  call mma_deallocate(mat)
  call mma_deallocate(vec)
  call mma_deallocate(exchange)
end subroutine ReferenceFock

subroutine TransformPairs(ierr)
  integer(kind=iwp),intent(out) :: ierr
  integer(kind=iwp) :: g,p,q,mu,nu,i,a,spin,np,nq,nlocal,nglobal,start
  real(kind=wp),allocatable :: projection(:,:),work(:,:)
  real(kind=wp) :: required
  ierr=1
  do g=1,ns
    np=blocks(g)%np
    nq=blocks(g)%nq
    nlocal=blocks(g)%nlocal
    nglobal=blocks(g)%nglobal
    required=real(np,wp)*real(nq+1,wp)
    if (cho) then
      required=required+real(nq,wp)*real(nglobal,wp)
    else
      required=required+real(nq,wp)**2
    end if
    if (.not.MemoryOK(required)) return
    call mma_allocate(projection,np,nq,label='symmetry AO to OV projection')
    projection=0.0_wp
    do q=1,nq
      i=blocks(g)%qi(q)
      a=blocks(g)%qa(q)
      spin=blocks(g)%spin(q)
      do p=1,np
        mu=blocks(g)%mu(p)
        nu=blocks(g)%nu(p)
        ! An AO pair contributes only to the corresponding irrep pair.
        if (((symorb(mu)/=symorb(a)).or.(symorb(nu)/=symorb(i))).and. &
            ((symorb(nu)/=symorb(a)).or.(symorb(mu)/=symorb(i)))) cycle
        projection(p,q)=coeff(mu,a,spin)*coeff(nu,i,spin)
        if (mu/=nu) projection(p,q)=projection(p,q)+coeff(nu,a,spin)*coeff(mu,i,spin)
      end do
    end do
    if (cho) then
      call mma_allocate(blocks(g)%projected,nq,nglobal,label='global symmetry OV vectors')
      blocks(g)%projected=0.0_wp
      start=blocks(g)%voff+1
      if ((np>0).and.(nq>0).and.(nlocal>0)) then
        call DGEMM_('T','N',nq,nlocal,np,1.0_wp,projection,np,blocks(g)%vectors,np, &
                    0.0_wp,blocks(g)%projected(1,start),nq)
      end if
      if (ROMP2_IsParallel().and.(nq*nglobal>0)) then
        call GADGOp(blocks(g)%projected(1,1),size(blocks(g)%projected),'+')
      end if
      call mma_deallocate(blocks(g)%vectors)
    else
      call mma_allocate(blocks(g)%kernel,nq,nq,label='symmetry OV integral kernel')
      call mma_allocate(work,np,1,label='symmetry metric projection')
      blocks(g)%kernel=0.0_wp
      if ((np>0).and.(nq>0)) then
        do q=1,nq
          if (mod(q-1,workers)/=rankid) cycle
          call DGEMM_('N','N',np,1,np,1.0_wp,blocks(g)%metric,np,projection(1,q),np,0.0_wp,work,np)
          call DGEMM_('T','N',nq,1,np,1.0_wp,projection,np,work,np,0.0_wp,blocks(g)%kernel(1,q),nq)
        end do
      end if
      if (ROMP2_IsParallel().and.(nq>0)) call GADGOp(blocks(g)%kernel(1,1),size(blocks(g)%kernel),'+')
      call mma_deallocate(work)
      call mma_deallocate(blocks(g)%metric)
    end if
    call mma_deallocate(projection)
  end do
  ierr=0
end subroutine TransformPairs

logical function IsOccupied(i,spin)
  integer(kind=iwp),intent(in) :: i,spin
  integer(kind=iwp) :: s
  s=symorb(i)
  IsOccupied=(i-off(s)>frozen(s)).and.(i-off(s)<=occ(s,spin))
end function IsOccupied

logical function IsVirtual(a,spin)
  integer(kind=iwp),intent(in) :: a,spin
  IsVirtual=a-off(symorb(a))>occ(symorb(a),spin)
end function IsVirtual

real(kind=wp) function MOIntegral(i,a,spin,j,b,other) result(value)
  integer(kind=iwp),intent(in) :: i,a,spin,j,b,other
  integer(kind=iwp) :: g,p,q
  real(kind=wp),external :: ddot_
  g=Mul(symorb(i),symorb(a))
  value=0.0_wp
  if (g/=Mul(symorb(j),symorb(b))) return
  p=qindex(i,a,spin)
  q=qindex(j,b,other)
  if (cho) then
    if (blocks(g)%nglobal>0) value=ddot_(blocks(g)%nglobal, &
        blocks(g)%projected(p,1),blocks(g)%nq,blocks(g)%projected(q,1),blocks(g)%nq)
  else
    value=blocks(g)%kernel(p,q)
  end if
end function MOIntegral

subroutine Doubles(energy,ierr)
  real(kind=wp),intent(out) :: energy(3)
  integer(kind=iwp),intent(out) :: ierr
  integer(kind=iwp) :: spin,other,i,j,a,b,component,paircount
  real(kind=wp) :: den,value
  energy=0.0_wp
  ierr=0
  do component=1,3
    spin=1
    other=1
    if (component==2) then
      spin=2
      other=2
    else if (component==3) then
      other=2
    end if
    paircount=0
    do i=1,nb
      if (.not.IsOccupied(i,spin)) cycle
      do j=1,nb
        if (.not.IsOccupied(j,other)) cycle
        if ((spin==other).and.(j>=i)) cycle
        paircount=paircount+1
        if (mod(paircount-1,workers)/=rankid) cycle
        do a=1,nb
          if (.not.IsVirtual(a,spin)) cycle
          do b=1,nb
            if (.not.IsVirtual(b,other)) cycle
            if ((spin==other).and.(b>=a)) cycle
            if (Mul(symorb(i),symorb(a))/=Mul(symorb(j),symorb(b))) cycle
            den=eps(i,spin)+eps(j,other)-eps(a,spin)-eps(b,other)
            if ((.not.ieee_is_finite(den)).or.(abs(den)<=1.0e-12_wp)) then
              ierr=1
              cycle
            end if
            value=MOIntegral(i,a,spin,j,b,other)
            if (spin==other) value=value-MOIntegral(i,b,spin,j,a,spin)
            energy(component)=energy(component)+value*value/den
          end do
        end do
      end do
    end do
  end do
  if (.not.all(ieee_is_finite(energy))) ierr=1
  message='Invalid denominator or nonfinite symmetry doubles energy.'
  call ROMP2_SyncError(ierr,message)
  if (ierr/=0) return
  if (ROMP2_IsParallel()) call GADGOp(energy(1),3,'+')
end subroutine Doubles
end module Open_Shell_Symmetry
