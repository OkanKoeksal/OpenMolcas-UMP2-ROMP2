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

module ROMP2_Reference
use Definitions, only: wp, iwp, u6
use UMP2_Global
use ROMP2_Cholesky, only: ROMP2_Cho_Fock
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
private
public :: ROMP2_Read, ROMP2_Fock
contains
subroutine ROMP2_Read(ierr,Message)
  use Para_Info, only: nProcs
  integer(kind=iwp), intent(out) :: ierr
  character(len=*), intent(out) :: Message
  integer(kind=iwp) :: nd(1),ns(1),nf(1),nv(1),nb(1),mult,nel,l
  logical(kind=iwp) :: found,df
  character(len=8) :: method
  ierr=1
  Message='ROMP2 requires serial C1 high-spin determinant RASSCF with conventional or Cholesky integrals.'
  if (nProcs/=1) return
  call Get_cArray('Relax Method',method,8)
  if ((method/='CASSCF  ').and.(method/='RASSCF  ')) return
  call DecideOnCholesky(DoCholesky)
  call DecideOnDF(df)
  if (df) then
    Message='ROMP2 RI/DF integrals are unsupported; use Cholesky or conventional integrals.'
    return
  end if
  call Get_iScalar('nSym',nSym)
  if (nSym/=1) return
  call Get_iArray('nBas',nb,1)
  call Get_iArray('nIsh',nd,1)
  call Get_iArray('nAsh',ns,1)
  call Get_iArray('nFro',nf,1)
  call Get_iArray('nDel',nv,1)
  call Get_iScalar('Multiplicity',mult)
  call Get_iScalar('nActel',nel)
  ! CAS(nOpen,nOpen), maximum spin: exactly one determinant.
  ! Frozen/deleted RASSCF orbitals are excluded in this first version.
  if ((nf(1)/=0).or.(nv(1)/=0)) return
  if ((ns(1)<0).or.(nd(1)<0).or.(nel/=ns(1)).or.(mult/=ns(1)+1)) return
  nBas=nb(1)
  nOrb=nBas
  nOccA=nd(1)+ns(1)
  nOccB=nd(1)
  if ((nBas<1).or.(nOccA>nBas)) return
  if (real(nBas,wp)**2>=real(huge(l),wp)/2.0_wp) return
  call qpg_dArray('RASSCF orbitals',found,l)
  if (.not.found) return
  if (l/=nBas*nBas) return
  call mma_allocate(CAlpha,nBas,nOrb,label='ROMP2 alpha orbitals')
  call mma_allocate(CBeta,nBas,nOrb,label='ROMP2 beta orbitals')
  call mma_allocate(EOrbA,nOrb,label='ROMP2 alpha energies')
  call mma_allocate(EOrbB,nOrb,label='ROMP2 beta energies')
  call Get_dArray('RASSCF orbitals',CAlpha,nBas*nBas)
  CBeta=CAlpha
  call Get_dScalar('Last energy',ESCF)
  if ((.not.all(ieee_is_finite(CAlpha))).or.(.not.ieee_is_finite(ESCF))) return
  ReferenceMethod='ROHF    '
  ierr=0
  Message=''
end subroutine ROMP2_Read

subroutine ROMP2_Fock(FA,FB,ierr,Message)
  use OneDat, only: AuxOne,sNoNuc,sNoOri
  use TwoDat, only: AuxTwo
  real(kind=wp), intent(out) :: FA(nBas,nBas),FB(nBas,nBas)
  integer(kind=iwp), intent(out) :: ierr
  character(len=*), intent(out) :: Message
  real(kind=wp), allocatable :: ERI(:,:,:,:),H(:,:),S(:,:),DA(:,:),DB(:,:),Buf(:),T(:,:),MA(:,:),MB(:,:)
  real(kind=wp) :: required,enuc,echeck,grad,oa,ob,pa,pb
  integer(kind=iwp) :: avail,np,rc,rc2,lu,comp,sym,opt,i,j,k,l,q,nmat,fs,fbasis(8),skip(8)
  logical(kind=iwp) :: exists,square
  logical :: opened,oneopened
  character(len=8) :: label
  ierr=1
  Message='ROMP2 Fock construction failed.'
  opened=.false.
  oneopened=.false.
  required=14.0_wp*real(nBas,wp)**2+131072.0_wp
  if (.not.DoCholesky) required=required+real(nBas,wp)**4
  call mma_maxDBLE(avail)
  if ((required>0.8_wp*real(avail,wp)).or.(required>real(huge(avail),wp)/2.0_wp)) then
    Message='Insufficient memory for ROMP2 reference Fock construction.'
    return
  end if
  np=nBas*(nBas+1)/2
  if (.not.DoCholesky) then
    call f_Inquire('ORDINT',exists)
    if (.not.exists) then
      Message='Conventional ROMP2 requires stored ORDINT; use SEWARD NoCholesky.'
      return
    end if
    if (AuxTwo%Opn) return
    call mma_allocate(ERI,nBas,nBas,nBas,nBas,label='ROMP2 AO ERI')
  end if
  call mma_allocate(Buf,np+4,label='ROMP2 integral buffer')
  call mma_allocate(H,nBas,nBas,label='ROMP2 core Hamiltonian')
  call mma_allocate(S,nBas,nBas,label='ROMP2 overlap')
  call mma_allocate(DA,nBas,nBas,label='ROMP2 density alpha')
  call mma_allocate(DB,nBas,nBas,label='ROMP2 density beta')
  call mma_allocate(T,nBas,nBas,label='ROMP2 work')
  call mma_allocate(MA,nBas,nBas,label='ROMP2 MO Fock alpha')
  call mma_allocate(MB,nBas,nBas,label='ROMP2 MO Fock beta')
  if (.not.AuxOne%Opn) then
    call OpnOne(rc,0,'ONEINT',41)
    oneopened=.true.
    if (rc/=0) goto 900
  end if
  opt=ibset(ibset(0,sNoOri),sNoNuc)
  comp=1
  sym=1
  label='OneHam  '
  call RdOne(rc,opt,label,comp,Buf,sym)
  if (rc/=0) goto 900
  call Unpack(Buf,H)
  label='Mltpl  0'
  call RdOne(rc,opt,label,comp,Buf,sym)
  if (rc/=0) goto 900
  call Unpack(Buf,S)
  if (oneopened) then
    call ClsOne(rc,0)
    oneopened=.false.
    if (rc/=0) goto 900
  end if
  T=matmul(S,CAlpha)
  MA=matmul(transpose(CAlpha),T)
  do i=1,nBas
    MA(i,i)=MA(i,i)-1.0_wp
  end do
  if ((.not.all(ieee_is_finite(MA))).or.(maxval(abs(MA))>1.0e-7_wp)) then
    Message='ROHF orbitals are not orthonormal in the ONEINT overlap metric.'
    goto 900
  end if
  if (.not.DoCholesky) then
    lu=43
    call OpnOrd(rc,0,'ORDINT',lu)
    opened=.true.
    if (rc/=0) goto 900
    call GetOrd(rc,square,fs,fbasis,skip)
    if (rc/=0) goto 900
    if ((fs/=1).or.(fbasis(1)/=nBas).or.(skip(1)/=0)) goto 900
    opt=1
    do i=1,nBas
      do j=1,i
        call RdOrd_(rc,opt,1,1,1,1,Buf,np+1,nmat)
        opt=2
        if ((rc/=0).or.(nmat/=1)) goto 900
        q=0
        do k=1,nBas
          do l=1,k
            q=q+1
            ERI(i,j,k,l)=Buf(q)
            ERI(j,i,k,l)=Buf(q)
            ERI(i,j,l,k)=Buf(q)
            ERI(j,i,l,k)=Buf(q)
          end do
        end do
      end do
    end do
    call ClsOrd(rc)
    opened=.false.
    if (rc/=0) goto 900
  end if
  DA=matmul(CAlpha(:,1:nOccA),transpose(CAlpha(:,1:nOccA)))
  DB=matmul(CBeta(:,1:nOccB),transpose(CBeta(:,1:nOccB)))
  FA=H
  FB=H
  if (DoCholesky) then
    call ROMP2_Cho_Fock(DA,DB,FA,FB,rc,Message)
    if (rc/=0) goto 900
  else
    do j=1,nBas
      do i=1,nBas
        do l=1,nBas
          do k=1,nBas
            FA(i,j)=FA(i,j)+(DA(k,l)+DB(k,l))*ERI(i,j,k,l)-DA(k,l)*ERI(i,k,j,l)
            FB(i,j)=FB(i,j)+(DA(k,l)+DB(k,l))*ERI(i,j,k,l)-DB(k,l)*ERI(i,k,j,l)
          end do
        end do
      end do
    end do
  end if
  if ((.not.all(ieee_is_finite(FA))).or.(.not.all(ieee_is_finite(FB)))) goto 900
  call Get_dScalar('PotNuc',enuc)
  echeck=enuc+0.5_wp*(sum(DA*(H+FA))+sum(DB*(H+FB)))
  if ((.not.ieee_is_finite(echeck)).or.(abs(echeck-ESCF)>1.0e-7_wp)) then
    write(u6,*) 'ROHF stored/reconstructed energies: ',ESCF,echeck
    Message='ROHF determinant energy does not match RASSCF; check reference and Hamiltonian.'
    goto 900
  end if
  MA=matmul(transpose(CAlpha),matmul(FA,CAlpha))
  MB=matmul(transpose(CAlpha),matmul(FB,CAlpha))
  grad=0.0_wp
  do i=1,nBas
    oa=0.0_wp
    ob=0.0_wp
    if (i<=nOccA) oa=1.0_wp
    if (i<=nOccB) ob=1.0_wp
    do j=1,nBas
      pa=0.0_wp
      pb=0.0_wp
      if (j<=nOccA) pa=1.0_wp
      if (j<=nOccB) pb=1.0_wp
      grad=max(grad,abs((oa-pa)*MA(i,j)+(ob-pb)*MB(i,j)))
    end do
  end do
  if (grad>1.0e-5_wp) then
    write(u6,*) 'ROHF stationarity residual: ',grad
    Message='ROHF orbital stationarity check failed; tighten RASSCF convergence.'
    goto 900
  end if
  ierr=0
  Message=''
900 continue
  if (opened) call ClsOrd(rc2)
  if (oneopened) call ClsOne(rc2,0)
  call mma_deallocate(ERI,safe='*')
  call mma_deallocate(Buf)
  call mma_deallocate(H)
  call mma_deallocate(S)
  call mma_deallocate(DA)
  call mma_deallocate(DB)
  call mma_deallocate(T)
  call mma_deallocate(MA)
  call mma_deallocate(MB)
contains
  subroutine Unpack(P,A)
    real(kind=wp),intent(in) :: P(:)
    real(kind=wp),intent(out) :: A(nBas,nBas)
    integer(kind=iwp) :: x,y,z
    z=0
    do x=1,nBas
      do y=1,x
        z=z+1
        A(x,y)=P(z)
        A(y,x)=P(z)
      end do
    end do
  end subroutine Unpack
end subroutine ROMP2_Fock
end module ROMP2_Reference
