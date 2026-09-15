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

subroutine UMP2_Energy(La,Lb,Eoa,Eva,Eob,Evb, &
                       nOA,nVA,nOB,nVB,nCho,DenTol, &
                       Eaa,Ebb,Eab,ierr)

use Definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: nOA,nVA,nOB,nVB,nCho

real(kind=wp), intent(in) :: La(nVA,nOA,nCho)
real(kind=wp), intent(in) :: Lb(nVB,nOB,nCho)
real(kind=wp), intent(in) :: Eoa(nOA),Eva(nVA)
real(kind=wp), intent(in) :: Eob(nOB),Evb(nVB)
real(kind=wp), intent(in) :: DenTol

real(kind=wp), intent(out) :: Eaa,Ebb,Eab
integer(kind=iwp), intent(out) :: ierr

integer(kind=iwp) :: i,j,a,b
real(kind=wp) :: Coul,Den

Eaa = 0.0_wp
Ebb = 0.0_wp
Eab = 0.0_wp
ierr = 0

! A positive denominator tolerance must be supplied.
if (DenTol <= 0.0_wp) then
  ierr = 1
  return
end if

call SameSpin(La,Eoa,Eva,nOA,nVA,Eaa)
if (ierr /= 0) then
  call ClearEnergies()
  return
end if

call SameSpin(Lb,Eob,Evb,nOB,nVB,Ebb)
if (ierr /= 0) then
  call ClearEnergies()
  return
end if

! Opposite-spin contribution: no exchange term.
do i=1,nOA
  do j=1,nOB
    do a=1,nVA
      do b=1,nVB

        Den = Eoa(i)+Eob(j)-Eva(a)-Evb(b)

        if (abs(Den) <= DenTol) then
          ierr = 2
          call ClearEnergies()
          return
        end if

        Coul = dot_product(La(a,i,:),Lb(b,j,:))
        Eab = Eab+Coul*Coul/Den

      end do
    end do
  end do
end do

contains

subroutine SameSpin(L,Eo,Ev,nO,nV,Ess)

  integer(kind=iwp), intent(in) :: nO,nV
  real(kind=wp), intent(in) :: L(nV,nO,nCho)
  real(kind=wp), intent(in) :: Eo(nO),Ev(nV)
  real(kind=wp), intent(out) :: Ess

  integer(kind=iwp) :: ii,jj,aa,bb
  real(kind=wp) :: Jint,Kint,Delta,Anti

  Ess = 0.0_wp

  ! No same-spin double excitations in these cases.
  if ((nO < 2) .or. (nV < 2)) return

  do ii=1,nO
    do jj=1,nO
      if (ii == jj) cycle

      do aa=1,nV
        do bb=1,nV
          if (aa == bb) cycle

          Delta = Eo(ii)+Eo(jj)-Ev(aa)-Ev(bb)

          if (abs(Delta) <= DenTol) then
            ierr = 2
            return
          end if

          Jint = dot_product(L(aa,ii,:),L(bb,jj,:))
          Kint = dot_product(L(bb,ii,:),L(aa,jj,:))
          Anti = Jint-Kint

          ! Full ordered sums require the factor 1/4.
          Ess = Ess+0.25_wp*Anti*Anti/Delta

        end do
      end do
    end do
  end do

end subroutine SameSpin

subroutine ClearEnergies()

  Eaa = 0.0_wp
  Ebb = 0.0_wp
  Eab = 0.0_wp

end subroutine ClearEnergies

end subroutine UMP2_Energy