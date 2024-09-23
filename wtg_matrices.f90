module wtg_matrices

!--------------------------------------------------------------
! Purpose:
!
! A collection of routines used to invert first- and second-
! order derivatives for the WTG schemes. Since these inversions
! are reusable, best to pull them all out into a single file
! 
! Author: Nathanael Wong
! Date: 2024 September
!--------------------------------------------------------------

contains

subroutine order1(nzm, z, ktrop, A)

  implicit none

  ! ======= inputs =======
  integer, intent(in) :: nzm ! number of model levels
  integer, intent(in) :: ktrop ! level of the tropopause
  real, intent(in) :: z(nzm+1) ! vertical grid

  ! ======= output =======
  real, intent(out) :: A(nzm,nzm) ! WTG large-scale pressure velocity in Pa/s on model levels.

  ! ======= local variables =======
  real :: dz(nzm) ! vertical grid spacing

  dz(1) = z(1)
  do k = 2,n
    dz(k) = z(k) - z(k-1)
  end do

  A = 0.
  A(1,1) = 1.
  do k = 2,(ktrop-1)
    A(k,k-1) = -1 / (dz(k-1)+dz(k))
    A(k,k+1) = 1 / (dz(k-1)+dz(k))
  end do
  do k = ktrop,nzm
    A(k,k) = 1.
  end do

end subroutine order1

subroutine order2(nzm, z, ktrop, A)

  implicit none

  ! ======= inputs =======
  integer, intent(in) :: nzm ! number of model levels
  integer, intent(in) :: ktrop ! level of the tropopause
  real, intent(in) :: z(nzm+1) ! vertical grid

  ! ======= output =======
  real, intent(out) :: A(nzm,nzm) ! WTG large-scale pressure velocity in Pa/s on model levels.

  ! ======= local variables =======
  real :: dz(nzm) ! vertical grid spacing

  dz(1) = z(1)
  do k = 2, nzm
    dz(k) = z(k) - z(k-1)
  end do

  A = 0.
  A(1,1) = 1.
  do k = 2,(ktrop-1)
    A(k,k-1) =  2 / (dz(k-1) + dz(k)) / dz(k-1)
    A(k,k)   = -4 / (dz(k-1) * dz(k))
    A(k,k+1) =  2 / (dz(k-1) + dz(k)) / dz(k)
  end do
  do k = ktrop,nzm
    A(k,k) = 1.
  end do

end subroutine order2

subroutine vec2eye(nzm, vec, A)

  implicit none

  ! ======= inputs =======
  integer, intent(in) :: nzm ! number of levels
  real, intent(in) :: vec(nzm) ! vertical grid

  ! ======= output =======
  real, intent(out) :: A(nzm,nzm) ! WTG large-scale pressure velocity in Pa/s on model levels.

  A = 0.
  do k = 1,nzm
    A(k,k) = vec(k)
  end do

end subroutine vec2eye

end module wtg_matrices