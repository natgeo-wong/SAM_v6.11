! This subroutine solves for perturbation omega based on WTG approximation as
! described in the appendix of Blossey, Bretherton & Wyant (2009):
!
!   http://dx.doi.org/10.3894/JAMES.2009.1.8
!
!  - original version: Peter Blossey (pblossey at gmail dot com), 
!                      September 2009 based on earlier implementaion in SAM.
!
! Note 1: this version handles the situation where the model top is
!   below the tropopause gracefully.  The boundary condition at the
!   tropopause (omega'=0) is enforced by adding additional levels
!   between the model top and the tropopause (whose height is assumed
!   to be 100 hPa if it is above the model top).
!
! Note 2: Below the tropopause, the temperature will be modified
!   by the WTG method through large-scale vertical motion.
!   Above the tropopause, the temperature should be nudged 
!   to the observed profile on a longish (~2 day?) timescale.
!
! Note 3: The wrapper routine allows the code to handle models
!   that index their pressure, temperature and moisture soundings
!   from either the surface-upwards or the top-downwards.  The driver
!   routine assumes that these soundings are indexed from the surface
!   upwards as in SAM, the model for which this routine was first
!   written.

subroutine wtg_jas2008(tabs_ref, ktrop)
implicit none

use grid, only: masterproc, nz, nzm, z
use vars, only: tabs0, rho
use vars, only: wtgtmp_o1d, wtgtmp_o2d, wtgtmp_v2I, wtgtmp_lhs, wtgtmp_rhs
use vars, only: o_wtg
use params, only: fcor, ggr,
use params, only: dowtg_timedependence, lambda_wtg, am_wtg, am_wtg_exp
use wtg_matrices

! ======= inputs =======
real, intent(in) :: tabs_ref(nzm) ! reference temperature sounding in K

! ======= output =======
integer, intent(out) :: ktrop ! index of interface just above the cold point.

! ======= local variables =======
integer :: k
real :: min_temp ! temporary variable used to find cold point of model sounding.
real :: coef_wtg ! coefficient on RHS of omega equation.

! ===== find index of cold point tropopause in vertical. =====
! reverse pressure coordinate, and find index
!   of cold point tropopause in the vertical.
ktrop = nzm+1 ! default is top of model/atmosphere (counting from surface)
min_temp = tabs0(nzm) 
do k = 1,nzm
  if(tabs0(k).lt.min_temp) then
    ktrop = k
    min_temp = tabs0(k)
  end if
end do
ktrop = ktrop + 1

! ========= Set up 1st-order, 2nd-order matrices ==========

order1(nzm, z, ktrop, wtgtmp_o1d)
order1(nzm, z, ktrop, wtgtmp_o2d)
order1(nzm, z, ktrop, wtgtmp_o2d)

! ========= Set up RHS ==========
tmp_rhs(:) = 0.
coef_wtg = (ggr*0.5*pi/lambda_wtg)**2 ! useful coefficient
do k = 2, ktrop-1
  ! RHS = k^2 * rho * g^2 * T' / T
  ! LHS will solve for omega, this is why rhs is g^2 and has no negative sign
  tmp_rhs(k) = coef_wtg * rho_model(k-1) &
               * (tabs_model(k-1) - tabs_ref(k-1)) / tabs_model(k-1)
end do

if(dowtg_timedependence)
  tmp_rhs = tmp_rhs + matmul(A,o_wtg) / dtn
  tmp_rhs(1) = 0
  do k = ktrop,nzm
    ! Setting all the right-hand side above the tropopause to zero
    tmp_rhs(k) = 0
  end do
end if

! Interpolate omega from the interfaces onto the model levels
!   from the surface up to the tropopause or top of model, 
!   whichever is lower.
omega_wtg(1:nzm) = &
      0.5*(tmp_rhs(1:nzm)+tmp_rhs(2:nzm+1))

if(ktrop.lt.nzm+1) then
  omega_wtg(ktrop-1:nzm) = 0.  ! set omega to zero above tropopause.
end if

end subroutine wtg_jas2008