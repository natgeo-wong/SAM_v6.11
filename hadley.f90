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

subroutine hadley(masterproc, nzm, nz, z, tabs_model, wmax, w1, w2, w3, w4, w5, whadley)

implicit none

! ======= inputs =======
logical, intent(in) :: masterproc ! For printing out to logs
integer, intent(in) :: nzm ! number of model levels
integer, intent(in) :: nz  ! number of vertical levels
real, intent(in) :: z(nz) ! pressure of model levels in Pa (domain-mean for LES)
real, intent(in) :: tabs_model(nzm) ! model temperature profile in K (domain-mean for LES)

real, intent(in) :: wmax ! pressure of model levels in Pa (domain-mean for LES)
real, intent(in) :: w1 ! pressure of model levels in Pa (domain-mean for LES)
real, intent(in) :: w2 ! pressure of model levels in Pa (domain-mean for LES)
real, intent(in) :: w3 ! pressure of model levels in Pa (domain-mean for LES)
real, intent(in) :: w4 ! pressure of model levels in Pa (domain-mean for LES)
real, intent(in) :: w5 ! pressure of model levels in Pa (domain-mean for LES)

! ======= output =======
real, intent(out) :: whadley(nzm) ! WTG large-scale pressure velocity in Pa/s on model levels.

! ======= local variables =======
integer :: k
integer :: ktrop
integer :: inum
! real :: min_temp ! temporary variable used to find cold point of model sounding.
real :: tabs_grad ! temporary variable used to find cold point of model sounding.
real :: ztrop ! Height of tropopause level (m)
real, parameter :: pi = 3.141592653589793 ! from MATLAB, format long.

if (z(nz) < 1.e4) then

  if(masterproc) then
    write(*,*) '******* ERROR in WTG Scheme of Raymond and Zeng (2005) *******'
    write(*,*) 'The model top is less than 10 km high. The WTG module of'
    write(*,*) 'Raymond and Zeng (2005) has not yet been configured to run on'
    write(*,*) 'a short domain.'
  end if
  call task_abort()

end if

! ! ===== find index of cold point tropopause in vertical. =====
! ! reverse pressure coordinate, and find index
! !   of cold point tropopause in the vertical.
! ktrop = nzm+1 ! default is top of model/atmosphere (counting from surface)
! min_temp = tabs_model(nzm) 
! do k = 1,nzm
!   if(tabs_model(k).lt.min_temp) then
!     ktrop = k
!     min_temp = tabs_model(k)
!   end if
! end do

! ===== find the tropopause using the 2 K/km method =====
! the lowest-temperature method fails, so trying to use the method where
! find the lowest level above 5 km where the gradient of the temperature
! falls below 2 K/km to determine the height of the troposphere
ktrop = nzm ! default is top of model/atmosphere (counting from surface)
tabs_grad = 0
do k = (nzm-1),2,-1
  tabs_grad = (tabs_model(k-1) - tabs_model(k+1)) / (z(k+1)-z(k-1)) * 1000
  if((tabs_grad.lt.2).AND.(z(k)>5000)) then
    ktrop = k
  end if
end do

! ===== prevents the Hadley Cell from growing too high =====
ztrop = z(ktrop)
if (ztrop>17500) then
  do k = 1,nzm
    if (z(k)<17500) then
      ktrop = k
    end if
  end do
end if
ztrop = z(ktrop)

do k = 1,ktrop
  whadley(k) = wmax * (sin(z(k) / ztrop * pi)     * w1 + &
                       sin(z(k) / ztrop * pi * 2) * w2 + &
                       sin(z(k) / ztrop * pi * 3) * w3 + &
                       sin(z(k) / ztrop * pi * 4) * w4 + &
                       sin(z(k) / ztrop * pi * 5) * w5)
end do

end subroutine hadley
