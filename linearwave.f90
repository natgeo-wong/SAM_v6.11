module linearwave

use vars
use params

implicit none

contains

subroutine wtg_linearwave

   tv_bg   = tg0   * (1+0.61*qg0)
   tv_wave = tabs0 * (1+0.61*qv0)
   dwdt    = 0.

   call calc_wtend(
      0.5*pi/lambda_wtg,w_wtg(1:nzm-2),tv_wave(1:nzm-2),tv_bg(1:nzm-2),
      rho(1:nzm-2),z(1:nzm-2),zi(1:nzm-1),dwdt(1:nzm-2),nzm-2
   )

   w_wtg = w_wtg * (1. - dt * am_wtg_time) + dwdt * dt

end subroutine wtg_linearwave

subroutine calc_wtend(wn, w_curr, tv_curr, tv_fullbg, rho_full, & 
                      z_full, z_half, wtend, nz)
!     ------------------------------ input arguments ------------------------------

      integer, intent(in) :: nz  ! number of midpoint levels, the number of interface levels is nz+1

      real, dimension(nz), intent(in) ::                 &
      tv_curr,           &       ! Cell center virtual temperature
      w_curr,            &       ! Cell center wave vertical velocity
      tv_fullbg,         &       ! Cell center background virtual temperature
      z_full,            &       ! Cell center height
      rho_full                   ! Cell center density
      real, dimension(nz+1), intent(in) ::                 &
      z_half                     ! Interface-level height, including the surface, which is at height zero

      real, intent (in) :: wn    !the specified horizontal wavenumber

!     ------------------------------ output arguments -----------------------------

      real, dimension(nz), intent(out) ::                 &
      wtend                      ! Cell center w tendency

!     ------------------------------ local variables ------------------------------
      real  N2top                ! Top-level buoyancy frequency squared

      real, dimension(1:nz+1) ::                                             &
      dz                         ! grid spacing between midpoint levels
      real, dimension(1:nz) ::                                             &
      rhs,               &       ! The right hand side for the Gauss Elimination
      aa,                &       ! aa, bb, cc are the three columns of the tridiagonal matrix
      bb,                &
      cc

      integer :: k               ! Vertical loop counter

!     ------------------------------ executable code ------------------------------
     
      N2top=ggr/tv_fullbg(nz)*((tv_fullbg(nz)-tv_fullbg(nz-1))/(z_full(nz)-z_full(nz-1))+ggr/cp)

!     compute grid spacing between midpoint levels
      do k=2,nz
         dz(k)=z_full(k)-z_full(k-1)
      end do
      dz(1)=2*z_full(1)
      dz(nz+1)=2*(z_half(nz+1)-z_full(nz))

!     Gaussian Elimination
      rhs=0.
      do k=1,nz
         rhs(k)=-rho_full(k)*ggr*wn*wn*(tv_curr(k)-tv_fullbg(k))/tv_fullbg(k)*dz(k)*dz(k+1)/2.
      end do

!     set up the tridiagonal matrix
      do k=1,nz
         aa(k)=dz(k+1)/(dz(k)+dz(k+1))
         bb(k)=-1.
         cc(k)=dz(k)/(dz(k)+dz(k+1))
      end do

      !symmetric lower BC
      aa(1)=0.
      bb(1)=-(2*dz(2)+dz(1))/(dz(1)+dz(2))

      !radiating upper BC
      aa(nz)=dz(nz+1)/dz(nz)
      bb(nz)=-1.*dz(nz+1)/dz(nz)
      rhs(nz)=rhs(nz)+rho_full(nz)*sqrt(N2top)*wn*w_curr(nz)*dz(nz+1)

      !Gaussian Elimination with no pivoting
      do k=1,nz-1
         bb(k+1)=bb(k+1)-aa(k+1)/bb(k)*cc(k)
         rhs(k+1)=rhs(k+1)-aa(k+1)/bb(k)*rhs(k)
      end do

      !backward substitution
      rhs(nz)=rhs(nz)/bb(nz)
      do k= nz-1,1,-1
         rhs(k)=(rhs(k)-cc(k)*rhs(k+1))/bb(k)
      end do

      do k=1,nz
         wtend(k)=rhs(k)/rho_full(k)
      end do

end subroutine calc_wtend

end module linearwave
