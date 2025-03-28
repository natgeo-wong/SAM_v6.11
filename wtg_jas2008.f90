subroutine wtg_jas2008

use vars
use params

implicit none

   integer :: k
   integer :: kk
   integer :: ktrop
   integer :: ktrop1
   real :: ztrop
   real :: ztrop1
   real :: min_temp
   real :: coef, buffer(nzm), buffer1(nzm)

   if (dowtg_timedependence) then

      ! ! ===== find index of cold point tropopause in vertical. =====
      ! ! reverse pressure coordinate, and find index
      ! !   of cold point tropopause in the vertical.
      ! ktrop = nzm+1 ! default is top of model/atmosphere (counting from surface)
      ! min_temp = tg0(nzm)
      ! do k = 1,nzm
      !    if(tg0(k).lt.min_temp) then
      !       ktrop = k
      !       min_temp = tg0(k)
      !    end if
      ! end do

      ! ! ===== find index of previous cold point tropopause in vertical. =====
      ! ! reverse pressure coordinate, and find index
      ! !   of cold point tropopause in the vertical.
      ! ktrop1 = nzm !
      ! do k = nzm,1,-1
      !    if(w_wtg(k).EQ.0) then
      !       ktrop1 = k
      !    end if
      ! end do

      ! if ((ktrop1.NE.1).AND.(ktrop1.NE.ktrop))

      !    ztrop = z(ktrop)
      !    ztrop1 = z(ktrop1)
      !    wwtgi = 0.

      !    do k = 2 : (ktrop-1)
      !       do kk = 1 : nzm
      !          wwtgi(k)
      !       end do
      !    end do

      !    do k = 1,nzm
      !       w_wtg(k) = wwtgi(k)
      !    end do

      ! end if

      ! ! At tropopause and above, quash any pre-existing vertical velocity from previous timestep to zero.
      ! do k = ktrop,nzm
      !    w_wtg(k) = 0
      ! end do

      ktrop = nzm-4
   
   else

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

   end if

   tv_lsbg = tg0   * (1. + 0.61*qg0)

   if(dompiensemble) then
      coef = 1. / nsubdomains
      do k = 1, nzm
         buffer(k) = tabs0(k) * (1. + 0.61*qv0(k) - qn0(k) - qp0(k))
      end do
      call task_sum_real8(buffer,buffer1,nzm)
      do k = 1, nzm
         tv_wave(k) = buffer1(k) * coef
      end do
   else
      tv_wave = tabs0 * (1. + 0.61*qv0 - qn0 - qp0)
   end if

   dwwtgdt = 0.

   call calc_wtend(0.5*pi/lambda_wtg, w_wtg(1:ktrop), dwwtgdt(1:ktrop), &
                     tabs0(1:ktrop), tg0(1:ktrop), tv_wave(1:ktrop), tv_lsbg(1:ktrop), &
                     rho(1:ktrop), z(1:ktrop), zi(1:ktrop+1), ktrop)

   if (dowtg_timedependence) then

      w_wtg(1:nzm) = (w_wtg(1:nzm) + dwwtgdt * dt) / (1. + dt * am_wtg_time)

   else

      w_wtg = 0.
      w_wtg(1:nzm) = dwwtgdt / am_wtg_time

   end if

   contains

   subroutine calc_wtend(wn, w_curr, wtend, &
                           ta_curr, tabg_curr, tv_curr, tv_fullbg, &
                           rho_full, z_full, z_half, nz)

      use params, only: dowtg_timedependence

   !     ------------------------------ input arguments ------------------------------

      integer, intent(in) :: nz  ! number of midpoint levels, the number of interface levels is nz+1

      real, dimension(nz), intent(in) ::                 &
      w_curr,            &       ! Cell center wave vertical velocity
      ta_curr,           &       ! Cell center pressure
      tabg_curr,         &       ! Cell center pressure
      tv_curr,           &       ! Cell center virtual temperature
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
      do k=1,nz-1
         aa(k)=dz(k+1)/(dz(k)+dz(k+1))
         bb(k)=-1.
         cc(k)=dz(k)/(dz(k)+dz(k+1))
      end do

      ! Symmetric boundary conditions?
      aa(1)=0.
      bb(1)=-(2*dz(2)+dz(1))/(dz(1)+dz(2))

      if (dowtg_timedependence) then

         !radiating upper BC
         aa(nz)=dz(nz+1)/dz(nz)
         bb(nz)=-1.*dz(nz+1)/dz(nz)
         rhs(nz)=rhs(nz)+rho_full(nz)*sqrt(N2top)*wn*w_curr(nz)*dz(nz+1)

      else

         ! dirichlet boundary conditions at tropopause (w' = 0)
         aa(nz)=0.
         bb(nz)=1.
         cc(nz)=0.
         rhs(nz)=0

      end if

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

end subroutine wtg_jas2008
