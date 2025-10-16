module simple_ocean

!------------------------------------------------------------
! Purpose:
!
! A collection of routines used to specify fixed 
! or compute interactive SSTs, like slab-ocean model, etc.
!
! Author: Marat Khairoutdinov
! Based on dynamic ocean impelemntation from the UW version of SAM.
!------------------------------------------------------------

use grid
use params, only: dodynamicocean, dosstislands
implicit none

public set_sst     ! set SST 
public sst_evolve ! evolve SST according to a model set by the ocean_type

CONTAINS


SUBROUTINE set_sst()

  use vars, only: sstxy,t00
  use params, only: tabs_s, delta_sst, ocean_type, &
                    nx_sinusoidalsst, dx_sinusoidalsst, ny_sinusoidalsst

! parameters of the sinusoidal SST destribution 
! along the X for Walker-type simulatons( ocean-type = 1):

  real(8) tmpx(nx), pii, lx, yy, ly
  integer i,j, it,jt

  select case (ocean_type)

    case(0) ! fixed constant SST

      sstxy = tabs_s - t00

    case(1) ! Sinusoidal distribution along the x-direction:

      lx = float(nx_gl)*dx/float(nx_sinusoidalsst)
      do i = 1,nx
          tmpx(i) = float(mod(rank,nsubdomains_x)*nx+i-1)*dx
      end do
      pii = atan2(0.d0,-1.d0)
      do j=1,ny
        do i=1,nx
          sstxy(i,j) = tabs_s-delta_sst*cos(2.*pii*(tmpx(i)/lx-dx_sinusoidalsst)) - t00
        end do
      end do
    
    case(2) ! Sinusoidal distribution along the y-direction:
      
      call task_rank_to_index(rank,it,jt)
      
      pii = atan2(0.d0,-1.d0)
      lx = float(ny_gl)*dy/float(ny_sinusoidalsst)
      do j=1,ny
        yy = dy*(j+jt-(ny_gl+YES3D-1)/2-1)
        do i=1,nx
          sstxy(i,j) = tabs_s+delta_sst*(2.*cos(pii*yy/lx)-1.) - t00
        end do
      end do 

    case(3) !  Aquaplanet Experiment (APE) "QOBS" (Neale and Hoskins 2001)

      call task_rank_to_index(rank,it,jt)

      pii = atan2(0.d0,-1.d0)
      lx = float(nx_gl)*dx
      ly = 0.5*float(ny_gl)*dy
      do j=1,ny
        yy = dy*(j+jt-(ny_gl+YES3D-1)/2-1)
        do i=1,nx
          sstxy(i,j) = tabs_s+delta_sst*(1-0.5*(sin(1.5*pii/180.*yy/40000000.*360.)**2 &
                                                +sin(1.5*pii/180.*yy/40000000.*360.)**4))-t00
        end do
      end do


    case default

      if(masterproc) then
          print*, 'unknown ocean type in set_sst. Exiting...'
          call task_abort
      end if

  end select

  if(dodynamicocean.and.dosstislands) call sst_islands_setmld()

end subroutine set_sst

SUBROUTINE sst_evolve
  use vars, only: sstxy, t00, fluxbt, fluxbq, rhow,qocean_xy
  use params, only: cp, lcond, tabs_s, ocean_type, dossthomo, dosstislands, &
                    depth_slab_ocean, Szero, deltaS, timesimpleocean
  use rad, only: swnsxy, lwnsxy

  real, parameter :: rhor = 1000. ! density of water (kg/m3)
  real, parameter :: cw = 4187.   ! Liquid Water heat capacity = 4187 J/kg/K
  real factor_cp, factor_lc, qoceanxy
  real tmpx(nx), lx, distsq
  real(8) sss(1),ssss(1)
  integer i, j, it, jt

    if(time.lt.timesimpleocean) return

    lx = float(nx_gl)*dx
    do i = 1,nx
      tmpx(i) = float(mod(rank,nsubdomains_x)*nx+i-1)*dx
    end do

    ! Define weight factors for the mixed layer heating due to
    ! the model's sensible and latent heat flux.
    factor_cp = rhow(1)*cp
    factor_lc = rhow(1)*lcond

    ! Use forward Euler to integrate the differential equation
    ! for the ocean mixed layer temperature: dT/dt = S - E.
    ! The source: CPT?GCSS WG4 idealized Walker-circulation 
    ! RCE Intercomparison proposed by C. Bretherton.
    if (.NOT.dosstislands) then
      do j=1,ny
        do i=1,nx
            qoceanxy       = Szero + deltaS*abs(2.*tmpx(i)/lx - 1)
            qocean_xy(i,j) = qocean_xy(i,j) + qoceanxy * dtfactor

            sstxy(i,j) = sstxy(i,j) &
                + dtn*(swnsxy(i,j)          & ! SW Radiative Heating
                - lwnsxy(i,j)               & ! LW Radiative Heating
                - factor_cp*fluxbt(i,j)     & ! Sensible Heat Flux
                - factor_lc*fluxbq(i,j)     & ! Latent Heat Flux
                + qoceanxy)            & ! Ocean Heating
                /(rhor*cw*depth_slab_ocean)        ! Convert W/m^2 Heating to K/s
        end do
      end do
    else
      ! Kuang Lab Addition
      ! Specify a small "island" where SST varies. SST elsewhere is held fixed.
      ! Added by Nathanael Wong on 2023/07/08
      call sst_islands()
    end if

    if(dossthomo) then
      sss = 0.
      do j=1,ny
        do i=1,nx
          sss(1) = sss(1) + sstxy(i,j)
        end do
      end do
      sss(1) = sss(1) / dble(nx*ny)
      if(dompi) then
          call task_sum_real8(sss,ssss,1)
          sss = ssss /float(nsubdomains)
      end if ! dompi
      if(ocean_type.eq.2) then
          tabs_s = sss(1) + t00
          call set_sst()
      else
          sstxy(:,:) = sss(1)
      end if
    end if

end subroutine sst_evolve

subroutine sst_islands_setmld

  ! Kuang Lab Addition
  ! Setup the mixed-layer depth for island archipelagoes
  ! Created submodule on 2023/07/10

  use vars, only: lsm_xy, mld_xy
  use params, only: sstislands_radius, sstislands_landmld, sstislands_oceanmld, &
                    sstislands_nrow, sstislands_ncol, sstislands_sep, &
                    readlsm, lsmfile

  real distsq, island_lon, island_lat
  integer i, j, it, jt, irow, icol
  integer lsm(1:nx_gl,1:ny_gl)

  call task_rank_to_index(rank,it,jt)

  if(readlsm) then

    if(masterproc) print*,'opening land-sea mask file: ',trim(lsmfile)
    open(11,file=trim(lsmfile),status='old',form='unformatted')
    if(masterproc) print*,'extracting land-sea mask data from file ...'
    read(11) lsm(1:nx_gl,1:ny_gl)
    if(masterproc) print*,'putting land-sea mask into respective task ...'
    lsm_xy(1:nx,1:ny) = lsm(1+it:nx+it,1+jt:ny+jt)
    close(11)

  else

    select case (sstisland_type)
      
      case(1) ! Round islands

        do icol = 1, sstislands_ncol
          do irow = 1, sstislands_nrow

            island_lon = nx_gl * dx / 2 + (icol-1 - (sstislands_ncol-1)*0.5) * sstislands_sep
            island_lat = ny_gl * dx / 2 + (irow-1 - (sstislands_nrow-1)*0.5) * sstislands_sep

            do j=1,ny
              do i=1,nx

                distsq = ((i+it-1)*dx - island_lon)**2 + ((j+jt-1)*dy - island_lat)**2
                if (distsq.lt.sstislands_radius**2) then
                  lsm_xy(i,j) = 1
                else
                  lsm_xy(i,j) = 0
                end if

              end do
            end do

          end do
        end do
      
      case(2) ! Islands channel

        island_lon = nx_gl * dx / 2
        do i=1,nx

          distsq = ((i+it-1)*dx - island_lon)**2
          if (distsq.lt.(sstislands_radius/2)**2) then
            lsm_xy(i,:) = 1
          else
            lsm_xy(i,:) = 0
          end if

        end do
      
      case default

        if(masterproc) then
            print*, 'unknown slab-island type in sst_islands_setmld. Exiting...'
            call task_abort
        end if

    end select

  end if

  do j = 1,ny
    do i = 1,nx
      if(lsm_xy(i,j).eq.1) then
        mld_xy(i,j) = sstislands_landmld
      else
        mld_xy(i,j) = sstislands_oceanmld
      end if
    end do
  end do

end subroutine sst_islands_setmld

subroutine sst_islands

  ! Kuang Lab Addition
  ! Specify an island archipelgao of islands where SST varies
  ! SST elsewhere held fixed, or varies according to MLD specified by sstislands_oceanmld
  ! Converted to submodule on 2023/07/09

  use vars, only: sstxy, fluxbt, fluxbq, rhow, qocean_xy, lsm_xy, mld_xy
  use params, only: cp, lcond, Szero, deltaS
  use rad, only: swnsxy, lwnsxy

  real, parameter :: rhor = 1000. ! density of water (kg/m3)
  real, parameter :: cw = 4187.   ! Liquid Water heat capacity = 4187 J/kg/K
  real factor_cp, factor_lc, qoceanxy
  real tmpx(nx), lx
  integer i, j


  ! Define weight factors for the mixed layer heating due to
  ! the model's sensible and latent heat flux.
  factor_cp = rhow(1)*cp
  factor_lc = rhow(1)*lcond

  do j=1,ny
    do i=1,nx

      if (.NOT.(mld_xy(i,j).EQ.0)) then
        if (lsm_xy(i,j).EQ.0) then
          qoceanxy       = Szero + deltaS*abs(2.*tmpx(i)/lx - 1)
          qocean_xy(i,j) = qocean_xy(i,j) + qoceanxy * dtfactor
        else
          qoceanxy       = 0
          qocean_xy(i,j) = 0
        end if
        sstxy(i,j) = sstxy(i,j) &
            + dtn * (swnsxy(i,j)        & ! SW Radiative Heating
            - lwnsxy(i,j)               & ! LW Radiative Heating
            - fluxbt(i,j) * factor_cp   & ! Sensible Heat Flux
            - fluxbq(i,j) * factor_lc   & ! Latent Heat Flux
            + qoceanxy)                 & ! Ocean Heating
            /(rhor*cw*mld_xy(i,j))        ! Convert W/m^2 Heating to K/s
      end if

    end do
  end do

end subroutine sst_islands

end module simple_ocean