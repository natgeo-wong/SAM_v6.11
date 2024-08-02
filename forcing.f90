
subroutine forcing
	
use vars
use params
use microphysics, only: micro_field, index_water_vapor, total_water, mklsadv
use simple_ocean, only: sst_evolve

implicit none


integer i,j,k,n,nn,m,iz,iday0,iday
real coef, radtend, dayy
real tt(nzm,2),qq(nzm,2),uu(nzm,2),vv(nzm,2),ww(nzm,2),tp(nzm,2),pp(nzm,2)
real tpm(nzm)
real ratio1, ratio2, ratio_t1, ratio_t2
logical zgrid, pgrid

! linear response perturbation (Song Qiyu, 2022)
real, save :: delt_t, delt_q    ! Layer by layer perturbation

! ktrop index for tropopause
integer :: ktrop

call t_startf ('forcing')


! if doseasons=.false. do perpetual forcing

if(doseasons) then
   dayy = day
else
   iday0 = day0
   iday = day
   dayy = day-iday
   dayy = iday0 + dayy
end if


! ---------------------------------------------------------------
! Large-scale sounding:

nn=1
do i=1,nsnd-1
   if(day.gt.daysnd(i)) then
      nn=i
   endif
end do

do n=1,2

   m = nn+n-1
   zgrid = .false.
   pgrid = .false.
   if(zsnd(2,m).gt.zsnd(1,m)) zgrid=.true.
   if(psnd(2,m).lt.psnd(1,m)) pgrid=.true.

   if((.not.zgrid).and.(.not.pgrid)) then
      if(masterproc) print*,'error in grid in snd'
      stop
   end if

   do iz = 1,nzm
   if(zgrid) then
      do i = 2,nzsnd
         if(z(iz).le.zsnd(i,m)) then
            coef = (z(iz)-zsnd(i-1,m))/(zsnd(i,m)-zsnd(i-1,m)) 
            tt(iz,n)=tsnd(i-1,m)+(tsnd(i,m)-tsnd(i-1,m))*coef
            if(pgrid) then
               pp(iz,n)=psnd(i-1,m)+(psnd(i,m)-psnd(i-1,m))*coef
               tt(iz,n)=tt(iz,n)/((1000./pp(iz,n))**(rgas/cp))
            else
               tt(iz,n)=tt(iz,n)/prespotb(iz)
            endif
            tp(iz,n)=tsnd(i-1,m)+(tsnd(i,m)-tsnd(i-1,m))*coef
            qq(iz,n)=qsnd(i-1,m)+(qsnd(i,m)-qsnd(i-1,m))*coef
            uu(iz,n)=usnd(i-1,m)+(usnd(i,m)-usnd(i-1,m))*coef
            vv(iz,n)=vsnd(i-1,m)+(vsnd(i,m)-vsnd(i-1,m))*coef
            goto 11
         endif
      end do
   else
      do i = 2,nzsnd
         if(pres(iz).ge.psnd(i,m)) then
            coef = (pres(iz)-psnd(i-1,m))/(psnd(i,m)-psnd(i-1,m))
            tt(iz,n)=tsnd(i-1,m)+(tsnd(i,m)-tsnd(i-1,m))*coef/prespotb(iz)
            tp(iz,n)=tsnd(i-1,m)+(tsnd(i,m)-tsnd(i-1,m))*coef
            qq(iz,n)=qsnd(i-1,m)+(qsnd(i,m)-qsnd(i-1,m))*coef
            uu(iz,n)=usnd(i-1,m)+(usnd(i,m)-usnd(i-1,m))*coef
            vv(iz,n)=vsnd(i-1,m)+(vsnd(i,m)-vsnd(i-1,m))*coef
            pp(iz,n)=psnd(i-1,m)+(psnd(i,m)-psnd(i-1,m))*coef
            goto 11
         endif
      end do
   end if

   call atmosphere(z(iz-1)/1000.,ratio1,ratio2,ratio_t1)
   call atmosphere(z(iz)/1000.,ratio1,ratio2,ratio_t2)

   tt(iz,n)=ratio_t2/ratio_t1*tt(iz-1,n)
!      qq(iz,n)=max(0.,2.*qq(iz-1,n)-qq(iz-2,n))
   qq(iz,n) = qq(iz-1,n)*exp(-(z(iz)-z(iz-1))/3000.)
   uu(iz,n)=uu(iz-1,n)
   vv(iz,n)=vv(iz-1,n)

   11 continue

end do ! iz 

end do ! n

coef=(day-daysnd(nn))/(daysnd(nn+1)-daysnd(nn))

do k=1,nzm
   tg0(k)=tt(k,1)+(tt(k,2)-tt(k,1))*coef
   tp0(k)=tp(k,1)+(tp(k,2)-tp(k,1))*coef
   qg0(k)=qq(k,1)+(qq(k,2)-qq(k,1))*coef
   qg0(k)=qg0(k)*1.e-3
   ! Note that ug0 and vg0 maybe reset if dolargescale is true)
   ug0(k)=uu(k,1)+(uu(k,2)-uu(k,1))*coef - ug
   vg0(k)=vv(k,1)+(vv(k,2)-vv(k,1))*coef - vg
   pg0(k)=pp(k,1)+(pp(k,2)-pp(k,1))*coef - vg
end do

! ---------------------------------------------------------------
! Initialize tendencies:

ttend(:) = 0.
qtend(:) = 0.

! ---------------------------------------------------------------
! Large-Scale Advection Forcing:

if(dolargescale.and.time.gt.timelargescale) then

   nn=1
   do i=1,nlsf-1
      if(day.gt.dayls(i)) nn=i
   end do

   do n=1,2

      m = nn+n-1

      if(zls(2,m).gt.zls(1,m)) then
         zgrid=.true.
      else if(pls(2,m).lt.pls(1,m)) then
         zgrid=.false.
      else
         if(masterproc) print*,'error in grid in lsf'
         stop
      end if
      do iz = 1,nzm
         if(zgrid) then
            do i = 2,nzlsf
               if(z(iz).le.zls(i,m)) then
                  coef = (z(iz)-zls(i-1,m))/(zls(i,m)-zls(i-1,m))
                  tt(iz,n)=dtls(i-1,m)+(dtls(i,m)-dtls(i-1,m))*coef
                  qq(iz,n)=dqls(i-1,m)+(dqls(i,m)-dqls(i-1,m))*coef
                  uu(iz,n)=ugls(i-1,m)+(ugls(i,m)-ugls(i-1,m))*coef
                  vv(iz,n)=vgls(i-1,m)+(vgls(i,m)-vgls(i-1,m))*coef
                  ww(iz,n)=wgls(i-1,m)+(wgls(i,m)-wgls(i-1,m))*coef
                  goto 12
               endif
            end do
         else
            do i = 2,nzlsf
               if(pres(iz).ge.pls(i,m)) then
                  coef = (pres(iz)-pls(i-1,m))/(pls(i,m)-pls(i-1,m))
                  tt(iz,n)=dtls(i-1,m)+(dtls(i,m)-dtls(i-1,m))*coef
                  qq(iz,n)=dqls(i-1,m)+(dqls(i,m)-dqls(i-1,m))*coef
                  uu(iz,n)=ugls(i-1,m)+(ugls(i,m)-ugls(i-1,m))*coef
                  vv(iz,n)=vgls(i-1,m)+(vgls(i,m)-vgls(i-1,m))*coef
                  ww(iz,n)=wgls(i-1,m)+(wgls(i,m)-wgls(i-1,m))*coef
                  goto 12
               endif
            end do
         end if
         tt(iz,n)=0.
         qq(iz,n)=0.
         uu(iz,n)=uu(iz-1,n)
         vv(iz,n)=vv(iz-1,n)
         ww(iz,n)=0.
         12 continue

      end do

   end do ! n

   ! linear response perturbation: layer by layer (Song Qiyu, 2022)
   if(dolayerperturb) then
      delt_t = 0.5/86400.
      delt_q = 1.e-3*0.2/86400.
      ! Apply perturbation forcing
      if (tperturbi.gt.0) then
         tt(tperturbi,:) = tt(tperturbi,:)+tperturbA*delt_t
      end if
      if (qperturbi.gt.0) then
         ! For height with small humidity, rescale humidity perturbation
         delt_q = min(delt_q,0.2*qg0(qperturbi)/7200.)
         qq(qperturbi,:) = qq(qperturbi,:)+qperturbA*delt_q
      end if
   end if

   coef=(day-dayls(nn))/(dayls(nn+1)-dayls(nn))
   dosubsidence = .false.
   do k=1,nzm
      ttend(k)=tt(k,1)+(tt(k,2)-tt(k,1))*coef
      qtend(k)=qq(k,1)+(qq(k,2)-qq(k,1))*coef
      ug0(k)=uu(k,1)+(uu(k,2)-uu(k,1))*coef - ug
      vg0(k)=vv(k,1)+(vv(k,2)-vv(k,1))*coef - vg
      wsub(k)=ww(k,1)+(ww(k,2)-ww(k,1))*coef
      dosubsidence = dosubsidence .or. wsub(k).ne.0.
      do j=1,ny
         do i=1,nx
            t(i,j,k)=t(i,j,k)+ttend(k) * dtn
            micro_field(i,j,k,index_water_vapor) = &
               max(0.,micro_field(i,j,k,index_water_vapor) + qtend(k) * dtn)
         end do
      end do
   end do 

   pres0 = pres0ls(nn)+(pres0ls(nn+1)-pres0ls(nn))*coef

   if(wgls_holds_omega) then
      ! convert omega (sitting in wsub) into large-scale vertical velocity.
      ! Note that omega was read in from SCAM IOP netcdf input file.
      do k = 1,nzm
         wsub(k) = -wsub(k)/rho(k)/ggr
      end do
   end if

   !-------------------------------------------------------------------------------
   ! Kuang Lab Addition
   ! Save reference copy of large-scale vertical velocity before modification
   ! by WTG or scaling techniques, similar to Blossey's version of SAM
   wsub_ref(1:nzm) = wsub(1:nzm)

   if(dodgw) then

      if(wtgscale_time.gt.0) then
         twtgmax = (nstop * dt - timelargescale) * wtgscale_time
         twtg = time-timelargescale
         if(twtg.gt.twtgmax) then
         am_wtg_time = am_wtg
         else
         am_wtg_time = am_wtg * twtgmax / twtg
         endif
      else
         am_wtg_time = am_wtg
      endif

      if (dowtg_blossey_etal_JAMES2009) call wtg_james2009(nzm, &
            100.*pres, tg0, qg0, tabs0, qv0, qn0+qp0, &
            fcor, lambda_wtg, am_wtg_time, am_wtg_exp, o_wtg, ktrop)
      if (dowtg_decompdgw) then
         call wtg_james2009(nzm, &
            100.*pres, tg0, qg0, tabs0, qv0, qn0+qp0, &
            fcor, lambda_wtg, am_wtg_time, am_wtg_exp, owtgr, ktrop)
         call wtg_decompdgw(masterproc, &
            nzm, nz, z, 100.*pg0, tg0, qg0, tabs0, qv0, qn0+qp0, &
            lambda_wtg, am_wtg_time, wtgscale_vertmodenum, wtgscale_vertmodescl, &
            o_wtg, wwtgc, ktrop)
      end if

      ! convert from omega in Pa/s to wsub in m/s
      w_wtg(1:nzm) = -o_wtg(1:nzm)/rho(1:nzm)/ggr
      if (dowtg_decompdgw) wwtgr(1:nzm) = -owtgr(1:nzm)/rho(1:nzm)/ggr

   end if

   if (dotgr) then

      if(wtgscale_time.gt.0) then
         twtgmax = (nstop * dt - timelargescale) * wtgscale_time
         twtg = time-timelargescale
         if(twtg.gt.twtgmax) then
         tau_wtg_time = tau_wtg
         else
         tau_wtg_time = tau_wtg * twtg / twtgmax
         endif
      else
         tau_wtg_time = tau_wtg
      endif

      do k = 1,nzm
         tpm(k) = tabs0(k) * prespot(k)
      end do

      if (dowtg_raymondzeng_QJRMS2005)   call wtg_qjrms2005(masterproc, nzm, nz, z, &
                              tp0, tpm, tabs0, tau_wtg_time, dowtgLBL, boundstatic, &
                              dthetadz_min, w_wtg, wwtgr)
      if (dowtg_hermanraymond_JAMES2014) call wtg_james2014(masterproc, nzm, nz, z, &
                              tp0, tpm, tabs0, tau_wtg_time, dowtgLBL, boundstatic, &
                              dthetadz_min, wtgscale_vertmodepwr, w_wtg, wwtgr, wwtgc)
      if (dowtg_decomptgr)               call wtg_decomptgr(masterproc, nzm, nz, z, &
                              tp0, tpm, tabs0, tau_wtg_time, &
                              wtgscale_vertmodenum, wtgscale_vertmodescl, &
                              dowtgLBL, boundstatic, dthetadz_min, w_wtg, wwtgr, wwtgc)

      ! convert from omega in Pa/s to wsub in m/s
      o_wtg(1:nzm) = -w_wtg(1:nzm)*rho(1:nzm)*ggr
      owtgr(1:nzm) = -wwtgr(1:nzm)*rho(1:nzm)*ggr

   end if

   if (dotgr.OR.dodgw) then

      ! add to reference large-scale vertical velocity.
      wsub(1:nzm) = wsub(1:nzm) + w_wtg(1:nzm)
      dosubsidence = .true.

   end if

   if (dohadley) then
      if(hadscale_time.gt.0) then
         thadmax = (nstop * dt - timelargescale) * hadscale_time
         thad = time - timelargescale
         if(thad.gt.thadmax) then
            whad = whadmax
         else
            whad = whadmax * thad / thadmax
         endif
      else
         whad = whadmax
      endif
      call hadley(masterproc, nzm, nz, z, tabs0, whad, zhadmax, whadley)
      if(.NOT.dodrivenequilibrium) then
         wsub(1:nzm) = wsub(1:nzm) + whadley(1:nzm)
         dosubsidence = .true.
      end if
   end if

   ! ---------------------------------------------------------------
   ! Initialize large-scale advection tendencies:

   ulsvadv(:)   = 0.
   vlsvadv(:)   = 0.
   qlsvadv(:)   = 0.
   tlsvadv(:)   = 0.
   mklsadv(:,:) = 0. ! large-scale microphysical tendencies

   if(dosubsidence) call subsidence()
   if(dodrivenequilibrium) call drivenequilibrium()

   ! normalize large-scale vertical momentum forcing
   ulsvadv(:) = ulsvadv(:) / float(nx*ny) 
   vlsvadv(:) = vlsvadv(:) / float(nx*ny) 
   mklsadv(1:nzm,index_water_vapor) = qlsvadv(1:nzm) * float(nx*ny)

end if 

!---------------------------------------------------------------------
! Prescribed Radiation Forcing:


if(doradforcing.and.time.gt.timelargescale) then

   nn=1
   do i=1,nrfc-1
      if(day.gt.dayrfc(i)) nn=i
   end do

   do n=1,2

      m = nn+n-1 

      if(prfc(2,m).gt.prfc(1,m)) then
         zgrid=.true.
      else if(prfc(2,m).lt.prfc(1,m)) then
         zgrid=.false.
      else
         if(masterproc) print*,'error in grid in rad'
         stop
      end if
      do iz = 1,nzm
         if(zgrid) then
            do i = 2,nzrfc
               if(z(iz).le.prfc(i,m)) then
                  tt(iz,n)=dtrfc(i-1,m)+(dtrfc(i,m)-dtrfc(i-1,m))/(prfc(i,m)-prfc(i-1,m)) &
                                                                     *(z(iz)-prfc(i-1,m))
                  goto 13
               endif
            end do
         else
            do i = 2,nzrfc
               if(pres(iz).ge.prfc(i,m)) then
                  tt(iz,n)=dtrfc(i-1,m)+(dtrfc(i,m)-dtrfc(i-1,m))/(prfc(i,m)-prfc(i-1,m)) &
                                                                     *(pres(iz)-prfc(i-1,m))
                  goto 13
               endif
            end do
         end if
         tt(iz,n)=0.
         13 continue
      end do

   end do ! n

   coef=(day-dayrfc(nn))/(dayrfc(nn+1)-dayrfc(nn))
   do k=1,nzm
      radtend=tt(k,1)+(tt(k,2)-tt(k,1))*coef
      radqrlw(k)=radtend*float(nx*ny)
      radqrsw(k)=0.
      do j=1,ny
         do i=1,nx
            t(i,j,k)=t(i,j,k)+radtend*dtn 
         end do
      end do
   end do

endif


!----------------------------------------------------------------------------
! Surface flux forcing:

if(dosfcforcing.and.time.gt.timelargescale) then

   nn=1
   do i=1,nsfc-1
      if(day.gt.daysfc(i)) nn=i
   end do

   coef=(day-daysfc(nn))/(daysfc(nn+1)-daysfc(nn))
   tabs_s=sstsfc(nn)+(sstsfc(nn+1)-sstsfc(nn))*coef
   fluxt0=(shsfc(nn)+(shsfc(nn+1)-shsfc(nn))*coef)/(rhow(1)*cp)
   fluxq0=(lhsfc(nn)+(lhsfc(nn+1)-lhsfc(nn))*coef)/(rhow(1)*lcond)
   tau0=tausfc(nn)+(tausfc(nn+1)-tausfc(nn))*coef

   do j=1,ny
      do i=1,nx
         sstxy(i,j) = tabs_s - t00
      end do
   end do

   if(dostatis) then
      sstobs = tabs_s  ! sst is not averaged over the sampling period
      lhobs = lhobs + fluxq0 * rhow(1)*lcond
      shobs = shobs + fluxt0 * rhow(1)*cp
   end if

endif

!----------------------------------------------------------------------------
! Temperature Tendency Forcing:
! Simple Radiative Tendencies taken from Pauluis & Garner [2006]

if(doradtendency.and.time.gt.timelargescale) then

  do k = 1,nzm
    do j=1,ny
      do i=1,nx
        if (tabs(i,j,k)>207.5) then
            t(i,j,k) = t(i,j,k) - dtn * troptend / 86400
        else
            t(i,j,k) = t(i,j,k) + dtn * (200 - tabs(i,j,k)) / (5*86400)
        end if
      end do
    end do
  end do

endif

!-------------------------------------------------------------------------------

if(.not.dosfcforcing.and.dodynamicocean) call sst_evolve()

!-------------------------------------------------------------------------------
call t_stopf ('forcing')


end subroutine forcing
