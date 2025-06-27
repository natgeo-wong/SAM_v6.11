
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

! mpiensemble mean
real :: coef_subdomain
real :: buffer(nzm,5), buffer1(nzm,5)
real, save :: pres_wtg(nzm), pres_wtg_calcw(nzm)
real, save :: prespot_wtg(nzm), prespot_wtg_calcw(nzm)

! wtg background
real :: buffer2(nzm,3), buffer3(nzm,3)

real :: tmp(nzm)

real :: wtgtimestep

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
   
   ! compute wtg background, since this may be different from initial profile
   if (docalcwtgbg.and.icycle.eq.1) then
      if (nstep.gt.nstartwtg.and.nstep.le.nstartwtg+nstepwtgbg) then
         t_wtgbg = t_wtgbg + dble(tabs0)
         q_wtgbg = q_wtgbg + dble(qv0)
         tp_wtgbg = tp_wtgbg + dble(tabs0*prespot)
      end if
      if (nstep.eq.nstartwtg+nstepwtgbg) then
         t_wtgbg=t_wtgbg/dble(nstepwtgbg)
         q_wtgbg=q_wtgbg/dble(nstepwtgbg)
         tp_wtgbg=tp_wtgbg/dble(nstepwtgbg)
         coef_subdomain = 1. / dble(nsubdomains)
         do k=1, nzm
            buffer2(k,1) = t_wtgbg(k)
            buffer2(k,2) = q_wtgbg(k)
            buffer2(k,3) = tp_wtgbg(k)
         end do
         call task_sum_real8(buffer2,buffer3,nzm*3)
         do k=1, nzm
            t_wtgbg(k)=buffer3(k,1) * coef_subdomain
            q_wtgbg(k)=buffer3(k,2) * coef_subdomain
            tp_wtgbg(k)=buffer3(k,3) * coef_subdomain
         end do
      end if
   end if
   ! if the above calculation is not done, use the initial profile
   if (.not.docalcwtgbg) then
      t_wtgbg = tg0
      q_wtgbg = qg0
      tp_wtgbg = tp0
   end if

   ! for mpiensemble run, use the ensemble mean profiles to calculate forcing for all members
   ! for dowtgtimestep, also need the profiles to be updated every nstepwtg steps
   if (dompiensemble.or.dowtgtimestep) then
      if (icycle.eq.1 &
              .and. nstep.le.nstartwtg+nstepwtgbg*merge(1,0,docalcwtgbg) &
              .and. mod(nstep-1-floor((timelargescale-1e-5)/dt)-1, nstepwtg).eq.0) then
         ! calculate wtg profiles to be used in subsidence_1d/3d
         ! this is only necessary if wsub_ref.neq.0
         ! otherwise, subsidence is not called
         ! if not dowtgtimestep, nstepwtg is set as 1
         coef_subdomain = 1. / dble(nsubdomains)
         do k = 1, nzm
            buffer(k,1) = tabs0(k)
            buffer(k,2) = qv0(k)
            buffer(k,3) = qn0(k) + qp0(k)
            buffer(k,4) = pres(k)
            buffer(k,5) = prespot(k)
         end do
         call task_sum_real8(buffer,buffer1,nzm*5)
         do k = 1, nzm
            ! these are from the fields at the end of last step
            ! will be used in advection tendency
            t_wtg(k) = buffer1(k,1) * coef_subdomain
            q_wtg(k) = buffer1(k,2) * coef_subdomain
            qcond_wtg(k) = buffer1(k,3) * coef_subdomain
            pres_wtg(k) = buffer1(k,4) * coef_subdomain
            prespot_wtg(k) = buffer1(k,5) * coef_subdomain
         end do
      end if

      ! need to keep the forcing identical in each subdomain/ensemble member
      ! so sync the mean profiles before calculation, only when icycle=1
      ! update the profiles for forcing calculation every nstepwtg steps
      ! these calculated profiles are then kept constant for the next nstepwtg steps
      ! if not dowtgtimestep, nstepwtg is set as 1
      if (icycle.eq.1 &
              .and. nstep.gt.nstartwtg+nstepwtgbg*merge(1,0,docalcwtgbg) &
              .and. mod(nstep-1-nstartwtg-nstepwtgbg*merge(1,0,docalcwtgbg),nstepwtg).eq.0) then
         coef_subdomain = 1. / dble(nsubdomains)
         do k = 1, nzm
            buffer(k,1) = tabs0(k)
            buffer(k,2) = qv0(k)
            buffer(k,3) = qn0(k) + qp0(k)
            buffer(k,4) = pres(k)
            buffer(k,5) = prespot(k)
         end do
         call task_sum_real8(buffer,buffer1,nzm*5)
         do k = 1, nzm
            ! these are from the fields at the end of last step
            ! will be used in advection tendency
            t_wtg(k) = buffer1(k,1) * coef_subdomain
            q_wtg(k) = buffer1(k,2) * coef_subdomain
            qcond_wtg(k) = buffer1(k,3) * coef_subdomain
            pres_wtg(k) = buffer1(k,4) * coef_subdomain
            prespot_wtg(k) = buffer1(k,5) * coef_subdomain
         end do

         ! use the mean profile to calculate w
         t_wtg_calcw = t_wtg
         q_wtg_calcw = q_wtg
         qcond_wtg_calcw = qcond_wtg
         pres_wtg_calcw = pres_wtg
         prespot_wtg_calcw = prespot_wtg
         if (doadvensnoise) then
            ! use the profile from certain member to calculate w
            ! make sure the calculation below is identical for each subdomain/ensemble member
            if(masterproc)   tmp = tabs0
            call task_bcast_real(0,tmp,nzm)
            t_wtg_calcw = tmp
   
            if(masterproc)   tmp = qv0
            call task_bcast_real(0,tmp,nzm)
            q_wtg_calcw = tmp
   
            if(masterproc)   tmp = qn0 + qp0
            call task_bcast_real(0,tmp,nzm)
            qcond_wtg_calcw = tmp
   
            if(masterproc)   tmp = pres
            call task_bcast_real(0,tmp,nzm)
            pres_wtg_calcw = tmp
   
            if(masterproc)   tmp = prespot
            call task_bcast_real(0,tmp,nzm)
            prespot_wtg_calcw = tmp
         end if ! doadvensnoise
      end if ! icycle.eq.1 and time to update
   else
      ! not using mpiensemble or wtgtimestep
      t_wtg = tabs0
      q_wtg = qv0
      qcond_wtg = qn0 + qp0
      pres_wtg = pres
      prespot_wtg = prespot

      t_wtg_calcw = t_wtg
      q_wtg_calcw = q_wtg
      qcond_wtg_calcw = qcond_wtg
      pres_wtg_calcw = pres_wtg
      prespot_wtg_calcw = prespot_wtg
   end if
   
   ! -----------------------------------------------------------
   ! time to calculate/update the wtg large-scale w:
   ! nstep should be after nstartwtg (wait time) and nstepwtgbg (if necessary)
   ! if not dompiensemble, no need to control sync
   ! otherwise, only update w at the beginning of each block of nstepwtg steps
   if (nstep.gt.nstartwtg+nstepwtgbg*merge(1,0,docalcwtgbg) &
           .and. (.not.(dompiensemble.or.dowtgtimestep) &
              .or. (icycle.eq.1 &
                 .and. mod(nstep-1-nstartwtg-nstepwtgbg*merge(1,0,docalcwtgbg),nstepwtg).eq.0))) then

      ! calculate wtg large-scle vertical velocity
      if(dodgw) then
   
         if(wtgscale_time.gt.0) then
            twtgmax = (nstop * dt - max(timelargescale, (nstartwtg+nstepwtgbg*merge(1,0,docalcwtgbg))*dt)) * wtgscale_time
            twtg = time-day0*86400. - max(timelargescale, (nstartwtg+nstepwtgbg*merge(1,0,docalcwtgbg))*dt)
            if(twtg.gt.twtgmax) then
            am_wtg_time = am_wtg
            else
            am_wtg_time = am_wtg * twtgmax / twtg
            endif
         else
            am_wtg_time = am_wtg
         endif
   
         if (dowtg_blossey_etal_JAMES2009) then
   
            call wtg_james2009(nzm, &
               100.*pres_wtg, t_wtgbg, q_wtgbg, t_wtg_calcw, q_wtg_calcw, qcond_wtg_calcw, &
               fcor, lambda_wtg, am_wtg_time, am_wtg_exp, o_wtg, ktrop)
            w_wtg(1:nzm) = -o_wtg(1:nzm)/rho(1:nzm)/ggr
   
         end if
   
         if (dowtg_kuang_JAS2008) then
            wtgtimestep = dtn
            if (dowtgtimestep) wtgtimestep=dt*dble(nstepwtg)
            call wtg_jas2008(nzm, wtgtimestep, z, zi, rho, t_wtgbg, q_wtgbg, t_wtg_calcw, &
               q_wtg_calcw, qcond_wtg_calcw, lambda_wtg, am_wtg_time, w_wtg, dwwtgdt)
            o_wtg(1:nzm) = -w_wtg(1:nzm)*rho(1:nzm)*ggr
   
         end if
   
         if (dowtg_decompdgw) then
   
            call wtg_james2009(nzm, &
               100.*pres_wtg, t_wtgbg, q_wtgbg, t_wtg_calcw, q_wtg_calcw, qcond_wtg_calcw, &
               fcor, lambda_wtg, am_wtg_time, am_wtg_exp, owtgr, ktrop)
            call wtg_decompdgw(masterproc, &
               nzm, nz, z, 100.*pg0, t_wtgbg, q_wtgbg, t_wtg_calcw, q_wtg_calcw, qcond_wtg_calcw, &
               lambda_wtg, am_wtg_time, wtgscale_vertmodenum, wtgscale_vertmodescl, &
               o_wtg, wwtgc, ktrop)
   
            w_wtg(1:nzm) = -o_wtg(1:nzm)/rho(1:nzm)/ggr
            wwtgr(1:nzm) = -owtgr(1:nzm)/rho(1:nzm)/ggr
            
         end if
   
      end if
   
      if (dotgr) then
   
         if(wtgscale_time.gt.0) then
            twtgmax = (nstop * dt - max(timelargescale, (nstartwtg+nstepwtgbg*merge(1,0,docalcwtgbg))*dt)) * wtgscale_time
            twtg = time-day0*86400. - max(timelargescale, (nstartwtg+nstepwtgbg*merge(1,0,docalcwtgbg))*dt)
            if(twtg.gt.twtgmax) then
            tau_wtg_time = tau_wtg
            else
            tau_wtg_time = tau_wtg * twtg / twtgmax
            endif
         else
            tau_wtg_time = tau_wtg
         endif
   
         do k = 1,nzm
            tpm(k) = t_wtg_calcw(k) * prespot_wtg_calcw(k)
         end do
   
         if (dowtg_raymondzeng_QJRMS2005)   call wtg_qjrms2005(masterproc, nzm, nz, z, &
                                 tp_wtgbg, tpm, t_wtg_calcw, tau_wtg_time, dowtgLBL, boundstatic, &
                                 dthetadz_min, w_wtg, wwtgr)
         if (dowtg_hermanraymond_JAMES2014) call wtg_james2014(masterproc, nzm, nz, z, &
                                 tp_wtgbg, tpm, t_wtg_calcw, tau_wtg_time, dowtgLBL, boundstatic, &
                                 dthetadz_min, wtgscale_vertmodepwr, w_wtg, wwtgr, wwtgc)
         if (dowtg_decomptgr)               call wtg_decomptgr(masterproc, nzm, nz, z, &
                                 tp_wtgbg, tpm, t_wtg_calcw, tau_wtg_time, &
                                 wtgscale_vertmodenum, wtgscale_vertmodescl, &
                                 dowtgLBL, boundstatic, dthetadz_min, w_wtg, wwtgr, wwtgc)
   
         ! convert from omega in Pa/s to wsub in m/s
         o_wtg(1:nzm) = -w_wtg(1:nzm)*rho(1:nzm)*ggr
         owtgr(1:nzm) = -wwtgr(1:nzm)*rho(1:nzm)*ggr
   
      end if
   
      if (dohadley) then
         if(hadscale_time.gt.0) then
            thadmax = (nstop * dt - max(timelargescale, (nstartwtg+nstepwtgbg*merge(1,0,docalcwtgbg))*dt)) * hadscale_time
            thad = time-day0*86400. - max(timelargescale, (nstartwtg+nstepwtgbg*merge(1,0,docalcwtgbg))*dt)
            if(thad.gt.thadmax) then
               whad = whadmax
            else
               whad = whadmax * thad / thadmax
            endif
         else
            whad = whadmax
         endif
         call hadley(masterproc, nzm, nz, z, t_wtg_calcw, whad, zhadmax, whadley)
      end if
   end if

   ! ---------------------------------------------------------------
   ! Initialize large-scale advection tendencies:
   ! subsidence should be done regardless of w is updated or not
   ! the calculated tendencies should remain identical in each of the nstepwtg blocks
   if (nstep.gt.nstartwtg+nstepwtgbg*merge(1,0,docalcwtgbg)) then
      if (dotgr.OR.dodgw) then
         ! add to reference large-scale vertical velocity.
         wsub(1:nzm) = wsub_ref(1:nzm) + w_wtg(1:nzm)
         dosubsidence = .true.
      end if

      if (dohadley) then
         wsub(1:nzm) = wsub_ref(1:nzm) + whadley(1:nzm)
         dosubsidence = .true.
      end if
   end if

   ulsvadv(:)   = 0.
   vlsvadv(:)   = 0.
   qlsvadv(:)   = 0.
   tlsvadv(:)   = 0.
   mklsadv(:,:) = 0. ! large-scale microphysical tendencies

   ! calculate large-scale advection tendencies
   ! apply t/q tendencies to fields
   ! include u/v tendencies in dudt/dvdt
   if(dosubsidence) then
      if(doadv3d) then
         call subsidence_3d()
      else
         call subsidence_1d()
      end if
   end if
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
               if(pres_wtg(iz).ge.prfc(i,m)) then
                  tt(iz,n)=dtrfc(i-1,m)+(dtrfc(i,m)-dtrfc(i-1,m))/(prfc(i,m)-prfc(i-1,m)) &
                                                                 *(pres_wtg(iz)-prfc(i-1,m))
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
