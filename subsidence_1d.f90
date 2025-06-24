
subroutine subsidence_1d()      
! in this case, domain mean U,V should be zero, so douvlsvadv should be effectively false 
use params, only: ggr, cp, dotqlsvadv, doadvinic, doadvbg, doadvensnoise
use vars
use microphysics, only: micro_field, index_water_vapor, nmicro_fields, mklsadv
implicit none

integer i,j,k,k1,k2
real rdz
real t_vtend, q_vtend

do k=1,nzm-1
   if(k.eq.1) then
      rdz=wsub(1)/(dz*adzw(2))
      k1 = 2
      k2 = 1
   else
      if(wsub(k).ge.0) then
         rdz=wsub(k)/(dz*adzw(k))
         k1 = k
         k2 = k-1 
      else
         rdz=wsub(k)/(dz*adzw(k+1))       
         k1 = k+1
         k2 = k
      end if
   end if
   t_vtend = 0.
   q_vtend = 0.
   if (dotqlsvadv) then
      if (doadvinic) then
         ! advect the initial profile as background
         t_vtend = - rdz * (tg0(k1)-tg0(k2)) - wsub(k) * ggr / cp
         q_vtend = - rdz * (qg0(k1)-qg0(k2))
      else if (doadvbg) then
         ! advect the wtg background - to be completed, placeholder now
         t_vtend = - rdz * (t_wtg(k1)-t_wtg(k2)) - wsub(k) * ggr / cp
         q_vtend = - rdz * (q_wtg(k1)-q_wtg(k2))
      else if (doadvensnoise) then
         ! advect the mean profile of each ensemble member
         t_vtend = - rdz * (tabs0(k1)-tabs0(k2)) - wsub(k) * ggr / cp
         q_vtend = - rdz * (qv0(k1)-qv0(k2))
      else
         ! advect ensemble mean profile
         t_vtend = - rdz * (t_wtg(k1)-t_wtg(k2)) - wsub(k) * ggr / cp
         q_vtend = - rdz * (q_wtg(k1)-q_wtg(k2))
      end if
      do j=1,ny
         do i=1,nx
            t(i,j,k) = t(i,j,k) + dtn * t_vtend
            micro_field(i,j,k,index_water_vapor) = max(0.,micro_field(i,j,k,index_water_vapor) & 
                    + dtn * q_vtend)
            ! currently not including the advection of condensates
         end do
      end do
   end if
   ! normalize tendencies
   mklsadv(k,index_water_vapor) = mklsadv(k,index_water_vapor) + q_vtend*float(nx*ny)
   ttend(k) = ttend(k) + t_vtend
   qtend(k) = qtend(k) + q_vtend
   tlsvadv(k) = tlsvadv(k) + t_vtend
   qlsvadv(k) = qlsvadv(k) + q_vtend

end do

end
