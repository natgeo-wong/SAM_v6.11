
subroutine subsidence_3d()      

        use params, only: ggr, cp, dotqlsvadv, douvlsvadv
use vars
use microphysics, only: micro_field, index_water_vapor, nmicro_fields, mklsadv
implicit none

integer i,j,k,k1,k2,n
real rdz, dq
real t_vtend, q_vtend
real t_tend(nx,ny,nzm), q_tend(nx,ny,nzm)

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
   do j=1,ny
      do i=1,nx
         if (douvlsvadv) then
            dudt(i,j,k,na) = dudt(i,j,k,na) - rdz*(u(i,j,k1)-u(i,j,k2))
            dvdt(i,j,k,na) = dvdt(i,j,k,na) - rdz*(v(i,j,k1)-v(i,j,k2))
            ! for calculating subdomain mean adv tendency for diag output
            ulsvadv(k) = ulsvadv(k) - rdz*(u(i,j,k1)-u(i,j,k2))
            vlsvadv(k) = vlsvadv(k) - rdz*(v(i,j,k1)-v(i,j,k2))
         end if
         
         if (dotqlsvadv) then
            t_tend(i,j,k) =  - rdz * (t(i,j,k1)-t(i,j,k2))
            q_tend(i,j,k) =  &
               - rdz * (micro_field(i,j,k1,index_water_vapor)-micro_field(i,j,k2,index_water_vapor))
            t(i,j,k) = t(i,j,k) + dtn * t_tend(i,j,k)
            micro_field(i,j,k,index_water_vapor) = max(0.,micro_field(i,j,k,index_water_vapor) &
                          + dtn * q_tend(i,j,k))
            ! for calculating subdomain mean adv tendency for diag output
            t_vtend = t_vtend + t_tend(i,j,k)
            q_vtend = q_vtend + q_tend(i,j,k)
            ! Apply large-scale vertical advection to all microphysics fields,
            !   not just water vapor or total water.  This resolves some issues
            !    when index_water_vapor refers to something other than total
            !    water (i.e., vapor+cloud).
            do n = 1,nmicro_fields
               if(n.ne.index_water_vapor) then
                  dq = - rdz * (micro_field(i,j,k1,n)-micro_field(i,j,k2,n))
                  micro_field(i,j,k,n) = micro_field(i,j,k,n) + dtn*dq
                  mklsadv(k,n) = mklsadv(k,n) + dq
               end if
            end do
         end if
      end do
   end do
   ! normalize tendencies
   mklsadv(k,index_water_vapor) = mklsadv(k,index_water_vapor) + q_vtend
   ulsvadv(k) = ulsvadv(k) / float(nx*ny)
   vlsvadv(k) = vlsvadv(k) / float(nx*ny)
   t_vtend = t_vtend / float(nx*ny)
   q_vtend = q_vtend / float(nx*ny)
   ttend(k) = ttend(k) + t_vtend
   qtend(k) = qtend(k) + q_vtend
   tlsvadv(k) = tlsvadv(k) + t_vtend
   qlsvadv(k) = qlsvadv(k) + q_vtend

end do

end
