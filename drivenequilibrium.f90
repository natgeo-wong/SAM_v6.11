
subroutine drivenequilibrium
	
use vars
use microphysics, only: micro_field, index_water_vapor, nmicro_fields, mklsadv
implicit none

integer i,j,k,k1,k2,n
real rdz, dq
real t_vtend, q_vtend
real t_tend(nx,ny,nzm), q_tend(nx,ny,nzm)

do k=2,nzm-1
  if(whadley(k).ge.0) then
     rdz=whadley(k)/(dz*adzw(k))	
     k1 = k
     k2 = k-1 
  else
     rdz=whadley(k)/(dz*adzw(k+1))       
     k1 = k+1
     k2 = k
  end if
  do j=1,ny
    do i=1,nx
      t_tend(i,j,k) =  - rdz * (tabs0(k1)-tabs0(k2)) - whadley * ggr / cp
      q_tend(i,j,k) =  - rdz * (qv0(k1)-qv0(k2))
    end do
  end do

end do
do k=2,nzm-1
  t_vtend = 0.
  q_vtend = 0.
  do j=1,ny
    do i=1,nx
      t(i,j,k) = t(i,j,k) + dtn * t_tend(i,j,k)
      micro_field(i,j,k,index_water_vapor) = max(0.,micro_field(i,j,k,index_water_vapor) &
                         + dtn * q_tend(i,j,k))
      t_vtend = t_vtend + t_tend(i,j,k)
      q_vtend = q_vtend + q_tend(i,j,k)
    end do
  end do
  t_vtend = t_vtend / float(nx*ny) - whadley(k) * ggr / cp
  q_vtend = q_vtend / float(nx*ny) 
  ttend(k) = ttend(k) + t_vtend
  qtend(k) = qtend(k) + q_vtend
  tlsvadv(k) = tlsvadv(k) + t_vtend
  qlsvadv(k) = qlsvadv(k) + q_vtend
end do

end
