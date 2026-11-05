! Read-only workbench diagnostics in a companion build of the untouched original
! sources. GNU ld intercepts WRITE_SHEAR_ACTIVITY from the driver, calls the
! original output routine, then exports state. No assignments to VPSC variables.
subroutine stock_report(step) bind(C,name='__wrap_write_shear_activity_')
  use iso_c_binding, only: c_int
  use vpsc8dim, only: time, sbarc, wgt, ngr, nmodes, nslmod, nsyst, nsm, ag, gamdot, gavmod, ihardlaw, nph
  use disden_rev, only: rho_f, rho_r
  use disden, only: rho_s, rho_deb
  implicit none
  integer(c_int), intent(in) :: step
  interface
    subroutine original_write(i) bind(C,name='__real_write_shear_activity_')
      import c_int
      integer(c_int), intent(in) :: i
    end subroutine
  end interface
  integer, save :: h,sys,index_file,rho_file,count,next=1
  real(8), save :: targets(16),rate
  integer :: u,i,j,g,k,offset
  real(8) :: strain, shares(nsyst(1)),den,rm(nslmod(1)),debris
  real(8) :: rot(3,3),a,b,c
  character(32) :: filename
  call original_write(step)
  if (step==1) then
    if (nph/=1) error stop 'The workbench observer requires a single native phase'
    open(newunit=u,file='stock_snapshots.in',status='old')
    read(u,*) rate,count
    read(u,*) targets(:count)
    close(u)
    open(newunit=h,file='MG_HISTORY.csv',status='replace')
    write(h,'(a)',advance='no') 'axial_strain,axial_stress_MPa'
    call columns(h,'family_',nmodes(1))
    open(newunit=sys,file='SYSTEM_ACTIVITY.csv',status='replace')
    write(sys,'(a)',advance='no') 'axial_strain'
    call columns(sys,'system_',nsyst(1))
    open(newunit=index_file,file='MG_SNAPSHOTS.csv',status='replace')
    write(index_file,'(a)') 'target_strain,actual_axial_strain,filename'
    if (ihardlaw==20.or.ihardlaw==23) then
      open(newunit=rho_file,file='STOCK_DENSITY.csv',status='replace')
      write(rho_file,'(a)',advance='no') 'axial_strain'
      do k=1,nslmod(1)
        write(rho_file,'(a,i0)',advance='no') ',rho_',k
      end do
      write(rho_file,'(a)') ',rho_debris,rho_total'
    end if
  end if
  strain=time*rate
  shares=0
  do g=1,ngr(1)
    shares=shares+wgt(g)*abs(gamdot(:nsyst(1),g))
  end do
  den=sum(shares)
  if (den>0) shares=shares/den
  write(h,'(*(es25.16e3,:,","))') strain,-sbarc(3,3),gavmod(:nmodes(1),1)
  write(sys,'(*(es25.16e3,:,","))') strain,shares
  if (ihardlaw==20.or.ihardlaw==23) then
    rm=0; debris=0
    do g=1,ngr(1)
      offset=0
      do k=1,nslmod(1)
        do i=offset+1,offset+nsm(k,1)
          if (ihardlaw==23) then
            rm(k)=rm(k)+wgt(g)*(rho_f(i,g)+rho_r(2*i-1,g)+rho_r(2*i,g))
          else
            rm(k)=rm(k)+wgt(g)*rho_s(i,g)
          end if
        end do
        offset=offset+nsm(k,1)
      end do
      if (ihardlaw==20) debris=debris+wgt(g)*rho_deb(g)
    end do
    write(rho_file,'(*(es25.16e3,:,","))') strain,rm,debris,sum(rm)+debris
  end if
  do while(next<=count)
    if (strain+1d-9<targets(next)) exit
    write(filename,'("MG_STATE_",i3.3,".csv")') next
    open(newunit=u,file=filename,status='replace')
    write(u,'(a)') 'id,phi1,Phi,phi2,volume_fraction'
    do g=1,ngr(1)
      rot=transpose(ag(:,:,g))
      call euler(1,a,b,c,rot)
      write(u,'(i0,4(",",es25.16e3))') g,a,b,c,wgt(g)
    end do
    close(u)
    write(index_file,'(2(es25.16e3,","),a)') targets(next),strain,trim(filename)
    next=next+1
  end do
contains
  subroutine columns(unit,prefix,n)
    integer, intent(in) :: unit,n
    character(*), intent(in) :: prefix
    integer :: j
    do j=1,n
      write(unit,'(a,i0)',advance='no') ','//prefix,j
    end do
    write(unit,'()')
  end subroutine
end subroutine
