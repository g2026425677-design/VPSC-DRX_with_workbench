! Zhou et al., Int. J. Plasticity 91 (2017) 48-76.
! DOI 10.1016/j.ijplas.2017.01.001. Optional IHARDLAW=40 (Cu).
! IHARDLAW=41: three-mode HCP extension of Eqs. 3-13; no extra strength law.
! SI units internally; CRSS handed to VPSC in MPa. Compile with
! -fdefault-real-8, as required by the original VPSC source.
module drx_state
  use iso_fortran_env, only: real64, int64
  use parameters, only: ngrmx,nmodmx,nsysmx
  implicit none
  real(real64) :: tau0, k1, k2, rhoc, nuc_c, mobility, temperature
  real(real64) :: burgers, chi, line_factor, rho_initial, diameter_initial
  real(real64) :: diameter_nucleus, q, dt_limit, strain_rate, mu, line_energy
  real(real64) :: box_radius_cubed, max_probability=0, volume_scale=1
  real(real64) :: density(ngrmx), radius(ngrmx), birth_strain(ngrmx)
  real(real64), allocatable :: mode_density(:,:), mode_shear(:,:)
  real(real64) :: mode_k1(nmodmx)=-1, mode_k2(nmodmx)=-1, mode_recovery_rate_m(nmodmx)=0
  real(real64) :: mode_rho_initial(nmodmx)=-1, mode_rho_nucleus(nmodmx)=-1
  logical :: has_nucleated(ngrmx)=.false.
  integer :: orientation_parent_id(ngrmx)=0, system_mode(nsysmx)=0
  integer :: generation(ngrmx), grain_id(ngrmx), parent_id(ngrmx)
  logical :: newborn(ngrmx), ready=.false., elastic_cached=.false.
  integer :: enabled=1, orientation_rule=0, snapshot_every=1000
  real(real64) :: orientation_spread=10
  integer :: created=0, removed=0, next_id=0, initial_count=0
  integer(int64) :: random_state=20261031_int64
  integer :: history_unit, cycles_unit
  integer :: mg_history_unit
  ! Law 41: HCP Mg extension. Law 40 retains the copper reproduction.
  integer :: material_kind=40, slip_exponent=10
  real(real64) :: mu_reference=16.6e9_real64, mu_slope=0.00053_real64
  real(real64) :: reference_temperature=300, reference_macro_rate=0.2_real64
  real(real64) :: strength_rate_m=0, recovery_rate_m=0
  real(real64) :: mobility_rate_m=0
  real(real64) :: mode_tau(nmodmx)=0, mode_b_ratio(nmodmx)=1
  real(real64) :: mode_chi_ratio(nmodmx)=1
  real(real64) :: system_tau(nsysmx)=0, system_b_ratio(nsysmx)=1
  real(real64) :: multiplicity(ngrmx)=1
  character(256) :: initial_size_file=''
  real(real64) :: rho_nucleus=-1
  real(real64) :: snapshot_strains(16)=-1
  integer :: next_snapshot=1, snapshot_index_unit
  real(real64) :: gos_shear_halfweight=-1
  real(real64) :: initial_strength_excess=0
  integer :: mg_observable_unit
  integer :: mg_cycle_unit
  integer :: mg_mode_unit, system_activity_unit
  integer :: births_by_generation(0:10)=0, consumed_by_generation(0:10)=0
contains
  subroutine numbered_columns(u,prefix,n)
    integer, intent(in) :: u,n
    character(*), intent(in) :: prefix
    integer :: i
    do i=1,n
      write(u,'(a,i0)',advance='no') ','//prefix,i
    end do
  end subroutine
  function uniform() result(u)
    ! Park-Miller using 64-bit intermediate: portable, deterministic stream.
    real(real64) :: u
    random_state=modulo(16807_int64*random_state,2147483647_int64)
    u=real(random_state,real64)/2147483647.0_real64
  end function

  function km_density(rho, shear, storage, recovery) result(r)
    ! Exact integration of Eq. 5 for constant k1,k2 over the shear increment.
    real(real64), intent(in) :: rho, shear
    real(real64), intent(in), optional :: storage, recovery
    real(real64) :: r, root, a, b
    if (shear==0) then
      r=rho
      return
    end if
    a=k1; b=k2
    if (present(storage)) a=storage
    if (present(recovery)) b=recovery
    if (b>0) then
      root=a/b+(sqrt(rho)-a/b)*exp(-0.5_real64*b*shear)
    else
      root=sqrt(rho)+0.5_real64*a*shear
    end if
    r=root*root
  end function

  function mode_resistance(im,g) result(tau)
    ! Eqs. 6-8: own-mode forest density, with debris and HP omitted as in paper.
    integer, intent(in) :: im,g
    real(real64) :: tau
    tau=mode_tau(im)+chi*mode_chi_ratio(im)*mu*burgers*mode_b_ratio(im)* &
      sqrt(mode_density(im,g))*1.e-6_real64
  end function
end module

subroutine drx_read(phase)
  use vpsc8dim
  use drx_state
  implicit none
  integer, intent(in) :: phase
  integer :: status, seed_value, im, is, last, nm
  character(256) :: config_file
  namelist /drx/ temperature,tau0,k1,k2,rhoc,nuc_c,mobility,burgers,chi, &
    line_factor,rho_initial,diameter_initial,diameter_nucleus,q,dt_limit, &
    strain_rate,seed_value,enabled,orientation_rule,orientation_spread,snapshot_every, &
    slip_exponent,mu_reference,mu_slope,reference_temperature,reference_macro_rate, &
    strength_rate_m,recovery_rate_m,mobility_rate_m,mode_tau,mode_b_ratio,mode_chi_ratio,initial_size_file, &
    rho_nucleus,snapshot_strains,gos_shear_halfweight,initial_strength_excess, &
    mode_k1,mode_k2,mode_recovery_rate_m,mode_rho_initial,mode_rho_nucleus
  if (phase/=1.or.nph/=1.or.nelem/=1) error stop 'DRX requires one phase and one element'
  ! SX record after IHARDLAW/IRATESENS specifies the namelist path.
  read(ur1,'(a)',iostat=status) config_file
  if (status/=0) error stop 'Missing DRX namelist filename in crystal file'
  ! All fitted parameters must be provided; no silent temperature interpolation.
  temperature=-1; tau0=-1; k1=-1; k2=-1; rhoc=-1; nuc_c=-1; mobility=-1
  burgers=2.56e-10_real64; chi=0.9_real64; line_factor=0.5_real64
  rho_initial=1.e9_real64; diameter_initial=78.e-6_real64
  diameter_nucleus=4.e-6_real64; q=0.9_real64; dt_limit=0.05_real64
  strain_rate=0.002_real64; seed_value=20261031
  material_kind=ihardlaw
  mode_tau=0; mode_b_ratio=1; initial_size_file=''
  open(newunit=status,file=trim(config_file),status='old',action='read')
  read(status,nml=drx)
  close(status)
  if (min(temperature,k1,rhoc,burgers,chi,rho_initial,diameter_initial, &
      diameter_nucleus,dt_limit,strain_rate)<=0.or.min(tau0,k2,nuc_c,mobility)<0) &
      error stop 'Invalid or missing DRX material parameter'
  if (orientation_rule<0.or.orientation_rule>1.or.snapshot_every<1) &
      error stop 'Invalid DRX output/orientation option'
  if ((enabled/=0.and.enabled/=1).or.line_factor<=0) error stop 'Invalid DRX switch or line factor'
  random_state=1_int64+modulo(int(seed_value,int64)-1_int64,2147483646_int64)
  ready=.false.
  next_snapshot=1
  if (rho_nucleus==0) error stop 'Nucleus density must be positive or default -1'
  if (initial_strength_excess<0) error stop 'Initial population extra resistance cannot be negative'
  if (material_kind==40) then
    nrs(1:nsyst(1),1)=20
  else
    nm=nmodes(1)
    if (material_kind/=41.or.nm<1.or.nm>nmodmx.or.ntwmod(1)/=0) error stop 'Law 41 requires slip modes'
    allocate(mode_density(nm,ngrmx),mode_shear(nm,ngrmx))
    if (minval(nsm(:nm,1))<1.or.minval(mode_tau(:nm))<0.or. &
        minval(mode_b_ratio(:nm))<=0.or.minval(mode_chi_ratio(:nm))<=0.or.slip_exponent<2.or. &
        reference_macro_rate<=0) error stop 'Invalid HCP mode or rate parameters'
    nrs(1:nsyst(1),1)=slip_exponent
    ! Legacy fields are accepted only at zero to prevent silent old-model runs.
    if (strength_rate_m/=0.or.initial_strength_excess/=0.or.mobility_rate_m/=0.or. &
        recovery_rate_m/=0.or.gos_shear_halfweight>0) &
        error stop 'Removed empirical Mg options: use the paper-revision parameter file'
    if (minval(mode_k1(:nm))<=0.or.minval(mode_k2(:nm))<0.or.minval(mode_rho_initial(:nm))<=0.or. &
        minval(mode_rho_nucleus(:nm))<=0) error stop 'Specify KM coefficients and densities for every mode'
    mode_k2(:nm)=mode_k2(:nm)*(strain_rate/reference_macro_rate)**(-mode_recovery_rate_m(:nm))
    rho_initial=sum(mode_rho_initial(:nm)); rho_nucleus=sum(mode_rho_nucleus(:nm))
    last=0
    do im=1,nmodes(1)
      do is=last+1,last+nsm(im,1)
        system_tau(is)=mode_tau(im)
        system_b_ratio(is)=mode_b_ratio(im)*mode_chi_ratio(im)
        system_mode(is)=im
      end do
      last=last+nsm(im,1)
    end do
  end if
  if (material_kind==40) then
    write(10,'(a)') 'Zhou et al. (2017) VPSC-DDRX, IHARDLAW=40'
  else
    write(10,'(a,i0)') 'Mode-resolved KM-DDRX, IHARDLAW=41. Active modes: ',nm
  end if
  write(10,nml=drx)
end subroutine

subroutine drx_prepare
  use vpsc8dim
  use drx_state
  implicit none
  integer :: g, saved_interaction, u, count_sizes
  real(real64) :: reference_rate, resistance, norm_rate, isotropic(6,6),lam,gmod
  integer :: j,k
  if (.not.ready) then
    if (ishape(1)/=0.or.iupdshp/=0.or.nneigh/=0.or.iflu/=0.or. &
        irecover/=0.or.isave/=0.or.ntwmod(1)/=0) &
        error stop 'DRX needs spherical grains, no neighbours, twins or restart'
    if (material_kind==40.and.(nsyst(1)/=12.or.nmodes(1)/=1)) &
        error stop 'Copper DRX requires 12 FCC systems in one mode'
    if (interaction/=3.or.abs(neffgrx-5.0_real64)>1.e-10_real64) &
        error stop 'Published DRX case requires interaction=3, meff=5'
    if (iupdhar/=1) error stop 'DRX requires hardening updates'
    if (abs(temp-temperature)>1.e-6_real64) error stop 'DRX temperature differs from process file'
    initial_count=ngr(1)
    ! In a one-phase isotropic elastic aggregate all rotations and weights
    ! give exactly the same elastic localization. Avoid rotating fourth-order
    ! elastic tensors at every plastic step; retain the original ELSC path
    ! for an anisotropic input. This changes no viscoplastic equation.
    lam=celccv(1,2,1); gmod=celccv(4,4,1); isotropic=0
    isotropic(1:3,1:3)=lam
    do j=1,3
      isotropic(j,j)=lam+2*gmod
      isotropic(j+3,j+3)=gmod
    end do
    elastic_cached=maxval(abs(celccv(:,:,1)-isotropic))<1.e-9_real64*maxval(abs(isotropic))
    box_radius_cubed=real(initial_count,real64)*(diameter_initial/2)**3
    do g=1,ngr(1)
      density(g)=rho_initial
      if (material_kind==41) then
        mode_density(:,g)=mode_rho_initial(:nmodes(1))
        mode_shear(:,g)=0
        has_nucleated(g)=.false.; orientation_parent_id(g)=0
      end if
      radius(g)=(wgt(g)*box_radius_cubed)**(1.0_real64/3)
      generation(g)=0; grain_id(g)=g; parent_id(g)=0; birth_strain(g)=0
    end do
    multiplicity(1:ngr(1))=1
    if (material_kind==41.and.len_trim(initial_size_file)>0) then
      open(newunit=u,file=trim(initial_size_file),status='old',action='read')
      read(u,*) count_sizes
      if (count_sizes/=ngr(1)) error stop 'Size file and texture count differ'
      do g=1,ngr(1)
        read(u,*) radius(g)
      end do
      close(u)
      radius(1:ngr(1))=radius(1:ngr(1))*0.5_real64 ! File holds diameters in m.
      if (minval(radius(1:ngr(1)))<=0) error stop 'Nonpositive initial size'
      ! EBSD area fractions are volume proxies; measured radii are independent.
      ! Each weighted orientation represents a (possibly fractional) grain cohort.
      box_radius_cubed=real(ngr(1),real64)/sum(wgt(1:ngr(1))/radius(1:ngr(1))**3)
      multiplicity(1:ngr(1))=wgt(1:ngr(1))*box_radius_cubed/radius(1:ngr(1))**3
    end if
    if (material_kind==41) then
      open(newunit=u,file='MG_INITIAL_GRAINS.csv',status='replace')
      write(u,'(a)') 'solver_id,diameter_um,volume_fraction,multiplicity,rho_m-2'
      do g=1,ngr(1)
        write(u,'(i0,4(",",es23.15))') g,2.e6_real64*radius(g),wgt(g),multiplicity(g),density(g)
      end do
      close(u)
    end if
    next_id=ngr(1); created=0; removed=0
    newborn(1:ngr(1))=.false.
    open(newunit=history_unit,file='DRX_HISTORY.csv',status='replace')
    write(history_unit,'(a)') 'strain,time_s,stress_MPa,drx_fraction,mean_diameter_number_um,'// &
      'mean_diameter_volume_um,rho_mean_m-2,grain_count,created,removed,max_nucleation_probability,volume_scale'
    open(newunit=cycles_unit,file='DRX_CYCLES.csv',status='replace')
    write(cycles_unit,'(a)') 'strain,generation_0,generation_1,generation_2,generation_3,'// &
      'generation_4,generation_5,generation_6,generation_7,generation_8,generation_9,generation_10plus'
    if (material_kind==41) then
      open(newunit=mg_history_unit,file='MG_HISTORY.csv',status='replace')
      write(mg_history_unit,'(a)',advance='no') 'axial_strain,axial_stress_MPa'
      call numbered_columns(mg_history_unit,'family_',nmodes(1))
      write(mg_history_unit,'(a)') ',cohort_multiplicity_sum'
      open(newunit=system_activity_unit,file='SYSTEM_ACTIVITY.csv',status='replace')
      write(system_activity_unit,'(a)',advance='no') 'axial_strain'
      call numbered_columns(system_activity_unit,'system_',nsyst(1))
      write(system_activity_unit,'()')
      open(newunit=mg_mode_unit,file='MG_MODE_HISTORY.csv',status='replace')
      write(mg_mode_unit,'(a)',advance='no') 'axial_strain'
      call numbered_columns(mg_mode_unit,'rho_',nmodes(1))
      write(mg_mode_unit,'(a)',advance='no') ',rho_total'
      call numbered_columns(mg_mode_unit,'crss_MPa_',nmodes(1))
      write(mg_mode_unit,'()')
      births_by_generation=0; consumed_by_generation=0
      births_by_generation(0)=ngr(1)
      open(newunit=mg_cycle_unit,file='MG_CYCLE_KINETICS.csv',status='replace')
      write(mg_cycle_unit,'(a)') 'axial_strain,generation,volume_fraction,rho_volume_mean_m-2,'// &
        'mean_shear_since_birth,nucleation_eligible_volume_fraction,births_cumulative,consumed_cumulative'
      if (snapshot_strains(1)>=0) then
        open(newunit=snapshot_index_unit,file='MG_SNAPSHOTS.csv',status='replace')
        write(snapshot_index_unit,'(a)') 'target_strain,actual_axial_strain,filename'
      end if
    end if
    ready=.true.
  end if
  if (abs(temp-temperature)>1.e-6_real64) error stop 'Only isothermal DRX is configured'
  ! Table 1: the denominator is Tm (1356 K), NOT Tm-300.
  if (material_kind==40) then
    mu=35.4e9_real64*(1-0.5_real64*(temp-300)/1356.0_real64)
  else
    mu=mu_reference*(1-mu_slope*(temp-reference_temperature))
    if (mu<=0) error stop 'Nonpositive shear modulus'
  end if
  line_energy=line_factor*mu*burgers**2
  reference_rate=sqrt(2.0_real64/3)*sqrt(sum(dbar(1:5)**2))
  if (material_kind==41) then
    ! Mg calibration's configured rate is the prescribed axial compression rate.
    ! Free lateral strains may differ in a textured aggregate, so its equivalent
    ! rate need not equal that axial rate. Keep Eq. 3's slip reference equivalent.
    if (abs(abs(dbarc(3,3))-strain_rate)>0.02_real64*strain_rate) &
      error stop 'Axial compression rate differs from Mg DRX configuration'
  else
    if (abs(reference_rate-strain_rate)>0.02_real64*strain_rate) &
      error stop 'Process strain rate differs from DRX configuration'
  end if
  norm_rate=sqrt(sum(dbar(1:6)**2))
  do g=1,ngr(1)
    resistance=tau0+chi*burgers*mu*sqrt(density(g))*1.e-6_real64
    crss(1:nsyst(1),g)=resistance
    if (material_kind==41) then
      do j=1,nsyst(1)
        crss(j,g)=mode_resistance(system_mode(j),g)
      end do
    end if
    gamd0s(1:nsyst(1),g)=reference_rate ! Eq. 3, macroscopic equivalent strain rate
    if (newborn(g)) then
      ! Initialize the new grain's moduli from its own orientation and CRSS.
      sgtry(:,g)=dbar(:)/norm_rate
      saved_interaction=interaction
      interaction=0
      call grain_stress(g,g,1,1)
      interaction=saved_interaction
      call grain_rate_and_moduli(0,g,g,1,1)
      newborn(g)=.false.
    end if
  end do
end subroutine

subroutine drx_update_elastic
  use drx_state, only: elastic_cached
  implicit none
  if (.not.elastic_cached) call elsc(1)
end subroutine

subroutine drx_harden
  use vpsc8dim
  use drx_state
  implicit none
  integer :: g,im,is,last
  real(real64) :: shear
  if (time_incr>dt_limit*(1+1.e-5_real64)) error stop 'DRX time increment exceeds configured limit'
  do g=1,ngr(1)
    shear=sum(abs(gamdot(1:nsyst(1),g)))*time_incr
    if (material_kind==41) gtotgr(g)=gtotgr(g)+shear
    if (material_kind==41) then
      last=0
      do im=1,nmodes(1)
        shear=sum(abs(gamdot(last+1:last+nsm(im,1),g)))*time_incr
        mode_shear(im,g)=mode_shear(im,g)+shear
        mode_density(im,g)=km_density(mode_density(im,g),shear,mode_k1(im),mode_k2(im))
        crss(last+1:last+nsm(im,1),g)=mode_resistance(im,g)
        last=last+nsm(im,1)
      end do
      density(g)=sum(mode_density(:,g)) ! Section 2.2: total density sums modes.
    else
      density(g)=km_density(density(g),shear)
      crss(1:nsyst(1),g)=tau0+chi*burgers*mu*sqrt(density(g))*1.e-6_real64
    end if
  end do
  if (enabled==1) call drx_population
end subroutine

subroutine drx_population
  use vpsc8dim
  use drx_state
  implicit none
  integer :: g, old_count, n, live, bucket
  real(real64) :: probability, rho_average, cube_sum, factor
  old_count=ngr(1)
  max_probability=0
  do g=1,old_count
    if (density(g)<rhoc) cycle
    probability=nuc_c*strain_rate**q*4*acos(-1.0_real64)*radius(g)**2*time_incr
    max_probability=max(max_probability,probability)
    if (probability>=1) error stop 'Eq.10 probability >=1: reduce time increment'
    if (uniform()<probability) call drx_nucleate(g)
  end do
  n=ngr(1)
  ! Append embryos at the fixed physical size, then normalize their weights.
  cube_sum=sum(multiplicity(1:n)*radius(1:n)**3)
  wgt(1:n)=multiplicity(1:n)*radius(1:n)**3/cube_sum
  rho_average=sum(wgt(1:n)*density(1:n))
  radius(1:n)=max(0.0_real64,radius(1:n)+ &
    mobility*line_energy*(rho_average-density(1:n))*time_incr)
  ! Remove consumed grains and pack every state array used by this branch.
  live=0
  do g=1,n
    if (radius(g)<=0) then
      removed=removed+1
      if (material_kind==41) then
        bucket=min(generation(g),10)
        consumed_by_generation(bucket)=consumed_by_generation(bucket)+1
      end if
      cycle
    end if
    live=live+1
    if (live/=g) call drx_copy_grain(g,live)
  end do
  if (live==0) error stop 'All grains consumed; reduce DRX time increment'
  ngr(1)=live
  cube_sum=sum(multiplicity(1:live)*radius(1:live)**3)
  volume_scale=cube_sum/box_radius_cubed
  wgt(1:live)=multiplicity(1:live)*radius(1:live)**3/cube_sum
  ! Restore the fixed aggregate volume. This closes the Appendix A1 mapping.
  factor=(box_radius_cubed/cube_sum)**(1.0_real64/3)
  radius(1:live)=radius(1:live)*factor
end subroutine

subroutine drx_nucleate(parent)
  use vpsc8dim
  use drx_state
  implicit none
  integer, intent(in) :: parent
  integer :: child,j,base,choice,pool_count
  real(real64) :: a(3,3),axis(3),rotation(3,3),cross(3,3),angle,z,azimuth
  real(real64) :: e1,e2,e3
  child=ngr(1)+1
  if (child>min(ngrmx,ngrpel)) error stop 'DRX grain capacity exceeded: increase NGRMX/NGRPEL'
  call drx_copy_grain(parent,child)
  base=parent
  if (material_kind==41) then
    has_nucleated(parent)=.true.
    has_nucleated(child)=.false.
    orientation_parent_id(child)=0
  end if
  if (orientation_rule==0) then
    e1=360*uniform(); e2=acos(2*uniform()-1)*180/acos(-1.0_real64); e3=360*uniform()
    call euler(2,e1,e2,e3,a)
    ag(:,:,child)=transpose(a)
  else
    ! Section 4.4: randomly select a surviving matrix grain that has nucleated.
    ! Cu legacy branch is unchanged; revised Mg uses the actual selection pool.
    if (material_kind==41) then
      pool_count=count(has_nucleated(1:ngr(1)))
      choice=1+int(uniform()*pool_count)
      do j=1,ngr(1)
        if (.not.has_nucleated(j)) cycle
        choice=choice-1
        if (choice==0) then
          base=j
          exit
        end if
      end do
      orientation_parent_id(child)=grain_id(base)
    end if
    z=2*uniform()-1; azimuth=2*acos(-1.0_real64)*uniform()
    axis=[sqrt(1-z*z)*cos(azimuth),sqrt(1-z*z)*sin(azimuth),z]
    angle=sqrt(-2*log(uniform()))*cos(2*acos(-1.0_real64)*uniform())* &
      orientation_spread*acos(-1.0_real64)/180
    cross=0
    cross(1,2)=-axis(3); cross(1,3)=axis(2); cross(2,1)=axis(3)
    cross(2,3)=-axis(1); cross(3,1)=-axis(2); cross(3,2)=axis(1)
    rotation=cos(angle)*xid3+(1-cos(angle))*spread(axis,2,3)*spread(axis,1,3)+sin(angle)*cross
    ag(:,:,child)=matmul(rotation,ag(:,:,base))
  end if
  density(child)=rho_initial
  if (material_kind==41.and.rho_nucleus>0) density(child)=rho_nucleus
  if (material_kind==41) then
    mode_density(:,child)=mode_rho_nucleus(:nmodes(1))
    mode_shear(:,child)=0
    density(child)=sum(mode_density(:,child))
  end if
  radius(child)=diameter_nucleus/2
  generation(child)=generation(parent)+1
  if (material_kind==41) then
    j=min(generation(child),10)
    births_by_generation(j)=births_by_generation(j)+1
  end if
  next_id=next_id+1; grain_id(child)=next_id; parent_id(child)=grain_id(parent)
  birth_strain(child)=epsvm; newborn(child)=.true.
  gtotgr(child)=0; epsvmgr(child)=0; workgr(child)=0
  gamdot(:,child)=0; dg_trans(:,child)=0
  fijgr(:,:,child)=xid3
  do j=1,3
    axisgr(0,j,child)=1
    axisgr(1:3,j,child)=xid3(:,j)
  end do
  neffgr(child)=neffgrx
  crss(:,child)=tau0+chi*burgers*mu*sqrt(density(child))*1.e-6_real64
  if (material_kind==41) then
    do j=1,nsyst(1)
      crss(j,child)=mode_resistance(system_mode(j),child)
    end do
  end if
  created=created+1; ngr(1)=child
end subroutine

subroutine drx_copy_grain(src,dst)
  use vpsc8dim
  use drx_state
  implicit none
  integer, intent(in) :: src,dst
  density(dst)=density(src); radius(dst)=radius(src)
  if (material_kind==41) then
    mode_density(:,dst)=mode_density(:,src); mode_shear(:,dst)=mode_shear(:,src)
    has_nucleated(dst)=has_nucleated(src); orientation_parent_id(dst)=orientation_parent_id(src)
  end if
  multiplicity(dst)=multiplicity(src)
  generation(dst)=generation(src); grain_id(dst)=grain_id(src)
  parent_id(dst)=parent_id(src); birth_strain(dst)=birth_strain(src); newborn(dst)=newborn(src)
  ag(:,:,dst)=ag(:,:,src); wgt(dst)=wgt(src)
  axisgr(:,:,dst)=axisgr(:,:,src); fijgr(:,:,dst)=fijgr(:,:,src)
  lijgr(:,:,dst)=lijgr(:,:,src); asgr(:,:,:,:,dst)=asgr(:,:,:,:,src)
  mastgr(:,:,dst)=mastgr(:,:,src)
  crss(:,dst)=crss(:,src); taue(:,dst)=taue(:,src)
  sg(:,dst)=sg(:,src); sg0(:,dst)=sg0(:,src); sgtry(:,dst)=sgtry(:,src)
  dg(:,dst)=dg(:,src); dg0(:,dst)=dg0(:,src); dg_trans(:,dst)=dg_trans(:,src)
  gtotgr(dst)=gtotgr(src); gamd0s(:,dst)=gamd0s(:,src)
  gamdot(:,dst)=gamdot(:,src); sch(:,:,dst)=sch(:,:,src)
  mgsec(:,:,dst)=mgsec(:,:,src); mgtg(:,:,dst)=mgtg(:,:,src)
  bg(:,:,dst)=bg(:,:,src); neffgr(dst)=neffgr(src)
  epsvmgr(dst)=epsvmgr(src); workgr(dst)=workgr(src)
  sgvm(dst)=sgvm(src); taylgr(dst)=taylgr(src)
  gavgr(dst)=gavgr(src); actgr(dst)=actgr(src)
  celcs(:,:,dst)=celcs(:,:,src); selcs(:,:,dst)=selcs(:,:,src)
  athcs(:,dst)=athcs(:,src); ag_el(:,:,dst)=ag_el(:,:,src)
  bg_el(:,:,dst)=bg_el(:,:,src); tg_el(:,dst)=tg_el(:,src)
  dg_el(:,dst)=dg_el(:,src); eelgr(:,dst)=eelgr(:,src)
  ktwsmx(dst)=0; ntwevents(dst)=0; twfrsy(:,dst)=0
  ! VPSC uses its neighbour table even when NNEIGH=0 (self rotation).
  neigh(0,dst)=dst; wneigh(0,dst)=1
end subroutine

subroutine drx_record(step)
  use vpsc8dim
  use drx_state
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, intent(in) :: step
  integer :: n,g,k,u,first,last
  character(80) :: state_filename
  real(real64) :: fractions(0:10),x,mean_number,mean_volume,rhobar,e1,e2,e3,a(3,3),activity(nmodes(1)),total_activity
  real(real64) :: rho_modes(nmodes(1)),tau_modes(nmodes(1)),system_activity(nsyst(1))
  real(real64) :: generation_rho(0:10),generation_shear(0:10),eligible_fraction(0:10)
  n=ngr(1)
  if (.not.ieee_is_finite(svm).or.any(.not.ieee_is_finite(density(1:n))).or. &
      any(.not.ieee_is_finite(wgt(1:n)))) error stop 'Nonfinite DRX state'
  if (abs(sum(wgt(1:n))-1)>1.e-9_real64.or.minval(wgt(1:n))<=0) &
    error stop 'DRX volume weights invalid'
  fractions=0
  generation_rho=0; generation_shear=0; eligible_fraction=0
  do g=1,n
    k=min(generation(g),10)
    fractions(k)=fractions(k)+wgt(g)
    if (material_kind==41) then
      generation_rho(k)=generation_rho(k)+wgt(g)*density(g)
      generation_shear(k)=generation_shear(k)+wgt(g)*gtotgr(g)
      if (density(g)>=rhoc) eligible_fraction(k)=eligible_fraction(k)+wgt(g)
    end if
  end do
  x=max(0.0_real64,min(1.0_real64,1-fractions(0)))
  mean_number=2.e6_real64*sum(radius(1:n))/n
  if (material_kind==41) mean_number=2.e6_real64*sum(radius(1:n)*multiplicity(1:n))/sum(multiplicity(1:n))
  mean_volume=2.e6_real64*sum(radius(1:n)*wgt(1:n))
  rhobar=sum(density(1:n)*wgt(1:n))
  write(history_unit,'(7(es23.15,","),3(i0,","),es23.15,",",es23.15)') &
    epsvm,time,svm,x,mean_number,mean_volume,rhobar,n,created,removed,max_probability,volume_scale
  write(cycles_unit,'(11(es23.15,","),es23.15)') epsvm,fractions
  if (material_kind==41) then
    if (any(.not.ieee_is_finite(mode_density(:,1:n))).or.minval(mode_density(:,1:n))<0) &
      error stop 'Invalid mode density'
    if (maxval(abs(density(1:n)-sum(mode_density(:,1:n),dim=1)))>1.e-12_real64*maxval(density(1:n))) &
      error stop 'Total density does not equal sum of modes'
    rho_modes=matmul(mode_density(:,1:n),wgt(1:n))
    tau_modes=0
    do g=1,n
      do k=1,nmodes(1)
        tau_modes(k)=tau_modes(k)+wgt(g)*mode_resistance(k,g)
      end do
    end do
    write(mg_mode_unit,'(*(es23.15,:,","))') time*strain_rate,rho_modes,rhobar,tau_modes
    flush(mg_mode_unit)
    do k=0,10
      if (fractions(k)>0) then
        write(mg_cycle_unit,'(es23.15,",",i0,4(",",es23.15),2(",",i0))') time*strain_rate,k, &
          fractions(k),generation_rho(k)/fractions(k),generation_shear(k)/fractions(k), &
          eligible_fraction(k),births_by_generation(k),consumed_by_generation(k)
      end if
    end do
    flush(mg_cycle_unit)
    activity=0; system_activity=0
    do g=1,n
      first=1
      do k=1,nmodes(1)
        last=first+nsm(k,1)-1
        activity(k)=activity(k)+wgt(g)*sum(abs(gamdot(first:last,g)))
        first=last+1
      end do
      system_activity=system_activity+wgt(g)*abs(gamdot(1:nsyst(1),g))
    end do
    total_activity=sum(activity)
    if (total_activity>0) then
      activity=activity/total_activity
      system_activity=system_activity/total_activity
    end if
    write(system_activity_unit,'(*(es23.15,:,","))') time*strain_rate,system_activity
    write(mg_history_unit,'(*(es23.15,:,","))') time*strain_rate,-sbarc(3,3), &
      activity,sum(multiplicity(1:n))
    flush(mg_history_unit)
    do while (next_snapshot<=size(snapshot_strains))
      if (snapshot_strains(next_snapshot)<0) exit
      if (time*strain_rate+1.e-9_real64<snapshot_strains(next_snapshot)) exit
      write(state_filename,'("MG_STATE_",i3.3,".csv")') next_snapshot
      call drx_write_state(trim(state_filename))
      write(snapshot_index_unit,'(2(es23.15,","),a)') &
        snapshot_strains(next_snapshot),time*strain_rate,trim(state_filename)
      flush(snapshot_index_unit)
      next_snapshot=next_snapshot+1
    end do
  end if
  flush(history_unit); flush(cycles_unit)
  if (mod(step-1,snapshot_every)/=0.and.step/=nsteps+1) return
  open(newunit=u,file='DRX_GRAINS.csv',status='replace')
  write(u,'(a)') 'id,parent_id,generation,birth_strain,diameter_um,volume_fraction,rho_m-2,phi1,Phi,phi2'
  do g=1,n
    a=transpose(ag(:,:,g)); call euler(1,e1,e2,e3,a)
    write(u,'(3(i0,","),6(es23.15,","),es23.15)') grain_id(g),parent_id(g),generation(g), &
      birth_strain(g),2.e6_real64*radius(g),wgt(g),density(g),e1,e2,e3
  end do
  close(u)
  if (material_kind==41) then
    open(newunit=u,file='MG_GRAIN_COHORTS.csv',status='replace')
    write(u,'(a)') 'id,multiplicity,radius_m'
    do g=1,n
      write(u,'(i0,2(",",es23.15))') grain_id(g),multiplicity(g),radius(g)
    end do
    close(u)
  end if
end subroutine

subroutine drx_write_state(filename)
  use vpsc8dim
  use drx_state
  implicit none
  character(*),intent(in) :: filename
  integer :: u,g
  real(real64) :: a(3,3),e1,e2,e3
  open(newunit=u,file=filename,status='replace')
  write(u,'(a)',advance='no') 'id,parent_id,generation,birth_strain,diameter_um,volume_fraction,rho_m-2,'// &
    'phi1,Phi,phi2,accumulated_shear,multiplicity'
  call numbered_columns(u,'rho_',nmodes(1))
  call numbered_columns(u,'shear_',nmodes(1))
  write(u,'(a)') ',orientation_base_id'
  do g=1,ngr(1)
    a=transpose(ag(:,:,g)); call euler(1,e1,e2,e3,a)
    write(u,'(3(i0,","),*(es23.15,","))',advance='no') grain_id(g),parent_id(g),generation(g), &
      birth_strain(g),2.e6_real64*radius(g),wgt(g),density(g),e1,e2,e3,gtotgr(g),multiplicity(g), &
      mode_density(:,g),mode_shear(:,g)
    write(u,'(i0)') orientation_parent_id(g)
  end do
  close(u)
end subroutine
