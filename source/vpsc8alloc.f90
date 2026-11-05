! Allocate legacy grain storage on the heap; preserve all lower bounds.
subroutine allocate_grain_storage
  use parameters
  use VPSC8DIM, only: alloc_VPSC8DIM_AG => AG
  use VPSC8DIM, only: alloc_VPSC8DIM_WGT => WGT
  use VPSC8DIM, only: alloc_VPSC8DIM_AXISGR => AXISGR
  use VPSC8DIM, only: alloc_VPSC8DIM_FIJGR => FIJGR
  use VPSC8DIM, only: alloc_VPSC8DIM_LIJGR => LIJGR
  use VPSC8DIM, only: alloc_VPSC8DIM_ASGR => ASGR
  use VPSC8DIM, only: alloc_VPSC8DIM_MASTGR => MASTGR
  use VPSC8DIM, only: alloc_VPSC8DIM_CRSS => CRSS
  use VPSC8DIM, only: alloc_VPSC8DIM_TAUE => TAUE
  use VPSC8DIM, only: alloc_VPSC8DIM_SG => SG
  use VPSC8DIM, only: alloc_VPSC8DIM_SG0 => SG0
  use VPSC8DIM, only: alloc_VPSC8DIM_SGTRY => SGTRY
  use VPSC8DIM, only: alloc_VPSC8DIM_DG => DG
  use VPSC8DIM, only: alloc_VPSC8DIM_DG0 => DG0
  use VPSC8DIM, only: alloc_VPSC8DIM_DG_TRANS => DG_TRANS
  use VPSC8DIM, only: alloc_VPSC8DIM_GTOTGR => GTOTGR
  use VPSC8DIM, only: alloc_VPSC8DIM_GAMD0S => GAMD0S
  use VPSC8DIM, only: alloc_VPSC8DIM_GAMDOT => GAMDOT
  use VPSC8DIM, only: alloc_VPSC8DIM_SCH => SCH
  use VPSC8DIM, only: alloc_VPSC8DIM_MGSEC => MGSEC
  use VPSC8DIM, only: alloc_VPSC8DIM_MGTG => MGTG
  use VPSC8DIM, only: alloc_VPSC8DIM_BG => BG
  use VPSC8DIM, only: alloc_VPSC8DIM_NEFFGR => NEFFGR
  use VPSC8DIM, only: alloc_VPSC8DIM_SGVM => SGVM
  use VPSC8DIM, only: alloc_VPSC8DIM_EPSVMGR => EPSVMGR
  use VPSC8DIM, only: alloc_VPSC8DIM_WORKGR => WORKGR
  use VPSC8DIM, only: alloc_VPSC8DIM_TAYLGR => TAYLGR
  use VPSC8DIM, only: alloc_VPSC8DIM_GAVGR => GAVGR
  use VPSC8DIM, only: alloc_VPSC8DIM_ACTGR => ACTGR
  use VPSC8DIM, only: alloc_VPSC8DIM_TWFRSY => TWFRSY
  use VPSC8DIM, only: alloc_VPSC8DIM_KTWSMX => KTWSMX
  use VPSC8DIM, only: alloc_VPSC8DIM_NTWEVENTS => NTWEVENTS
  use VPSC8DIM, only: alloc_VPSC8DIM_KTWMMX => KTWMMX
  use VPSC8DIM, only: alloc_VPSC8DIM_ITWVFT => ITWVFT
  use VPSC8DIM, only: alloc_VPSC8DIM_TWFRGR => TWFRGR
  use VPSC8DIM, only: alloc_VPSC8DIM_IPTSGR => IPTSGR
  use VPSC8DIM, only: alloc_VPSC8DIM_KPTMGR => KPTMGR
  use VPSC8DIM, only: alloc_VPSC8DIM_WNEIGH => WNEIGH
  use VPSC8DIM, only: alloc_VPSC8DIM_NEIGH => NEIGH
  use VPSC8DIM, only: alloc_VPSC8DIM_CELCS => CELCS
  use VPSC8DIM, only: alloc_VPSC8DIM_SELCS => SELCS
  use VPSC8DIM, only: alloc_VPSC8DIM_ATHCS => ATHCS
  use VPSC8DIM, only: alloc_VPSC8DIM_EELGR => EELGR
  use VPSC8DIM, only: alloc_VPSC8DIM_AG_EL => AG_EL
  use VPSC8DIM, only: alloc_VPSC8DIM_BG_EL => BG_EL
  use VPSC8DIM, only: alloc_VPSC8DIM_TG_EL => TG_EL
  use VPSC8DIM, only: alloc_VPSC8DIM_DG_EL => DG_EL
  use VPSC8DIM, only: alloc_VPSC8DIM_SEQ2 => SEQ2
  use VPSC8DIM, only: alloc_VPSC8DIM_DEQ2 => DEQ2
  use VPSC8DIM, only: alloc_VPSC8DIM_ASO => ASO
  use VPSC8DIM, only: alloc_VPSC8DIM_ESO => ESO
  use VPSC8DIM, only: alloc_VPSC8DIM_SECMOM5 => SECMOM5
  use VPSC8DIM, only: alloc_VPSC8DIM_SECMOM5D => SECMOM5D
  use VPSC8DIM, only: alloc_VPSC8DIM_BETFLU => BETFLU
  use VPSC8DIM, only: alloc_VPSC8DIM_CHIFLU => CHIFLU
  use CUB_COMP, only: alloc_CUB_COMP_igrtype => igrtype
  use DIFFRACT, only: alloc_DIFFRACT_IGRSET => IGRSET
  use DISDEN, only: alloc_DISDEN_RHO_S => RHO_S
  use DISDEN, only: alloc_DISDEN_RHO_DEB => RHO_DEB
  use DISDEN, only: alloc_DISDEN_RHO_FOR => RHO_FOR
  use DISDEN_REV, only: alloc_DISDEN_REV_xlambda => xlambda
  use DISDEN_REV, only: alloc_DISDEN_REV_xp => xp
  use DISDEN_REV, only: alloc_DISDEN_REV_rho_f => rho_f
  use DISDEN_REV, only: alloc_DISDEN_REV_rho_r => rho_r
  use DISDEN_REV, only: alloc_DISDEN_REV_rho_prev => rho_prev
  use DISDEN_REV, only: alloc_DISDEN_REV_rho_deb => rho_deb
  use DISDEN_REV, only: alloc_DISDEN_REV_shear_tag => shear_tag
  use DISDEN_REV, only: alloc_DISDEN_REV_deltaD => deltaD
  use DISDEN_REV, only: alloc_DISDEN_REV_xp_g => xp_g
  use DISDEN_REV, only: alloc_DISDEN_REV_drho_tot => drho_tot
  use DISDEN_REV, only: alloc_DISDEN_REV_debs => debs
  use DISDEN_REV, only: alloc_DISDEN_REV_rho_ps => rho_ps
  use EXTRAPSO, only: alloc_EXTRAPSO_ASOOLD => ASOOLD
  use EXTRAPSO, only: alloc_EXTRAPSO_ESOOLD => ESOOLD
  use IRR_VARS, only: alloc_IRR_VARS_rho_irr => rho_irr
  use IRR_VARS, only: alloc_IRR_VARS_xN_i => xN_i
  use IRR_VARS, only: alloc_IRR_VARS_xr_i => xr_i
  use IRR_VARS, only: alloc_IRR_VARS_xN_v => xN_v
  use IRR_VARS, only: alloc_IRR_VARS_xr_v => xr_v
  implicit none
  integer :: status
  if (.not.allocated(alloc_VPSC8DIM_AG)) then
    allocate(alloc_VPSC8DIM_AG(3,3,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:AG'
    alloc_VPSC8DIM_AG=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_WGT)) then
    allocate(alloc_VPSC8DIM_WGT(NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:WGT'
    alloc_VPSC8DIM_WGT=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_AXISGR)) then
    allocate(alloc_VPSC8DIM_AXISGR(0:3,3,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:AXISGR'
    alloc_VPSC8DIM_AXISGR=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_FIJGR)) then
    allocate(alloc_VPSC8DIM_FIJGR(3,3,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:FIJGR'
    alloc_VPSC8DIM_FIJGR=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_LIJGR)) then
    allocate(alloc_VPSC8DIM_LIJGR(3,3,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:LIJGR'
    alloc_VPSC8DIM_LIJGR=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_ASGR)) then
    allocate(alloc_VPSC8DIM_ASGR(3,3,3,3,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:ASGR'
    alloc_VPSC8DIM_ASGR=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_MASTGR)) then
    allocate(alloc_VPSC8DIM_MASTGR(MDIM,MDIM,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:MASTGR'
    alloc_VPSC8DIM_MASTGR=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_CRSS)) then
    allocate(alloc_VPSC8DIM_CRSS(NSYSMX,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:CRSS'
    alloc_VPSC8DIM_CRSS=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_TAUE)) then
    allocate(alloc_VPSC8DIM_TAUE(NSYSMX,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:TAUE'
    alloc_VPSC8DIM_TAUE=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_SG)) then
    allocate(alloc_VPSC8DIM_SG(6,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:SG'
    alloc_VPSC8DIM_SG=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_SG0)) then
    allocate(alloc_VPSC8DIM_SG0(6,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:SG0'
    alloc_VPSC8DIM_SG0=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_SGTRY)) then
    allocate(alloc_VPSC8DIM_SGTRY(6,NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:SGTRY'
    alloc_VPSC8DIM_SGTRY=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_DG)) then
    allocate(alloc_VPSC8DIM_DG(6,NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:DG'
    alloc_VPSC8DIM_DG=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_DG0)) then
    allocate(alloc_VPSC8DIM_DG0(6,NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:DG0'
    alloc_VPSC8DIM_DG0=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_DG_TRANS)) then
    allocate(alloc_VPSC8DIM_DG_TRANS(6,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:DG_TRANS'
    alloc_VPSC8DIM_DG_TRANS=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_GTOTGR)) then
    allocate(alloc_VPSC8DIM_GTOTGR(NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:GTOTGR'
    alloc_VPSC8DIM_GTOTGR=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_GAMD0S)) then
    allocate(alloc_VPSC8DIM_GAMD0S(NSYSMX,NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:GAMD0S'
    alloc_VPSC8DIM_GAMD0S=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_GAMDOT)) then
    allocate(alloc_VPSC8DIM_GAMDOT(NSYSMX,NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:GAMDOT'
    alloc_VPSC8DIM_GAMDOT=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_SCH)) then
    allocate(alloc_VPSC8DIM_SCH(5,NSYSMX,NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:SCH'
    alloc_VPSC8DIM_SCH=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_MGSEC)) then
    allocate(alloc_VPSC8DIM_MGSEC(MDIM,MDIM,NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:MGSEC'
    alloc_VPSC8DIM_MGSEC=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_MGTG)) then
    allocate(alloc_VPSC8DIM_MGTG(MDIM,MDIM,NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:MGTG'
    alloc_VPSC8DIM_MGTG=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_BG)) then
    allocate(alloc_VPSC8DIM_BG(MDIM,MDIM,NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:BG'
    alloc_VPSC8DIM_BG=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_NEFFGR)) then
    allocate(alloc_VPSC8DIM_NEFFGR(NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:NEFFGR'
    alloc_VPSC8DIM_NEFFGR=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_SGVM)) then
    allocate(alloc_VPSC8DIM_SGVM(NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:SGVM'
    alloc_VPSC8DIM_SGVM=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_EPSVMGR)) then
    allocate(alloc_VPSC8DIM_EPSVMGR(NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:EPSVMGR'
    alloc_VPSC8DIM_EPSVMGR=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_WORKGR)) then
    allocate(alloc_VPSC8DIM_WORKGR(NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:WORKGR'
    alloc_VPSC8DIM_WORKGR=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_TAYLGR)) then
    allocate(alloc_VPSC8DIM_TAYLGR(NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:TAYLGR'
    alloc_VPSC8DIM_TAYLGR=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_GAVGR)) then
    allocate(alloc_VPSC8DIM_GAVGR(NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:GAVGR'
    alloc_VPSC8DIM_GAVGR=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_ACTGR)) then
    allocate(alloc_VPSC8DIM_ACTGR(NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:ACTGR'
    alloc_VPSC8DIM_ACTGR=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_TWFRSY)) then
    allocate(alloc_VPSC8DIM_TWFRSY(NTWSMX,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:TWFRSY'
    alloc_VPSC8DIM_TWFRSY=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_KTWSMX)) then
    allocate(alloc_VPSC8DIM_KTWSMX(NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:KTWSMX'
    alloc_VPSC8DIM_KTWSMX=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_NTWEVENTS)) then
    allocate(alloc_VPSC8DIM_NTWEVENTS(NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:NTWEVENTS'
    alloc_VPSC8DIM_NTWEVENTS=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_KTWMMX)) then
    allocate(alloc_VPSC8DIM_KTWMMX(NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:KTWMMX'
    alloc_VPSC8DIM_KTWMMX=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_ITWVFT)) then
    allocate(alloc_VPSC8DIM_ITWVFT(2,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:ITWVFT'
    alloc_VPSC8DIM_ITWVFT=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_TWFRGR)) then
    allocate(alloc_VPSC8DIM_TWFRGR(0:NTWMMX,NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:TWFRGR'
    alloc_VPSC8DIM_TWFRGR=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_IPTSGR)) then
    allocate(alloc_VPSC8DIM_IPTSGR(NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:IPTSGR'
    alloc_VPSC8DIM_IPTSGR=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_KPTMGR)) then
    allocate(alloc_VPSC8DIM_KPTMGR(NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:KPTMGR'
    alloc_VPSC8DIM_KPTMGR=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_WNEIGH)) then
    allocate(alloc_VPSC8DIM_WNEIGH(0:NNEIMX,NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:WNEIGH'
    alloc_VPSC8DIM_WNEIGH=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_NEIGH)) then
    allocate(alloc_VPSC8DIM_NEIGH(0:NNEIMX,NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:NEIGH'
    alloc_VPSC8DIM_NEIGH=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_CELCS)) then
    allocate(alloc_VPSC8DIM_CELCS(6,6,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:CELCS'
    alloc_VPSC8DIM_CELCS=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_SELCS)) then
    allocate(alloc_VPSC8DIM_SELCS(6,6,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:SELCS'
    alloc_VPSC8DIM_SELCS=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_ATHCS)) then
    allocate(alloc_VPSC8DIM_ATHCS(6,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:ATHCS'
    alloc_VPSC8DIM_ATHCS=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_EELGR)) then
    allocate(alloc_VPSC8DIM_EELGR(6,NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:EELGR'
    alloc_VPSC8DIM_EELGR=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_AG_EL)) then
    allocate(alloc_VPSC8DIM_AG_EL(6,6,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:AG_EL'
    alloc_VPSC8DIM_AG_EL=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_BG_EL)) then
    allocate(alloc_VPSC8DIM_BG_EL(6,6,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:BG_EL'
    alloc_VPSC8DIM_BG_EL=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_TG_EL)) then
    allocate(alloc_VPSC8DIM_TG_EL(6,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:TG_EL'
    alloc_VPSC8DIM_TG_EL=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_DG_EL)) then
    allocate(alloc_VPSC8DIM_DG_EL(6,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:DG_EL'
    alloc_VPSC8DIM_DG_EL=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_SEQ2)) then
    allocate(alloc_VPSC8DIM_SEQ2(NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:SEQ2'
    alloc_VPSC8DIM_SEQ2=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_DEQ2)) then
    allocate(alloc_VPSC8DIM_DEQ2(NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:DEQ2'
    alloc_VPSC8DIM_DEQ2=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_ASO)) then
    allocate(alloc_VPSC8DIM_ASO(NSYSMX,NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:ASO'
    alloc_VPSC8DIM_ASO=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_ESO)) then
    allocate(alloc_VPSC8DIM_ESO(NSYSMX,NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:ESO'
    alloc_VPSC8DIM_ESO=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_SECMOM5)) then
    allocate(alloc_VPSC8DIM_SECMOM5(5,5,NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:SECMOM5'
    alloc_VPSC8DIM_SECMOM5=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_SECMOM5D)) then
    allocate(alloc_VPSC8DIM_SECMOM5D(5,5,NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:SECMOM5D'
    alloc_VPSC8DIM_SECMOM5D=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_BETFLU)) then
    allocate(alloc_VPSC8DIM_BETFLU(5,5,NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:BETFLU'
    alloc_VPSC8DIM_BETFLU=0
  end if
  if (.not.allocated(alloc_VPSC8DIM_CHIFLU)) then
    allocate(alloc_VPSC8DIM_CHIFLU(5,5,NGRPEL),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array VPSC8DIM:CHIFLU'
    alloc_VPSC8DIM_CHIFLU=0
  end if
  if (.not.allocated(alloc_CUB_COMP_igrtype)) then
    allocate(alloc_CUB_COMP_igrtype(NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array CUB_COMP:igrtype'
    alloc_CUB_COMP_igrtype=0
  end if
  if (.not.allocated(alloc_DIFFRACT_IGRSET)) then
    allocate(alloc_DIFFRACT_IGRSET(NDIFMX,NGRPEL,NPHMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array DIFFRACT:IGRSET'
    alloc_DIFFRACT_IGRSET=0
  end if
  if (.not.allocated(alloc_DISDEN_RHO_S)) then
    allocate(alloc_DISDEN_RHO_S(NSYSMX,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array DISDEN:RHO_S'
    alloc_DISDEN_RHO_S=0
  end if
  if (.not.allocated(alloc_DISDEN_RHO_DEB)) then
    allocate(alloc_DISDEN_RHO_DEB(NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array DISDEN:RHO_DEB'
    alloc_DISDEN_RHO_DEB=0
  end if
  if (.not.allocated(alloc_DISDEN_RHO_FOR)) then
    allocate(alloc_DISDEN_RHO_FOR(NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array DISDEN:RHO_FOR'
    alloc_DISDEN_RHO_FOR=0
  end if
  if (.not.allocated(alloc_DISDEN_REV_xlambda)) then
    allocate(alloc_DISDEN_REV_xlambda(NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array DISDEN_REV:xlambda'
    alloc_DISDEN_REV_xlambda=0
  end if
  if (.not.allocated(alloc_DISDEN_REV_xp)) then
    allocate(alloc_DISDEN_REV_xp(nsysmx,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array DISDEN_REV:xp'
    alloc_DISDEN_REV_xp=0
  end if
  if (.not.allocated(alloc_DISDEN_REV_rho_f)) then
    allocate(alloc_DISDEN_REV_rho_f(NSYSMX,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array DISDEN_REV:rho_f'
    alloc_DISDEN_REV_rho_f=0
  end if
  if (.not.allocated(alloc_DISDEN_REV_rho_r)) then
    allocate(alloc_DISDEN_REV_rho_r(2*NSYSMX,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array DISDEN_REV:rho_r'
    alloc_DISDEN_REV_rho_r=0
  end if
  if (.not.allocated(alloc_DISDEN_REV_rho_prev)) then
    allocate(alloc_DISDEN_REV_rho_prev(NSYSMX,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array DISDEN_REV:rho_prev'
    alloc_DISDEN_REV_rho_prev=0
  end if
  if (.not.allocated(alloc_DISDEN_REV_rho_deb)) then
    allocate(alloc_DISDEN_REV_rho_deb(nsysmx,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array DISDEN_REV:rho_deb'
    alloc_DISDEN_REV_rho_deb=0
  end if
  if (.not.allocated(alloc_DISDEN_REV_shear_tag)) then
    allocate(alloc_DISDEN_REV_shear_tag(NSYSMX,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array DISDEN_REV:shear_tag'
    alloc_DISDEN_REV_shear_tag=0
  end if
  if (.not.allocated(alloc_DISDEN_REV_deltaD)) then
    allocate(alloc_DISDEN_REV_deltaD(NSYSMX,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array DISDEN_REV:deltaD'
    alloc_DISDEN_REV_deltaD=0
  end if
  if (.not.allocated(alloc_DISDEN_REV_xp_g)) then
    allocate(alloc_DISDEN_REV_xp_g(NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array DISDEN_REV:xp_g'
    alloc_DISDEN_REV_xp_g=0
  end if
  if (.not.allocated(alloc_DISDEN_REV_drho_tot)) then
    allocate(alloc_DISDEN_REV_drho_tot(NSYSMX,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array DISDEN_REV:drho_tot'
    alloc_DISDEN_REV_drho_tot=0
  end if
  if (.not.allocated(alloc_DISDEN_REV_debs)) then
    allocate(alloc_DISDEN_REV_debs(NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array DISDEN_REV:debs'
    alloc_DISDEN_REV_debs=0
  end if
  if (.not.allocated(alloc_DISDEN_REV_rho_ps)) then
    allocate(alloc_DISDEN_REV_rho_ps(NSYSMX,NGRMX),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array DISDEN_REV:rho_ps'
    alloc_DISDEN_REV_rho_ps=0
  end if
  if (.not.allocated(alloc_EXTRAPSO_ASOOLD)) then
    allocate(alloc_EXTRAPSO_ASOOLD(NSYSMX,NGRPEL,3),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array EXTRAPSO:ASOOLD'
    alloc_EXTRAPSO_ASOOLD=0
  end if
  if (.not.allocated(alloc_EXTRAPSO_ESOOLD)) then
    allocate(alloc_EXTRAPSO_ESOOLD(NSYSMX,NGRPEL,3),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array EXTRAPSO:ESOOLD'
    alloc_EXTRAPSO_ESOOLD=0
  end if
  if (.not.allocated(alloc_IRR_VARS_rho_irr)) then
    allocate(alloc_IRR_VARS_rho_irr(18,ngrmx),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array IRR_VARS:rho_irr'
    alloc_IRR_VARS_rho_irr=0
  end if
  if (.not.allocated(alloc_IRR_VARS_xN_i)) then
    allocate(alloc_IRR_VARS_xN_i(4,ngrmx),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array IRR_VARS:xN_i'
    alloc_IRR_VARS_xN_i=0
  end if
  if (.not.allocated(alloc_IRR_VARS_xr_i)) then
    allocate(alloc_IRR_VARS_xr_i(4,ngrmx),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array IRR_VARS:xr_i'
    alloc_IRR_VARS_xr_i=0
  end if
  if (.not.allocated(alloc_IRR_VARS_xN_v)) then
    allocate(alloc_IRR_VARS_xN_v(4,ngrmx),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array IRR_VARS:xN_v'
    alloc_IRR_VARS_xN_v=0
  end if
  if (.not.allocated(alloc_IRR_VARS_xr_v)) then
    allocate(alloc_IRR_VARS_xr_v(4,ngrmx),stat=status)
    if (status/=0) error stop 'Cannot allocate grain array IRR_VARS:xr_v'
    alloc_IRR_VARS_xr_v=0
  end if
end subroutine
