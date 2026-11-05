
C *****************************************************************************
C     VPSC8DIM.FOR   -->   VERSION OF 2022/06/16
C
C     REPLACES VPSC8.DIM OF 2021/06/09 AND OTHER INCLUDE STATEMENTS DEFINING 
C     COMMON STATEMENTS & PARAMETERS TO BE INCLUDED DURING COMPILATION OF VPSC8
C *****************************************************************************

      MODULE PARAMETERS

      PARAMETER(NELMX =1)     ! MAXIMUM # OF ELEMENTS
      PARAMETER(NPHPEL=4)     ! MAXIMUM # OF PHASES PER ELEMENT
      PARAMETER(NPHMX =4)     ! MAXIMUM # OF PHASES OVER ALL ELEMENTS
      PARAMETER(NGRPEL=180000)   ! MAXIMUM # OF GRAINS PER ELEMENT (DRX)
      PARAMETER(NGRMX =180000)   ! MAXIMUM # OF GRAINS OVER ALL PHASES & ELEMS
      PARAMETER(NMODMX=8)     ! MAXIMUM # OF ACTIVE SL+TW MODES IN ANY PHASE
      PARAMETER(NTWMMX=3)     ! MAXIMUM # OF ACTIVE TWIN  MODES IN ANY PHASE
      PARAMETER(NSYSMX=48)    ! MAXIMUM # OF ACTIVE SL+TW SYSTEMS IN ANY PHASE
      PARAMETER(NTWSMX=12)    ! MAXIMUM # OF ACTIVE TWIN  SYSTEMS IN ANY PHASE
      PARAMETER(NNEIMX=7)     ! MAXIMUM # OF NEIGHBOURS
      PARAMETER(NDIFMX=20)    ! MAXIMUM # OF DIFFRACTING PLANES IN ANY PHASE
	  
      PARAMETER(MDIM=5)       ! MDIM=5 used for VP runs and MDIM=6 for EVP runs

      END MODULE PARAMETERS 
C *****************************************************************************  

      MODULE VPSC8DIM

      USE PARAMETERS

C --> ex COMMON /FILES/ /IOUNITS/ /MISCEL/
      CHARACTER*80  FILETEXT,FILECRYS,FILEAXES,FILEHIST,FILEDIFF
      CHARACTER*130 PROSA,LABEL(10)
      INTEGER  UR0,UR1,UR2,UR3,UR4,UR5,UR6,UR7,UW1,UW2,UW3,UW4,UW5
      INTEGER  STRAIN_CONTROL,FREEZE_SHAPE,JRAN
      REAL*8   MGSEC,MGTG,MBARTG,LBARTG,MASTGR,MASTPH,ImSINVPH
      REAL*8   LIJBARc,LIJGR,LIJPH,MU_MODE,NEFFGR,NEFFGRX
      LOGICAL  CONVERGED,initializeF
      REAL PI
      DIMENSION XID6(6,6),XID5(5,5),XID3(3,3)      ! identity matrix
      DIMENSION AUX5(5),AUX6(6),AUX33(3,3),AUX55(5,5),AUX66(6,6),
     #          AUX3333(3,3,3,3)

C --> ex COMMON/AVERAGE/
      REAL SVM,DVM,EPSACU,EPSVM,EIJINCR,EVMINCR
      DIMENSION SAV(6),SAVc(3,3),DAV(6),EAV(6),
     #       STDEVS(6),STDEVD(6),STDEVE(6),
     #       SASTAV(6),SASTBAR(6),DASTBAR(6), 
     #       SAVPH(6,NPHMX),SAVPHDEV(6,NPHMX),
     #       MBARTG(MDIM,MDIM),LBARTG(MDIM,MDIM),DBAR0(6),
     #       EPSTOTc(3,3),EPSTOTv(6)
	  
C --> ex COMMON/DATACRYST/
      REAL GAMD0
      DIMENSION DNCA(3,NSYSMX,NPHPEL),DBCA(3,NSYSMX,NPHPEL),
     #  SCHCA(5,NSYSMX,NPHPEL),HARD(NSYSMX,NSYSMX,NPHPEL),
     #  NRS(NSYSMX,NPHPEL),ISENSE(NSYSMX,NPHPEL),ITWTYPE(NSYSMX,NPHPEL)

C --> ex COMMON/DATAGRAIN/
      DIMENSION AG(:,:,:),WGT(:),NGR(0:NPHMX)
     #       ,AXISGR(:,:,:),FIJGR(:,:,:),LIJGR(:,:,:),
     #       ASGR(:,:,:,:,:),MASTGR(:,:,:)

C --> ex COMMON/DATAMODES/
      DIMENSION TWTHRES(2,NTWMMX,NPHPEL),TWSH(NTWMMX,NPHPEL),
     #       NMODES(NPHPEL),NSLMOD(NPHPEL),NTWMOD(NPHPEL),
     #       NSYST(NPHPEL),NTWSYS(NPHPEL),NSLSYS(NPHPEL),
     #       NSM(NMODMX,NPHPEL),ISECTW(NSYSMX,NPHPEL)
	 
C --> ex COMMON/DATAPHASE/
      INTEGER NPH,NELEM,IPHBOT,IPHTOP
      DIMENSION WPH(NPHMX),EULERPH(3,0:NPHMX),
     #       AXISPH(0:3,3,0:NPHMX),AVAX(3,NPHMX),SDAX(3,NPHMX),
     #       FIJPH(3,3,0:NPHMX),LIJPH(3,3,NPHMX),
     #       ASPH(3,3,3,3,NPHMX),MASTPH(MDIM,MDIM,NPHMX),
     #       ICRYSYMPH(NPHMX),
     #       ISHAPE(0:NPHMX),FREEZE_SHAPE(0:NPHMX),IFRAG(0:NPHMX),
     #       CRIT_SHP(NPHMX)

C --> ex COMMON/GRAINPROPS/
      DIMENSION CRSS(:,:),TAUE(:,:),
     #       SG(:,:),SG0(:,:),SGTRY(:,:),
     #       DG(:,:),DG0(:,:),
     #       DG_TRANSCc(3,3),DG_TRANSC(6),DG_TRANS(:,:),
     #       GTOTGR(:),GAMD0S(:,:),
     #       GAMDOT(:,:),SCH(:,:,:),
     #       MGSEC(:,:,:),MGTG(:,:,:),
     #       BG(:,:,:),NEFFGR(:)

C --> ex COMMON/POLAR_DEC/
      DIMENSION FIJ_PD(3,3),ROT_PD(3,3),ROT_PD_t(3,3),RCOROT(3,3)

C --> ex COMMON/RUN_COND/
      REAL ERR,ERRS,ERRD,ERRM,ERRSO
      INTEGER ITMAXEXT,ITMAXINT,ITMAXSO,IRECOVER,ISAVE,NPROC,
     #       IVGVAR,NWRITE,NWRITEX,INTERACTION,ICUBCOM,ICS,
     #       NUNIQUE,NRSMIN,IHARDLAW,IRATESENS,ITWINLAW,
     #       IUPDORI,IUPDSHP,IUPDHAR,IRSVAR,IFLU,
     #       jxrs,jxrsini,jxrsfin,jxrstep
      DIMENSION IDIFF(NPHMX),NDIFF(NPHMX)
      LOGICAL  RDCrun

C --> ex COMMON/RXL1/
       REAL TAYLAV,WORKAV,WRATE1,WRATE2
       DIMENSION SGVM(:),EPSVMGR(:),WORKGR(:),
     #       TAYLGR(:)

C --> ex COMMON/STAT_SHEARS/
      DIMENSION GAVMOD(NMODMX,NPHMX),GAVPH(NPHMX),GAVGR(:),
     #       ACTGR(:),ACTPH(NPHMX)

C --> ex COMMON/STATS_TWIN/
      DIMENSION EFTWFR(NTWMMX,NPHMX),TWFRPH(NTWMMX,NPHMX),
     #       TWFRSY(:,:),KTWSMX(:),NTWEVENTS(:),
     #       KTWMMX(:),PRITW(NPHMX),SECTW(NPHMX),
     #       SCH_STAT(-5:15,NTWMMX,NPHMX),VAR_STAT(6,NTWMMX,NPHMX),
     #       VAR_VF(6,NTWMMX,NPHMX),ITWVFT(:,:)
      DIMENSION TWFRGR(:,:),IPTSGR(:),KPTMGR(:)

C --> ex COMMON/TEST_COND/
      INTEGER ICTRL,NSTEPS
      REAL DBARCTRL,TIME,TIME_INCR,TEMP,TEMP_INCR,TEMP_INI,TEMP_FIN
      DIMENSION LIJBARc(3,3),DBAR(6),DBARv(6),DBARc(3,3),
     #       SBAR(6),SBARv(6),SBARc(3,3), WBARc(3,3),
     #       ILBAR(3,3),ISBARv(6),IDBARv(6),IDBARb(6),ISBARb(6)

C --> ex COMMON/NEIGHB/
      INTEGER NNEIGH
      DIMENSION WNEIGH(:,:),NEIGH(:,:)

C --> ex COMMON/THERMOELAST/
      DIMENSION MU_MODE(NMODMX,NPHPEL),
     #      CELCCv(6,6,NPHPEL),CELCC(6,6,NPHPEL),CELCS(:,:,:),       
     #      SELCCv(6,6,NPHPEL),SELCC(6,6,NPHPEL),SELCS(:,:,:), 
     #      ATHCCv(6,NPHPEL)  ,ATHCC(6,NPHPEL)  ,ATHCS(:,:), 
     #      CELAV(6,6),CELAVv(6,6) ,SELAV(6,6),SELAVv(6,6),
     #      ATHAV(6),ATHAVv(6),ATHAVc(3,3),
     #      EELAV(6),EELAVc(3,3),EELAVv(6),
     #      ETHAV(6),ETHAVv(6),ETHAVc(3,3),EELGR(:,:),
     #      AG_EL(:,:,:),BG_EL(:,:,:),TG_EL(:,:),
     #      DG_EL(:,:)

C --> ex COMMON/FLUCT/
      REAL SDSEQINTER,SDDEQINTER,SDSEQINTRA,SDDEQINTRA,UTILDE
      DIMENSION SEQ2(:),DEQ2(:),
     #      ASO(:,:),ESO(:,:),
     #      SECMOM5(:,:,:),SECMOM5D(:,:,:),
     #      BETFLU(:,:,:),CHIFLU(:,:,:),
     #      FSPH(5,5,NPHMX),ImSinvPH(5,5,NPHMX)      ! used inside GET_THEFLU


C     Grain storage is allocated to avoid the Windows PE image limit.
      ALLOCATABLE :: AG
      ALLOCATABLE :: WGT
      ALLOCATABLE :: AXISGR
      ALLOCATABLE :: FIJGR
      ALLOCATABLE :: LIJGR
      ALLOCATABLE :: ASGR
      ALLOCATABLE :: MASTGR
      ALLOCATABLE :: CRSS
      ALLOCATABLE :: TAUE
      ALLOCATABLE :: SG
      ALLOCATABLE :: SG0
      ALLOCATABLE :: SGTRY
      ALLOCATABLE :: DG
      ALLOCATABLE :: DG0
      ALLOCATABLE :: DG_TRANS
      ALLOCATABLE :: GTOTGR
      ALLOCATABLE :: GAMD0S
      ALLOCATABLE :: GAMDOT
      ALLOCATABLE :: SCH
      ALLOCATABLE :: MGSEC
      ALLOCATABLE :: MGTG
      ALLOCATABLE :: BG
      ALLOCATABLE :: NEFFGR
      ALLOCATABLE :: SGVM
      ALLOCATABLE :: EPSVMGR
      ALLOCATABLE :: WORKGR
      ALLOCATABLE :: TAYLGR
      ALLOCATABLE :: GAVGR
      ALLOCATABLE :: ACTGR
      ALLOCATABLE :: TWFRSY
      ALLOCATABLE :: KTWSMX
      ALLOCATABLE :: NTWEVENTS
      ALLOCATABLE :: KTWMMX
      ALLOCATABLE :: ITWVFT
      ALLOCATABLE :: TWFRGR
      ALLOCATABLE :: IPTSGR
      ALLOCATABLE :: KPTMGR
      ALLOCATABLE :: WNEIGH
      ALLOCATABLE :: NEIGH
      ALLOCATABLE :: CELCS
      ALLOCATABLE :: SELCS
      ALLOCATABLE :: ATHCS
      ALLOCATABLE :: EELGR
      ALLOCATABLE :: AG_EL
      ALLOCATABLE :: BG_EL
      ALLOCATABLE :: TG_EL
      ALLOCATABLE :: DG_EL
      ALLOCATABLE :: SEQ2
      ALLOCATABLE :: DEQ2
      ALLOCATABLE :: ASO
      ALLOCATABLE :: ESO
      ALLOCATABLE :: SECMOM5
      ALLOCATABLE :: SECMOM5D
      ALLOCATABLE :: BETFLU
      ALLOCATABLE :: CHIFLU
      END MODULE VPSC8DIM
C *****************************************************************************

      MODULE C4GAFLU

C --> ex COMMON/FLUCT2/ stifness in ellipsoid axes for passing to GET_THE_FLU
      DIMENSION C4GA_FLU(3,3,3,3)

      END MODULE C4GAFLU
C *****************************************************************************

      MODULE CHANGE_BASIS

      DIMENSION B(3,3,6)

      END MODULE CHANGE_BASIS
C *****************************************************************************

      MODULE CRSS_VOCE

      USE PARAMETERS

      DIMENSION TAU(NSYSMX,0:1,NPHPEL),THET(NSYSMX,0:1,NPHPEL),
     #            HPFAC(NSYSMX,NPHPEL),GRSZE(NPHPEL)

      END MODULE CRSS_VOCE
C *****************************************************************************

      MODULE CRYSTAL_SYM

      INTEGER   nsymop
      DIMENSION h(3,3,24),cvec(3,3)

      END MODULE CRYSTAL_SYM 
C *****************************************************************************

      MODULE CUB_COMP

      USE PARAMETERS

      PARAMETER (NIDORMX=300)
      PARAMETER (NIDMODMX=6)
      CHARACTER*3 IDLABEL
C --> ex COMMON/CUCO/ & /CUB/
      DIMENSION widmod(0:6),igrtype(:),idlabel(0:6)
      DIMENSION aidort(3,3,NIDORMX),itype(NIDORMX),normod(NIDMODMX)
      INTEGER NIDMOD,NOR
	  

C     Grain storage is allocated to avoid the Windows PE image limit.
      ALLOCATABLE :: igrtype
      END MODULE CUB_COMP
C *****************************************************************************

      MODULE DIFFRACT

      USE PARAMETERS

C --> ex COMMON/DIFFRACT/ 
      REAL ANGDETECTOR
      DIMENSION DIFF_VS(3,NDIFMX,NPHMX), DIFF_VC(3,24,NDIFMX,NPHMX),
     #      NFAMILY(NDIFMX,NPHMX), WGTSETPH(NDIFMX,NPHMX),
     #      IGRSET(:,:,:), NGRSETPH(NDIFMX,NPHMX),
     #      PARA_W(NDIFMX), EPS_W(NDIFMX,NPHMX)


C     Grain storage is allocated to avoid the Windows PE image limit.
      ALLOCATABLE :: IGRSET
      END MODULE DIFFRACT
C *****************************************************************************

      MODULE DISDEN

      USE PARAMETERS

C --> ex COMMON/DisDen/
       REAL SH_MOD,BOLTZ,GRSZE,CHIfor,CHIsub,K1gener,K2recom,
     #      RHO_initW,TEMPref,EDOTref
       DIMENSION BURG(NMODMX),ACTENER(NMODMX),EDOT0(NMODMX),
     # K1gener(NMODMX),K2recom(NMODMX),DRAG(NMODMX),HPcoef(NMODMX),
     # ATAU(NMODMX),Atau0(NMODMX),Atau1(NMODMX),Atau2(NMODMX),
     # FDEB(NMODMX),Fdeb0(NMODMX),Fdeb1(NMODMX),Fdeb2(NMODMX),

     # TAU_NUCL(NMODMX),T0nucl(NMODMX),T1nucl(NMODMX),T2nucl(NMODMX),
     # TAU_PROP(NMODMX),T0prop(NMODMX),T1prop(NMODMX),T2prop(NMODMX),
     # CTWSL(NTWMMX,NMODMX),CTS1(NTWMMX,NMODMX),CTS2(NTWMMX,NMODMX),

     # RHO_initS(NMODMX),RHOsat(NMODMX),
     # DELRHOS(NSYSMX),DELRHOW(NSYSMX),DELGAM(NSYSMX),
     # RHO_S(:,:),RHO_DEB(:),RHO_FOR(:),
     # RHO_FOR_PH(NMODMX,NPHMX),RHO_DEB_PH(NPHMX),RHO_TOT_PH(NPHMX)


C     Grain storage is allocated to avoid the Windows PE image limit.
      ALLOCATABLE :: RHO_S
      ALLOCATABLE :: RHO_DEB
      ALLOCATABLE :: RHO_FOR
      END MODULE DISDEN
C *****************************************************************************

      MODULE DISDEN_REV

      USE PARAMETERS

C --> ex COMMON/DisDen_rev/
      integer shear_tag
      real xD, qBS, mREV
      dimension xMu(NMODMX),xB(NMODMX),xBs(nsysmx)
     #   ,tau0(NMODMX),xf_self(NMODMX),xf_lat(NMODMX),xK(NMODMX)
     #   ,xlambda(:),xp(:,:),fBSm(nmodmx),fBSs(nsysmx)
     #   ,xRho_min(NMODMX),xRho_max(NMODMX),xRho_infl(NMODMX)
     #   ,rho_f(:,:),rho_r(:,:),taud(nsysmx)
     #   ,rho_prev(:,:),rho_deb(:,:),tauini(NSYSMX)
     #   ,shear_tag(:,:),rhos(nmodmx),rho_mod(nmodmx)
     #   ,for_rhos(nmodmx),rev_rhos(nmodmx) 
     #   ,for_rhosl(nmodmx),rev_rhosl(nmodmx)
     #   ,AMSS(NSYSMX,NSYSMX),deltaD(:,:),xp_g(:) 
     #   ,drho_tot(:,:),debs(:),rho_ps(:,:)   
     #   ,cts(NSYSMX,NSYSMX),tautwint(nsysmx)                      


C     Grain storage is allocated to avoid the Windows PE image limit.
      ALLOCATABLE :: xlambda
      ALLOCATABLE :: xp
      ALLOCATABLE :: rho_f
      ALLOCATABLE :: rho_r
      ALLOCATABLE :: rho_prev
      ALLOCATABLE :: rho_deb
      ALLOCATABLE :: shear_tag
      ALLOCATABLE :: deltaD
      ALLOCATABLE :: xp_g
      ALLOCATABLE :: drho_tot
      ALLOCATABLE :: debs
      ALLOCATABLE :: rho_ps
      END MODULE DISDEN_REV
C *****************************************************************************

      MODULE ESH123

      PARAMETER (ngaumx=40,ngaumx2=1600)

C --> ex COMMON/ESHELBY 1 2 3/integration points, weights, aspect ratio limits

      dimension ngaussph(12),ngaussth(12)
      dimension alpha(12,3,ngaumx2),aa1(12,6,ngaumx2),
     #          aww(12,3,ngaumx2),aaww1(12,6,ngaumx2),ww(12,ngaumx2)
      dimension punti(10,11),pesi(10,11),dte(0:10),
     #          puntigl(ngaumx),pesigl(ngaumx)

      END MODULE ESH123
C *****************************************************************************

      MODULE EXTRAPSO

      USE PARAMETERS

C --> ex COMMON/EXTRAPSO/
      dimension ASOOLD(:,:,:),ESOOLD(:,:,:),
     #          XMRSOLD(3)    


C     Grain storage is allocated to avoid the Windows PE image limit.
      ALLOCATABLE :: ASOOLD
      ALLOCATABLE :: ESOOLD
      END MODULE EXTRAPSO
C *****************************************************************************

      MODULE IRR_VARS

      USE PARAMETERS

C --> ex COMMON/IRR_VARS/

      REAL dpa_dose,dpa_rate,f_recomb,f_cl,
     & xNmax_a,xNmax_c,dose_max_a,dose_crit_c,dose_max_c,A_,
     & bmag_a,bmag_c,t_count,B_,rho_ref

      DIMENSION rho_irr(:,:),xN_i(:,:),
     & xr_i(:,:),xN_v(:,:),xr_v(:,:),
     & rho0(18),EpsGlobal(3,3),rad_mean(nphmx)


C     Grain storage is allocated to avoid the Windows PE image limit.
      ALLOCATABLE :: rho_irr
      ALLOCATABLE :: xN_i
      ALLOCATABLE :: xr_i
      ALLOCATABLE :: xN_v
      ALLOCATABLE :: xr_v
      END MODULE IRR_VARS
C *****************************************************************************

      MODULE MTSSAVE

      real KOVB3,MU0,D0,T0,TAUa,TAUi,TAUeini,TH0,KAP,MU
      real G0i,ED0i,QQi,PPi,G0e,ED0e,QQe,PPe,G0esat,EDesat0,TAUesat0
      real PSI,RHO,Cp1,Cp2,Cp3,Cp

      END MODULE MTSSAVE
C *****************************************************************************

      MODULE NEFFSAVE

      REAL RDCmin,RDCmax
      INTEGER KGmin,KGmax
	  
      END MODULE NEFFSAVE
C *****************************************************************************

      MODULE NLNR

      USE PARAMETERS

C --> ex COMMON/NLNR/

      INTEGER NSYSTX
      DIMENSION F0NR(5),XMASTX(5,5),SCHX(5,NSYSMX),TAUX(NSYSMX),
     #             GAMD0X(NSYSMX),NRSX(NSYSMX),ISENSEX(NSYSMX)

      END MODULE NLNR
C *****************************************************************************

      MODULE PCYSSAVE

      PARAMETER(NPROBE=72)

      DIMENSION STRAINR(5,NPROBE),STRESS(5,NPROBE),
     #          DPCYS(2,-2:NPROBE+3),SPCYS(2,-2:NPROBE+3)
      INTEGER NPART,NPROBEX
      REAL    WRATEREF,STNORM

      END MODULE PCYSSAVE
C *****************************************************************************

      MODULE PCYSTABLE

      PARAMETER(NPRBMX=20000)

      DIMENSION ACTION(0:12,0:6,0:6,0:6,5),
     #          REACTION(0:12,0:6,0:6,0:6,5),
     #          THGRID(4,NPRBMX),ITHGRID(4,NPRBMX)
      INTEGER ITABLE,NPROBE,NPROBERED

      END MODULE PCYSTABLE
C *****************************************************************************

