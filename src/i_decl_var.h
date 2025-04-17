! Declaration file for all the variables used in the program
! Structure
! Group
!   - site; species; climate; parameters
! Inputs are ordered by the order they are provided for the program
! Other variables are ordered by group and alphabeticaly

!***************************************
! INPUT

! Global configuration
!integer, parameter :: rk = selected_real_kind(8)      ! 8 byte real

! Site data ----------------------------
real(kind=8) :: lat                             ! site latitude
integer :: soil_class                           ! soil parameters for soil class
real(kind=8) :: asw                             ! available soil water
real(kind=8) :: asw_max                         ! maximum available soil water
real(kind=8) :: asw_min                         ! minimum available soil water
real(kind=8) :: psi_soil                        ! soil water potential
integer :: year_i                               ! initial year when the simulations starts
integer :: month_i                              ! initial month when the simulation starts
real(kind=8) :: altitude                        ! altitude of the site location, m
real(kind=8) :: soil_carbon_i                   ! initial soil carbon
real(kind=8) :: soil_nitrogen_i                 ! initial soil nitrogen
real(kind=8) :: soil_phosphorus_i               ! initial soil phosphorus
real(kind=8) :: bb_cbp			                ! bark beetle spatial scaling
real(kind=8) :: bb_cst			                ! bark beetle temporal scaling

! Species data -------------------------
integer, dimension(n_sp) :: year_p              ! year when species was planted
integer, dimension(n_sp) :: month_p             ! month when species was planted
real(kind=8), dimension(n_sp) :: fertility      ! initial site fertility rating for a given species
real(kind=8), dimension(n_sp) :: biom_foliage_i ! initial foliage biomass for a species
real(kind=8), dimension(n_sp) :: biom_root_i    ! initial root biomass for a species
real(kind=8), dimension(n_sp) :: biom_stem_i    ! initial stem biomass for a species
real(kind=8), dimension(n_sp) :: stems_n_i      ! initial stand stocking for a species
real(kind=8), dimension(n_sp) :: n_sapling      ! initial density after planting/regenerating
real(kind=8), dimension(n_sp) :: sapling_wf     ! foliage biomass of saplings (g/sapling)
real(kind=8), dimension(n_sp) :: sapling_wr     ! root biomass of saplings (g/sapling)
real(kind=8), dimension(n_sp) :: sapling_ws     ! stem biomass of saplings (g/sapling)
real(kind=8), dimension(n_sp) :: deadwood_carbon_i               ! initial deadwood carbon
real(kind=8), dimension(n_sp) :: deadwood_nitrogen_i             ! initial deadwood nitrogen
real(kind=8), dimension(n_sp) :: litter_carbon_i                 ! initial litter carbon
real(kind=8), dimension(n_sp) :: litter_nitrogen_i               ! initial litter nitrogen
real(kind=8), dimension(n_sp) :: n_sapling_i    ! number of seedlings/saplings
real(kind=8), dimension(n_sp) :: sapling_wf_i   ! foliage mass of saplings (g DW)
real(kind=8), dimension(n_sp) :: sapling_wr_i   ! root mass of saplings (g DW)
real(kind=8), dimension(n_sp) :: sapling_ws_i   ! stem mass of saplings (g DW)
real(kind=c_double), dimension(n_m, n_sp) :: dynamic_fertility	 ! downregulate CO2 fertilization

! Climate ------------------------------
real(kind=8), dimension(n_m) :: tmp_min         ! minimum monthly temperature
real(kind=8), dimension(n_m) :: tmp_max         ! maximum monthly temperature
real(kind=8), dimension(n_m) :: tmp_ave			! mean monthly temperature
real(kind=8) :: tmp_year, prcp_year				! yearly average temperature
real(kind=8), dimension(12) :: tmp_yasso_year	! yearly average temperature
real(kind=8), dimension(n_m) :: prcp            ! monthly precipitation sum
real(kind=8), dimension(n_m) :: rh	            ! monthly relative humidity
real(kind=8), dimension(n_m) :: snow_prcp       ! monthly snow precipitation sum
!real(kind=8), dimension(n_m) :: snowmelt_f     ! potential monthly snowmelt
real(kind=8), dimension(n_m) :: solar_rad       ! mean daily incident solar radiation
real(kind=8), dimension(n_m) :: frost_days      ! number of frost days per month
real(kind=8), dimension(n_m) :: vpd_day
real(kind=8), dimension(n_m) :: co2             ! atmospheric CO2
real(kind=8), dimension(n_m) :: d13catm         ! added d13C of atmospheric CO2 (per mil)
real(kind=8), dimension(n_m) :: n_deposition    ! nitrogen deposition (t/ha)
real(kind=8), dimension(n_m) :: weibull_a    	! scale of wind distribution 
real(kind=8), dimension(n_m) :: weibull_k    	! shape of wind distribution 
real(kind=8), dimension(n_m) :: thinn_adjust, lag, previous_density, current_density    ! thinning efect on critical wind speed

! Parameters ---------------------------

! Biomass partitioning and turnover
real(kind=8), dimension(n_sp) :: pFS2           ! Foliage:stem partitioning ratio at D=2 cm
real(kind=8), dimension(n_sp) :: pFS20          ! Foliage:stem partitioning ratio at D=20 cm
real(kind=8), dimension(n_sp) :: aWs            ! Constant in the stem mass v. diam. relationship
real(kind=8), dimension(n_sp) :: nWs            ! Power in the stem mass v. diam. relationship
real(kind=8), dimension(n_sp) :: pRn            ! Minimum fraction of NPP to roots
real(kind=8), dimension(n_sp) :: pRx            ! Maximum fraction of NPP to roots
real(kind=8), dimension(n_sp) :: gammaF1        ! Coefficients in monthly litterfall rate
real(kind=8), dimension(n_sp) :: gammaF0        ! Coefficients in monthly litterfall rate
real(kind=8), dimension(n_sp) :: tgammaF        ! Coefficients in monthly litterfall rate
real(kind=8), dimension(n_sp) :: gammaR         ! Average monthly root turnover rate
integer, dimension(n_sp) :: leafgrow            ! If deciduous, leaves are produced at end of this month
integer, dimension(n_sp) :: leaffall            ! If deciduous, leaves all fall at start of this month
integer, dimension(n_sp) :: leafgrow_base            ! If deciduous, leaves are produced at end of this month
integer, dimension(n_sp) :: leaffall_base           ! If deciduous, leaves all fall at start of this month

! NPP & conductance modifiers
real(kind=8), dimension(n_sp) :: Tmin           ! Minimum temperature for growth
real(kind=8), dimension(n_sp) :: Topt           ! Optimum temperature for growth
real(kind=8), dimension(n_sp) :: Tmax           ! Maximum temperature for growth
real(kind=8), dimension(n_sp) :: kF             ! Days production lost per frost day
real(kind=8), dimension(n_sp) :: SWconst0       ! Moisture ratio deficit for fq = 0.5
real(kind=8), dimension(n_sp) :: SWpower0       ! Power of moisture ratio deficit
real(kind=8), dimension(n_sp) :: fCalpha700     ! Assimilation enhancement factor at 700 ppm
real(kind=8), dimension(n_sp) :: fCg700         ! Canopy conductance enhancement factor at 700 ppm
real(kind=8), dimension(n_sp) :: m0             ! Value of 'm' when FR = 0
real(kind=8), dimension(n_sp) :: fN0            ! Value of 'fNutr' when FR = 0
real(kind=8), dimension(n_sp) :: fNn            ! Power of (1-FR) in 'fNutr'
real(kind=8), dimension(n_sp) :: MaxAge         ! Maximum stand age used in age modifier
real(kind=8), dimension(n_sp) :: nAge           ! Power of relative age in function for f_age
real(kind=8), dimension(n_sp) :: rAge           ! Relative age to give f_age = 0.5

! Stem mortality & self-thinning
real(kind=8), dimension(n_sp) :: gammaN1        ! Mortality rate for large t
real(kind=8), dimension(n_sp) :: gammaN0        ! Seedling mortality rate (t = 0)
real(kind=8), dimension(n_sp) :: tgammaN        ! Age at which mortality rate has median value
real(kind=8), dimension(n_sp) :: ngammaN        ! Shape of mortality response
real(kind=8), dimension(n_sp) :: wSx1000        ! Max. stem mass per tree @ 1000 trees/hectare
real(kind=8), dimension(n_sp) :: thinPower      ! Power in self-thinning rule
real(kind=8), dimension(n_sp) :: mF             ! Fraction mean single-tree foliage biomass lost per dead tree
real(kind=8), dimension(n_sp) :: mR             ! Fraction mean single-tree root biomass lost per dead tree
real(kind=8), dimension(n_sp) :: mS             ! Fraction mean single-tree stem biomass lost per dead tree
real(kind=8) :: beta_repar             			! Beta parameter in Jandl et al. (2020)
real(kind=8), dimension(n_sp) :: dead_volume    ! Volume of dead trees
integer :: nyear, dist_start			  		! Number of years to compute the growth effciiency parameter, year of disturbance activation
real(kind=8), dimension(n_sp) :: greff_par		! Growth efficiency for bark beetle module   
real(kind=8) :: greff_spruce					! Growth efficiency for bark beetle module   


! Canopy structure and processes
real(kind=8), dimension(n_sp) :: SLA0           ! Specific leaf area at age 0
real(kind=8), dimension(n_sp) :: SLA1           ! Specific leaf area for mature leaves
real(kind=8), dimension(n_sp) :: tSLA           ! Age at which specific leaf area = (SLA0+SLA1)/2
real(kind=8), dimension(n_sp) :: k              ! Extinction coefficient for absorption of PAR by canopy
real(kind=8), dimension(n_sp) :: fullCanAge     ! Age at canopy closure
real(kind=8), dimension(n_sp) :: MaxIntcptn     ! Maximum proportion of rainfall evaporated from canopy
real(kind=8), dimension(n_sp) :: LAImaxIntcptn  ! LAI for maximum rainfall interception
real(kind=8), dimension(n_sp) :: cVPD           ! 'DF LAI for 50% reduction of VPD in canopy
real(kind=8), dimension(n_sp) :: alphaCx        ! Canopy quantum efficiency
real(kind=8), dimension(n_sp) :: y              ! Ratio NPP/GPP
real(kind=8), dimension(n_sp) :: MinCond        ! Minimum canopy conductance
real(kind=8), dimension(n_sp) :: MaxCond        ! Maximum canopy conductance
real(kind=8), dimension(n_sp) :: LAIgcx         ! LAI for maximum canopy conductance
real(kind=8), dimension(n_sp) :: CoeffCond      ! Defines stomatal response to VPD
real(kind=8), dimension(n_sp) :: BLcond         ! Canopy boundary layer conductance
real(kind=8), dimension(n_sp) :: RGcGW          ! The ratio of diffusivities of CO2 and water vapour in air
real(kind=8), dimension(n_sp) :: D13CTissueDif  ! d13C difference of modelled tissue and new photosynthate
real(kind=8), dimension(n_sp) :: aFracDiffu     ! Fractionation against 13C in diffusion
real(kind=8), dimension(n_sp) :: bFracRubi      ! Enzymatic fractionation by Rubisco

! Maintenance respiration
real(kind=8), dimension(n_sp) :: pssN1          ! Share of sapwood for large t
real(kind=8), dimension(n_sp) :: pssN0          ! Share of sapwood for seedlings
real(kind=8), dimension(n_sp) :: tpssN          ! Age in which share of sapwood has median value
real(kind=8), dimension(n_sp) :: npssN          ! Shape response of share of sapwood
real(kind=8), dimension(n_sp) :: fr_ratio       ! Ratio of fine root to foliage biomass
real(kind=8), dimension(n_sp) :: fr_wn          ! Modifier of fine root to foliage ratio

! Wood and stand properties
real(kind=8), dimension(n_sp) :: fracBB0        ! Branch and bark fraction at age 0
real(kind=8), dimension(n_sp) :: fracBB1        ! Branch and bark fraction for mature stands
real(kind=8), dimension(n_sp) :: tBB            ! Age at which fracBB = (fracBB0+fracBB1)/2
real(kind=8), dimension(n_sp) :: rho0           ! Minimum basic density - for young trees
real(kind=8), dimension(n_sp) :: rho1           ! Maximum basic density - for older trees
real(kind=8), dimension(n_sp) :: tRho           ! Age at which rho = (rhoMin+rhoMax)/2
real(kind=8), dimension(n_sp) :: stem_bark_corr ! correction of stem bark biomass for soil models
integer, dimension(n_sp) :: CrownShape          !***DF crown shape of a given species; 1=cone, 2=ellipsoid, 3=half-ellipsoid, 4=rectangular

! Height and Wolume
real(kind=8), dimension(n_sp) :: aH, nHB, nHC
real(kind=8), dimension(n_sp) :: aV, nVB, nVH, nVBH
real(kind=8), dimension(n_sp) :: aK, nKB, nKH, nKC, nKrh
real(kind=8), dimension(n_sp) :: aHL, nHLB, nHLL, nHLC, nHLrh

! Delta 13
real(kind=8), dimension(n_sp) :: Qa, Qb
real(kind=8), dimension(n_sp) :: gDM_mol, molPAR_MJ

! Bias correction
real(kind=8), dimension(n_sp) :: Dscale0, DscaleB, Dscalerh, Dscalet, DscaleC
real(kind=8), dimension(n_sp) :: Dshape0, DshapeB, Dshaperh, Dshapet, DshapeC
real(kind=8), dimension(n_sp) :: Dlocation0, DlocationB, Dlocationrh, Dlocationt, DlocationC
real(kind=8), dimension(n_sp) :: wsscale0, wsscaleB, wsscalerh, wsscalet, wsscaleC
real(kind=8), dimension(n_sp) :: wsshape0, wsshapeB, wsshaperh, wsshapet, wsshapeC
real(kind=8), dimension(n_sp) :: wslocation0, wslocationB, wslocationrh, wslocationt, wslocationC

! ICBM/2N soil model
real(kind=8), dimension(n_sp) :: klmax, krmax, komax, hc, FR_mod
real(kind=8), dimension(n_sp) :: qir, qil, qh, qbc
real(kind=8), dimension(n_sp) :: el, er
real(kind=8), dimension(n_sp) :: Ncf, Ncr, Ncs
real(kind=8), dimension(n_sp) :: Yl_Coutflux, Yr_Coutflux, Litter_Coutflux, Litter_C, Yr_C, Yl_C, &
								 Yb_C, Yb_input, Yr_input, Yl_input, fr_cr_ratio, &
								 Yb_input_aux, Yr_input_aux, Yl_input_aux
real(kind=8), dimension(n_sp) :: Yl_Noutflux, Yr_Noutflux, Litter_Noutflux, Yr_N, Yl_N
real(kind=8), dimension(n_sp) :: humification_Litter, humification_l, humification_r, humification_N_l, humification_N_r 
real(kind=8), dimension(n_sp) :: Un, Up
real(kind=8) :: O_Coutflux, O_Noutflux, O_C, O_N, TotalCarbo, TotalNitro 
real(kind=8) :: FR, FRin, Nav
real(kind=8) :: kl, kr, ko
real(kind=8) :: sand_prop, leaching_org, leaching_min

! Yasso soil model
integer:: step_aux
real(kind=8),dimension(35) :: parameters_yasso ! parameters
real(kind=8),dimension(15,n_sp) :: parameters_AWEN ! parameters
real (kind=8), dimension(5) :: folAWENH, stAWENH, crbAWENH, huAWENH
real (kind=8), dimension(n_m,n_sp,4,5) :: soilC
real (kind=8) :: litter_input, tstep, size_eff, leaching
real(kind=8), dimension(5), parameter :: iniAWEN = (/0.3d0, 0.03d0, 0.05d0, 0.62d0, 0.d0/)

! Soil temperature
real (kind=8) :: t_soil, t_soil_surf

! Nutrient distribution
real (kind=8), dimension(n_sp) :: root_dist		! root distribution parameter
real (kind=8), dimension(n_sp) :: SRA			! specific root area
real (kind=8), dimension(n_sp) :: root_area		! total root area
real (kind=8) :: n_dist							! nitorgen distribution parameter

!***************************************
! DERIVED VARIABLES

! Helpers ------------------------------
integer :: i = 1                                ! indexing for species
integer :: ii = 1                               ! indexing for month (row of climatic data)
integer :: month = 1
integer :: b_n = 2                              ! how many times to iterate for biass correction
integer :: n = 1                                ! count for bias correction
logical :: b_cor = .TRUE.                            ! if something has changed and wee need to correct bias
logical, dimension(n_sp) :: final_harvest		! If current management operation is final harvest

! Climatic variables -------------
real(kind=8), dimension(12) :: adjSolarZenithAngle
real(kind=8), dimension(12) :: day_length


! Stand variables ----------------
real(kind=8), dimension(n_m, n_sp) :: age     ! Age of each species and month
real(kind=8), dimension(n_m, n_sp) :: age_m   ! Age of each species used for calculating modifiers (one month less than s_age)
real(kind=8), dimension(n_sp) :: stems_n
real(kind=8), dimension(n_sp) :: stems_n_ha     ! potential number of stems per ha

real(kind=8), dimension(n_sp) :: basal_area     ! stand level basal area
real(kind=8), dimension(n_sp) :: basal_area_prop    ! proportion of basal area
real(kind=8), dimension(n_sp) :: dbh, hbh            ! average tree DBH, cm
real(kind=8), dimension(n_sp) :: dcb, hcb            ! average tree diameter at crown base, cm
real(kind=8), dimension(n_sp) :: h_vol	        ! heartwood volume
real(kind=8), dimension(n_sp) :: height, biom_sapwood_ratio, biom_vol         ! average tree height, m
real(kind=8) :: Height_max
! Crop tree variables
integer, dimension(n_sp) :: crop_trees, crop_trees_aux
real(kind=8), dimension(n_sp) :: dbh_crop
real(kind=8), dimension(n_sp) :: vol_crop
real(kind=8), dimension(n_sp) :: basal_area_crop
real(kind=8), dimension(n_sp) :: biom_stem_crop
real(kind=8), dimension(n_sp) :: height_crop
real(kind=8), dimension(n_sp) :: dbh_harv
real(kind=8), dimension(n_sp) :: stem_harv
real(kind=8), dimension(n_sp) :: biom_stem_harv
real(kind=8), dimension(n_sp) :: height_harv

real(kind=8), dimension(n_sp) :: crown_length   !***DF mean live-crown length (m) of a species
real(kind=8), dimension(n_sp) :: crown_width    ! ***DF mean crown diameter (m)

real(kind=8), dimension(n_sp) :: volume
real(kind=8), dimension(n_sp) :: volume_mai
real(kind=8), dimension(n_sp) :: volume_old
real(kind=8), dimension(n_sp) :: volume_cum
real(kind=8), dimension(n_sp) :: volume_change
real(kind=8), dimension(n_sp) :: harvesting		! Harvested volume
real(kind=8), dimension(n_sp) :: salvage		! Salvaged volume
real(kind=8), dimension(n_sp) :: residues		! Harvesting residues - assumed to stay on site
real(kind=8), dimension(n_sp) :: root_residues		! Roots pool of harvested trees - assumed to stay on site
real(kind=8), dimension(n_sp) :: foliage_residues 	! Foliage of harvested trees - assumed to stay on site

real(kind=8), dimension(n_sp) :: competition_total

real(kind=8), dimension(n_m, n_sp) :: SLA       ! Specific leaf area
real(kind=8), dimension(n_m, n_sp) :: fracBB    ! Fraction of stem biomass as branch and bark
real(kind=8), dimension(n_m, n_sp) :: wood_density  ! Whole-tree basic density


! Canopy variables ---------------
real(kind=8), dimension(n_sp) :: LAI            ! Canopy LAI (mean annual LAI if output time step is annual, and final year LAI if step is whole rotation)
real(kind=8), dimension(n_sp) :: lai_total      ! total competition of the forest
real(kind=8), dimension(n_sp) :: elai           ! effective LAI Ding et al. (2014)
real(kind=8), dimension(n_sp) :: LAI_per        ! species specific proportion of lai
real(kind=8), dimension(n_sp) :: lai_above      ! leaf area above the given species
real(kind=8), dimension(n_sp) :: canopy_vol_frac
real(kind=8), dimension(n_sp) :: lai_sa_ratio !the ratio of mean tree leaf area (m2) to crownSA (m2)
integer, dimension(n_sp) :: layer_id


! Stocks variables ---------------
real(kind=8), dimension(n_sp) :: biom_foliage   ! Foliage biomass
real(kind=8), dimension(n_sp) :: biom_foliage_debt
real(kind=8), dimension(n_sp) :: biom_root      ! Root biomass
real(kind=8), dimension(n_sp) :: biom_stem      ! Stem biomass, including branches and bark
real(kind=8), dimension(n_sp) :: biom_tree      ! average tree stem mass
real(kind=8), dimension(n_sp) :: biom_tree_max  ! Max. mean tree stem mass at current stocking

real(kind=8), dimension(n_sp) :: biom_incr_foliage
real(kind=8), dimension(n_sp) :: biom_incr_root
real(kind=8), dimension(n_sp) :: biom_incr_stem

real(kind=8), dimension(n_sp) :: biom_loss_foliage  ! Litter fall
real(kind=8), dimension(n_sp) :: biom_loss_root

real(kind=8), dimension(n_sp) :: biom_sapwood, height_sapwood, biom_fineroot

! Modifiers ----------------------
real(kind=8), dimension(n_m, n_sp) :: f_age     ! Age related modifier
real(kind=8), dimension(n_m, n_sp) :: f_tmp     ! Temperature modifier
real(kind=8), dimension(n_m, n_sp) :: f_tmp_gc  ! gc canopy conductance modifier as in Feikema et al 2010 FEM 260,663–678
real(kind=8), dimension(n_m, n_sp) :: f_frost   ! Frost modifier
real(kind=8), dimension(n_m, n_sp) :: f_calpha  !
real(kind=8), dimension(n_m, n_sp) :: f_cg      !

real(kind=8), dimension(n_sp) :: f_vpd
real(kind=8), dimension(n_sp) :: f_sw
real(kind=8), dimension(n_sp) :: f_nutr
real(kind=8), dimension(n_sp) :: f_phys

real(kind=8), dimension(n_m, n_sp) :: f_inund, istress  !Inundation modifier 

! Production ---------------------
real(kind=8), dimension(n_sp) :: pfsConst       ! Derived from pFS2 and PFS20
real(kind=8), dimension(n_sp) :: pfsPower       ! Derived from pFS2 and PFS20

real(kind=8), dimension(n_sp) :: pFS
real(kind=8), dimension(n_sp) :: fi             !***DF the proportion of above canopy PAR absorbed by each species
real(kind=8), dimension(n_sp) :: lambda_h       !Constant to account for horizontal canopy heterogeneity such as gaps between trees and the change in zenith angle (and shading) with latitude and season (see Equations 2 and 5 of Forrester et al., 2014, Forest Ecosystems, 1:17)
real(kind=8), dimension(n_sp) :: lambda_v       !Constant to partition light between species and to account for vertical canopy heterogeneity (see Equations 2 and 3 of Forrester et al., 2014, Forest Ecosystems, 1:17)

real(kind=8), dimension(n_sp) :: npp_fract_root
real(kind=8), dimension(n_sp) :: npp_fract_stem
real(kind=8), dimension(n_sp) :: npp_fract_foliage

real(kind=8), dimension(n_sp) :: apar           ! RADint
real(kind=8), dimension(n_sp) :: aero_resist    ! # 'DF aerodynamic resistance within the canopy at the height of the given species (s m-1)
real(kind=8), dimension(n_sp) :: VPD_sp         ! # 'DF VPD around the crowns of the given species
real(kind=8), dimension(n_sp) :: alpha_c        ! Canopy quantum efficiency after modifiers
real(kind=8), dimension(n_sp) :: epsilon    ! Light-use efficiency based on GPP
real(kind=8), dimension(n_sp) :: epsilon_gpp    ! Light-use efficiency based on GPP
real(kind=8), dimension(n_sp) :: epsilon_npp    ! Light-use efficiency based on NPP
real(kind=8), dimension(n_sp) :: epsilon_biom_stem !Light-use efficiency based on stem biomass (increment in WS)
real(kind=8), dimension(n_sp) :: GPP
real(kind=8), dimension(n_sp) :: NPP
real(kind=8), dimension(n_sp) :: NPP_f          ! the full NPP before substraction of depth
real(kind=8), dimension(n_sp) :: gC
real(kind=8), dimension(n_sp) :: conduct_canopy
real(kind=8), dimension(n_sp) :: m
real(kind=8), dimension(n_sp) :: ra
real(kind=8), dimension(n_m,n_sp) :: rsap, icount


! Mortality ----------------------
real(kind=8), dimension(n_sp) :: mort_stress     ! Number of trees that died due to stress-related mortality
real(kind=8), dimension(n_sp) :: mort_thinn      ! Number of trees that died due to density-dependent mortality

real(kind=8), dimension(n_m, n_sp) :: gammaN
real(kind=8), dimension(n_m, n_sp) :: gammaF


integer, dimension(n_sp) :: t_n ! currnet thinnign number
real(kind=8), dimension(n_sp) :: mort_manag, dbh_removal, cwd_removal ! mortality due to management
real(kind=8) :: target_ba, target_vol, removal ! target basal area/volume


! Water use ----------------------
real(kind=8), dimension(n_sp) :: SWconst         ! soil parameters for soil class
real(kind=8), dimension(n_sp) :: SWpower         ! soil parameters for soil class
real(kind=8) :: Irrig
real(kind=8) :: poolFractn                       ! Determines fraction of excess water that remains on site
real(kind=8) :: water_runoff_polled              ! current stored runoff

real(kind=8) :: snow_water	                 ! current stored snow
real(kind=8) :: snowmelt	                 ! monthly snowmelt
real(kind=8) :: irrig_supl
real(kind=8) :: prcp_runoff
real(kind=8) :: excessSW

real(kind=8) :: conduct_soil
real(kind=8) :: minvtot, surfsw, zmin, grad_lim, acro_por, az, z1, z2, maxhei, wtd, wtp_d

! Transpiration
real(kind=8), dimension(n_sp) :: transp_veg       ! Traspiration from the forest
real(kind=8) :: evapotra_soil
real(kind=8) :: transp_total

real(kind=8), dimension(n_sp) :: prcp_interc_fract
real(kind=8), dimension(n_sp) :: prcp_interc

real(kind=8) :: prcp_interc_total                 ! total rain interception

real(kind=8) :: evapo_transp
real(kind=8) :: f_transp_scale                    !***DF scales GPP and NPP down if evapotranspiration is greater than ASW

real(kind=8), dimension(n_sp) :: WUE
real(kind=8), dimension(n_sp) :: WUE_transp       !***DF

real(kind=8), dimension(n_sp) :: inund_height   ! Inundation height tolerance 
real(kind=8), dimension(n_sp) :: inund_dur  ! Inundation duration tolerance 

! Hydrological parameters -------

real(kind=c_double), dimension(3,n_sp) :: par_plant	! Hydraulic traits
real(kind=c_double), dimension(2,n_sp) :: par_cost	! Cost parameters

! Wood Delta --------------------
real(kind=8), dimension(n_sp) :: fCalphax
real(kind=8), dimension(n_sp) :: fCg0

real(kind=8) :: air_pressure                    !Air pressure (kPa)
real(kind=8), dimension(n_sp) :: GPP_molsec     !GPP per second (mol/m2 s)
real(kind=8), dimension(n_sp) :: Gw_mol         !Canopy conductance for water vapor in mol/m2s
real(kind=8), dimension(n_sp) :: Gc_mol         !Canopy conductance for CO2 in mol/m2s
real(kind=8), dimension(n_sp) :: canopy_cover
real(kind=8), dimension(n_sp) :: InterCi        !intercellular CO2 concentration
real(kind=8), dimension(n_sp) :: D13CNewPS
real(kind=8), dimension(n_sp) :: D13CTissue


! Weibull -----------------------
real(kind=8), dimension(15, n_sp) :: bias_scale

! Wind disturbances -----------------------
real(kind=8), dimension(n_sp) :: Creg		! Tree pulling coefficient
real(kind=8), dimension(n_sp) :: MOE		! Modulus of elasticity
real(kind=8), dimension(n_sp) :: MOR		! Modulus of rupture
real(kind=8), dimension(n_sp) :: wet_fact	! Wet factor for wet stem weight calculation
real(kind=8), dimension(n_sp) :: crown_area_fact! Crown area factor
real(kind=8), dimension(n_sp) :: Cdrag, Ndrag
real(kind=8), dimension(n_sp) :: cws, cws10		!Critical wind speed
real(kind=8), dimension(n_sp) :: wind_damage	!Damage probability
real(kind=8), dimension(n_sp) :: wind_mortality	!Damage probability
real(kind=8), dimension(n_sp) :: biom_stem_wind, biom_root_wind, biom_foliage_wind, crown_width_aux, crown_length_aux	!Wind damage
real(kind=8) :: rnd_nr 						! Random number

! Bark beetle disturbances -----------------------
real(kind=8):: infestation_share		! Landsacape infestation index
real(kind=8):: background_infestation		! Background infestation index
real(kind=8):: bb_mort		! Landsacape infestation index
real(kind=8):: wind_lit, lit_dbh, stress_lit, thinn_lit		! Landsacape infestation index
real(kind=8), dimension(n_sp) :: c_st					! Temporal scaling factor
real(kind=8), dimension(n_sp) :: c_sb					! Spatial scaling factor
real(kind=8), dimension(n_sp) :: background_pressure	! Background beetle pressure
real(kind=8):: w_bio					!Biophysical weight on BB pressure index
real(kind=8):: w_wind					!Wind weight on BB pressure index
real(kind=8):: w_age					!Age weight on BB pressure index
real(kind=8):: w_host					!Host weight on BB pressure index
real(kind=8):: w_vita					!Tree vitality weight on BB pressure index
real(kind=8):: w_dens					!Tree density weight on BB pressure index
real(kind=8):: previous_damage, maxbiom_sp		! Previous damage parameters
real(kind=8), dimension(n_sp) :: biom_stem_bb, biom_root_bb, biom_foliage_bb	!Wind damage
real(kind=8):: wtmin					! Minimum winter temperature
integer:: spruce_idx									! Spruce indicator
integer:: d_type										! Drought susceptibility type
integer:: vita_type										! Vitality modifier type

! Salvage logging --------------------------------
real(kind=8), dimension(n_sp) :: salvage_wind 			! Share of wind damage salvaged
real(kind=8), dimension(n_sp) :: salvage_biotic 		! Share of bark beetle damage salvaged
real(kind=8), dimension(n_sp) :: salvage_fire 			! Share of fire damage salvaged

! Settings ----------------------
integer :: light_model                          ! 1 - 3PGpjs; 2 - 3PGmix
integer :: transp_model                         ! 1 - 3PGpjs; 2 - 3PGmix; 3 - Phydro
integer :: phys_model                           ! 1 - 3PGpjs; 2 - 3PGmix
integer :: height_model                         ! 1 - linear; 2-non-linear
integer :: correct_bias                         ! 0 - no; 1 - 3PGmix
integer :: calculate_d13c                       ! 0 - no; 1 - 3PGmix
integer :: soil_model	                        ! 0 - no; 1 - ICBM/2N; 3 - Yasso
integer :: gpp_model	                        ! 1 - 3PGpjs/mix; 2 - p-model; 3 - phydro
integer :: method_jmaxlim                       ! 1 - wang17; 2 - smith19: 3 - none
integer :: nut_lim                           	! 0 - no nutrient limitation; 1 - with nutrient limitation
integer :: maint_resp                           ! 1 - fixed NPP/GPP; 2 - calculate Rm
integer :: canopy_cond                          ! 1 - jarvis; 2 - p-model, 3 - phydro
integer :: mortality_model                      ! 1 - 3PGpjs; 2 - Brandl; 3 - LPJ-GUESS; 4 - iLand
integer :: wind_dist                         	! 0 - no disturbance; 1 - activate disturbance
integer :: beetle_dist                         	! 0 - no disturbance; 1 - activate disturbance
integer :: d_sus_type                         	! 1 - LandClim; 2 - growth efficiency based
