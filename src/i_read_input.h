! Reading the mapin input
! Structure
! Group
!   - site; species; climate; parameters
! Inputs are ordered by the order they are provided for the program
! Other variables are ordered by group and alphabeticaly


! Site data ----------------------------
lat           = siteInputs(1)
altitude      = int( siteInputs(2) )
soil_class    = int( siteInputs(3) )
aSW           = siteInputs(4)
asw_min       = siteInputs(5)
asw_max       = siteInputs(6)
year_i        = int( siteInputs(7) )
month_i       = int( siteInputs(8) )
soil_carbon_i   = siteInputs(9)
soil_nitrogen_i = siteInputs(10)
soil_phosphorus_i = siteInputs(11)
bb_cbp = siteInputs(12)
bb_cst = siteInputs(13)
background_infestation = siteInputs(14)

! Species data -------------------------
year_p      = int( speciesInputs(:,1) )
month_p     = int( speciesInputs(:,2) )
fertility   = speciesInputs(:,3)
stems_n_i   = speciesInputs(:,4)
biom_stem_i = speciesInputs(:,5)
biom_root_i = speciesInputs(:,6)
biom_foliage_i = speciesInputs(:,7)
deadwood_carbon_i   = speciesInputs(:,8)
deadwood_nitrogen_i = speciesInputs(:,9)
litter_carbon_i   = speciesInputs(:,10)
litter_nitrogen_i = speciesInputs(:,11)
n_sapling_i   = speciesInputs(:,12)
sapling_wf_i = speciesInputs(:,13)
sapling_wr_i = speciesInputs(:,14)
sapling_ws_i = speciesInputs(:,15)
cwd_removal = speciesInputs(:,16)
salvage_wind = speciesInputs(:,17)
salvage_biotic = speciesInputs(:,18)
salvage_fire = speciesInputs(:,19)

! Climate ------------------------------
tmp_min     = forcingInputs(:,1)
tmp_max     = forcingInputs(:,2)
tmp_ave     = forcingInputs(:,3)
prcp        = forcingInputs(:,4)
rh          = forcingInputs(:,5)
solar_rad   = forcingInputs(:,6)
frost_days  = forcingInputs(:,7)
vpd_day     = forcingInputs(:,8)
co2         = forcingInputs(:,9)
d13catm     = forcingInputs(:,10)
N_deposition    = forcingInputs(:,11)
weibull_a    = forcingInputs(:,12)
weibull_k    = forcingInputs(:,13)

! Settings ----------------------------
light_model = settings(1)
transp_model = settings(2)
phys_model = settings(3)
height_model = settings(4)
correct_bias = settings(5)
calculate_d13c = settings(6)
soil_model = settings(7)
gpp_model = settings(8)
method_jmaxlim = settings(9)
nut_lim = settings(10)
maint_resp = settings(11)
canopy_cond = settings(12)
mortality_model = settings(13)
wind_dist = settings(14)
beetle_dist = settings(15)
d_sus_type = settings(16)
dist_start = settings(17)
spruce_idx = settings(18)