! Initialize variables

! Output ------------------------------
output(:,:,:,:) = 0.d0

! Climatic variables -------------

! Stand variables ---------------
age(:,:) = 0.d0
age_m(:,:) = 0.d0
stems_n(:) = 0.d0
basal_area(:) = 0.d0
dbh(:) = 0.d0
height(:) = 0.d0
volume(:) = 0.d0
volume_mai(:) = 0.d0
volume_old(:) = 0.d0
volume_cum(:) = 0.d0
volume_change(:) = 0.d0
dbh_crop(:) = 0.d0
vol_crop(:) = 0.d0
basal_area_crop(:) = 0.d0
biom_stem_crop(:) = 0.d0
height_crop(:) = 0.d0
crop_trees(:) = 0
crop_trees_aux(:) = 0
dbh_harv(:) = 0.d0
stem_harv(:) = 0.d0
biom_stem_harv(:) = 0.d0
height_harv(:) = 0.d0

! Canopy variables ---------------
layer_id(:) = 0.d0
sla(:,:) = 0.d0
lai(:) = 0.d0
canopy_cover(:) = 0.d0
canopy_vol_frac(:) = 0.d0
lai_above(:) = 0.d0
lambda_v(:) = 0.d0
lambda_h(:) = 0.d0
aero_resist(:) = 0.d0
vpd_sp(:) = 0.d0
lai_sa_ratio(:) = 0.d0

! Stocks variables ---------------
biom_foliage(:) = 0.d0
biom_root(:) = 0.d0
biom_stem(:) = 0.d0
biom_sapwood(:) = 0.d0
stems_n(:) = 0.d0
biom_tree(:) = 0.d0
biom_loss_foliage(:) = 0.d0
biom_loss_root(:) = 0.d0
biom_incr_foliage(:) = 0.d0
biom_incr_root(:) = 0.d0
biom_incr_stem(:) = 0.d0
biom_foliage_debt(:) = 0.d0

! Modifiers ---------------
f_age(:,:) = 0.d0
f_vpd(:) = 0.d0
f_tmp(:,:) = 0.d0
f_tmp_gc(:,:) = 0.d0
f_frost(:,:) = 0.d0
f_sw(:) = 0.d0
f_nutr(:) = 0.d0
f_calpha(:,:) = 0.d0
f_cg(:,:) = 0.d0
f_phys(:) = 0.d0
f_transp_scale = 0.d0

! Production ---------------
gpp(:) = 0.d0
npp(:) = 0.d0
npp_f(:) = 0.d0
apar(:) = 0.d0
fi(:) = 0.d0
alpha_c(:) = 0.d0
epsilon(:) = 0.d0
epsilon_gpp(:) = 0.d0
epsilon_npp(:) = 0.d0
epsilon_biom_stem(:) = 0.d0
npp_fract_stem(:) = 0.d0
npp_fract_foliage(:) = 0.d0
npp_fract_root(:) = 0.d0
pFS(:) = 0.d0
ra(:) = 0.d0
harvesting(:) = 0.d0
salvage(:) = 0.d0

! Water use ---------------
conduct_canopy(:) = 0.d0
conduct_soil = 0.d0
evapotra_soil = 0.d0
prcp_interc(:) = 0.d0
prcp_interc_fract(:) = 0.d0
prcp_runoff = 0.d0
irrig_supl = 0.d0
wue(:) = 0.d0
wue_transp(:) = 0.d0
evapo_transp = 0.d0
transp_veg(:) = 0.d0

transp_total = 0.d0

! Mortality ---------------
biom_tree_max(:) = 0.d0
mort_stress(:) = 0.d0
mort_thinn(:) = 0.d0
mort_manag(:) = 0.d0
t_n(:) = 1
cws(:) = 0.d0
wind_damage(:) = 0.d0
wind_mortality(:) = 0.d0
w_bio = 1.d0 / 6.d0
w_age = 0.125d0
w_host = 0.125d0
w_vita = 0.125d0
w_dens = 0.125d0
w_wind = 1.d0 / 2.d0
previous_damage = 0.d0
biom_stem_bb(:) = 0.d0
biom_root_bb(:) = 0.d0
biom_foliage_bb(:) = 0.d0
biom_stem_wind(:) = 0.d0
biom_root_wind(:) = 0.d0
biom_foliage_wind(:) = 0.d0
thinn_adjust(:) = 1.d0
lag(:) = 10.d0
infestation_share=0.d0

! Wood Delta ------------------
Gc_mol(:) = 0.d0
Gw_mol(:) = 0.d0
D13CNewPS(:) = 0.d0
D13CTissue(:) = 0.d0
InterCi(:) = 0.d0 / 1000000.d0

! Weibull ---------------------
bias_scale(:,:) = 0.d0

! Soil ------------------
kr = 0.d0
kl = 0.d0
ko = 0.d0
Yl_C(:) = 0.d0
Yr_C(:) = 0.d0
Yl_Coutflux(:) = 0.d0
Yr_Coutflux(:) = 0.d0
Litter_Coutflux(:) = 0.d0
Litter_C(:) = 0.d0
humification_Litter(:) = 0.d0
humification_r(:) = 0.d0
humification_l(:) = 0.d0
humification_N_r(:) = 0.d0
humification_N_l(:) = 0.d0
O_Coutflux = 0.d0
Yl_N(:) = 0.d0
Yr_N(:) = 0.d0
Yl_Noutflux(:) = 0.d0
Yr_Noutflux(:) = 0.d0
O_Noutflux = 0.d0
TotalCarbo = 0.d0
soilC(:,:,:,:) = 0.d0
psi_soil = 0.d0
nyear = 0
Yr_input(:) = 0.d0 
Yb_input(:) = 0.d0
Yl_input(:) = 0.d0
