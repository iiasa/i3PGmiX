! Writing the main output

! Climatic variables -------------
output(ii,:,1,1) = tmp_min(ii)
output(ii,:,1,2) = tmp_max(ii)
output(ii,:,1,3) = tmp_ave(ii)
output(ii,:,1,4) = frost_days(ii)
output(ii,:,1,5) = solar_rad(ii)
output(ii,:,1,6) = day_length(month)
output(ii,:,1,7) = prcp(ii)
output(ii,:,1,8) = vpd_day(ii)
output(ii,:,1,9) = co2(ii)
output(ii,:,1,10) = d13catm(ii)
output(ii,:,1,11) = snow_prcp(ii)
output(ii,:,1,12) = t_soil
output(ii,:,1,13) = t_soil_surf
output(ii,:,1,14) = psi_soil

! Stand variables ---------------
output(ii,:,2,1) = age(ii,:)
output(ii,:,2,2) = stems_n(:)
output(ii,:,2,3) = basal_area(:)
output(ii,:,2,4) = basal_area_prop(:)
output(ii,:,2,5) = dbh(:)
output(ii,:,2,6) = height(:)
output(ii,:,2,7) = bias_scale(15,:) ! relative height
output(ii,:,2,8) = crown_length(:)
output(ii,:,2,9) = crown_width(:)
output(ii,:,2,10) = volume(:)
output(ii,:,2,11) = volume_mai(:)
output(ii,:,2,12) = volume_change(:)
output(ii,:,2,13) = volume_cum(:)
output(ii,:,2,14) = harvesting(:)
output(ii,:,2,15) = dead_volume(:)


! Canopy variables ---------------
output(ii,:,3,1) = layer_id(:)
output(ii,:,3,2) = sla(ii,:)
output(ii,:,3,3) = lai(:)
output(ii,:,3,4) = lai_above(:)
output(ii,:,3,5) = lai_sa_ratio(:)
output(ii,:,3,6) = canopy_vol_frac(:)
output(ii,:,3,7) = canopy_cover(:)
output(ii,:,3,8) = lambda_v(:)
output(ii,:,3,9) = lambda_h(:)
output(ii,:,3,10) = aero_resist(:)
output(ii,:,3,11) = vpd_sp(:)


! Stocks variables ---------------
output(ii,:,4,1) = biom_stem(:)
output(ii,:,4,2) = biom_root(:)
output(ii,:,4,3) = biom_foliage(:)
output(ii,:,4,4) = biom_tree(:)
output(ii,:,4,5) = wood_density(ii,:)
output(ii,:,4,6) = fracBB(ii,:)
output(ii,:,4,7) = biom_loss_foliage(:)
output(ii,:,4,8) = biom_loss_root(:)
output(ii,:,4,9) = biom_incr_foliage(:)
output(ii,:,4,10) = biom_incr_root(:)
output(ii,:,4,11) = biom_incr_stem(:)
output(ii,:,4,12) = biom_sapwood(:)
output(ii,:,4,13) =  biom_fineroot(:)
																

! Modifiers ---------------
output(ii,:,5,1) = f_age(ii,:)
output(ii,:,5,2) = f_vpd(:)
output(ii,:,5,3) = f_tmp(ii,:)
output(ii,:,5,4) = f_tmp_gc(ii,:)
output(ii,:,5,5) = f_frost(ii,:)
output(ii,:,5,6) = f_sw(:)
output(ii,:,5,7) = f_nutr(:)
output(ii,:,5,8) = f_calpha(ii,:)
output(ii,:,5,9) = f_cg(ii,:)
output(ii,:,5,10) = f_phys(:)
output(ii,:,5,11) = gammaF(ii,:)
output(ii,:,5,12) = f_transp_scale
output(ii,:,5,13) = fertility(:)	

! Production ---------------
output(ii,:,6,1) = gpp(:)
output(ii,:,6,2) = npp_f(:)
output(ii,:,6,3) = apar(:)
output(ii,:,6,4) = fi(:)
output(ii,:,6,5) = alpha_c(:)
output(ii,:,6,6) = epsilon_gpp(:)
output(ii,:,6,7) = epsilon_npp(:)
output(ii,:,6,8) = epsilon_biom_stem(:)
output(ii,:,6,9) = npp_fract_stem(:)
output(ii,:,6,10) = npp_fract_foliage(:)
output(ii,:,6,11) = npp_fract_root(:)
output(ii,:,6,12) = pFS(:)
output(ii,:,6,13) = biom_foliage_debt(:)
if (maint_resp .eq. int(2)) then
output(ii,:,6,14) = ra(:)
end if
									 
! Water use ---------------
output(ii,:,7,1) = conduct_canopy(:)
output(ii,:,7,2) = conduct_soil
output(ii,:,7,3) = evapotra_soil
output(ii,:,7,4) = prcp_interc(:)
output(ii,:,7,5) = prcp_interc_fract(:)
output(ii,:,7,6) = prcp_runoff
output(ii,:,7,7) = irrig_supl
output(ii,:,7,8) = wue(:)
output(ii,:,7,9) = wue_transp(:)
output(ii,:,7,10) = evapo_transp
output(ii,:,7,11) = transp_veg(:)
output(ii,:,7,12) = asw
output(ii,:,7,13) = water_runoff_polled
output(ii,:,7,14) = snow_water
output(ii,:,7,15) = snowmelt	
						  
! Mortality ---------------
output(ii,:,8,1) = biom_tree_max(:)
output(ii,:,8,2) = gammaN(ii,:)
output(ii,:,8,3) = mort_thinn(:)
output(ii,:,8,4) = mort_stress(:)
output(ii,:,8,5) = mort_manag(:)

! Crop trees and harvested trees ---------------
output(ii,:,8,6) = dbh_crop(:)
output(ii,:,8,7) = height_crop(:) 
output(ii,:,8,8) = basal_area_crop(:)
output(ii,:,8,9) = biom_stem_crop(:)
output(ii,:,8,10) = vol_crop(:)
output(ii,:,8,11) = dbh_harv(:)
output(ii,:,8,12) = height_harv(:) 
output(ii,:,8,13) = biom_stem_harv(:)
output(ii,:,8,14) = stem_harv(:)


! Wood Delta ------------------
output(ii,:,9,1) = Gc_mol(:)
output(ii,:,9,2) = Gw_mol(:)
output(ii,:,9,3) = D13CNewPS(:)
output(ii,:,9,4) = D13CTissue(:)
output(ii,:,9,5) = InterCi(:) * 1000000.d0


! Weibull ---------------------
output(ii,:,10,1) = bias_scale(1,:)
output(ii,:,10,2) = bias_scale(2,:)
output(ii,:,10,3) = bias_scale(3,:)
output(ii,:,10,4) = bias_scale(4,:)
output(ii,:,10,5) = bias_scale(5,:)
output(ii,:,10,6) = bias_scale(6,:)
output(ii,:,10,7) = bias_scale(7,:)
output(ii,:,10,8) = bias_scale(8,:)
output(ii,:,10,9) = bias_scale(9,:)
output(ii,:,10,10) = bias_scale(10,:)
output(ii,:,10,11) = bias_scale(11,:)
output(ii,:,10,12) = bias_scale(12,:)
output(ii,:,10,13) = bias_scale(13,:)
output(ii,:,10,14) = bias_scale(14,:)

! Soil ---------------------
if (soil_model .eq. int(1)) then
output(ii,:,11,1) = Yl_C(:)
output(ii,:,11,2) = Yr_C(:)
output(ii,:,11,3) = O_C
output(ii,:,11,4) = O_N
output(ii,:,11,5) = (Yl_Coutflux(:) + Yr_Coutflux(:))
output(ii,:,11,6) = O_Coutflux
output(ii,:,11,7) = TotalCarbo
output(ii,:,11,8) = TotalNitro
output(ii,:,11,9) = Nav
output(ii,:,11,10) = Un(:)
output(ii,:,11,11) = Yl_input
output(ii,:,11,12) = Yr_input
output(ii,:,11,13) = kr
else if (soil_model .eq. int(2)) then
output(ii,:,11,1) = Yl_input(:)
output(ii,:,11,2) = Yr_input(:)
output(ii,:,11,3) = Yb_input(:)
output(ii,:,11,4) = sum(soilC(ii,i,1,:))
output(ii,:,11,5) = (Yl_Coutflux(:) + Yr_Coutflux(:))
output(ii,:,11,6) = O_Coutflux
output(ii,:,11,7) = O_C
output(ii,:,11,8) = TotalNitro
output(ii,:,11,9) = Nav
output(ii,:,11,10) = Un(:)
output(ii,:,11,11) = sum(soilC(ii,:,:,1))
output(ii,:,11,12) = sum(soilC(ii,:,:,2))
output(ii,:,11,13) = sum(soilC(ii,:,:,3))
output(ii,:,11,14) = sum(soilC(ii,:,:,4))
output(ii,:,11,15) = sum(soilC(ii,:,:,5))
end if

! Disturbances ---------------

if (wind_dist .eq. int(1)) then
output(ii,:,12,1) = cws(:)
output(ii,:,12,2) = cws10(:)
output(ii,:,12,3) = wind_mortality(:) * (biom_stem(:) / stems_n(:)) * (1.d0 - fracBB(ii,:)) / wood_density(ii,:)
output(ii,:,12,4) = 1.2d0 * wind_mortality(:) / output(ii-1,:,2,2)
end if

if (beetle_dist .eq. int(1)) then
output(ii,:,12,5) = infestation_share
output(ii,:,12,6) = biom_stem_bb(spruce_idx) 
output(ii,:,12,7) = biom_stem_bb(spruce_idx) * (1.d0 - fracBB(ii,spruce_idx)) / wood_density(ii,spruce_idx)
end if
