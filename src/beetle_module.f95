module beetle_disturbance

	use beetle_subroutines

	! LandClim bark beetle model from Temperli et al. (2013) and extension by Marie et al. (2023):
	!
	! Temperli, C., Bugmann, H., & Elkin, C. (2013). Cross‐scale interactions among bark beetles, 
	! climate change, and wind disturbances: A landscape modeling approach. 
	! Ecological Monographs, 83(3), 383-402.
	!
	! Marie, G., Jeong, J., Jactel, H., Petter, G., Cailleret, M., McGrath, M., ... & Luyssaert, S. (2023). 
	! Simulating Bark Beetle Outbreak Dynamics and their Influence on Carbon Balance Estimates with 
	! ORCHIDEE r7791. EGUsphere, 2023, 1-35.

	implicit none
	
	public :: bb_landscape_risk, bb_mortality
	
	contains

	! Estimate the landscape risk index based on the susceptibility and bark beetle pressure index
	subroutine bb_landscape_risk(	infestation_share, &
									bb_mort, &
									age, & 
									f_sw , d_type, & 
									wind_damage, sp_biom_stem, & 
									sp_ba_share, & 
									w_bio, w_age, w_host, w_vita, w_dens, w_wind, &
									greff, greff_sp, vita_type, &
									previous_damage, max_sp_biomass, c_st, & 
									lat, tmean, wtmin, rad, lai, k, day_length, c_bp, & 
									background_pressure	)

		! input 
		integer, intent(in):: d_type, vita_type
		real(kind=8), intent (in):: age, w_wind, & 
						wind_damage, sp_biom_stem, sp_ba_share, w_bio, previous_damage, &
						max_sp_biomass, c_st, c_bp, background_pressure, lat, greff, greff_sp
						 

		real(kind=8), intent(in):: k, wtmin
		real(kind=8), dimension(12), intent(in):: tmean, rad, lai, day_length, f_sw	

		! output 
		real(kind=8), intent(out):: infestation_share, bb_mort

		! local
		real(kind=8):: susceptibility_idx, bpi_idx, risk_idx, beta_scale, beta_s1, beta_s2, rX, rY, &
						w_age, w_host, w_vita, w_dens, age_sus, host_sus, wind_sus, drought_sus, &
						dens_sus, c_bp_aux, c_st_aux, lag_idx
						integer:: epidemic

		! Check for epidemic phase
		epidemic = 0
		!call get_weighting(min(bpi_idx,1.d0), w_age, w_host, w_vita, w_dens, w_bio, w_wind, epidemic)
		call get_weighting(wind_damage, w_age, w_host, w_vita, w_dens, w_bio, w_wind, epidemic)

		if (epidemic > 0) then
			c_bp_aux = max(c_bp, 1.d0)
			c_st_aux = max(c_st, 1.d0)
		else
			c_bp_aux = c_bp
			c_st_aux = c_st
		end if

		! Get compound susceptibility
		call local_susceptibility( 	susceptibility_idx, &
									age_sus, host_sus, wind_sus, drought_sus, dens_sus, &
									age, & 
									f_sw , d_type, & 
									greff, greff_sp, vita_type, &
									wind_damage, sp_biom_stem, & 
									sp_ba_share, & 
									w_bio, w_wind) 

		
		! Get BB pressure index
		call bb_pressure(	bpi_idx, susceptibility_idx, &
							previous_damage, max_sp_biomass, c_st_aux,  & 
							lat, tmean, rad, lai, k, day_length, c_bp_aux, & 
							background_pressure, lag_idx)

		! Get average risk index
		risk_idx  = min(susceptibility_idx,1.d0)*min(bpi_idx,1.d0)

		! Scale of the beta function for the stochastic component
		beta_scale = 10000.d0

		! Stochastic component - beta random variable
		beta_s1 = 1.d0 +(risk_idx*beta_scale)
		beta_s2 = 1.d0 +((1.d0-risk_idx)*beta_scale)

		!call rand_gamma_scalar(beta_s1,rX)
		!call rand_gamma_scalar(beta_s2,rY)
		call  random_gamma(beta_s1,1.d0,rX)
		call  random_gamma(beta_s2,1.d0,rY)
		!rX = gengam(beta_s1, beta_scale)
		!rY = gengam(beta_s2, beta_scale)

		! Calculate infestation share in landscape
		infestation_share  = rX/(rX+rY)

		if(risk_idx < 0.0025d0 .or. wtmin < -10.d0) infestation_share  = 0.d0 ! Limits damage in areas with no wind risk

		call bb_mortality (bb_mort, infestation_share, age_sus, drought_sus, dens_sus, wind_sus, &
							w_age, w_vita, w_dens, bpi_idx, c_st, epidemic)
							infestation_share = gen_effect(calc_rdd(lat, tmean, rad, lai, k, day_length))  
	end subroutine bb_landscape_risk
	
	! Mortality rate of infested cohorts
	subroutine bb_mortality( bb_mort, infestation_share, age_idx, drought_idx, dens_idx, wind_idx, &
							w_age, w_vita, w_dens, bpi_idx, c_st, epidemic )

		!input
		integer, intent(in):: epidemic
		real(kind=8), intent(in):: infestation_share, age_idx, drought_idx, dens_idx, wind_idx, &
									w_age, w_vita, w_dens, bpi_idx, c_st
		
		!output
		real(kind=8), intent(out):: bb_mort

		!local
		real(kind=8):: bio_idx

		!bio_idx = (age_idx * w_age + drought_idx *w_vita + dens_idx * w_dens) / (w_age + w_vita + w_dens)
		bio_idx = max( (age_idx * w_age + drought_idx *w_vita) / (w_age + w_vita), wind_idx )
		!bio_idx = max((dens_idx * w_dens + drought_idx *w_vita) / (w_dens + w_vita), wind_idx )
		!bio_idx = (age_idx * w_age + drought_idx *w_vita) / (w_age + w_vita)
		!bio_idx = max(0.5d0 *  (age_idx + drought_idx), wind_idx ) 

		!if (epidemic > 0) then
		!	bb_mort =  0.5d0 *( bio_idx + wind_idx ) * infestation_share
		!else
			bb_mort = 0.5d0 * ( bio_idx + bpi_idx ) * infestation_share 
		!end if

	end subroutine bb_mortality
	
	! Beetle pressure index 
	subroutine bb_pressure(	bdi, susceptibility_idx, &
							previous_damage, max_sp_biomass, c_st, & ! Lag effect
							lat, tmean, rad, lai, k, day_length, c_bp, & ! Generation effect
							background_pressure, lag_idx) ! Background infestation
		
		!input
		real(kind=8), intent(in):: previous_damage, max_sp_biomass, c_bp, &
									c_st, lat, background_pressure, k, susceptibility_idx
		real(kind=8), dimension(12), intent(in):: tmean, rad, lai, day_length

		!output
		real(kind=8), intent(out):: bdi, lag_idx

		!local
		real(kind=8) :: gen_idx, bdi_base, rel_gdd !lag_idx, 

		! Get previous damage effect
		lag_idx = lag_effect(previous_damage, max_sp_biomass, c_st)

		! Get generation effect
		rel_gdd = calc_rdd(lat, tmean, rad, lai, k, day_length)
		gen_idx = gen_effect(rel_gdd)


		! Bark beetle pressure index
		bdi_base = 0.5d0 * (lag_idx + gen_idx)  * susceptibility_idx * c_bp

		! Return background beetle pressure if beetle pressure is lower than background beetle pressure
		bdi = max(bdi_base, background_pressure)

		if (rel_gdd < 1.d0) then
			bdi = 0.d0
		end if

	end subroutine bb_pressure

	! Beetle generations Index (G, Eq. A10 in Appendix A)
	function gen_effect(rel_gdd) result(bb_gen_idx)
		
		! input
		!real(kind=8), dimension(12), intent(in) :: tmean, rad, day_length, lai, rel_gdd
		!real(kind=8):: lat, k
		real(kind=8), intent(in):: rel_gdd

		! local
		real(kind=8) :: bb_gen_idx, pr, PM

		!rel_gdd = calc_rdd(lat, tmean, rad, lai, k, day_length)

		pr = 3.307963d0; !Parameters for logistic function using original hazard scores in Fig. 7
		pM = 1.980938d0

		bb_gen_idx = (1.d0/(1.d0 + exp(-pr*(rel_gdd-pM)))) !Index based on hazard rating by Baier et al. 2007 and Netherer and Nopp-Mayr 2005.

	end function gen_effect


	! Impact of previous damage on current damage
	function lag_effect(killed_biomass, max_host_biomass, c_st) result(bb_dam_idx)
		
		! input
		real(kind=8), intent(in) :: killed_biomass, max_host_biomass, c_st

		! local
		real(kind=8) :: bb_dam_idx

		if(killed_biomass / max_host_biomass > 1.d0) then
			bb_dam_idx = 1.d0
		else
			bb_dam_idx = killed_biomass / (max_host_biomass * c_st)
		end if

	end function lag_effect


	! Estimates local susceptibility to BB damage Eq. (2) in Temperli et al. (2013)
	subroutine local_susceptibility(bb_susceptibility, & !total
									age_sus, host_sus, wind_sus, drought_sus, dens_sus, & ! components
									age, & ! age related
									f_sw , d_type, & ! drought related
									greff, greff_sp, vita_type, & ! drought mode
									wind_damage, sp_biom_stem, & ! wind related
									sp_ba_share, & ! host related
									w_bio, w_wind) ! weights
									

		!input 
		integer, intent(in):: d_type, vita_type							
		real(kind=8), intent(in):: age, wind_damage, sp_biom_stem, sp_ba_share, &
									 w_bio, w_wind, greff, greff_sp
		real(kind=8), dimension(12), intent(in):: f_sw

		!output
		real(kind=8), intent(out):: bb_susceptibility

		!local
		real(kind=8), intent(out):: age_sus, host_sus, drought_sus, wind_sus, dens_sus

		age_sus = age_effect(age)
		host_sus = host_effect(sp_ba_share)
		wind_sus = min(wind_damage / 0.1d0 , 1.d0) !wind_effect(wind_damage, sp_biom_stem) ! gets directly the share of damage from wind module

		if (vita_type == int(1)) then 
			drought_sus = drought_effect(f_sw, d_type)
		else
			drought_sus = vitality_effect( greff , greff_sp ) 
		end if

		!dens_sus = dens_effect(stems_n, dbh, thinPower, aWS, nWS, wSx1000)

		bb_susceptibility = (age_sus + host_sus + drought_sus) * w_bio + wind_sus * w_wind
	
	end subroutine local_susceptibility

	! Estimate weighting scheme for the susceptibility, according to ORCHIDEE implementation (Marie et al. 2023)
	subroutine get_weighting(risk , w_age, w_host, w_vita, w_dens, w_bio, w_wind, epidemic )

		! input
		real(kind=8):: risk

		!output
		real(kind=8):: w_age, w_host, w_vita, w_dens, w_bio, w_wind

		!local
		integer, intent(out):: epidemic ! BPI value where epidemic start and weighting changes
		real(kind=8):: epidemic_threshold = 0.1d0 ! Risk index value where epidemic start and weighting changes
		real(kind=8):: r1, r2

		! Needs to be defined (not shown in the paper), for now assigns a 0.5 weight for low wind damage and 0.8 for high wind damage
		r1 = 4.74710078
		r2 = 2.886096695
		r1 = 18.48574061
		r2 = 9.648420753
	

		! Check outbreak phase
		if (risk > epidemic_threshold) then
			epidemic = 1
		else
			epidemic = 0
		end if

		! Tree density effect weight
		w_dens = 0.067d0

		! Host effect weight
		w_host = 0.067d0

		! Age effect weight
		w_age = 0.067d0

		if (epidemic > 0) then

			! Stress weighting
			w_vita = 1.d0 / ( 1.d0 + exp( r1 * risk/0.3d0 - r2 )) * ( 1.d0 - (w_dens + w_host + w_age) ) ! original

			w_bio = (1.d0 / ( 1.d0 + exp( r1 * risk/0.3d0 - r2 )) * ( 1.d0 - (w_dens + w_host + w_age) ) + 0.2d0) / 3.d0

			! Breeding substrate
			w_wind = 1.d0 - ( w_vita + w_dens + w_host + w_age) ! original


			w_wind = 1.d0 - w_bio

		else
			w_wind = 0.d0 ! original

			w_vita = 1.d0 - (w_dens + w_host + w_age) ! original
			w_bio = 1.d0 / 3.d0 

			w_wind = 0.5d0 ! original
			w_bio = 1.d0 / 6.d0 

		end if

	end subroutine get_weighting


	! Age-induced susceptibility (S_a; Eq. A4 in Appendix A):
	function age_effect(age) result(age_idx)
		
		! input
		real(kind=8), intent(in) :: age

		! local
		real(kind=8) :: age_idx

		age_idx = (0.2d0 + 1.d0-0.2d0) / (1.d0 + exp(-0.1094542d0 * (age-70.d0)))

	end function age_effect


	! Drought effect on susceptibility
	function drought_effect( f_sw , type) result(drought_idx)

		!input
		real(kind=8), dimension(12), intent(in) :: f_sw 
		integer, intent(in) :: type

		!local
		real(kind=8):: logFun_M, logFun_r, asw_maxsus, drought_idx, drought_tol, LowerAsymDrS, ParamCDrS, di

		di = (7.d0 - sum(f_sw(4:10))) / 7.d0 ! Drought index: soil water modfier between April and October (growing season)

		! Spruce drought tolerance factor (LandClim default)
		drought_tol = 0.2d0 

		! lower Asymptote of logistic function that relates droughtIndex to drought-induced beetle susceptibility (LandClim default)
		LowerAsymDrS = 0.d0

		! Parameter to define minimum susceptibility in polynomial variant (LandClim default)
		ParamCDrS = 0.d0

		! value at which SI_drought is maximal (LandClim default)
		asw_maxsus = 0.2d0

		!logFun_M = (1.d0 - asw_maxsus**2.d0) * drought_tol/2.d0 
		logFun_M = 4.8d0 * drought_tol / 2.d0 ! Since the soil water modifier in 3PG ranges between 0 and 1 we scale the (1.d0 - asw_maxsus**2.d0) term to yield the same curve in LandClim

		! Default
		! Maximum slope in logistic function: Parameterised so that S_d = 0.99 at droughtIndex = M*2
		if (type .eq. int(1)) then
			logFun_r = -log(1.d0/0.99d0 - 1.d0)/logFun_M
			drought_idx = (1.d0/(1.d0+exp(-logFun_r*(di-logFun_M))))
		! negLAsym
		!Strategy: Assume negative susceptibilityIndexDrought at low droughtIndex values --> set lower
		else if (type .eq. int(2)) then
			logFun_r = -log( (1.d0-0.99d0) / 0.99d0 - LowerAsymDrS)/logFun_M
			drought_idx = LowerAsymDrS + (1.d0 - LowerAsymDrS) / (1.d0+exp(-logFun_r*(di-logFun_M)))
		!polyn
		! Strategy: Subtract unimodal function that does not go below zero --> prob dens. function of beta distribution
		else if (type .eq. int(3)) then
			logFun_r = -log(1.d0/0.99d0 - 1.d0)/logFun_M
			drought_idx = (1.d0/(1.d0+exp(-logFun_r*(di-logFun_M)))) - &
							ParamCDrS * di**1.5d0 * (1.d0-di)**49
		end if

	end function drought_effect

	! Vitality effect on susceptibility - maybe can be used in place of the drought effect, uses greff concept from LPJ-GUESS
	function vitality_effect( greff , greff_sp ) result(vital_idx)

		!input
		real(kind=8), intent(in) :: greff
		real(kind=8), intent(in) :: greff_sp

		!local
		real(kind=8):: vital_idx, si
		real(kind=8) :: greff_aux

		! Get carbon balance status
		greff_aux = greff / greff_sp

		if (greff_aux > 1.d0) then
			greff_aux = 1.d0
		end if

		!Get vitality index
		!si = (7.d0 - greff_aux) / 7.d0 

		! Uses same relationship as iLand SI
		!vital_idx = 1.d0 - (0.85d0 * si + 0.15d0)
		vital_idx = 1.d0 - (0.85d0 * greff_aux + 0.15d0)

	end function vitality_effect


	! Beetle host share-induced susceptibility (S_s; Eq. A5 in Appendix A)
	function host_effect(sp_ba_share) result(host_idx)

		!input
		real(kind=8):: sp_ba_share ! spruce basal area share

		!local
		real(kind=8):: host_idx

		host_idx = (0.08d0 + (1.d0-0.08d0)/ (1.d0 + 0.2674846d0* &
					exp(-6.8455438d0*(sp_ba_share-0.375d0)))**(1.d0/0.3002734d0))

	end function host_effect

	! Density effect from Marie et al. (2023)
	function dens_effect(stems_n, dbh, thinPower, aWS, nWS, wSx1000) result(dens_idx)

		! Need to figure out the parameter values, not provided in the publication

		!input
		real(kind=8):: stems_n, dbh, thinPower, aWS, nWS, wSx1000 

		!local
		real(kind=8):: dens_idx, max_stems, rdi, par_a, a1, a2, a3

		! Needs to be defined (not shown in the paper)
		a1 = 0.d0
		a2 = -12.d0
		a3 = 0.6d0

		! Intercept of self-thinning curve (log-scale)
		par_a = log(1000.d0) + thinPower * log( (wSx1000 / aWS)**(1.d0/nWS) )

		! Maximum density at current dg
		max_stems = exp( par_a - thinPower * log(dbh) )

		! Relative density index
		rdi = min(stems_n / max_stems, 1.d0)

		! Density score
		dens_idx = a1 + (1.d0 - a1) / (1.d0 + exp(a2 * (rdi - a3)))

	end function dens_effect

	! Wind throw-induced susceptibility (S_w)
	function wind_effect(wind_damage, sp_biom_stem) result(wind_idx)
	
		!input
		real(kind=8), intent(in):: wind_damage, sp_biom_stem ! windthrow damage

		!local
		real(kind=8):: wind_idx

		if (wind_damage > sp_biom_stem) then
			wind_idx = 1.d0
		else
			wind_idx = wind_damage / sp_biom_stem
		end if

	end function wind_effect

end module beetle_disturbance


