module wind_disturbance

	implicit none
	
	public :: get_cws_tmc, get_damage_probability, get_thinning_effect !,get_gust
	
	contains
	
	function get_crown_wind_speed(cw,cl,h,stems,w_speed,cdrag, ndrag, lai, n)
	
		integer, intent(in) :: n
		real(kind=8), intent(in) :: lai
		real(kind=8), dimension(n), intent(in) :: cw, cl, h, stems, w_speed, cdrag, ndrag
		real(kind=8):: porosity, cdl, surface_drag_coefficient, element_drag_coefficient, &
					   kaman_constant, caf 
		real(kind=8), dimension(n) :: lambda, dzero, zzero, u_factor, u_crown,  lambda_drag, &
									 gamma, w_speed_corr
		real(kind = 8) :: get_crown_wind_speed(n)
		
		! Bounds on wind for drag coeff
		w_speed_corr(:) = w_speed(:)

		where (w_speed_corr(:) < 10.d0)
			w_speed_corr(:) = 10.d0
		end where
		where (w_speed_corr(:) > 25.d0)
			w_speed_corr(:) = 25.d0
		end where

		! frontal area index
		caf = 0.74d0 ! Assuming constant due to the low variation
    
		porosity = 0.5d0
		lambda(:) = (cw(:) * cl(:)) * caf * porosity / (10000.d0 / sum(stems(:)))
		lambda(:) = (cw(:) * cl(:)) * cdrag(:)*(w_speed_corr(:))**(-ndrag(:)) / (10000.d0 / sum(stems(:)))
		
		cdl = 7.5d0
		dzero(:) = h(:) * ( 1.d0 - (1.d0-exp(-sqrt(cdl*lambda(:))))/(sqrt(cdl*lambda(:))))

		surface_drag_coefficient = 0.003d0
		element_drag_coefficient = 0.3d0
		kaman_constant = 0.4d0

		lambda_drag(:) = min(lambda(:), 0.6d0)   	
		
		gamma(:) = 1.d0/sqrt(surface_drag_coefficient + element_drag_coefficient*lambda_drag(:)/2.d0)

		zzero = (h(:) - dzero(:))*exp(-kaman_constant*gamma(:) + 0.193d0)

		u_factor(:) = log(10.d0/zzero(:)) / log((h(:)*1.05d0-dzero(:))/zzero(:))  
		
		! Based on Yi 2008
		!u_factor(:) =  exp(0.5d0 * lai * (((10.d0 + dzero(:))/zzero(:))/(h(:)*1.05)-1.d0) ) / &
		!				exp(0.5d0 * lai * ((dzero(:)/zzero(:))/(h(:)*1.05d0) - 1.d0)) 

		u_crown(:) = w_speed(:) * u_factor(:)

		get_crown_wind_speed(:) = u_crown(:)
	
	end function get_crown_wind_speed
	
	
	! !Correct gust speed above canopy 
	! function get_gust(crown_width,crown_length,height,crown_area_factor,wind_speed,n)
	
	! 	! Get wind gust speed above the canopy - based on iLand implementation
	! 	! n = number of species
	! 	! crown_width = crown width(m)
	! 	! crown_length = crown_length (m)
	! 	! height = tree height (m)
	! 	! crown_area_factor
	! 	! wind_speed = wind gust speed (m/s)
	! 	! n = number of tree species in the stand
	
	! 	integer, intent(in) :: n
	! 	real(kind=8), dimension(n) :: u_crown, crown_width, crown_length,height, crown_area_factor, wind_speed
	! 	real(kind=8) :: corr_gust
	! 	real(kind=8) :: topo_mod, gust_factor, r_num, r_num_aux
	! 	real(kind = 8)  ::  get_gust(n)
		
	! 	! get deviate for wind gust calculation
	! 	r_num_aux = rand()
		
	! 	!random speed modifier
	! 	if (r_num_aux <= 0.5d0) then
	! 		gust_factor = 1.d0 + rand()
	! 	else
	! 		gust_factor = 1.d0 - rand()
	! 	end if
		
	! 	!topographic exposure modifier
	! 	topo_mod = 1.d0
		
	! 	!corrected gust speed		
	! 	corr_gust = wind_speed * topo_mod * gust_factor
		
	! 	!calculate crown wind speed		
	! 	u_crown(:) = get_crown_wind_speed(crown_width,crown_length,height,crown_area_factor,corr_gust, n)

	! 	!return wind speed at the crown
	! 	get_gust(:) = u_crown(:)
		
	! end function get_gust
	
	function calc_dlf(h, dbh, cr_depth, cr_width, moe, stem_weight, crown_density, breaking_moment, tc, n)
		
		! calculate deflection
		! h = tree height (m)
		! dbh = diameter at breast height (cm)
		! cr_depth = crown length (m)
		! cr_width = crown width (m)
		! moe = modulus of elasticity (MPa)
		! stem_weight = weight of wet stem (t of biomass) 
		! crown_density (unitless)
		! breaking_moment = breaking moment ( )
		! tc =  turning coefficient 
		! n =  number of species in the stand
		
		integer, intent(in) :: n 
		real(kind=8) :: an, bn, cn, grav, snow_density, snow_depth, hh, pi
		real(kind=8), dimension(n), intent(in) :: h, cr_depth, cr_width, moe, &
												  stem_weight, breaking_moment, &
												  tc, dbh, crown_density
		real(kind=8), dimension(n) :: bm, wind_action_point, i_value, lever_arm, &
									  r_val, x, deflection_1, deflection_2, dbase, &
									  wind_force, DLF_calc
		real(kind=8) :: calc_dlf(n)
		
		!define constants and parameters for computations
		an = -1.4d0  
		bn = -0.4d0
		cn = 0.6d0
		grav = 9.81d0
		snow_density = 150.d0
		snow_depth = 0.d0 !for now not included
		pi = 3.141592654d0
		
		hh = an*bn*cn
		
		!calculate bending moment
		bm(:) = tc(:) * sqrt(breaking_moment(:) / tc(:))
		
		!calculate wind action point
		wind_action_point(:) = h(:) - (cr_depth(:)/2.d0)
		where (wind_action_point(:) < 1.3*2.d0 ) 
			wind_action_point(:) = 1.3*2.d0
		end where
		
		!diameter at tree base
		dbase(:) = (dbh(:) / (wind_action_point(:) - 1.3d0)**0.333) * wind_action_point(:)**0.333d0 
		
		!second area moment of Inertia calculated at tree base
		i_value(:) = pi * ((dbase(:))**4.d0 / 64.d0)
		
		!get force of wind
		wind_force(:) = bm(:) / (h(:) - cr_depth(:)/2.d0)
		
		!get lever arm
		lever_arm(:) = (h(:) - cr_depth(:)/2.d0)
		
		!ratio of lever arm to total tree height
		r_val(:) = (h(:) - lever_arm(:)) / h(:)
		
		!get deflection at mid canopy
		x(:) = cr_depth(:)/2.d0
		deflection_1(:) =  ((wind_force(:)*lever_arm(:)**3.d0) / (moe(:) * i_value(:) *hh)) * &
							(an*(x(:)/lever_arm(:))**cn - r_val(:) * cn*(x(:)/lever_arm(:))**bn + &
							r_val(:) * bn*cn*(x(:)/lever_arm(:)) - an*cn*(x(:)/lever_arm(:)) + &
							an*bn - r_val(:) *an*cn)
		
		!get deflection at 3/4 of tree height
		x(:) = 3.d0/4.d0 * h(:)													  
		deflection_2(:) =  ((wind_force(:)*lever_arm(:)**3.d0) / (moe(:) * i_value(:) *hh)) * &
							(an*(x(:)/lever_arm(:))**cn - r_val(:) * cn*(x(:)/lever_arm(:))**bn + &
							r_val(:) * bn*cn*(x(:)/lever_arm(:)) - an*cn*(x(:)/lever_arm(:)) + &
							an*bn - r_val(:) *an*cn)
		
		!get deflection
		DLF_calc(:) =  1.d0 / (1.d0 - ((deflection_1(:) * ((((pi * cr_depth(:) * (cr_width(:)/2.d0)**2.d0)/3.d0) * &
						crown_density(:)) + (pi * (cr_width(:)/2.d0)**2.d0 * snow_depth * snow_density)) * grav) + &
						deflection_2(:) * stem_weight(:)  * grav) / bm(:))
		
		where(DLF_calc(:) < 1.d0)
			DLF_calc(:) = 1.136d0
		end where
		where(DLF_calc(:) > 2.d0)
			DLF_calc(:) = 1.136d0
		end where
		
		calc_dlf(:) = DLF_calc(:)
	
	end function calc_dlf
	
	
	function get_cws_tmc(h, dbh, cr_depth, cr_width, moe, mor, biom_stem, &
						wet_factor, competition, nha, c_reg, snow_water, corr_comp, thinn_adj, n)
	
		integer, intent(in) :: n
		logical, intent(in) :: corr_comp
		integer :: i
		real(kind=8), dimension(n), intent(in) :: h, cr_width, cr_depth, moe, mor, biom_stem, &
												  nha, c_reg, dbh, competition, wet_factor, thinn_adj
		real(kind=8), intent(in) ::  snow_water
		real(kind=8) :: snow_depth, tree_heights_inside_forest, pi, dlf, maximum_gap_size, &
						ro_default, soil_group, snow_weight
		real(kind=8), dimension(n) :: density, edge_gap_factor, s_h, rooting, gap_size, dist_edge, &
									  wind_action_point, dbase, breaking_moment, overturning_moment,&
									  tc, dlf_used, cws_u, cws_b, stem_weight, fknot, &
									  ci_value, cws, crown_density
		real(kind=8) :: get_cws_tmc(n)
		
		ro_default = 1.2226d0
		tree_heights_inside_forest  = 9.d0
		maximum_gap_size = 10.d0
		dlf = 1.136d0
		pi = 3.141592654d0
		crown_density(:) = 0.5d0
	
		!get average spacing
		density(:) = sqrt(10000.d0 / sum(nha(:)))
		
		soil_group = 1.d0
		rooting(:) = 3.d0
		gap_size(:) = 0.d0
		dist_edge(:) = tree_heights_inside_forest * h(:)
		fknot(:) = 1.d0
		
		stem_weight(:) = biom_stem(:) * wet_factor(:)
		snow_weight = snow_water / 1000.d0
		
		!calculate edge factor
		s_h(:) = density(:) / h(:)
		where (s_h(:) < 0.075d0)
			s_h(:) = 0.075d0
		end where
		where (s_h(:) > 0.45d0)
			s_h(:) = 0.45d0
		end where
		 
		edge_gap_factor(:) = (2.7193d0 * s_h(:) -0.061d0) + (-1.273d0 * s_h(:) + 0.9701d0) * &
							 ((1.1127d0 * s_h(:) + 0.0311d0)**(tree_heights_inside_forest)) + &
                        (-1.273d0 * s_h(:) + 0.9701d0) *(((1.1127d0 * s_h(:) + 0.0311d0)**(dist_edge/h(:))) - &
						((1.1127d0 * s_h(:) + 0.0311d0)**(tree_heights_inside_forest))) * h(:) / (2.7193d0 * s_h(:) -0.061d0) + &
						(-1.273d0 * s_h(:) + 0.9701d0) * ((1.1127d0 * s_h(:) + 0.0311d0)**(tree_heights_inside_forest))

		
		!calculate moments
		wind_action_point(:) = h(:) - (cr_depth(:)/2.d0)
		where (wind_action_point(:) < (1.3*2.d0) ) 
			wind_action_point(:) = 1.3d0*2.d0
		end where
		
		dbase(:) = (dbh(:) / (wind_action_point(:) - 1.3d0)**0.333d0) * wind_action_point(:)**0.333d0 
  
		breaking_moment(:) = mor(:) * fknot(:) * pi * ((dbase(:)/100.d0)**3.d0) / 32.d0 
		overturning_moment(:) = c_reg(:) * stem_weight(:)
		
		ci_value(:) = log(competition(:))
			
		if (corr_comp .eqv. .TRUE. ) then
			tc(:) = 3.86d0*ci_value(:) + 124.252d0*(dbh(:)/100.d0)**2.d0 * h(:) - 17.285d0*(dbh(:)/100.d0)**2.d0 * h(:) *ci_value(:)
		else
			tc(:) = 111.604d0*(dbh(:)/100.d0)**2.d0 * h(:) 
		end if
		
		!edge_gap_factor(:) = 1.05d0
		dlf_used(:) = calc_dlf(h, dbh, cr_depth, cr_width, moe, stem_weight, crown_density, breaking_moment, tc, n)
		!dlf_used(:) = 1.136d0
			
		cws_u(:) = sqrt(overturning_moment(:) / (tc(:)*dlf_used(:)*thinn_adj(:)*edge_gap_factor(:)))
		cws_b(:) = sqrt(breaking_moment(:) / (tc(:)*dlf_used(:)*thinn_adj(:)*edge_gap_factor(:)))
		 
		do i=1, n
			cws(i) = min(cws_u(i), cws_b(i))
		end do
	
		get_cws_tmc(:) = cws(:)
		!get_cws_tmc(:) = stem_weight(:)
		
	end function get_cws_tmc
	
	function get_critical_wind_speed(tmp_soil, gap_length, h, dbh, biom_stem, wet_factor, creg, competition, n)
	
		integer, intent(in) :: n 
		integer :: i
		real(kind=8), intent(in) :: tmp_soil, gap_length
		real(kind=8), dimension(n), intent(in) :: h, biom_stem, wet_factor,creg, competition, dbh
		real(kind=8) :: pi
		real(kind=8), dimension(n) :: rel_gap, f_gap, f_knot, cws_u, cws_b, MOR, ci, tc, f_edge, cws, biom_weight
		real(kind=8) :: get_critical_wind_speed(n)
		
		!relative gap length
		rel_gap(:) = min(gap_length / h(:), 10.d0)
		
		!formulation from Peltola et al. (1999), based on Gardiner et al. (1997) 0.2158... = scale to situation with gapsize = 0
		f_gap = (0.001d0+0.001d0*rel_gap**0.562d0)/(0.00465d0) 

		!calculate the wet stem weight 
		biom_weight(:) = biom_stem(:) * wet_factor(:)

		!calculate the turning moment coefficient from tree level ForestGales
		ci(:) = log(competition(:)) 
		
		tc(:) = 4.42d0+122.1d0*dbh(:)*h(:)-0.141d0*ci(:)-14.6d0*dbh(:)*h(:)*ci(:)
    
		! now derive the critital wind speeds for uprooting and breakage
	
		!account for the presence of knots
		f_knot(:) = 1.d0
		pi = 3.141592654d0
		
		!calculate cws for uprooting and breakage

		! Turning moments at stand edges are significantly higher at stand edges compared to conditions
		! well inside the forest. Data from Gardiner(1997) and recalculations from Byrne (2011) indicate,
		! that the maximum turning moment at the edge is about 5 times as high as "well inside" the forest.
		f_edge(:) = 5.d0;

		cws_u(:) = sqrt(biom_stem(:) * creg(:) / (tc(:)*f_gap(:)*f_edge(:)))
		cws_b(:) = sqrt(MOR(:)*dbh(:)**3.d0*f_knot(:)*pi/(32.d0*tc(:)*f_gap(:)*f_edge(:)))
		
		do i=1, n
			cws(i) = min(cws_u(i), cws_b(i))
		end do
		
		get_critical_wind_speed(:) = cws(:)

	end function get_critical_wind_speed
	
	function get_damage_probability(a,k, u_crit, n)
	
		integer, intent(in) :: n
		real(kind=8), intent(in) :: a, k
		real(kind=8), dimension(n), intent(in) :: u_crit
		real(kind=8) :: u_a, u_c1, u_c2, u_c3, u_c4
		real(kind=8), dimension(n) :: u, u_c, u_aux, aep
		real(kind=8) :: get_damage_probability(n)
		
		u_a = 5.d0
		u_c1 = -0.5903d0
		u_c2 = 4.4345d0
		u_c3 = -11.8633d0
		u_c4 = 13.569d0
		

		! Check if wind speed is based on weibull distribution or gust
		if (k == 0.d0) then
			u(:) = a
		else
			u_c(:) = u_c1 * (k)**3.d0 + u_c2 * (k)**2.d0 + u_c3 * k + u_c4
		
			u(:) = (u_c(:) * a)**2.d0		
			
		end if
		
		u_aux(:) = u(:) / u_a

		aep(:) = 1.d0 - exp(-1.d0 * exp(-1.d0 * (u_crit(:)**2.d0 - u(:)) / u_aux(:)))
		
		get_damage_probability(:) = aep(:)

	end function get_damage_probability

	function get_damage_wind(a,k, u_crit, stems, rnd) result( damage )
	
		real(kind=8), intent(in) :: a, k, u_crit, stems, rnd
		real(kind=8) :: u_a, u_c1, u_c2, u_c3, u_c4, u, u_c, u_aux, gust, rf, gust_corr
		real(kind=8) :: damage
		
		! Gumbel distribution parameters
		u_a = 5.d0
		u_c1 = -0.5903d0
		u_c2 = 4.4345d0
		u_c3 = -11.8633d0
		u_c4 = 13.569d0
		
		! Damage smoothing parameter
		rf = 12.d0

		gust = 0.d0
		gust_corr = 1.d0

		! Check if wind speed is based on weibull distribution or gust
		if (k <= 0.d0) then
			gust = a
			gust_corr = 1.d0 / (1.d0 + exp(-(gust-25.d0)/6.d0))
		else
			u_c = u_c1 * (k)**3.d0 + u_c2 * (k)**2.d0 + u_c3 * k + u_c4
		
			u = (u_c * a)**2.d0
			
			u_aux = u / u_a
			
			! Random draw from the Gumbel distribution
			gust = sqrt( u - u_aux * log(-log(rnd)) )
			gust_corr = 1.d0 / (1.d0 + exp(-(gust-25.d0)/6.d0))
			
		end if

		
		
		! Get damage according to Chen et al. (2018) Simulating damage for wind storms in the land surface model ORCHIDEE-CAN.  Eq.9
		! Added an additional factor to reduce damage for low wind speeds (< 20m/s)
		damage = stems * ((1.d0 / (1.d0 + exp(- ( gust * gust_corr - u_crit) / rf ))) - (1.d0 / (1.d0 + exp(u_crit/rf))))

	end function get_damage_wind
	
	! Get effect of thinning on wind risk
	function get_thinning_effect(curr_n,prev_n, lag, n) 
	
		integer, intent(in)::n
		integer:: i
		real(kind=8), dimension(n), intent(in) :: curr_n,prev_n, lag
		real(kind=8), dimension(n) :: spacing_current, spacing_before
		real(kind=8) :: get_thinning_effect(n)
	
		get_thinning_effect(:) = 1.d0
		spacing_current(:) = sqrt(10000.d0/max(curr_n(:),1.d0))
		spacing_before(:) = sqrt(10000.d0/max(prev_n(:),1.d0))

		! now it is calculating per species separately, but under the assumption of homogeneous mixing,
		! the average value is used for the TMR accorection for all species
		do i=1, n
			if (lag(i) >= 5.d0) then
				get_thinning_effect(i) = 1.d0
			else
				get_thinning_effect(i) = max((0.99d0 * spacing_current(i)/spacing_before(i)) * &
										 (1.d0 - lag(i)/5.d0) + (lag(i)/5.d0), 1.d0)
			end if
		end do

	end function get_thinning_effect
	
end module wind_disturbance

