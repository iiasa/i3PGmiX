module mod_3PG

    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_bool
    
    use allocation
    use beetle_disturbance
    use bias_correction
	use gpp_calc
    use light_int
    use management
    use mod_decl_const
    use mortality
    use phenology
    use ra_calc
    use soil_cnp_subroutines
    use soil_cnp
    use transpiration
    use utils
    use wind_disturbance

    implicit none
    private
    public :: s_3PG_f

contains

    subroutine s_3PG_f ( siteInputs, speciesInputs, forcingInputs, managementInputs, pars_i, pars_b, &
        n_sp, n_m, n_man, t_t, settings, output) bind(C, name = "s_3PG_f_")

        implicit none

        !********************************************************************************************
        !   Declaration

        ! Number of species and month
        integer(kind=c_int), intent(in) :: n_m
        integer(kind=c_int), intent(in) :: n_sp
        integer(kind=c_int), intent(in) :: n_man ! number of management interactiosn
        integer(kind=c_int), dimension(n_sp), intent(in) :: t_t ! number of management interactiosn
        integer(kind=c_int), dimension(18), intent(in) :: settings    ! settings for the models

        ! Initial, forcing, parameters
        real(kind=c_double), dimension(14), intent(in) :: siteInputs
        real(kind=c_double), dimension(n_sp,19), intent(in) :: speciesInputs
        real(kind=c_double), dimension(n_man,8,n_sp), intent(in) :: managementInputs
        real(kind=c_double), dimension(n_m,13), intent(in) :: forcingInputs
        real(kind=c_double), dimension(165,n_sp), intent(in) :: pars_i
        real(kind=c_double), dimension(30,n_sp), intent(in) :: pars_b
		real(kind=c_double), dimension(22) :: out_aux
		real(kind=c_double), dimension(10) :: out_aux_h
        real(kind=c_double), dimension(n_sp,22) :: out_rpmodel
		real(kind=c_double), dimension(n_sp,10) :: out_rphmodel
        		
        ! Output array
        real(kind=c_double), dimension(n_m,n_sp,12,15), intent(inout) :: output

        ! Variables, Parameters, Constants
        include 'i_decl_var.h'

        include 'i_read_input.h'
        include 'i_read_param.h'
        include 'i_read_param_sizeDist.h'
        include 'i_read_param_icbm.h'
        include 'i_read_param_yasso.h'	
		include 'i_read_param_wind.h'
		include 'i_read_param_hydro.h'


        ! Initialization
        include 'i_init_var.h'
		
        !*************************************************************************************
        ! INITIALISATION (Age independent)
		

        ! Day-length calculations
        adjSolarZenithAngle(:) = f_get_solarangle( Lat )

        day_length(:) = 86400.d0 * f_get_daylength( Lat ) !Seconds

        ! CO2 modifiers helpers
        fCalphax(:) = fCalpha700(:) / (2.d0 - fCalpha700(:))
        fCg0(:) = fCg700(:) / (2.d0 * fCg700(:) - 1.d0)
		
		
        ! Temperature --------
        do i = 1, n_sp
            ! calculate temperature response function to apply to alphaCx
            f_tmp(:,i) = ((tmp_ave(:) - Tmin(i)) / (Topt(i) - Tmin(i))) * &
                ((Tmax(i) - tmp_ave(:)) / (Tmax(i) - Topt(i))) ** ((Tmax(i) - Topt(i)) / (Topt(i) - Tmin(i)))

            where( tmp_ave(:) <= Tmin(i) .or. tmp_ave(:) >= Tmax(i) )
                f_tmp(:,i) = 0.d0
            end where

            ! calculate temperature response function to apply to gc (uses mean of Tx and Tav instead of Tav, Feikema et al 2010)
            f_tmp_gc(:,i) = (((tmp_ave(:) + tmp_max(:)) / 2 - Tmin(i)) / (Topt(i) - Tmin(i))) * &
                ((Tmax(i) - (tmp_ave(:) + tmp_max(:)) / 2) / (Tmax(i) - Topt(i))) ** ((Tmax(i) - Topt(i)) / (Topt(i) - Tmin(i)))

            where( (tmp_ave(:) + tmp_max(:)) / 2 <= Tmin(i) .or. (tmp_ave(:) + tmp_max(:)) / 2 >= Tmax(i) )
                f_tmp_gc(:,i) = 0.d0
            end where

            ! frost modifier
            f_frost(:,i) = 1.d0 - kF(i) * ( frost_days(:) / 30.d0)

            ! CO2 modifiers
            f_calpha(:,i) = fCalphax(i) * co2(:) / (350.d0 * (fCalphax(i) - 1.d0) + co2(:))
            f_cg(:,i) = fCg0(i) / (1.d0 + (fCg0(i) - 1.d0) * co2(:) / 350.d0)

        end do


        ! air pressure
        air_pressure = 101.3d0 * Exp(-1.d0 * altitude / 8200.d0)


        ! SOIL WATER --------
        ! Assign the SWconst and SWpower parameters for this soil class
        if ( soil_class > 0.d0 ) then
            ! Standard soil type
            SWconst(:) = 0.8d0 - 0.10d0 * soil_class
            SWpower(:) = 11.d0 - 2.d0 * soil_class
        elseIf ( soil_class < 0.d0 ) then
            ! Use supplied parameters
            SWconst(:) = SWconst0(:)
            SWpower(:) = SWpower0(:)
        else
            ! No soil-water effects
            SWconst(:) = 999
            SWpower(:) = SWpower0(:)
        end if

        ! Initial ASW must be between min and max ASW
        if (asw_min > asw_max) then
            asw_min = asw_max
        end if

        ASW = max( min( ASW, asw_max ), asw_min )

        ! Silvicultural events are currently not active
        Irrig = 0.d0
        water_runoff_polled = 0.d0
        poolFractn = 0.d0
        poolFractn = max(0.d0, min(1.d0, poolFractn))

        ! Snow precipitation and snowmelt from MWBM
		where ( tmp_ave(:) <= -1.d0 )
			 snow_prcp(:) = prcp(:)
		end where
		where ( tmp_ave(:) >= 3.d0 )
			 snow_prcp(:) = 0.d0
		end where
		where ( tmp_ave(:) > -1.d0 .and. tmp_ave(:) < 3.d0)
			 snow_prcp(:) = prcp(:) * ( 3.d0 - tmp_ave(:) ) / 4.d0
		end where

		prcp(:) = prcp(:) - snow_prcp(:)

		snow_water = 0.d0

		where ( prcp(:) < 0.d0 )
			 prcp(:) = 0.d0
		end where

        ! NUTRITIONS --------
        ! Check fN(FR) for no effect: fNn = 0 ==> fN(FR)=1 for all FR
        where( fNn(:) == 0.d0 ) fN0(:) = 1.d0


        ! Partitioning  --------
        pfsPower(:) = Log( pFS20(:) / pFS2(:) ) / Log( 20.d0 / 2.d0 )
        pfsConst(:) = pFS2(:) / 2.d0 ** pfsPower(:)



        ! INITIALISATION (Age dependent)---------------------
        ! Calculate the species specific modifiers
        do i = 1, n_sp
            age(:,i) = 12.d0 * ( year_i - year_p(i) ) + month_i - month_p(i) - 1.d0 !
            age(:,i) =  ( age(:,i) + int( (/(i, i=1, n_m)/) ) ) / 12.d0 ! translate to years
            age_m(:,i) =  age(:,i) - 1.d0/12.d0
            age_m(1,i) =  age(1,i)

            SLA(:,i) = f_exp( n_m, age_m(:,i), SLA0(i), SLA1(i), tSLA(i), 2.d0)
            fracBB(:,i) = f_exp( n_m, age_m(:,i), fracBB0(i), fracBB1(i), tBB(i), 1.d0)
            wood_density(:,i) = f_exp( n_m, age_m(:,i), rho0(i), rho1(i), tRho(i), 1.d0)
            gammaN(:,i) = f_exp( n_m, age(:,i), gammaN0(i), gammaN1(i), tgammaN(i), ngammaN(i))
			
			if (maint_resp .eq. int(2)) then
				rsap(:,i) = f_exp( n_m, age(:,i), pssN0(i), pssN1(i), tpssN(i), npssN(i))			
				where ( rsap(:,i) > 1.d0 )
					rsap(:,i) = 1.d0
				end where
			end if 
			
            gammaF(:,i) = f_exp_foliage( n_m, age_m(:,i), gammaF1(i), gammaF0(i), tgammaF(i))


            ! age modifier
            if (nAge(i) == 0.d0) then
                f_age(:,i) = 1.d0
            else
                ! I'm not declaring relative age, but directly put it inside
                f_age(:,i) = 1.d0 / (1.d0 + ( (age_m(:,i) / MaxAge(i) ) / rAge(i)) ** nAge(i))
            end if

            !age(:,i) = age(:,i) + 1.d0 ! that how the VBA works

        end do


        ! INITIALISATION (Stand)---------------------
        ii = 1
        month = month_i
		
		if (soil_model .gt. int(0)) then
			tmp_year = sum(tmp_ave(max((ii - month + 1),1):(ii + 12 - month)))/12.d0
			tmp_yasso_year(:) = tmp_ave(max((ii - month + 1),1):(ii + 12 - month))
			prcp_year = sum(prcp(max((ii - month + 1),1):(ii + 12 - month)))
		end if

        where (age(ii,:) >= 0.d0 )
          stems_n(:) = stems_n_i(:)
          biom_stem(:) = biom_stem_i(:)
          biom_foliage(:) = biom_foliage_i(:)
          biom_root(:) = biom_root_i(:)
          n_sapling(:) = n_sapling_i(:)
          sapling_wf(:) = sapling_wf_i(:)
          sapling_wr(:) = sapling_wr_i(:)
          sapling_ws(:) = sapling_ws_i(:)
        end where

        ! Initial soil characteristics

        O_C = soil_carbon_i
        O_N = soil_nitrogen_i

		! Initialize soil accounting variables
		TotalCarbo = soil_carbon_i + sum(litter_carbon_i(:)) + sum(deadwood_carbon_i(:))
		TotalNitro = soil_nitrogen_i + sum(litter_nitrogen_i(:)) + sum(deadwood_nitrogen_i(:))

		! Initial soil characteristics for Yasso
		!folAWENH(:)  = 0.d0
		!crbAWENH(:)  = 0.d0
		!stAWENH(:) = 0.d0
		!Yb_input(:) = 0.d0
		!Yr_input(:) = 0.d0
		!Yl_input(:) = 0.d0
		tstep = 1.d0
		size_eff = 0.d0
		leaching = 0.d0

		do i = 1,4
			soilC(ii,:,:,i) = 0.d0
		end do

		soilC(ii,:,1:3,5) =  0.d0
		
		litter_input = 0.d0

        do i = 1, n_sp
			Yr_C(i) = deadwood_carbon_i(i)
			Yr_N(i) = deadwood_nitrogen_i(i)
			Yl_C(i) = litter_carbon_i(i)
			Yl_N(i) = litter_nitrogen_i(i)

			!Initialize leaf and root litter
			litter_input = Yl_C(i)
			call compAWENH(litter_input,folAWENH,parameters_AWEN(1:5,i))
			soilC(ii,i,1,:) = folAWENH(:)
			
			!Initialize branches and coarse roots
			litter_input = Yr_C(i) * 0.2d0
			call compAWENH(litter_input,crbAWENH,iniAWEN)
			soilC(ii,i,2,:) = crbAWENH(:)

			! Initialize CWD
			litter_input = Yr_C(i) * 0.8d0
			call compAWENH(litter_input,stAWENH,iniAWEN)
			soilC(ii,i,3,:) = stAWENH(:)
		end do

		! Initialize SOC
		soilC(ii,1,4,5) =  soil_carbon_i

		! Initialize soil temperature
		call soil_temperature_toy(t_soil_surf, tmp_ave(ii), month, sum(lai(:)), snow_water)
		
        ! Check if this is the dormant period or previous/following period is dormant
        ! to allocate foliage if needed, etc.
        do i = 1, n_sp
            ! if this is a dormant month
            if( f_dormant(month, leafgrow(i), leaffall(i)) .eqv. .TRUE. ) then
                biom_foliage_debt(i)= biom_foliage(i)
                biom_foliage(i) = 0.d0
            end if
        end do

        ! Initial stand characteristics
        where (age(ii,:) >= 0.d0 )
          biom_tree(:) = biom_stem(:) * 1000.d0 / stems_n(:)  ! kg/tree
          dbh(:) = ( biom_tree(:) / aWs(:)) ** (1.d0 / nWs(:))
          basal_area(:) = dbh(:) ** 2.d0 / 4.d0 * Pi * stems_n(:) / 10000.d0
          lai(:) =  biom_foliage(:) * SLA(ii,:) * 0.1d0
        end where

        competition_total(:) = sum( wood_density(ii,:) * basal_area(:) )
        
        if( height_model .eq. 1 ) then
            height(:) = aH(:) * dbh(:) ** nHB(:) * competition_total(:) ** nHC(:)
        else if ( height_model .eq. 2 ) then
            height(:) = 1.3d0 + aH(:) * Exp(1.d0)**(-nHB(:)/dbh(:)) + nHC(:) * competition_total(:) * dbh(:)
        end if

        ! Correct the bias
        do n = 1, b_n
            competition_total(:) = sum( wood_density(ii,:) * basal_area(:) )

            call s_sizeDist_correct(n_sp, age(ii,:), stems_n(:), biom_tree(:), competition_total(:), lai(:), &
                correct_bias, height_model,  pars_i(69:85,:), pars_b, aWs(:), nWs(:), pfsPower(:), pfsConst(:), &
                dbh(:), basal_area(:), height(:), crown_length(:), crown_width(:), pFS(:), bias_scale(:,:) )
        end do
        
        Height_max = maxval( height(:), mask=lai(:)>0.d0 )

        ! Volume and Volume increment
        volume(:) = biom_stem(:) * (1.d0 - fracBB(ii,:)) / wood_density(ii,:)

        where( aV(:) > 0 ) volume(:) = aV(:) * dbh(:) ** nVB(:) * height(:) ** nVH(:) * &
                (dbh(:) * dbh(:) * height(:)) ** nVBH(:) * stems_n(:)

        volume_cum(:) = volume(:)
        volume_old(:) = volume(:)
		harvesting(:) = 0.d0
		residues(:) = 0.d0
		root_residues(:) = 0.d0
		foliage_residues(:) = 0.d0
		
		where(cwd_removal(:) >= 1.d0) cwd_removal(:) = 1.d0

        volume_mai(:) = volume_cum(:) / age(ii,:)

        basal_area_prop(:) = basal_area(:) / sum( basal_area(:) )

        ! crop trees
        ! Check if crop trees are defined
        do i=1, n_sp
            if (managementInputs(t_n(i),8,i) > 0.d0) then
                crop_trees(i) = NINT(managementInputs(t_n(i),8,i))
            end if
        end do 

        if (sum(crop_trees(:)) > 0.d0) then
            
            call    get_crop_trees(dbh_crop(:), vol_crop(:), basal_area_crop(:), biom_stem_crop(:), height_crop(:), &
                                 dbh(:), age(ii,:), height(:), stems_n(:), basal_area(:), volume(:), &
                                 volume_cum(:) * wood_density(ii,:), &
                                 competition_total(:), aWs(:), nWs(:), aH(:), nHB(:), nHC(:), &
                                 aV(:), nVB(:), nVH(:), nVBH(:), Dscale0(:), DscaleB(:), Dscalerh(:), Dscalet(:), &
                                 DscaleC(:), Dshape0(:), DshapeB(:), Dshaperh(:), Dshapet(:), DshapeC(:), Dlocation0(:), &
                                 DlocationB(:), Dlocationrh(:), Dlocationt(:), DlocationC(:), crop_trees(:), fracBB(ii,:), &
                                 wood_density(ii,:), height_model, 5.d0, n_sp) 
        end if

		! Base phenology
		leaffall_base(:) = leaffall(:)
		leafgrow_base(:) = leafgrow(:)

        ! INITIALISATION (Write output)---------------------
        include 'i_write_out.h'



        !*************************************************************************************
        ! Monthly simulations

        do ii = 2, n_m
		
            ! month update
            month = month + int(1)

            if (month > 12) then
                month = int(1)
                if (soil_model .gt. int(0)) then
                    tmp_year = sum(tmp_ave(max((ii - month + 1),2):(ii + 12 - month)))/12.d0
                    ! Climate for Yasso modifiers
                    if (soil_model == int(2)) then
                        tmp_yasso_year(:) = tmp_ave(max((ii - month + 1),2):(ii + 12 - month))
                        prcp_year = sum(prcp(max((ii - month + 1),2):(ii + 12 - month)))
                    end if
                end if
            end if
			
			
            ! Add new cohort ----------------------------------------------------------------------
            where (age(ii,:) .eq. 0.d0 )
              stems_n(:) = stems_n_i(:)
              biom_stem(:) = biom_stem_i(:)
              biom_foliage(:) = biom_foliage_i(:)
              biom_root(:) = biom_root_i(:)
            end where

            if( any(age(ii,:) .eq. 0.d0) ) then
              b_cor = .TRUE.
            end if


            ! Test for dormancy ----------------------------------------------------------------------

            ! If this is first month after dormancy we need to make potential LAI, so the
            ! PAR absorbption can be applied, otherwise it will be zero.
            ! In the end of the month we will re-calculate it based on the actual values
            do i = 1, n_sp
                if( f_dormant(month, leafgrow(i), leaffall(i)) .eqv. .FALSE. ) then
                    if( f_dormant(month-1, leafgrow(i), leaffall(i)) .eqv. .TRUE. ) then
                        lai(i) = biom_foliage_debt(i) * SLA(ii,i) * 0.1d0
                        b_cor = .TRUE.
                    end if
                end if

                ! If this is first dormant month, then we set WF to 0 and move everything to the dept
                if( f_dormant(month, leafgrow(i), leaffall(i)) .eqv. .TRUE. ) then
                    if( f_dormant(month-1, leafgrow(i), leaffall(i)) .eqv. .FALSE. ) then
                        biom_foliage_debt(i)= biom_foliage(i)
                        biom_foliage(i) = 0.d0
                        lai(i) =  0.d0
                        b_cor = .TRUE.
                    end if
                end if

            end do

            !****** We shall call this only if the any of the above is TRUE
            if ( b_cor .eqv. .TRUE. ) then
                do n = 1, b_n
                    competition_total(:) = sum( wood_density(ii,:) * basal_area(:) )

                    call s_sizeDist_correct(n_sp, age(ii,:), stems_n(:), biom_tree(:), competition_total(:), lai(:), &
                        correct_bias, height_model,  pars_i(69:85,:), pars_b, aWs(:), nWs(:), pfsPower(:), pfsConst(:), &
                        dbh(:), basal_area(:), height(:), crown_length(:), crown_width(:), pFS(:), bias_scale(:,:) )
                end do
                b_cor = .FALSE.
            end if

            !Radiation and assimilation ----------------------------------------------------------------------
            if ( light_model .eq. int(1) .or. age(ii,1) < fullCanAge(1)) then
                call s_light_3pgpjs ( n_sp, age_m(ii,:), fullCanAge(:), k(:), lai(:), &
                    solar_rad(ii), daysInMonth(month), &
                    canopy_cover(:), apar(:) )

                VPD_sp(:) = vpd_day(ii)

            else if ( light_model .eq. int(2) .and. age(ii,1) >= fullCanAge(1)) then

                ! Calculate the absorbed PAR. If this is first month, then it will be only potential
                call s_light_3pgmix ( n_sp, height(:), crown_length(:), crown_width(:), lai(:), stems_n(:), &
                    solar_rad(ii), CrownShape(:), k(:), adjSolarZenithAngle(month), daysInMonth(month), &
                    apar(:), lai_above(:), fi(:), lambda_v(:), lambda_h(:), canopy_vol_frac(:), layer_id(:), lai_sa_ratio(:))

                VPD_sp(:) = vpd_day(ii) * Exp(lai_above(:) * (-Log(2.d0)) / cVPD(:))

            end if


            ! Determine the various environmental modifiers which were not calculated before
            ! calculate VPD modifier
            ! Get within-canopy climatic conditions this is exponential function
            Height_max = maxval( height(:), mask=lai(:)>0.d0 )

            ! but since BLcond is a vector we can't use the expF
            aero_resist(:) = (1.d0 / BLcond(:)) + (5.d0 * sum( lai(:) ) - (1.d0 / BLcond(:))) * &
                Exp(-ln2 * ( height(:) / (Height_max / 2.d0)) ** 2.d0)
            ! if this is the highest tree
            where( height(:) == Height_max)
                aero_resist(:) = 1.d0 / BLcond(:)
            end where
            ! Check for dormancy
            where( lai(:) .eq. 0.d0)
                aero_resist(:) = 0.d0
            end where

            f_vpd(:) = Exp( -CoeffCond(:) * VPD_sp(:))

            ! soil water modifier
            f_sw(:) = 1.d0 / (1.d0 + ((1.d0 -  ASW / asw_max) / SWconst(:)) ** SWpower(:))

            ! soil nutrition modifier
            f_nutr(:) = 1.d0 - (1.d0 - fN0(:)) * (1.d0 - fertility(:)) ** fNn(:)
            where( fNn(:) == 0.d0 ) f_nutr(:) = 1.d0


            !Inundation modifier 
            !icount is the number of month, the WTP was above the species tolerance
            !level. Compare icount with inund_dur, i.e. the number of months,
            !the species can tolerate inundation.
            ! if (anoxic_io .eq. int(1)) then
            !     do i = 1, n_sp
            !         where (icount(:,i) .LT. inund_dur(i))
            !             f_inund(:,i) = (icount(:,i) ) / inund_dur(i)
            !         end where
            !     end do
            ! end if
			
            ! calculate physiological modifier applied to conductance and alphaCx.
            if ( phys_model .eq. int(1) ) then

                f_phys(:) = min( f_vpd(:), f_sw(:) ) * f_age(ii,:)
                f_tmp_gc(ii,:) = 1.d0

            else if ( phys_model .eq. int(2) ) then

                f_phys(:) = f_vpd(:) * f_sw(:) * f_age(ii,:)

            end if

			if (gpp_model .eq. int(2)) then
				alpha_c(:) = alphaCx(:) * f_nutr(:) * f_tmp(ii,:) * f_frost(ii,:) * f_sw(:) * f_age(ii,:)
				!alpha_c(:) = alphaCx(:) * f_nutr(:) * f_tmp(ii,:) * f_frost(ii,:) * f_phys(:)
			else if (gpp_model .eq. int(3)) then
				alpha_c(:) = alphaCx(:) * f_nutr(:) * f_tmp(ii,:) * f_frost(ii,:) * f_age(ii,:)
			else
				alpha_c(:) = alphaCx(:) * f_nutr(:) * f_tmp(ii,:) * f_frost(ii,:) * f_calpha(ii,:) * f_phys(:)
			end if

            ! if (anoxic_io .eq. int(1)) then
            !     alpha_c(:) = alpha_c(:) * ( 1.d0 - f_inund(ii,:))
            ! end if


			! Call rpmodel for GPP computation
			if (gpp_model .eq. int(2)) then
				do i = 1, n_sp
					out_aux(:) = 0.d0
					apar(i) = apar(i) * 2.30d0 ! Convert to mol/m2

					if ( apar(i) > 0.d0 ) then
						call rpmodel( tmp_ave(ii), VPD_sp(i)*100.d0, co2(ii), 1.d0, apar(i) , -1.d0, &
							 altitude, alpha_c(i), 146.d0, f_sw(i), 1.d0, 0.d0, 0.733d0, .FALSE., &
							 int(method_jmaxlim), .FALSE., .FALSE., out_aux)						 							 							 							 
					end if

					 out_rpmodel(i,:) = out_aux(:)
				end do
            ! Phydro GPP model
			else if (gpp_model .eq. int(3)) then
				do i = 1, n_sp
					out_aux_h(:) = 0.d0
					apar(i) = apar(i) * molPAR_MJ(i) ! Convert to mol/m2 according to internal parameters
					
                    ! Soil water potential
					call calc_psi_soil(asw, asw_max, soil_class, psi_soil)

					if ( apar(i) > 0.d0 ) then
							 
						call rphmodel( tmp_ave(ii), apar(i) * 11.57407d0 / daysInMonth(month) , VPD_sp(i)*100.d0, co2(ii), altitude, 1.d0, &
						  	 alpha_c(i), psi_soil, 0.02d0, par_plant, par_cost, out_aux_h)							 
							 							 							 
					end if

					 out_rphmodel(i,:) = out_aux_h(:)
				end do
            end if

            ! Calculate assimilation before the water ballance is done

            where( lai(:) == 0.d0 ) alpha_c(:) = 0.d0
            epsilon(:) = gDM_mol * molPAR_MJ * alpha_c(:)

			if (gpp_model .eq. int(2)) then
				GPP(:) = out_rpmodel(:,1) / 100.d0 * (1.d0 / dmC) ! convert rpmodel output to tDW/ha
			else if (gpp_model .eq. int(3)) then
				GPP(:) = out_rphmodel(:,1) * daysInMonth(month) / 100.d0 * (1.d0 / dmC) ! convert rp hydro model output to tDW/ha				
			else
				GPP(:) = epsilon(:) * apar(:) / 100        ! tDM/ha (apar is MJ/m^2)
			end if

			if (maint_resp .eq. int(2)) then

				
				!soil temperature based on swat			
				call soil_temperature_swat(t_soil, tmp_ave(ii), tmp_year, ASW, snow_water, &
										sum(lai(:)), month, soil_class)


            ! Get autotrophic respiration
                call   get_raut(ra(:), rsap(ii,:), pssN0(:), pssN1(:), tpssN(:), npssN(:), age(ii,:), dbh(:), height(:), &
                            basal_area(:), crown_length(:), stems_n(:), f_phys(:), m(:), biom_root(:), biom_foliage(:), &
                            biom_foliage_debt(:), biom_fineroot(:), biom_stem(:), biom_sapwood(:), fr_ratio(:), Ncr(:), &
                            Ncs(:), out_rpmodel(:,22), out_rphmodel(:,4), t_soil, tmp_ave(ii), month, gpp_model, n_sp )                                     
						
                !Get NPP
				NPP(:) = (1.d0 - 0.25d0) * (GPP(:) - ra(:)) 
			else
				NPP(:) = GPP(:) * y(:) ! assumes respiratory rate is constant
			end if
	

			where(NPP(:) < 0.d0)
				NPP(:) = 0.d0
			end where

            ! Water Balance ----------------------------------------------------------------------
            ! Calculate each species proportion
            lai_total(:) = sum( lai(:) )
            lai_per(:) = lai(:) / lai_total(:)
            where( lai_total(:) .eq. 0.d0 ) lai_per(:) = 0.d0

            ! Calculate conductance
            gC(:) = MaxCond(:)
            where( lai_total(:) <= LAIgcx(:) )
                gC(:) = MinCond(:) + (MaxCond(:) - MinCond(:)) * lai_total(:) / LAIgcx(:)
            end where

			if (gpp_model .eq. int(2)) then

                ! Get canopy conductance
				if (canopy_cond .eq. int(2)) then

					! Canopy conductance based on the upscaling of the P-model output
					conduct_canopy(:) = out_rpmodel(:,11) * (1.d0 - exp( - k(:) * lai(:))) / k(:)  / & 
						(day_length(month) * daysInMonth(month)) 

				else
                    ! Default conductance
					conduct_canopy(:) =  gC(:) * lai_per(:) * f_phys(:) * f_tmp_gc(ii,:) * f_cg(ii,:)
				end if

				where ( lai(:) == 0.d0 )
					 conduct_canopy(:) = 0.d0
				end where

			else if (gpp_model .eq. int(3)) then
			
				if (canopy_cond .eq. int(2)) then	
                    
                    ! Canopy conductance based on the upscaling of the Phydro model output
					conduct_canopy(:) = out_rphmodel(:,3) * (1.d0 - exp( - k(:) * lai(:))) / k(:)  

				else
                    ! Default conductance
					conduct_canopy(:) =  gC(:) * lai_per(:) * f_phys(:) * f_tmp_gc(ii,:) * f_cg(ii,:)
				end if

				where ( lai(:) == 0.d0 )
					 conduct_canopy(:) = 0.d0
				end where
			else
				conduct_canopy(:) = gC(:) * lai_per(:) * f_phys(:) * f_tmp_gc(ii,:) * f_cg(ii,:)
			end if

            conduct_soil = MaxSoilCond * ASW / asw_max


            ! Calculate transpiration
            if ( transp_model .eq. int(1) ) then

                call s_transpiration_3pgpjs( n_sp, solar_rad(ii), day_length(month), VPD_sp(:), BLcond(:), &
                    conduct_canopy(:), daysInMonth(month), Qa, Qb, &
                    transp_veg(:))
                evapotra_soil = 0.d0

            else if ( transp_model .eq. int(2) ) then

                call s_transpiration_3pgmix( n_sp, solar_rad(ii), vpd_day(ii), day_length(month), daysInMonth(month), &
                    lai(:), fi(:), VPD_sp(:), aero_resist(:), conduct_canopy(:), conduct_soil, Qa, Qb, &
                    transp_veg(:), evapotra_soil)
					
			else if (gpp_model .eq. int(3) .and. transp_model .eq. int(3) ) then
			
				transp_veg(:) = out_rphmodel(:,5) * daysInMonth(month)
				evapotra_soil = 0.d0
				call soil_evap( solar_rad(ii), vpd_day(ii), day_length(month), daysInMonth(month), &
				            	sum(lai(:)), sum(fi(:)), conduct_soil, Qa(1), Qb(1), evapotra_soil)
            end if
			

            transp_total = sum( transp_veg(:) ) + evapotra_soil


            ! rainfall interception
            prcp_interc_fract(:) = MaxIntcptn(:)
            where (LAImaxIntcptn(:) > 0.d0)
                prcp_interc_fract(:) = MaxIntcptn(:) * min(1.d0, lai_total(:) / LAImaxIntcptn(:)) * LAI_per(:)
            end where

            prcp_interc(:) = prcp(ii) * prcp_interc_fract(:)
            prcp_interc_total = sum( prcp_interc(:) )

		   ! Do snow water balance
			snow_water = snow_water + snow_prcp(ii)
			if (tmp_ave(ii) > 5.5) then
				snowmelt = snow_water
			else
				snowmelt = min( snow_water * (0.5 * (tmp_ave(ii) + 1.d0) / 4.d0 ), snow_water)
				if ( snowmelt < 0.d0 ) then
					snowmelt = 0.d0
				end if
			end if
			snow_water = snow_water - snowmelt

            ! Do soil water balance Need to constrain irrigation only to the growing season
            ASW = ASW + prcp(ii) + (100.d0 * Irrig / 12.0d0) + water_runoff_polled + snowmelt
            evapo_transp = min( ASW, transp_total + prcp_interc_total)  !ET can not exceed ASW
            excessSW = max(ASW - evapo_transp - asw_max, 0.d0)
            ASW = ASW - evapo_transp - excessSW
            water_runoff_polled = poolFractn * excessSW
            prcp_runoff = (1.d0 - poolFractn) * excessSW
			ASW = max(min(ASW , asw_max), asw_min)

            ! if (peat_io .eq. int(1)) then
            !     !GRADIENT, with which water content decreases with WTP
            !     !Granberg used a gradient for the peat
            !     z1 = (acro_por - surfsw) / grad_lim
            !     z2 = (surfsw - minvtot) / (zmin - grad_lim)
            !     az = (acro_por - minvtot) / zmin
                
            !     !Estimating the water over surface 
            !     if (ASW .GT. z1 * acro_por) then 
            !         wtp_d = ASW - z1 * acro_por
            !         if (wtp_d .GT. maxhei) wtp_d = maxhei
            !     else 
            !         wtd = SQRT(3 * (acro_por * z1 - ASW)/ (2 * az))
            !         if(wtd .GT. zmin) then
            !             wtd = 3 * (acro_por * z1 - ASW) / (2 * (acro_por - minvtot))
            !         end if 
            !         wtp_d = - wtd
            !     end if 
            ! end if 

            !INUNDATION STRESS
            !Here we take a amount of water standing and comapre to the 
            !stress tolerance level of each species, and number of months 
            ! do i = 1, n_sp
            !     wtp = wtp_d
            !     if (wtp .GT. inund_height(i)) then
            !         icount(:,i) = icount(:,i) + 1.d0
            !     end if
            ! end do    

            if (ASW < asw_min) then
                irrig_supl = asw_min - ASW
                ASW = asw_min
            end if


            if ( ( transp_total + prcp_interc_total ) == 0 ) then
                !this might be close to 0 if the only existing species is dormant during this month
                ! (it will include the soil evaporation if Apply3PGpjswaterbalance = no)
                f_transp_scale = 1.
            else
                f_transp_scale = evapo_transp / (transp_total + prcp_interc_total)  !scales NPP and GPP
            end if

            ! correct for actual ET
            GPP = GPP * f_transp_scale
            NPP = NPP * f_transp_scale
            NPP_f = NPP


            if ( transp_total > 0 .and. f_transp_scale < 1 ) then
                ! a different scaler is required for transpiration because all of the scaling needs
                ! to be done to the transpiration and not to the RainIntcpth, which occurs regardless of the growth
                transp_veg(:) = (evapo_transp - prcp_interc_total) / transp_total * transp_veg(:)
                evapotra_soil = (evapo_transp - prcp_interc_total) / transp_total * evapotra_soil
            end if


            ! NEED TO CROSS CHECK THIS PART, DON'T FULLY AGREE WITH IT
            if ( evapo_transp /= 0.d0 .and. n_sp == 1 ) then
                ! in case ET is zero! Also, for mixtures it is not possible to calculate WUE based on
                ! ET because the soil evaporation cannot simply be divided between species.
                WUE(:) = 100.d0 * NPP(:) / evapo_transp
            else
                WUE(:) = 0.d0
            end if

            WUE_transp(:) = 0.d0
            where ( transp_veg(:) > 0.d0 )
                WUE_transp(:) = 100.d0 * NPP(:) / transp_veg(:)
            end where


            if ( calculate_d13c .eq. int(1) ) then
                ! d13C module ----------------------------------------------------------------------
                ! Calculating d13C - This is based on Wei et al. 2014 (Plant, Cell and Environment 37, 82-100)
                ! and Wei et al. 2014 (Forest Ecology and Management 313, 69-82). This is simply calculated from
                ! other variables and has no influence on any processes

                !convert GPP (currently in tDM/ha/month) to GPP in mol/m2/s.
                GPP_molsec(:) = GPP(:) * 100.d0 / ( daysInMonth(month) * 24.0d0 * 3600.0d0 * gDM_mol)

                !Canopy conductance for water vapour in mol/m2s, unit conversion (CanCond is m/s)
                Gw_mol(:) = conduct_canopy(:) * 44.6d0 * (273.15d0 / (273.15d0 + tmp_ave(ii) ) ) * (air_pressure / 101.3d0)

                ! Canopy conductance for CO2 in mol/m2s
                ! This calculation needs to consider the area covered by leaves as opposed to the total ground area of the stand.
                ! The explanation that Wei et al. provided for adding the "/Maximum(0.0000001, CanCover)" is
                ! that 3PG is a big leaf leaf model for conductance and the leaf area is assumed to be evenly distributed
                ! across the land area. So GwMol is divided by Maximum(0.0000001, CanCover) to convert the conductance
                ! to the area covered by the leaves only, which is smaller than the land area if the canopy has not
                ! closed. If the original light model has been selected then a CanCover value has already been calculated
                ! although Wei et al. also warn against using d13C calculations in stands with CanCover < 1.
                ! If the new light model has been selected then CanCover still needs to be calculated.

                canopy_cover(:) = (stems_n(:) * (crown_width(:) + 0.25d0) ** 2.d0) / 10000.d0
				
                where( canopy_cover(:) > 1.d0) canopy_cover(:) = 1.d0

                Gc_mol(:) = Gw_mol(:) * RGcGW(:) / max(0.0000001d0, canopy_cover(:))

                !Calculating monthly average intercellular CO2 concentration. Ci = Ca - A/g
                InterCi(:) = CO2(ii) * 0.000001d0 - GPP_molsec(:) / Gc_mol(:)

                !Calculating monthly d13C of new photosynthate, = d13Catm- a-(b-a) (ci/ca)
                D13CNewPS(:) = d13Catm(ii) - aFracDiffu(:) - (bFracRubi(:) - aFracDiffu(:)) * (InterCi(:) / (CO2(ii) * 0.000001d0))
                D13CTissue(:) = D13CNewPS(:) + D13CTissueDif(:)

                ! correct for dormancy
                where( Gc_mol(:) .eq. 0.d0 )
                    InterCi(:) = 0.d0
                    D13CNewPS(:) = 0.d0
                    D13CTissue(:) = 0.d0
                end where

            end if


            ! Biomass increment and loss module ----------------------------------------------
            ! determine biomass increments and losses

            ! Get allocation coeficients
            call allocation_3pg(m0(:), fertility(:), pRx(:), pRn(:), f_phys(:), pFS(:), m(:), &
                                npp_fract_root(:), npp_fract_stem(:), npp_fract_foliage(:), n_sp) 


            do i = 1, n_sp

                !  Dormant period -----------
                if ( f_dormant(month, leafgrow(i), leaffall(i)) .eqv. .TRUE. ) then

                    ! There is no increment. But if this is the first dormant period then there is litterfall
                    if ( f_dormant(month-1, leafgrow(i), leaffall(i))  .eqv. .TRUE. ) then
                        biom_loss_foliage(i) = 0.d0
                    else
                        biom_loss_foliage(i) = biom_foliage_debt(i)
                    end if

                    biom_loss_root(i) = 0.d0


                    ! No changes during dormant period
                    biom_incr_foliage(i) = 0.d0
                    biom_incr_root(i) = 0.d0
                    biom_incr_stem(i) = 0.d0

                else

                    ! if there is some leaves to be growth put first NPP to the leaf growth
                    ! if there is enough NPP then growth all the leaves, otherwise wait for next period
                    if( biom_foliage(i) == 0.d0 ) then
                        biom_foliage(i) = biom_foliage_debt(i)
                    end if

                    if( NPP(i) >= biom_foliage_debt(i) ) then
                        !if there is enough NPP
                        NPP(i) = NPP(i) - biom_foliage_debt(i)
                        biom_foliage_debt(i) = 0.d0
                    else
                        ! IF there is not enough NPP to regrow the leaves we regrow part and wait for
                        biom_foliage_debt(i) = biom_foliage_debt(i) - NPP(i)
                        NPP(i) = 0.d0
                    end if

                    ! Calculate biomass loss
                    biom_loss_foliage(i) = gammaF(ii, i) * biom_foliage(i)
                    biom_loss_root(i) = gammaR(i) * biom_root(i)


                    ! Calculate biomass increments
                    biom_incr_foliage(i) = NPP(i) * npp_fract_foliage(i)
                    biom_incr_root(i) = NPP(i) * npp_fract_root(i)
                    biom_incr_stem(i) = NPP(i) * npp_fract_stem(i)


                    ! Calculate end-of-month biomass
                    biom_foliage(i) = biom_foliage(i) + biom_incr_foliage(i) - biom_loss_foliage(i)
                    biom_root(i) = biom_root(i) + biom_incr_root(i) - biom_loss_root(i)
                    biom_stem(i) = biom_stem(i) + biom_incr_stem(i)

                end if

            end do

            ! Correct the bias
            biom_tree(:) = biom_stem(:) * 1000.d0 / stems_n(:)  ! kg/tree
            where( stems_n(:) .eq. 0.d0 ) biom_tree(:) = 0.d0

            lai(:) =  biom_foliage(:) * SLA(ii,:) * 0.1d0

            do n = 1, b_n
                competition_total(:) = sum( wood_density(ii,:) * basal_area(:) )

                call s_sizeDist_correct(n_sp, age(ii,:), stems_n(:), biom_tree(:), competition_total(:), lai(:), &
                    correct_bias, height_model,  pars_i(69:85,:), pars_b, aWs(:), nWs(:), pfsPower(:), pfsConst(:), &
                    dbh(:), basal_area(:), height(:), crown_length(:), crown_width(:), pFS(:), bias_scale(:,:) )
            end do

            ! Volume and Volume increment
            ! This is done before thinning and mortality part
            volume(:) = biom_stem(:) * (1.d0 - fracBB(ii,:)) / wood_density(ii,:)
            where( aV(:) > 0.d0 ) volume(:) = aV(:) * dbh(:) ** nVB(:) * height(:) ** nVH(:) * &
                    (dbh(:) * dbh(:) * height(:)) ** nVBH(:) * stems_n(:)

            volume_change(:) = volume(:) - volume_old(:)
            where( lai(:) .eq. 0.d0 ) volume_change(:) = 0.d0
            where( volume_change(:) .le. 0.d0 ) volume_change(:) = 0.d0

            volume_cum(:) = volume_cum(:) + volume_change(:)
            volume_old(:) = volume(:)

            volume_mai(:) = volume_cum(:) / age(ii,:)


            ! Management -------------------------------------------------------------------------

		   ! Initialize accounting variables
			harvesting(:) = 0.d0
            dbh_harv(:) = 0.d0
            height_harv(:) = 0.d0
            stem_harv(:) = 0.d0
            biom_stem_harv(:) = 0.d0
			residues(:) = 0.d0
			root_residues(:) = 0.d0
			foliage_residues(:) = 0.d0
            if (wind_dist .eq. int(1)) lag(:) = lag(:) + 1.d0/12.d0

			do i = 1, n_sp
			
				target_ba = 0.d0
				target_vol = 0.d0
				removal = 0.d0
				mort_manag(i) = 0.d0
                final_harvest(i) = .FALSE.

                if( t_t(i) > 0 ) then

                    if(t_n(i) <= t_t(i)) then

                        if( age(ii,i) >= managementInputs(t_n(i),1,i) ) then

                            if (wind_dist .eq. int(1)) then
                                lag(i) = 0.d0
                                previous_density(i) = stems_n(i)
                            end if

                            ! Check if crop trees are defined
                            if (managementInputs(t_n(i),8,i) > 0.d0) then
                                crop_trees(i) = NINT(managementInputs(t_n(i),8,i))
                            end if

							! Check if thinning is to be performed based on basal area removal (2nd priority)
							if( managementInputs(t_n(i),6,i) >= 0.d0 .and. managementInputs(t_n(i),7,i) < 0.d0 ) then

								! Check if target basal area is in percentage or absolute value
								if ( managementInputs(t_n(i),6,i) < 1.d0) then
									target_ba = basal_area(i)*( 1.d0 - managementInputs(t_n(i),6,i) )
								else if ( managementInputs(t_n(i),6,i) == 1.d0) then
									target_ba = 0.d0
								else
									target_ba = managementInputs(t_n(i),6,i)
								end if


								if (target_ba > 0.d0 .and. basal_area(i) > target_ba) then

                                    ! Compute stem biomass removal to reach target basal area
                                    removal = (200.d0 * sqrt(target_ba / (stems_n(i) * pi))) ** nWs(i) * &
                                        aWs(i) * stems_n(i) / 1000

                                    mort_manag(i) = (stems_n(i) - ceiling(removal / (biom_tree(i) / 1000.d0))) / &
                                            stems_n(i) / managementInputs(t_n(i),3,i)

                                    if (mort_manag(i) < 0.d0) then
                                        mort_manag(i) = 0.d0
                                    else if (mort_manag(i) > 1.d0) then
                                        mort_manag(i) = 1.d0
                                    end if

                                    final_harvest(i) = .FALSE.

								else if ( target_ba == 0.d0 ) then

                                    final_harvest(i) = .TRUE.

								end if

                                call remove_trees(mort_manag(i), stem_harv(i), stems_n(i), biom_foliage_debt(i), &
                                    biom_foliage(i), biom_stem_harv(i), biom_stem(i), biom_root(i), harvesting(i), &
                                    residues(i), root_residues(i), foliage_residues(i), managementInputs(t_n(i),5,i), &
                                    managementInputs(t_n(i),4,i), managementInputs(t_n(i),3,i), wood_density(ii,i), &
                                    fracBB(ii,i), f_dormant(month, leafgrow(i), leaffall(i)), b_cor, final_harvest(i))
                                
							! Check if management is performed based on removed volume (1st priority)
							else if ( managementInputs(t_n(i),7,i) >= 0.d0 ) then

								! Check if target volume is in percentage or absolute value
								if ( managementInputs(t_n(i),7,i) < 1.d0) then
									target_vol = volume(i) * managementInputs(t_n(i),7,i) 
								else if ( managementInputs(t_n(i),7,i) == 1.d0) then
									target_vol = 0.d0
								else
									target_vol = managementInputs(t_n(i),7,i)
								end if
								
								if (target_vol > 0.d0 .and. volume(i) > target_vol) then

                                    ! Compute stem biomass removal to reach target volume
                                    removal = biom_stem(i) - target_vol * wood_density(ii,i) / (1.d0 - fracBB(ii,i)) 
                                    
                                    mort_manag(i) = (stems_n(i) - ceiling(removal / (biom_tree(i) / 1000.d0))) / &
                                            stems_n(i) / managementInputs(t_n(i),3,i)

                                    if (mort_manag(i) < 0.d0) then
                                        mort_manag(i) = 0.d0
                                    else if (mort_manag(i) > 1.d0) then
                                        mort_manag(i) = 1.d0
                                    end if

                                    final_harvest(i) = .FALSE.

								else if ( target_vol == 0.d0 ) then
				
                                    final_harvest(i) = .TRUE.

								end if

                                call remove_trees(mort_manag(i), stem_harv(i), stems_n(i), biom_foliage_debt(i), &
                                    biom_foliage(i), biom_stem_harv(i), biom_stem(i), biom_root(i), harvesting(i), &
                                    residues(i), root_residues(i), foliage_residues(i), managementInputs(t_n(i),5,i), &
                                    managementInputs(t_n(i),4,i), managementInputs(t_n(i),3,i), wood_density(ii,i), &
                                    fracBB(ii,i), f_dormant(month, leafgrow(i), leaffall(i)), b_cor, final_harvest(i))
								
							! Otherwise use stem count
							else

								if( stems_n(i) > managementInputs(t_n(i),2,i) ) then

									mort_manag(i) = (stems_n(i) - managementInputs(t_n(i),2,i) ) / stems_n(i)

									!if the stand is thinned from above, then the ratios (F, R and S) of stem,
									! foliage and roots to be removed relative to the mean tree in the stand
									! will be >1. If the product of this ratio and delN is > 1 then the new
									! WF, WR or WS will be < 0, which is impossible. Therefore, make sure this is >= 0.

									if( maxval( mort_manag(i) * managementInputs( t_n(i),3:5,i)) >= 1.d0 ) then

                                        final_harvest(i) = .TRUE.

									else

                                        final_harvest(i) = .FALSE.

									end if

                                    call remove_trees(mort_manag(i), stem_harv(i), stems_n(i), biom_foliage_debt(i), &
                                        biom_foliage(i), biom_stem_harv(i), biom_stem(i), biom_root(i), harvesting(i), &
                                        residues(i), root_residues(i), foliage_residues(i), managementInputs(t_n(i),5,i), &
                                        managementInputs(t_n(i),4,i), managementInputs(t_n(i),3,i), wood_density(ii,i), &
                                        fracBB(ii,i), f_dormant(month, leafgrow(i), leaffall(i)), b_cor, final_harvest(i))

								end if

							end if

						   ! Add tree density related stand replacement

							if ( stems_n(i) <= 15.d0 .and. stems_n(i) > 0.d0 .or. dbh(i) > 80.d0) then

                                call remove_trees(mort_manag(i), stem_harv(i), stems_n(i), biom_foliage_debt(i), &
                                    biom_foliage(i), biom_stem_harv(i), biom_stem(i), biom_root(i), harvesting(i), &
                                    residues(i), root_residues(i), foliage_residues(i), managementInputs(t_n(i),5,i), &
                                    managementInputs(t_n(i),4,i), managementInputs(t_n(i),3,i), wood_density(ii,i), &
                                    fracBB(ii,i), f_dormant(month, leafgrow(i), leaffall(i)), b_cor, .TRUE.)
								
								!Update management schedule if clearcut is to be performed at a later period
								if (minval(managementInputs(:,2,i))==0.d0) then
									do while (managementInputs(t_n(i),2,i) > 0.d0 .and. t_n(i) < 200)
										t_n(i) = t_n(i) + 1
									end do
                                else
                                    do while (managementInputs(t_n(i),1,i) > 30.d0 .and. t_n(i) < 200)
										t_n(i) = t_n(i) + 1
									end do
								end if
							end if 

							! If final harvest was conducted replant/regenerate
							if (stems_n(i) < 1) then

								! Recalculate age-dependent parameters
                                call replant_stand(age, age_m, SLA, fracBB, wood_density, gammaN, gammaF, &
                                    f_age, stems_n, biom_foliage_debt, biom_foliage, biom_stem, &
                                    biom_root, biom_tree, sapling_wf, sapling_ws, sapling_wr, &
                                    n_sapling, leafgrow, leaffall,  MaxAge, rAge, nAge, &
                                    SLA0, SLA1, tSLA, fracBB0, fracBB1, tBB,  rho0, rho1, tRho, &
                                    gammaN0, gammaN1, tgammaN, ngammaN, gammaF1, gammaF0, tgammaF, &
                                    crop_trees, dbh_crop, height_crop, basal_area_crop, biom_stem_crop, &
                                    vol_crop, volume_cum, basal_area_prop, n_m, n_sp, ii, i, month) 

							end if

							t_n(i) = t_n(i) + 1
                            if (wind_dist .eq. int(1)) current_density(i) = stems_n(i)
                            
                        end if

                    end if

                end if

            end do

            ! Correct the bias
            if ( b_cor .eqv. .TRUE. ) then

                biom_tree(:) = biom_stem(:) * 1000.d0 / stems_n(:)  ! kg/tree
                where( stems_n(:) .eq. 0.d0 ) biom_tree(:) = 0.d0
                lai(:) =  biom_foliage(:) * SLA(ii,:) * 0.1d0

                do n = 1, b_n
                    competition_total(:) = sum( wood_density(ii,:) * basal_area(:) )

                    call s_sizeDist_correct(n_sp, age(ii,:), stems_n(:), biom_tree(:), competition_total(:), lai(:), &
                        correct_bias, height_model,  pars_i(69:85,:), pars_b, aWs(:), nWs(:), pfsPower(:), pfsConst(:), &
                        dbh(:), basal_area(:), height(:), crown_length(:), crown_width(:), pFS(:), bias_scale(:,:) )
                end do

                ! Adjust the old wolume after thinning
                volume(:) = biom_stem(:) * (1.d0 - fracBB(ii,:)) / wood_density(ii,:)
                where( aV(:) > 0 ) volume(:) = aV(:) * dbh(:) ** nVB(:) * height(:) ** nVH(:) * &
                    (dbh(:) * dbh(:) * height(:)) ** nVBH(:) * stems_n(:)

                volume_old(:) = volume(:)

                b_cor = .FALSE.
            end if

            ! Estimate crop tree parameters
            if (sum(crop_trees(:)) > 0.d0 .and. month == 12 ) then
                crop_trees_aux(:) = crop_trees(:)

                where (dbh(:) < 5.d0)
                    crop_trees_aux(:) = 0
                end where

                call    get_crop_trees(dbh_crop(:), vol_crop(:), basal_area_crop(:), biom_stem_crop(:), height_crop(:), &
                                     dbh(:), age(ii,:), height(:), stems_n(:), basal_area(:), volume(:), &
                                     volume_cum(:)*wood_density(ii,:), &
                                     competition_total(:), aWs(:), nWs(:), aH(:), nHB(:), nHC(:), &
                                     aV(:), nVB(:), nVH(:), nVBH(:), Dscale0(:), DscaleB(:), Dscalerh(:), Dscalet(:), &
                                     DscaleC(:), Dshape0(:), DshapeB(:), Dshaperh(:), Dshapet(:), DshapeC(:), Dlocation0(:), &
                                     DlocationB(:), Dlocationrh(:), Dlocationt(:), DlocationC(:), crop_trees_aux(:), fracBB(ii,:), &
                                     wood_density(ii,:), height_model, 5.d0, n_sp) 
            end if

            
            ! get dimensions of harvested trees
            call get_harvested_trees(dbh_harv(:), height_harv(:), aWs(:), nWs(:), aH(:), nHB(:), nHC(:), &
                                    competition_total(:), biom_stem_harv(:), stem_harv(:), n_sp, height_model)

            ! Mortality --------------------------------------------------------------------------

            ! Stress related ------------------
            mort_stress(:) = 0.d0
            dead_volume(:) = 0.d0
            do i = 1, n_sp
                if( f_dormant(month, leafgrow(i), leaffall(i)) .eqv. .FALSE.) then

                    if ( mortality_model == int(1) ) then

                        mort_stress(i) = gammaN(ii,i) / 12.d0 /100.d0
                        !mort_stress(i) = gammaN(ii,i) * stems_n(i) / 12.d0 /100.d0
                        !mort_stress(i) = ceiling( mort_stress(i) )
                        ! mort_stress(i) = min( mort_stress(i), stems_n(i)) ! Mortality can't be more than available

                        ! biom_foliage(i) = biom_foliage(i) - mF(i) * mort_stress(i) * (biom_foliage(i) / stems_n(i))
                        ! biom_root(i) = biom_root(i) - mR(i) * mort_stress(i) * (biom_root(i) / stems_n(i))
                        ! biom_stem(i) = biom_stem(i) - mS(i) * mort_stress(i) * (biom_stem(i) / stems_n(i))
                        ! stems_n(i) = stems_n(i) - mort_stress(i)

                        call remove_mortality(mort_stress(i), stems_n(i), biom_foliage(i), biom_root(i), biom_stem(i), &
                                                dead_volume(i), fracBB(ii,i), wood_density(ii,i), mF(i), mR(i), mS(i))

                        b_cor = .TRUE.
					!end if

				    ! Empirical survival probability curves by Brandl et al. (2020)
                    else if (  mortality_model == int(2) .and. month == 7) then

                        mort_stress(i) = mortality_brandl(gammaN(ii,i), age(ii,i), tmp_max(ii-6:ii+5), tmp_min(ii-6:ii+5), &
                                                           tmp_ave(ii-6:ii+5), prcp(ii-6:ii+5), month)

                        call remove_mortality(mort_stress(i), stems_n(i), biom_foliage(i), biom_root(i), biom_stem(i), &
                                                dead_volume(i), fracBB(ii,i), wood_density(ii,i), mF(i), mR(i), mS(i))

						b_cor = .TRUE.

                    ! LPJ-GUESS implementation 
                    else if ( mortality_model == int(3) .and. month == 7) then

                        if (ii <= 12) then
                            nyear = 1

                            mort_stress(i) = mortality_lpj(gammaN1(i), output((ii-nyear*6 + 1):ii,i,6,2), &
                            output((ii-nyear*6 + 1):ii,i,3,3), dmC, nyear, 2.d0)
                        else if (ii > 12 .and. ii <= 24) then
                            nyear = 1

                            mort_stress(i) = mortality_lpj(gammaN1(i), output((ii-nyear*12 + 1):ii,i,6,2), &
                            output((ii-nyear*12 + 1):ii,i,3,3), dmC, nyear, 1.d0)
                        else if (ii > 24 .and. ii <= 36) then
                            nyear = 2

                            mort_stress(i) = mortality_lpj(gammaN1(i), output((ii-nyear*12 + 1):ii,i,6,2), &
                            output((ii-nyear*12 + 1):ii,i,3,3), dmC, nyear, 1.d0)
                        else 
                            nyear = 3

                            mort_stress(i) = mortality_lpj(gammaN1(i), output((ii-nyear*12 + 1):ii,i,6,2), &
                            output((ii-nyear*12 + 1):ii,i,3,3), dmC, nyear, 1.d0)
                        endif

                        call remove_mortality(mort_stress(i), stems_n(i), biom_foliage(i), biom_root(i), biom_stem(i), &
                                                dead_volume(i), fracBB(ii,i), wood_density(ii,i), mF(i), mR(i), mS(i))

                        b_cor = .TRUE.

                    else if ( mortality_model == int(4) .and. month == 7) then

                        if (ii <= 12) then
                            mort_stress(i) = mortality_iland(output((ii-6):(ii+5),i,6,2), output((ii-6):(ii+5),i,4,7), &
                                                             output((ii-6):(ii+5),i,4,8), gammaN(ii,i))
                        else  
                            mort_stress(i) = mortality_iland(output((ii-11):ii,i,6,2), output((ii-11):ii,i,4,7), &
                                                             output((ii-11):ii,i,4,8), gammaN(ii,i))
                        end if

                        call remove_mortality(mort_stress(i), stems_n(i), biom_foliage(i), biom_root(i), biom_stem(i), &
                                                dead_volume(i), fracBB(ii,i), wood_density(ii,i), mF(i), mR(i), mS(i))

                        b_cor = .TRUE.
					end if

				else
					mort_stress(i) = 0.d0
				end if

            end do


            ! Correct the bias
            if ( b_cor .eqv. .TRUE. ) then

                biom_tree(:) = biom_stem(:) * 1000.d0 / stems_n(:)  ! kg/tree
                where( stems_n(:) .eq. 0.d0 ) biom_tree(:) = 0.d0
                lai(:) =  biom_foliage(:) * SLA(ii,:) * 0.1d0

                do n = 1, b_n
                    competition_total(:) = sum( wood_density(ii,:) * basal_area(:) )

                    call s_sizeDist_correct(n_sp, age(ii,:), stems_n(:), biom_tree(:), competition_total(:), lai(:), &
                        correct_bias, height_model,  pars_i(69:85,:), pars_b, aWs(:), nWs(:), pfsPower(:), pfsConst(:), &
                        dbh(:), basal_area(:), height(:), crown_length(:), crown_width(:), pFS(:), bias_scale(:,:) )
                end do

                b_cor = .FALSE.
            end if

            ! Self-thinning related ------------------
            mort_thinn(:) = 0.d0
            basal_area_prop(:) = basal_area(:) / sum( basal_area(:) )
            !basal_area_prop(:) if basal_area_prop(:) > 0 and basal_area_prop(:) < 0.01 put 0.01

			where( basal_area_prop(:) <0.01d0 )
			    basal_area_prop(:) = 0.01d0
			end where
            stems_n_ha(:) = stems_n(:) / basal_area_prop(:)

            biom_tree_max(:) = wSx1000(:) * (1000.d0 / stems_n_ha(:)) ** thinPower(:)

            do i = 1, n_sp

                if( f_dormant(month, leafgrow(i), leaffall(i)) .eqv. .FALSE.) then

                    if ( biom_tree_max(i) < biom_tree(i) ) then

                        mort_thinn(i) = f_get_mortality( stems_n_ha(i), biom_stem(i) / basal_area_prop(i) , &
                        mS(i), wSx1000(i), thinPower(i) ) * basal_area_prop(i)
						
                        !if( stems_n(i) < 1.d0 ) mort_thinn(i) = stems_n(i)
                        !mort_thinn(i) = ceiling( mort_thinn(i) )

                        if( mort_thinn(i) < stems_n(i) ) then

                            biom_foliage(i) = biom_foliage(i) - mF(i) * mort_thinn(i) * (biom_foliage(i) / stems_n(i))
                            biom_root(i) = biom_root(i) - mR(i) * mort_thinn(i) * (biom_root(i) / stems_n(i))
                            biom_stem(i) = biom_stem(i) - mS(i) * mort_thinn(i) * (biom_stem(i) / stems_n(i))
                            stems_n(i) = stems_n(i) - mort_thinn(i)
                        else

                            biom_foliage(i) = 0.d0
                            biom_root(i) = 0.d0
                            biom_stem(i) = 0.d0
                            stems_n(i) = 0.d0
                        end if

                        b_cor = .TRUE.

                    end if

                else
                    mort_thinn(i) = 0.d0

                end if
            end do


            ! Correct the bias
            if ( b_cor .eqv. .TRUE. ) then

                biom_tree(:) = biom_stem(:) * 1000.d0 / stems_n(:)  ! kg/tree
                where( stems_n(:) .eq. 0.d0 ) biom_tree(:) = 0.d0
                lai(:) =  biom_foliage(:) * SLA(ii,:) * 0.1d0

                do n = 1, b_n
                    competition_total(:) = sum( wood_density(ii,:) * basal_area(:) )

                    call s_sizeDist_correct(n_sp, age(ii,:), stems_n(:), biom_tree(:), competition_total(:), lai(:), &
                        correct_bias, height_model,  pars_i(69:85,:), pars_b, aWs(:), nWs(:), pfsPower(:), pfsConst(:), &
                        dbh(:), basal_area(:), height(:), crown_length(:), crown_width(:), pFS(:), bias_scale(:,:) )
                end do

                b_cor = .FALSE.
            end if
			
			
			! Wind disturbance ------------------
			
			if (wind_dist .eq. int(1) .and. ii >= dist_start) then

				wind_mortality(:) = 0.d0
                biom_stem_wind(:) = 0.d0
                biom_root_wind(:) = 0.d0
                biom_foliage_wind(:) = 0.d0

				if ( month == 12) then
				
					!get_gust(crown_width(:),crown_length(:),height(:),crown_area_fact(:),wind_speed,n_sp)
                    !get critical wind speed
                    do i = 1, n_sp
                        if (crown_width(i)==0.d0) then
                            crown_width_aux(i) = maxval(output((ii-12):ii,i,2,9))
                        else
                            crown_width_aux(i) = crown_width(i)
                        end if
                    end do

                    thinn_adjust(:) = get_thinning_effect(current_density, previous_density, lag, n_sp) 

					cws(:) = get_cws_tmc(height(:), dbh(:), crown_length(:), crown_width_aux(:), MOE(:) * 10.d0**6.d0,&
                             MOR(:) * 10.d0**6.d0, biom_tree(:), wet_fact(:), competition_total(:), stems_n(:), Creg(:), &
                             snow_water, .FALSE. , thinn_adjust(:), n_sp)

                    ! Elevate CWS to anemometer height
                    cws10(:) = get_crown_wind_speed(crown_width_aux(:),crown_length(:),height(:),stems_n(:),cws(:), &
                                                    Cdrag(:),Ndrag(:), sum(lai(:)), n_sp)                             
                    ! If deciduous in dormant period assume 30% higher cws and no drag
                    where(crown_width(:) <= 0.1d0)
                        cws10(:) = cws(:) * 1.3d0
                    end where

					! Define return period (10 years)
                    !wind_damage(:) = get_damage_probability(weibull_a(ii),weibull_k(ii), cws(:), n_sp)
                    wind_damage(:) = 0.d0        

                    if (weibull_k(ii) <= 0.d0 .and. ii == dist_start) then
                        wind_damage(:) = 1.d0
                    else if (weibull_k(ii) > 0.d0) then
                        ! wind_damage(:) = min(get_damage_probability(weibull_a(ii),weibull_k(ii), min(cws10(:),30.d0), n_sp), &
                        !                     1.d0/15.d0)
                        wind_damage(:) = 1.d0/10.d0
                    end if
                        
					where ( wind_damage(:) > 1.d0)
						wind_damage(:) = 1.d0
					end where
					
					where ( wind_damage(:) < 0.d0)
						wind_damage(:) = 0.d0
					end where
					
					call random_number(rnd_nr)					
				
					if ( rnd_nr < MAXVAL(wind_damage(:))) then

                        call random_number(rnd_nr)

						do i = 1, n_sp
							if (age(ii,i) > 30.d0 .and. height(i) > 15.d0) then
							
								! Get damage - upscales wind damage probability to share of damaged trees.
								
                                wind_mortality(i) = get_damage_wind(weibull_a(ii),weibull_k(ii), cws10(i), stems_n(i), rnd_nr) 
								
								wind_mortality(i) = min( wind_mortality(i), stems_n(i)) !* 0.01d0*(-log(1.d0-rnd_nr))/0.15d0 ! Mortality can't be more than stock

                                biom_stem_wind(i) = wind_mortality(i) * (1.2d0 * biom_stem(i) / stems_n(i))
                                biom_root_wind(i) = wind_mortality(i) * (1.2d0 * biom_root(i) / stems_n(i))
                                biom_foliage_wind(i) = wind_mortality(i) * (1.2d0 * biom_foliage(i) / stems_n(i))

								biom_foliage(i) = biom_foliage(i) - biom_foliage_wind(i) 
								biom_root(i) = biom_root(i) - biom_root_wind(i) 
								biom_stem(i) = biom_stem(i) - biom_stem_wind(i)
								stems_n(i) = stems_n(i) - wind_mortality(i)
								dead_volume(i) = dead_volume(i) + biom_stem_wind(i) * (1.d0 - fracBB(ii,i)) / wood_density(ii,i)
									
								! Add wind related stand replacement
								if ( stems_n(i) <= 15.d0 .or. biom_stem(i) <= 0.d0) then

                                    call remove_trees(mort_manag(i), stem_harv(i), stems_n(i), biom_foliage_debt(i), &
                                        biom_foliage(i), biom_stem_harv(i), biom_stem(i), biom_root(i), harvesting(i), &
                                        residues(i), root_residues(i), foliage_residues(i), managementInputs(t_n(i),5,i), &
                                        managementInputs(t_n(i),4,i), managementInputs(t_n(i),3,i), wood_density(ii,i), &
                                        fracBB(ii,i), f_dormant(month, leafgrow(i), leaffall(i)), b_cor, .TRUE.)
									
									! Update management schedule if clearcut is to be performed at a later period
									if (minval(managementInputs(:,2,i))==0.d0) then
										do while (managementInputs(t_n(i),2,i) > 0.d0 .and. t_n(i) < int(200))
											t_n(i) = t_n(i) + 1
										end do
									end if

                                    call replant_stand(age, age_m, SLA, fracBB, wood_density, gammaN, gammaF, &
                                        f_age, stems_n, biom_foliage_debt, biom_foliage, biom_stem, &
                                        biom_root, biom_tree, sapling_wf, sapling_ws, sapling_wr, &
                                        n_sapling, leafgrow, leaffall,  MaxAge, rAge, nAge, &
                                        SLA0, SLA1, tSLA, fracBB0, fracBB1, tBB,  rho0, rho1, tRho, &
                                        gammaN0, gammaN1, tgammaN, ngammaN, gammaF1, gammaF0, tgammaF, &
                                        crop_trees, dbh_crop, height_crop, basal_area_crop, biom_stem_crop, &
                                        vol_crop, volume_cum, basal_area_prop, n_m, n_sp, ii, i, month) 
							

								end if 	
																							
								
							end if
						end do
						
					end if 
				
				end if
			end if


            ! Bark beetle disturbance
            if (beetle_dist .eq. int(1) .and. ii > dist_start .and. spruce_idx > 0) then
                !spruce_idx = 1
                bb_mort = 0.d0
                wind_lit = 0.d0
                thinn_lit = 0.d0
                biom_stem_bb(:) = 0.d0
                biom_root_bb(:) = 0.d0
                biom_foliage_bb(:) = 0.d0
                greff_spruce = 80.d0 ! set as intermediate, in LPJ = 40.d0

                if (ii <= 12) then
                    wtmin = 0.d0
                    wtmin = (sum(tmp_min(1:2))/2.d0)
                    if (wtmin < -10.d0) then
                        background_infestation = 0.d0
                    end if
                    
                else
                    wtmin = 0.d0
                    wtmin = (sum(tmp_min((ii-12):(ii-10)))/3.d0)
                    if (wtmin < -10.d0) then
                        background_infestation = 0.d0
                    end if
                end if

                if (dbh(spruce_idx) > 10.d0 .and. age(ii, spruce_idx) > 30.d0 .and. month == 12) then
                    
                    if (ii <= 12) then
                        nyear = 1
                        previous_damage = 0.d0
                        maxbiom_sp =  output(ii- 1 ,spruce_idx,4,1)
                        wind_lit = 0.d0
                        greff_par = calc_greff(output((ii-nyear*12 + 1):ii,spruce_idx,6,2), &
                                    output((ii-nyear*12 + 1):ii,spruce_idx,3,3), dmC, nyear, 1.d0)

                    else if (ii > 12 .and. ii <= 24) then
                        nyear = 2
                        previous_damage = output(ii-12,spruce_idx,12,6)  * (1.d0-salvage_biotic(spruce_idx))
                        maxbiom_sp = output(ii-12 ,spruce_idx,4,1) + previous_damage
                        wind_lit = output(ii-12,spruce_idx,12,4) * (1.d0-salvage_wind(spruce_idx))
                        
                        lit_dbh = ( mS(spruce_idx) *  output(ii-12,spruce_idx,4,4) / aWs(spruce_idx)) ** (1.d0 / nWs(spruce_idx))

                        greff_par = calc_greff(output((ii-nyear*12 + 1):ii,spruce_idx,6,2), &
                                    output((ii-nyear*12 + 1):ii,spruce_idx,3,3), dmC, nyear, 1.d0) 

                        ! check if natural mortality is added to beetle substrate
                        if ( lit_dbh > 10.d0 ) then
                            stress_lit = mS(spruce_idx)*output(ii-12,spruce_idx,8,3) / output(ii- 13,spruce_idx,2,2)
                            thinn_lit = mS(spruce_idx)*output(ii-12,spruce_idx,8,4) / output(ii- 13,spruce_idx,2,2)
                            wind_lit  =  wind_lit + stress_lit + thinn_lit
                        end if
                    else
                        nyear = 3
                        previous_damage = output(ii-12,spruce_idx,12,6) * (1.d0-salvage_biotic(spruce_idx))
                        maxbiom_sp = output(ii-12 ,spruce_idx,4,1) + previous_damage

                        wind_lit = (output(ii-12,spruce_idx,12,4) + output(ii-24,spruce_idx,12,4)) * (1.d0-salvage_wind(spruce_idx))

                        lit_dbh = ( mS(spruce_idx) *  output(ii-12,spruce_idx,4,4) / aWs(spruce_idx)) ** (1.d0 / nWs(spruce_idx))

                        greff_par = calc_greff(output((ii-nyear*12 + 1):ii,spruce_idx,6,2), &
                                    output((ii-nyear*12 + 1):ii,spruce_idx,3,3), dmC, nyear, 1.d0)

                        ! check if natural mortality is added to beetle substrate
                        if ( lit_dbh > 10.d0 ) then
                            stress_lit = mS(spruce_idx)*output(ii-24,spruce_idx,8,4) / output(ii- 25,spruce_idx,2,2)
                            thinn_lit = mS(spruce_idx)*output(ii-24,spruce_idx,8,3) / output(ii- 25,spruce_idx,2,2)
                            wind_lit  =  wind_lit + stress_lit + thinn_lit
                        end if

                    end if

                    call 	bb_landscape_risk(	infestation_share, &
                                                bb_mort, &
                                                age(ii, spruce_idx), & 
                                                output((ii-11):ii,spruce_idx,5,6) , int(1), & 
                                                wind_lit, output(ii-1,spruce_idx,4,1), & 
                                                basal_area_prop(spruce_idx), & 
                                                w_bio, w_age, w_host, w_vita, w_dens, w_wind, &
                                                greff_par(spruce_idx), greff_spruce, d_sus_type, &
                                                previous_damage, maxbiom_sp, bb_cst, & 
                                                lat, tmp_ave((ii-11):ii), wtmin, solar_rad((ii-11):ii), &
                                                output((ii-11):ii,spruce_idx,3,3), &
                                                k(spruce_idx), day_length(:), bb_cbp, & 
                                                background_infestation )

                    ! Get damage
                    biom_stem_bb(spruce_idx) = bb_mort * biom_stem(spruce_idx)
                    biom_root_bb(spruce_idx) = bb_mort * biom_root(spruce_idx)
                    biom_foliage_bb(spruce_idx) = bb_mort * biom_foliage(spruce_idx)


                    biom_foliage(spruce_idx) = biom_foliage(spruce_idx) - biom_foliage_bb(spruce_idx)
                    biom_root(spruce_idx) = biom_root(spruce_idx) - biom_root_bb(spruce_idx)
                    biom_stem(spruce_idx) = biom_stem(spruce_idx) - biom_stem_bb(spruce_idx) 
                    stems_n(spruce_idx) = stems_n(spruce_idx) - bb_mort * stems_n(spruce_idx) / 1.2d0
                    dead_volume(spruce_idx) = dead_volume(spruce_idx) + biom_stem_bb(spruce_idx) * &
                                            (1.d0 - fracBB(ii,spruce_idx)) / wood_density(ii,spruce_idx)
                  
                    if ( stems_n(spruce_idx) <= 15.d0 .or. biom_stem(spruce_idx) <= 0.d0) then

                        call remove_trees(mort_manag(spruce_idx), stem_harv(spruce_idx), stems_n(spruce_idx), &
                            biom_foliage_debt(spruce_idx), biom_foliage(spruce_idx), biom_stem_harv(spruce_idx), &
                            biom_stem(spruce_idx), biom_root(spruce_idx), harvesting(spruce_idx), residues(spruce_idx), &
                            root_residues(spruce_idx), foliage_residues(spruce_idx), &
                            managementInputs(t_n(spruce_idx),5,spruce_idx), managementInputs(t_n(spruce_idx),4,spruce_idx), & 
                            managementInputs(t_n(spruce_idx),3,spruce_idx), wood_density(ii,spruce_idx), &
                            fracBB(ii,spruce_idx), .FALSE., b_cor, .TRUE.)
                        
                        ! Update management schedule if clearcut is to be performed at a later period
                        if (minval(managementInputs(:,2,spruce_idx))==0.d0) then
                            do while (managementInputs(t_n(spruce_idx),2,spruce_idx) > 0.d0 .and. t_n(spruce_idx) < int(200))
                                t_n(spruce_idx) = t_n(spruce_idx) + 1
                            end do
                        end if

                        call replant_stand(age, age_m, SLA, fracBB, wood_density, gammaN, gammaF, &
                            f_age, stems_n, biom_foliage_debt, biom_foliage, biom_stem, &
                            biom_root, biom_tree, sapling_wf, sapling_ws, sapling_wr, &
                            n_sapling, leafgrow, leaffall,  MaxAge, rAge, nAge, &
                            SLA0, SLA1, tSLA, fracBB0, fracBB1, tBB,  rho0, rho1, tRho, &
                            gammaN0, gammaN1, tgammaN, ngammaN, gammaF1, gammaF0, tgammaF, &
                            crop_trees, dbh_crop, height_crop, basal_area_crop, biom_stem_crop, &
                            vol_crop, volume_cum, basal_area_prop, n_m, n_sp, ii, spruce_idx, month) 
                        

                    end if 	 	 	

                end if
            else
                biom_stem_bb(:) = 0.d0
                biom_root_bb(:) = 0.d0
                biom_foliage_bb(:) = 0.d0
            end if
			
			! Inputs for soil models
			fr_cr_ratio(:) = biom_foliage(:) * fr_ratio(:) * fr_wn(:) /  biom_root(:)
            stem_bark_corr(:) = min( 0.06d0, fracBB(ii,:))

            dead_volume(:) = (mS(:) * (mort_stress(:) + mort_thinn(:)) * (biom_stem(:) / stems_n(:))) * &
                cwd_removal(:) / wood_density(ii,:)

            residues(:) = residues(:) + dead_volume(:)
			
			if (soil_model .ne. int(2)) then

                !Recalcitrant pool
				Yr_input(:) = (mS(:) * (mort_stress(:) + mort_thinn(:)) * (biom_stem(:) / stems_n(:))  * dmC) * &
						(1.d0-cwd_removal(:)) + (mR(:) * (mort_stress(:) + mort_thinn(:)) * &
						(biom_root(:) / stems_n(:)) * (1.d0-fr_cr_ratio(:))  * dmC) + &                         
						(residues(:) * (1.d0-cwd_removal(:)) + root_residues(:) * (1.d0-fr_cr_ratio(:))) * &
                        wood_density(ii,:) * dmC 

                !Labile pool
				Yl_input(:) = (biom_loss_foliage(:) + biom_loss_root(:)) * dmC + &
						(mF(:) * (mort_stress(:) + mort_thinn(:)) * (biom_foliage(:) / stems_n(:)) * dmC) + &
						(mR(:) * (mort_stress(:) + mort_thinn(:)) * (biom_root(:) / stems_n(:)) * fr_cr_ratio(:)  * dmC) + & 
						(foliage_residues(:) * (1.d0-cwd_removal(:)) + &
                        root_residues(:) * wood_density(ii,:) * fr_cr_ratio(:)) * dmC 
                   
                ! Add wind disturbance damage to litter
                if (wind_dist .eq. int(1)) then

                    ! Recalcitrant pool
                    Yr_input(:) = Yr_input(:) + biom_stem_wind(:) * (1.d0 - salvage_wind(:)) *  dmC + &
                                biom_root_wind(:) * (1.d0-fr_cr_ratio(:)) * dmC

                    !Labile pool
                    Yl_input(:) = Yl_input(:) + biom_foliage_wind(:) * dmC + biom_root_wind(:) * (fr_cr_ratio(:)) * dmC

                end if  
                ! Add wind disturbance damage to litter
                if (beetle_dist .eq. int(1)) then

                    ! Recalcitrant pool
                    Yr_input(:) = Yr_input(:) + biom_stem_bb(:)*(1.d0 - salvage_biotic(:)) *  dmC + &
                                    biom_root_bb(:) * (1.d0-fr_cr_ratio(:)) * dmC

                    !Labile pool
                    Yl_input(:) = Yl_input(:) + biom_foliage_bb(:) * dmC + biom_root_bb(:) * (fr_cr_ratio(:)) * dmC

                end if 
			else
                ! CDW
				Yr_input(:) = (mS(:) * (mort_stress(:) + mort_thinn(:) ) * (biom_stem(:) / stems_n(:))  * dmC) * &
						(1.d0-cwd_removal(:)) * (1.d0 - fracBB(ii,:) + stem_bark_corr(:)) 		

				! Branches and roots		
				Yb_input(:) = (mS(:) * (mort_stress(:) + mort_thinn(:)) * (biom_stem(:) / stems_n(:))  * dmC) * &
						(1.d0-cwd_removal(:)) * (fracBB(ii,:) - stem_bark_corr(:)) + (mR(:) * (mort_stress(:) + mort_thinn(:) ) * &
						(biom_root(:) / stems_n(:)) * (1.d0-fr_cr_ratio(:))  * dmC) + & 
						(residues(:) * (1.d0-cwd_removal(:)) + root_residues(:) * (1.d0-fr_cr_ratio(:))) * wood_density(ii,:) * dmC 							
				! Labile		
				Yl_input(:) = (biom_loss_foliage(:) + biom_loss_root(:)) * dmC + &
						(mF(:) * (mort_stress(:) + mort_thinn(:)) * (1.d0-cwd_removal(:)) * (biom_foliage(:) / stems_n(:)) * dmC) + &
						(mR(:) * (mort_stress(:) + mort_thinn(:)) * (biom_root(:) / stems_n(:)) * fr_cr_ratio(:)  * dmC) + & 
						(foliage_residues(:) * (1.d0-cwd_removal(:)) + root_residues(:) * fr_cr_ratio(:) * wood_density(ii,:)) * dmC 
                
                ! Add wind disturbance damage to litter
                if (wind_dist .eq. int(1)) then
                    ! CWD
                    Yr_input(:) = Yr_input(:) + (biom_stem_wind(:) * (1.d0 - salvage_wind(:))) * &
                                    (1.d0 - fracBB(ii,:) + stem_bark_corr(:)) *  dmC 

                    ! Branches and roots
                    Yb_input(:) = Yb_input(:) + biom_root_wind(:) * (1.d0-fr_cr_ratio(:)) *dmC + & 
                                biom_stem_wind(:) * (1.d0 - salvage_wind(:)) * (fracBB(ii,:) - stem_bark_corr(:)) *  dmC

                    ! Litter pool
                    Yl_input(:) = Yl_input(:) + biom_foliage_wind(:) * dmC + biom_root_wind(:) * (fr_cr_ratio(:)) * dmC
                end if 
                ! Add bark beetle disturbance damage to litter
                if (beetle_dist .eq. int(1)) then

                    ! CWD
                    Yr_input(:) = Yr_input(:) + biom_stem_bb(:)*(1.d0 - salvage_biotic(:)) * &
                                    (1.d0 - fracBB(ii,:) + stem_bark_corr(:)) *  dmC 

                    ! Branches and roots
                    Yb_input(:) = Yb_input(:) + biom_root_bb(:) * (1.d0-fr_cr_ratio(:)) *dmC + & 
                                biom_stem_bb(:)*(1.d0 - salvage_biotic(:)) * (fracBB(ii,:) - stem_bark_corr(:)) *  dmC

                    ! Litter pool
                    Yl_input(:) = Yl_input(:) + biom_foliage_bb(:) * dmC + biom_root_bb(:) * (fr_cr_ratio(:)) * dmC
                end if 
			end if


            if (wind_dist .eq. int(1) .or. beetle_dist .eq. int(1)) then
                salvage(:) = (biom_stem_bb(:) * salvage_biotic(:) + biom_stem_wind(:) * salvage_wind(:)) * &
                             (1.d0 - fracBB(ii,:)) / wood_density(ii,:)
            end if        
										
							
            ! Soil C and N modules

            ! ICBM/2N soil module -------------------------------------------------------------------------------
			if ( soil_model .eq. int(1) ) then
			
			
			call	run_icbm(krmax, klmax, komax(1), Yl_C, Yr_C, O_C, Yl_N, Yr_N, O_N, Litter_C, &
					 hc, el, er, qbc, qh, qil, qir, Yl_input, Yr_input, O_Coutflux, &
					 Yl_Coutflux, Yr_Coutflux, Nav, TotalCarbo, TotalNitro, ((biom_loss_foliage(:) + foliage_residues(:)) * dmC), &
					 tmp_ave(ii), prcp(ii), excessSW, ASW, asw_max, soil_class, N_deposition(ii), dmC, n_sp)

			
                ! Get N availability and demand
                     
				Un(:) = 0.d0
                Up(:) = 0.d0

                ! Get nitrogen demand
			    do i=1, n_sp

					Un(i) =  biom_incr_foliage(i) * Ncf(i) / 100.d0 + biom_incr_root(i) * 0.67d0 * Ncr(i) / 100.d0 + &
						(biom_incr_stem(i) + biom_incr_root(i) * 0.33d0) * Ncs(i) / 100.d0

			    end do
                
                ! Estimate fertility
                if ( nut_lim .eq. int(1)) then	
                    root_area(:) = SRA(:) * biom_fineroot(:)
                    call compute_fertility(Nav, Un(:), Up(:), root_dist(:), root_area(:), 0.031d0, fertility(:), n_sp)
                end if

				
			end if

			! End of ICBM/2N soil module ---------------------------------------------------------------------


			! Yasso soil module, based on PREBAS implementation (https://github.com/ForModLabUHel/Rprebasso)-----

			if ( soil_model .eq. int(2) ) then

				litter_input = 0.d0

                hc(:) = parameters_yasso(28)
                
                ! If nutrient balance is necessary, a monthly time step is used, otherwise use yearly computations 
                if (nut_lim .eq. int(1)) then 

                    do i=1,n_sp

                        Yr_C(i) = Yr_input(i)
                        Yl_C(i) = Yl_input(i)
                        Yb_C(i) = Yb_input(i)

                        !Call Yasso for leaf and root litter
                        litter_input = Yl_C(i)
                        call compAWENH(litter_input,folAWENH,parameters_AWEN(1:5,i))
                        call mod5c(parameters_yasso,tstep,tmp_ave(ii),prcp(ii),soilC(ii-1,i,1,:),&
                                    folAWENH,size_eff,leaching,soilC(ii,i,1,:),0)
                        !call mod5c(parameters_yasso,tstep,tmp_yasso_year,prcp_year,soilC(ii-1,i,1,:),folAWENH,size_eff,leaching,soilC(ii,i,1,:),0)
                        
                        !Call Yasso for branches and coarse roots
                        litter_input = Yb_C(i)
                        call compAWENH(litter_input,crbAWENH,parameters_AWEN(6:10,i))
                        call mod5c(parameters_yasso,tstep,tmp_ave(ii),prcp(ii),soilC(ii-1,i,2,:),&
                                    crbAWENH,size_eff,leaching,soilC(ii,i,2,:),0)
                        !call mod5c(parameters_yasso,tstep,tmp_yasso_year,prcp_year,soilC(ii-1,i,1,:),folAWENH,size_eff,leaching,soilC(ii,i,1,:),0)

                        ! Call Yasso for CWD
                        litter_input = Yr_C(i)
                        call compAWENH(litter_input,stAWENH,parameters_AWEN(11:15,i))
                        call mod5c(parameters_yasso,tstep,tmp_ave(ii),prcp(ii),soilC(ii-1,i,3,:),&
                                    stAWENH,dbh(i),leaching,soilC(ii,i,3,:),0)
                        !call mod5c(parameters_yasso,tstep,tmp_yasso_year,prcp_year,soilC(ii-1,i,2,:),stAWENH,dbh(i),leaching,soilC(ii,i,2,:),0)

                    end do

                    ! Call Yasso for Humus layer
                    huAWENH(:) = 0.d0
                    call mod5c(parameters_yasso,tstep,tmp_ave(ii),prcp(ii),soilC(ii-1,1,4,:),&
                                huAWENH,size_eff,leaching,soilC(ii,1,4,:),0)
                    !call mod5c(parameters_yasso,tstep,tmp_yasso_year,prcp_year,soilC(ii-1,1,4,:),huAWENH,size_eff,leaching,soilC(ii,1,4,:),0)
                    call yasso_outflux (soilC, Yr_C, Yl_C, O_C, Yr_input, Yb_input, Yl_input, Yl_Coutflux, Yr_Coutflux, &
                                        O_Coutflux, n_sp, n_m, ii)

                    call yasso_nav (Yr_C, Yl_C, O_C, Yl_Coutflux, Yr_Coutflux, O_Coutflux, humification_N_l, humification_N_r, & 
                                    Yr_Noutflux, Yl_Noutflux, O_Noutflux, Yr_N, Yl_N, O_N, hc, qh, el, qbc, qir, qil, er, &
                                    Yr_input, Yl_input, Yb_input, soil_class, excessSW, ASW, N_deposition(ii), Nav, n_sp)

                    ! Get N availability and demand

                    Un(:) = 0.d0
                    Up(:) = 0.d0
        
                    ! Get nitrogen demand
                    do i=1, n_sp
                        Un(i) =  biom_incr_foliage(i) * Ncf(i) / 100.d0 + biom_incr_root(i) * 0.67d0 * Ncr(i) / 100.d0 + &
                            (biom_incr_stem(i) + biom_incr_root(i) * 0.33d0) * Ncs(i) / 100.d0
        
                    end do

                    ! Get total soil, litter and deadwood carbon
                    TotalCarbo = sum(soilC(ii,:,:,:))

                    !Estimate fertility
                    if ( nut_lim .eq. int(1)) then	
                        root_area(:) = SRA(:) * biom_fineroot(:)
                        call compute_fertility(Nav, Un(:), Up(:), root_dist(:), root_area(:), 0.031d0, fertility(:), n_sp)
                    end if

                else

                    if (month == 12) then

                        if (ii <= 12) then
                            step_aux = 11
                        else
                            step_aux = 12
                        end if

                        do i=1,n_sp

                            ! accounting variables
                            Yb_input_aux(i) = sum( output((ii-11):(ii-1),i,11,3) ) + Yb_input(i)
                            Yr_input_aux(i) = sum( output((ii-11):(ii-1),i,11,2) ) + Yr_input(i)
                            Yl_input_aux(i) = sum( output((ii-11):(ii-1),i,11,1) ) + Yl_input(i)

                            !Call Yasso for leaf and root litter
                            litter_input = Yl_input_aux(i)
                            call compAWENH(litter_input,folAWENH,parameters_AWEN(1:5,i))
                            call mod5c_year(parameters_yasso,tstep,tmp_ave((ii-11):ii),sum(prcp((ii-11):ii)), &
                                            soilC(ii-step_aux,i,1,:),folAWENH,size_eff,leaching,soilC(ii,i,1,:),0)
                            
                            !Call Yasso for branches and coarse roots
                            litter_input = Yb_input_aux(i)
                            call compAWENH(litter_input,crbAWENH,parameters_AWEN(6:10,i))
                            call mod5c_year(parameters_yasso,tstep,tmp_ave((ii-11):ii),sum(prcp((ii-11):ii)), &
                                            soilC(ii-step_aux,i,2,:), crbAWENH,size_eff,leaching,soilC(ii,i,2,:),0)

                            ! Call Yasso for CWD
                            litter_input = Yr_input_aux(i)
                            call compAWENH(litter_input,stAWENH,parameters_AWEN(11:15,i))
                            call mod5c_year(parameters_yasso,tstep,tmp_ave((ii-11):ii), sum(prcp((ii-11):ii)), &
                                            soilC(ii-step_aux,i,3,:), stAWENH,dbh(i),leaching,soilC(ii,i,3,:),0)

                        end do     
                        
                        huAWENH(:) = 0.d0
                        call mod5c_year(parameters_yasso,tstep,tmp_ave((ii-11):ii),sum(prcp((ii-11):ii)), &
                                        soilC(ii-step_aux,1,4,:), huAWENH,size_eff,leaching,soilC(ii,1,4,:),0)

                        call yasso_outflux_year(soilC, Yr_C, Yl_C, O_C, Yr_input_aux(:), Yb_input_aux(:), &
                                                Yl_input_aux(:), Yl_Coutflux, Yr_Coutflux, O_Coutflux, n_sp, &
                                                n_m, ii, step_aux)

                    end if

                end if
			end if

					
            ! Additional calculations ------------------
            basal_area_prop(:) = basal_area(:) / sum( basal_area(:) )

            ! Efficiency
            epsilon_gpp(:) = 100.d0 * gpp(:) / apar(:)
            epsilon_npp(:) = 100.d0 * npp_f(:) / apar(:)
            epsilon_biom_stem(:) = 100.d0 * biom_incr_stem(:) / apar(:)

            where( apar(:) .eq. 0.d0 )
                epsilon_gpp(:) = 0.d0
                epsilon_npp(:) = 0.d0
                epsilon_biom_stem(:) = 0.d0
            end where
            
            ! Save end of the month results
            include 'i_write_out.h'

        end do

    end subroutine s_3PG_f    

end module mod_3PG
