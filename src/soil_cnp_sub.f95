module soil_cnp_subroutines


IMPLICIT NONE
public :: compute_leaching, soil_temperature_swat, soil_temperature_toy, &
			calc_psi_soil, compute_fertility

contains

	!Subroutines for leaching------- ----------------------------------------------------
	subroutine compute_leaching( soil_class, leaching_org, leaching_min, excessSW, ASW )

		integer, intent(in) :: soil_class
		real (kind=8), intent(in) :: excessSW, ASW
		real (kind=8), intent(inout) :: leaching_org, leaching_min
		real (kind=8) :: sand_prop, leach_org, leach_min 

		if (soil_class == 1) then
			sand_prop = 0.85d0
		else if (soil_class == 2) then
			sand_prop = 0.725d0
		else
			sand_prop = 0.5d0
		end if

		leaching_org = 	0.000183d0 * excessSW * (0.01d0 + 0.04d0 * sand_prop)
		leaching_min = 	excessSW / ASW
		
	end subroutine compute_leaching


	! Subroutine to calculate nitrogen gas losses----------------------------------------------------
	subroutine nitrogen_gas(n_sp, n_inorg_T, n_gas, nminl_gas, minl_n_tot, immob_n_tot )
		! Input parameters
		integer, intent(in) :: n_sp
		! In/Out parameters
		real(kind=8), dimension(n_sp), intent(inout) :: n_inorg_T, n_gas

		real(kind=8), intent(in) :: nminl_gas
		real(kind=8), dimension(n_sp), intent(in) :: minl_n_tot, immob_n_tot



		! Local variables
		integer :: i
		real(kind=8), dimension(n_sp)::  n_minl_tot_pr

		do i=1, n_sp
			n_minl_tot_pr(i) = MAX(minl_n_tot(i) - immob_n_tot(i), 0.0)  ! calculate total mineralization of N
			n_gas(i) = nminl_gas * n_minl_tot_pr(i)  ! calculate fraction of net mineralization of N that is lost as gas
			n_inorg_T(i) = n_inorg_T(i) - n_gas(i)  ! subtract N lost as gas from inorganic N pool
		end do

	end subroutine nitrogen_gas


	!Subroutines for soil temperature ----------------------------------------------------
	subroutine soil_temperature_swat(t_soil, t_avg, t_year_avg, asw, snow, lai, m, sc)
	!Calculates the soil temperature based on the SWAT model implementation,
	!assumes temperature at 50 cm depth 

		!Input
		integer, intent(in) :: m, sc
		real(kind=8), intent(in) :: t_avg, t_year_avg, asw, snow, lai
		real(kind=8), intent(out) :: t_soil
		real(kind=8) :: t_soil_upper, t_soil_lower
	
		!Parameters
		real(kind=8) :: lambda, df, zd, z, dd, dd_max, z_tot, phi, t_bare, theta, bdod, t_soil_surf !, bcv, bcv_snow, bcv_veg
		!lag parameter
		lambda = 0.0d0
		
		!soil depth and profile depth (mm)
		z_tot = 1000.d0
		
		!get bare ground temperature
		call soil_temperature_toy(t_soil_surf, t_avg, m, lai, snow)
		

		! soil water correction factor
		if (sc == 1) then
			theta = (50.d0 + asw)
			bdod = 1.5d0
		else if (sc == 2) then
			theta = (80.d0 + asw)
			bdod = 1.45d0
		else if (sc == 3) then
			theta = (150.d0 + asw)
			bdod = 1.42d0
		else 
			theta = (200.d0 + asw)
			bdod = 1.4d0
		end if
		
		!maximum damping depth (mm)
		dd_max = 1000.d0 + 2500.d0 * bdod / ( bdod + 686.d0 * ( exp( -5.63d0 * bdod ) ) )

		phi = theta / ( ( 0.356d0 - 0.144d0 * bdod ) * z_tot )
		
		!damping depth
		dd = dd_max * exp(  log( 500.d0 / dd_max ) * ( ( 1.d0 - phi ) / ( 1.d0 + phi ) )**2.d0 )
		
		!layer depth
		z = 200.d0
		
		!ratio of depth to damping depth
		zd = z / dd
		
		!depth factor
		df = zd / ( zd + exp(-0.867d0  - 2.078d0 * zd ))
		
		!soil temperature update
		t_soil_upper = lambda * t_soil + ( 1.d0 - lambda ) * ( df * ( t_year_avg - t_soil_surf ) + t_soil_surf )
		
		!layer depth
		z = 700.d0
		
		!ratio of depth to damping depth
		zd = z / dd
		
		!depth factor
		df = zd / ( zd + exp(-0.867d0  - 2.078d0 * zd ))
		
		!soil temperature update
		t_soil_lower = lambda * t_soil + ( 1.d0 - lambda ) * ( df * ( t_year_avg - t_soil_surf ) + t_soil_surf )
		
		t_soil = 0.7d0 * t_soil_upper + 0.3d0 * t_soil_lower

	end subroutine soil_temperature_swat
	
	subroutine soil_temperature_toy(t_soil, t_avg, m, lai, snow)

		!Input
		integer, intent(in) :: m
		real(kind=8), intent(in) :: t_avg, lai, snow
		real(kind=8), intent(out) :: t_soil
		real(kind=8) :: bcv_snow, corr, t_soil_out
		

		! Calculates soil surface temperature based on Toy et al. (1979)
		if (m == 12 .or. m == 1 .or. m == 2) then
			t_soil_out = t_avg*0.656d0 + 2.3967d0
		else if (m == 3 .or. m == 4 .or. m == 5) then
			t_soil_out = t_avg*1.052d0 + 1.0239d0
		else if (m == 6 .or. m == 7 .or. m == 8) then
			t_soil_out = t_avg*0.856d0 + 6.3928d0
		else
			t_soil_out = t_avg*1.023d0 + 1.2856d0
		end if
		
		!snow correction
		bcv_snow = snow / ( snow + exp( 6.055d0 * snow - 0.3002d0 * snow ) )
		
		corr = min( (1.d0 - bcv_snow) ,Exp( -0.06d0 * lai ) )
		
		!correct for LAI/snow
		if (t_soil_out > t_avg) then
			t_soil = t_soil_out * corr
		end if

	end subroutine soil_temperature_toy


	! Calculate soil matric potential
	subroutine calc_psi_soil(asw, maxasw, sc, psi)
		
		real(kind=8), intent(in) :: asw, maxasw
		integer, intent(in) :: sc
			
		real(kind=8) :: psi, a, n
		! Parameters from Schaap, M. G., & Van Genuchten, M. T. (2006). 
			
		if (sc .eq. int(1)) then
			a = 0.02630268d0	
			n = 2.233572223d0	
		else if (sc .eq. int(2)) then
			a = 0.040738028d0	
			n = 1.191242008d0	
		else if (sc .eq. int(3)) then
			a = 0.012022644d0		
			n = 1.377209469d0	
		else
			a = 0.011220185d0	
			n = 1.300169578d0	
		endif
			
		! Parameters from Lu, N., Godt, J. W., & Wu, D. T. (2010). 
		! if (sc .eq. int(1)) then
			! a = 0.3d0
			! n = 3.d0
		! else if (sc .eq. int(2) .or. sc .eq. int(3)) then
			! a = 0.05d0
			! n = 2.5d0
		! else
			! a = 0.01d0
			! n = 1.8d0
		! endif
		
		psi = -(((1.d0/ (asw/maxasw))**(n/(n-1)) - 1.d0)**(1.d0/n) / a) / 1000.d0
		
		if (psi < -10.d0) then
			psi = -10.d0
		else if (psi > -0.033d0) then
			psi = -0.033d0
		end if
			
	end subroutine calc_psi_soil


	! Nutrient competition

	subroutine compute_fertility(Nav, Un, Up, root_dist, root_area, n_dist, fertility, n_sp)
		
		integer, intent(in) :: n_sp
		integer :: i
		real(kind=8), intent(in) :: Nav,  n_dist
		real(kind=8), dimension(n_sp), intent(in) :: Un, Up, root_dist, root_area 
		real(kind=8), dimension(n_sp), intent(out) :: fertility
		real(kind=8), dimension(n_sp) :: Nav_sp, Nav_ratio, fr_n_sp, fr_p_sp
		
		!Partition available nitrogen according to  root distribution, N distribution and specific root area/length (1 m of soil)
		!Root distribution is computed based on a power law (1-gamma^z) in accordance with 4C model implementation and nitrogen distribution
		!is computed based on an exponential function (1-exp(-kz)) based on the ICP soil survey shares of nitrogen at 30 cm and 100 cm soil depth (k=0.031).
		!Both the product of the function are integrated up to 1m soil depth to calculate the ratio of available nutrients for each species

		!Nav_ratio(:) = root_area(:) * ((root_dist(:)**100.d0 * n_dist**100.d0 + 100.d0 * (log(root_dist(:)) + log(n_dist)) - 1.d0) / &
		!			 (log(root_dist(:)) + log(n_dist)) + (1.d0 - root_dist(:)**100.d0)/log(root_dist(:)) + &
		!			 (1.d0 - n_dist**100.d0)/log(n_dist))

		Nav_ratio(:) = root_area(:) * ( 100.d0 + (-1.d0 + exp(-100.d0*n_dist))/n_dist + (1-exp(-100*n_dist)*root_dist(:)**100.d0)/ &
						(n_dist - log(root_dist(:))) + (1.d0 - root_dist(:)**100.d0)/log(root_dist(:)) )
		
		Nav_ratio(:) = Nav_ratio(:) / sum(Nav_ratio(:))
		
		!Define available N for each species
		Nav_sp(:) = Nav * Nav_ratio(:)
		
		
		!Get fertility rating
		do i=1, n_sp	
			!Nitrogen driven fertility ratio
			if (Un(i) == 0.d0) then
				fr_n_sp(i) = 1.d0
			else
				fr_n_sp(i) = Nav_sp(i) / Un(i)
			end if
			
			if (fr_n_sp(i) > 1.d0) then
				fr_n_sp(i) = 1.d0
			else if (fr_n_sp(i) < 0.d0) then
				fr_n_sp(i) = 0.d0
			else if (isnan(fr_n_sp(i))) then
				fr_n_sp(i) = 1.d0
			end if	
			
			fertility(i) = fr_n_sp(i)
			
		end do
		
	end subroutine compute_fertility



!-------------------------------------------------------------------------------------

end module soil_cnp_subroutines