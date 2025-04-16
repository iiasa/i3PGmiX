module gpp_calc

use gpp_subroutines

! P-model GPP calculation according to the rpmodel package (Stocker et al. 2020 - https://github.com/geco-bern/rpmodel) and
! P-hydro model  GPP calculation according to Joshi et al. (2023) - https://github.com/jaideep777/phydro

IMPLICIT NONE
public :: rpmodel, rphmodel, default_gpp

contains

subroutine default_gpp (gpp, eps , apar, n)

	INTEGER, INTENT(IN) ::   n
	REAL (kind = 8), DIMENSION(n), INTENT(IN) ::   eps , apar 
	REAL (kind = 8), DIMENSION(n), INTENT(OUT) :: gpp
	
	gpp(:) = eps(:) * apar(:) / 100.d0

end subroutine default_gpp

subroutine rpmodel ( tc, vpd, co2, fapar, ppfd, patm, elev, kphio, beta, soilm, meanalpha, apar_soilm, bpar_soilm, c4, &
                        method_jmaxlim, do_ftemp_kphio, do_soilmstress, out)

	REAL (kind = 8), INTENT(IN) ::   tc,vpd, co2, fapar, ppfd, patm, elev, kphio, beta, soilm, meanalpha, apar_soilm, bpar_soilm
	INTEGER, INTENT(IN) ::   method_jmaxlim
	LOGICAL, INTENT(IN) ::   c4, do_ftemp_kphio, do_soilmstress
	REAL (kind = 8) :: c_molmass, kPo, kTo, rd_to_vcmax, ns, ns25, ns_star, ftemp25_inst_vcmax, vcmax25_unitiabs, &
			   ftemp_inst_rd_v,rd_unitiabs, iabs, gpp, vcmax, vcmax25, rd, fact_jmaxlim, jmax, ftemp25_inst_jmax, &
			   jmax25, a_j, a_c, assim, gs, ca, ci, gammastar, iwue, kmm, soilmstress, beta_v, kphio_v, patm_v
	REAL (kind = 8), DIMENSION(5) :: out_optchi, out_lue_vcmax
	REAL (kind = 8), DIMENSION(22), INTENT(OUT) :: out
	!REAL (kind = 8) :: rpmodel(16)

	out(:) = 0.d0

	if (patm == -1.d0) then
		patm_v = calc_patm( elev )
	else 
		patm_v = patm
	end if
	  
	if (c4 .eqv. .TRUE.) then
		beta_v = 146.d0/9.d0
	else
		beta_v = 146.d0
	end if


	!---- Fixed parameters--------------------------------------------------------
	c_molmass = 12.0107d0  ! molecular mass of carbon (g)
	  
	kPo = 101325.d0       ! standard atmosphere, Pa (Allen, 1973)
	kTo = 25.d0           ! base temperature, deg C (Prentice, unpublished)
	rd_to_vcmax = 0.015d0  ! Ratio of Rdark to Vcmax25, number from Atkin et al., 2015 for C3 herbaceous


	!---- Temperature dependence of quantum yield efficiency----------------------
	!! 'do_ftemp_kphio' is not actually a stress function, but is the temperature-dependency of
	!! the quantum yield efficiency after Bernacchi et al., 2003 PCE

	if (c4 .eqv. .TRUE.) then
		kphio_v = 1.d0
	else
		if (do_ftemp_kphio .eqv. .TRUE.) then
			if (do_soilmstress .eqv. .TRUE.) then
				kphio_v = 0.087182d0
			else
				kphio_v = 0.087182d0
			end if
		else
			kphio_v = kphio
		end if

	end if

	if (do_ftemp_kphio .eqv. .TRUE.) then
		kphio_v = ftemp_kphio( tc, c4 ) * kphio_v
	else 
		if (c4 .eqv. .TRUE.) then 
			kphio_v = ftemp_kphio( 15.d0, c4 ) * kphio_v
		end if
	end if


	!---- soil moisture stress as a function of soil moisture and mean alpha -----
	if (do_soilmstress .eqv. .TRUE.) then
		soilmstress = calc_soilmstress( soilm, apar_soilm, bpar_soilm )
		soilmstress = soilm
	else 
		soilmstress = 1.d0
	end if

	!---- Photosynthesis parameters depending on temperature, pressure, and CO2. -
	!! ambient CO2 partial pression (Pa)
	ca = co2_to_ca( co2, patm_v )

	!! photorespiratory compensation point - Gamma-star (Pa)
	gammastar = calc_gammastar( tc, patm_v )

	!! Michaelis-Menten coef. (Pa)
	kmm = calc_kmm( tc, patm_v )

	!! viscosity correction factor = viscosity( temp, press )/viscosity( 25 degC, 1013.25 Pa)
	ns      = viscosity_h2o( tc, patm_v )  ! Pa sc4, 1.0,
	ns25    = viscosity_h2o( kTo, kPo )  ! Pa s
	ns_star = ns / ns25  ! (unitless)

	!!----Optimal ci -------------------------------------------------------------
	!! The heart of the P-model: calculate ci:ca ratio (chi) and additional terms
	out_optchi = optimal_chi( kmm, gammastar, ns_star, ca, vpd, beta_v, c4 )
	  
	!! leaf-internal CO2 partial pressure (Pa)
	ci = out_optchi(2) * ca

	!---- Corrolary preditions ---------------------------------------------------
	!! intrinsic water use efficiency (in Pa)
	iwue = ( ca - ci ) / 1.6d0

	!---- Vcmax and light use efficiency -----------------------------------------
	! Jmax limitation comes in only at this step
	if (c4 .eqv. .TRUE.) then

		out_lue_vcmax = lue_vcmax_c4(kphio_v,c_molmass,soilmstress)

	else if (method_jmaxlim==1) then

		!! apply correction by Jmax limitation
		out_lue_vcmax = lue_vcmax_wang17( out_optchi, kphio_v, c_molmass, soilmstress)

	else if (method_jmaxlim==2) then

		out_lue_vcmax = lue_vcmax_smith19(out_optchi, kphio_v, c_molmass, soilmstress)

	else if (method_jmaxlim==3) then

		out_lue_vcmax = lue_vcmax_none(out_optchi, kphio_v, c_molmass, soilmstress)

	else 

		ERROR STOP "rpmodel(): argument method_jmaxlim not idetified."

	end if

	!---- Corrolary preditions ---------------------------------------------------
	! Vcmax25 (vcmax normalized to 25 deg C)
	ftemp25_inst_vcmax  = ftemp_inst_vcmax( tc )
	vcmax25_unitiabs  = out_lue_vcmax(3) / ftemp25_inst_vcmax

	!! Dark respiration at growth temperature
	ftemp_inst_rd_v = ftemp_inst_rd( tc )
	rd_unitiabs  = rd_to_vcmax * (ftemp_inst_rd_v / ftemp25_inst_vcmax) * out_lue_vcmax(3)

	!---- Quantities that scale linearly with absorbed light ---------------------
	iabs = fapar * ppfd

	! Gross Primary Productivity
	gpp = iabs * out_lue_vcmax(2)   ! in g C m-2 s-1

	! Vcmax per unit ground area is the product of the intrinsic quantum
	! efficiency, the absorbed PAR, and 'n'
	vcmax = iabs * out_lue_vcmax(3)

	!! (vcmax normalized to 25 deg C)
	vcmax25 = iabs * vcmax25_unitiabs

	!! Dark respiration
	rd = iabs * rd_unitiabs

	! Jmax using again A_J = A_C, derive the "Jmax limitation factor" 
	fact_jmaxlim = vcmax * (ci + 2.d0 * gammastar) / (kphio_v * iabs * (ci + kmm))
	  
	! use definition of Jmax limitation factor (L in Eq. 13) and solve for Jmax.
	jmax = 4.d0 * kphio_v * iabs / sqrt( (1.d0/fact_jmaxlim)**2.d0 - 1.d0 )
	  
	! !! Alternatively, Jmax can be calculated from Eq. F10 in Stocker et al., 2020
	! kc <- 0.41
	! jmax_alt <- 4.0 * kphio_v * iabs * sqrt((out_optchi$mj / kc)**(2/3) - 1.0)
	! fact_jmaxlim_alt <- 1.0 / sqrt(1 + (4.0 * kphio * iabs / jmax_alt)**2)
	  
	ftemp25_inst_jmax = ftemp_inst_jmax( tc )
	jmax25 = jmax / ftemp25_inst_jmax

	!! Test: at this stage, verify if A_J = A_C
	if (c4 .eqv. .TRUE.) then
		a_j = kphio_v * iabs * out_optchi(4) * fact_jmaxlim
		a_c = vcmax * out_optchi(3)
	else 
		a_j = kphio_v * iabs * (ci - gammastar)/(ci + 2.d0 * gammastar) * fact_jmaxlim
		a_c = vcmax * (ci - gammastar) / (ci + kmm)
	end if
	  
	if (abs(a_j - a_c) > 0.001d0) then
		! Print "rpmodel(): light and Rubisco-limited assimilation rates are not identical"
	end if

	! Assimilation is not returned because it should not be confused with what 
	! is usually measured should use instantaneous assimilation for comparison to
	! measurements. This is returned by inst_rpmodel().
	if (a_j < a_c) then
		assim = a_j
	else
		assim = a_c
	end if

	if (abs(assim - gpp / c_molmass) > 0.001d0) then
		! Print "rpmodel(): Assimilation and GPP are not identical"
	end if	

	!! average stomatal conductance
	gs = 1.6e6 * assim / ((ca - ci) / (( 1.0e-6 ) * patm_v)) * (8.31432d0 * (tc + 273.15d0)) / patm_v ! convert to m/month
  
	!! construct list for output
	out(1) = gpp
	out(2) = ca
	out(3) = gammastar
	out(4) = ns_star
	out(5) = out_optchi(2)
	out(6) = out_optchi(1)
	out(7) = out_optchi(4)
	out(8) = out_optchi(3)
	out(9) = ci
	out(10) = 1.6d0 * assim / (ca - ci)
	out(11) = gs 
	out(12) = vcmax
	out(13) = vcmax25
	out(14) = jmax
	out(15) = jmax25
	out(16) = out_optchi(5)
	out(17) = out_lue_vcmax(1)
	out(18) = out_lue_vcmax(2)
	out(19) = out_lue_vcmax(3)
	out(20) = out_lue_vcmax(4)
	out(21) = out_lue_vcmax(5)
	out(22) = rd * c_molmass

	!rpmodel(:) = out(:)

end subroutine rpmodel


subroutine rphmodel( tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, rdark, par_plant, par_cost, out)

	REAL (kind = 8), INTENT(IN) ::   tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, rdark
	REAL (kind = 8), DIMENSION(3), INTENT(IN) ::  par_plant
	REAL (kind = 8), DIMENSION(2), INTENT(IN) ::  par_cost
	integer :: iflag


	REAL (kind = 8) :: kPo, kTo, rd_to_vcmax, ns, ns25, ns_star, ftemp25_inst_vcmax, vcmax25_unitiabs, &
			   ftemp_inst_rd_v,rd_unitiabs, iabs, gpp, vcmax, vcmax25, rd, fact_jmaxlim, jmax, ftemp25_inst_jmax, &
			   jmax25, a_j, a_c, assim, ca, ci, gammastar, iwue, kmm, soilmstress, beta_v, kphio_v, patm_v, &
			   viscosity_water, density_water, chi, delta, u, dpsi, gs, J, g, c_molmass, dpsi_bounds, gamma, kappa, &
			   phi0, xzero, fzero
			   
	REAL (kind = 8), DIMENSION(10) :: out
	
	c_molmass = 12.0107d0  ! molecular mass of carbon (g)
	
	patm_v = calc_patm( elv )
    kmm = calc_kmm(tc, patm_v)
    gammastar = calc_gammastar(tc, patm_v)
    phi0 = kphio*ftemp_kphio(tc, .FALSE.)
    iabs = ppfd*fapar
    ca = co2*patm_v*1e-6
    delta = rdark	
	
	viscosity_water = viscosity_h2o(tc, patm_v)
    density_water = density_h2o(tc, patm_v)

	dpsi_bounds = calc_dpsi_bound(psi_soil, ca, gammastar, patm_v, vpd, par_plant, viscosity_water, &
										 density_water, par_cost,kmm, phi0, iabs, delta)
	

	!u = uniroot(f_f, dpsi_bounds*0.001d0, dpsi_bounds*0.99d0,0.0001d0)
	!dpsi_bounds = 2.d0
	call root_scalar('zhang',f_f,dpsi_bounds*0.001d0, dpsi_bounds*0.99d0,xzero, fzero,iflag)
	
	dpsi = xzero
	!dpsi = 1.d0 !0.01140841d0
	
	chi = calc_x_from_dpsi(dpsi, psi_soil, par_plant, gammastar, patm_v, vpd, viscosity_water, &
										density_water, kmm, ca, delta, par_cost)
		
	gs = calc_gs(dpsi, psi_soil, par_plant, vpd, patm_v, viscosity_water, density_water)
	  
	J = calc_J(gs, chi, gammastar, ca, kmm, delta, patm_v)	
	
	Jmax = calc_jmax_from_J(J, phi0, iabs)
	  
	ca = ca/patm_v*1e6
	kmm = kmm/patm_v*1e6
	g = gammastar/patm_v*1e6
	Vcmax = (J/4.d0)*(chi*ca + kmm)/(chi*ca + 2.d0*g)
	!Vcmax = J*(chi*ca + kmm)/(chi*ca + 2.d0*g)
  
	assim = gs * ca*(1.d0-chi)
	
	gpp = assim * c_molmass * 0.0864d0
	
	gamma = gammastar / ca
	kappa = kmm / ca
	
	rd_unitiabs  = 0.015d0 * (ftemp_inst_rd( tc ) / ftemp_inst_vcmax( tc )) * Vcmax * 0.0864d0 * c_molmass 
	
	!rd = iabs * rd_unitiabs * c_molmass 
	
	out(1) = gpp 
    out(2) = dpsi
    out(3) = gs * (8.31432d0 * (tc + 273.15d0)) / patm_v ! convert to m/s
    out(4) = rd
    out(5) = transpiration(dpsi, psi_soil, par_plant, vpd, patm_v, viscosity_water, density_water) * 1556.52d0 ! convert to mm/day
    out(6) = chi
    out(7) = gammastar
    out(8) = Vcmax
    out(9) = J
    out(10) = Jmax
	
	contains
	
	REAL (kind = 8) FUNCTION f_f(dpsi)
	
		real(kind=8), intent(in) :: dpsi
		
		f_f = dFdx(dpsi, psi_soil, gammastar, ca, kmm, delta, patm_v, vpd, &
						phi0, iabs, par_plant,viscosity_water, density_water, par_cost)

	end function f_f
	
	
end subroutine rphmodel

end module gpp_calc