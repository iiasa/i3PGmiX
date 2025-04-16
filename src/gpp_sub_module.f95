module gpp_subroutines

use root_module

IMPLICIT NONE
public :: ftemp_arrh, calc_patm, density_h2o, QUADP, QUADM, omega, calc_kmm, &
			calc_gammastar, ftemp_inst_rd, ftemp_kphio, calc_soilmstress, &
			ftemp_inst_vcmax, ftemp_inst_jmax, viscosity_h2o, co2_to_ca, optimal_chi, &
			lue_vcmax_wang17, lue_vcmax_smith19, lue_vcmax_none, lue_vcmax_c4, &
			calc_gs, calc_gsprime, calc_x_from_dpsi, calc_J, calc_djmax_dJ, &
			scale_conductivity, P, Pprime, Pprimeprime, integral_P_ana, gammainc, &
			calc_dJ_dchi, calc_jmax_from_J, f1, f2, dFdx, chi_jmax_limited, &
			calc_dpsi_bound



contains

!-----------------------------------------------------------------------------------------
!                               P-Model auxiliary functions
!-----------------------------------------------------------------------------------------

! Calculates the Arrhenius-type temperature response
! tk Temperature (Kelvin)
! dha Activation energy (J mol-1)
! tkref Reference temperature (Kelvin)
REAL (kind = 8) FUNCTION ftemp_arrh( tk, dha )
	REAL (kind = 8), INTENT(IN) :: tk, dha
	REAL (kind = 8), PARAMETER :: tkref = 298.15d0
	REAL (kind = 8), PARAMETER :: kR = 8.3145d0     ! Universal gas constant, J/mol/K
	REAL (kind = 8) :: ftemp 

	! Note that the following forms are equivalent:
	! ftemp = exp( dha * (tk - 298.15) / (298.15 * kR * tk) )
	! ftemp = exp( dha * (tc - 25.0)/(298.15 * kR * (tc + 273.15)) )
	! ftemp = exp( (dha/kR) * (1/298.15 - 1/tk) ) 
	  
	ftemp = exp( dha * (tk - tkref) / (tkref * kR * tk) )
	  
	ftemp_arrh = ftemp
	
END FUNCTION ftemp_arrh


! Calculates atmospheric pressure
! elev elevation in m a.s.l.
REAL (kind = 8) FUNCTION calc_patm( elv )

	REAL (kind = 8), INTENT(IN) :: elv
	REAL (kind = 8), PARAMETER :: patm0 = 101325d0   ! standard atmosphere a.s.l.
	REAL (kind = 8), PARAMETER :: kTo = 298.15d0    ! base temperature, K (Prentice, unpublished)
	REAL (kind = 8), PARAMETER :: kL  = 0.0065d0    ! gravitational acceleration, m/s^2 (Allen, 1973)
	REAL (kind = 8), PARAMETER :: kG  = 9.80665d0   ! gravitational acceleration, m/s^2 (Allen, 1973)
	REAL (kind = 8), PARAMETER :: kR  = 8.3145d0   ! universal gas constant, J/mol/K (Allen, 1973)
	REAL (kind = 8), PARAMETER :: kMa = 0.028963d0  ! molecular weight of dry air, kg/mol (Tsilingiris, 2008)
	REAL (kind=8) :: out

	! Convert elevation to pressure, Pa:
	out = patm0*(1.d0 - kL*elv/kTo)**(kG*kMa/(kR*kL))
	  
	calc_patm = out
	
END FUNCTION calc_patm


! Density of water
! tc air temperature deg Celsius
! p atmospheric pressure Pa
REAL (kind = 8) FUNCTION density_h2o( tc,p )

	REAL (kind = 8), INTENT(IN) :: tc, p
	REAL (kind = 8)  ::  my_lambda, po, vinf, pbar, v, rho

	! Calculate lambda, (bar cm^3)/g:
	my_lambda = 1788.316d0 + 21.55053d0*tc -0.4695911d0*tc**2.d0 + &
			(3.096363e-3)*tc**3.d0 -(7.341182e-6)*tc**4.d0
	  
	! Calculate po, bar
	po = 5918.499d0 + 58.05267d0*tc -1.1253317d0*tc**2.d0 + &
		(6.6123869e-3)*tc**3.d0 -(1.4661625e-5)*tc**4.d0
	  
	  ! Calculate vinf, cm^3/g
	vinf = 0.6980547d0 -(7.435626e-4)*tc + (3.704258e-5)*tc**2.d0 + &
		-(6.315724e-7)*tc**3.d0 + (9.829576e-9)*tc**4.d0 -(1.197269e-10)*tc**5.d0 + &
		(1.005461e-12)*tc**6.d0 -(5.437898e-15)*tc**7.d0 + (1.69946e-17)*tc**8.d0 + &
		-(2.295063e-20)*tc**9.d0
	  
	! Convert pressure to bars (1 bar <- 100000 Pa)
	pbar = (1e-5)*p
	  
	! Calculate the specific volume (cm^3 g^-1):
	v = vinf + my_lambda/(po + pbar)
	  
	! Convert to density (g cm^-3) -> 1000 g/kg; 1000000 cm^3/m^3 -> kg/m^3:
	rho = (1e3/v)
	  
	density_h2o = rho

END FUNCTION density_h2o


!Larger quadratic root
! a, b, c polynomial coefficients
REAL (kind = 8) FUNCTION QUADP( A, B ,C )

IMPLICIT NONE

	REAL (kind = 8), INTENT(IN) :: A, B, C
	REAL (kind = 8):: out

	if ((B**2.d0 - 4.d0*A*C) < 0.d0) then
		out = 0.d0
	else
		if (A==0.d0) then
			if (B==0.d0) then
				out = 0.d0
			else
				out = -C/B
			end if
		else
			out = (- B + sqrt(B**2.d0 - 4.d0*A*C)) / (2.d0*A)
		end if
	end if

	QUADP = out

END FUNCTION QUADP



!Minor quadratic root
! a, b, c polynomial coefficients
REAL (kind = 8) FUNCTION QUADM( A, B ,C )

IMPLICIT NONE

	REAL (kind = 8), INTENT(IN) :: A, B, C
	REAL (kind = 8):: out

	if ((B**2.d0 - 4.d0*A*C) < 0.d0) then
		out = 0.d0
	else
	  
		if (A==0.d0) then
			if (B==0.d0) then
				out = 0.0
			else
				out = -C/B
			end if
		else
			out = (- B - sqrt(B**2.d0 - 4.d0*A*C)) / (2.d0*A)
		end if
	end if

	QUADM = out

END FUNCTION QUADM

REAL (kind = 8) FUNCTION omega ( theta, c_cost, m )

REAL (kind = 8), INTENT(IN) :: theta, c_cost, m 
REAL (kind=8) :: cm, v, capP, aquad, bquad, cquad, m_star, omega_v

    cm = 4.d0 * c_cost / m                        ! simplification term for omega calculation
    v  = 1.d0/(cm * (1.d0 - theta * cm)) - 4.d0 * theta ! simplification term for omega calculation
    
    ! account for non-linearities at low m values
    capP = (((1.d0/1.4d0) - 0.7d0)**2.d0 / (1.d0-theta)) + 3.4d0
    aquad = -1.d0
    bquad = capP
    cquad = -(capP * theta)
    m_star = (4.d0 * c_cost) / QUADP(cquad,bquad,aquad) !polyroot( aquad, bquad, cquad )

    if (m < m_star) then
    	omega_v = -( 1.d0 - (2.d0 * theta) ) - sqrt( (1.d0 - theta) * v)
    else
	omega_v = -( 1.d0 - (2.d0 * theta))  + sqrt( (1.d0 - theta) * v)
    end if
    
    omega = omega_v

END FUNCTION omega

REAL (kind = 8) FUNCTION calc_kmm( tc, patm )

	REAL (kind = 8), INTENT(IN) :: tc, patm

	! Declare constants
	REAL (kind = 8), PARAMETER :: dhac = 79430.d0     ! (J/mol) Activation energy, Bernacchi et al. (2001) 
	REAL (kind = 8), PARAMETER :: dhao = 36380.d0      ! (J/mol) Activation energy, Bernacchi et al. (2001)
	REAL (kind = 8), PARAMETER :: kco = 2.09476e5  ! (ppm) O2 partial pressure, Standard Atmosphere
	REAL (kind = 8), PARAMETER :: kc25 = 39.97d0   ! Pa, value based on Bernacchi et al. (2001)
	REAL (kind = 8), PARAMETER :: ko25 = 27480.d0   ! Pa, value based on Bernacchi et al. (2001)

	REAL (kind=8) :: tk, po, kmm, kc, ko

	! conversion to Kelvin
	tk = tc + 273.15d0
	  
	kc = ftemp_arrh( tk, dhac )
	kc = kc * kc25 
	ko = ftemp_arrh( tk, dhao )
	ko = ko * ko25
	  
	po  = kco * (1e-6) * patm         ! O2 partial pressure
	kmm = kc * (1.d0 + po/ko)

	calc_kmm = kmm
	  
END FUNCTION calc_kmm


! Calculates the CO2 compensation point
! tc temperature in deg. Celsius
! patm atmospheric pressure in Pa
REAL (kind = 8) FUNCTION calc_gammastar ( tc, patm )

	REAL (kind = 8), INTENT(IN) :: tc, patm
	REAL (kind = 8), PARAMETER :: dha = 37830.d0 ! (J/mol) Activation energy, Bernacchi et al. (2001)
	REAL (kind = 8), PARAMETER :: gs25_0 = 4.332d0
	REAL (kind=8) :: gammastar

	gammastar = gs25_0 * patm / calc_patm(0.d0) * ftemp_arrh( (tc + 273.15d0), dha )

	calc_gammastar = gammastar

END FUNCTION calc_gammastar


! Calculates the temperature response of dark respiration
! tc temperature in deg Celsius
REAL (kind = 8) FUNCTION ftemp_inst_rd( tc )

	REAL (kind = 8), INTENT(IN) :: tc
	REAL (kind = 8), PARAMETER :: apar = 0.1012d0	!loal parameter
	REAL (kind = 8), PARAMETER :: bpar = 0.0005d0	!loal parameter
	REAL (kind = 8) :: fr

	fr = exp( apar * (tc - 25.d0) - bpar * (tc**2.d0 - 25.d0**2.d0) )
	ftemp_inst_rd = fr

END FUNCTION ftemp_inst_rd


!Calculates the temperature dependence of the quantum yield efficiency
! tc air temperature deg Celsius
REAL (kind = 8) FUNCTION ftemp_kphio ( tc, c4 )

	REAL (kind = 8), INTENT(IN) :: tc
	REAL (kind = 8) :: ftemp
	LOGICAL :: c4

	if (c4 .eqv. .TRUE.) then
		ftemp = -0.064d0 + 0.03d0 * tc - 0.000464d0 * tc**2.d0     
	else
		ftemp = 0.352d0 + 0.022d0 * tc - 3.4e-4 * tc**2.d0
	end if

	if (ftemp < 0.d0) then
		ftemp = 0.d0
	end if

	ftemp_kphio = ftemp

END FUNCTION ftemp_kphio


!Calculates an empirical soil moisture stress factor
REAL (kind = 8) FUNCTION calc_soilmstress ( soilm, apar_soilm, bpar_soilm )

	REAL (kind = 8), INTENT(IN) :: soilm, apar_soilm, bpar_soilm
	REAL (kind = 8), PARAMETER ::  meanalpha = 1.d0
	REAL (kind = 8), PARAMETER ::  x0 = 0.d0
	REAL (kind = 8), PARAMETER ::  x1 = 0.6d0
	REAL (kind = 8) ::  y0, beta, outstress

	!apar_soilm = 0.0
	!bpar_soilm = 0.685

	y0 = (apar_soilm + bpar_soilm * meanalpha)
	beta = (1.d0 - y0) / (x0 - x1)**2
	outstress = 1.d0 - beta * ( soilm - x1 )**2

	! bound between 0 and 1
	outstress = min(max(outstress, 0.d0), 1.d0)
	! and set to 1.0 above soil moisture threshold x1.
	if (soilm > x1) then
		outstress = 1.d0
	end if

	calc_soilmstress = outstress

END FUNCTION calc_soilmstress



!Calculates the instantaneous temperature response of Vcmax
! tcleaf Leaf temperature
! tcgrowth Growth temperature
! tcref Reference temperature
REAL (kind = 8) FUNCTION ftemp_inst_vcmax ( tcleaf )

	REAL (kind = 8), INTENT(IN) :: tcleaf
	REAL (kind = 8), PARAMETER :: tcref = 25.d0
	REAL (kind = 8), PARAMETER :: Ha    = 71513.d0  ! activation energy (J/mol)
	REAL (kind = 8), PARAMETER :: Hd    = 200000.d0 ! deactivation energy (J/mol)
	REAL (kind = 8), PARAMETER :: Rgas  = 8.3145d0 ! universal gas constant (J/mol/K)
	REAL (kind = 8), PARAMETER :: a_ent = 668.39d0 ! offset of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K)
	REAL (kind = 8), PARAMETER :: b_ent = 1.07d0   ! slope of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K^2)
	REAL (kind = 8) :: tcgrowth, tkref, tkleaf, dent, fva, fvb, fv

	tcgrowth = tcleaf
	tkref = tcref + 273.15  ! to Kelvin
	tkleaf = tcleaf + 273.15 ! conversion of temperature to Kelvin

	! calculate entropy following Kattge and Knorr (2007)
	dent = a_ent - b_ent * tcgrowth   ! 'tcgrowth' corresponds to 'tmean' in Nicks, 'tc25' is 'to' in Nick's
	fva = ftemp_arrh( tkleaf, Ha)
	fvb = (1.d0 + exp((tkref * dent - Hd)/(Rgas * tkref))) / (1.d0 + exp( (tkleaf * dent - Hd)/(Rgas * tkleaf) ))
	fv  = fva * fvb

	ftemp_inst_vcmax = fv

END FUNCTION ftemp_inst_vcmax



!Calculates the instantaneous temperature response of Jmax
! tcleaf Leaf temperature
! tcgrowth Growth temperature
! tcref Reference temperature
REAL (kind = 8) FUNCTION ftemp_inst_jmax( tcleaf )

	REAL (kind = 8), INTENT(IN) :: tcleaf
	REAL (kind = 8), PARAMETER :: tcref = 25.d0
	REAL (kind = 8), PARAMETER :: Ha = 49884.d0  ! activation energy (J/mol)
	REAL (kind = 8), PARAMETER :: Hd = 200000.d0 ! deactivation energy (J/mol)
	REAL (kind = 8), PARAMETER :: Rgas = 8.3145d0 ! universal gas constant (J/mol/K)
	REAL (kind = 8), PARAMETER :: a_ent = 659.7d0 ! offset and slope (a / b) of entropy vs. temperature relationship 
	REAL (kind = 8), PARAMETER :: b_ent = 0.75d0   ! offset and slope (a / b) of entropy vs. temperature relationship 
	REAL (kind = 8) :: tcgrowth, tkref, tkleaf, dent, fva, fvb, fv

	tcgrowth = tcleaf
	tkref = tcref + 273.15d0
	tkleaf = tcleaf + 273.15d0

	! calculate entropy following Kattge & Knorr (2007), negative slope and 
	! y-axis intersect is when expressed as a function of temperature in
	! degrees Celsius, not Kelvin !!!
	  
	! 'tcgrowth' corresponds to 'tmean' in Nicks, 'tc25' is 'to' in Nick's
	dent = a_ent - b_ent * tcgrowth   
	fva = ftemp_arrh( tkleaf, Ha )
	fvb = (1.d0 + exp( (tkref * dent - Hd)/(Rgas * tkref) ) ) / (1.d0 + exp( (tkleaf * dent - Hd)/(Rgas * tkleaf) ) )
	fv  = fva * fvb

	ftemp_inst_jmax = fv

END FUNCTION ftemp_inst_jmax


! Viscosity of water
! tc air temperaturwe deg Celsius
! p atmospheric pressure Pa
REAL (kind = 8) FUNCTION viscosity_h2o( tc,p )

	REAL (kind = 8), INTENT(IN) :: tc, p
	REAL (kind = 8), PARAMETER :: tcref = 25.d0
	REAL (kind = 8), PARAMETER :: tk_ast = 647.096d0    ! Kelvin
	REAL (kind = 8), PARAMETER :: rho_ast = 322.d0      ! kg/m^3
	REAL (kind = 8), PARAMETER :: mu_ast = 1e-6       ! Pa s
	REAL (kind = 8) :: rho, tbar, tbarx, tbar2, tbar3, rbar, mu, mu0, mu1, ctbar, coef1, coef2, mu_bar
	REAL (kind = 8), DIMENSION(7,6) :: h_array
	INTEGER :: i,j

	! Get the density of water, kg/m^3
	rho = density_h2o(tc, p)
	  
	! Calculate dimensionless parameters:
	tbar  = (tc + 273.15d0)/tk_ast
	tbarx = tbar**(0.5d0)
	tbar2 = tbar**2.d0
	tbar3 = tbar**3.d0
	rbar  = rho/rho_ast
  
	! Calculate mu0 (Eq. 11 & Table 2, Huber et al., 2009):
	mu0 = 1.67752d0 + 2.20462d0/tbar + 0.6366564d0/tbar2 - 0.241605d0/tbar3
	mu0 = 1e2*tbarx/mu0

	! Create Table 3, Huber et al. (2009):
	h_array =reshape((/0.520094d0, 0.0850895d0, -1.08374d0, -0.289555d0, 0.d0, 0.d0, &
					0.222531d0, 0.999115d0, 1.88797d0, 1.26613d0, 0.d0, 0.120573d0, &
					-0.281378d0, -0.906851d0, -0.772479d0, -0.489837d0, -0.25704d0, 0.d0, &
					0.161913d0,  0.257399d0, 0.d0, 0.d0, 0.d0, 0.d0, &
					-0.0325372d0, 0.d0, 0.d0, 0.0698452d0, 0.d0, 0.d0, &
					0.d0, 0.d0, 0.d0, 0.d0, 0.00872102d0, 0.d0, &
					0.d0, 0.d0, 0.d0, -0.00435673d0, 0.d0, -0.000593264d0/), &
					shape(h_array), order=(/2,1/))


	! Calculate mu1 (Eq. 12 & Table 3, Huber et al., 2009):
	mu1 = 0.d0
	ctbar = (1.d0/tbar) - 1.d0
	! print(paste("ctbar",ctbar))
	! for i in xrange(6):
	do i=1, 6
		coef1 = ctbar**(i-1.d0)  
		coef2 = 0.d0
		do j=1,7
			coef2 = coef2 + h_array(j,i) * (rbar - 1.d0)**(j-1)
		end do
		mu1 = mu1 + coef1 * coef2
	end do
	mu1 = exp( rbar * mu1 )
  
	! Calculate mu_bar (Eq. 2, Huber et al., 2009)
	!   assumes mu2 = 1
	mu_bar = mu0 * mu1
	  
	! Calculate mu (Eq. 1, Huber et al., 2009)
	mu = mu_bar * mu_ast    ! Pa s

	viscosity_h2o = mu

END FUNCTION viscosity_h2o



!CO2 partial pressure
! co2 Atmospheric CO2 concentration (ppm)
! patm Atmospheric pressure (Pa)
REAL (kind = 8) FUNCTION co2_to_ca( co2,patm )

	REAL (kind = 8), INTENT(IN) :: co2, patm
	REAL (kind = 8) :: ca
	 
	! Pa, atms. CO2
	ca   = ( 1.0e-6 ) * co2 * patm
	co2_to_ca = ca

END FUNCTION co2_to_ca


!Calculate optimal chi
! kmm : Pa, Michaelis-Menten coeff.
! gammastar
! ns_star  : (unitless) viscosity correction factor for water
! ca
! vpd : Pa, vapor pressure deficit
! beta
! c4

FUNCTION optimal_chi ( kmm, gammastar, ns_star, ca, vpd, beta, c4 )

	REAL (kind = 8), INTENT(IN) :: kmm, gammastar, ns_star, ca, vpd, beta
	REAL (kind = 8), DIMENSION(5) :: out
	LOGICAL, INTENT(IN) :: c4
	REAL (kind = 8)  ::  xi, chi, gamma, kappa, mj, mc, mjoc, vpd_aux
	REAL (kind = 8)  ::  optimal_chi(5)

	if (vpd < 0.d0) then
		vpd_aux = 0.d0
	else
		vpd_aux = vpd
	end if 

	! leaf-internal-to-ambient CO2 partial pressure (ci/ca) ratio
	xi  = sqrt( (beta * ( kmm + gammastar ) ) / ( 1.6 * ns_star ) )
	chi = gammastar / ca + ( 1.d0 - gammastar / ca ) * xi / ( xi + sqrt(vpd_aux) )

	if (c4 .eqv. .TRUE.) then
		out(1) = xi
		out(2) = chi
		out(3) = 1.d0
		out(4) = 1.d0
		out(5) = 1.d0
	else 
		!! alternative variables
		gamma = gammastar / ca
		kappa = kmm / ca
		
		! use chi for calculating mj
		mj = (chi - gamma) / (chi + 2.d0 * gamma)
		
		! mc
		mc = (chi - gamma) / (chi + kappa)
		
		! mj:mv
		mjoc = (chi + kappa) / (chi + 2.d0 * gamma)
		
		! format output list
		out(1) = xi
		out(2) = chi
		out(3) = mc
		out(4) = mj
		out(5) = mjoc
	end if

	optimal_chi(:) = out(:)

END FUNCTION optimal_chi



! Light use efficiency based on Wang et al. 2017
! out_optchi
! kphio
! c_molmass
! soilmstress
FUNCTION lue_vcmax_wang17 ( out_optchi, kphio, c_molmass, soilmstress )

	REAL (kind = 8), INTENT(IN) :: kphio, c_molmass, soilmstress
	REAL (kind = 8), DIMENSION(5), INTENT(IN) :: out_optchi
	REAL (kind = 8), PARAMETER :: kc= 0.41
	REAL (kind=8) :: len, tmp, mprime, lue, vcmax_unitiabs, omega, omega_star
	REAL (kind = 8), DIMENSION(5) :: out
	REAL (kind = 8) :: lue_vcmax_wang17(5)

	! ## Following eq. 17 in Stocker et al., 2020 GMD
	! tmp <- 1.0 - (kc / out_optchi$mj)**(2.0/3.0)
	! mprime <- ifelse(tmp > 0, out_optchi$mj * sqrt(tmp), NA)  # avoid square root of negative number

	! original code - is equivalent to eq. 17-based formulation above
	tmp = out_optchi(4)**2.d0 - kc**(2.d0/3.d0) * (out_optchi(4)**(4.d0/3.d0))

	if (tmp > 0.d0) then
		mprime = sqrt(tmp)
	else
		mprime = 0.d0
	end if

	   
	! Light use efficiency (gpp per unit absorbed light)
	lue = kphio * mprime * c_molmass * soilmstress
		
	! Vcmax normalised per unit absorbed PPFD (assuming iabs=1), with Jmax limitation
	vcmax_unitiabs = soilmstress * kphio * mprime / out_optchi(3)
		
	! complement for non-smith19
	omega      = 0.d0
	omega_star = 0.d0

	out(1) = mprime
	out(2) = lue
	out(3) = vcmax_unitiabs
	out(4) = omega
	out(5) = omega_star

	lue_vcmax_wang17(:) = out(:)

END FUNCTION lue_vcmax_wang17



! Light use efficiency based on Smith et al. 2019
! out_optchi
! kphio
! c_molmass
! soilmstress
FUNCTION lue_vcmax_smith19( out_optchi, kphio, c_molmass, soilmstress )

	REAL (kind = 8), INTENT(IN) :: kphio, c_molmass, soilmstress
	REAL (kind = 8), DIMENSION(5), INTENT(IN) :: out_optchi
	REAL (kind = 8), PARAMETER :: theta = 0.85d0    ! should be calibratable?
	REAL (kind = 8), PARAMETER :: c_cost = 0.05336251d0
	REAL (kind=8) :: len, tmp, mprime, lue, vcmax_unitiabs, omega_v, omega_star
	REAL (kind = 8), DIMENSION(5) :: out
	REAL (kind = 8) :: lue_vcmax_smith19(5)

	! factors derived as in Smith et al., 2019
	omega_v = omega( theta , c_cost , out_optchi(4) )          ! Eq. S4
	omega_star = 1.d0 + omega_v - sqrt( (1.d0 + omega_v)**2.d0 - (4.d0 * theta * omega_v) )       ! Eq. 18
	  
	! Effect of Jmax limitation
	mprime = out_optchi(4) * omega_star / (8.d0 * theta)
	  
	! Light use efficiency (gpp per unit absorbed light)
	lue = kphio * mprime * c_molmass * soilmstress
	  
	! calculate Vcmax per unit aborbed light
	vcmax_unitiabs  = kphio * out_optchi(5) * omega_star / (8.d0 * theta) * soilmstress   ! Eq. 19
	  
	out(1) = mprime
	out(2) = lue
	out(3) = vcmax_unitiabs
	out(4) = omega_v
	out(5) = omega_star
	  
	lue_vcmax_smith19(:) = out(:)

END FUNCTION lue_vcmax_smith19



! Find roots of polynomial function
!REAL (kind = 8) FUNCTION polyroot ( a,b,c )

!END FUNCTION polyroot



! Do not include effect of Jmax limitation
! out_optchi
! kphio
! c_molmass
! soilmstress

FUNCTION lue_vcmax_none( out_optchi, kphio, c_molmass, soilmstress ) 

	REAL (kind = 8), INTENT(IN) :: kphio, c_molmass, soilmstress
	REAL (kind = 8), DIMENSION(5), INTENT(IN) :: out_optchi
	REAL (kind=8) :: lue, vcmax_unitiabs, omega_v, omega_star
	REAL (kind = 8), DIMENSION(5) :: out
	REAL (kind = 8) :: lue_vcmax_none(5)

	! Light use efficiency (gpp per unit absorbed light)
	lue = kphio * out_optchi(4) * c_molmass * soilmstress
	  
	! calculate Vcmax per unit aborbed light
	vcmax_unitiabs  = kphio * out_optchi(5) * soilmstress   ! Eq. 19

	omega_v = 0.d0
	omega_star = 0.d0

	out(1) = 0.d0
	out(2) = lue
	out(3) = vcmax_unitiabs
	out(4) = omega_v
	out(5) = omega_star
	  
	lue_vcmax_none(:) = out(:)

END FUNCTION lue_vcmax_none



! LUE for c4
! out_optchi
! kphio
! c_molmass
! soilmstress
FUNCTION lue_vcmax_c4( kphio, c_molmass, soilmstress )

	REAL (kind = 8), INTENT(IN) :: kphio, c_molmass, soilmstress
	REAL (kind=8) :: lue, vcmax_unitiabs, omega_v, omega_star
	REAL (kind = 8), DIMENSION(5) :: out
	REAL (kind = 8) :: lue_vcmax_c4(5)

	! Light use efficiency (gpp per unit absorbed light)
	lue = kphio * c_molmass * soilmstress
	  
	! calculate Vcmax per unit aborbed light
	vcmax_unitiabs  = kphio * soilmstress   ! Eq. 19

	omega_v = 0.d0
	omega_star = 0.d0
	  
	out(2) = lue
	out(3) = vcmax_unitiabs
	out(4) = omega_v
	out(5) = omega_star
	  
	lue_vcmax_c4(:) = out(:)

END FUNCTION lue_vcmax_c4


!-----------------------------------------------------------------------------------------
!                               P hydro-Model auxiliary functions
!-----------------------------------------------------------------------------------------

!Hydraulic vulnerability curve
!
! Calculates the conductivity of the hydraulic pathway (stem or leaf) as a function 
! of water potential, using a 2-parameter Weibull function to describe the hydraulic 
! vulnerability curve:
! P(psi) = (1/2)^((psi/psi_50)^b)
!
! @param psi numeric, Water potential (Mpa)
! @param psi50 numeric, Water potential at which 50\% conductivity is lost (Mpa)
! @param b numeric, Slope of the vulnerability curve
!
! @return Fractional conductivity

REAL (kind = 8) FUNCTION P(psi, psi50, b)
	real(kind = 8), intent(in) :: psi, psi50, b
	P = 0.5d0**((psi/psi50)**b)
end function P

REAL (kind = 8) FUNCTION Pprime(psi, psi50, b)
	real(kind = 8), intent(in) :: psi, psi50, b
	Pprime = log(0.5d0)*(0.5d0**((psi/psi50)**b))*b*(psi/psi50)**(b-1.d0)/psi50
end function Pprime

REAL (kind = 8) FUNCTION Pprimeprime(psi, psi50, b)
	real(kind = 8), intent(in) :: psi, psi50, b
	Pprimeprime = log(0.5d0)*b*(psi/psi50)**(b-1.0)/psi50*(log(0.5d0)*(0.5d0**((psi/psi50)**b))*b*(psi/psi50)**(b-1.d0)/psi50) + &
				  log(0.5d0)*(0.5d0**((psi/psi50)**b))/psi50**2.d0*b*(b-1.d0)*(psi/psi50)**(b-2.d0)
end function Pprimeprime


! Integral of the vulnerabiity curve
! integral of P(psi)dpsi  (Mpa)
REAL (kind = 8) FUNCTION integral_P_ana(dpsi, psi_soil, psi50, b)

	real(kind = 8), intent(in) :: dpsi, psi_soil, psi50, b
	
	real(kind=8) :: ps, pl, l2
	
	ps = psi_soil/psi50
	pl = (psi_soil-dpsi)/psi50
	l2 = log(2.d0) 

	integral_P_ana = -(psi50/b)*(l2**(-1.d0/b))*(gammainc(l2*pl**b, 1.d0/b) - gammainc(l2*ps**b, 1.d0/b))
  
end function integral_P_ana

REAL (kind = 8) FUNCTION integral_P_approx (dpsi, psi_soil, psi50, b)

	real(kind = 8), intent(in) :: dpsi, psi_soil, psi50, b
	
	integral_P_approx = -dpsi * 0.5d0**(( (psi_soil-dpsi/2) /psi50)**b) 
	
end function integral_P_approx


! Returns conductivity in mol/m2/s/Mpa
REAL (kind = 8) FUNCTION scale_conductivity(conductivity, viscosity_water, density_water)

	real(kind = 8), intent(in) :: conductivity, viscosity_water, density_water
	real(kind = 8) :: mol_h20_per_kg_h20, K2, K3, K4
	
	mol_h20_per_kg_h20 = 55.5d0
	
	! Flow rate in m3/m2/s/Pa
	K2 = conductivity / viscosity_water
	K3 = K2 * density_water * mol_h20_per_kg_h20
	K4 = K3*1e6

	scale_conductivity = K4

end function scale_conductivity


REAL (kind = 8) FUNCTION gammainc(x, a)
    real(kind = 8), intent(in) :: x, a
    real(kind = 8) :: xam, gim, ga, s, r, gin, t0
    integer :: k

	xam = -x + a*log(x)

	! Computation of the (upper) incomplete gamma function
	gim = 0.d0
	  
	if (x == 0.d0) then
		ga = gamma(a)
		gim = ga
	else if (x <= 1.d0 + a) then
		s = 1.d0/a
		r = s
		do k=1, 60  
			r = r * x/(a+k)
			s = s+r
			if (abs(r/s) < 1e-15) exit
		end do
		gin = exp(xam) * s
		ga = gamma(a)
		gim = ga - gin
	else if (x > 1.0 + a) then
		t0 = 0.d0
		do k = 60, 1, -1 
		  t0 = (k-a)/(1.d0 + k/(x+t0))
		end do
		gim = exp(xam)/(x+t0)
	end if
  
	gammainc = gim

end function gammainc



REAL (kind = 8) FUNCTION uniroot(func, a, b, tol) 
	
	real(kind=8), intent(in) :: tol 
	real(kind=8) :: x, dx, fa, fb, fx,  a, b
	real(kind=8), external :: func
	
	integer :: max_iter, iter

	max_iter = 1000
	iter = 0
	dx = abs(b - a)
	x = (a + b) / 2.d0

	fa = func(a)
	fb = func(b)
	fx = func(x)

	do while (dx > tol .and. iter < max_iter)
		if (sign(1.d0, fa) /= sign(1.d0, fx)) then
			b = x
			fb = fx
		else
			a = x
			fa = fx
		end if
		dx = dx / 2.d0
		x = (a + b) / 2.d0
		fx = func(x)
		iter = iter + 1
	end do

	uniroot = x
	
end function uniroot



! Stomatal conductance and Transpiration
!
! Calculates regulated stomatal conducatnce given the leaf water potential, 
! plant hydraulic traits, and the environment.
!
! @param dpsi numeric, soil-to-leaf water potential difference (\eqn{\psi_s-\psi_l}), Mpa
! @param psi_soil numeric, soil water potential, Mpa
! @param par_plant numeric, A list of plant hydraulic parameters, which must include ...... 
! @param par_env numeric, A list of environmental parameters
!
! @return numeric, stomatal conductance in mol/m2/s

REAL (kind = 8) FUNCTION calc_gs(dpsi, psi_soil, par_plant, vpd, patm, viscosity_water, density_water)

	real(kind=8), intent(in) :: dpsi, psi_soil, vpd, patm, viscosity_water, density_water
	real(kind=8), dimension(3), intent(in) :: par_plant
	real(kind=8) :: K, D

	K = scale_conductivity(par_plant(1), viscosity_water, density_water)
	D = vpd/patm
	calc_gs = K/1.6d0/D * (- 1.d0) * integral_P_ana(dpsi, psi_soil, par_plant(2), par_plant(3))

end function calc_gs  

REAL (kind = 8) FUNCTION calc_gsprime(dpsi, psi_soil, par_plant, vpd, patm, viscosity_water, density_water)	

	real(kind=8), intent(in) :: dpsi, psi_soil, vpd, patm, viscosity_water, density_water
	real(kind=8), dimension(3), intent(in) :: par_plant
	real(kind=8) :: K, D

	K = scale_conductivity(par_plant(1), viscosity_water, density_water)
	D = vpd/patm
	calc_gsprime = K/1.6d0/D * P(psi_soil-dpsi, par_plant(2), par_plant(3))

end function calc_gsprime  


REAL (kind = 8) FUNCTION calc_jmax_from_J(J, phi0, iabs)
	
	real(kind=8), intent(in) :: J, phi0, iabs
	real(kind=8) :: p

   	!p = phi0 * iabs
    !calc_jmax_from_J=4.d0*p/((4.d0*p/J)**2.d0-1.d0)**(1.d0/2.d0)
	!calc_jmax_from_J=4.d0*p/((p/J)**2.d0-1.d0)**(1.d0/2.d0)
	p = 4.d0*phi0 * iabs
	calc_jmax_from_J=p/((p/J)**2.d0-1.d0)**(1.d0/2.d0)
	
end function calc_jmax_from_J

REAL (kind = 8) FUNCTION calc_djmax_dJ(J, phi0, iabs)

	real(kind=8), intent(in) :: J, phi0, iabs
	real(kind=8) :: p, sq, psq
	
	!p = phi0 * iabs
	!calc_djmax_dJ = (4.d0*p)**3.d0/((4.d0*p)**2.d0-J**2.d0)**(3.d0/2.d0)
	!calc_djmax_dJ = 4.d0*p**3.d0/J**3.d0/((p/J)**2.d0-1.d0)**(3.d0/2.d0)
	p = 4.d0 * phi0 * iabs
	sq = (p**2.d0-J**2.d0)**(1.d0/2.d0)
	calc_djmax_dJ = (p/sq)**3.d0
	
end function calc_djmax_dJ

REAL (kind = 8) FUNCTION calc_dJ_dchi(gs, x, gammastar, ca_v, kmm, patm, delta)

	real(kind=8), intent(in) :: gs, x, gammastar, ca_v, kmm, patm, delta
	real(kind=8) :: g, k, ca, d

	g = gammastar/ca_v
	k = kmm/ca_v
	ca = ca_v/patm*1e6
	d = delta

	calc_dJ_dchi = 4.d0*gs*ca * ((d*(2*g*(k + 1.d0) + k*(2.d0*x - 1.d0) + x**2.d0) - &
				((x-g)**2.d0+3.d0*g*(1.d0-g)))/(d*(k + x) + g - x)**2.d0)
	
	!calc_dJ_dchi = gs*ca * ((d*(2.d0*g*(k + 1.d0) + k*(2.d0*x - 1.d0) + x**2.d0) - &
	!			((x-g)**2.d0+3.d0*g*(1.d0-g)))/(d*(k + x) + g - x)**2.d0)


end function calc_dJ_dchi


REAL (kind = 8) FUNCTION calc_x_from_dpsi(dpsi, psi_soil, par_plant, gammastar, patm, vpd, viscosity_water, &
										density_water, kmm, ca_v, delta, par_cost)

	real(kind=8), intent(in) :: dpsi, psi_soil, gammastar, patm, kmm, ca_v, delta, vpd, viscosity_water, density_water
	real(kind=8), dimension(3), intent(in) :: par_plant
	real(kind=8), dimension(2), intent(in) :: par_cost
	
	real(kind=8) :: gstar, Km, ca, br, y, K, D, gsprime, x

	gstar = gammastar/patm*1e6
	Km = kmm/patm*1e6
	ca = ca_v/patm*1e6
	br = delta
	y = par_cost(2)
	
	K = scale_conductivity(par_plant(1), viscosity_water, density_water)
    D = (vpd/patm)
	gsprime = K/1.6d0/D*P(psi_soil-dpsi, par_plant(2), par_plant(3))
	
	x = (-2.d0*ca*dpsi*(gstar + br*Km)*y + &
		ca**2.d0*((3.d0 - 2.d0*br)*gstar + br*Km)*gsprime +  &
		-sqrt(2.d0)*sqrt( ca**2.d0*dpsi*((-3.d0 + 2.d0*br)*gstar - br*Km)*((-1.d0 + br)*ca + gstar + & 
		br*Km)*y*(-2.d0*dpsi*y + (ca + 2.d0*gstar)*gsprime)))/ &
		(ca**2.d0*(2.d0*(-1.d0 + br)*dpsi*y + ((3.d0 - 2.d0*br)*gstar + br*Km)*gsprime)) 
		
	if (x<((gstar + br*Km)/(ca - br*ca))) then
		x=(gstar + br*Km)/(ca - br*ca)+1e-12
	end if
	
	calc_x_from_dpsi = x
	
end function calc_x_from_dpsi  


REAL (kind = 8) FUNCTION calc_J(gs, x, gammastar, ca_v, kmm, delta, patm)

	real(kind=8), intent(in) :: gs, x, gammastar, ca_v, kmm, delta, patm
	real(kind=8) :: g, k, ca, d

	g = gammastar/ca_v
	k = kmm/ca_v
	ca = ca_v/patm*1e6
	d = delta
	calc_J=4.d0*gs*ca*(1.d0-x)*(x+2.d0*g)/(x*(1.d0-d)-(g+d*k)) !- documentation version is different from the package
	!calc_J=gs*ca*(1.d0-x)*(x+2.d0*g)/(x*(1.d0-d)-(g+d*k))	
  
end function calc_J


REAL (kind = 8) FUNCTION f1(dpsi, psi_soil, viscosity_water, density_water, par_plant, vpd, &
							patm, gammastar, kmm, ca_v, delta, par_cost, phi0, iabs)

	real(kind=8), intent(in) :: dpsi, psi_soil, viscosity_water, density_water, vpd, patm, &
								gammastar, kmm, ca_v, delta,  phi0, iabs	
	real(kind=8), dimension(3), intent(in) :: par_plant
	real(kind=8), dimension(2), intent(in) :: par_cost
	
	real(kind=8) :: gs, x, ajmax

    gs = calc_gs(dpsi, psi_soil, par_plant, vpd, patm, viscosity_water, density_water)
    x = calc_x_from_dpsi(dpsi, psi_soil, par_plant, gammastar, patm, vpd, viscosity_water, &
										density_water, kmm, ca_v, delta, par_cost)
    ajmax=calc_J(gs, x, gammastar, ca_v, kmm, delta, patm)-4.d0*phi0*iabs
	
    f1=ajmax
	
end function f1


REAL (kind = 8) FUNCTION f2(dpsi,psi_soil,par_plant,viscosity_water, density_water, ca_v, gammastar, patm, vpd, par_cost)

	real(kind=8), intent(in) :: dpsi,psi_soil,viscosity_water, density_water, ca_v, gammastar, patm, vpd
	real(kind=8), dimension(2), intent(in) :: par_cost
	real(kind=8), dimension(3), intent(in) :: par_plant
	real(kind=8) :: K, D, gstar, y, gsprime

	K = scale_conductivity(par_plant(1), viscosity_water, density_water)
    D = vpd/patm
	gstar = gammastar/patm*1e6
	y = par_cost(2)
	
	gsprime = K/1.6d0/D*P(psi_soil-dpsi, par_plant(2), par_plant(3))
 
	f2 = -2.d0*dpsi*y + (ca_v + 2.d0*gstar)* gsprime
	
end function f2


REAL (kind = 8) FUNCTION dFdx(dpsi, psi_soil, gammastar, ca_v, kmm, delta, patm, vpd, &
						phi0, iabs, par_plant,viscosity_water, density_water, par_cost)

	real(kind=8), intent(in) :: dpsi, psi_soil, gammastar, ca_v, kmm, delta, patm, vpd, phi0, iabs, &
								viscosity_water, density_water
	real(kind=8), dimension(3), intent(in) :: par_plant
	real(kind=8), dimension(2), intent(in) :: par_cost
	
	real(kind=8) :: gs, gsprime, x, J, ca, g, djmax_dJ, dJ_dchi

	gs = calc_gs(dpsi, psi_soil, par_plant, vpd, patm, viscosity_water, density_water)
	gsprime = calc_gsprime(dpsi, psi_soil, par_plant, vpd, patm, viscosity_water, density_water)
	  
	x =  calc_x_from_dpsi(dpsi, psi_soil, par_plant, gammastar, patm, vpd, viscosity_water, &
										density_water, kmm, ca_v, delta, par_cost)
	  
	J = calc_J(gs, x, gammastar, ca_v, kmm, delta, patm)
	  
	ca = ca_v/patm*1e6
	g = gammastar/ca
	  
	djmax_dJ = calc_djmax_dJ(J, phi0, iabs)
	dJ_dchi = calc_dJ_dchi(gs, x, gammastar, ca_v, kmm, patm, delta)
	  
	dFdx = -gs*ca - par_cost(1) * djmax_dJ * dJ_dchi
  
end function dFdx

! Analytical chi in the case of strong Jmax limitation
REAL (kind = 8) FUNCTION chi_jmax_limited(gammastar, kmm, ca, delta, par_cost)

	real(kind = 8), intent(in) :: gammastar, kmm, ca, delta
	real(kind = 8), dimension(2) :: par_cost
	
	real(kind=8) :: g, k, b, a
	
	g = gammastar/ca
	k = kmm/ca
	b = delta
	a = par_cost(1)
	  
	chi_jmax_limited = (2.d0*sqrt(-a*(4.d0*a + b - 1.d0)*(-3.d0*g + 2.d0*b*g - b*k)* &
		(-1.d0 + b + g + b*k)) - (4.d0*a + b - 1.d0)*(b*k + g))/((b - 1.d0)*(4.d0*a + b - 1.d0))
  
end function chi_jmax_limited


! Canopy transpiration
REAL (kind = 8) FUNCTION transpiration(dpsi, psi_soil, par_plant, vpd, patm,viscosity_water, density_water)

	real(kind = 8), intent(in) :: dpsi, psi_soil, vpd, patm, viscosity_water, density_water
	real(kind = 8), dimension(3) :: par_plant
	
	real(kind=8) :: k, d
	
	k = scale_conductivity(par_plant(1), viscosity_water, density_water)
	d = vpd/patm
	transpiration = -k * integral_P_ana(dpsi, psi_soil, par_plant(2), par_plant(3))
  
end function transpiration

REAL (kind = 8) FUNCTION calc_dpsi_bound(psi_soil, ca_v, gammastar, patm, vpd, par_plant, viscosity_water, &
										 density_water, par_cost,kmm, phi0, iabs, delta)

	real(kind = 8), intent(in) :: psi_soil, ca_v, gammastar, patm, vpd, viscosity_water, density_water, &
									kmm, phi0, iabs, delta										 
	real(kind = 8), dimension(3), intent(in) ::  par_plant
	real(kind = 8), dimension(2), intent(in) ::  par_cost

	real(kind=8) :: gstar, K, y, Pox, Ppox, Pppox, iabs_bound, ca, a, b, c, del, approx_O2, exact, use_bound, &
					xzero, fzero
	integer iflag


	gstar = gammastar/patm*1e6
	
	ca = ca_v/patm*1e6
	y = par_cost(2)
	K = scale_conductivity(par_plant(1), viscosity_water, density_water)
	K = K/(1.6d0*vpd/patm)
	Pox = P(psi_soil, par_plant(2), par_plant(3))
	Ppox = Pprime(psi_soil, par_plant(2), par_plant(3))
	Pppox = Pprimeprime(psi_soil, par_plant(2), par_plant(3))
   
	a = (ca + 2*gstar)*K*Pppox*4.d0/8.d0
	b = -(2.d0*y + (ca + 2.d0*gstar)*K*Ppox)
	c = (ca + 2.d0*gstar)*K*Pox
	del = b**2.d0-4.d0*a*c

	approx_O2 = (-b-sqrt(del))/2.d0/a
	
	call root_scalar('zhang',f_f2,0.d0,10.d0,xzero, fzero,iflag)
	
	!exact = xzero

	use_bound = xzero
	
	call root_scalar('zhang',f_f1,use_bound*0.001d0, use_bound*0.99d0, xzero, fzero,iflag)
	
	calc_dpsi_bound = xzero 
	
	contains
	
	REAL (kind = 8) FUNCTION f_f2 (dpsi)
	
		real(kind = 8), intent(in) :: dpsi

		f_f2 = f2(dpsi,psi_soil,par_plant,viscosity_water, density_water, ca, gammastar, patm, vpd, par_cost)
	
	end function f_f2
	
	
	REAL (kind = 8) FUNCTION f_f1 (dpsi)
	
		real(kind = 8), intent(in) :: dpsi

		f_f1 = f1(dpsi, psi_soil, viscosity_water, density_water, par_plant, vpd, &
							patm, gammastar, kmm, ca_v, delta, par_cost, phi0, iabs)
	
	end function f_f1
	
	  
END FUNCTION calc_dpsi_bound

end module gpp_subroutines