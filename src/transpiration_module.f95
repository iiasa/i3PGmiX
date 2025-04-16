module transpiration

use mod_decl_const


IMPLICIT NONE
public :: s_transpiration_3pgpjs, s_transpiration_3pgmix, soil_evap

contains

	subroutine s_transpiration_3pgpjs ( n_sp, solar_rad, day_length, VPD_sp, BLcond, conduct_canopy, days_in_month, Qa, Qb, &
		transp_veg)

	implicit none

	! input
	integer, intent(in) :: n_sp ! number of species
	real(kind=8), intent(in) :: solar_rad
	real(kind=8), intent(in) ::  day_length
	real(kind=8), dimension(n_sp), intent(in) :: VPD_sp
	real(kind=8), dimension(n_sp), intent(in) :: BLcond
	real(kind=8), dimension(n_sp), intent(in) :: conduct_canopy
	integer, intent(in) :: days_in_month
	real(kind=8), dimension(n_sp), intent(in) :: Qa, Qb

	! output
	real(kind=8), dimension(n_sp), intent(out) :: transp_veg

	! derived variables
	real(kind=8), dimension(n_sp) :: netRad
	real(kind=8), dimension(n_sp) :: defTerm
	real(kind=8), dimension(n_sp) :: div


	if( sum(VPD_sp(:)) == 0.d0 ) then
		transp_veg(:) = 0.d0

	else
		netRad(:) = (Qa + Qb * (solar_rad * 10.d0 ** 6.d0 / day_length))
		!netRad(:) = max(netRad(:), 0.d0) ! net radiation can't be negative
		!SolarRad in MJ/m2/day ---> * 10^6 J/m2/day ---> /day_length converts to only daytime period ---> W/m2
		defTerm(:) = rhoAir * lambda * (VPDconv * VPD_sp(:)) * BLcond(:)
		div(:) = conduct_canopy(:) * (1.d0 + e20) + BLcond(:)

		transp_veg(:) = days_in_month * conduct_canopy(:) * (e20 * netRad(:) + defTerm(:)) / div(:) / lambda * day_length
		! in J/m2/s then the "/lambda*h" converts to kg/m2/day and the days in month then coverts this to kg/m2/month
		transp_veg(:) = max(0.d0, transp_veg(:)) ! transpiration can't be negative

	end if

	end subroutine s_transpiration_3pgpjs


	subroutine s_transpiration_3pgmix ( n_sp, solar_rad, vpd_day, day_length, days_in_month, lai, fi, VPD_sp, &
		aero_resist, conduct_canopy, conduct_soil, Qa, Qb, &
		transp_veg, evapotra_soil)

	implicit none

	! input
	integer, intent(in) :: n_sp ! number of species

	real(kind=8), intent(in) :: solar_rad
	real(kind=8), intent(in) :: vpd_day
	real(kind=8), intent(in) ::  day_length
	integer, intent(in) :: days_in_month
	real(kind=8), dimension(n_sp), intent(in) :: lai
	real(kind=8), dimension(n_sp), intent(in) :: fi
	real(kind=8), dimension(n_sp), intent(in) :: VPD_sp
	real(kind=8), dimension(n_sp), intent(in) :: aero_resist
	real(kind=8), dimension(n_sp), intent(in) :: conduct_canopy
	real(kind=8), intent(in) ::  conduct_soil
	real(kind=8), dimension(n_sp), intent(in) :: Qa, Qb

	! output
	real(kind=8), dimension(n_sp), intent(out) :: transp_veg
	real(kind=8), intent(out) :: evapotra_soil

	! derived variables
	real(kind=8), dimension(n_sp) :: netRad
	real(kind=8), dimension(n_sp) :: defTerm
	real(kind=8), dimension(n_sp) :: div
	real(kind=8) :: lai_total ! here is is a number, while in the main subroutine it is a vector
	real(kind=8) :: netRad_so
	real(kind=8) :: defTerm_so
	real(kind=8) :: div_so ! ending `so` mean soil

	! Species level calculations ---
	! the within canopy aero_resist and VPDspecies have been calculated using information from the light submodel
	! and from the calculation of the modifiers. The netrad for each species is calculated
	! using the fi (proportion of PAR absorbed by the given species) and is calculated by the light submodel.

	if( sum(lai(:)) == 0.d0 ) then
		transp_veg(:) = 0.d0
	else
		netRad(:) = (Qa + Qb * (solar_rad * 10.d0 ** 6.d0 / day_length)) * fi(:)
		!netRad(:) = max(netRad(:), 0.d0) ! net radiation can't be negative
		!SolarRad in MJ/m2/day ---> * 10^6 J/m2/day ---> /day_length converts to only daytime period ---> W/m2
		defTerm(:) = rhoAir * lambda * (VPDconv * VPD_sp(:)) / aero_resist(:)
		div(:) = conduct_canopy(:) * (1.d0 + e20) + 1.d0 / aero_resist(:)

		transp_veg(:) = days_in_month * conduct_canopy(:) * (e20 * netRad(:) + defTerm(:)) / div(:) / lambda * day_length
		! in J/m2/s then the "/lambda*h" converts to kg/m2/day and the days in month then coverts this to kg/m2/month

		transp_veg(:) = max(transp_veg(:), 0.d0)
		
		where( lai(:) == 0.d0 )
			transp_veg(:) = 0.d0
		end where

	end if

	! now get the soil evaporation (soil aero_resist = 5 * lai_total, and VPD of soil = VPD * Exp(lai_total * -Log(2) / 5))
	lai_total = sum( LAI(:) )

	if( lai_total > 0 ) then
		defTerm_so = rhoAir * lambda * (VPDconv * (vpd_day * Exp(lai_total * (-ln2) / 5.d0))) / (5.d0 * lai_total)
		div_so = conduct_soil * (1.d0 + e20) + 1.d0 / (5.d0 * lai_total)
	else
		!defTerm_so = 0.d0
		defTerm_so = rhoAir * lambda * (VPDconv * (vpd_day * Exp(lai_total * (-ln2) / 5.d0)))
		div_so = conduct_soil * (1.d0 + e20) + 1.d0
	end if

	netRad_so = (Qa(1) + Qb(1) * (solar_rad * 10.d0 ** 6.d0 / day_length)) * (1.d0 - sum( fi(:) ) )
	!SolarRad in MJ/m2/day ---> * 10^6 J/m2/day ---> /day_length converts to only daytime period ---> W/m2

	evapotra_soil = days_in_month * conduct_soil * (e20 * netRad_so + defTerm_so) / div_so / lambda * day_length
	!in J/m2/s then the "/lambda*h" converts to kg/m2/day and the days in month then coverts this to kg/m2/month

	end subroutine s_transpiration_3pgmix

	! Soil evaporation
	subroutine soil_evap( solar_rad, vpd_day, day_length, days_in_month,lai_total,  fi, &
			conduct_soil, Qa, Qb, evapotra_soil)

	implicit none

	! input

	real(kind=8), intent(in) :: solar_rad
	real(kind=8), intent(in) :: lai_total
	real(kind=8), intent(in) :: vpd_day
	real(kind=8), intent(in) ::  day_length
	integer, intent(in) :: days_in_month
	real(kind=8), intent(in) :: fi
	real(kind=8), intent(in) ::  conduct_soil
	real(kind=8), intent(in) :: Qa, Qb

	! output
	real(kind=8), intent(out) :: evapotra_soil

	! derived variables

	real(kind=8) :: netRad_so
	real(kind=8) :: defTerm_so
	real(kind=8) :: div_so ! ending `so` mean soil


	if( lai_total > 0 ) then
		defTerm_so = rhoAir * lambda * (VPDconv * (vpd_day * Exp(lai_total * (-ln2) / 5.d0))) / (5.d0 * lai_total)
		div_so = conduct_soil * (1.d0 + e20) + 1.d0 / (5.d0 * lai_total)
	else
		!defTerm_so = 0.d0
		defTerm_so = rhoAir * lambda * (VPDconv * (vpd_day * Exp(lai_total * (-ln2) / 5.d0)))
		div_so = conduct_soil * (1.d0 + e20) + 1.d0
	end if

	netRad_so = (Qa + Qb * (solar_rad * 10.d0 ** 6.d0 / day_length)) * (1.d0 - fi )
	!SolarRad in MJ/m2/day ---> * 10^6 J/m2/day ---> /day_length converts to only daytime period ---> W/m2

	evapotra_soil = days_in_month * conduct_soil * (e20 * netRad_so + defTerm_so) / div_so / lambda * day_length
	!in J/m2/s then the "/lambda*h" converts to kg/m2/day and the days in month then coverts this to kg/m2/month

	end subroutine soil_evap

end module transpiration