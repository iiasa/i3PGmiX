module ra_calc

use utils
use soil_cnp_subroutines
use mod_decl_const

! Autotrophic respiration module
! Computes maintenance respiration based on the fine root and sapwood biomass

IMPLICIT NONE
public :: get_raut

contains

	! autotrophic respiration from sapwood and fine roots
	subroutine get_raut(ra, rsap, pssN0, pssN1, tpssN, npssN, age, dbh, height, basal_area, crown_length, &
						stems_n, f_phys, m, biom_root, biom_foliage, biom_foliage_debt, biom_fineroot, biom_stem, &
						biom_sapwood, fr_ratio, Nfr, Nsw, rd_p, rd_ph, t_soil, tmp_ave, month, gpp_model, n_sp )

		! Input
		integer, intent(in):: n_sp, gpp_model, month
		real(kind=8), dimension(n_sp), intent(in):: pssN0, pssN1, tpssN, npssN, age, dbh, height, stems_n, &
													basal_area, crown_length, f_phys, m, biom_root, &
													biom_foliage, biom_foliage_debt, biom_stem, fr_ratio, &
													Nfr, Nsw, rd_p, rd_ph
		
		real(kind=8), intent(in):: t_soil, tmp_ave

		! Output
		real(kind=8), dimension(n_sp), intent(out):: ra, biom_sapwood, biom_fineroot, rsap 

		! Local
		integer :: i
		real(kind=8), dimension(n_sp):: height_sapwood, biom_sapwood_ratio, fr_wn, bmfr

		! Sapwood share                       
		rsap(:) =  min(pssN1(:) + (pssN0(:) -  pssN1(:)) * Exp(-ln2 * ( age(:) / tpssN(:)) ** npssN(:)) , 1.d0)

		! Sapwood hieght according to 4C implementation
		height_sapwood(:) = ((2.d0/3.d0) * (height(:) - crown_length(:)) + (1.d0/3.d0) * height(:))

		! Share of sapwood/heartwood
		biom_sapwood_ratio(:) = f_heartwood(age(:), dbh(:), basal_area(:), stems_n(:), rsap(:), &
									height(:), crown_length(:), height_sapwood(:), n_sp )  

		! fine root modifier
		fr_wn(:) = 1.3d0 / (1.d0  + 0.3d0 * f_phys(:) * m(:))

		biom_fineroot(:) = biom_foliage(:) * fr_ratio(:) * fr_wn(:)

		do i=1, n_sp
			if (biom_foliage(i) == 0.d0) then
				bmfr(i) = biom_fineroot(i)
			else
				bmfr(i) = biom_foliage_debt(i) * fr_ratio(i) * fr_wn(i)
			end if
		end do

		! Sapwood biomass                                    
		biom_sapwood(:) = biom_stem(:) * biom_sapwood_ratio(:) + (biom_root(:) - bmfr(:))
				
		where (crown_length(:) == 0.d0 .or. height(:) == 0.d0) 
			biom_sapwood(:) = biom_stem(:) * rsap(:) + (biom_root(:) - biom_fineroot(:))
		end where


		!Maintenance respiration - temperature effect with acclimation based on the 3D-CMCC implementation from Collalti et al. (2016)
		ra(:) = (biom_fineroot(:)*Nfr(:)/100.d0) * 0.218d0 * daysInMonth(month) * (1.d0 / dmC) * &
			((3.22d0 - t_soil/21.7391d0)**((t_soil-20.d0)/10.d0)) * (10.d0 ** (-0.00794d0 * (t_soil - 20.d0))) + &  
			(biom_sapwood(:)*Nsw(:)/100.d0) * 0.218d0 * daysInMonth(month) * (1.d0 / dmC) * &
			2.d0**((tmp_ave - 20.d0)/10.d0) * (10.d0 ** (-0.00794d0 * (tmp_ave - 20.d0))) 


		! Add foliage maintenance respiration
		if (gpp_model .eq. int(1)) then
			! In case not given by the GPP computation, use the GOLTILWA+ approach (Gracia & Sabaté, 2003)
			ra(:) = ra(:) + biom_foliage(:) * (55.5d0 / 4700.d0) * daysInMonth(month) * 2.d0 ** ((tmp_ave-25.d0)/10.d0)
		else if (gpp_model .eq. int(2)) then
			ra(:) = ra(:) + (rd_p) / 100.d0 * (1.d0 / dmC) ! Foliage respiration from the P-model
		else if (gpp_model .eq. int(3)) then
			ra(:) = ra(:) + (rd_ph * daysInMonth(month)) / 100.d0 * (1.d0 / dmC) ! Foliage respiration from the Phydro-model
		end if


	end subroutine get_raut


end module ra_calc