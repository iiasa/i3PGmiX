module mortality

IMPLICIT NONE
public :: f_get_mortality, mortality_brandl, mortality_lpj, mortality_iland, remove_mortality

contains

    function f_get_mortality(stems_n, WS, mS, wSx1000, thinPower) result(mort_n)
        ! calculate the mortality

        !input
        real(kind=8), intent(in) :: stems_n, WS, mS, wSx1000, thinPower

        ! output
        real(kind=8) :: mort_n

        ! local
        real(kind=8), parameter :: accuracy = 1.d0 / 1000.d0
        integer :: i
        real(kind=8) :: fN,dfN,dN,n,x1,x2


        n = stems_n / 1000.d0
        x1 = 1000.d0 * mS * WS / stems_n
        i = 0

        do
            i = i + 1

            if (n <= 0.d0) exit !added in 3PG+

            x2 = wSx1000 * n ** (1.d0 - thinPower)
            fN = x2 - x1 * n - (1.d0 - mS) * WS
            dfN = (1.d0 - thinPower) * x2 / n - x1
            dN = -fN / dfN
            n = n + dN

            if (abs(dN) <= accuracy .Or. i >= 5) exit

        end do

        mort_n = stems_n - 1000.d0 * n

    end function f_get_mortality

    ! Empirical survival probability curves from Brandl et al. (2020) fitted to ICP data.
    ! Mortality is given by the difference in survival probability from one year to the next,
    ! based on species and climate.
    function mortality_brandl(species_group, age, tmax, tmin, tmean, pr, month) result(mort_n)

        !input
        integer, intent(in):: month
        real(kind=8), intent(in) :: species_group, age
        real(kind=8), dimension(12), intent(in) :: tmax, tmin, tmean, pr

        !local variables
        real(kind=8):: beta_repar   

        ! output
        real(kind=8):: mort_n

        ! Beech
        if (species_group == 1.0) then
            beta_repar = 3494.69016d0 / ( Exp(0.108d0 * tmax(month)) * Exp(-0.05d0 * tmin(month-6)) )
            !mort_n = exp(-(age / beta_repar )**1.10627d0) - exp(-( (age+1) / beta_repar )**1.10627d0)  
            mort_n = get_P_survival(1.10627d0, beta_repar, age ) - get_P_survival(1.10627d0, beta_repar, age + 1.d0 )
            
        end if

        ! Douglas Fir
        if (species_group == 2.0) then
            beta_repar = 843.8713d0 / ( Exp(0.188d0 * 1.d0/12.d0 * sum(tmean(month-6:month+5))) * &
                Exp(- 0.001 * 1/12 * sum(pr(month-6:month+5))) * Exp(0.671d0 * 0.745d0) )
            !mort_n = exp(( age / beta_repar )**1.05654d0) - exp(( (age+1) / beta_repar )**1.05654d0) 
            mort_n = get_P_survival(1.05654d0, beta_repar, age ) - get_P_survival(1.05654d0, beta_repar, age + 1.d0 )
        end if

        ! Fir
        if (species_group == 3.0) then
            beta_repar = 438.78081d0 / Exp(0.188d0 * 1.d0/12.d0 * sum(tmean(month-6:month+5)))
            !mort_n = exp(( age / beta_repar )**1.05654d0) - exp(( (age+1) / beta_repar )**1.05654d0) 
            mort_n = get_P_survival(1.05654d0, beta_repar, age ) - get_P_survival(1.05654d0, beta_repar, age + 1.d0 )
        end if

        ! Oak
        if (species_group == 4.0) then
            beta_repar = 2177.64653d0 / Exp(0.095d0 * tmax(month))
            !mort_n = exp(( age / beta_repar )**1.05865d0) - exp(( (age+1) / beta_repar )**1.05865d0) 
            mort_n = get_P_survival(1.05865d0, beta_repar, age ) - get_P_survival(1.05865d0, beta_repar, age + 1.d0 )
        end if

        ! Pine
        if (species_group == 5.0) then
            beta_repar = 681.97977d0 / ( Exp(0.077d0 * tmax(month)) * &
                    Exp(-0.002d0 * sum(pr(month-1:month+1)) ) )
            !mort_n = exp(( age / beta_repar )**0.55878d0) - exp(( (age+1) / beta_repar )**0.55878d0) 
            mort_n = get_P_survival(0.55878d0, beta_repar, age ) - get_P_survival(0.55878d0, beta_repar, age + 1.d0 )
        end if

        ! Spruce
        if (species_group == 6.0) then
            beta_repar = 607.28609d0 / ( Exp(0.058d0 * tmax(month)) * &
                    Exp(-0.0002d0 * sum(pr(month-1:month+1)) ) * Exp(0.653d0 * 0.801d0) )
            !mort_n = exp(( age / beta_repar )**1.29305d0) - exp(( (age+1) / beta_repar )**1.29305d0) 
            mort_n = get_P_survival(1.29305d0, beta_repar, age ) - get_P_survival(1.29305d0, beta_repar, age + 1.d0 )
        end if

    end function mortality_brandl

    !LPJ-GUESS mortality implementation (Smith et al. 2014).
    function mortality_lpj(greff_sp, npp, lai, dmC, nyears, npp_corr) result(mort_n)

        !input
        integer, intent(in):: nyears
        real(kind=8), intent(in) :: greff_sp, dmC, npp_corr
        real(kind=8), dimension(nyears*12), intent(in) :: npp, lai

        !output
        real(kind=8) :: mort_n

        ! local
        integer :: i, n, cnt
        real(kind=8) :: greff, lai_avg

        ! Numbwer of months
        n = nyears * 12

        ! Get average LAI, excluding the dormant period
        cnt = 0
        lai_avg = 0.d0
        do i=1, n
            if (lai(i) > 0.d0) then
                lai_avg = lai_avg + lai(i)
                cnt = cnt + 1
            end if
        end do

        ! Average LAI over the period
        lai_avg = lai_avg / cnt

        ! Average growth efficiency over the period
        greff = npp_corr * 1000.d0 * (sum(npp) * dmC / (10.d0 * nyears)) / lai_avg !(sum(lai)/n)  ! ! convert to gC/m²/year

        ! Mortality share
        ! European version
        mort_n = 0.1d0 / (1 + (greff/greff_sp)**5.d0)

        !Default version
        !mort_n = 0.05d0 / (1 + 35.d0 * greff)

    end function mortality_lpj

    ! Growth efficiency parameter
    function calc_greff(npp, lai, dmC, nyears, npp_corr) result(greff)

        !input
        integer, intent(in):: nyears
        real(kind=8), intent(in) :: dmC, npp_corr
        real(kind=8), dimension(nyears*12), intent(in) :: npp, lai

        !output
        real(kind=8) :: mort_n

        ! local
        integer :: i, n, cnt
        real(kind=8) :: greff, lai_avg

        ! Numbwer of months
        n = nyears * 12

        ! Get average LAI, excluding the dormant period
        cnt = 0
        lai_avg = 0.d0
        do i=1, n
            if (lai(i) > 0.d0) then
                lai_avg = lai_avg + lai(i)
                cnt = cnt + 1
            end if
        end do

        ! Average LAI over the period
        lai_avg = lai_avg / cnt

        ! Average growth efficiency over the period
        greff = npp_corr * 1000.d0 * (sum(npp) * dmC / (10.d0 * nyears)) / lai_avg !(sum(lai)/n)  ! ! convert to gC/m²/year

    end function calc_greff

    ! iLand mortality implementation
    function mortality_iland(npp, wf_loss, wr_loss, bs) result(mort_n)

        !input
        real(kind=8), intent(in):: bs
        real(kind=8), dimension(12), intent(in) :: npp, wf_loss, wr_loss

        !output
        real(kind=8) :: mort_n

        ! local
        real(kind=8) :: si

        ! get stress index
        si = max(1.d0 - sum(npp) / (sum(wf_loss) + sum(wr_loss)), 0.d0)

        ! get mortality
        mort_n = 1.d0 - exp(-bs * si)

    end function mortality_iland

    ! Update stand after mortality routine
    subroutine remove_mortality(mort_stress, stems_n, biom_foliage, biom_root, biom_stem, dead_volume, &
                                fracBB, wood_density, mF, mR, mS)

        !input
        real(kind=8), intent(in):: fracBB, wood_density, mF, mR, mS

        !output
        real(kind=8), intent(inout) :: mort_stress, stems_n, biom_foliage, biom_root, &
                                        biom_stem, dead_volume

        mort_stress = mort_stress * stems_n

        mort_stress = min( mort_stress, stems_n) ! Mortality can't be more than available

        biom_foliage = biom_foliage - mF * mort_stress * (biom_foliage / stems_n)
        biom_root = biom_root - mR * mort_stress * (biom_root / stems_n)
        biom_stem = biom_stem - mS * mort_stress * (biom_stem / stems_n)
        stems_n = stems_n - mort_stress
        dead_volume = mS * mort_stress * (biom_stem / stems_n) * &
                        (1.d0 - fracBB) / wood_density

    end subroutine remove_mortality

    ! Get survival probability based on Weibull function
    function get_P_survival(alpha, beta, age) result(mort_prob)
        ! calculate the mortality

        !input
        real(kind=8), intent(in) :: alpha, beta, age

        ! output
        real(kind=8) :: mort_prob

        mort_prob = exp(-(age/beta)**alpha)

    end function get_P_survival


end module mortality