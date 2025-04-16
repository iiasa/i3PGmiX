module utils

use mod_decl_const

IMPLICIT NONE
public :: f_exp, f_exp_foliage, f_gamma_dist, f_orderId, p_min_max, f_rsap, f_elai, f_heartwood, f_greff

contains

    function f_exp(n_m, x, g0, gx, tg, ng) result( out )

        implicit none

        ! input
        integer, intent(in) :: n_m
        real(kind=8), dimension(n_m), intent(in) :: x
        real(kind=8), intent(in) :: g0, gx, tg, ng

        ! output
        real(kind=8), dimension(n_m) :: out

        out(:) = gx

        if ( tg /= 0.d0 ) then
            out(:) = gx + (g0 - gx) * Exp(-ln2 * ( x(:) / tg) ** ng)
        end if

    end function f_exp


    function f_exp_foliage(n_m, x, f1, f0, tg) result( out )

        implicit none

        ! input
        integer, intent(in) :: n_m
        real(kind=8), dimension(n_m), intent(in) :: x
        real(kind=8), intent(in) :: f1, f0, tg

        ! output
        real(kind=8), dimension(n_m) :: out

        ! local
        real(kind=8) :: kg

        if( tg * f1 == 0.d0 ) then
            out(:) = f1
        else
            kg = 12.d0 * Log(1.d0 + f1 / f0) / tg
            out(:) = f1 * f0 / (f0 + (f1 - f0) * Exp(-kg * x))
        end if

    end function f_exp_foliage


    function f_gamma_dist( x, n ) result( out )

        implicit none

        ! input
        integer, intent(in) :: n
        real(kind=8), dimension(n), intent(in) :: x

        ! output
        real(kind=8), dimension(n) :: out

        out = x ** (x - 0.5d0) * 2.718282d0 ** (-x) * (2.d0 * Pi) ** (0.5d0) * &
            (1.d0 + 1.d0 / (12.d0 * x) + 1.d0 / (288.d0 * x ** 2.d0) - 139.d0 / (51840.d0 * x ** 3.d0) - &
            571.d0 / (2488320.d0 * x ** 4.d0))

    end function f_gamma_dist


    function f_orderId(x) result(id)
        ! Returns the indices that would sort an array.

        implicit none

        ! input
        real(kind=8), intent(in) :: x(:)       ! array of numbers

        ! output
        integer :: id( size(x) )            ! indices into the array 'x' that sort it

        ! local
        integer :: i, n, imin, temp1        ! helpers
        real(kind=8) :: temp2
        real(kind=8) :: x2( size(x) )

        x2 = x
        n = size(x)

        do i = 1, n
            id(i) = i
        end do

        do i = 1, n-1
            ! find ith smallest in 'a'
            imin = minloc(x2(i:),1) + i - 1
            ! swap to position i in 'a' and 'b', if not already there
            if (imin /= i) then
                temp2 = x2(i); x2(i) = x2(imin); x2(imin) = temp2
                temp1 = id(i); id(i) = id(imin); id(imin) = temp1
            end if
        end do
    end function f_orderId


    function p_min_max ( x, mn, mx, n ) result( out )
        ! correct the values to be within the minimum and maximum range

        implicit none

        ! input
        integer, intent(in) :: n
        real(kind=8), intent(in) :: mn, mx
        real(kind=8), dimension(n) :: x

        ! output
        real(kind=8), dimension(n) :: out

        where( x(:) > mx) x(:) = mx
        where( x(:) < mn) x(:) = mn

        out = x

    end function p_min_max


    function f_rsap ( age, pss_par, nsp, hb ) result( r )
        ! Get sapwood fraction based on 4C routine

        implicit none

        ! input
        integer, intent(in) :: nsp
        real(kind=8), dimension(nsp), intent(in) :: age, pss_par, hb

        ! output
        real(kind=8), dimension(nsp) :: r

        integer :: i, j
        integer, dimension(nsp) :: n, taumax
        real(kind=8), dimension(nsp) :: pss_aux

        n(:) = floor(age(:))
        r(:)=0.d0
        pss_aux(:) = 1.d0 - pss_par(:)

        do i = 1, nsp
            if (n(i) < 5.d0) then
                r(i) = 1.d0
            else
                if (hb(i) < 1.3d0) then
                    taumax(i)=n(i)-INT(hb(i)/1.3d0*5.d0)
                else
                    taumax(i)=n(i)-5.d0
                end if
            end if

            do j = 0, (taumax(i)-1)
                r(i) = r(i) + exp((j * log(pss_aux(i)))) * (2.d0*(taumax(i) - j) - 1.d0)
            end do
        end do

        r(:) = r(:) / taumax(:) ** 2.d0

        where( r(:) > 1.d0)
            r(:) = 1.d0
        end where

    end function f_rsap

    function f_heartwood (age, dbh, ba, stems, rsap, height, crown_length, height_sap, nsp ) result( r )
        ! Get sapwood/heartwood volume share

        implicit none

        ! input
        integer, intent(in) :: nsp
        real(kind=8), dimension(nsp), intent(in) :: age, dbh, ba, stems, rsap,  height, crown_length, height_sap

        ! local
        real(kind=8), dimension(nsp) :: h_vol, sapwood_ratio, dcb, dsb, dmid, acb, asb, amid, s_vol, hb, r

    	! Bole height
	    hb(:) = height(:) - crown_length(:)

        where (crown_length(:) <= 0.d0)
            hb(:) = height(:) * 0.5d0
        end where

        ! Diameter at crown base
        dcb(:) = dbh(:) * (1.37 - hb(:)) / height(:) + dbh(:)

        ! Diameter at mid section
        dmid(:) = dbh(:) * (1.37 - hb(:) * 0.5d0) / height(:) + dbh(:)

        ! Diameter at stem base
        dsb(:) = dbh(:) * 1.37d0 / height(:) + dbh(:)

        ! Area at crown base
        acb(:) = Pi * dcb(:)**2 / 40000.d0 * stems(:) - ba(:) * rsap(:) 

        where (acb(:) < 0.d0)
            acb(:) = 0.d0
        end where
        
        ! Area in mid section
        amid(:) = Pi * dmid(:)**2 / 40000.d0 * stems(:) * (1-rsap(:))

        !Area at stem base
        asb(:) = Pi * dsb(:)**2 / 40000.d0 * stems(:) * (1-rsap(:)) 

        ! Heartwood volume
        h_vol(:) =  (asb(:) + 4*amid(:) + acb(:))/6.d0 * height_sap(:)

        ! Based on form factor
        h_vol(:) = ba(:) * (1.d0-rsap(:)) * height_sap(:) * 0.42d0

        where(age(:) <= 5 .or. height(:) < 1.37d0)
            h_vol(:) = 0.d0
        end where

        ! Sapwood volume
        s_vol(:) = ba(:) * rsap(:) * height_sap(:)

        where (acb(:) == 0.d0)
            s_vol(:) = Pi * dcb(:)**2 / 40000.d0 * stems(:) * height_sap(:)
        end where

        ! Ratio
        r(:) = s_vol(:) / (h_vol(:) + s_vol(:))

    end function f_heartwood


    function f_elai ( lai, n ) result( out )
        ! get elai based on lai 

        implicit none

        ! input
        real(kind=8), dimension(n), intent(in) :: lai
        integer, intent(in) :: n
        real(kind=8), dimension(n) :: elai

        ! output
        real(kind=8), dimension(n) :: out
        integer :: i

        do i = 1, n
            if (lai(i) < 2.d0) then
                elai(i) = lai(i)
            else if (lai(i) >= 2.d0 .and. lai(i) < 4.d0) then
                elai(i) = 2.d0
            else
                elai(i) = lai(i) / 2.d0
            end if
        end do

        out = elai

    end function f_elai
    
    ! LPJ-GUESS stress index, used for the bark beetle module
    function f_greff(greff_sp, npp, lai, dmC, nyears, npp_corr) result(greff_out)

        !input
        integer, intent(in):: nyears
        real(kind=8), intent(in) :: greff_sp, dmC, npp_corr
        real(kind=8), dimension(nyears*12), intent(in) :: npp, lai

        !output
        real(kind=8) :: mort_n

        ! local
        integer :: i, n, cnt
        real(kind=8) :: greff, lai_avg, greff_out

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
        greff = npp_corr * 1000.d0 * (sum(npp) * dmC / (10.d0 * nyears)) / lai_avg ! convert to gC/m²/year

        greff_out = greff / greff_sp

    end function f_greff

    ! iLand stress index, used for the bark beetle module
    function f_si(npp , wf_loss, wr_loss) result(si)

        !input
        real(kind=8), dimension(12), intent(in):: npp , wf_loss, wr_loss

        !local
        real(kind=8):: si

        si = max(1.d0 - sum(npp) / (sum(wf_loss) + sum(wr_loss)), 0.d0)

    end function f_si

end module utils