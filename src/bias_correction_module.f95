module bias_correction

use mod_decl_const
use utils

IMPLICIT NONE
public :: s_sizeDist_correct

contains

    subroutine s_sizeDist_correct (n_sp, age, stems_n, biom_tree, competition_total, lai, &
        correct_bias, height_model, pars_s, pars_b, aWs, nWs, pfsPower, pfsConst, &
        dbh, basal_area, height, crown_length, crown_width, pFS, bias_scale)

        ! Diameter distributions are used to correct for bias when calculating pFS from mean dbh, and ws distributions are
        ! used to correct for bias when calculating mean dbh from mean ws. This bias is caused by Jensen's inequality and is
        ! corrected using the approach described by Duursma and Robinson (2003) FEM 186, 373-380, which uses the CV of the
        ! distributions and the exponent of the relationship between predicted and predictor variables.

        ! The default is to ignore the bias. The alternative is to correct for it by using empirically derived weibull distributions
        ! from the weibull parameters provided by the user. If the weibull distribution does not vary then just provide scale0 and shape0.

        implicit none

        ! input
        integer, intent(in) :: n_sp ! number of species
        real(kind=8), dimension(n_sp), intent(in) :: age
        real(kind=8), dimension(n_sp), intent(in) :: stems_n
        real(kind=8), dimension(n_sp), intent(in) :: biom_tree
        real(kind=8), dimension(n_sp), intent(in) :: competition_total
        real(kind=8), dimension(n_sp), intent(in) :: lai

        ! parameters
        integer, intent(in) :: correct_bias ! if the distribution shall be fitted
        integer, intent(in) :: height_model ! which heigh equation
        real(kind=8), dimension(17, n_sp), intent(in) :: pars_s ! parameters for bias
        real(kind=8), dimension(30, n_sp), intent(in) :: pars_b ! parameters for bias
        real(kind=8), dimension(n_sp), intent(in) :: aWs, nWs
        real(kind=8), dimension(n_sp), intent(in) :: pfsPower, pfsConst

        ! output
        real(kind=8), dimension(n_sp), intent(inout) :: dbh
        real(kind=8), dimension(n_sp), intent(inout) :: basal_area
        real(kind=8), dimension(n_sp), intent(inout) :: height
        real(kind=8), dimension(n_sp), intent(out) :: crown_length
        real(kind=8), dimension(n_sp), intent(out) :: crown_width
        real(kind=8), dimension(n_sp), intent(out) :: pFS
        real(kind=8), dimension(15, n_sp), intent(out) :: bias_scale


        ! Variables and parameters
        real(kind=8), dimension(n_sp) :: lai_total
        real(kind=8), dimension(n_sp) :: height_rel

        real(kind=8), dimension(n_sp) :: aH, nHB, nHC
        real(kind=8), dimension(n_sp) :: aV, nVB, nVH, nVBH
        real(kind=8), dimension(n_sp) :: aK, nKB, nKH, nKC, nKrh
        real(kind=8), dimension(n_sp) :: aHL, nHLB, nHLL, nHLC, nHLrh
        real(kind=8), dimension(n_sp) :: Dscale0, DscaleB, Dscalerh, Dscalet, DscaleC
        real(kind=8), dimension(n_sp) :: Dshape0, DshapeB, Dshaperh, Dshapet, DshapeC
        real(kind=8), dimension(n_sp) :: Dlocation0, DlocationB, Dlocationrh, Dlocationt, DlocationC
        real(kind=8), dimension(n_sp) :: wsscale0, wsscaleB, wsscalerh, wsscalet, wsscaleC
        real(kind=8), dimension(n_sp) :: wsshape0, wsshapeB, wsshaperh, wsshapet, wsshapeC
        real(kind=8), dimension(n_sp) :: wslocation0, wslocationB, wslocationrh, wslocationt, wslocationC

        ! Additional variables for calculation distribution
        real(kind=8), dimension(n_sp) :: DWeibullScale, DWeibullShape, DWeibullLocation
        real(kind=8), dimension(n_sp) :: wsWeibullScale, wsWeibullShape, wsWeibullLocation
        real(kind=8), dimension(n_sp) :: Ex, Varx, CVdbhDistribution, CVwsDistribution
        real(kind=8), dimension(n_sp) :: DrelBiaspFS, DrelBiasheight, DrelBiasBasArea, DrelBiasLCL, DrelBiasCrowndiameter
        real(kind=8), dimension(n_sp) :: wsrelBias

        real(kind=8), dimension(n_sp) :: dlocation, wslocation
        real(kind=8), dimension(n_sp) :: DWeibullShape_gamma, wsWeibullShape_gamma

        include 'i_read_param_sizeDist.h'
        include 'i_read_param_sub.h'

        bias_scale(:,:) = 0.d0

        ! LAI
        lai_total = sum( lai(:) )

        ! Calculate the relative height
        ! height(:) = aH(:) * dbh(:) ** nHB(:) * competition_total(:) ** nHC(:)
        height_rel(:) = height(:) / ( sum( height(:) * stems_n(:) ) / sum( stems_n(:) ) )


        ! Check where all the locations are provided
        dlocation(:) = 1.d0
        wslocation(:) = 1.d0

        where( Dlocation0(:)==0.d0 .and. DlocationB(:)==0.d0 .and. Dlocationrh(:)==0.d0 .and. &
            Dlocationt(:)==0.d0 .and. DlocationC (:)==0.d0 ) dlocation(:) = 0.d0

        where( wslocation0(:)==0.d0 .and. wslocationB(:)==0.d0 .and. wslocationrh(:)==0.d0 .and. &
            wslocationt(:)==0.d0 .and. wslocationC (:)==0.d0 ) wslocation(:) = 0.d0

        if (correct_bias .eq. 1 ) then

            ! Calculate the DW scale -------------------
            DWeibullScale(:) = Exp( Dscale0(:) + DscaleB(:) * Log(dbh(:)) + Dscalerh(:) * Log(height_rel(:)) + &
                Dscalet(:) * Log(age(:)) + DscaleC(:) * Log(competition_total(:)))

            DWeibullShape(:) = Exp( Dshape0(:) + DshapeB(:) * Log( dbh(:) ) + Dshaperh(:) * Log(height_rel(:)) + &
                Dshapet(:) * Log(age(:)) + DshapeC(:) * Log(competition_total(:)))

            DWeibullShape_gamma(:) = f_gamma_dist(1.d0 + 1.d0 / DWeibullShape(:), n_sp)

            DWeibullLocation(:) = Exp( Dlocation0(:) + DlocationB(:) * Log(dbh(:)) + &
                    Dlocationrh(:) * Log(height_rel(:)) + Dlocationt(:) * Log(age(:)) + &
                    DlocationC(:) * Log(competition_total(:)))

            where( dlocation(:) == 0.d0 )
                DWeibullLocation(:) = NINT(dbh(:)) / 1.d0 - 1.d0 - DWeibullScale(:) * DWeibullShape_gamma(:)
            end where

            where( DWeibullLocation(:) < 0.01d0 ) DWeibullLocation(:) = 0.01d0


            Ex(:) = DWeibullLocation(:) + DWeibullScale(:) * DWeibullShape_gamma(:)
            !now convert the Ex from weibull scale to actual scale of diameter units in cm
            Varx(:) = DWeibullScale(:) ** 2.d0 * (f_gamma_dist(1.d0 + 2.d0 / DWeibullShape(:), n_sp) -  DWeibullShape_gamma ** 2.d0)
            CVdbhDistribution(:) = Varx(:) ** 0.5d0 / Ex(:)

            ! calculate the bias
            DrelBiaspFS(:) = 0.5d0 * (pfsPower(:) * (pfsPower(:) - 1.d0)) * CVdbhDistribution(:) ** 2.d0
            DrelBiasheight(:) = 0.5d0 * (nHB(:) * (nHB(:) - 1.d0)) * CVdbhDistribution(:) ** 2.d0
            DrelBiasBasArea(:) = 0.5d0 * (2.d0 * (2.d0 - 1.d0)) * CVdbhDistribution(:) ** 2.d0
            DrelBiasLCL(:) = 0.5d0 * (nHLB(:) * (nHLB(:) - 1.d0)) * CVdbhDistribution(:) ** 2.d0
            DrelBiasCrowndiameter(:) = 0.5d0 * (nKB(:) * (nKB(:) - 1.d0)) * CVdbhDistribution(:) ** 2.d0

            ! prevent unrealisticly large bias, by restricting it to within + or - 50%
            DrelBiaspFS(:) = p_min_max( DrelBiaspFS(:), -0.5d0, 0.5d0, n_sp)
            DrelBiasheight(:) = p_min_max( DrelBiasheight(:), -0.5d0, 0.5d0, n_sp)
            DrelBiasBasArea(:) = p_min_max( DrelBiasBasArea(:), -0.5d0, 0.5d0, n_sp)
            DrelBiasLCL(:) = p_min_max( DrelBiasLCL(:), -0.5d0, 0.5d0, n_sp)
            DrelBiasCrowndiameter(:) = p_min_max( DrelBiasCrowndiameter(:), -0.5d0, 0.5d0, n_sp)


            ! Calculate the biom_stem scale -------------------
            wsWeibullScale(:) = Exp( wsscale0(:) + wsscaleB(:) * Log(dbh(:)) + wsscalerh(:) * Log(height_rel(:)) + &
                wsscalet(:) * Log(age(:)) + wsscaleC(:) * Log(competition_total(:)))

            wsWeibullShape(:) = Exp( wsshape0(:) + wsshapeB(:) * Log(dbh(:)) + wsshaperh(:) * Log(height_rel(:)) + &
                wsshapet(:) * Log(age(:)) + wsshapeC(:) * Log(competition_total(:)))
            wsWeibullShape_gamma = f_gamma_dist(1.d0 + 1.d0 / wsWeibullShape(:), n_sp)

            wsWeibullLocation(:) = Exp( wslocation0(:) + wslocationB(:) * Log(dbh(:)) + &
                    wslocationrh(:) * Log(height_rel(:)) + wslocationt(:) * Log(age(:)) + &
                    wslocationC(:) * Log(competition_total(:)))

            where( wslocation(:) == 0.d0 )
                wsWeibullLocation(:) = NINT(biom_tree(:)) / 10.d0 - 1.d0 - wsWeibullScale(:) * wsWeibullShape_gamma(:)
            end where

            where( wsWeibullLocation(:) < 0.01d0 ) wsWeibullLocation(:) = 0.01d0


            Ex(:) = wsWeibullLocation(:) + wsWeibullScale(:) * wsWeibullShape_gamma
            !now convert the Ex from weibull scale to actual scale of diameter units in cm
            Varx(:) = wsWeibullScale(:) ** 2.d0 * (f_gamma_dist(1.d0 + 2.d0 / wsWeibullShape(:), n_sp) - &
                wsWeibullShape_gamma ** 2.d0)
            CVwsDistribution(:) = Varx(:) ** 0.5d0 / Ex(:)

            wsrelBias(:) = 0.5d0 * (1.d0 / nWs(:) * (1.d0 / nWs(:) - 1.d0)) * CVwsDistribution(:) ** 2.d0 !DF the nWS is replaced with 1/nWs because the equation is inverted to predict dbh from ws, instead of ws from dbh

            wsrelBias(:) = p_min_max( wsrelBias(:), -0.5d0, 0.5d0, n_sp)

        else
            DrelBiaspFS(:) = 0.d0
            DrelBiasBasArea(:) = 0.d0
            DrelBiasheight(:) = 0.d0
            DrelBiasLCL(:) = 0.d0
            DrelBiasCrowndiameter(:) = 0.d0
            wsrelBias(:) = 0.d0
        end if

        ! Correct for trees that have age 0 or are thinned (e.g. n_trees = 0)
        where( age(:) .eq. 0.d0 .or. stems_n(:) .eq. 0.d0 ) 
            DrelBiaspFS(:) = 0.d0
            DrelBiasBasArea(:) = 0.d0
            DrelBiasheight(:) = 0.d0
            DrelBiasLCL(:) = 0.d0
            DrelBiasCrowndiameter(:) = 0.d0
            wsrelBias(:) = 0.d0
        end where

        ! Avoid issues when planting mixed stands at initial age
        where( isnan(DrelBiaspFS(:)) ) 
            DrelBiaspFS(:) = 0.d0
        end where

        where( isnan(DrelBiasBasArea(:)) ) !causes problems for planting mixed stands at initial age
            DrelBiasBasArea(:) = 0.d0
        end where

        where( isnan(DrelBiasheight(:)) ) !causes problems for planting mixed stands at initial age
            DrelBiasheight(:) = 0.d0
        end where

        where( isnan(DrelBiasLCL(:)) ) !causes problems for planting mixed stands at initial age
            DrelBiasLCL(:) = 0.d0
        end where

        where( isnan(DrelBiasCrowndiameter(:)) ) !causes problems for planting mixed stands at initial age
            DrelBiasCrowndiameter(:) = 0.d0
        end where
        
        where( isnan(wsrelBias(:)) ) !causes problems for planting mixed stands at initial age
            wsrelBias(:) = 0.d0
        end where
        

        ! Correct for bias ------------------
        dbh(:) = (biom_tree(:) / aWs(:)) ** (1.d0 / nWs(:)) * (1.d0 + wsrelBias(:))
        basal_area(:) = ( dbh(:) ** 2.d0 / 4.d0 * Pi * stems_n(:) / 10000.d0) * (1.d0 + DrelBiasBasArea(:))

        if( height_model .eq. 1 ) then

            height(:) = ( aH(:) * dbh(:) ** nHB(:) * competition_total(:) ** nHC(:)) * (1.d0 + DrelBiasheight(:))

            crown_length(:) = ( aHL(:) * dbh(:) ** nHLB(:) * lai_total ** nHLL(:) * competition_total(:) ** nHLC(:) * &
                height_rel(:) ** nHLrh(:)) * (1.d0 + DrelBiasLCL(:))

        else if ( height_model .eq. 2 ) then

            height(:) = 1.3d0 + aH(:) * exp(1.d0)**(-nHB(:)/dbh(:)) + nHC(:) * competition_total(:) * dbh(:)
            crown_length(:) = 1.3d0 + aHL(:) * exp(1.d0)**(-nHLB(:)/dbh(:)) + nHLC(:) * competition_total(:) * dbh(:)
        end if

        crown_width(:) = min(( aK(:) * dbh(:) ** nKB(:) * height(:) ** nKH(:) * competition_total(:) ** nKC(:) * &
            height_rel(:) ** nKrh(:)) * (1.d0 + DrelBiasCrowndiameter(:)), 10.d0)
        where( lai(:) .eq. 0.d0 ) crown_width(:) = 0.d0


        pFS(:) = ( pfsConst(:) * dbh(:) ** pfsPower(:)) * (1.d0 + DrelBiaspFS(:))


        ! check that the height and LCL allometric equations have not predicted that height - LCL < 0
        ! and if so reduce LCL so that height - LCL = 0 (assumes height allometry is more reliable than LCL allometry)
        where ( crown_length(:) > height(:) )
            crown_length(:) = height(:)
        end where

        ! output the matrix of biasses
        bias_scale(1,:) = DWeibullScale(:)
        bias_scale(2,:) = DWeibullShape(:)
        bias_scale(3,:) = DWeibullLocation(:)
        bias_scale(4,:) = wsWeibullScale(:)
        bias_scale(5,:) = wsWeibullShape(:)
        bias_scale(6,:) = wsWeibullLocation(:)
        bias_scale(7,:) = CVdbhDistribution(:)
        bias_scale(8,:) = CVwsDistribution(:)
        bias_scale(9,:) = wsrelBias(:)
        bias_scale(10,:) = DrelBiaspFS(:)
        bias_scale(11,:) = DrelBiasheight(:)
        bias_scale(12,:) = DrelBiasBasArea(:)
        bias_scale(13,:) = DrelBiasLCL(:)
        bias_scale(14,:) = DrelBiasCrowndiameter(:)
        bias_scale(15,:) = height_rel


    end subroutine s_sizeDist_correct

end module bias_correction