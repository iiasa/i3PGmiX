module management

use utils
use phenology
use mod_decl_const

IMPLICIT NONE
public :: get_crop_trees, get_harvested_trees, remove_trees, replant_stand

contains
    ! This routine restarts the age related parameters and accounting variables after
    ! final harvestings and restarts the stand with established saplings, based on the
    ! species input parameters. In case more than one species is present and only one
    ! of them is harvested, the number of sappling is proportional to the basal area
    ! proportion of the species.
    subroutine replant_stand(age, age_m, SLA, fracBB, wood_density, gammaN, gammaF, &
                            f_age, stems_n, biom_foliage_debt, biom_foliage, biom_stem, &
                            biom_root, biom_tree, sapling_wf, sapling_ws, sapling_wr, &
                            n_sapling, leafgrow, leaffall,  MaxAge, rAge, nAge, &
                            SLA0, SLA1, tSLA, fracBB0, fracBB1, tBB,  rho0, rho1, tRho, &
                            gammaN0, gammaN1, tgammaN, ngammaN, gammaF1, gammaF0, tgammaF, &
                            crop_trees, dbh_crop, height_crop, basal_area_crop, biom_stem_crop, &
                            vol_crop, volume_cum, basal_area_prop, n_m, n_sp, ii, i, month) 

        ! input 
        integer, intent(in):: n_m, ii, i, month, n_sp
        integer, dimension(n_sp), intent(in):: leafgrow, leaffall, crop_trees
        real(kind=8), dimension(n_sp), intent(in):: sapling_wf, sapling_ws, sapling_wr, n_sapling, &
                                   MaxAge, rAge, nAge, basal_area_prop, SLA0, SLA1, tSLA, fracBB0, &
                                   fracBB1, tBB,  rho0, rho1, tRho, gammaN0, gammaN1, tgammaN, ngammaN, &
                                   gammaF1, gammaF0, tgammaF

        ! output
        real(kind=8), dimension(n_sp), intent(out):: stems_n, biom_foliage_debt, biom_foliage, biom_stem, &
                                    biom_root, biom_tree, dbh_crop, height_crop, basal_area_crop, &
                                    biom_stem_crop, vol_crop, volume_cum 

        real(kind=8), dimension(n_m,n_sp), intent(out):: age, age_m, SLA, fracBB, wood_density, gammaN, gammaF, &
                                    f_age

            ! Recompute age related parameters
            age(ii:n_m,i) = age(ii:n_m,i) - age(ii,i)
            age_m(ii:n_m,i) =  age(ii:n_m,i) !- 1.d0/12.d0
            !age_m(ii,i) =  age(ii,i)

            ! Specific leaf area
            SLA(ii:n_m,i) = f_exp( (n_m-ii+1), age_m(ii:n_m,i), SLA0(i), SLA1(i), tSLA(i), 2.d0)

            ! Branch and bark fraction
            fracBB(ii:n_m,i) = f_exp( (n_m-ii+1), age_m(ii:n_m,i), fracBB0(i), fracBB1(i), tBB(i), 1.d0)

            ! Wood density
            wood_density(ii:n_m,i) = f_exp( (n_m-ii+1), age_m(ii:n_m,i), rho0(i), rho1(i), tRho(i), 1.d0)

            ! Mortality rate
            gammaN(ii:n_m,i) = f_exp( (n_m-ii+1), age(ii:n_m,i), gammaN0(i), gammaN1(i), tgammaN(i), ngammaN(i))
            gammaF(ii:n_m,i) = f_exp_foliage( (n_m-ii+1), age_m(ii:n_m,i), gammaF1(i), gammaF0(i), tgammaF(i))

            ! restart cummulative harvested volume
            volume_cum(i) = 0.d0
            
            ! Crop trees
            if ( crop_trees(i) > 0.d0) then
                dbh_crop(i) = 0.d0
                height_crop(i) = 0.d0
                basal_area_crop(i) = 0.d0
                biom_stem_crop(i) = 0.d0
                vol_crop(i) = 0.d0
            end if

            ! age modifier
            if (nAge(i) == 0.d0) then
                f_age(:,i) = 1.d0
            else
            ! I'm not declaring relative age, but directly put it inside
                f_age(ii:n_m,i) = 1.d0 / (1.d0 + ( (age_m(ii:n_m,i) / MaxAge(i) ) / rAge(i)) ** nAge(i))
            end if

            stems_n(i) = n_sapling(i) * basal_area_prop(i)

            if( f_dormant(month, leafgrow(i), leaffall(i)) .eqv. .TRUE.) then
                biom_foliage_debt(i) = sapling_wf(i)*stems_n(i)/10.d0**6.d0
                biom_foliage(i) = 0.d0
            else
                biom_foliage(i) = sapling_wf(i)*stems_n(i)/10.d0**6.d0
            end if


            biom_stem(i) = sapling_ws(i)*stems_n(i)/10.d0**6.d0
            biom_root(i) = sapling_wr(i)*stems_n(i)/10.d0**6.d0
            biom_tree(i) = biom_stem(i) * 1000.d0 / stems_n(i)  ! kg/tree
 
    end subroutine replant_stand

    ! Remove trees from thinnings and harvesting operations and update stand structure
    subroutine remove_trees(mort_manag, stem_harv, stems_n, biom_foliage_debt, biom_foliage, biom_stem_harv, &
                            biom_stem, biom_root, harvesting, residues, root_residues, foliage_residues, hF, &
                            hR, hS, wood_density, fracBB, dormant, b_cor, final_harvesting)

        ! Input
        logical, intent(in) :: dormant, final_harvesting
        real(kind=8), intent(in):: mort_manag, hF, hR, hS, wood_density, fracBB
        
        ! Output
        logical, intent(out) :: b_cor
        real(kind=8), intent(out):: stem_harv, stems_n, biom_foliage_debt, biom_foliage, biom_stem_harv, &
                                    biom_stem, biom_root, harvesting, residues, root_residues, foliage_residues

        ! Check if this is a partial or final harvesting operations
        if (final_harvesting .eqv. .TRUE.) then

            stem_harv = stem_harv + stems_n
            biom_stem_harv = biom_stem_harv + biom_stem

            ! Get harvesting volume
            harvesting = harvesting + biom_stem * (1.d0 - fracBB) / wood_density

            ! Get harvesting residues (bark + branches)
            residues = residues + biom_stem * fracBB / wood_density
            
            ! Get root volume of harvested trees
            root_residues = root_residues + biom_root / wood_density

            ! Get foliage biomass of harvested trees
            foliage_residues =  foliage_residues + biom_foliage

            ! Update stand structure
            if( dormant .eqv. .TRUE. ) then
                biom_foliage_debt = 0.d0
            else
                biom_foliage = 0.d0
            end if

            biom_root = 0.d0
            biom_stem = 0.d0
            stems_n = 0.d0
        else
            ! Update remaining and harvested stems
            stem_harv = stems_n * mort_manag
            stems_n = stems_n * (1.d0 - mort_manag)

            ! Remove foliage biomass of removed trees
            if( dormant .eqv. .TRUE.) then
                biom_foliage_debt = biom_foliage_debt * (1.d0 - mort_manag * hF)
            else
                biom_foliage = biom_foliage * (1.d0 - mort_manag * hF)
            end if
            
            ! Update root biomass
            biom_root = biom_root  * (1.d0 - mort_manag * hR )

            ! Define harvested biomass
            biom_stem_harv = biom_stem  * mort_manag * hS

            ! Update stem biomass
            biom_stem = biom_stem  * (1.d0 - mort_manag * hS )
           
            ! Get Harvesting volume
            harvesting = biom_stem * (mort_manag * hS) * (1.d0 - fracBB) / wood_density

            ! Get harvesting residues (bark + branches)
            residues = biom_stem * (mort_manag * hS) * fracBB / wood_density

            ! Get root volume of harvested trees
            root_residues = biom_root * (mort_manag * hR ) / wood_density

            ! Get foliage biomass of harvested trees
            foliage_residues =  biom_foliage * (mort_manag * hF)

            b_cor = .TRUE.
        end if


    end subroutine remove_trees


    ! This is used to calculate the stand variables for crop trees as opposed to the whole stand.
    ! Currently 5 cm dbh classes are defined for the computation.
    ! The computation was also modified compared to the original version, to use the cummulative Weibull 
    ! distribution and avoid looping to find minimum and maximum DBH classes. Furthermore, to ensure consistency
    ! with the model, the rationale was modified to calculate the share of stem biomass of crop trees compared
    ! to the stem biomass of the whole stand, and the computation of further parameters follow from the stem biomass
    subroutine get_crop_trees(dbh_crop, vol_crop, basal_area_crop, biom_stem_crop, height_crop, dbh, age, &
                            height, stems_n, basal_area, volume, biom_stem, competition_total, &
                            aWs, nWs, aH, nHB, nHC, aV, nVB, nVH, nVBH, &
                            Dscale0, DscaleB, Dscalerh, Dscalet, DscaleC, &
                            Dshape0, DshapeB, Dshaperh, Dshapet, DshapeC, &
                            Dlocation0, DlocationB, Dlocationrh, Dlocationt, DlocationC, &
                            crop_trees, fracBB, wood_density, height_model, step, n_sp) 
        ! input
        integer, intent(in):: n_sp, height_model
        real(kind=8), intent(in) :: step ! DBH size class
        integer, dimension(n_sp), intent(in):: crop_trees
        real(kind=8), dimension(n_sp), intent(in):: aWs, nWs, aH, nHB, nHC, aV, nVB, nVH, nVBH, &
                                                    dbh, age, height, stems_n, competition_total, &
                                                    basal_area, volume, biom_stem, fracBB, wood_density
        
        ! Distribution parameters
        real(kind=8), dimension(n_sp), intent(in) :: Dscale0, DscaleB, Dscalerh, Dscalet, DscaleC
        real(kind=8), dimension(n_sp), intent(in) :: Dshape0, DshapeB, Dshaperh, Dshapet, DshapeC
        real(kind=8), dimension(n_sp), intent(in) :: Dlocation0, DlocationB, Dlocationrh, Dlocationt, DlocationC

        ! output
        real(kind=8), dimension(n_sp), intent(out):: dbh_crop, vol_crop, basal_area_crop, biom_stem_crop, height_crop

        ! local variables
        integer:: croptreecounter, i   ! counts the number of crop trees in the diameter classes already considered when determining the smallest diameter class containing crop trees
        real(kind=8):: sizeclass, Nclass, sizeclass_min, dbh_class, n_trees_class, n_low, n_up, biom_stem_crop_total
        real(kind=8), dimension(n_sp):: height_rel, HeightTemp, DWeibullShape, DWeibullScale, DWeibullLocation, &
                                        dlocation, DWeibullShape_gamma   

        ! Calculate the relative height
        height_rel(:) = height(:) / ( sum( height(:) * stems_n(:) ) / sum( stems_n(:) ) )

        where( Dlocation0(:)==0.d0 .and. DlocationB(:)==0.d0 .and. Dlocationrh(:)==0.d0 .and. &
            Dlocationt(:)==0.d0 .and. DlocationC (:)==0.d0 )   dlocation(:) = 0.d0

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


        do i=1, n_sp
        
            !only continue with this if crop trees have been requested and if the location parameter (for weibull
            !distributions) is more than 0.01. This is because location cannot be negative when bias is being
            !corrected but if it is restricted to be >=0 then unrealistic crop tree values will be predicted
            !until the location parameter would have been >0 without the restriction.
            if (crop_trees(i) > 0.d0 .and. DWeibullLocation(i) > 0.01d0) then

            !Find the largest size class with at least 1 tree, but first find the smallest size class with at least 1 tree.
                ! sizeclass = 1.d0

                ! if (sizeclass < DWeibullLocation(i)) then
                !     Nclass = 0.d0
                ! else
                !     Nclass = (DWeibullShape(i) / DWeibullScale(i) * (((sizeclass - DWeibullLocation(i)) / &
                !               DWeibullScale(i)) ** (DWeibullShape(i) - 1.d0)) * exp(-((sizeclass - DWeibullLocation(i)) / &
                !               DWeibullScale(i)) ** DWeibullShape(i))) * stems_n(i)
                ! end if

                ! if (Nclass >= 1.d0) then

                !     do while (Nclass >= 1.d0)

                !         sizeclass = sizeclass + 1.d0

                !         if (sizeclass < DWeibullLocation(i)) then
                !             Nclass = 0.d0
                !         else
                !             Nclass = (DWeibullShape(i) / DWeibullScale(i) * (((sizeclass - DWeibullLocation(i)) / &
                !                       DWeibullScale(i)) ** (DWeibullShape(i) - 1.d0)) * exp(-((sizeclass - DWeibullLocation(i)) / &
                !                       DWeibullScale(i)) ** DWeibullShape(i))) * stems_n(i)
                !         end if

                !     end do

                !     sizeclass = sizeclass - 1.d0
            
                ! else !firstly find the minimum tree size

                    ! sizeclass = 1.d0

                    ! if (sizeclass < DWeibullLocation(i)) then
                    !     Nclass = 0.d0
                    ! else
                    !     Nclass = (DWeibullShape(i) / DWeibullScale(i) * (((sizeclass - DWeibullLocation(i)) / &
                    !               DWeibullScale(i)) ** (DWeibullShape(i) - 1.d0)) * exp(-((sizeclass - DWeibullLocation(i)) / &
                    !               DWeibullScale(i)) ** DWeibullShape(i))) * stems_n(i)
                    ! end if

                    ! do while (Nclass < 1.d0)

                    !     sizeclass = sizeclass + 1.d0

                    !     if (sizeclass < DWeibullLocation(i)) then
                    !         Nclass = 0.d0
                    !     else
                    !         Nclass = (DWeibullShape(i) / DWeibullScale(i) * (((sizeclass - DWeibullLocation(i)) / &
                    !                  DWeibullScale(i)) ** (DWeibullShape(i) - 1.d0)) * exp(-((sizeclass - DWeibullLocation(i)) / &
                    !                  DWeibullScale(i)) ** DWeibullShape(i))) * stems_n(i)
                    !     end if

                    ! end do

                    ! !now find the largest size class
                    ! do while (Nclass >= 1.d0)

                    !     sizeclass = sizeclass + 1.d0

                    !     if (sizeclass < DWeibullLocation(i)) then
                    !         Nclass = 0.d0
                    !     else
                    !         Nclass = (DWeibullShape(i) / DWeibullScale(i) * (((sizeclass - DWeibullLocation(i)) / &
                    !                   DWeibullScale(i)) ** (DWeibullShape(i) - 1.d0)) * exp(-((sizeclass - DWeibullLocation(i)) / &
                    !                   DWeibullScale(i)) ** DWeibullShape(i))) * stems_n(i)
                    !     end if

                    ! end do

                    !sizeclass = sizeclass - 1.d0

                ! end if

                !Now loop back until the number of trees counted is equal to X (note that only variables
                !that are predicted from B (diameter) can be considered here)
                dbh_crop(i) = 0.d0
                height_crop(i) = 0.d0
                basal_area_crop(i) = 0.d0
                biom_stem_crop(i) = 0.d0
                vol_crop(i) = 0.d0

                if (stems_n(i) <= crop_trees(i)) then !all trees are then crop trees

                    dbh_crop(i) = dbh(i)
                    height_crop(i) = height(i)
                    basal_area_crop(i) = basal_area(i)
                    biom_stem_crop(i) = biom_stem(i)
                    vol_crop(i) = volume(i)

                else

                    ! croptreecounter = 0

                    ! do while (croptreecounter < crop_trees(i) .or. sizeclass >= 1.d0)
                        
                    !     if (sizeclass < DWeibullLocation(i)) then
                    !         Nclass = 0.d0
                    !     else
                    !         Nclass = (DWeibullShape(i) / DWeibullScale(i) * (((sizeclass - DWeibullLocation(i)) / &
                    !                   DWeibullScale(i)) ** (DWeibullShape(i) - 1.d0)) * exp(-((sizeclass - DWeibullLocation(i)) / &
                    !                   DWeibullScale(i)) ** DWeibullShape(i))) * stems_n(i)
                    !     end if

                    !     croptreecounter = croptreecounter + NINT(Nclass)
                        
                    !     avDBHCrop(i) = avDBHCrop(i) + Nclass * sizeclass
                    !     BasAreaCrop(i) = BasAreaCrop(i) + Nclass * sizeclass **2.d0 * Pi / 40000.d0
                    !     WSCrop(i) = WSCrop(i) + Nclass * aWs(i) * sizeclass ** nWs(i) / 1000.d0
                    !     HeightTemp(i) = (aH(i) * sizeclass ** nHB(i) * competition_total(i) ** nHC(i))
                    !     HeightCrop(i) = HeightCrop(i) + Nclass * HeightTemp(i)
                    !     StandVolCrop(i) = StandVolCrop(i) + Nclass * (aV(i) * sizeclass ** nVB(i) * HeightTemp(i) ** nVH(i) * &
                    !                      (sizeclass **2.d0 * HeightTemp(i)) ** nVBH(i))
                    
                    !     sizeclass = sizeclass - 1.d0

                    ! end do

                    sizeclass = NINT(DWeibullScale(i) * ( (-log( 1.d0 - (stems_n(i) - 1.d0) / &
                                stems_n(i))) ** (1.d0 / DWeibullShape(i)) ) + DWeibullLocation(i) - step)


                    do while (sizeclass > 0.d0)
            
                        ! if (sizeclass < DWeibullLocation(i) .or. croptreecounter < crop_trees(i)) then
                        !     Nclass = 0.d0
                        ! else
                        !     Nclass = (DWeibullShape(i) / DWeibullScale(i) * (((sizeclass - DWeibullLocation(i)) / &
                        !             DWeibullScale(i)) ** (DWeibullShape(i) - 1.d0)) * exp(-((sizeclass - DWeibullLocation(i)) / &
                        !             DWeibullScale(i)) ** DWeibullShape(i))) * stems_n(i)
                        ! end if
                        n_trees_class = 0.d0
                        
                        dbh_class = sizeclass + step/2.d0

                        if (sizeclass < DWeibullLocation(i)) then
                            n_trees_class = 0.d0
                        else 
                        ! Get number of trees in the DBH class
                            n_low = 1.d0 - exp(-((sizeclass - DWeibullLocation(i)) / DWeibullScale(i) )**DWeibullShape(i))
                            n_up = 1.d0 - exp(-((sizeclass + step - DWeibullLocation(i)) / DWeibullScale(i) )**DWeibullShape(i))
                            n_trees_class = NINT( stems_n(i) * (n_up - n_low) )
                        end if
                        
                        biom_stem_crop(i) = biom_stem_crop(i) + n_trees_class * aWs(i) * dbh_class ** nWs(i) / 1000.d0
                                            
                        sizeclass = sizeclass - step

                    end do

                    biom_stem_crop_total =  biom_stem_crop(i)

                    croptreecounter = 0

                    dbh_crop(i) = 0.d0
                    height_crop(i) = 0.d0
                    basal_area_crop(i) = 0.d0
                    biom_stem_crop(i) = 0.d0
                    vol_crop(i) = 0.d0

                    sizeclass_min = NINT(DWeibullScale(i) * ( (-log( 1.d0 - (stems_n(i) - crop_trees(i)) / &
                                    stems_n(i))) ** (1.d0 / DWeibullShape(i)) ) + DWeibullLocation(i) - step)
        
                    sizeclass = NINT(DWeibullScale(i) * ( (-log( 1.d0 - (stems_n(i) - 1.d0) / &
                                    stems_n(i))) ** (1.d0 / DWeibullShape(i)) ) + DWeibullLocation(i) - step)

                    do while (sizeclass > sizeclass_min)
                        
                        ! if (sizeclass < DWeibullLocation(i) .or. croptreecounter < crop_trees(i)) then
                        !     Nclass = 0.d0
                        ! else
                        !     Nclass = (DWeibullShape(i) / DWeibullScale(i) * (((sizeclass - DWeibullLocation(i)) / &
                        !               DWeibullScale(i)) ** (DWeibullShape(i) - 1.d0)) * exp(-((sizeclass - DWeibullLocation(i)) / &
                        !               DWeibullScale(i)) ** DWeibullShape(i))) * stems_n(i)
                        ! end if
                        ! n_trees_class = 0.d0
                        
                        dbh_class = sizeclass + step/2.d0

                        if (sizeclass < DWeibullLocation(i)) then
                            n_trees_class = 0.d0
                        else 
                        ! Get number of trees in the DBH class
                            n_low = 1.d0 - exp(-((sizeclass - DWeibullLocation(i)) / DWeibullScale(i) )**DWeibullShape(i))
                            n_up = 1.d0 - exp(-((sizeclass + step - DWeibullLocation(i)) / DWeibullScale(i) )**DWeibullShape(i))
                            n_trees_class = NINT( stems_n(i) * (n_up - n_low) )
                        end if

                        croptreecounter = croptreecounter + n_trees_class
                        
                        ! avDBHCrop(i) = avDBHCrop(i) + n_trees_class * dbh_class
                        ! BasAreaCrop(i) = BasAreaCrop(i) + n_trees_class * dbh_class **2.d0 * Pi / 40000.d0
                        biom_stem_crop(i) = biom_stem_crop(i) + n_trees_class * aWs(i) * dbh_class ** nWs(i) / 1000.d0

                        ! if( height_model .eq. 1 ) then
                        !     HeightTemp(i) = aH(i) * dbh_class ** nHB(i) * competition_total(i) ** nHC(i)
                        ! else if ( height_model .eq. 2 ) then
                        !     HeightTemp(i) = 1.3d0 + aH(i) * Exp(1.d0)**(-nHB(i)/dbh_class) + nHC(i) * &
                        !                     competition_total(i) * dbh_class
                        ! end if

                        ! HeightCrop(i) = HeightCrop(i) + n_trees_class * HeightTemp(i)
                                            
                        sizeclass = sizeclass - step

                    end do
            
                end if

            ! remove any extra trees that were in the smallest sizeclass with crop trees
                if (croptreecounter > crop_trees(i)) then

                    sizeclass = NINT(DWeibullScale(i) * ( (-log( 1.d0 - (stems_n(i) - crop_trees(i)) / &
                                stems_n(i))) ** (1.d0 / DWeibullShape(i)) ) + DWeibullLocation(i) - step)

                    dbh_class = sizeclass + step/2.d0

                    n_trees_class = croptreecounter - crop_trees(i)

                !     !do while (croptreecounter /= crop_trees(i))

                        croptreecounter = croptreecounter - n_trees_class
                        ! avDBHCrop(i) = avDBHCrop(i) - n_trees_class * dbh_class
                        ! BasAreaCrop(i) = BasAreaCrop(i) - n_trees_class * dbh_class ** 2.d0 * Pi / 40000.d0
                        biom_stem_crop(i) = biom_stem_crop(i) - n_trees_class * aWs(i) * dbh_class ** nWs(i) / 1000.d0
                        ! if( height_model .eq. 1 ) then
                        !     HeightTemp(i) = aH(i) * dbh_class ** nHB(i) * competition_total(i) ** nHC(i)
                        ! else if ( height_model .eq. 2 ) then
                        !     HeightTemp(i) = 1.3d0 + aH(i) * Exp(1.d0)**(-nHB(i)/dbh_class) + nHC(i) * &
                        !                     competition_total(i) * dbh_class
                        ! end if

                        ! HeightCrop(i) = HeightCrop(i) - n_trees_class * HeightTemp(i)

                !     !end do

                end if

                biom_stem_crop(i) = biom_stem(i) * biom_stem_crop(i) / biom_stem_crop_total 
                dbh_crop(i) = ( ( biom_stem(i) * 1000.d0 / crop_trees(i) ) / aWs(i)) ** ( 1.d0 / nWs(i) ) 
                basal_area_crop(i) = crop_trees(i) * dbh_crop(i) ** 2.d0 * Pi / 40000.d0
                if( height_model .eq. 1 ) then
                    height_crop(i) = aH(i) * dbh_crop(i) ** nHB(i) * competition_total(i) ** nHC(i)
                else if ( height_model .eq. 2 ) then
                    height_crop(i) = 1.3d0 + aH(i) * Exp(1.d0)**(-nHB(i)/dbh_crop(i)) + nHC(i) * &
                                    competition_total(i) * dbh_crop(i)
                end if

                !convert DBH and height to means rather than sums
                !avDBHCrop(i) = avDBHCrop(i) / crop_trees(i)
                !HeightCrop(i) = HeightCrop(i) / crop_trees(i)

                !get crop tree volume
                vol_crop(i) = vol_crop(i) + (biom_stem_crop(i) * (1.d0-fracBB(i)) / wood_density(i))
            ! else 
                ! dbh_crop(i) = 0.d0
                ! height_crop(i) = 0.d0
                ! basal_area_crop(i) = 0.d0
                ! biom_stem_crop(i) = 0.d0
                ! vol_crop(i) = 0.d0
            end if !if crop_trees(i) was not > 0 then dont do any calculations for this species
            
        end do 

    end subroutine get_crop_trees

    ! Routine to compute the dbh and height of harvested trees, used for defining assortments in
    ! a post-processing step.
    subroutine get_harvested_trees(dbh_harv, height_harv, aWs, nWs, aH, nHB, nHC, competition_total, &
                                    biom_stem_harv, stems_n, n_sp, height_model)
        
        ! input
        integer, intent(in):: n_sp, height_model

        real(kind=8), dimension(n_sp), intent(in):: aWs, nWs, aH, nHB, nHC, stems_n, competition_total, &
                                                    biom_stem_harv

        ! local
        integer:: i                                           

        ! output
        real(kind=8), dimension(n_sp), intent(out):: dbh_harv, height_harv

        do i=1, n_sp
            if (biom_stem_harv(i) > 0.d0) then

                ! get average dbh
                dbh_harv(i) = ( ( biom_stem_harv(i) * 1000.d0 /stems_n(i) ) / aWs(i)) ** ( 1.d0 / nWs(i) ) 

                ! get average height
                if( height_model .eq. 1 ) then
                    height_harv(i) = aH(i) * dbh_harv(i) ** nHB(i) * competition_total(i) ** nHC(i)
                else if ( height_model .eq. 2 ) then
                    height_harv(i) = 1.3d0 + aH(i) * Exp(1.d0)**(-nHB(i)/dbh_harv(i)) + nHC(i) * &
                                    competition_total(i) * dbh_harv(i)
                end if

            end if
        end do

    end subroutine get_harvested_trees
end module management