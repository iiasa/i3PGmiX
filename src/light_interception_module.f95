module light_int

use mod_decl_const
use utils

IMPLICIT NONE
public :: s_light_3pgpjs, s_light_3pgmix, f_get_daylength, f_get_layer, f_get_layer_sum, f_get_solarangle

contains

	function f_get_daylength( Lat ) result( day_length )
		! Day-length calculations

		implicit none

		! input
		real(kind=8), intent(in) :: Lat

		! output
		real(kind=8), dimension(12) :: day_length

		! local
		real(kind=8) :: SLAt, cLat
		real(kind=8), dimension(12) :: sinDec, cosH0


		SLAt = sin(Pi * Lat / 180.d0)
		cLat = cos(Pi * Lat / 180.d0)
		sinDec(:) = 0.4d0 * sin(0.0172d0 * (dayOfYear(:) - 80.d0) )
		cosH0(:) = -sinDec(:) * SLAt / (cLat * sqrt(1.d0 - (sinDec(:)) ** 2.d0))

		day_length(:) = Acos(cosH0(:)) / Pi

		where( cosH0 > 1.d0 ) day_length = 0.d0
		where( cosH0 < -1.d0 ) day_length = 1.d0

	end function f_get_daylength

	function f_get_solarangle( Lat ) result( solarangle )

		implicit none

		! input
		real(kind=8), intent(in) :: Lat

		! output
		real(kind=8), dimension(12) :: solarangle

		! local
		real(kind=8) :: secondxaxisintercept, firstxaxisintercept
		real(kind=8), dimension(12) :: gamma, declinationangle, szaprep, solarzenithangle


		secondxaxisintercept = 0.0018d0 * Lat ** 3.d0 - 0.0031d0 * Lat ** 2.d0 + 2.3826d0 * Lat + 266.62d0
		firstxaxisintercept = -0.0018d0 * Lat ** 3.d0 + 0.0021d0 * Lat ** 2.d0 - 2.3459d0 * Lat + 80.097d0

		gamma(:) = 2.d0 * Pi / 365.d0 * ( dayOfYear(:) - 1.d0)

		declinationangle(:) = 0.006918d0 - (0.399912d0 * Cos(gamma(:))) + 0.070257d0 * Sin(gamma(:)) - &
			0.006758d0 * Cos(2.d0 * gamma(:)) + 0.000907d0 * Sin(2.d0 * gamma(:)) - 0.002697d0 * Cos(3.d0 * gamma(:)) + &
			0.00148d0 * Sin(3.d0 * gamma(:))

		szaprep(:) = Sin(Pi / 180.d0 * Lat * ( -1.d0) ) * Sin(declinationangle(:)) + &
			Cos(Pi / 180.d0 * Lat * (-1.d0) ) * Cos(declinationangle(:))
		solarzenithangle(:) = 180.d0 / Pi * (Atan(-szaprep(:) / ((-szaprep(:) * szaprep(:) + 1.d0) ** 0.5d0)) + 2.d0 * Atan(1.d0))

		solarangle(:) = solarzenithangle(:)

		if ( Lat >= 0.d0 .and. Lat <= 23.4d0) Then
			!the zenith angle only needs to be adjusted if the lat is between about -23.4 and 23.4
			where( dayOfYear(:) > secondxaxisintercept .or. dayOfYear(:) < firstxaxisintercept )
				solarangle(:) = -1.d0 * solarzenithangle(:)
			end where
		end if

		if (  Lat >= -23.4d0 .and. Lat < 0.d0 ) Then
			!the zenith angle only needs to be adjusted if the lat is between about -23.4 and 23.4
			where( dayOfYear(:) > firstxaxisintercept .and. dayOfYear(:) < secondxaxisintercept )
				solarangle(:) = -1.d0 * solarzenithangle(:)
			end where
		end if

	end function f_get_solarangle


	function f_get_layer ( n_sp, height, Heightcrown) result(layer_id)
		! function to allocate each tree to the layer based on height and crown heigh
		! First layer (1) is the highest
		! According to Forrester, D.I., Guisasola, R., Tang, X. et al. For. Ecosyst. (2014) 1: 17.
		! Calculations based on example https://it.mathworks.com/matlabcentral/answers/366626-overlapping-time-intervals

		implicit none

		integer, intent(in) :: n_sp ! number of species
		real(kind=8), dimension(n_sp), intent(in) :: height, Heightcrown

		! output
		integer, dimension(n_sp) :: layer_id ! array of layer id

		! local
		real(kind=8), dimension( n_sp*2 ) :: Height_all
		integer, dimension( n_sp*2 ) :: Height_ind
		integer, dimension( n_sp*2 ) :: ones,  ones_sum ! vector of 1, 0, -1 for calculation
		real(kind=8), allocatable, dimension(:) :: Height_layer ! maximum height of each layer

		integer :: i
		integer :: n_l

		! Sort all height and crown heigh
		Height_all = [Heightcrown(:), height(:)] ! put height and crown beginning into vector
		Height_ind = f_orderId(Height_all) ! sort the array

		! Assign index order for further calculations
		ones(:) = -1; ones(1:n_sp) = 1
		ones = ones(Height_ind)

	!   cumulative sum
		ones_sum(1) = ones(1)
		!if( n_sp > 1 ) then
			do i = 2, n_sp*2
				ones_sum(i) = ones_sum(i-1) + ones(i)
			end do
		!end if

		! Max height of each layer
		n_l = count(ones_sum == 0)
		allocate( Height_layer(n_l) )
		Height_layer(:) = 0
		Height_layer = Height_all(PACK(Height_ind, ones_sum == 0))

		! Assign layer to each species
		layer_id(:) = 1
		if( n_l > 1 ) then
			do i = 1, n_l-1
				where ( height(:) > Height_layer(i) ) layer_id(:) = i+1
			end do
		end if

		deallocate( Height_layer )

		! revert the order, so highest trees are 1 layer and lowest is n
		layer_id(:) = maxval( layer_id(:) ) - layer_id(:) + 1

	end function f_get_layer


	function f_get_layer_sum ( n_sp, nLayers, x, layer_id) result (y)
		! function to sum any array x, based on the vector of layers id

		implicit none

		! input
		integer, intent(in) :: n_sp, nLayers ! number of species and layers
		real(kind=8), dimension(n_sp), intent(in) :: x
		integer, dimension(n_sp), intent(in) :: layer_id

		! output
		real(kind=8), dimension(n_sp) :: y

		! local
		integer :: i = 1

		y(:) = 0.d0

		do i = 1, nLayers
			where ( layer_id(:) == i )
				y(:) = sum(x(:), mask=layer_id(:)==i)
			end where
		end do

	end function f_get_layer_sum


	subroutine s_light_3pgpjs ( n_sp, age, fullCanAge, k, lai, solar_rad, days_in_month, &
		canopy_cover, apar )

		implicit none

		! input
		integer, intent(in) :: n_sp ! number of species
		real(kind=8), dimension(n_sp), intent(in) :: age
		real(kind=8), dimension(n_sp), intent(in) :: fullCanAge     ! Age at canopy closure
		real(kind=8), dimension(n_sp), intent(in) :: k
		real(kind=8), dimension(n_sp), intent(in) :: lai
		real(kind=8), intent(in) :: solar_rad
		integer, intent(in) :: days_in_month

		! output
		real(kind=8), dimension(n_sp), intent(out) :: canopy_cover
		real(kind=8), dimension(n_sp), intent(out) :: apar

		! Additional variables for calculation distribution
		real(kind=8) :: RADt ! Total available radiation
		real(kind=8), dimension(n_sp) :: lightIntcptn


		canopy_cover(:) = 1.d0
		where (fullCanAge(:) > 0.d0 .and. age(:) < fullCanAge(:) )
			canopy_cover(:) = (age(:) + 0.01d0) / fullCanAge(:)
		end where

		lightIntcptn = (1.d0 - (Exp(-k * lai / canopy_cover)))

		RADt = solar_rad * days_in_month ! MJ m-2 month-1
		apar(:) = RADt * lightIntcptn(:) * canopy_cover(:)

	end subroutine s_light_3pgpjs


	subroutine s_light_3pgmix ( n_sp, height, crown_length, crown_width, lai, stems_n, solar_rad, &
		CrownShape, k, solarAngle,days_in_month, &
		apar, lai_above, fi, lambda_v, lambda_h, canopy_vol_frac, layer_id, lai_sa_ratio)

		! Subroutine calculate the apar for the mixed species forest
		! It first allocate each species to a specific layer based on height and crown length
		! and then distribute the light between those layers

		! If LAI is equal to 0, this is an indicator that the species is currently in the dormant period


		implicit none

		! input
		integer, intent(in) :: n_sp ! number of species
		real(kind=8), dimension(n_sp) :: height ! i'm not putting the intent(in) here as we modify those variables later
		real(kind=8), dimension(n_sp) :: crown_length ! i'm not putting the intent(in) here as we modify those variables later
		real(kind=8), dimension(n_sp), intent(in) :: crown_width
		real(kind=8), dimension(n_sp), intent(in) :: lai
		real(kind=8), dimension(n_sp), intent(in) :: stems_n
		real(kind=8), intent(in) :: solar_rad
		integer, dimension(n_sp), intent(in) :: CrownShape   !***DF crown shape of a given species; 1=cone, 2=ellipsoid, 3=half-ellipsoid, 4=rectangular
		real(kind=8), dimension(n_sp), intent(in) :: k
		real(kind=8), intent(in) :: solarAngle
		integer, intent(in) :: days_in_month

		! output
		real(kind=8), dimension(n_sp), intent(out) :: apar
		real(kind=8), dimension(n_sp), intent(out) :: lai_above !leaf area above the given species
		real(kind=8), dimension(n_sp), intent(out) :: fi !***DF the proportion of above canopy apar absorbed by each species
		real(kind=8), dimension(n_sp), intent(out) :: lambda_v       !Constant to partition light between species and to account for vertical canopy heterogeneity (see Equations 2 and 3 of Forrester et al., 2014, Forest Ecosystems, 1:17)
		real(kind=8), dimension(n_sp), intent(out) :: lambda_h         !Constant to account for horizontal canopy heterogeneity such as gaps between trees and the change in zenith angle (and shading) with latitude and season (see Equations 2 and 5 of Forrester et al., 2014, Forest Ecosystems, 1:17)
		real(kind=8), dimension(n_sp), intent(out) :: canopy_vol_frac !Fraction of canopy space (between lowest crown crown height to tallest height) filled by crowns
		integer, dimension(n_sp), intent(out) :: layer_id
		real(kind=8), dimension(n_sp), intent(out) :: lai_sa_ratio !the ratio of mean tree leaf area (m2) to crownSA (m2)

		! Additional variables for calculation distribution
		integer :: i
		real(kind=8), dimension(n_sp) :: Heightmidcrown    !mean height of the middle of the crown (height - height to crown base)/2 + height to crown base       !***DF
		real(kind=8), dimension(n_sp) :: Heightcrown ! height of the crown begining
		real(kind=8), dimension(n_sp) :: CrownSA  !mean crown surface area (m2) of a species
		real(kind=8), dimension(n_sp) :: Crownvolume   !***DF the crown volume of a given species
		integer :: nLayers ! number of layers
		real(kind=8), dimension(n_sp) :: Height_max_l
		real(kind=8), dimension(n_sp) :: Heightcrown_min_l
		real(kind=8), dimension(n_sp) :: Heightmidcrown_l ! maximum and minimum height of layer
		real(kind=8), dimension(n_sp) :: Heightmidcrown_r !ratio of the mid height of the crown of a given species to the mid height of a canopy layer
		real(kind=8), dimension(n_sp) :: kL_l          !sum of k x L for all species within the given layer
		real(kind=8), dimension(n_sp) :: lambdaV_l     ! sum of lambda_v per layer
		real(kind=8), dimension(n_sp) :: kLSweightedave   !calculates the contribution each species makes to the sum of all kLS products in a given layer (see Equation 6 of Forrester et al., 2014, Forest Ecosystems, 1:17)
		real(kind=8), dimension(n_sp) :: aparl  !The absorbed apar for the given  layer
		real(kind=8) :: RADt ! Total available radiation
		real(kind=8), dimension(n_sp) :: LAI_l ! Layer LAI

		! initialization
		CrownSA(:) = 0.d0
		Crownvolume(:) = 0.d0
		Height_max_l(:) = 0.d0
		Heightcrown_min_l(:) = 0.d0
		aparl(:) = 0.d0
		apar(:) = 0.d0
		lai_above(:) = 0.d0

		!Calculate the mid crown height, crown surface and volume
		! check if species is dormant
		! where( lai(:) == 0 )
		!    height(:) = 0.d0
		!    crown_length(:) = 0.d0
		! end where

		Heightcrown(:) = height(:) - crown_length(:)
		Heightmidcrown(:) = height(:) - crown_length(:) / 2


		! Calculate the crown area and volume
		! We only do it for species that have LAI, otherwise it stays 0 as was initialized above
		do i = 1, n_sp
			if( lai(i) > 0.d0 ) then
				if( CrownShape(i) == int(1) ) then !cone shaped
					CrownSA(i) = Pi * ((crown_width(i) / 2.d0) ** 2.d0) + Pi * crown_width(i) / 2.d0 * &
						(((crown_width(i) / 2.d0) ** 2.d0) + crown_length(i) ** 2.d0) ** 0.5d0
					Crownvolume(i) = Pi * crown_width(i) * crown_width(i) * crown_length(i) / 12.d0
				else if( CrownShape(i) == int(2) ) then !ellipsoid
					CrownSA(i) = 4.d0 * Pi * ((((crown_width(i) / 2.d0) ** 1.6075d0) * ((crown_width(i) / 2.d0) ** 1.6075d0) + &
						((crown_width(i) / 2.d0) ** 1.6075d0) * ((crown_length(i) / 2.d0) ** 1.6075d0) + &
						((crown_width(i) / 2.d0) ** 1.6075d0) * ((crown_length(i) / 2.d0) ** 1.6075d0)) / 3.d0) ** (1.d0 / 1.6075d0)
					Crownvolume(i) = Pi * crown_width(i) * crown_width(i) * crown_length(i) * 4.d0 / 24.d0
				else if( CrownShape(i) == int(3) ) then !half-ellipsoid
					CrownSA(i) = Pi * ((crown_width(i) / 2.d0) ** 2.d0) + (4.d0 * Pi * ((((crown_width(i) / 2.d0) ** 1.6075d0) * &
						((crown_width(i) / 2.d0) ** 1.6075d0) + ((crown_width(i) / 2.d0) ** 1.6075d0) * &
						((crown_length(i)) ** 1.6075d0) + ((crown_width(i) / 2.d0) ** 1.6075d0) * &
						((crown_length(i)) ** 1.6075d0)) / 3.d0) ** (1 / 1.6075d0)) / 2.d0
					Crownvolume(i) = Pi * crown_width(i) * crown_width(i) * crown_length(i) * 4.d0 / 24.d0
				else if( CrownShape(i) == int(4) ) then !rectangular
					CrownSA(i) = crown_width(i) * crown_width(i) * 2.d0 + crown_width(i) * crown_length(i) * 4.d0
					Crownvolume(i) = crown_width(i) * crown_width(i) * crown_length(i)
				end if
			end if
		end do


		!calculate the ratio of tree leaf area to crown surface area restrict kLS to 1
		lai_sa_ratio(:) = lai(:) * 10000.d0 / stems_n(:) / CrownSA(:)
		where ( lai(:) == 0.d0 ) lai_sa_ratio(:) = 0.d0


		! separate trees into layers
		layer_id(:) = f_get_layer(n_sp, height(:), Heightcrown(:) )
		!where ( lai(:) == 0.0d0 ) layer_id(:) = -1.d0
		nLayers = maxval( layer_id(:) )


		! Now calculate the proportion of the canopy space that is filled by the crowns. The canopy space is the
		! volume between the top and bottom of a layer that is filled by crowns in that layer.
		! We calculate it only for the trees that have LAI and are in that particular year. Thus the tree can be in that
		! layer, but currently will not have LAI
		do i = 1, nLayers
			where ( layer_id(:) == i )
				Height_max_l(:) = maxval(height(:), mask=layer_id(:) .eq. i .and. lai(:) .ne. 0.d0)
				Heightcrown_min_l(:) = minval(Heightcrown(:), mask=layer_id(:) .eq. i .and. lai(:) .ne. 0.d0)
			end where
		end do


		! sum the canopy volume fraction per layer and save it at each species
		canopy_vol_frac(:) = Crownvolume(:) * stems_n(:) / ( (Height_max_l(:) - Heightcrown_min_l(:)) * 10000.d0)
		canopy_vol_frac(:) = f_get_layer_sum(n_sp, nLayers, canopy_vol_frac(:), layer_id(:))

		! if the canopy volume fraction is < 0.01 (very small seedlings) then it is outside the range of the model there is no need for lambda_h so, make canopy_vol_frac = 0.01
		!where( canopy_vol_frac(:) < 0.01d0 ) canopy_vol_frac(:) = 0.01d0

		Heightmidcrown_l(:) = Heightcrown_min_l(:) + ( Height_max_l(:) - Heightcrown_min_l(:) ) / 2.d0

		!determine the ratio between the mid height of the given species and the mid height of the layer.
		Heightmidcrown_r(:) = Heightmidcrown(:) / Heightmidcrown_l(:)

		! Calculate the sum of kL for all species in a layer
		kL_l(:) =  k(:) * lai(:)
		kL_l(:) = f_get_layer_sum(n_sp, nLayers, kL_l(:), layer_id(:))


		! Constant to partition light between species and to account for vertical canopy heterogeneity
		! (see Equations 2 and 3 of Forrester et al., 2014, Forest Ecosystems, 1:17)
		lambda_v(:) = 0.012306d0 + 0.2366090d0 * k(:) * LAI(:) / kL_l(:) + 0.029118d0 * Heightmidcrown_r(:) + &
			0.608381d0 * k(:) * LAI(:) / kL_l(:) * Heightmidcrown_r(:)

		! check for dormant
		where ( lai(:) == 0.d0 )
			lambda_v(:) = 0.d0
		end where

		! make sure the sum of all lambda_v = 1
		lambdaV_l(:) = f_get_layer_sum(n_sp, nLayers, lambda_v(:), layer_id(:))

		where( lambdaV_l(:) .ne. 0.d0 )
			lambda_v(:) = lambda_v(:) / lambdaV_l(:)
		end where



		! Calculate the weighted kLS based on kL/sumkL
		kLSweightedave(:) = k(:) * lai_sa_ratio(:) * k(:) * lai(:) / kL_l(:)
		kLSweightedave(:) = f_get_layer_sum( n_sp, nLayers, kLSweightedave(:), layer_id(:))
		! the kLS should not be greater than 1 (based on the data used to fit the light model in Forrester et al. 2014)
		! This is because when there is a high k then LS is likely to be small.
		where( kLSweightedave(:) > 1.d0) kLSweightedave(:) = 1.d0

		!Constant to account for horizontal canopy heterogeneity such as gaps between trees and the change in zenith angle (and shading) with latitude and season (see Equations 2 and 5 of Forrester et al., 2014, Forest Ecosystems, 1:17)
		lambda_h(:) = 0.8285d0 + ((1.09498d0 - 0.781928d0 * kLSweightedave(:)) * 0.1d0 ** (canopy_vol_frac(:))) - &
			0.6714096d0 * 0.1d0 ** (canopy_vol_frac(:))
		if ( solarAngle > 30.d0 ) then
			lambda_h(:) = lambda_h(:) + 0.00097d0 * 1.08259d0 ** solarAngle
		end if
		! check for dormant
		where ( lai(:) == 0.0d0 )
			lambda_h(:) = 0.0d0
		end where


		RADt = solar_rad * days_in_month ! MJ m-2 month-1
		do i = 1, nLayers
			where ( layer_id(:) == i )
				aparl(:) = RADt * (1.d0 - 2.71828182845905d0 ** (-kL_l(:)))
			end where
			RADt = RADt - maxval(aparl(:), mask=layer_id(:)==i ) ! subtract the layer RAD from total
		end do

		! ***DF this used to have month in it but this whole sub is run each month so month is now redundant here.
		apar(:) = aparl(:) * lambda_h(:) * lambda_v(:)

		! The proportion of above canopy apar absorbed by each species. This is used for net radiation calculations in the gettranspiration sub
		fi(:) = apar(:) / (solar_rad * days_in_month)

		! calculate the LAI above the given species for within canopy VPD calculations
		LAI_l = f_get_layer_sum(n_sp, nLayers, LAI(:), layer_id(:))

		! now calculate the LAI of all layers above and part of the current layer if the species
		! is in the lower half of the layer then also take the proportion of the LAI above
		! the proportion is based on the Relative height of the mid crown

		do i = 1, n_sp
			lai_above(i) = sum( lai(:), mask = layer_id(:) < layer_id(i) )
			if ( Heightmidcrown_r(i) < 0.9999999999999d0 ) then
				lai_above(i) =  lai_above(i) + sum( LAI(:), mask = layer_id(:) == layer_id(i) ) * ( 1.d0-Heightmidcrown_r(i) )
			end if
		end do

	end subroutine s_light_3pgmix

end module light_int