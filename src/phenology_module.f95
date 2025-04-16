module phenology

IMPLICIT NONE
public :: f_dormant!, f_dormant_dyn

contains

	function f_dormant(month, leafgrow, leaffall) result( out )

		implicit none

		! input
		integer, intent(in) :: month, leafgrow, leaffall

		! output
		logical :: out

		out = .FALSE.
		
		! This is called if the leafgrow parameter is not 0, and hence the species is Deciduous
		! This is true if "currentmonth" is part of the dormant season
		if ( leafgrow > leaffall ) then
			! check which hemisphere
			if  ( month >= leaffall .and. month <= leafgrow) then ! growing at winter
				out = .TRUE.
			end if
		else if ( leafgrow < leaffall ) then
			if ( month < leafgrow .or. month >= leaffall ) then ! growing at summer
				out = .TRUE.
			end if
		end if

	end function f_dormant


	! function f_dormant_dyn(month, leafgrow, leaffall, gdd) result( out )

	! 	implicit none

	! 	! input
	! 	integer, intent(in) :: month, leafgrow, leaffall
	!	real(kind=8):: gdd

	! 	! output
	! 	logical :: out

	! 	out = .FALSE.
		
	! 	! This is called if the leafgrow parameter is not 0, and hence the species is Deciduous
	! 	! This is true if "currentmonth" is part of the dormant season
	! 	if ( leafgrow > leaffall ) then
	! 		! check which hemisphere
	! 		if  ( month >= leaffall .and. month <= leafgrow) then ! growing at winter
	! 			out = .TRUE.
	! 		end if
	! 	else if ( leafgrow < leaffall ) then
	! 		if ( month < leafgrow .or. month >= leaffall ) then ! growing at summer
	! 			out = .TRUE.
	! 		end if
	! 	end if

	! end function f_dormant_dyn


	! climate sensitive phenology
	!do i = 1, n_sp
	!	if (month .eq. leaffall_base(i) - int(1) .and. tmp_min(ii) < 4.3d0) then
	!		leaffall(i) = leaffall_base(i) - int(1)
	!	else if (month .eq. leaffall_base(i) .and. tmp_min(ii) >= 4.3d0) then
	!		leaffall(i) = leaffall_base(i) + int(1)
	!	end if
		
	!	if (month .eq. leafgrow_base(i) - int(1) .and. tmp_min(ii) >= 4.3d0) then
	!		leafgrow(i) = leafgrow_base(i) - int(1)
	!	else if (month .eq. leafgrow_base(i) .and. tmp_min(ii) < 4.3d0) then
	!		leafgrow(i) = leafgrow_base(i) + int(1)
	!	end if
	
	!end do

end module phenology