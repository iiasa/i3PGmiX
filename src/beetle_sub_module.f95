module beetle_subroutines

	use mod_decl_const

	implicit none
	
	public :: rand_gamma_scalar, calc_rdd
	
	contains

	subroutine rand_normal(x, mean,stdev) 

		real(kind=8):: mean,stdev,c,r,theta
		real(kind=8), dimension(2):: temp
		real(kind=8), intent(out):: x

		if(stdev <= 0.0d0) then
		   !Write(*,*) "Standard Deviation must be +ve"
		   x=mean
		else
		   call random_number(temp)
		   r=(-2.0d0*log(temp(1)))**0.5d0
		   theta = 2.0d0*Pi*temp(2)
		   x= mean+stdev*r*sin(theta)
		end if
	  end subroutine rand_normal

	! Returns a psuedorandom scalar drawn from the gamma distribution.
	subroutine rand_gamma0(a, first, fn_val)
		
		!
		! The shape parameter a >= 1.
		!
		! [1] Marsaglia, G., & Tsang, W. W. (2000). A Simple Method for Generating
		! Gamma Variables. ACM Transactions on Mathematical Software (TOMS), 26(3),
		! 363–372.
		real(kind=8), intent(in) :: a
		logical, intent(in) :: first
		real(kind=8), intent(out) :: fn_val
		real(kind=8), save :: c, d
		real(kind=8) :: U, v, x
		if (a < 1) then 
			stop
		end if 
		if (first) then
			d = a - 1.d0/3.d0
			c = 1.d0/sqrt(9.d0*d)
		end if
		do
			do
				call rand_normal(x, 0.d0, 1.d0)
				v = (1.d0 + c*x)**3.d0
				if (v > 0.d0) exit
			end do
			call random_number(U)
			! Note: the number 0.0331 below is exact, see [1].
			if (U < 1.d0 - 0.0331d0*x**4.d0) then
				fn_val = d*v
				exit
			else if (log(U) < x**2.d0/2.d0 + d*(1.d0 - v + log(v))) then
				fn_val = d*v
				exit
			end if
		end do
		end subroutine
	
	! Gamma random number
	subroutine rand_gamma_scalar(a, x)
		real(kind=8), intent(in) :: a
		real(kind=8), intent(out) :: x

		call rand_gamma0(a, .true., x)
	end subroutine

	! Get month of diapause and diapause day
	subroutine get_diapause_month(diapause_month, month_share, lat, day_length_h)

		!input
		real(kind=8), intent(in) :: lat
		real(kind=8), dimension(12), intent(in) :: day_length_h

		!output
		integer, intent(out):: diapause_month
		real(kind=8), intent(out) :: month_share

		!local
		integer:: index, start_day, end_day, i, dim

		index = 0
		! start in April
		do i = 4, 10
			if (day_length_h(i) > 14.5d0) then
				index = i
				exit
			endif
		end do

		! Set diapause month
		diapause_month = index

		! Get DOY limits from the diapause month
		start_day = startDay(diapause_month)
		end_day = endDay(diapause_month)
		dim = daysInMonth(diapause_month)

		! Number of days until diapause
		month_share = diapause_day(dim, start_day, end_day, lat) / dim

	end subroutine get_diapause_month

	! Get  diapause day
	function diapause_day(days_in_month , start_day, end_day, lat) result(diap_day)

		!input
		integer, intent(in):: days_in_month, start_day, end_day
		real(kind=8), intent(in) :: lat

		!local
		integer:: cnt, i, diap_day
		integer, dimension(days_in_month):: doy, day_length, lower, upper, mid
		real(kind=8):: dl

		cnt = 1
		do i = start_day, end_day
			doy(cnt) =  i
			cnt = cnt + 1
		end do

		day_length = get_day_len_doy(doy, Lat, days_in_month)

		! diap_day = 1
		! dl = day_length(diap_day)
		! do while(dl < 14.5d0 .or. diap_day < days_in_month)
		! 	diap_day = diap_day + 1 
		! 	dl = day_length(diap_day)
		! end do
		

		!  ! Initialize bounds
		lower = 1
		upper = days_in_month
		diap_day = 0
		
		! Binary search
		do i = 1, days_in_month

			dl = day_length(i)

			if (dl >= 14.5d0) exit

		end do

		diap_day = i

	end function diapause_day

	! Get day lengths of diapause month
	function get_day_len_doy(doy, Lat, days_in_month) result(day_len)

		! input
		integer, intent(in) :: days_in_month
		integer, dimension(days_in_month), intent(in) :: doy
		real(kind=8), intent(in) :: Lat

		! output
		real(kind=8), dimension(days_in_month) :: day_len

		! local
		real(kind=8) :: SLAt, cLat
		real(kind=8), dimension(days_in_month) :: sinDec, cosH0


		SLAt = sin(Pi * Lat / 180.d0)
		cLat = cos(Pi * Lat / 180.d0)
		sinDec(:) = 0.4d0 * sin(0.0172d0 * (doy(:) - 80.d0) )
		cosH0(:) = -sinDec(:) * SLAt / (cLat * sqrt(1.d0 - (sinDec(:)) ** 2.d0))

		day_len(:) = Acos(cosH0(:)) / Pi

	end function get_day_len_doy


	! Get relative degree days 
	! effective temperature sum based on monthly data, according to LandClim implementation
	function calc_rdd(lat, tmean, rad, lai, k, day_length) result(rdd)
		
		!input
		real(kind=8):: k, lat
		real(kind=8), dimension(12):: tmean, day_length, rad, lai

		!local
		integer:: diapause_month, i, ndays
		real(kind=8):: tbark, teff, topt, tmax, dt_l, dt_u, alpha, beta, &
					 gamma, teff_sum, month_share, par_k, rdd
		real(kind=8), dimension(12):: day_length_h

		! Parameters
		topt = 30.4d0
		tmax = 40.9958913d0
		dt_l = 8.3d0
		dt_u = 38.9d0
		alpha = 0.02876507d0
		beta = 3.5922336
		gamma = 1.24657367
		par_k = 557.d0

		teff_sum = 0.d0

		! Convert to hours
		day_length_h(:) = day_length(:) / 3600.d0

		! Get month when photoperiod > 14.5
		call get_diapause_month(diapause_month, month_share, lat, day_length_h)

		do i=1, diapause_month
			
			tbark = calc_tbark(tmean(i), rad(i), lai(i), k) 

			teff_sum = 0.d0

			
			! Get effective temperature
			if (tbark >= dt_l .and. tbark <= dt_u) then

				! If photoperiod > 14.5h within this month, correct the number of days used in the Teff sum
				if (i == diapause_month) then
					ndays = daysInMonth(i) * month_share
				else
					ndays = daysInMonth(i)
				end if

				! Get month Teff
				teff = ndays * (topt - dt_l) * exp(alpha * tbark) - exp(alpha * tmax - (tmax - tbark) / beta ) - gamma	
					
			else
				teff = 0.d0
			end if

			teff_sum = teff_sum + teff

		end do

		! Relative degree days
		rdd = teff_sum / par_k

	end function calc_rdd

	! Calculate bark temperature based on iLand implementation
	function calc_tbark(tmean, rad, lai, k) result(tbark)

		!input
		real(kind=8):: tmean, rad, lai, k

		!local
		real(kind=8):: tbark, rad_corr

		! Get radiation in Wh/m2/day, corrected by canopy interception
		rad_corr = rad / 0.0036d0 * exp(-k * lai)

		! Compute bark temperature based on empirical relation from iLand (Seidl & Rammer 2016)
		tbark = -0.173d0 + 0.0008518d0 * rad_corr + 1.054d0 * tmean

	end function calc_tbark


 ! Random number generator
    function gengam ( a, r ) result(rangam)
        
        real(kind=8):: a
        real(kind=8):: rangam
        real(kind=8):: r
     
        rangam = sgamma( r ) / a

    end function gengam

    function sgamma(a) result(sgamm)

        real(kind=8):: a
        real(kind=8), parameter :: a1 =  0.3333333E+00
        real(kind=8), parameter :: a2 = -0.2500030E+00
        real(kind=8), parameter :: a3 =  0.2000062E+00
        real(kind=8), parameter :: a4 = -0.1662921E+00
        real(kind=8), parameter :: a5 =  0.1423657E+00
        real(kind=8), parameter :: a6 = -0.1367177E+00
        real(kind=8), parameter :: a7 =  0.1233795E+00
        real(kind=8):: b
        real(kind=8):: c
        real(kind=8):: d
        real(kind=8):: e
        real(kind=8), parameter :: e1 = 1.0E+00
        real(kind=8), parameter :: e2 = 0.4999897E+00
        real(kind=8), parameter :: e3 = 0.1668290E+00
        real(kind=8), parameter :: e4 = 0.0407753E+00
        real(kind=8), parameter :: e5 = 0.0102930E+00
        real(kind=8):: p
        real(kind=8):: q
        real(kind=8):: q0
        real(kind=8), parameter :: q1 =  0.04166669E+00
        real(kind=8), parameter :: q2 =  0.02083148E+00
        real(kind=8), parameter :: q3 =  0.00801191E+00
        real(kind=8), parameter :: q4 =  0.00144121E+00
        real(kind=8), parameter :: q5 = -0.00007388E+00
        real(kind=8), parameter :: q6 =  0.00024511E+00
        real(kind=8), parameter :: q7 =  0.00024240E+00
        real(kind=8):: r
        real(kind=8):: r4_uni_01
        real(kind=8):: s
        real(kind=8):: s2
        !real(kind=8):: sexpo
        real(kind=8):: si
        real(kind=8):: sgamm
        !real(kind=8):: snorm
        real(kind=8), parameter :: sqrt32 = 5.656854E+00
        real(kind=8):: t
        real(kind=8):: u
        real(kind=8):: v
        real(kind=8):: w
        real(kind=8):: x

		if ( 1.0E+00 <= a ) then

			s2 = a - 0.5E+00
			s = sqrt ( s2 )
			d = sqrt32 - 12.0E+00 * s
		!
		!  Immediate acceptance.
		!
			t = snorm( )
			x = s + 0.5E+00 * t
			sgamm = x * x
		
			if ( 0.0E+00 <= t ) then
			  	return
			end if
		!
		!  Squeeze acceptance.
		!
			call random_number(u)
			if ( d * u <= t * t * t ) then
			  	return
			end if
		
			r = 1.0E+00 / a
			q0 = (((((( q7 &
			  * r + q6 ) &
			  * r + q5 ) &
			  * r + q4 ) &
			  * r + q3 ) &
			  * r + q2 ) &
			  * r + q1 ) &
			  * r
		!
		!  Approximation depending on size of parameter A.
		!
			if ( 13.022E+00 < a ) then
				b = 1.77E+00
				si = 0.75E+00
				c = 0.1515E+00 / s
			else if ( 3.686E+00 < a ) then
				b = 1.654E+00 + 0.0076E+00 * s2
				si = 1.68E+00 / s + 0.275E+00
				c = 0.062E+00 / s + 0.024E+00
			else
				b = 0.463E+00 + s + 0.178E+00 * s2
				si = 1.235E+00
				c = 0.195E+00 / s - 0.079E+00 + 0.16E+00 * s
			end if
		!
		!  Quotient test.
		!
			if ( 0.0E+00 < x ) then
		
				v = 0.5E+00 * t / s
			
				if ( 0.25E+00 < abs ( v ) ) then
					q = q0 - s * t + 0.25E+00 * t * t + 2.0E+00 * s2 * log ( 1.0E+00 + v )
				else
					q = q0 + 0.5E+00 * t * t * (((((( a7 &
					* v + a6 ) &
					* v + a5 ) &
					* v + a4 ) &
					* v + a3 ) &
					* v + a2 ) &
					* v + a1 ) &
					* v
				end if
			
				if ( log ( 1.0E+00 - u ) <= q ) then
					return
				end if
		
		end if
		
			do
		
				e = sexpo( )
				call random_number(u)
				u = 2.0E+00 * u - 1.0E+00
			
				if ( 0.0E+00 <= u ) then
						t = b + abs ( si * e )
				else
						t = b - abs ( si * e )
				end if
			!
			!  Possible rejection.
			!
				if ( t < -0.7187449E+00 ) then
						cycle
				end if
			!
			!  Calculate V and quotient Q.
			!
				v = 0.5E+00 * t / s
			
				if ( 0.25E+00 < abs ( v ) ) then
					q = q0 - s * t + 0.25E+00 * t * t + 2.0E+00 * s2 * log ( 1.0E+00 + v )
				else
					q = q0 + 0.5E+00 * t * t * (((((( a7 &
					* v + a6 ) &
					* v + a5 ) &
					* v + a4 ) &
					* v + a3 ) &
					* v + a2 ) &
					* v + a1 ) &
					*  v
				end if
		!
			!  Hat acceptance.
			!
				if ( q <= 0.0E+00 ) then
					cycle
				end if
			
				if ( 0.5E+00 < q ) then
					w = exp ( q ) - 1.0E+00
				else
					w = (((( e5 * q + e4 ) * q + e3 ) * q + e2 ) * q + e1 ) * q
				end if
		!
		!  May have to sample again.
		!
				if ( c * abs ( u ) <= w * exp ( e - 0.5E+00 * t * t ) ) then
					exit
				end if
		
			end do
		
			x = s + 0.5E+00 * t
			sgamm = x * x
		
			return
		!
		!  Method for A < 1.
		!
		  else
		
			b = 1.0E+00 + 0.3678794E+00 * a
		
			do
				e = sexpo()
				call random_number(u)
				p = b * u
			
				if ( p < 1.0E+00 ) then
			
					sgamm = exp ( log ( p ) / a )
			
					if ( sgamm <= e ) then
						return
					end if
			
					cycle
			
				end if
			
				sgamm = - log ( ( b - p ) / a )
			
				if ( ( 1.0E+00 - a ) * log ( sgamm ) <= e ) then
					exit
				end if
		
			end do
		
		  end if
		
		return
		                                                 

    end function sgamma


    function sexpo( ) result(sexp)

        
          real(kind=8):: a
          integer:: i
          real(kind=8):: q(8)
          real(kind=8):: runi
          real(kind=8):: sexp
          real(kind=8):: u
          real(kind=8):: umin
          real(kind=8):: ustar
        
          save q
        
          data q / &
               0.6931472E+00, &
               0.9333737E+00, &
               0.9888778E+00, &
               0.9984959E+00, &
               0.9998293E+00, &
               0.9999833E+00, &
               0.9999986E+00, &
               0.9999999E+00 /
        
          a = 0.0E+00
          call random_number(u)
        
          do
        
            u = u + u
        
            if ( 1.0E+00 < u ) then
              exit
            end if
        
            a = a + q(1)
        
          end do
        
          u = u - 1.0E+00
        
          if ( u <= q(1) ) then
            sexp = a + u
          end if
        
          i = 1
          call random_number(ustar)
          umin = ustar
        
          do
        
            call random_number(ustar)
            umin = min ( umin, ustar )
            i = i + 1
        
            if ( u <= q(i) ) then
              exit
            end if
        
          end do
        
          sexp = a + umin * q(1)
        
          
        end function sexpo


        function snorm( ) result(snor)

            !*****************************************************************************80
            !
            !! SNORM samples the standard normal distribution.
            !
            !  Discussion:
            !
            !    This procedure corresponds to algorithm FL, with M = 5, in the reference.
            !
            !  Licensing:
            !
            !    This code is distributed under the GNU LGPL license.
            !
            !  Modified:
            !
            !    31 March 2013
            !
            !  Author:
            !
            !    Original FORTRAN77 version by Barry Brown, James Lovato.
            !    FORTRAN90 version by John Burkardt.
            !
            !  Reference:
            !
            !    Joachim Ahrens, Ulrich Dieter,
            !    Extensions of Forsythe's Method for Random
            !    Sampling from the Normal Distribution,
            !    Mathematics of Computation,
            !    Volume 27, Number 124, October 1973, page 927-937.
            !
            !  Parameters:
            !
            !    Output, real SNORM, a random deviate from the distribution.
            !
              implicit none
            
              real(kind=8):: a(32)
              real(kind=8):: aa
              real(kind=8):: d(31)
              real(kind=8):: h(31)
              integer:: i
              real(kind=8):: runi
              real(kind=8):: s
              real(kind=8):: snor
              real(kind=8):: t(31)
              real(kind=8):: tt
              real(kind=8):: u
              real(kind=8):: ustar
              real(kind=8):: w
              real(kind=8):: y
            
              data a / &
                    0.0000000E+00, 0.3917609E-01, 0.7841241E-01, 0.1177699E+00, &
                    0.1573107E+00, 0.1970991E+00, 0.2372021E+00, 0.2776904E+00, &
                    0.3186394E+00, 0.3601299E+00, 0.4022501E+00, 0.4450965E+00, &
                    0.4887764E+00, 0.5334097E+00, 0.5791322E+00, 0.6260990E+00, &
                    0.6744898E+00, 0.7245144E+00, 0.7764218E+00, 0.8305109E+00, &
                    0.8871466E+00, 0.9467818E+00, 1.009990E+00,  1.077516E+00, &
                    1.150349E+00,  1.229859E+00,  1.318011E+00,  1.417797E+00, &
                    1.534121E+00,  1.675940E+00,  1.862732E+00,  2.153875E+00 /
            
              data d / &
                    0.0000000E+00, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00, &
                    0.0000000E+00, 0.2636843E+00, 0.2425085E+00, 0.2255674E+00, &
                    0.2116342E+00, 0.1999243E+00, 0.1899108E+00, 0.1812252E+00, &
                    0.1736014E+00, 0.1668419E+00, 0.1607967E+00, 0.1553497E+00, &
                    0.1504094E+00, 0.1459026E+00, 0.1417700E+00, 0.1379632E+00, &
                    0.1344418E+00, 0.1311722E+00, 0.1281260E+00, 0.1252791E+00, &
                    0.1226109E+00, 0.1201036E+00, 0.1177417E+00, 0.1155119E+00, &
                    0.1134023E+00, 0.1114027E+00, 0.1095039E+00 /
            
              data h / &
                    0.3920617E-01, 0.3932705E-01, 0.3950999E-01, 0.3975703E-01, &
                    0.4007093E-01, 0.4045533E-01, 0.4091481E-01, 0.4145507E-01, &
                    0.4208311E-01, 0.4280748E-01, 0.4363863E-01, 0.4458932E-01, &
                    0.4567523E-01, 0.4691571E-01, 0.4833487E-01, 0.4996298E-01, &
                    0.5183859E-01, 0.5401138E-01, 0.5654656E-01, 0.5953130E-01, &
                    0.6308489E-01, 0.6737503E-01, 0.7264544E-01, 0.7926471E-01, &
                    0.8781922E-01, 0.9930398E-01, 0.1155599E+00, 0.1404344E+00, &
                    0.1836142E+00, 0.2790016E+00, 0.7010474E+00 /
            
              data t / &
                    0.7673828E-03, 0.2306870E-02, 0.3860618E-02, 0.5438454E-02, &
                    0.7050699E-02, 0.8708396E-02, 0.1042357E-01, 0.1220953E-01, &
                    0.1408125E-01, 0.1605579E-01, 0.1815290E-01, 0.2039573E-01, &
                    0.2281177E-01, 0.2543407E-01, 0.2830296E-01, 0.3146822E-01, &
                    0.3499233E-01, 0.3895483E-01, 0.4345878E-01, 0.4864035E-01, &
                    0.5468334E-01, 0.6184222E-01, 0.7047983E-01, 0.8113195E-01, &
                    0.9462444E-01, 0.1123001E+00, 0.1364980E+00, 0.1716886E+00, &
                    0.2276241E+00, 0.3304980E+00, 0.5847031E+00 /
            
              call random_number(u)
              if ( u <= 0.5E+00 ) then
                s = 0.0E+00
              else
                s = 1.0E+00
              end if
              u = 2.0E+00 * u - s
              u = 32.0E+00 * u
              i = int ( u )
              if ( i == 32 ) then
                i = 31
              end if
            !
            !  Center
            !
              if ( i /= 0 ) then
            
                ustar = u - real ( i )
                aa = a(i)
            
                do
            
                  if ( t(i) < ustar ) then
            
                    w = ( ustar - t(i) ) * h(i)
            
                    y = aa + w
            
                    if ( s /= 1.0E+00 ) then
                      snor = y
                    else
                      snor = -y
                    end if
            
                    return
            
                  end if
            
                  call random_number(u)
                  w = u * ( a(i+1) - aa )
                  tt = ( 0.5E+00 * w + aa ) * w
            
                  do
            
                    if ( tt < ustar ) then
                      y = aa + w
                      if ( s /= 1.0E+00 ) then
                        snor = y
                      else
                        snor = -y
                      end if
                      return
                    end if
            
                    call random_number(u)
            
                    if ( ustar < u ) then
                      exit
                    end if
            
                    tt = u
                    call random_number(ustar)
            
                  end do
            
                  call random_number(ustar)
            
                end do
            !
            !  Tail
            !
              else
            
                i = 6
                aa = a(32)
            
                do
            
                  u = u + u
            
                  if ( 1.0E+00 <= u ) then
                    exit
                  end if
            
                  aa = aa + d(i)
                  i = i + 1
            
                end do
            
                u = u - 1.0E+00
                w = u * d(i)
                tt = ( 0.5E+00 * w + aa ) * w
            
                do
            
                    call random_number(ustar)
            
                  if ( tt < ustar ) then
                    y = aa + w
                    if ( s /= 1.0E+00 ) then
                      snor = y
                    else
                      snor = -y
                    end if
                    return
                  end if
            
                  call random_number(u)
            
                  if ( u <= ustar ) then
                    tt = u
                  else
                    call random_number(u)
                    w = u * d(i)
                    tt = ( 0.5E+00 * w + aa ) * w
                  end if
            
                end do
            
              end if
            
            end function snorm

			subroutine random_stdgamma_alpha_ge_1(alpha,x)
				real(kind=8),intent(in) :: alpha
				real(kind=8),intent(out) :: x
				real(kind=8) :: y,t,u1,u2
				do
				   call random_stduniform(u1)
				   y = -log(u1)
				   t = (y/exp(y-1))**(alpha-1)
				   call random_stduniform(u2)
				   if(u2 <= t) then
					  x = alpha*y
					  exit
				   end if
				end do
			 end subroutine random_stdgamma_alpha_ge_1

			 subroutine random_stdgamma(alpha,x)
				real(kind=8),intent(in) :: alpha
				real(kind=8),intent(out) :: x
				real(kind=8) :: g,u
				if(alpha<=0) then
				   stop "alpha<=0"
				else if(alpha<1) then
				   call random_stdgamma_alpha_ge_1(alpha+1.0,g)
				   call random_stduniform(u)
				   x = g*u**(1.0/alpha)
				else
				   call  random_stdgamma_alpha_ge_1(alpha,x)
				end if
			 end subroutine random_stdgamma

			 subroutine random_gamma(alpha,beta,x)
				real(kind=8),intent(in) :: alpha,beta
				real(kind=8),intent(out) :: x
				call random_stdgamma(alpha,x)
				x = x*beta
			 end subroutine random_gamma

			 subroutine random_stduniform(u)
				implicit none
				real(kind=8),intent(out) :: u
				real(kind=8) :: r
				call random_number(r)
				u = 1 - r
			 end subroutine random_stduniform


end module beetle_subroutines


