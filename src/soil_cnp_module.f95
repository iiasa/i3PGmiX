module soil_cnp

use soil_cnp_subroutines

IMPLICIT NONE
public :: run_icbm, mod5c, mod5c_year, compAWENH, yasso_outflux, yasso_outflux_year, &
            yasso_nav

contains

    !Subroutines for ICBM/2N--------------------------------------------------------------
    subroutine run_icbm(krmax, klmax, komax, Yl_C, Yr_C, O_C, Yl_N, Yr_N, O_N, Litter_C, &
                        hc, el, er, qbc, qh, qil, qir, Yl_input, Yr_input, O_Coutflux, &
                        Yl_Coutflux, Yr_Coutflux, Nav, TotalCarbo, TotalNitro, fol_litter, &
                        temp, pr, excessSW, asw, aswmax, soil_class, n_depo, dmC, n_sp)
        implicit none		 
        integer, intent(in):: n_sp, soil_class
        integer :: i
        real(kind=8), dimension(n_sp), intent(in):: krmax, klmax, qil, qir, fol_litter, &
                                                    hc, el, er, qbc, qh, Yl_input, Yr_input
        real(kind=8), intent(out):: O_C, O_N, Nav, TotalCarbo, TotalNitro, O_Coutflux	
        real(kind=8), dimension(n_sp), intent(out):: Yl_C, Yr_C, Yl_N, Yr_N, Litter_C, &
                                                    Yl_Coutflux, Yr_Coutflux											
        real(kind=8), intent(in):: komax, dmC, temp, pr, n_depo, excessSW, asw, aswmax						
        real(kind=8), dimension(n_sp) :: Litter_Coutflux, humification_Litter, &
                                        humification_l, humification_r, Yl_Noutflux, humification_N_l, &
                                        humification_N_r, Yr_Noutflux
        real(kind=8) :: leaching_org, leaching_min, water_ratio, kl, kr, ko, O_Noutflux								
        
        water_ratio = asw / aswmax
        
        do i=1, n_sp
            !Adjust decomposition rates based on iLand implementation (Seidl et al. 2012)
            kr = krmax(i) * (1.d0/(1.d0+30.d0*exp(-8.5d0*(water_ratio)))) * exp(5.508032845d0 - 308.56d0/(273.15d0+temp-227.13d0)) 
            kl = klmax(i) * (1.d0/(1.d0+30.d0*exp(-8.5d0*(water_ratio)))) * exp(5.508032845d0 - 308.56d0/(273.15d0+temp-227.13d0)) 

            !Compute monthly out fluxes
            Yl_Coutflux(i) = kl * (1.d0 - hc(i)) * Yl_C(i)  
            Yr_Coutflux(i) = kr * (1.d0 - hc(i)) * Yr_C(i) 
            Litter_Coutflux(i) = kl * (1.d0 - hc(i)) * Litter_C(i) 
            humification_Litter(i) = kl * hc(i) * Litter_C(i)
            humification_l(i) = kl * hc(i) * Yl_C(i)
            humification_r(i) = kr * hc(i) * Yr_C(i)
            
            !Compute N out fluxes
            Yl_Noutflux(i) = kl * ((1.d0 - hc(i)) / (1.d0 - el(i))) * (Yl_N(i) - el(i) * (Yl_C(i) / qbc(i)))
            humification_N_l(i) = kl * hc(i) * (Yl_N(i) / qh(i))
            humification_N_r(i) = kr * hc(i) * (Yr_N(i) / qh(i))
            Yr_Noutflux(i) = kr * ((1.d0 - hc(i)) / (1.d0 - er(i))) * (Yr_N(i) - er(i) * (Yr_C(i) / qbc(i)))

            !Now calculate the end-of-month carbon and nitrogen pools				
            Yr_C(i) = Yr_C(i) + Yr_input(i) - Yr_Coutflux(i) - humification_r(i)  
            Yl_C(i) = Yl_C(i) + Yl_input(i) - Yl_Coutflux(i) - humification_l(i) 
            Litter_C(i) = Litter_C(i) + fol_litter(i) - Litter_Coutflux(i) - humification_Litter(i)

            Yr_N(i) = Yr_N(i) + ((Yr_input(i)) / (2.d0 * qir(i))) - Yr_Noutflux(i) - humification_N_r(i)
            Yl_N(i) = Yl_N(i) + ((Yl_input(i)) / (2.d0 * qil(i))) - Yl_Noutflux(i) - humification_N_l(i)
        end do

        !Repeat for the old carbon and nitrogen pools
        ko = komax * (1.d0/(1.d0+30.d0*exp(-8.5d0*(water_ratio)))) * exp(5.508032845d0 - 308.56d0/(273.15d0+temp-227.13d0)) 

        O_Coutflux = ko * O_C  
        O_Noutflux = ko * O_N

        O_C = O_C + sum(humification_l(:)) + sum(humification_r(:)) - O_Coutflux
        O_N = O_N + sum(humification_N_r(:)) + sum(humification_N_l(:)) - O_Noutflux

        TotalCarbo = sum(Yr_C(:) + Yl_C(:)) + O_C
        TotalNitro = sum(Yr_N(:) + Yl_N(:)) + O_N

        !Calculate N leaching
        call compute_leaching( soil_class, leaching_org, leaching_min, excessSW, asw )
        
        !Get available soil N
        Nav = sum(Yr_Noutflux(:)) * (1.d0-leaching_org)  + sum(+ Yl_Noutflux(:)) * (1.d0-leaching_org) + &
                        O_Noutflux * (1.d0 - leaching_min)  + n_depo


    end subroutine run_icbm
    ! End subroutines for ICBM/2N---------------------------------------------------------

   


    !Subroutines for Yasso ---------------------------------------------------------------

    SUBROUTINE mod5c(theta,time,temp,prec,init,b,d,leac,xt,steadystate_pred)
    IMPLICIT NONE
        !********************************************* &
        ! GENERAL DESCRIPTION FOR ALL THE MEASUREMENTS
        !********************************************* &
        ! returns the model prediction xt for the given parameters
        ! 1-16 matrix A entries: 4*alpha, 12*p

        ! 17-21 Leaching parameters: w1,...,w5 IGNORED IN THIS FUNCTION

        ! 22-23 Temperature-dependence parameters for AWE fractions: beta_1, beta_2

        ! 24-25 Temperature-dependence parameters for N fraction: beta_N1, beta_N2

        ! 26-27 Temperature-dependence parameters for H fraction: beta_H1, beta_H2

        ! 28-30 Precipitation-dependence parameters for AWE, N and H fraction: gamma, gamma_N, gamma_H

        ! 31-32 Humus decomposition parameters: p_H, alpha_H (Note the order!)

        ! 33-35 Woody parameters: theta_1, theta_2, r

        REAL (kind=8),DIMENSION(35),INTENT(IN) :: theta ! parameters
        REAL (kind=8),INTENT(IN) :: time,d,leac ! time,size,leaching
        !REAL (kind=8),DIMENSION(12),INTENT(IN) :: temp ! monthly mean temperatures
        REAL (kind=8),INTENT(IN) :: temp ! monthly mean temperatures
        REAL (kind=8),INTENT(IN) :: prec ! annual/monthly precipitation
        REAL (kind=8),DIMENSION(5),INTENT(IN) :: init ! initial state
        REAL (kind=8),DIMENSION(5),INTENT(IN) :: b ! infall
        REAL (kind=8),DIMENSION(5),INTENT(OUT) :: xt ! the result i.e. x(t)
        INTEGER, INTENT(IN) :: steadystate_pred
        ! LOGICAL,OPTIONAL,INTENT(IN) :: steadystate_pred ! set to true if ignore 'time' and compute solution
        ! in steady-state conditions (which sould give equal solution as if time is set large enough)
        REAL (kind=8),DIMENSION(5,5) :: A,At,mexpAt
        INTEGER :: i
        REAL (kind=8),PARAMETER :: pi = 3.141592653589793d0
        REAL (kind=8) :: tem,temN,temH,size_dep
        REAL (kind=8),DIMENSION(5) :: z1,z2
        REAL (kind=8),PARAMETER :: tol = 10E-12
        LOGICAL :: ss_pred

        ! IF(PRESENT(steadystate_pred)) THEN
            ! ss_pred = steadystate_pred
        ! ENDIF
        IF(steadystate_pred == 1) THEN
            ss_pred = .TRUE.
        ELSE
            ss_pred = .FALSE.
        ENDIF

        !#########################################################################
        ! Compute the coefficient matrix A for the differential equation

        ! temperature annual cycle approximation
        ! te(1) = climate(1)+4*climate(3)*(1/SQRT(2.0)-1)/pi
        ! te(2) = climate(1)-4*climate(3)/SQRT(2.0)/pi
        ! te(3) = climate(1)+4*climate(3)*(1-1/SQRT(2.0))/pi
        ! te(4) = climate(1)+4*climate(3)/SQRT(2.0)/pi

        ! DO i = 1,4 ! Average temperature dependence
        !     tem = tem+EXP(theta(22)*te(i)+theta(23)*te(i)**2.0)/4.0 ! Gaussian
        !     temN = temN+EXP(theta(24)*te(i)+theta(25)*te(i)**2.0)/4.0
        !     temH = temH+EXP(theta(26)*te(i)+theta(27)*te(i)**2.0)/4.0
        ! END DO

        ! Set up climate dependence factors
        tem = 0.d0
        temN = 0.d0
        temH = 0.d0

        ! Monthly temperature dependence
        tem = tem+EXP(theta(22)*temp+theta(23)*temp**2.d0)
        temN = temN+EXP(theta(24)*temp+theta(25)*temp**2.d0)
        temH = temH+EXP(theta(26)*temp+theta(27)*temp**2.d0)
        
        ! Precipitation dependence
        tem = tem*(1.d0-EXP(theta(28)*prec/83.3333d0))
        temN = temN*(1.d0-EXP(theta(29)*prec/83.3333d0))
        temH = temH*(1.d0-EXP(theta(30)*prec/83.3333d0))
        
        ! Monthly temperature dependence
        !DO i = 1,12
        !    tem = tem+EXP(theta(22)*temp(i)+theta(23)*temp(i)**2.d0)
        !    temN = temN+EXP(theta(24)*temp(i)+theta(25)*temp(i)**2.d0)
        !    temH = temH+EXP(theta(26)*temp(i)+theta(27)*temp(i)**2.d0)
        !END DO
        
        !tem = tem*(1.d0-EXP(theta(28)*prec/1000.d0))/12
        !temN = temN*(1.d0-EXP(theta(29)*prec/1000.d0))/12
        !temH = temH*(1.d0-EXP(theta(30)*prec/1000.d0))/12

        ! Size class dependence -- no effect if d == 0.0
        size_dep = MIN(1.d0,(1.d0+theta(33)*d+theta(34)*d**2.d0)**(-ABS(theta(35))))

        ! check rare case where no decomposition happens for some compartments
        ! (basically, if no rain)
        IF (tem <= tol) THEN
            xt = init + b*time
            return
        END IF

        ! Calculating matrix A (will work ok despite the sign of alphas)
        DO i = 1,3
            A(i,i) = -ABS(theta(i))*tem*size_dep
        END DO
        A(4,4) = -ABS(theta(4))*temN*size_dep

        A(1,2) = theta(5)*ABS(A(2,2))
        A(1,3) = theta(6)*ABS(A(3,3))
        A(1,4) = theta(7)*ABS(A(4,4))
        A(1,5) = 0.0 ! no mass flows from H -> AWEN
        A(2,1) = theta(8)*ABS(A(1,1))
        A(2,3) = theta(9)*ABS(A(3,3))
        A(2,4) = theta(10)*ABS(A(4,4))
        A(2,5) = 0.0
        A(3,1) = theta(11)*ABS(A(1,1))
        A(3,2) = theta(12)*ABS(A(2,2))
        A(3,4) = theta(13)*ABS(A(4,4))
        A(3,5) = 0.0
        A(4,1) = theta(14)*ABS(A(1,1))
        A(4,2) = theta(15)*ABS(A(2,2))
        A(4,3) = theta(16)*ABS(A(3,3))
        A(4,5) = 0.0
        A(5,5) = -ABS(theta(32))*temH ! no size effect in humus
        DO i = 1,4
            A(5,i) = theta(31)*ABS(A(i,i)) ! mass flows AWEN -> H (size effect is present here)
        END DO

        ! Leaching (no leaching for humus)
        DO i = 1,4
            ! A(i,i) = A(i,i)+leac*climate(2)/1000.0
            A(i,i) = A(i,i)+leac*prec/1000.d0
        END DO

        !#########################################################################
        ! Solve the differential equation x'(t) = A(theta)*x(t) + b, x(0) = init

        IF(ss_pred) THEN
            ! Solve DE directly in steady state conditions (time = infinity)
            ! using the formula 0 = x'(t) = A*x + b => x = -A**-1*b
            CALL solve(-A, b, xt)
        ELSE
            ! Solve DE in given time
            z1 = MATMUL(A,init) + b * 12.d0 ! correct to the full monthly litter input
            At = A/12.d0 !A*time !At = A*t
            CALL matrixexp(At,mexpAt)
            z2 = MATMUL(mexpAt,z1) - b * 12.d0 ! correct to the full monthly litter input
            CALL solve(A,z2,xt) ! now it can be assumed A is non-singular
        ENDIF

    END SUBROUTINE mod5c


    SUBROUTINE mod5c_year(theta,time,temp,prec,init,b,d,leac,xt,steadystate_pred)
        IMPLICIT NONE
        !********************************************* &
        ! GENERAL DESCRIPTION FOR ALL THE MEASUREMENTS
        !********************************************* &
        ! returns the model prediction xt for the given parameters
        ! 1-16 matrix A entries: 4*alpha, 12*p

        ! 17-21 Leaching parameters: w1,...,w5 IGNORED IN THIS FUNCTION

        ! 22-23 Temperature-dependence parameters for AWE fractions: beta_1, beta_2

        ! 24-25 Temperature-dependence parameters for N fraction: beta_N1, beta_N2

        ! 26-27 Temperature-dependence parameters for H fraction: beta_H1, beta_H2

        ! 28-30 Precipitation-dependence parameters for AWE, N and H fraction: gamma, gamma_N, gamma_H

        ! 31-32 Humus decomposition parameters: p_H, alpha_H (Note the order!)

        ! 33-35 Woody parameters: theta_1, theta_2, r

        REAL (kind=8),DIMENSION(35),INTENT(IN) :: theta ! parameters
        REAL (kind=8),INTENT(IN) :: time,d,leac ! time,size,leaching
        REAL (kind=8),DIMENSION(12),INTENT(IN) :: temp ! monthly mean temperatures
        REAL (kind=8),INTENT(IN) :: prec ! annual precipitation
        REAL (kind=8),DIMENSION(5),INTENT(IN) :: init ! initial state
        REAL (kind=8),DIMENSION(5),INTENT(IN) :: b ! infall
        REAL (kind=8),DIMENSION(5),INTENT(OUT) :: xt ! the result i.e. x(t)
        INTEGER, INTENT(IN) :: steadystate_pred
        ! LOGICAL,OPTIONAL,INTENT(IN) :: steadystate_pred ! set to true if ignore 'time' and compute solution
        ! in steady-state conditions (which sould give equal solution as if time is set large enough)
        REAL (kind=8),DIMENSION(5,5) :: A,At,mexpAt
        INTEGER :: i
        REAL (kind=8),PARAMETER :: pi = 3.141592653589793
        REAL (kind=8) :: tem,temN,temH,size_dep
        REAL (kind=8),DIMENSION(5) :: z1,z2
        REAL (kind=8),PARAMETER :: tol = 1E-12
        LOGICAL :: ss_pred

        ! IF(PRESENT(steadystate_pred)) THEN
            ! ss_pred = steadystate_pred
        ! ENDIF
        IF(steadystate_pred == 1) THEN
            ss_pred = .TRUE.
        ELSE
            ss_pred = .FALSE.
        ENDIF

        !#########################################################################
        ! Compute the coefficient matrix A for the differential equation

        ! temperature annual cycle approximation
        ! te(1) = climate(1)+4*climate(3)*(1/SQRT(2.0)-1)/pi
        ! te(2) = climate(1)-4*climate(3)/SQRT(2.0)/pi
        ! te(3) = climate(1)+4*climate(3)*(1-1/SQRT(2.0))/pi
        ! te(4) = climate(1)+4*climate(3)/SQRT(2.0)/pi

        ! DO i = 1,4 ! Average temperature dependence
        !     tem = tem+EXP(theta(22)*te(i)+theta(23)*te(i)**2.0)/4.0 ! Gaussian
        !     temN = temN+EXP(theta(24)*te(i)+theta(25)*te(i)**2.0)/4.0
        !     temH = temH+EXP(theta(26)*te(i)+theta(27)*te(i)**2.0)/4.0
        ! END DO

        ! Set up climate dependence factors
        tem = 0.0
        temN = 0.0
        temH = 0.0

        ! Monthly temperature dependence
        DO i = 1,12
            tem = tem+EXP(theta(22)*temp(i)+theta(23)*temp(i)**2.0)
            temN = temN+EXP(theta(24)*temp(i)+theta(25)*temp(i)**2.0)
            temH = temH+EXP(theta(26)*temp(i)+theta(27)*temp(i)**2.0)
         END DO

        ! Precipitation dependence
        tem = tem*(1.0-EXP(theta(28)*prec/1000.0))/12
        temN = temN*(1.0-EXP(theta(29)*prec/1000.0))/12
        temH = temH*(1.0-EXP(theta(30)*prec/1000.0))/12

        ! Size class dependence -- no effect if d == 0.0
        size_dep = MIN(1.0,(1.0+theta(33)*d+theta(34)*d**2.0)**(-ABS(theta(35))))

        ! check rare case where no decomposition happens for some compartments
        ! (basically, if no rain)
        IF (tem <= tol) THEN
            xt = init + b*time
            return
        END IF

        ! Calculating matrix A (will work ok despite the sign of alphas)
        DO i = 1,3
            A(i,i) = -ABS(theta(i))*tem*size_dep
        END DO
        A(4,4) = -ABS(theta(4))*temN*size_dep

        A(1,2) = theta(5)*ABS(A(2,2))
        A(1,3) = theta(6)*ABS(A(3,3))
        A(1,4) = theta(7)*ABS(A(4,4))
        A(1,5) = 0.0 ! no mass flows from H -> AWEN
        A(2,1) = theta(8)*ABS(A(1,1))
        A(2,3) = theta(9)*ABS(A(3,3))
        A(2,4) = theta(10)*ABS(A(4,4))
        A(2,5) = 0.0
        A(3,1) = theta(11)*ABS(A(1,1))
        A(3,2) = theta(12)*ABS(A(2,2))
        A(3,4) = theta(13)*ABS(A(4,4))
        A(3,5) = 0.0
        A(4,1) = theta(14)*ABS(A(1,1))
        A(4,2) = theta(15)*ABS(A(2,2))
        A(4,3) = theta(16)*ABS(A(3,3))
        A(4,5) = 0.0
        A(5,5) = -ABS(theta(32))*temH ! no size effect in humus
        DO i = 1,4
            A(5,i) = theta(31)*ABS(A(i,i)) ! mass flows AWEN -> H (size effect is present here)
        END DO

        ! Leaching (no leaching for humus)
        DO i = 1,4
            ! A(i,i) = A(i,i)+leac*climate(2)/1000.0
            A(i,i) = A(i,i)+leac*prec/1000.0
        END DO

        !#########################################################################
        ! Solve the differential equation x'(t) = A(theta)*x(t) + b, x(0) = init

        IF(ss_pred) THEN
            ! Solve DE directly in steady state conditions (time = infinity)
            ! using the formula 0 = x'(t) = A*x + b => x = -A**-1*b
            CALL solve(-A, b, xt)
        ELSE
            ! Solve DE in given time
            z1 = MATMUL(A,init) + b
            At = A*time !At = A*t
            CALL matrixexp(At,mexpAt)
            z2 = MATMUL(mexpAt,z1) - b
            CALL solve(A,z2,xt) ! now it can be assumed A is non-singular
        ENDIF                                   
    
    END SUBROUTINE mod5c_year



                                                        



    !#########################################################################
    ! Functions for solving the diff. equation, adapted for the Yasso case
    SUBROUTINE matrixexp(A,B)
        IMPLICIT NONE
        ! Approximated matrix exponential using Taylor series with scaling & squaring
        ! Accurate enough for the Yasso case
        INTEGER,PARAMETER :: n = 5
        REAL (kind=8),DIMENSION(n,n),INTENT(IN) :: A
        REAL (kind=8),DIMENSION(n,n),INTENT(OUT) :: B
        REAL (kind=8),DIMENSION(n,n) :: C,D
        REAL (kind=8) :: p,normiter
        INTEGER :: i,q,j
        q = 10 ! #terms in Taylor
        B = 0.d0
        DO i = 1,n
            B(i,i) = 1.d0
        END DO
        normiter = 2.d0 ! Amount of scaling & squaring
        j = 1
        CALL matrixnorm(A, p)
        DO
            IF (p<normiter) THEN
                EXIT
            END IF
            normiter = normiter*2.d0
            j = j+1
        END DO
        !write(*,*) normiter
        C = A/normiter ! scale
        B = B+C
        D = C
        DO i = 2,q ! compute Taylor expansion
            D = MATMUL(C,D)/REAL(i)
            B = B+D
        END DO
        DO i = 1,j ! square
            B = MATMUL(B,B)
        END DO
    END SUBROUTINE matrixexp

    SUBROUTINE matrixnorm(A,B)
        !returns elementwise (i.e. Frobenius) norm of a square matrix
        IMPLICIT NONE
        INTEGER,PARAMETER :: n = 5
        REAL (kind=8),DIMENSION(n,n),INTENT(IN) :: A
        REAL (kind=8),INTENT(OUT) :: b
        INTEGER :: i
        b = 0.d0
        DO i = 1,n
            b = b+SUM(A(:,i)**2.d0)
        END DO
        b = SQRT(b)
    END SUBROUTINE matrixnorm


    SUBROUTINE solve(A, b, x)
        ! Solve linear system A*x = b
        IMPLICIT NONE
        INTEGER,PARAMETER :: n = 5
        REAL (kind=8),DIMENSION(n,n),INTENT(IN) :: A
        REAL (kind=8),DIMENSION(n),INTENT(IN) :: b
        REAL (kind=8),DIMENSION(n),INTENT(OUT) :: x
        REAL (kind=8),DIMENSION(n,n) :: U
        REAL (kind=8),DIMENSION(n) :: c
        INTEGER :: i

        ! transform the problem to upper diagonal form
        CALL pgauss(A, b, U, c)

        ! solve U*x = c via back substitution
        x(n) = c(n)/U(n,n)
        DO i = n-1,1,-1
            x(i) = (c(i) - DOT_PRODUCT(U(i,i+1:n),x(i+1:n)))/U(i,i)
        END DO
    END SUBROUTINE solve

    SUBROUTINE pgauss(A, b, U, c)
        ! Transform the lin. system to upper diagonal form using gaussian elimination
        ! with pivoting
        IMPLICIT NONE
        INTEGER,PARAMETER :: n = 5
        REAL (kind=8),DIMENSION(n,n),INTENT(IN) :: A
        REAL (kind=8),DIMENSION(n),INTENT(IN) :: b
        REAL (kind=8),DIMENSION(n,n),INTENT(OUT) :: U
        REAL (kind=8),DIMENSION(n),INTENT(OUT) :: c
        INTEGER :: k, j
        REAL,PARAMETER :: tol = 1E-12

        U = A
        c = b
        DO k = 1,n-1
            CALL pivot(U,c,k) ! do pivoting (though may not be necessary in our case)
            IF (ABS(U(k,k)) <= tol) THEN
                write(*,*) 'Warning!!! Matrix is singular to working precision!'
            END IF
            U(k+1:n,k) = U(k+1:n,k)/U(k,k)
            DO j = k+1,n
                U(j,k+1:n) = U(j,k+1:n) - U(j,k)*U(k,k+1:n)
            END DO
            c(k+1:n) = c(k+1:n) - c(k)*U(k+1:n,k)
        END DO
    END SUBROUTINE pgauss

    SUBROUTINE pivot(A, b, k)
        ! perform pivoting to matrix A and vector b at row k
        IMPLICIT NONE
        INTEGER,PARAMETER :: n = 5
        REAL (kind=8),DIMENSION(n,n),INTENT(INOUT) :: A
        REAL (kind=8),DIMENSION(n),INTENT(INOUT) :: b
        INTEGER,INTENT(IN) :: k
        INTEGER :: q, pk

        !write(*,*) 'Pivot elements are: ', A(k:n,k)
        q = MAXLOC(ABS(A(k:n,k)),1)
        !write(*,*) q
        IF (q > 1) THEN
            pk = k-1+q
            A(k:pk:pk-k,:) = A(pk:k:k-pk,:)
            b(k:pk:pk-k) = b(pk:k:k-pk)
        END IF
        !write(*,*) 'Pivot elements are: ', A(k:n,k)
    END SUBROUTINE pivot

    SUBROUTINE compAWENH(Lit,AWENH,parsAWEN)
        IMPLICIT NONE
        REAL (kind=8),DIMENSION(5),INTENT(OUT) :: AWENH
        REAL (kind=8),DIMENSION(5),INTENT(IN) :: parsAWEN
        REAL (kind=8),INTENT(IN) :: Lit
        AWENH(1) = parsAWEN(1)*Lit
        AWENH(2) = parsAWEN(2)*Lit
        AWENH(3) = parsAWEN(3)*Lit
        AWENH(4) = parsAWEN(4)*Lit
        AWENH(5) = 0.d0
    END SUBROUTINE compAWENH

    ! Get stock changes from Yasso
    subroutine yasso_outflux (soilC, Yr_C, Yl_C, O_C, Yr_input, Yb_input, Yl_input, Yl_Coutflux, &
                            Yr_Coutflux, O_Coutflux, n_sp, nm, ii)

        integer, intent(in):: nm, n_sp, ii
        integer :: i
        real(kind=8), dimension(nm,n_sp,4,5), intent(inout) :: soilC
        real(kind=8), dimension(n_sp), intent(in) :: Yr_input, Yb_input, Yl_input
        real(kind=8), dimension(n_sp), intent(out) :: Yl_Coutflux, Yr_Coutflux, Yr_C, Yl_C
        real(kind=8), intent(out) :: O_C, O_Coutflux

        do i =1, n_sp

            Yl_Coutflux(i) = sum(soilC(ii-1,i,1,1:4)) - sum(soilC(ii,i,1,1:4)) + Yl_input(i)  
            Yr_Coutflux(i) = sum(soilC(ii-1,i,2:3,1:4)) - sum(soilC(ii,i,2:3,1:4)) + Yr_input(i) + Yb_input(i)
            
            Yr_C(i) = sum(soilC(ii,i,2:3,1:4))
            Yl_C(i) = sum(soilC(ii,i,1,1:4))

        end do

        O_Coutflux = soilC(ii-1,1,4,5) - soilC(ii,1,4,5)

        soilC(ii,1,4,5) = soilC(ii,1,4,5) + sum(soilC(ii,:,1:3,5))
        soilC(ii,:,1:3,5) = 0.d0

        O_C = soilC(ii,1,4,5)

    end subroutine  yasso_outflux 

        ! Get stock changes from Yasso at yearly time step
    subroutine yasso_outflux_year (soilC, Yr_C, Yl_C, O_C, Yr_input, Yb_input, Yl_input, Yl_Coutflux, &
                                    Yr_Coutflux, O_Coutflux, n_sp, nm, ii, step_aux)

        integer, intent(in):: nm, n_sp, ii, step_aux
        integer :: i
        real(kind=8), dimension(nm,n_sp,4,5), intent(inout) :: soilC
        real(kind=8), dimension(n_sp), intent(in) :: Yr_input, Yb_input, Yl_input
        real(kind=8), dimension(n_sp), intent(out) :: Yl_Coutflux, Yr_Coutflux, Yr_C, Yl_C
        real(kind=8), intent(out) :: O_C, O_Coutflux

        do i =1, n_sp

            Yl_Coutflux(i) = sum(soilC(ii-step_aux,i,1,1:4)) - sum(soilC(ii,i,1,1:4)) + Yl_input(i)  
            Yr_Coutflux(i) = sum(soilC(ii-step_aux,i,2:3,1:4)) - sum(soilC(ii,i,2:3,1:4)) + Yr_input(i) + Yb_input(i)

            Yr_C(i) = sum(soilC(ii,i,2:3,1:4))
            Yl_C(i) = sum(soilC(ii,i,1,1:4))

        end do

        O_Coutflux = soilC(ii-step_aux,1,4,5) - soilC(ii,1,4,5)

        soilC(ii,1,4,5) = soilC(ii,1,4,5) + sum(soilC(ii,:,1:3,5))
        soilC(ii,:,1:3,5) = 0.d0

        O_C = soilC(ii,1,4,5)

    end subroutine  yasso_outflux_year 

    ! Get available nitrogen based on stoichiometry pof carbon decomposition
    subroutine yasso_nav (Yr_C, Yl_C, O_C, Yl_Coutflux, Yr_Coutflux, O_Coutflux, humification_N_l, humification_N_r, & 
                        Yr_Noutflux, Yl_Noutflux, O_Noutflux, Yr_N, Yl_N, O_N, hc, qh, el, qbc, qir, qil, er, &
                        Yr_input, Yl_input, Yb_input, soil_class, excessSW, asw, n_depo, Nav, n_sp)
    
        integer, intent(in):: n_sp, soil_class
        integer :: i
        real(kind=8), dimension(n_sp), intent(in) :: hc, qh, el, qbc, er, Yr_input, Yl_input, Yb_input, qir, qil, &
                                                    Yl_Coutflux, Yr_Coutflux,  Yr_C, Yl_C
        real(kind=8), intent(in) :: excessSW, asw, n_depo, O_C, O_Coutflux
        real(kind=8), dimension(n_sp), intent(out) :: Yl_Noutflux, Yr_Noutflux, humification_N_l, &
                                                        humification_N_r, Yr_N, Yl_N
        real(kind=8), intent(out) :: Nav,  O_N,  O_Noutflux
        real(kind=8) :: leaching_min, leaching_org

        do i =1, n_sp
            
            humification_N_l(i) = Yl_Coutflux(i)/(Yl_C(i) + 0.01d0) * hc(i) * (Yl_N(i) / qh(i))
            Yl_Noutflux(i) = max((Yl_Coutflux(i)/(Yl_C(i) + 0.01d0)) * (1-hc(i)) / (1 - el(i)) * &
                    (Yl_N(i) - el(i) * (Yl_C(i) / qbc(i))),0.d0)

            humification_N_r(i) = Yr_Coutflux(i)/(Yr_C(i) + 0.01d0) * hc(i) * (Yr_N(i) / qh(i))
            Yr_Noutflux(i) = max((Yr_Coutflux(i)/(Yr_C(i) + 0.01d0)) * (1-hc(i)) / (1 - er(i)) * &
                        (Yr_N(i) - er(i) * (Yr_C(i) / qbc(i))),0.d0)

            !Now calculate the end-of-month carbon and nitrogen pools					
            Yr_N(i) = Yr_N(i) + ((Yr_input(i) + Yb_input(i)) / (2.d0 * qir(i))) - Yr_Noutflux(i) - humification_N_r(i)
            Yl_N(i) = Yl_N(i) + ((Yl_input(i)) / (2.d0 * qil(i))) - Yl_Noutflux(i) - humification_N_l(i)
    
        end do

        O_Noutflux = (O_Coutflux/O_C) * O_N

        O_N = O_N - O_Noutflux + sum(humification_N_r(:)) + sum(humification_N_l(:))
                
        call compute_leaching( soil_class, leaching_org, leaching_min, excessSW, asw )
    
        Nav = sum(Yr_Noutflux(:)) * (1.d0-leaching_org)  + sum(+ Yl_Noutflux(:)) * (1.d0-leaching_org) + &
                    O_Noutflux * (1.d0 - leaching_min)  + n_depo		
        
    end subroutine  yasso_nav 
    ! End subroutines for Yasso -----------------------------------------------------------


end module soil_cnp