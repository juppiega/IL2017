module IL2017_mod
    use IL2017_coeff_mod
    implicit none
    private
    public IL2017

    type modelStruct
        real(kind = 8) ::            P10, P11, P20, P21, P22, &
                                     P30, P31, P32, P33, P40, P41, P42, P43, P44, P50, P51, &
                                     P52, P53, P60, P61, P62, P63, P71,  &
                                     mP10,mP20,mP30,mP40,mP50,mP60,mP70,mP11,mP31,mP51, &
                                     mP21, mP41, mP32, mP52, yv, dv, lv, Z, F, FA

        real(kind = 8) :: aeInt(4)
    end type

contains

subroutine IL2017(day_of_year, altitude, latitude, longitude, solar_time, F30, F30_average,&
                  ae_integrals, rho, Temperature, components)
    implicit none
    ! INPUTS
    real(kind = 8), intent(in) :: day_of_year, altitude, latitude, longitude, solar_time, F30, F30_average
    real(kind = 8), intent(in) :: ae_integrals(:)
    ! OUTPUTS
    real(kind = 8), intent(out) :: rho, Temperature, components(8)
    ! LOCAL VARIABLES
    type(modelStruct) :: S
    real(kind = 8) :: Tex, T0, dT, O_lb, N2_lb, He_lb, Ar_lb, O2_lb ! Boundary conditions

    S = compute_variables_for_fit(day_of_year, altitude, latitude, longitude, solar_time, &
                                  F30, F30_average, ae_integrals)

    Tex = evalTex(S, coeff(TexInd))
    T0 = evalT0(S, T0Coeffs)
    dT = evalDT(S, dTCoeffs)

    O_lb = evalMajorSpecies(S, coeff(OInd))
    N2_lb = evalMajorSpecies(S, coeff(N2Ind))
    He_lb = evalMajorSpecies(S, coeff(HeInd))
    Ar_lb = evalMajorSpecies(S, coeff(ArInd))
    O2_lb = exp(coeff(O2Ind(1)))

    call computeRhoAndTemp(T0, dT, Tex, S%Z, O_lb, N2_lb, He_lb, Ar_lb, O2_lb,&
                           rho, Temperature, components(1:5))

    components(6) = Tex
    components(7) = T0
    components(8) = dT

end subroutine

function evalMajorSpecies(S, speciesCoeff)
    implicit none
    type(modelStruct), intent(in) :: S
    real(kind = 8), intent(in) :: speciesCoeff(:)
    real(kind = 8) :: evalMajorSpecies

    evalMajorSpecies = exp(speciesCoeff(1) + G_majorTex(speciesCoeff, S));

end function

function evalTex(S, TexCoeff)
    implicit none
    type(modelStruct), intent(in) :: S
    real(kind = 8), intent(in) :: TexCoeff(:)
    real(kind = 8) :: evalTex

    evalTex = TexCoeff(1) * (1 + G_majorTex(TexCoeff, S))

end function

function evalT0(S, T0Coeff)
    implicit none
    type(modelStruct), intent(in) :: S
    real(kind = 8), intent(in) :: T0Coeff(:)
    real(kind = 8) :: evalT0

    evalT0 = T0Coeff(1) * (1 + G_quiet(T0Coeff(2:size(T0Coeff)), S))

end function

function evaldT(S, dTCoeff)
    implicit none
    type(modelStruct), intent(in) :: S
    real(kind = 8), intent(in) :: dTCoeff(:)
    real(kind = 8) :: evaldT

    evaldT = dTCoeff(1) * (1 + G_quiet(dTCoeff(2:size(dTCoeff)), S))

end function

function G_majorTex(a, S)
    implicit none
    type(modelStruct), intent(in) :: S
    real(kind = 8), intent(in) :: a(:)
    integer :: k
    real(kind = 8) :: G_majorTex

    k = 1; ! Counter, which helps at counting the terms
    G_majorTex = G_quiet(a(k+1:k+111), S) + G_storm(a(k+112:size(a)), S);

end function

function G_quiet(a, S)
    implicit none
    type(modelStruct), intent(in) :: S
    real(kind = 8), intent(in) :: a(:)
    real(kind = 8) :: G_quiet, latitudeTerm, solarTerm, annual, diurnal, semidiurnal, &
                                   terdiurnal, quaterdiurnal, geomagnetic, geom_symmetric, geom_yearly,&
                                   geom_lst, AE_base, longitudinal
    integer :: k, dPh, numInts, dPy
    integer(kind = 4) :: mexPrintf, i
    real(kind = 8) :: pi
    pi = 4.0 * atan(1.0D0)

    k = 0;
    latitudeTerm = a(k+1)*S%P10 + a(k+2)*S%P20 + a(k+3)*S%P30 + a(k+4)*S%P40 + a(k+5)*S%P50 + a(k+6)*S%P60 + &
                     a(k+7)*S%FA*S%P10 + a(k+8)*S%FA*S%P20 + a(k+9)*S%FA*S%P30 + a(k+10)*S%FA*S%P40 + &
                     a(k+11)*S%F*S%P10 + a(k+12)*S%F*S%P20 + a(k+13)*S%F*S%P30 + a(k+14)*S%F*S%P40;
    k = k + 14;

    solarTerm = a(k+1)*S%F + a(k+2)*S%F**2 + a(k+3)*S%FA + a(k+4)*S%FA**2 + a(k+5)*S%F*S%FA;
    k = k + 5;

    annual = (a(k+1) + a(k+2)*S%P20 + a(k+3)*S%P40)*cos(S%yv-pi*a(k+4))*(1+a(k+5)*S%FA+a(k+6)*S%FA**2) + &
           (a(k+7) + a(k+8)*S%P20)*cos(2*(S%yv-pi*a(k+9)))*(1+a(k+10)*S%FA+a(k+11)*S%FA**2) + &
           (a(k+12) + a(k+13)*S%P20)*cos(3*(S%yv-pi*a(k+14)))*(1+a(k+15)*S%FA) + &
           a(k+16)*cos(4*(S%yv-pi*a(k+17)))*(1+a(k+18)*S%FA);
    k = k + 18;

    annual = annual + (a(k+1)*S%P10 + a(k+2)*S%P30 + a(k+3)*S%P50)*cos(S%yv-pi*a(k+4))*(1+a(k+5)*S%FA+a(k+6)*S%FA**2)+&
           (a(k+7)*S%P10 + a(k+8)*S%P30)*cos(2*(S%yv-pi*a(k+9)))*(1+a(k+10)*S%FA+a(k+11)*S%FA**2) + &
           a(k+12)*S%P10*cos(3*(S%yv-pi*a(k+13)))*(1+a(k+14)*S%FA);
    k = k + 14;

    dPh = k + 11;
    diurnal = ((a(k+1)*S%P11 + a(k+2)*S%P31 + a(k+3)*S%P51 + a(k+4)*S%P71 + a(k+5)*S%FA + a(k+6)*S%FA**2 + &
                 a(k+7)*(S%F - S%FA)) + (a(k+8)*S%P11 + a(k+9)*S%P21 + a(k+10)*S%P31)*(cos(S%yv-pi*a(dPh))))*cos(S%dv) + &
                ((a(k+12)*S%P11 + a(k+13)*S%P31 + a(k+14)*S%P51 + a(k+15)*S%P71 + a(k+16)*S%FA + a(k+17)*S%FA**2 + &
                a(k+18)*(S%F - S%FA)) + (a(k+19)*S%P11 + a(k+20)*S%P21 + a(k+21)*S%P31)*(cos(S%yv-pi*a(dPh))))*sin(S%dv);
    k = k + 21;

    semidiurnal = (a(k+1)*S%P22 + a(k+2)*S%P32 + a(k+3)*S%P52 + a(k+4)*S%FA + a(k+5)*S%FA**2 + a(k+6)*(S%F-S%FA) + &
                (a(k+7)*S%P32 + a(k+8)*S%P52)*cos(S%yv-pi*a(dPh)))*cos(2*S%dv) + &
                (a(k+9)*S%P22 + a(k+10)*S%P32 + a(k+11)*S%P52 + a(k+12)*S%FA + a(k+13)*S%FA**2 + a(k+14)*(S%F-S%FA) + &
                (a(k+15)*S%P32 + a(k+16)*S%P52)*cos(S%yv-pi*a(dPh)))*sin(2*S%dv);
    k = k + 16;

    terdiurnal = (a(k+1)*S%P33 + a(k+2)*S%P53 + (a(k+3)*S%P43 + a(k+4)*S%P63)*cos(S%yv-pi*a(dPh)))*cos(3*S%dv) + &
                   (a(k+5)*S%P33 + a(k+6)*S%P53 + (a(k+7)*S%P43 + a(k+8)*S%P63)*cos(S%yv-pi*a(dPh)))*sin(3*S%dv);
    k = k + 8;

    quaterdiurnal = a(k+1)*S%P44*cos(4*S%dv) + a(k+2)*S%P44*sin(4*S%dv);
    k = k + 2;

    dPy = k+7;
    longitudinal = (1.0 + a(k+1)*S%FA)*(a(k+2)*S%P21+a(k+3)*S%P41+a(k+4)*S%P61 + (a(k+5)*S%P11+a(k+6)*S%P31)*cos(S%yv-pi*a(dPy)))*&
                    cos(S%lv)+&
                     (1.0 + a(k+8)*S%FA)*(a(k+9)*S%P21+a(k+10)*S%P41+a(k+11)*S%P61 + (a(k+12)*S%P11+a(k+13)*S%P31)*&
                    cos(S%yv-pi*a(dPy)))*sin(S%lv);
    k = k + 13;

    G_quiet = latitudeTerm + solarTerm + annual + diurnal + semidiurnal + terdiurnal + quaterdiurnal + longitudinal
end function

function G_storm(a, S)
    implicit none
    type(modelStruct), intent(in) :: S
    real(kind = 8), intent(in) :: a(:)
    real(kind = 8) :: G_storm

    G_storm = 0.0

end function

subroutine computeRhoAndTemp(T0, dT0, Tex, Z, OlbDens, N2lbDens, HelbDens, ArlbDens, O2lbDens,&
                      rho, T, components)
    implicit none
    real(kind = 8), intent(in) :: Tex, dT0, T0, Z, olbDens, N2lbDens, HelbDens, &
                                  ArlbDens, O2lbDens
    real(kind = 8), intent(inout) :: T, rho, components(5)
    real(kind = 8) :: sigma, gamma_O, gamma_N2, gamma_He, gamma_Ar, gamma_O2,f_O,f_N2,f_He,f_Ar,f_O2
    real(kind = 8) :: OnumDens, N2numDens, HeNumDens, ArNumDens, O2NumDens
    real(kind = 8), parameter :: u2kg = 1.660538921E-27, k = 1.38064852E-23, g = 9.418

    sigma = dT0 / (Tex - T0);

    T = Tex - (Tex - T0) * exp(-sigma * (Z));

    gamma_O = 16 * u2kg * g / (sigma*1E-3 * k * Tex);
    f_O = (T0 / T)**(1+gamma_O) * exp(-sigma * (Z) * gamma_O);
    OnumDens = OlbDens*f_O; ! [1/cm^3]

    gamma_N2 = 28 * u2kg * g / (sigma*1E-3 * k * Tex);
    f_N2 = (T0 / T)**(1+gamma_N2) * exp(-sigma * (Z) * gamma_N2);
    N2numDens = N2lbDens*f_N2; ! [1/cm^3]

    gamma_He = 4 * u2kg * g / (sigma*1E-3 * k * Tex);
    f_He = (T0 / T)**(1+gamma_He-0.38) * exp(-sigma * (Z) * gamma_He);
    HeNumDens = HelbDens*f_He; ! [1/cm^3]

    gamma_Ar = 40 * u2kg * g / (sigma*1E-3 * k * Tex);
    f_Ar = (T0 / T)**(1+gamma_Ar) * exp(-sigma * (Z) * gamma_Ar);
    ArNumDens = ArlbDens*f_Ar; ! [1/cm^3]

    gamma_O2 = 32 * u2kg * g / (sigma*1E-3 * k * Tex);
    f_O2 = (T0 / T)**(1+gamma_O2) * exp(-sigma * (Z) * gamma_O2);
    O2NumDens = O2lbDens*f_O2; ! [1/cm^3]

    rho = (16*OnumDens + 28*N2numDens + 4*HeNumDens + 40*ArNumDens + 32*O2NumDens) * u2kg * 1E6; ! [kg/m^3]

    components(1) = OnumDens
    components(2) = N2numDens
    components(3) = HenumDens
    components(4) = ArnumDens
    components(5) = O2numDens

end subroutine

function compute_variables_for_fit(day_of_year, altitude, latitude, longitude, solar_time, &
                                  F30, F30_average, ae_integrals) result(S)
    implicit none
     ! INPUTS
    real(kind = 8), intent(in) :: day_of_year, altitude, latitude, longitude, solar_time, F30, F30_average
    real(kind = 8), intent(in) :: ae_integrals(:)
    ! LOCAL VARIABLES
    type(modelStruct) :: S
    real(kind = 8) :: x, P(0:7), pi
    real(kind = 8), parameter :: R = 6356770, z0 = 130D3
    pi = 4.0 * atan(1.0D0)

    x = sin(latitude * pi / 180.0)

    P = legendre(0, x)
    S%P10 = P(1)
    S%P20 = P(2)
    S%P30 = P(3)
    S%P40 = P(4)
    S%P50 = P(5)
    S%P60 = P(6)

    P = legendre(1, x)
    S%P11 = P(1)
    S%P21 = P(2)
    S%P31 = P(3)
    S%P41 = P(4)
    S%P51 = P(5)
    S%P61 = P(6)
    S%P71 = P(7)

    P = legendre(2, x)
    S%P22 = P(2)
    S%P32 = P(3)
    S%P42 = P(4)
    S%P52 = P(5)
    S%P62 = P(6)

    P = legendre(3, x)
    S%P33 = P(3)
    S%P43 = P(4)
    S%P53 = P(5)
    S%P63 = P(6)

    P = legendre(4, x)
    S%P44 = P(4)

    S%yv = 2 * pi * (day_of_year - 1) / 365
    S%dv = 2 * pi * (solar_time) / 24
    S%lv = longitude * pi / 180

    S%F = F30
    S%FA = F30_average

    S%aeInt = ae_integrals

    S%Z = (R + z0) * (altitude - z0 * 1D-3) / (R + altitude*1000)

end function

function legendre(m,X)
    implicit none
    ! INPUTS
    integer, intent(in) :: m
    real(kind = 8), intent(in) :: X
    ! FUNCTION
    real(kind = 8) :: legendre(0:7)
    ! LOCAL VARIABLES
    integer(kind = 4), parameter :: numPoints = 1
    real(kind = 8) :: inputX(1), legendre_result(1,0:7)

    inputX(1) = X
    call pm_polynomial_value(numPoints, 7, m, inputX, legendre_result)
    legendre = legendre_result(1,:)
    !print *, legendre_result

end function

subroutine pm_polynomial_value ( mm, n, m, x, cx )

!*****************************************************************************80
!
!! PM_POLYNOMIAL_VALUE evaluates the Legendre polynomials Pm(n,m,x).
!
!  Differential equation:
!
!    (1-X*X) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0
!
!  First terms:
!
!    M = 0  ( = Legendre polynomials of first kind P(N,X) )
!
!    Pm(0,0,x) =    1
!    Pm(1,0,x) =    1 X
!    Pm(2,0,x) = (  3 X^2 -   1)/2
!    Pm(3,0,x) = (  5 X^3 -   3 X)/2
!    Pm(4,0,x) = ( 35 X^4 -  30 X^2 +   3)/8
!    Pm(5,0,x) = ( 63 X^5 -  70 X^3 +  15 X)/8
!    Pm(6,0,x) = (231 X^6 - 315 X^4 + 105 X^2 -  5)/16
!    Pm(7,0,x) = (429 X^7 - 693 X^5 + 315 X^3 - 35 X)/16
!
!    M = 1
!
!    Pm(0,1,x) =   0
!    Pm(1,1,x) =   1 * SQRT(1-X^2)
!    Pm(2,1,x) =   3 * SQRT(1-X^2) * X
!    Pm(3,1,x) = 1.5 * SQRT(1-X^2) * (5*X^2-1)
!    Pm(4,1,x) = 2.5 * SQRT(1-X^2) * (7*X^3-3*X)
!
!    M = 2
!
!    Pm(0,2,x) =   0
!    Pm(1,2,x) =   0
!    Pm(2,2,x) =   3 * (1-X^2)
!    Pm(3,2,x) =  15 * (1-X^2) * X
!    Pm(4,2,x) = 7.5 * (1-X^2) * (7*X^2-1)
!
!    M = 3
!
!    Pm(0,3,x) =   0
!    Pm(1,3,x) =   0
!    Pm(2,3,x) =   0
!    Pm(3,3,x) =  15 * (1-X^2)^1.5
!    Pm(4,3,x) = 105 * (1-X^2)^1.5 * X
!
!    M = 4
!
!    Pm(0,4,x) =   0
!    Pm(1,4,x) =   0
!    Pm(2,4,x) =   0
!    Pm(3,4,x) =   0
!    Pm(4,4,x) = 105 * (1-X^2)^2
!
!  Recursion:
!
!    if N < M:
!      Pm(N,M,x) = 0
!    if N = M:
!      Pm(N,M,x) = (2*M-1)!! * (1-X*X)^(M/2) where N!! means the product of
!      all the odd integers less than or equal to N.
!    if N = M+1:
!      Pm(N,M,x) = X*(2*M+1)*Pm(M,M,x)
!    if M+1 < N:
!      Pm(N,M,x) = ( X*(2*N-1)*Pm(N-1,M,x) - (N+M-1)*Pm(N-2,M,x) )/(N-M)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MM, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the maximum first index of the Legendre
!    function, which must be at least 0.
!
!    Input, integer ( kind = 4 ) M, the second index of the Legendre function,
!    which must be at least 0, and no greater than N.
!
!    Input, real ( kind = 8 ) X(MM), the point at which the function is to be
!    evaluated.
!
!    Output, real ( kind = 8 ) CX(MM,0:N), the function values.
!
  implicit none

  integer ( kind = 4 ) mm
  integer ( kind = 4 ) n

  real ( kind = 8 ) cx(mm,0:n)
  real ( kind = 8 ) fact
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  real ( kind = 8 ) x(mm)

  cx(1:mm,0:n) = 0.0D+00
!
!  J = M is the first nonzero function.
!
  if ( m <= n ) then
    cx(1:mm,m) = 1.0D+00

    fact = 1.0D+00
    do j = 1, m
      cx(1:mm,m) = - cx(1:mm,m) * fact * sqrt ( 1.0D+00 - x(1:mm)**2 )
      fact = fact + 2.0D+00
    end do

  end if
!
!  J = M + 1 is the second nonzero function.
!
  if ( m + 1 <= n ) then
    cx(1:mm,m+1) = x(1:mm) * real ( 2 * m + 1, kind = 8 ) * cx(1:mm,m)
  end if
!
!  Now we use a three term recurrence.
!
  do j = m + 2, n
    cx(1:mm,j) = ( real ( 2 * j     - 1, kind = 8 ) * x(1:mm) * cx(1:mm,j-1) &
                 + real (   - j - m + 1, kind = 8 ) *           cx(1:mm,j-2) ) &
                 / real (     j - m,     kind = 8 )
  end do

  return
end subroutine

end module IL2017_mod





program test_il
    use IL2017_mod

    implicit none
    real(kind = 8) :: day_of_year, altitude, latitude, longitude, solar_time, F30, F30_average
    real(kind = 8) :: ae_integrals(4)
    real(kind = 8) :: rho, Temperature, components(8)

    day_of_year = 100
    altitude = 400
    latitude = 0
    longitude = 25
    solar_time = 12
    F30 = 100
    F30_average = 100
    ae_integrals(:) = [0,0,0,0]

    call IL2017(day_of_year, altitude, latitude, longitude, solar_time, F30, F30_average,&
                  ae_integrals, rho, Temperature, components)

    print *, rho
    print *, Temperature
    print *, components
end program test_il
