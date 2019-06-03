module IL2017_mod
    use IL2017_coeff_mod
    implicit none
    private
    public IL2017

    ! TYPE MODELSTRUCT:
    ! An internal data type used to store the independent variables (spherical harmonics, geopotential height and solar fluxes).
    type modelStruct
        ! PNM and mPNM are the associated Legendre functions in geographic and magnetic coordinates, respectively.
        ! All parameters are scalars
        real(kind = 8) ::            P10, P11, P20, P21, P22, &
                                     P30, P31, P32, P33, P40, P41, P42, P43, P44, P50, P51, &
                                     P52, P53, P60, P61, P62, P63, P71,  &
                                     mP10,mP20,mP30,mP40,mP50,mP60, mP11,mP21,mP31, &
                                     mP41, mP51, mP61, mP22, mP32, mP42, mP52, mP62, &
                                     yv, dv, dv_mag, lv, Z, F, FA
    end type

contains

! subroutine IL2017(day_of_year, altitude, latitude, longitude, UTC_hour, F30, F30_average,&
!                  ae_integrals, rho, temperature, components)
! PURPOSE: The main subroutine of the IL2017 model.
! INPUTS:
!       real*8    day_of_year     : Day of the year (1.0 - 366.0). Can include a fraction.
!       real*8    altitude        : Altitude (WGS84, km)
!       real*8    latitude        : Geographic latitude (WGS84, degrees)
!       real*8    longitude       : Geographic longitude (WGS84, degrees)
!       real*8    UTC_hour        : Number of hours since 00 UTC. Please include a fraction for accuracy.
!       real*8    F30             : Top-of-the-atmosphere daily 30 cm radio flux on the previous day (sfu).
!       real*8    F30_average     : Arithmetic mean of the top-of-the-atmosphere 30 cm radio flux over the preceding 81 days (sfu).
!       real*8    ae_integrals(:) : Exponentially weighted averages of the AE index:
!                                   ae_integrals(1): e-folding time (tau) of 6.9622 hours (for Exospheric temperature)
!                                   ae_integrals(2): tau = 17.9177 hours (Atomic Oxygen)
!                                   ae_integrals(3): tau = 6.1054 hours  (Molecular Nitrogen)
!                                   ae_integrals(4): tau = 5.7896 hours  (Helium)
! OUTPUTS:
!       real*8    rho             : Total (neutral) mass density (kg/m^3)
!       real*8    temperature     : Neutral temperature (K)
!       real*8    components(8)   : Number densities and further parameters
!                                   components(1): O  number density (1/cm^3)
!                                   components(2): N2 number density (1/cm^3)
!                                   components(3): He number density (1/cm^3)
!                                   components(4): Ar number density (1/cm^3) (NOT IMPLEMENTED)
!                                   components(5): O2 number density (1/cm^3)
!                                   components(6): Exospheric temperature (K)
!                                   components(7): 130-km temperature (K)
!                                   components(8): 130-km temperature gradient (K/km)
subroutine IL2017(day_of_year, altitude, latitude, longitude, UTC_hour, F30, F30_average,&
                  ae_integrals, rho, temperature, components)
    implicit none
    ! INPUTS
    real(kind = 8), intent(in) :: day_of_year, altitude, latitude, longitude, UTC_hour, F30, F30_average
    real(kind = 8), intent(in) :: ae_integrals(:)
    ! OUTPUTS
    real(kind = 8), intent(out) :: rho, temperature, components(8)
    ! LOCAL VARIABLES
    type(modelStruct) :: S ! Type for the processed input parameters
    real(kind = 8) :: Tex, T0, dT, O_lb, N2_lb, He_lb, Ar_lb, O2_lb ! Boundary conditions

    ! Input parameter post processin (compute spherical harmonics)
    S = compute_variables_for_fit(day_of_year, altitude, latitude, longitude, UTC_hour, F30, F30_average)

    ! Compute the boundary temperatures.
    Tex = evalTex(S, coeff(TexInd), ae_integrals(1))
    T0 = evalT0(S, T0Coeffs)
    dT = evalDT(S, dTCoeffs)

    ! Compute the densities at the lower boundary.
    O_lb = evalMajorSpecies(S, coeff(OInd), ae_integrals(2))
    N2_lb = evalMajorSpecies(S, coeff(N2Ind), ae_integrals(3))
    He_lb = evalMajorSpecies(S, coeff(HeInd), ae_integrals(4))
    Ar_lb = 0.0 ! Not implemented
    O2_lb = exp(coeff(O2Ind(1)))

    ! Hydrostatic integration of the neutral densities. Computes the final outputs.
    call computeRhoAndTemp(T0, dT, Tex, S%Z, O_lb, N2_lb, He_lb, Ar_lb, O2_lb,&
                           rho, temperature, components(1:5))

    ! Output the boundary temperatures.
    components(6) = Tex
    components(7) = T0
    components(8) = dT

end subroutine

! function evalMajorSpecies(S, speciesCoeff, aeIntegral)
! PURPOSE: Use the model equations (G) to compute the number density of a neutral species at the lower boundary.
function evalMajorSpecies(S, speciesCoeff, aeIntegral)
    implicit none
    type(modelStruct), intent(in) :: S
    real(kind = 8), intent(in) :: speciesCoeff(:), aeIntegral
    real(kind = 8) :: evalMajorSpecies

    evalMajorSpecies = exp(speciesCoeff(1) + G_majorTex(speciesCoeff, S, aeIntegral));

end function

! function evalTex(S, TexCoeff, aeIntegral)
! PURPOSE: Use the model equations to compute the exospheric temperature.
function evalTex(S, TexCoeff, aeIntegral)
    implicit none
    type(modelStruct), intent(in) :: S
    real(kind = 8), intent(in) :: TexCoeff(:), aeIntegral
    real(kind = 8) :: evalTex

    evalTex = TexCoeff(1) * (1 + G_majorTex(TexCoeff, S, aeIntegral))

end function

! function evalT0(S, T0Coeff)
! PURPOSE: Compute the 130-km temperature.
function evalT0(S, T0Coeff)
    implicit none
    type(modelStruct), intent(in) :: S
    real(kind = 8), intent(in) :: T0Coeff(:)
    real(kind = 8) :: evalT0

    evalT0 = T0Coeff(1) * (1 + G_quiet(T0Coeff(2:size(T0Coeff)), S))

end function

! function evaldT(S, dTCoeff)
! PURPOSE: Compute the 130-km temperature gradient.
function evaldT(S, dTCoeff)
    implicit none
    type(modelStruct), intent(in) :: S
    real(kind = 8), intent(in) :: dTCoeff(:)
    real(kind = 8) :: evaldT

    evaldT = dTCoeff(1) * (1 + G_quiet(dTCoeff(2:size(dTCoeff)), S))

end function

! function G_majorTex(a, S, aeIntegral)
! PURPOSE: Upper-level function of the model equation. Sum the quiet and the storm effects.
function G_majorTex(a, S, aeIntegral)
    implicit none
    type(modelStruct), intent(in) :: S
    real(kind = 8), intent(in) :: a(:), aeIntegral
    integer :: k
    real(kind = 8) :: G_majorTex

    k = 1; ! Counter, which helps at counting terms.
    G_majorTex = G_quiet(a(k+1:k+111), S) + G_storm(a(k+112:size(a)), S, aeIntegral);

end function

! function G_quiet(a, S)
! PURPOSE: The quiet-time model equation.
function G_quiet(a, S)
    implicit none

    type(modelStruct), intent(in) :: S
    real(kind = 8), intent(in) :: a(:)
    real(kind = 8) :: G_quiet, latitudeTerm, solarTerm, annual, diurnal, semidiurnal, &
                                   terdiurnal, quaterdiurnal, longitudinal
    integer :: k, dPh, dPy
    real(kind = 8) :: pi

    pi = 4.0 * atan(1.0)

    k = 0; ! Counter.

    latitudeTerm = a(k+1)*S%P10 + a(k+2)*S%P20 + a(k+3)*S%P30 + a(k+4)*S%P40 + a(k+5)*S%P50 + a(k+6)*S%P60 + &
                     a(k+7)*S%FA*S%P10 + a(k+8)*S%FA*S%P20 + a(k+9)*S%FA*S%P30 + a(k+10)*S%FA*S%P40 + &
                     a(k+11)*S%F*S%P10 + a(k+12)*S%F*S%P20 + a(k+13)*S%F*S%P30 + a(k+14)*S%F*S%P40;
    k = k + 14;

    solarTerm = a(k+1)*(S%F-S%FA) + a(k+2)*(S%F-S%FA)**2 + a(k+3)*(S%FA-80) + a(k+4)*(S%FA-80)**2 + a(k+5)*(S%F-S%FA)*(S%FA-80);
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
    diurnal = ((a(k+1)*S%P11 + a(k+2)*S%P31 + a(k+3)*S%P51 + a(k+4)*S%P71 + a(k+5)*(S%FA-80) + a(k+6)*(S%FA-80)**2 + &
                 a(k+7)*(S%F - S%FA)) + (a(k+8)*S%P11 + a(k+9)*S%P21 + a(k+10)*S%P31)*(cos(S%yv-pi*a(dPh))))*cos(S%dv) + &
                ((a(k+12)*S%P11 + a(k+13)*S%P31 + a(k+14)*S%P51 + a(k+15)*S%P71 + a(k+16)*(S%FA-80) + a(k+17)*(S%FA-80)**2 + &
                a(k+18)*(S%F - S%FA)) + (a(k+19)*S%P11 + a(k+20)*S%P21 + a(k+21)*S%P31)*(cos(S%yv-pi*a(dPh))))*sin(S%dv);
    k = k + 21;

    semidiurnal = (a(k+1)*S%P22 + a(k+2)*S%P32 + a(k+3)*S%P52 + a(k+4)*(S%FA-80) + a(k+5)*(S%FA-80)**2 + a(k+6)*(S%F-S%FA) + &
                (a(k+7)*S%P32 + a(k+8)*S%P52)*cos(S%yv-pi*a(dPh)))*cos(2*S%dv) + &
                (a(k+9)*S%P22 + a(k+10)*S%P32 + a(k+11)*S%P52 + a(k+12)*(S%FA-80) + a(k+13)*(S%FA-80)**2 + a(k+14)*(S%F-S%FA) + &
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

! function G_storm(a, S)
! PURPOSE: The storm-time model equation.
function G_storm(a, S, aeIntegral)
    implicit none
    type(modelStruct), intent(in) :: S
    real(kind = 8), intent(in) :: a(:), aeIntegral
    real(kind = 8) :: G_storm
    real(kind = 8) :: an, order_0, order_1, order_2
    integer :: k
    k = 1

    an = a(k+8);
    order_0 = (a(k+1)+a(k+2)*S%mP20+a(k+3)*S%mP40+a(k+4)*S%mP60 + (a(k+5)*S%mP10+a(k+6)*S%mP30+a(k+7)*S%mP50)*&
              cos(S%yv-an))*(1+a(k+9)*S%F)*aeIntegral;
    k = k + 9;!10

    order_1 = (a(k+1)*S%mP11+a(k+2)*S%mP31+a(k+3)*S%mP51 + (a(k+4)*S%mP21+a(k+5)*S%mP41+a(k+6)*S%mP61)*&
                cos(S%yv-an))*(1+a(k+7)*S%F)*aeIntegral*cos(S%dv_mag-a(k+8));
    k = k + 8;!18

    order_2 = (a(k+1)*S%mP22+a(k+2)*S%mP42+a(k+3)*S%mP62 + (a(k+4)*S%mP32+a(k+5)*S%mP52)*&
              cos(S%yv-an))*(1+a(k+6)*S%F)*aeIntegral*cos(2*(S%dv_mag-a(k+7)));
    k = k + 7;!25

    G_storm = order_0 + order_1 + order_2;

end function

! subroutine computeRhoAndTemp(T0, dT0, Tex, Z, OlbDens, N2lbDens, HelbDens, ArlbDens, O2lbDens,&
!                      rho, T, components)
! PURPOSE: Perform hydrostatic integration to compute the neutral densities at the correct altitude.
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

! subroutine convertToMagneticCoordinates(latitude_deg, longitude_deg, altitude, magLat, magLon)
! PURPOSE: Compute the magnetic latitude and longitude for a given geographic position.
subroutine convertToMagneticCoordinates(latitude_deg, longitude_deg, magLat, magLon)
    implicit none
    real(kind = 8), intent(in) :: latitude_deg, longitude_deg
    real(kind = 8), intent(out) :: magLat, magLon
    real(kind = 8) :: pi, latitude, longitude, N, x, y, z, mag_x, mag_y, mag_z, r
    real(kind = 8), parameter :: a = 6378137, f = 1/298.257223563, b = a*(1 - f), e2 = 1 - (b/a)**2

    pi = 4.0 * atan(1.0)

    latitude = latitude_deg*pi/180; longitude = longitude_deg*pi/180;

    ! Convert to ECEF (see e.g. https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_geodetic_to_ECEF_coordinates)
    N = a / sqrt(1 - e2*sin(latitude)**2);
    x = N * cos(latitude) * cos(longitude);
    y = N * cos(latitude) * sin(longitude);
    z = N * (1 - e2) * sin(latitude);

    ! Rotate coordinates to IGRF (e.g. Russel, C. T. (1971), Geophysical coordinate transformations, Cosmic Elecrodynamics, 2, 184-196)
    mag_x = 0.33907*x - 0.91964*y - 0.19826*z
    mag_y = 0.93826*x + 0.34594*y
    mag_z = 0.06859*x + 0.18602*y + 0.98015*z

    r = sqrt(mag_x*mag_x + mag_y*mag_y + mag_z*mag_z)
    magLat = pi/2 - acos(mag_z / r)
    magLat = magLat * 180 / pi
    magLon = 180 * atan2(mag_y, mag_x) / pi

end subroutine

! function computeDeclination (doy) result(declination)
! PURPOSE: Compute declination of the Sun.
! Calculation from:
! U.S. Naval Observatory; U.K. Hydrographic Office, H.M. Nautical Almanac Office (2008),
! The Astronomical Almanac for the Year 2010. U.S. Govt. Printing Office. Section C.
! Year 2000 is assumed
function computeDeclination (doy) result(declination)
    implicit none
    real(kind = 8), intent(in) :: doy
    real(kind = 8) :: declination
    real(kind = 8) :: n, L, g, lambda, eps, pi

    pi = 4.0 * atan(1.0D0)

    n = doy - 1.5;
    L = mod(280.46 + 0.9856474*n, 360D0);
    g = mod(357.528 + 0.9856003*n, 360D0);
    lambda = L + 1.915*sin(g*pi/180) + 0.020*sin(2*g*pi/180);

    eps = 23.439 * pi / 180;
    declination = asin(sin(eps) * sin(lambda*pi/180)) * 180 / pi;

end function

! function computeEquationOfTime (doy) result(eqOfTime)
! PURPOSE: Compute the equation of time.
! Calculation from:
! U.S. Naval Observatory; U.K. Hydrographic Office, H.M. Nautical Almanac Office (2008),
! The Astronomical Almanac for the Year 2010. U.S. Govt. Printing Office. Section C.
! Year 2000 is assumed
function computeEquationOfTime (doy) result(eqOfTime)
    implicit none
    real(kind = 8), intent(in) :: doy
    real(kind = 8) :: eqOfTime
    real(kind = 8) :: n, g, pi

    pi = 4.0 * atan(1.0D0)

    n = doy - 1.5;
    g = mod(357.528 + 0.9856003*n, 360D0) * pi / 180;

    eqOfTime = (-7.659*sin(g) + 9.863*sin(2*g + 3.5932))/60; ! hours

end

! function computeMagneticLocalTime(magLon, doy, UTC_hour) result (magneticLocalTime)
! PURPOSE: Compute the magnetic local time used in the storm equation.
function computeMagneticLocalTime(magLon, doy, UTC_hour) result (magneticLocalTime)
    implicit none
    real(kind = 8), intent(in) :: magLon, doy, UTC_hour
    real(kind = 8) :: subsolarLat, subsolarLon, magneticSubSolarLat, magneticSubSolarLon, magneticLocalTime

    subSolarLat = computeDeclination(doy);
    subSolarLon = 15 * (12 - computeEquationOfTime(doy) - UTC_hour);

    call convertToMagneticCoordinates(subSolarLat, subSolarLon, magneticSubSolarLat, magneticSubSolarLon);

    magneticLocalTime = 12 + (magLon - magneticSubSolarLon) / 15;

end function

! function compute_variables_for_fit(day_of_year, altitude, latitude, longitude, UTC_hour,&
!                                  F30, F30_average) result(S)
! PURPOSE: Further process the input parameters.
function compute_variables_for_fit(day_of_year, altitude, latitude, longitude, UTC_hour,&
                                  F30, F30_average) result(S)
    implicit none
     ! INPUTS
    real(kind = 8), intent(in) :: day_of_year, altitude, latitude, longitude, UTC_hour, F30, F30_average
    ! OUTPUT
    type(modelStruct) :: S
    ! LOCAL VARIABLES
    real(kind = 8) :: x, x_mag, P(0:7), pi, magLat, magLon, MLT, LST
    real(kind = 8), parameter :: R = 6356770, z0 = 130D3

    pi = 4.0 * atan(1.0D0)

    x = sin(latitude * pi / 180.0)
    call convertToMagneticCoordinates(latitude, longitude, magLat, magLon)
    x_mag = sin(magLat * pi / 180.0)

    ! Compute the associated Legendre functions in geographic coordinates.
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
    S%P52 = P(5)

    P = legendre(3, x)
    S%P33 = P(3)
    S%P43 = P(4)
    S%P53 = P(5)
    S%P63 = P(6)

    P = legendre(4, x)
    S%P44 = P(4)

    ! Same for the magnetic coordinates.
    P = legendre(0, x_mag)
    S%mP10 = P(1)
    S%mP20 = P(2)
    S%mP30 = P(3)
    S%mP40 = P(4)
    S%mP50 = P(5)
    S%mP60 = P(6)

    P = legendre(1, x_mag)
    S%mP11 = P(1)
    S%mP21 = P(2)
    S%mP31 = P(3)
    S%mP41 = P(4)
    S%mP51 = P(5)
    S%mP61 = P(6)

    P = legendre(2, x_mag)
    S%mP22 = P(2)
    S%mP32 = P(3)
    S%mP42 = P(4)
    S%mP52 = P(5)
    S%mP62 = P(6)

    ! Variations in the "zonal" direction (yv: yearly-, dv: diurnal-, lv: longitudinal variation)
    S%yv = 2 * pi * (day_of_year - 1) / 365
    LST = UTC_hour + longitude/15 ! Local solar time
    S%dv = 2 * pi * LST / 24
    S%lv = longitude * pi / 180

    ! Magnetic diurnal variation.
    MLT = computeMagneticLocalTime(magLon, day_of_year, UTC_hour)
    S%dv_mag = 2 * pi * MLT / 24

    S%F = F30
    S%FA = F30_average

    ! Geopotential height.
    S%Z = (R + z0) * (altitude - z0 * 1D-3) / (R + altitude*1000)

end function

! Factorial function for the Legendre computation.
function fac(k)
    implicit none
    integer :: k
    integer(kind = 8) :: fac
    integer :: j

    fac = 1
    do j = 2, k
        fac = fac * j
    end do

end function

! function legendre(m,X)
! PURPOSE: Compute the fully-normalized associated Legendre functions for order m (and degree 0 to 7).
function legendre(m,X)
    implicit none
    ! INPUTS
    integer, intent(in) :: m
    real(kind = 8), intent(in) :: X
    ! OUTPUT
    real(kind = 8) :: legendre(0:7)
    ! LOCAL VARIABLES
    integer(kind = 4), parameter :: numPoints = 1, maxDegree = 7
    integer :: n
    real(kind = 8) :: inputX(1), legendre_result(1,0:7)

    inputX(1) = X
    ! Evaluate the polynomials at X.
    call pm_polynomial_value(numPoints, maxDegree, int(m,4), inputX, legendre_result)
    ! Normalize to produce the fully-normalized functions.
    do n = m, 7
        legendre_result(1,n) = (-1)**m * sqrt((n+0.5D0) * fac(n-m) / fac(n+m)) * legendre_result(1,n)
    end do
    legendre = legendre_result(1,:)

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
