! Test program for the IL2017 model
! The output should be:
! Mass density [kg/m^3]   :    6.1668469943140263E-012
! Neutral temperature [K] :    1321.6949700126679
program IL2017_test
    use IL2017_mod

    implicit none
    integer :: day_of_year
    real(kind = 8) :: altitude, latitude, longitude, UTC_hour, F30, F30_average
    real(kind = 8) :: ae_integrals(4)
    real(kind = 8) :: rho, temperature, components(8)

    day_of_year = 105
    altitude = 400
    latitude = 10
    longitude = -60
    UTC_hour = 16
    F30 = 100
    F30_average = 100
    ae_integrals(:) = [600D0,600D0,600D0,600D0]

    call IL2017(day_of_year, altitude, latitude, longitude, UTC_hour, F30, F30_average,&
                  ae_integrals, rho, temperature, components)

    print *, 'Mass density [kg/m^3]   : ', rho
    print *, 'Neutral temperature [K] : ', temperature
end program IL2017_test
