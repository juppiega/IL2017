function IL2017_matlab_example()
% IL2017_matlab_example    An example illustrating the call to IL2017
%     Plots the neutral mass density at 400 km and the exospheric temperature 
%     by calling the IL2017 model through a MATLAB MEX-function. 
%
%     Please start by compiling the MEX-file:
%     mex FCFLAGS="\$FCFLAGS -O3" -output IL2017_mex IL2017_coeff_mod.f90 IL2017_mod.f90  IL2017_mex.F90
%
%     This creates the function IL2017_mex:
%     [rho, temperature, components] = ...
%                    IL2017_mex(day_of_year, altitude, latitude, longitude,...
%                               UTC_hour, F30, F30_average, ae_integrals);
%
%     The example demonstrates, how you would call the function in a MATLAB
%     script.

% Make sure we've compiled all the files:
mex FCFLAGS="\$FCFLAGS -O3" -output IL2017_mex IL2017_coeff_mod.f90 IL2017_mod.f90  IL2017_mex.F90

% Let's start by defining the values that remain constant:
day_of_year = 171;  % Day of year (range from 1 to 366)
altitude = 400;     % Kilometers (WGS84)
longitude = 0;     % Geographic, degrees (WGS84)
F30 = 80;          % 30 cm radio flux for previous day (sfu)
F30_average = 80;  % 30 cm flux for the prevous 81 days (sfu)
ae_integrals = [20,20,20,20]; % Geomagnetic activity -- the exponentially weighed AE integrals

% Varying parameters
UTC_time = 0:0.2:24; % Time since 00 UTC (hours)
latitude = -90:2:90; % Geographic, degrees (WGS84)

% The MEX-file does not take vector inputs, so multiple values have to be
% computed through for-loops.
density = zeros(length(latitude), length(UTC_time));
exospheric_T = zeros(length(latitude), length(UTC_time));
for j = 1:length(UTC_time)
    for i = 1:length(latitude)
        % Call the MEX-function
        [density(i,j), ~, components] = IL2017_mex(day_of_year, altitude, latitude(i), longitude, UTC_time(j), F30, F30_average, ae_integrals);
        % Exospheric temperature is part of the components vector (size:
        % 8), which includes the secondary outputs; the number densities
        % and the temperature parameters.
        exospheric_T(i,j) = components(6);
    end
end

% Plot the results
figure;
[X,Y] = meshgrid(UTC_time, latitude);
surf(X, Y, density, 'edgecolor', 'none');
view(2)
axis tight
colorbar
title('IL2017 Mass density [kg/m^3] @ 400 km')
xlabel('UTC time [h]')
ylabel('Geographic latitude')

figure;
[X,Y] = meshgrid(UTC_time, latitude);
surf(X, Y, exospheric_T, 'edgecolor', 'none');
view(2)
axis tight
colorbar
title('IL2017 Exospheric temperature [K]')
xlabel('UTC time [h]')
ylabel('Geographic latitude')

end