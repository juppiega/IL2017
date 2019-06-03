function compute_AE_integrals()
% compute_AE_integrals    Illustrates, how to compute the AE-integrals for the IL2017 model.
%     Exospheric temperature, O, N2 and He equations require exponentially
%     weighed AE-integrals, each with a different e-folding time (length of
%     system memory):
%     

% Read in the AE timeseries
load AE_supermag.mat % Contains variables called "timestamps" and "AE"

% e-folding times (in hours) for exospheric temperature, O, N2 and He, respectively:
e_foldings = [6.9622, 17.9177, 6.1054, 5.7896];

% Make sure the timeseries is continuous and interpolate to 1 min:
t = (timestamps - timestamps(1))*24*60; 
t1min = t(1):t(end);
AE_continuous = interp1(t, AE, t1min, 'linear', 0);

AE_integrals = zeros(length(AE_continuous), 4);

% Compute the integrals one-by-one:
for i = 1:length(e_foldings)
    % Number of e-foldings, after which the AE-index has no impact on the thermosphere:
    num_e_foldings = 5; % The exact definition calls for an indefinite integral (infinite history), but a cutoff is used for numerical evaluation.
    t_exp = -round(e_foldings(i) * num_e_foldings * 60) : 1 : 0; % For the exponential weighing function.
    weight_function = exp(t_exp / (e_foldings(i) * 60));
    
    % Compute the weighted AE index using cross-correlation:
    integral = xcorr(AE_continuous, weight_function);
    
    % Remove the other branch given by xcorr:
    num_zero_pad = max(length(AE_continuous) - length(weight_function), 0);
    conserve_ind = num_zero_pad + (1:length(AE_continuous));
    
    % Normalize, such that the integrals represent a weighted mean:
    AE_integrals(:,i) = integral(conserve_ind) / sum(weight_function);
end

% Back to actual timestamps
timestamps = t1min/1440 + timestamps(1);

save('AE_integrals.mat','timestamps','AE_integrals')

end 

