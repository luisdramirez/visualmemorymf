

%%% Fit perceived contrast measures with the normalization model for both perception and visual working memory conditions

function [r2] = fitNormalizationModel_contrastMatch(free_params, fixed_params)

% free parameter
c50_est = free_params(1);
n_est = free_params(2); % Slope
Wi_est = free_params(3); %Suppression weight

% fixed parameters
C_Test = fixed_params{1}; % True center contrast
C_baseline = fixed_params{2}'; % Baseline center contrast estimates
y_data = fixed_params{3}'; % Variable center contrast estimates
surroundContrast = fixed_params{4}; % Surround grating contrast

% Suppression ratio - takes into account the baseline before fitting
C_match = y_data; % - C_baseline;
C_true = C_Test; %- C_baseline;

% compute response to matching data - no surround condition to get estimates of alpha and n 
% Rm = (C_Test.^n_est) ./ (1 + a_est .* (C_Test.^n_est));
Rm = (C_Test.^n_est) ./ ((c50_est.^n_est) + (C_Test.^n_est));

r2_m = 1 - (sum((C_baseline - Rm).^2) / sum((C_baseline - mean(C_baseline)).^2));

% compute the inhibitory weigth necessary to equate Rm = Rt
% Rt = (C_Test.^n_est) ./ (1 + (a_est .* (C_Test.^n_est)) + real((C_Surround*Wi_est).^n_est));

% Rt = (C_Test.^n_est) ./ ((c50_est.^n_est) + (C_Test.^n_est) + ((C_Surround*Wi_est).^n_est));

Rt = (C_Test.^n_est) ./ ((c50_est.^n_est) + (C_Test.^n_est) + (Wi_est*(surroundContrast.^n_est)));


r2_t = 1 - (sum((C_match - Rt).^2) / sum((C_match - mean(C_match)).^2));

% r2 = sum((Rt - C_match).^2);

r2=-sum([r2_t r2_m]);

