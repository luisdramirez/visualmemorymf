%% gaussian function

function r2 = mygauss_allContrasts(params, y_data, x)

numContrasts = 5;

% only allow mean to vary
if length(params) == numContrasts
    mu = params(1:numContrasts); % means of the distributions
    A = ones(1,numContrasts); % amplitude - scaling
    sigma = ones(1,numContrasts); % width distribution
    sigma = params(numContrasts+1:end); % spread

% Allow both mean and width to vary
elseif length(params) == numContrasts*2
    mu = params(1:numContrasts); % means of the distributions
    A = ones(1,numContrasts); % amplitude - scaling
    sigma = params(numContrasts+1:end); % spread
end

for n = 1:length(mu)
    
    output_hat = A(n)*exp(-(x-mu(n)).^2/(2*(sigma(n)^2)));
    r2_temp(n) = 1-sum((y_data(n,:)-output_hat).^2)/sum((y_data(n,:)- mean(y_data(n,:))).^2);

end

r2 = -sum(r2_temp); % sum the fits for all surround conditions to minimize
