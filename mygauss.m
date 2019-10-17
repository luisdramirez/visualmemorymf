%% gaussian function

function sse = mygauss(params, y_data, x)

mu = params(1); % means of the distributions
A = params(3); % amplitude - scaling
sigma = params(2); % width distribution

    
    output_hat = A*exp(-(x-mu).^2/(2*(sigma^2)));
    sse_temp = sum((y_data-output_hat).^2);
   % r2_temp = 1-sum((y_data-output_hat).^2)/sum((y_data- mean(y_data)).^2);


%r2 = -r2_temp; % sum the fits for all surround conditions to minimize
sse = sse_temp;