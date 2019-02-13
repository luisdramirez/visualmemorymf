%%% Maximize r2 for gaussian fit
function [gauss, r2] = gaussfunc(data,params)

mean = params(1); sd = params(2); amp = params(3);

gauss = amp*exp(-0.5*((data-mean/sd)^2);

if sum(~isreal(estThresholds)) > 0
    r2 = Inf;
end
end