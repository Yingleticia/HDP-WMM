%--------- Weibull cdf ---------%
% alpha (v): b in matalb --> vector
% beta (w): a^b in matlab --> vector
% x: scalar

function y = weibull_cdf(x, alpha, beta)

temp1 = exp(-(x.^alpha)./beta);

y = 1 - temp1;
end