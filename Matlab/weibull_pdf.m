%--------- Weibull pdf ---------%
% alpha (v): b in matalb --> vector
% beta (w): a^b in matlab --> vector
% x: scalar

function y = weibull_pdf(x, alpha, beta)

temp1 = alpha./beta;
temp2 = x.^(alpha-1);
temp3 = exp(-(x.^alpha)./beta);

y = temp1 .* temp2 .* temp3;
end
