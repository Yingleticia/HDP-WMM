%--------- New dish sampler for x_ji ---------%
% only for sampling a new dish when a new table is proposed
% step 1: sample alpha using a discretizing method
% step 2: sample beta from its Inv-Gamma posterior
% x: x_ji
% delta: delta_ji (1: censored; 0: not censored)

function [alpha, beta] = NewDishSampler_x_ji(x, delta, d, eta, rho, numDiscrete)

pp = zeros(1, numDiscrete+1);
yy = linspace(0, rho, numDiscrete+1); % numDiscrete intervals
fun = @(alpha)f_x_knew_int(x, delta, d, eta, alpha, rho);

for i=1:numDiscrete
    pp(i+1) = integral(fun, yy(i), yy(i+1));
end

pp = pp/sum(pp);
CDF = cumsum(pp);
r = unifrnd(0, 1);
idx = find(CDF > r,1);
alpha = mean([yy(idx-1), yy(idx)]); % inverse CDF sampling
% debug
if isnan(alpha)
    disp('new alpha sample error (for x_ji)')
    disp(['pp',num2str(pp)])
    disp(['CDF',num2str(CDF)])
    disp(['r',num2str(r)])
    disp(['idx',num2str(idx)])
    disp(['alpha',num2str(mean([yy(idx-1), yy(idx)]))])
end
% Inv-Gamma posterior
aa = d + 1 - delta; % shape
bb = eta + x^alpha; % scale
beta = 1/gamrnd(aa, 1/bb); % Inv-Gamma(aa,bb)
% debug
if isnan(beta)
    disp('new beta sample error (for x_ji)')
    disp(['aa',num2str(aa)])
    disp(['bb',num2str(bb)])
    disp(['beta',num2str(beta)])
end
end

function f = f_x_knew_int(x, delta, d, eta, alpha, rho)

% temp3 = (alpha .* (x.^(alpha-1))).^(1-delta);
% temp4 = (eta + x.^alpha).^(-d-1+delta);
% f =  temp3 .* temp4;

temp1 = (1-delta) .* log(d);
temp2 = d*log(eta) - log(rho);
temp3 = (1-delta) .* (log(alpha) + (alpha-1) .* log(x));
temp4 = (-d-1+delta) .* log(eta + x.^alpha);
f = exp(temp1 + temp2 + temp3 + temp4);

end