%--------- New dish sampler for x_jt ---------%
% sampling a new dish for a table
% step 1: sample alpha using a discretizing method
% step 2: sample beta from its Inv-Gamma posterior
% dat_jt: all obs at table t in restaurant j

function [alpha, beta] = NewDishSampler_x_jt(dat_jt, d, eta, rho, numDiscrete)

X_t = dat_jt(:,1);
X_t_o = dat_jt(dat_jt(:,2)==0,1); % not censored
if isempty(X_t_o)
    X_t_o = 1;
    n_t_o = 0;
else
   n_t_o = length(X_t_o); 
end

pp = zeros(1, numDiscrete+1);
yy = linspace(0, rho, numDiscrete+1); % numDiscrete intervals
fun = @(alpha)f_x_jt_int(X_t, X_t_o, n_t_o, d, eta, alpha, rho);
% debug: figure; fun_yy = fun(yy); plot(fun_yy)

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
    disp('new alpha sample error (for x_jt)')
    disp(['pp',num2str(pp)])
    disp(['CDF',num2str(CDF)])
    disp(['r',num2str(r)])
    disp(['idx',num2str(idx)])
    disp(['alpha',num2str(mean([yy(idx-1), yy(idx)]))])
end

% Inv-Gamma posterior
aa = d + n_t_o; % shape
bb = eta + sum(X_t.^alpha); % scale
beta = 1/gamrnd(aa, 1/bb);
% debug
if isnan(beta)
    disp('new beta sample error (for x_ji)')
    disp(['aa',num2str(aa)])
    disp(['bb',num2str(bb)])
    disp(['beta',num2str(beta)])
end

end

function f = f_x_jt_int(X_t, X_t_o, n_t_o, d, eta, alpha, rho)

% temp1 = alpha .* n_t_o;
% temp2 = prod(X_t_o).^(alpha-1);
% temp3 = (eta + sum(X_t.^alpha)).^(-d-n_t_o);
% f =  temp1 .* temp2 .* temp3;

temp0 = gammaln(d + n_t_o) - gammaln(d) + d*log(eta) - log(rho);
temp1 = n_t_o .* log(alpha);
temp2 = (alpha-1) .* sum(log(X_t_o));
temp3 = (-d-n_t_o) .* log(eta + sum(X_t.^alpha));
f = exp(temp0 + temp1 + temp2 + temp3);

end

