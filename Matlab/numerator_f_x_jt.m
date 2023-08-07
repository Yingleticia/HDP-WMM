%--------- numerator for f_x_jt ---------%
% dat_jt: all obs at table t in restaurant j
% dat_jt_minus_k: all obs with dish k except table jt


function f = numerator_f_x_jt(dat_jt_minus_k, dat_jt, d, eta, rho)

fun = @(alpha)numerator_int(dat_jt_minus_k, alpha, d, eta, dat_jt);
f = integral(fun, 0, rho); % , 'AbsTol', 1e-600

% f = int(fun, 0, rho);
% f = integral(fun, 0, rho);
% f = quadgk(fun, 0, rho);

end


function f = numerator_int(dat_jt_minus_k, alpha, d, eta, dat_jt)

X = dat_jt_minus_k(:,1);
X_o = dat_jt_minus_k(dat_jt_minus_k(:,2)==0,1); % not censored
if isempty(X_o)
    X_o = 1;
    n_o = 0;
else
   n_o = length(X_o); 
end

X_t = dat_jt(:,1);
X_t_o = dat_jt(dat_jt(:,2)==0,1); % not censored
if isempty(X_t_o)
    X_t_o = 1;
    n_t_o = 0;
else
   n_t_o = length(X_t_o); 
end

% temp1 = gamma(d + n_t_o + n_o)/gamma(d + n_o);
% temp2 = alpha.^(n_o + n_t_o);
% temp3 = (prod(X_o) .* prod(X_t_o)).^(alpha-1);
% temp4 = (eta + sum(X.^alpha) + sum(X_t.^alpha)).^(-n_o-d-n_t_o);
% f = temp1 .* temp2 .* temp3 .* temp4;


temp1 = gammaln(d + n_t_o + n_o) - gammaln(d + n_o); % gammaln(x) = log(gamma(x))
temp2 = (n_o + n_t_o) .* log(alpha);
temp3 = (alpha-1) .* (sum(log(X_o)) + sum(log(X_t_o)));
temp4 = (-n_o-d-n_t_o) .* log(eta + sum(X.^alpha) + sum(X_t.^alpha));
f = exp(temp1 + temp2 + temp3 + temp4);


end


