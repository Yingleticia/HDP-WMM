%--------- f_x_jt ---------%
% f(x_jt|x^-jt, k^new)
% dat_jt: all obs at table t in restaurant j

function f = f_x_jt_knew(dat_jt, d, eta, rho)
fun = @(alpha)f_x_jt_knew_int(dat_jt, d, eta, rho, alpha);
% f = int(fun, 0, rho);
f = integral(fun, 0, rho);
end

function f = f_x_jt_knew_int(dat_jt, d, eta, rho, alpha)

X_t = dat_jt(:,1);
X_t_o = dat_jt(dat_jt(:,2)==0,1); % not censored
if isempty(X_t_o)
    X_t_o = 1;
    n_t_o = 0;
else
   n_t_o = length(X_t_o); 
end

% temp1 = (gamma(d + n_t_o) / gamma(d)) * ((eta^d) / rho);
% temp2 = alpha .^ n_t_o;
% temp3 = prod(X_t_o).^(alpha-1);
% temp4 = (eta + sum(X_t.^alpha)).^(-d-n_t_o);
% f = temp1 .* temp2 .* temp3 .* temp4;

temp1 = gammaln(d + n_t_o) - gammaln(d) + d*log(eta) - log(rho);
temp2 = n_t_o .* log(alpha);
temp3 = (alpha-1) .* sum(log(X_t_o));
temp4 = (-d-n_t_o) .* log(eta + sum(X_t.^alpha));
f = exp(temp1 + temp2 + temp3 + temp4);

end