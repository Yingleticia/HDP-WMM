%--------- f_x_ji ---------%
% f(x_ji|x^-ji, z_ji=k^new)
% delta: delta_ji (1: censored; 0: not censored) 

function f = f_x_ji_knew(x, delta, d, eta, rho)
fun = @(alpha)f_x_knew_int(x, delta, d, eta, rho, alpha);
% f = int(fun, 0, rho);
f = integral(fun, 0, rho);
end

function f = f_x_knew_int(x, delta, d, eta, rho, alpha)

% temp1 = d^(1-delta);
% temp2 = (eta^d)/rho;
% temp3 = (alpha .* (x.^(alpha-1))).^(1-delta);
% temp4 = (eta + x.^alpha).^(-d-1+delta);
% f = temp1 .* temp2 .* temp3 .* temp4;

temp1 = (1-delta) .* log(d);
temp2 = d*log(eta) - log(rho);
temp3 = (1-delta) .* (log(alpha) + (alpha-1) .* log(x));
temp4 = (-d-1+delta) .* log(eta + x.^alpha);
f = exp(temp1 + temp2 + temp3 + temp4);

end