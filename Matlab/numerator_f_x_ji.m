%--------- numerator for f_x_ji ---------%
% data: all obs with dish k except ji
% x: x_ji
% delta: delta_ji (1: censored; 0: not censored)


function f = numerator_f_x_ji(data, x, delta, d, eta, rho)

fun = @(alpha)numerator_int(data, alpha, d, eta, x, delta);
f = integral(fun, 0, rho);

% f = int(fun, 0, rho);
% f = integral(fun, 0, rho);
% f = quadgk(fun, 0, rho);

end


function f = numerator_int(data, alpha, d, eta, x, delta)

X = data(:,1);
X_o = data(data(:,2)==0,1); % not censored
if isempty(X_o)
    X_o = 1;
    n_o = 0;
else
   n_o = length(X_o); 
end

% temp1 = (n_o + d).^(1-delta); 
% temp2 = alpha.^(n_o + 1 - delta);
% temp3 = (prod(X_o) .* x.^(1-delta)).^(alpha-1);
% temp4 = (eta + sum(X.^alpha) + x.^alpha).^(-n_o-d-1+delta);
% f = temp1 .* temp2 .* temp3 .* temp4;

temp1 = (1-delta) .* log(n_o + d); 
temp2 = (n_o + 1 - delta) .* log(alpha);
temp3 = (alpha-1) .* (sum(log(X_o)) + (1-delta).*log(x));
temp4 = (-n_o-d-1+delta) .* log(eta + sum(X.^alpha) + x.^alpha);
f = exp(temp1 + temp2 + temp3 + temp4);


end


