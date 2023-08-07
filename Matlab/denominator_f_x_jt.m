%--------- denominator for f_x_jt ---------%
% dat_jt_minus_k: all obs with dish k except table jt

function f = denominator_f_x_jt(dat_jt_minus_k, d, eta, rho)

fun = @(alpha)denominator_int(dat_jt_minus_k, alpha, d, eta);
f = integral(fun, 0, rho); %, 'AbsTol', 1e-300

% f = int(fun, 0, rho);
% f = integral(fun, 0, rho);
% f = quadgk(fun, 0, rho);

end


function f = denominator_int(dat_jt_minus_k, alpha, d, eta)

X = dat_jt_minus_k(:,1);
X_o = dat_jt_minus_k(dat_jt_minus_k(:,2)==0,1); % not censored
if isempty(X_o)
    X_o = 1;
    n_o = 0;
else
   n_o = length(X_o); 
end

% temp1 = alpha.^n_o;
% temp2 = prod(X_o).^(alpha-1);
% temp3 = (eta+sum(X.^alpha)).^(-n_o-d);
% f = temp1 .* temp2 .* temp3;

temp1 = n_o .* log(alpha);
temp2 = (alpha-1) .* sum(log(X_o));
temp3 = (-n_o-d) .* log(eta+sum(X.^alpha));
f = exp(temp1 + temp2 + temp3);

end

