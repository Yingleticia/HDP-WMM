%--------- denominator for f_x_ji ---------%
% data: all obs with dish k except x_ji

function f = denominator_f_x_ji(data, d, eta, rho)

fun = @(alpha)denominator_int(data, alpha, d, eta);
f = integral(fun, 0, rho);

% f = int(fun, 0, rho);
% f = integral(fun, 0, rho);
% f = quadgk(fun, 0, rho);

end


function f = denominator_int(data, alpha, d, eta)

X = data(:,1);
X_o = data(data(:,2)==0,1); % not censored
if isempty(X_o)
    X_o = 1;
    n_o = 0;
else
   n_o = length(X_o); 
end

temp1 = n_o .* log(alpha);
temp2 = (alpha-1) .* sum(log(X_o));
temp3 = (-n_o-d) .* log(eta+sum(X.^alpha));
f = exp(temp1 + temp2 + temp3);

% f = (alpha.^n_o) .* (prod(X_o).^(alpha-1)) .* ((eta+sum(X.^alpha)).^(-n_o-d));


end

