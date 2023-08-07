%--------- f_x_ji ---------%
% f(x_ji|x^-ji, z_ji=k)
% data: all obs with dish k except ji
% x: x_ji
% delta: delta_ji (1: censored; 0: not censored)


function f = f_x_ji(data, x, delta, d, eta, rho)
if size(data,1)==0
    f = 0;
else
    denominator = denominator_f_x_ji(data, d, eta, rho);
    numerator = numerator_f_x_ji(data, x, delta, d, eta, rho);
    if denominator~=0 && ~isnan(numerator)&& ~isnan(denominator)
%         f = numerator/denominator;
        f = exp(log(numerator)-log(denominator));
    else 
        f = 0;
    end
end


end
