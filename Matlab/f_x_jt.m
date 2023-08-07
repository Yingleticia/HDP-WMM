%--------- f_x_jt ---------%
% f(x_jt|x^-jt, k_jt=k)
% dat_jt: all obs at table t in restaurant j
% dat_jt_minus_k: all obs with dish k except table jt


function f = f_x_jt(dat_jt_minus_k, dat_jt, d, eta, rho)

if size(dat_jt_minus_k,1)==0
    f = 0;
else
    denominator = denominator_f_x_jt(dat_jt_minus_k, d, eta, rho);
    numerator = numerator_f_x_jt(dat_jt_minus_k, dat_jt, d, eta, rho);
    if denominator~=0 && ~isnan(numerator)&& ~isnan(denominator)
%         f = numerator/denominator;
        f = exp(log(numerator)-log(denominator));
    else 
        f = 0;
    end
end
% debug
% disp([numerator, denominator])


end
