%--------- p(xji|t^-ji,t_ji = t^new, k) ---------% (Not used!)
% coef = [m1, ... ,m_K, gamma]
% pp: likelihood
% data: all data
% dat: all data except x_ji
% x: x_ji
% delta: delta_ji (1: censored; 0: not censored)
% f: probability of choosing a new table
% prob: probability vector of choosing a dish k for this new table

function [f, prob] = p_xji_tminusji_tnew_k(data, i, delta, d, eta, rho, J_group, gamma)
x = data(i,1);
K = max(data(:,5)); % number of dishes
dat = data;
dat(i,:)=[];

coef = zeros(1, K+1);
pp = zeros(1, K+1);

for kk = 1:K
    dat_temp = dat(dat(:,5)==kk,:);
    pp(kk) = f_x_ji(dat_temp, x, delta, d, eta, rho);
    
    % numbr of tables serving dish k
    for j = 1:J_group
        temp = dat((dat(:,5)==kk).*(dat(:,3)==j)==1,8);
        coef(kk) = coef(kk) + length(unique(temp));
    end
        
end

coef(end) = gamma;
coef = coef/sum(coef);
pp(end) = f_x_ji_knew(x, delta, d, eta, rho);

prob = coef .* pp;
f = sum(prob);

end