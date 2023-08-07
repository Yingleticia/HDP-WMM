%--------- Dish sampler ---------%
% for sampling dish for each table at each restaurant
% p(k_jt = k|t, k^-jt)
% data: all data
% j: restaurant j
% t: table t 
% dat_jt: all obs at table t in restaurant j
% dat_jt_minus: data\dat_jt
% dat_jt_minus_k: all obs with dish k except table jt
% coef = [m1, ... ,m_K, gamma]
% pp: likelihood

function k_new = dishSampler(data, j, t, d, eta, rho, J_group, gamma)

K = max(data(:,5)); % number of dishes
dat_jt = data((data(:,3)==j) .* (data(:,8)==t)==1,:);
dat_jt_minus = data;
dat_jt_minus((data(:,3)==j) .* (data(:,8)==t)==1,:)=[];

coef = zeros(1,K+1);
pp = zeros(1,K+1);

for k = 1:K
    dat_jt_minus_k = dat_jt_minus(dat_jt_minus(:,5)==k,:);
    pp(k) = f_x_jt(dat_jt_minus_k, dat_jt, d, eta, rho);
    
    % number of tables serving dish k
    for jj = 1:J_group
        temp = dat_jt_minus((dat_jt_minus(:,5)==k).*(dat_jt_minus(:,3)==jj)==1,8);
        if ~isempty(temp)
            coef(k) = coef(k) + length(unique(temp));
        end
    end
end

coef(end) = gamma;
pp(end) = f_x_jt_knew(dat_jt, d, eta, rho);
% debug
% disp(pp)
prob = pp .* coef;
prob = prob/sum(prob);
% debug
% disp(prob)
k_new = randsample(K+1, 1, true, prob);

end