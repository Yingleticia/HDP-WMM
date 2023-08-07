%--------- Table sampler ---------%
% p(t_ji = t|t^-ji,k)
% data: all data
% dat: all data except x_ji
% j: restaurant j
% i: customer i
% data_j: all data in restaurant j
% num = [n_j1, ... ,n_jT, alpha_0]
% pp_t: likelihood for tables
% coef_k = [m1, ... ,m_K, gamma]
% pp_k: likelihood for dishes

function [t_i_new, prob_new] = tableSampler(data, i, d, eta, rho, alpha_0, J_group, gamma)
x_i = data(i, 1); % obs
delta_i = data(i, 2); % censor indicator
j = data(i, 3); % restaurant
t_i = data(i, 8); % table

dat = data;
dat(i,:) = [];
data_j = data(data(:,3)==j,:); % all data in restaurant j
T = max(data_j(:,8)); % number of tables in this restaurant
K = max(data(:,5)); % number of dishes

num = zeros(1,T+1); % number of customers in table t except x_ji
pp_t = zeros(1,T+1);

pp_k = zeros(1,K+1);
coef_k = zeros(1, K+1);
for kk = 1:K
    dat_temp = dat(dat(:,5)==kk,:);
    pp_k(kk) = f_x_ji(dat_temp, x_i, delta_i, d, eta, rho);
    
    % numbr of tables serving dish k
    for j = 1:J_group
        temp = dat((dat(:,5)==kk).*(dat(:,3)==j)==1,8);
        if ~isempty(temp)
            coef_k(kk) = coef_k(kk) + length(unique(temp));
        end
    end
        
end
coef_k(end) = gamma;
coef_k = coef_k/sum(coef_k);
pp_k(end) = f_x_ji_knew(x_i, delta_i, d, eta, rho);
prob_new = coef_k .* pp_k;

for t = 1:T
    tb = data_j(:,8);
    
    if t_i == t
        num(t) = sum(tb==t)-1;
    else
        num(t) = sum(tb==t);
    end
    
    if num(t)==0
        pp_t(t) = 0;
    else
        k_jt = data_j(find(tb==t,1),5); % dish index for table t
        pp_t(t) = pp_k(k_jt);
    end
end

pp_t(end) = sum(prob_new);
num(end) = alpha_0;
prob = pp_t .* num;
prob = prob/sum(prob);

prob_new = prob_new/sum(prob_new);
t_i_new = randsample(T+1, 1, true, prob);

end