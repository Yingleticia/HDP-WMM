%--------- Dish sampler for DP-WMM ---------%
% data: all data
% dat: all data except x_i
% i: customer i
% coef_k = [m1, ... ,m_K, gamma]
% pp_k: likelihood for dishes

function [k_i_new, K] = DP_dishSampler(data, i, d, eta, rho, gamma)
x_i = data(i, 1); % obs
delta_i = data(i, 2); % censor indicator

dat = data;
dat(i,:) = [];

K = max(dat(:,5)); % number of dishes

pp_k = zeros(1,K+1);
coef_k = zeros(1, K+1);
for kk = 1:K
    dat_temp = dat(dat(:,5)==kk,:);
    % number of customers having dish k
    coef_k(kk) = size(dat_temp,1);
    
    % likelihood
    if coef_k(kk)~=0 % in case the cluster has no obs
        alpha_k = dat_temp(1,6);
        beta_k = dat_temp(1,7); 
        if delta_i == 1 % censored
            pp_k(kk) = 1 - weibull_cdf(x_i, alpha_k, beta_k);
        elseif delta_i == 0 % not censored
            pp_k(kk) = weibull_pdf(x_i, alpha_k, beta_k);
        end
    end
end
coef_k(end) = gamma;
coef_k = coef_k/sum(coef_k);
pp_k(end) = f_x_ji_knew(x_i, delta_i, d, eta, rho);
prob_new = coef_k .* pp_k;
prob_new = prob_new/sum(prob_new);

k_i_new = randsample(K+1, 1, true, prob_new);

end