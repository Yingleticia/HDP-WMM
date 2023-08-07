%--------- construct G_0 ---------%
% DP(gamma+m, gamma/(gamma+m)H(.|rho,eta) + 1/(gamma+m)sum_{k=1}^K m_k delta_{theta_k})
%%% posterior samples %%%
% data: all obs
% coef = [m1, ... ,m_K, gamma]

function [G0_omega, G0_alpha, G0_beta, K] = construct_G0(data, d, eta, rho, gamma, epsilon, J_group)

dish = unique(data(:,5));
K = length(dish);

% validation
if K ~=max(dish)
    disp('dish estimates have errors')
end

coef = zeros(1, K+1);
alpha_K = zeros(1, K);
beta_K = zeros(1, K);

for kk = 1:K
    alpha_K(kk) = data(find(data(:,5)==kk,1),6);
    beta_K(kk) = data(find(data(:,5)==kk,1),7);
    % number of tables serving dish k
    for jj = 1:J_group
        temp = data((data(:,5)==kk).*(data(:,3)==jj)==1,8);
        if ~isempty(temp)
            coef(kk) = coef(kk) + length(unique(temp));
        end
    end
end
coef(end) = gamma;

% verification %
alpha_temp = unique(data(:,6)); 
beta_temp = unique(data(:,7)); 
alpha_K_sort = sort(alpha_K);
alpha_temp = sort(alpha_temp); 
alpha_temp = alpha_temp';
beta_K_sort = sort(beta_K);
beta_temp = sort(beta_temp); 
beta_temp = beta_temp';
if prod(alpha_K_sort==alpha_temp)~=1
    disp('alpha estimates do not match')
end
if prod(beta_K_sort==beta_temp)~=1
    disp('beta estimates do not match')
end

m = sum(coef(1:K));
R_temp = log(epsilon) / (log(gamma+m)-log(1+gamma+m));
R = min(4000, ceil(R_temp));
G0_omega = zeros(1 ,R); % weights
G0_alpha = zeros(1, R); % alpha
G0_beta = zeros(1, R); % beta
 
% weights %
z_l = betarnd(ones(1,R-1),(gamma+m)*ones(1,R-1));
for r = 1:(R-1)
    G0_omega(r) = z_l(r) * prod(1-z_l(1:(r-1)));
end
G0_omega(R) = 1 - sum(G0_omega(1:(R-1)));

% parameters: alpha + beta %
prob = coef/sum(coef);
k_idx = randsample(K+1, R, true, prob);

for r = 1:R
    if k_idx(r) > K
        G0_alpha(r) = unifrnd(0 ,rho); % uniform distribution
        G0_beta(r) = 1/gamrnd(d, 1/eta); % Inv-Gamma distribution
    else
        G0_alpha(r) = alpha_K(k_idx(r));
        G0_beta(r) = beta_K(k_idx(r));
    end
end

    

end