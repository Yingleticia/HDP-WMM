%--------- construct G ---------%
% DP(gamma+N, gamma/(gamma+N)H(.|rho,eta) + 1/(gamma+N)sum_{k=1}^K m_k delta_{theta_k})
%%% posterior samples %%%
% data: all obs
% coef = [m1, ... ,m_K, gamma] --> m_k is the number of customers having dish k

function [G_omega, G_alpha, G_beta] = DP_construct_G(data, d, eta, rho, gamma, epsilon)

dish = unique(data(:,5));
K = length(dish);

% validation
if K ~=max(dish)
    disp('dish estimates have errors: K~=max(unique(dish indicators))')
end

coef = zeros(1, K+1);
alpha_K = zeros(1, K);
beta_K = zeros(1, K);

for kk = 1:K
    data_temp = data(data(:,5)==kk,:);
    coef(kk) = size(data_temp,1);
    % validation
    if coef(kk) == 0
       disp('dish estimates have errors: some dish does not have customers')
    end
    alpha_K(kk) = data_temp(1,6);
    beta_K(kk) = data_temp(1,7);
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

N = size(data,1);
if N ~= sum(coef(1:K))
    disp('Sum of all dish counts do not match with N')
end

R_temp = log(epsilon) / (log(gamma+N)-log(1+gamma+N));
R = min(4000, ceil(R_temp));
G_omega = zeros(1 ,R); % weights
G_alpha = zeros(1, R); % alpha
G_beta = zeros(1, R); % beta

% weights %
z_l = betarnd(ones(1,R-1),(gamma+N)*ones(1,R-1));
for r = 1:(R-1)
    G_omega(r) = z_l(r) * prod(1-z_l(1:(r-1)));
end
G_omega(R) = 1 - sum(G_omega(1:(R-1)));

% parameters: alpha + beta %
prob = coef/sum(coef);
k_idx = randsample(K+1, R, true, prob);

for r = 1:R
    if k_idx(r) > K
        G_alpha(r) = unifrnd(0 ,rho); % uniform distribution
        G_beta(r) = 1/gamrnd(d, 1/eta); % Inv-Gamma distribution
    else
        G_alpha(r) = alpha_K(k_idx(r));
        G_beta(r) = beta_K(k_idx(r));
    end
end




end