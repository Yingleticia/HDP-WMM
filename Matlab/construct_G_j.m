%--------- construct G_j ---------%
% DP(gamma+m, gamma/(gamma+m)H(.|rho,eta) + 1/(gamma+m)sum_{k=1}^K m_k delta_{theta_k})
%%% posterior samples %%%
% data_j: all obs in restaurant j
% G0 weights: G0_omega
% G0 parameters: G0_alpha + G0_beta

function [G_j_omega, G_j_alpha, G_j_beta] = construct_G_j(data_j, G0_omega, G0_alpha, G0_beta, alpha_0, epsilon)

dish = unique(data_j(:,5));
n_j = size(data_j, 1);
K = length(dish);
coef = zeros(1, K+1);
alpha_K = zeros(1, K);
beta_K = zeros(1, K);
for kk = 1:K
    alpha_K(kk) = data_j(find(data_j(:,5)==dish(kk),1),6);
    beta_K(kk) = data_j(find(data_j(:,5)==dish(kk),1),7);
    coef(kk) = sum(data_j(:,5)==dish(kk));
end
coef(end) = alpha_0;

R_temp = log(epsilon) / (log(alpha_0+n_j)-log(1+alpha_0+n_j));
R = min(4000, ceil(R_temp));
G_j_omega = zeros(1 ,R); % weights
G_j_alpha = zeros(1, R); % alpha
G_j_beta = zeros(1, R); % beta

% weights %
z_l = betarnd(ones(1,R-1),(alpha_0+n_j)*ones(1,R-1));
for r = 1:(R-1)
    G_j_omega(r) = z_l(r) * prod(1-z_l(1:(r-1)));
end
G_j_omega(R) = 1 - sum(G_j_omega(1:(R-1)));

% parameters: alpha + beta %
prob = coef/sum(coef);
k_idx = randsample(K+1, R, true, prob);

for r = 1:R
    if k_idx(r) > K
        u = randsample(length(G0_omega), 1, true, G0_omega);
        G_j_alpha(r) = G0_alpha(u); 
        G_j_beta(r) = G0_beta(u); 
    else
        G_j_alpha(r) = alpha_K(k_idx(r));
        G_j_beta(r) = beta_K(k_idx(r));
    end
end





end