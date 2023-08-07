

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Survival function inference && Bayesian clustering datasets generating  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;close all;

tic 
%======================= Parameter setting =======================%
seed = 3; % 1; 2; 3
label = 'censor'; % 'censor' ; 'group';     'both'; 'none'
censor_rate = 0.0; % 0.0, 0.1, 0.2, 0.3, 0.4, 0.5
J_group = 2;

non_parallel = 1; % only for 'censor':  J=2; r_c=0.0;
%=================================================================%

fig_num = 0;
SameRange = 1; % 1: same grid values for all groups; 2: group-specific range
tail = 0; % add to the max value for tails
CI_method = 1; % 1: t-distribution; 2: quantile
CI_coef = 0.05;
epsilon = 0.001; % for DP approximation
bound = 0; % plot bound in the figures
hzd = 2; % 1: using posterior samples; 2: pdf/cdf
if hzd == 2
    bound = 0; 
end
%% Load results
if strcmp(label,'group') % sensitivity analysis on group number
    folder = ['results/SimulationStudy/',label,'/seed_',num2str(seed),'/',...
        label,'_',num2str(J_group),'/'];
elseif strcmp(label,'censor') % sensitivity analysis on censoring rate    
    folder = ['results/SimulationStudy/',label,'/seed_',num2str(seed),'/',...
        label,'_',num2str(censor_rate*10),'/'];
    if non_parallel == 1
        folder = ['results/SimulationStudy/',label,'/seed_',num2str(seed),'/',...
            label,'_',num2str(censor_rate*10),'/no-parallel/'];
    end
end

%------------- HDP-WMM --------------%
filename = [folder,'HDP_seed_',num2str(seed),'_censor_',num2str(censor_rate*10),'.mat'];
load(filename)
% 'L','X_max','J_group','weight_all','mu_all','sigma2_all','B',...
%    'seed','d','censor_rate','data','idx_true','Num_J','N','K0',...
%    'data_record_HDP','gamma_record_HDP','alpha_0_record_HDP',...
%    'rho_record_HDP','eta_record_HDP','computation_time'
disp(['HDP-WMM (proposed): Sampling computation time is ',num2str(computation_time),' hours'])

%------------- DP-WMM (All data) --------------%
filename = [folder,'DP_all_seed_',num2str(seed),'_censor_',num2str(censor_rate*10),'.mat'];
load(filename)
% 'L','X_max','J_group','weight_all','mu_all','sigma2_all','B',...
%    'seed','d','censor_rate','data','idx_true','Num_J','N','K0',...
%    'data_record_DP_all','gamma_record_DP_all','rho_record_DP_all',...
%    'eta_record_DP_all','computation_time'
disp(['DP-WMM (All data): Sampling computation time is ',num2str(computation_time),' hours'])

%------------- DP-WMM (Each group) --------------%
filename = [folder,'DP_each_seed_',num2str(seed),'_censor_',num2str(censor_rate*10),'.mat'];
load(filename)
% 'L','X_max','J_group','weight_all','mu_all','sigma2_all','B',...
%    'seed','d','censor_rate','data','idx_true','Num_J','N','K0',...
%    'results_all','computation_time'
%  results_all: group X (data, gamma, rho, eta)
disp(['DP-WMM (Each group): Sampling computation time is ',num2str(computation_time),' hours'])


data = data_record_HDP(:,:,1);

%% Model fitting preparation: 
% model: (1) true (2) HDP-WMM (3) DP-WMM (all data) (4) DP-WMM (each group)
% functions: (1) pdf (2) cdf (3) hazard (4) cumulative hazard

% 1: survival time (1-2)
% 2: censor indicator (0: not censored; 1: censored)
% 3: group (restaurant)
% 4: original survival time
% 5: dish index
% 6: component parameter alpha
% 7: component parameter beta
% 8: table index

%----------- prepare grid values -----------%
X_upper = zeros(1, J_group);
xx = zeros(L, J_group);
yy = zeros(L, J_group);

for j = 1:J_group
    if SameRange == 1
        x_max_temp = X_max;
    else
        x_max_temp = max(data(data(:,3)==j,4));
    end
    X_upper(j) = x_max_temp + tail;
    xx_temp = linspace(0,X_upper(j),L); % grid values
    xx_temp = xx_temp';
    xx(:,j) = xx_temp;
    yy(:,j) = 1 + xx_temp/X_max; 
end

%----------- save true and estimated functions -----------%
% (1) true (2) HDP-WMM (3) DP-WMM (all data) (4) DP-WMM (each group)
PDF_all = cell(1,4);
CDF_all = cell(1,4);
Survival_all = cell(1,4);
hazard_all = cell(1,4);
Hazard_all = cell(1,4);
ClusterDraws_all = cell(1,4);

PDF_B_all = cell(1,4);
CDF_B_all = cell(1,4);
Survival_B_all = cell(1,4);
hazard_B_all = cell(1,4);
Hazard_B_all = cell(1,4);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------             True                ---------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PDF_true = zeros(L, J_group);
CDF_true = zeros(L, J_group);
Survival_true = zeros(L, J_group);
hazard_true = zeros(L, J_group);
Hazard_true = zeros(L, J_group);

if non_parallel == 1
    mu_all = [0.1, 1.0, 1.4];
    sigma2_all = [0.04, 0.02, 0.02];
    weight_all = [0.6, 0.4, 0.0;... % group 1
                  0.3, 0.0, 0.7;... % group 2
                  0.2, 0.8, 0.0;... % group 3
                  0.4, 0.2, 0.4;... % group 4
                  1.0, 0.0, 0.0];   % group 5    
end
for j = 1:J_group
    weight_j = weight_all(j,:);
    weight_j = weight_j';
    PDF_true(:,j) = lognpdf(xx(:,j),mu_all,sqrt(sigma2_all)) * weight_j; 
    CDF_true(:,j) = logncdf(xx(:,j),mu_all,sqrt(sigma2_all)) * weight_j; 

    Survival_true(:,j) = 1 - CDF_true(:,j);
    hazard_true(:,j) = PDF_true(:,j) ./ (1 - CDF_true(:,j)); % hazard = pdf/(1-cdf)
    Hazard_true(:,j) = -log(1 - CDF_true(:,j)); % Hazard = -log(1-cdf)
        
end

%----- save estimates for metrics and figures -----%
PDF_all{1} = PDF_true;
CDF_all{1} = CDF_true;
Survival_all{1} = Survival_true;
hazard_all{1} = hazard_true;
Hazard_all{1} = Hazard_true;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------       Estimation (HDP-WMM)      ---------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PDF_est_B = zeros(L, J_group, B); 
CDF_est_B = zeros(L, J_group, B); 
Survival_est_B = zeros(L, J_group, B);
hazard_est_B = zeros(L, J_group, B); 
Hazard_est_B = zeros(L, J_group, B);  

PDF_est = cell(1, J_group); % mean; lower; upper
CDF_est = cell(1, J_group); % mean; lower; upper
Survival_est = cell(1, J_group); % mean; lower; upper
hazard_est = cell(1, J_group); % mean; lower; upper
Hazard_est = cell(1, J_group); % mean; lower; upper

ClusterDraws = nan(N, B); % for Bayesian clustering
K_est = zeros(1, B); % number of components
K_common_est = zeros(1, B); % number of components that are shared
T_est = zeros(J_group, B); % number of tables

for b = 1:B
    
    % b-th posterior samples
    data_b = data_record_HDP(:,:,b);
    gamma_b = max(gamma_record_HDP);
    alpha_0_b = alpha_0_record_HDP(b);
    rho_b = rho_record_HDP(b);
    eta_b = eta_record_HDP(b);
    ClusterDraws(:,b) = data_b(:,5);
    
    %----------- construct G_0 -----------%
    [G0_omega, G0_alpha, G0_beta, K] = construct_G0(data_b, d, eta_b, rho_b, ...
        gamma_b, epsilon, J_group);
    K_est(b) = K;
    dish_common = 1:K;
    dish_common = dish_common';
    for j = 1:J_group
        data_j = data_b(data_b(:,3)==j,:);
        dish_temp = unique(data_j(:,5));
        dish_common = intersect(dish_common, dish_temp);
    end
    K_common_est(b) = length(dish_common);
     
    %----------- construct G_j -----------%
    for j = 1:J_group
        data_j = data_b(data_b(:,3)==j,:);
        T_est(j,b) = max(data_j(:,8));
        [G_j_omega, G_j_alpha, G_j_beta] = construct_G_j(data_j, G0_omega,...
            G0_alpha, G0_beta, alpha_0_b, epsilon);
        %----------- function inference -----------%
        temp1 = weibull_cdf(1, G_j_alpha, G_j_beta);
        temp = 1 - temp1; % temp2 - temp1; 1 - temp1
        idd = temp==0;
        if sum(idd)~=0
          disp(['The number of components that have Pr(Y<1)=1: ',num2str(sum(idd))])
          disp(['The sum of the corresponding weights: ',num2str(sum(G_j_omega(idd)))])
        end
        for i = 1:L
            x_i = yy(i,j);
            % pdf
            pdf_temp = (1/X_max) * weibull_pdf(x_i, G_j_alpha, G_j_beta)./temp;
            pdf_temp(idd) = 0;
            PDF_est_B(i,j,b) = sum(G_j_omega .* pdf_temp);
            % cdf
            cdf_temp = (weibull_cdf(x_i, G_j_alpha,G_j_beta) - temp1)./temp;
            cdf_temp(idd) = 0;
            CDF_est_B(i,j,b) = sum(G_j_omega .* cdf_temp);
            % survival
            Survival_est_B(i,j,b) = 1 - CDF_est_B(i,j,b);
            % hazard
%             hazard_est(i,j,b) = sum(G_j_omega .* (pdf_temp./(1-cdf_temp)));
            hazard_est_B(i,j,b) = PDF_est_B(i,j,b)/(1-CDF_est_B(i,j,b));
            % Hazard
            Hazard_est_B(i,j,b) = -log(1-CDF_est_B(i,j,b));
        end
    end
   
end


%%% compute posterior mean and confidence interval %%%
%%% PDF %%%
PDF_est_mean = mean(PDF_est_B,3);
% PDF_est_median = median(PDF_est_B,3);
if CI_method==1
    PDF_est_std = std(PDF_est_B,0,3);
    PDF_est_upper = PDF_est_mean + norminv(1 - CI_coef/2) * PDF_est_std/sqrt(B-1);
    PDF_est_lower = PDF_est_mean - norminv(1 - CI_coef/2) * PDF_est_std/sqrt(B-1);
else
    PDF_est_upper = quantile(PDF_est_B,CI_coef/2,3);
    PDF_est_lower = quantile(PDF_est_B,1 - CI_coef/2,3);
end
for jj = 1:J_group
    temp = [PDF_est_mean(:,jj),PDF_est_lower(:,jj),PDF_est_upper(:,jj)];
    PDF_est{jj} = temp;
end

%%% CDF %%%
CDF_est_mean = mean(CDF_est_B,3);
% CDF_est_median = median(CDF_est_B,3);
if CI_method==1
    CDF_est_std = std(CDF_est_B,0,3);
    CDF_est_upper = CDF_est_mean + norminv(1 - CI_coef/2) * CDF_est_std/sqrt(B-1);
    CDF_est_lower = CDF_est_mean - norminv(1 - CI_coef/2) * CDF_est_std/sqrt(B-1);
else
    CDF_est_upper = quantile(CDF_est_B,CI_coef/2,3);
    CDF_est_lower = quantile(CDF_est_B,1 - CI_coef/2,3);
end
for jj = 1:J_group
    temp = [CDF_est_mean(:,jj),CDF_est_lower(:,jj),CDF_est_upper(:,jj)];
    CDF_est{jj} = temp;
end

%%% Survival %%%
Survival_est_mean = mean(Survival_est_B,3);
% Survival_est_median = median(Survival_est_B,3);
if CI_method==1
    Survival_est_std = std(Survival_est_B,0,3);
    Survival_est_upper = Survival_est_mean + norminv(1 - CI_coef/2) * Survival_est_std/sqrt(B-1);
    Survival_est_lower = Survival_est_mean - norminv(1 - CI_coef/2) * Survival_est_std/sqrt(B-1);
else
    Survival_est_upper = quantile(Survival_est_B,CI_coef/2,3);
    Survival_est_lower = quantile(Survival_est_B,1 - CI_coef/2,3);
end
for jj = 1:J_group
    temp = [Survival_est_mean(:,jj),Survival_est_lower(:,jj),Survival_est_upper(:,jj)];
    Survival_est{jj} = temp;
end

%%% hazard %%%
hazard_est_mean = mean(hazard_est_B,3);
% hazard_est_median = median(hazard_est_B,3);
if CI_method==1
    hazard_est_std = std(hazard_est_B,0,3);
    hazard_est_upper = hazard_est_mean + norminv(1 - CI_coef/2) * hazard_est_std/sqrt(B-1);
    hazard_est_lower = hazard_est_mean - norminv(1 - CI_coef/2) * hazard_est_std/sqrt(B-1);
else
    hazard_est_upper = quantile(hazard_est_B,CI_coef/2,3);
    hazard_est_lower = quantile(hazard_est_B,1 - CI_coef/2,3);
end
for jj = 1:J_group
    temp = [hazard_est_mean(:,jj),hazard_est_lower(:,jj),hazard_est_upper(:,jj)];
    hazard_est{jj} = temp;
end

%%% Hazard %%%
Hazard_est_mean = mean(Hazard_est_B,3);
% Hazard_est_median = median(Hazard_est_B,3);
if CI_method==1
    Hazard_est_std = std(Hazard_est_B,0,3);
    Hazard_est_upper = Hazard_est_mean + norminv(1 - CI_coef/2) * Hazard_est_std/sqrt(B-1);
    Hazard_est_lower = Hazard_est_mean - norminv(1 - CI_coef/2) * Hazard_est_std/sqrt(B-1);
else
    Hazard_est_upper = quantile(Hazard_est_B,CI_coef/2,3);
    Hazard_est_lower = quantile(Hazard_est_B,1 - CI_coef/2,3);
end
for jj = 1:J_group
    temp = [Hazard_est_mean(:,jj),Hazard_est_lower(:,jj),Hazard_est_upper(:,jj)];
    Hazard_est{jj} = temp;
end

%----- save estimates for metrics and figures -----%
PDF_all{2} = PDF_est;
CDF_all{2} = CDF_est;
Survival_all{2} = Survival_est;
hazard_all{2} = hazard_est;
Hazard_all{2} = Hazard_est;
ClusterDraws_all{2} = ClusterDraws;


PDF_B_all{2} = PDF_est_B;
CDF_B_all{2} = CDF_est_B;
Survival_B_all{2} = Survival_est_B;
hazard_B_all{2} = hazard_est_B;
Hazard_B_all{2} = Hazard_est_B;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------     Estimation (DP-WMM-All)     ---------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PDF_est_B = zeros(L, J_group, B); 
CDF_est_B = zeros(L, J_group, B); 
Survival_est_B = zeros(L, J_group, B);
hazard_est_B = zeros(L, J_group, B); 
Hazard_est_B = zeros(L, J_group, B);  

PDF_est = cell(1, J_group); % mean; lower; upper
CDF_est = cell(1, J_group); % mean; lower; upper
Survival_est = cell(1, J_group); % mean; lower; upper
hazard_est = cell(1, J_group); % mean; lower; upper
Hazard_est = cell(1, J_group); % mean; lower; upper

ClusterDraws = nan(N, B); % for Bayesian clustering

for b = 1:B
    
    %----------- b-th posterior sample -----------%
    data_b = data_record_DP_all(:,:,b);
    gamma_b = max(gamma_record_DP_all);
    rho_b = rho_record_DP_all(b);
    eta_b = eta_record_DP_all(b);
    ClusterDraws(:,b) = data_b(:,5);
    
    %----------- construct G -----------%
    [G_omega, G_alpha, G_beta] = DP_construct_G(data_b, d, eta_b,...
        rho_b, gamma_b, epsilon);
    
    %----------- function inference -----------%
    temp1 = weibull_cdf(1, G_alpha, G_beta);
    temp = 1 - temp1;
    idd = temp==0;
    if sum(idd)~=0
      disp(['The number of components that have Pr(Y<1)=1: ',num2str(sum(idd))])
      disp(['The sum of the corresponding weights: ',num2str(sum(G_omega(idd)))])
    end
    for j = 1:J_group
        for i = 1:L
            x_i = yy(i,j);
            % pdf
            pdf_temp = (1/X_max) * weibull_pdf(x_i, G_alpha, G_beta)./temp;
            pdf_temp(idd) = 0;
            PDF_est_B(i,j,b) = sum(G_omega .* pdf_temp);
            % cdf
            cdf_temp = (weibull_cdf(x_i, G_alpha,G_beta) - temp1)./temp;
            cdf_temp(idd) = 0;
            CDF_est_B(i,j,b) = sum(G_omega .* cdf_temp);
            % survival
            Survival_est_B(i,j,b) = 1 - CDF_est_B(i,j,b);
            % hazard
            hazard_est_B(i,j,b) = PDF_est_B(i,j,b)/Survival_est_B(i,j,b);
            % Hazard
            Hazard_est_B(i,j,b) = -log(Survival_est_B(i,j,b));
        end
    end
     
end

%%% compute posterior mean and confidence interval %%%
%%% PDF %%%
PDF_est_mean = mean(PDF_est_B,3);
if CI_method==1
    PDF_est_std = std(PDF_est_B,0,3);
    PDF_est_upper = PDF_est_mean + norminv(1 - CI_coef/2) * PDF_est_std/sqrt(B-1);
    PDF_est_lower = PDF_est_mean - norminv(1 - CI_coef/2) * PDF_est_std/sqrt(B-1);
else
    PDF_est_upper = quantile(PDF_est_B,CI_coef/2,3);
    PDF_est_lower = quantile(PDF_est_B,1 - CI_coef/2,3);
end
for jj = 1:J_group
    temp = [PDF_est_mean(:,jj),PDF_est_lower(:,jj),PDF_est_upper(:,jj)];
    PDF_est{jj} = temp;
end

%%% CDF %%%
CDF_est_mean = mean(CDF_est_B,3);
% CDF_est_median = median(CDF_est_B,3);
if CI_method==1
    CDF_est_std = std(CDF_est_B,0,3);
    CDF_est_upper = CDF_est_mean + norminv(1 - CI_coef/2) * CDF_est_std/sqrt(B-1);
    CDF_est_lower = CDF_est_mean - norminv(1 - CI_coef/2) * CDF_est_std/sqrt(B-1);
else
    CDF_est_upper = quantile(CDF_est_B,CI_coef/2,3);
    CDF_est_lower = quantile(CDF_est_B,1 - CI_coef/2,3);
end
for jj = 1:J_group
    temp = [CDF_est_mean(:,jj),CDF_est_lower(:,jj),CDF_est_upper(:,jj)];
    CDF_est{jj} = temp;
end

%%% Survival %%%
Survival_est_mean = mean(Survival_est_B,3);
% Survival_est_median = median(Survival_est_B,3);
if CI_method==1
    Survival_est_std = std(Survival_est_B,0,3);
    Survival_est_upper = Survival_est_mean + norminv(1 - CI_coef/2) * Survival_est_std/sqrt(B-1);
    Survival_est_lower = Survival_est_mean - norminv(1 - CI_coef/2) * Survival_est_std/sqrt(B-1);
else
    Survival_est_upper = quantile(Survival_est_B,CI_coef/2,3);
    Survival_est_lower = quantile(Survival_est_B,1 - CI_coef/2,3);
end
for jj = 1:J_group
    temp = [Survival_est_mean(:,jj),Survival_est_lower(:,jj),Survival_est_upper(:,jj)];
    Survival_est{jj} = temp;
end

%%% hazard %%%
hazard_est_mean = mean(hazard_est_B,3);
% hazard_est_median = median(hazard_est_B,3);
if CI_method==1
    hazard_est_std = std(hazard_est_B,0,3);
    hazard_est_upper = hazard_est_mean + norminv(1 - CI_coef/2) * hazard_est_std/sqrt(B-1);
    hazard_est_lower = hazard_est_mean - norminv(1 - CI_coef/2) * hazard_est_std/sqrt(B-1);
else
    hazard_est_upper = quantile(hazard_est_B,CI_coef/2,3);
    hazard_est_lower = quantile(hazard_est_B,1 - CI_coef/2,3);
end
for jj = 1:J_group
    temp = [hazard_est_mean(:,jj),hazard_est_lower(:,jj),hazard_est_upper(:,jj)];
    hazard_est{jj} = temp;
end

%%% Hazard %%%
Hazard_est_mean = mean(Hazard_est_B,3);
% Hazard_est_median = median(Hazard_est_B,3);
if CI_method==1
    Hazard_est_std = std(Hazard_est_B,0,3);
    Hazard_est_upper = Hazard_est_mean + norminv(1 - CI_coef/2) * Hazard_est_std/sqrt(B-1);
    Hazard_est_lower = Hazard_est_mean - norminv(1 - CI_coef/2) * Hazard_est_std/sqrt(B-1);
else
    Hazard_est_upper = quantile(Hazard_est_B,CI_coef/2,3);
    Hazard_est_lower = quantile(Hazard_est_B,1 - CI_coef/2,3);
end
for jj = 1:J_group
    temp = [Hazard_est_mean(:,jj),Hazard_est_lower(:,jj),Hazard_est_upper(:,jj)];
    Hazard_est{jj} = temp;
end

%----- save estimates for metrics and figures -----%
PDF_all{3} = PDF_est;
CDF_all{3} = CDF_est;
Survival_all{3} = Survival_est;
hazard_all{3} = hazard_est;
Hazard_all{3} = Hazard_est;
ClusterDraws_all{3} = ClusterDraws;


PDF_B_all{3} = PDF_est_B;
CDF_B_all{3} = CDF_est_B;
Survival_B_all{3} = Survival_est_B;
hazard_B_all{3} = hazard_est_B;
Hazard_B_all{3} = Hazard_est_B;

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------     Estimation (DP-WMM-Each)    ---------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PDF_est_B = zeros(L, J_group, B); 
CDF_est_B = zeros(L, J_group, B); 
Survival_est_B = zeros(L, J_group, B);
hazard_est_B = zeros(L, J_group, B); 
Hazard_est_B = zeros(L, J_group, B);  

PDF_est = cell(1, J_group); % mean; lower; upper
CDF_est = cell(1, J_group); % mean; lower; upper
Survival_est = cell(1, J_group); % mean; lower; upper
hazard_est = cell(1, J_group); % mean; lower; upper
Hazard_est = cell(1, J_group); % mean; lower; upper

ClusterDraws = nan(N, B); % for Bayesian clustering
dish_K = 0; % current largest dish indicator -> make sure there is no overlapping among groups

for j = 1:J_group
    
    % cell(J_group,4):  group X (data, gamma, rho, eta)
    data_record_DP_j = results_all{j,1}; 
    gamma_record_DP_j = results_all{j,2}; 
    rho_record_DP_j = results_all{j,3}; 
    eta_record_DP_j = results_all{j,4}; 
    idx_j = (sum(Num_J(1:(j-1)))+1):sum(Num_J(1:j));
    
    for b = 1:B
        %----------- b-th posterior sample -----------%
        data_b = data_record_DP_j(:,:,b);
        gamma_b = max(gamma_record_DP_j);
        rho_b = rho_record_DP_j(b);
        eta_b = eta_record_DP_j(b);
        ClusterDraws(idx_j,b) = data_b(:,5) + dish_K;
        
        %----------- construct G -----------%
        [G_omega, G_alpha, G_beta] = DP_construct_G(data_b, d, eta_b,...
            rho_b, gamma_b, epsilon);
        
        %----------- function inference -----------%
        temp1 = weibull_cdf(1, G_alpha, G_beta);
        temp = 1 - temp1;
        idd = temp==0;
        if sum(idd)~=0
          disp(['The number of components that have Pr(Y<1)=1: ',num2str(sum(idd))])
          disp(['The sum of the corresponding weights: ',num2str(sum(G_omega(idd)))])
        end
        for i = 1:L
            x_i = yy(i,j);
            % pdf
            pdf_temp = (1/X_max) * weibull_pdf(x_i, G_alpha, G_beta)./temp;
            pdf_temp(idd) = 0;
            PDF_est_B(i,j,b) = sum(G_omega .* pdf_temp);
            % cdf
            cdf_temp = (weibull_cdf(x_i, G_alpha,G_beta) - temp1)./temp;
            cdf_temp(idd) = 0;
            CDF_est_B(i,j,b) = sum(G_omega .* cdf_temp);
            % survival
            Survival_est_B(i,j,b) = 1 - CDF_est_B(i,j,b);
            % hazard
            hazard_est_B(i,j,b) = PDF_est_B(i,j,b)/Survival_est_B(i,j,b);
            % Hazard
            Hazard_est_B(i,j,b) = -log(Survival_est_B(i,j,b));
        end
        
    end
    
    dish_K_temp = max(ClusterDraws(idx_j,:),[],'all');
    dish_K = dish_K_temp + dish_K;
end

%%% compute posterior mean and confidence interval %%%
%%% PDF %%%
PDF_est_mean = mean(PDF_est_B,3);
if CI_method==1
    PDF_est_std = std(PDF_est_B,0,3);
    PDF_est_upper = PDF_est_mean + norminv(1 - CI_coef/2) * PDF_est_std/sqrt(B-1);
    PDF_est_lower = PDF_est_mean - norminv(1 - CI_coef/2) * PDF_est_std/sqrt(B-1);
else
    PDF_est_upper = quantile(PDF_est_B,CI_coef/2,3);
    PDF_est_lower = quantile(PDF_est_B,1 - CI_coef/2,3);
end
for jj = 1:J_group
    temp = [PDF_est_mean(:,jj),PDF_est_lower(:,jj),PDF_est_upper(:,jj)];
    PDF_est{jj} = temp;
end

%%% CDF %%%
CDF_est_mean = mean(CDF_est_B,3);
% CDF_est_median = median(CDF_est_B,3);
if CI_method==1
    CDF_est_std = std(CDF_est_B,0,3);
    CDF_est_upper = CDF_est_mean + norminv(1 - CI_coef/2) * CDF_est_std/sqrt(B-1);
    CDF_est_lower = CDF_est_mean - norminv(1 - CI_coef/2) * CDF_est_std/sqrt(B-1);
else
    CDF_est_upper = quantile(CDF_est_B,CI_coef/2,3);
    CDF_est_lower = quantile(CDF_est_B,1 - CI_coef/2,3);
end
for jj = 1:J_group
    temp = [CDF_est_mean(:,jj),CDF_est_lower(:,jj),CDF_est_upper(:,jj)];
    CDF_est{jj} = temp;
end

%%% Survival %%%
Survival_est_mean = mean(Survival_est_B,3);
% Survival_est_median = median(Survival_est_B,3);
if CI_method==1
    Survival_est_std = std(Survival_est_B,0,3);
    Survival_est_upper = Survival_est_mean + norminv(1 - CI_coef/2) * Survival_est_std/sqrt(B-1);
    Survival_est_lower = Survival_est_mean - norminv(1 - CI_coef/2) * Survival_est_std/sqrt(B-1);
else
    Survival_est_upper = quantile(Survival_est_B,CI_coef/2,3);
    Survival_est_lower = quantile(Survival_est_B,1 - CI_coef/2,3);
end
for jj = 1:J_group
    temp = [Survival_est_mean(:,jj),Survival_est_lower(:,jj),Survival_est_upper(:,jj)];
    Survival_est{jj} = temp;
end

%%% hazard %%%
hazard_est_mean = mean(hazard_est_B,3);
% hazard_est_median = median(hazard_est_B,3);
if CI_method==1
    hazard_est_std = std(hazard_est_B,0,3);
    hazard_est_upper = hazard_est_mean + norminv(1 - CI_coef/2) * hazard_est_std/sqrt(B-1);
    hazard_est_lower = hazard_est_mean - norminv(1 - CI_coef/2) * hazard_est_std/sqrt(B-1);
else
    hazard_est_upper = quantile(hazard_est_B,CI_coef/2,3);
    hazard_est_lower = quantile(hazard_est_B,1 - CI_coef/2,3);
end
for jj = 1:J_group
    temp = [hazard_est_mean(:,jj),hazard_est_lower(:,jj),hazard_est_upper(:,jj)];
    hazard_est{jj} = temp;
end

%%% Hazard %%%
Hazard_est_mean = mean(Hazard_est_B,3);
% Hazard_est_median = median(Hazard_est_B,3);
if CI_method==1
    Hazard_est_std = std(Hazard_est_B,0,3);
    Hazard_est_upper = Hazard_est_mean + norminv(1 - CI_coef/2) * Hazard_est_std/sqrt(B-1);
    Hazard_est_lower = Hazard_est_mean - norminv(1 - CI_coef/2) * Hazard_est_std/sqrt(B-1);
else
    Hazard_est_upper = quantile(Hazard_est_B,CI_coef/2,3);
    Hazard_est_lower = quantile(Hazard_est_B,1 - CI_coef/2,3);
end
for jj = 1:J_group
    temp = [Hazard_est_mean(:,jj),Hazard_est_lower(:,jj),Hazard_est_upper(:,jj)];
    Hazard_est{jj} = temp;
end

%----- save estimates for metrics and figures -----%
PDF_all{4} = PDF_est;
CDF_all{4} = CDF_est;
Survival_all{4} = Survival_est;
hazard_all{4} = hazard_est;
Hazard_all{4} = Hazard_est;
ClusterDraws_all{4} = ClusterDraws;


PDF_B_all{4} = PDF_est_B;
CDF_B_all{4} = CDF_est_B;
Survival_B_all{4} = Survival_est_B;
hazard_B_all{4} = hazard_est_B;
Hazard_B_all{4} = Hazard_est_B;

%% Bayesian clustering 
% true labels
idx_true = idx_true';
idx_true = array2table(idx_true);
filename = [folder,'idx_true.xlsx'];
writetable(idx_true,filename,'WriteVariableNames',false)

% estimates
% HDP-WMM
ClusterDraws = ClusterDraws_all{2};
ClusterDraws = ClusterDraws';
ClusterDraws = array2table(ClusterDraws);
filename = [folder,'ClusterDraws_HDP_WMM.xlsx'];
writetable(ClusterDraws,filename,'WriteVariableNames',0)

% DP-WMM-All
ClusterDraws = ClusterDraws_all{3};
ClusterDraws = ClusterDraws';
ClusterDraws = array2table(ClusterDraws);
filename = [folder,'ClusterDraws_DP_WMM_All.xlsx'];
writetable(ClusterDraws,filename,'WriteVariableNames',0)

% DP-WMM-Each
ClusterDraws = ClusterDraws_all{4};
ClusterDraws = ClusterDraws';
ClusterDraws = array2table(ClusterDraws);
filename = [folder,'ClusterDraws_DP_WMM_Each.xlsx'];
writetable(ClusterDraws,filename,'WriteVariableNames',0)


%% Figures 

nbins = 10;
if bound == 1
    legend_name={'Empirical','True','HDP-WMM','HDP-WMM(95% CIs)',...
        'DP-WMM-All','DP-WMM-All(95% CIs)','DP-WMM-Each','DP-WMM-Each(95% CIs)'};
else
    legend_name={'Empirical','True','HDP-WMM','DP-WMM-All','DP-WMM-Each'};
end

line_style = {'-','-.',':','--'};

for j = 1:J_group
    
    fig_num = fig_num + 1;
    figure(fig_num);
    color_order = get(gca,'colororder');
    tiledlayout(1,4,'TileSpacing','compact'); % loose; compact; tight; none
    sur = data(data(:,3)==j,1);
    sur_original = data(data(:,3)==j,4);
    censored = data(data(:,3)==j,2);
    
    %%% PDF %%%
    nexttile(1)
    %-------- Empirical --------%
    histogram(sur_original,nbins,'Normalization','probability','FaceAlpha',0.3, 'FaceColor',color_order(1,:));
    xlim([0 X_upper(j)]) % xlim([1 2])
    set(gca,'FontSize',16,'FontName', 'times')
    xlabel('(a) Density function')
%     title('Density function','FontSize',20,'FontName', 'times', 'FontWeight','normal')
    hold on;
    %-------- Different models --------%
    for mm  = 1:size(PDF_all,2)
        temp = PDF_all{mm};
        if mm == 1
            plot(xx(:,j),temp(:,j),'LineWidth', 2,'linestyle',line_style{mm},'Color',color_order(mm,:))
        else
            temp_j = temp{j};
            plot(xx(:,j),temp_j(:,1),'LineWidth', 2,'linestyle',line_style{mm},'Color',color_order(mm,:))
            if bound == 1 % confidence interval
                x = xx(:,j);
                p_low = temp_j(:,2);
                p_up = temp_j(:,3);
                x = x';
                p_low = p_low';
                p_up = p_up';
                h_temp=fill([x,fliplr(x)],[p_low,fliplr(p_up)],color_order(mm,:));
                set(h_temp,'edgealpha',0,'facealpha',0.1)
            end 
        end
    end
    legend(legend_name,'FontSize',14,'FontName', 'times','Location','best')   
    legend('boxoff')
 
       
    %%% Survival %%%
    nexttile(2)
    %-------- Different models --------%
    for mm  = 1:size(Survival_all,2)
        temp = Survival_all{mm};
        if mm == 1
            plot(xx(:,j),temp(:,j),'LineWidth', 2,'linestyle',line_style{mm},'Color',color_order(mm,:)); hold on;
        else
            temp_j = temp{j};
            plot(xx(:,j),temp_j(:,1),'LineWidth', 2,'linestyle',line_style{mm},'Color',color_order(mm,:))
            if bound == 1 % confidence interval
                x = xx(:,j);
                p_low = temp_j(:,2);
                p_up = temp_j(:,3);
                x = x';
                p_low = p_low';
                p_up = p_up';
                h_temp=fill([x,fliplr(x)],[p_low,fliplr(p_up)],color_order(mm,:));
                set(h_temp,'edgealpha',0,'facealpha',0.1)
            end 
        end
    end
    set(gca,'FontSize',16,'FontName', 'times')
    xlabel('(b) Survival function')
     
    
    %%% hazard %%%
    nexttile(3)
    %-------- Different models --------%
    for mm  = 1:size(hazard_all,2)
        temp = hazard_all{mm};
        if mm == 1
            plot(xx(:,j),temp(:,j),'LineWidth', 2,'linestyle',line_style{mm},'Color',color_order(mm,:)); hold on;
        else
            if hzd == 1
                temp_j = temp{j};
                plot(xx(:,j),temp_j(:,1),'LineWidth', 2,'linestyle',line_style{mm},'Color',color_order(mm,:))
            else
                temp = PDF_all{mm}; 
                temp_j = temp{j};
                pp = temp_j(:,1);
                temp = Survival_all{mm}; 
                temp_j = temp{j};
                ss = temp_j(:,1);                
                plot(xx(:,j),pp./ss,'LineWidth', 2,'linestyle',line_style{mm},'Color',color_order(mm,:))
            end
            if bound == 1 % confidence interval
                x = xx(:,j);
                p_low = temp_j(:,2);
                p_up = temp_j(:,3);
                x = x';
                p_low = p_low';
                p_up = p_up';
                h_temp=fill([x,fliplr(x)],[p_low,fliplr(p_up)],color_order(mm,:));
                set(h_temp,'edgealpha',0,'facealpha',0.1)
            end 
        end
    end
    set(gca,'FontSize',16,'FontName', 'times')
    xlabel('(c) Hazard rate function')
    
       
    %%% Hazard %%%
    nexttile(4)
    %-------- Different models --------%
    for mm  = 1:size(Hazard_all,2)
        temp = Hazard_all{mm};
        if mm == 1
            plot(xx(:,j),temp(:,j),'LineWidth', 2,'linestyle',line_style{mm},'Color',color_order(mm,:)); hold on;
        else
            if hzd == 1
                temp_j = temp{j};
                plot(xx(:,j),temp_j(:,1),'LineWidth', 2,'linestyle',line_style{mm},'Color',color_order(mm,:))
            else
                temp = Survival_all{mm}; 
                temp_j = temp{j};
                ss = temp_j(:,1); 
                plot(xx(:,j),-log(ss),'LineWidth', 2,'linestyle',line_style{mm},'Color',color_order(mm,:))
            end
            if bound == 1 % confidence interval
                x = xx(:,j);
                p_low = temp_j(:,2);
                p_up = temp_j(:,3);
                x = x';
                p_low = p_low';
                p_up = p_up';
                h_temp=fill([x,fliplr(x)],[p_low,fliplr(p_up)],color_order(mm,:));
                set(h_temp,'edgealpha',0,'facealpha',0.1)
            end 
        end
    end
    set(gca,'FontSize',16,'FontName', 'times')
    xlabel('(d) Cumulative hazard function')

    
end

%% Metrics: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------ Model fitting ------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions: (1) pdf (2) survival (3) hazard (4) cumulative hazard
% fitness: true vs estimation (Euclidian distance + MSE + MPE)

method = 1; % 1: true-posterior_mean; 2: mean(true-sample_b)
sur_low = 0.0001; % define valid values for MPE
MSE_fit = zeros(4, 3, J_group); % Function X Model X Groups
MPE_fit = zeros(4, 3, J_group); % Function X Model X Groups
MAE_fit = zeros(4, 3, J_group); % Function X Model X Groups
SAE_fit = zeros(4, 3, J_group); % Function X Model X Groups
% EuclidianDist_fit = zeros(4, 3, J_group); % Function X Model X Groups

for j = 1:J_group
    
    % PDF %
    if method == 1 % posterior_mean
        temp = PDF_all{1};
        fun_true = temp(:,j);
        for mm = 2:size(PDF_all,2)
            temp = PDF_all{mm};
            temp_j = temp{j};
            fun_est = temp_j(:,1);
            MSE_fit(1,mm-1,j) = mean((fun_true-fun_est).^2);
            MAE_fit(1,mm-1,j) = mean(abs(fun_true-fun_est));
            SAE_fit(1,mm-1,j) = sum(abs(fun_true-fun_est));
%             idd_valid = fun_true~=0;
            idd_valid = fun_true>sur_low;
            MPE_fit(1,mm-1,j) = mean(abs(fun_true(idd_valid)-fun_est(idd_valid))...
                ./fun_true(idd_valid));
%             EuclidianDist_fit(1,mm-1,j) = sqrt(sum((fun_true-fun_est).^2));
        end
    else % all samples
%         temp = PDF_all{1};
%         fun_true = temp(:,j);
%         fun_true = repmat(fun_true,1,B);
%         for mm = 2:size(PDF_all,2)
%             PDF_est_B = PDF_B_all{mm};
%             fun_est_B = PDF_est_B(:,j,:);
%             fun_est_B = reshape(fun_est_B,[L,B]);
%             
%         end
        
        
    end
        
    % Survival %
    temp = Survival_all{1};
    fun_true = temp(:,j);
    for mm = 2:size(Survival_all,2)
        temp = Survival_all{mm};
        temp_j = temp{j};
        fun_est = temp_j(:,1);
        MSE_fit(2,mm-1,j) = mean((fun_true-fun_est).^2);
        MAE_fit(2,mm-1,j) = mean(abs(fun_true-fun_est));
        SAE_fit(2,mm-1,j) = sum(abs(fun_true-fun_est));
%         idd_valid = fun_true~=0;
        idd_valid = fun_true>sur_low;
        MPE_fit(2,mm-1,j) = mean(abs(fun_true(idd_valid)-fun_est(idd_valid))...
            ./fun_true(idd_valid));
%         EuclidianDist_fit(2,mm-1,j) = sqrt(sum((fun_true-fun_est).^2));
    end
    
    % hazard %
    temp = hazard_all{1};
    fun_true = temp(:,j);
    for mm = 2:size(hazard_all,2)
        if hzd == 1
            temp = hazard_all{mm};
            temp_j = temp{j};
            fun_est = temp_j(:,1);
        else
            temp = PDF_all{mm};
            temp_j = temp{j};
            pp = temp_j(:,1);
            temp = Survival_all{mm};
            temp_j = temp{j};
            ss = temp_j(:,1);
            fun_est = pp./ss;
        end
        MSE_fit(3,mm-1,j) = mean((fun_true-fun_est).^2);
        MAE_fit(3,mm-1,j) = mean(abs(fun_true-fun_est));
        SAE_fit(3,mm-1,j) = sum(abs(fun_true-fun_est));
%         idd_valid = fun_true~=0;
        idd_valid = fun_true>sur_low;
        MPE_fit(3,mm-1,j) = mean(abs(fun_true(idd_valid)-fun_est(idd_valid))...
            ./fun_true(idd_valid));
%         EuclidianDist_fit(3,mm-1,j) = sqrt(sum((fun_true-fun_est).^2));
    end
    
    % Hazard %
    temp = Hazard_all{1};
    fun_true = temp(:,j);
    for mm = 2:size(Hazard_all,2)
        if hzd == 1
            temp = Hazard_all{mm};
            temp_j = temp{j};
            fun_est = temp_j(:,1);
        else
            temp = Survival_all{mm};
            temp_j = temp{j};
            ss = temp_j(:,1);
            fun_est = -log(ss);
        end
        MSE_fit(4,mm-1,j) = mean((fun_true-fun_est).^2);
        MAE_fit(4,mm-1,j) = mean(abs(fun_true-fun_est));
        SAE_fit(4,mm-1,j) = sum(abs(fun_true-fun_est));
%         idd_valid = fun_true~=0;
        idd_valid = fun_true>sur_low;
        MPE_fit(4,mm-1,j) = mean(abs(fun_true(idd_valid)-fun_est(idd_valid))...
            ./fun_true(idd_valid));
%         EuclidianDist_fit(4,mm-1,j) = sqrt(sum((fun_true-fun_est).^2));
    end
    
end

% disp('===================================== MSE (fitting) =====================================')
% tbl_sum = zeros(size(MSE_fit(:,:,1)));
% for j = 1:J_group
%     disp(['-------------------------- Group ',num2str(j),' --------------------------'])
%     tbl = MSE_fit(:,:,j);
%     tbl_sum = tbl_sum + tbl;
%     tbl = array2table(tbl,'VariableNames',{'HDP-WMM','DP-WMM-All','DP-WMM-Each'});
%     Function = {'PDF';'Survival';'hazard';'Hazard'};
%     Function = cell2table(Function);
%     T = [Function,tbl];
%     disp(T)
% end
% 
% disp('-------------------------- All Groups --------------------------')
% tbl_sum = array2table(tbl_sum,'VariableNames',{'HDP-WMM','DP-WMM-All','DP-WMM-Each'});
% Function = {'PDF';'Survival';'hazard';'Hazard'};
% Function = cell2table(Function);
% T = [Function,tbl_sum];
% disp(T)



% disp('===================================== MPE (fitting) =====================================')
% tbl_sum = zeros(size(MPE_fit(:,:,1)));
% for j = 1:J_group
%     disp(['-------------------------- Group ',num2str(j),' --------------------------'])
%     tbl = MPE_fit(:,:,j);
%     tbl_sum = tbl_sum + tbl;
%     tbl = array2table(tbl,'VariableNames',{'HDP-WMM','DP-WMM-All','DP-WMM-Each'});
%     Function = {'PDF';'Survival';'hazard';'Hazard'};
%     Function = cell2table(Function);
%     T = [Function,tbl];
%     disp(T)
% end
% 
% disp('-------------------------- All Groups --------------------------')
% tbl_sum = array2table(tbl_sum,'VariableNames',{'HDP-WMM','DP-WMM-All','DP-WMM-Each'});
% Function = {'PDF';'Survival';'hazard';'Hazard'};
% Function = cell2table(Function);
% T = [Function,tbl_sum];
% disp(T)


disp('===================================== MAE (fitting) =====================================')
tbl_sum = zeros(size(MAE_fit(:,:,1)));
for j = 1:J_group
    disp(['-------------------------- Group ',num2str(j),' --------------------------'])
    tbl = MAE_fit(:,:,j);
    tbl_sum = tbl_sum + tbl;
    tbl = array2table(tbl,'VariableNames',{'HDP-WMM','DP-WMM-All','DP-WMM-Each'});
    Function = {'PDF';'Survival';'hazard';'Hazard'};
    Function = cell2table(Function);
    T = [Function,tbl];
    disp(T)
end

disp('-------------------------- All Groups --------------------------')
tbl_sum = array2table(tbl_sum,'VariableNames',{'HDP-WMM','DP-WMM-All','DP-WMM-Each'});
Function = {'PDF';'Survival';'hazard';'Hazard'};
Function = cell2table(Function);
T = [Function,tbl_sum];
disp(T)



% disp('===================================== SAE (fitting) =====================================')
% tbl_sum = zeros(size(SAE_fit(:,:,1)));
% for j = 1:J_group
%     disp(['-------------------------- Group ',num2str(j),' --------------------------'])
%     tbl = SAE_fit(:,:,j);
%     tbl_sum = tbl_sum + tbl;
%     tbl = array2table(tbl,'VariableNames',{'HDP-WMM','DP-WMM-All','DP-WMM-Each'});
%     Function = {'PDF';'Survival';'hazard';'Hazard'};
%     Function = cell2table(Function);
%     T = [Function,tbl];
%     disp(T)
% end
% 
% disp('-------------------------- All Groups --------------------------')
% tbl_sum = array2table(tbl_sum,'VariableNames',{'HDP-WMM','DP-WMM-All','DP-WMM-Each'});
% Function = {'PDF';'Survival';'hazard';'Hazard'};
% Function = cell2table(Function);
% T = [Function,tbl_sum];
% disp(T)




% disp('===================================== Euclidian Distance (fitting) =====================================')
% tbl_sum = zeros(size(EuclidianDist_fit(:,:,1)));
% for j = 1:J_group
%     disp(['-------------------------- Group ',num2str(j),' --------------------------'])
%     tbl = EuclidianDist_fit(:,:,j);
%     tbl_sum = tbl_sum + tbl;
%     tbl = array2table(tbl,'VariableNames',{'HDP-WMM','DP-WMM-All','DP-WMM-Each'});
%     Function = {'PDF';'Survival';'hazard';'Hazard'};
%     Function = cell2table(Function);
%     T = [Function,tbl];
%     disp(T)
% end
% 
% disp('-------------------------- All Groups --------------------------')
% tbl_sum = array2table(tbl_sum,'VariableNames',{'HDP-WMM','DP-WMM-All','DP-WMM-Each'});
% Function = {'PDF';'Survival';'hazard';'Hazard'};
% Function = cell2table(Function);
% T = [Function,tbl_sum];
% disp(T)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------ Group difference (only for J=2) ------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% group difference: (1) hazard ratio (2) survival difference 
% metric: percentage of (1) HR>1 (2) survival diff>0

% if J_group == 2
%     
%     percentage = zeros(2, 3); % difference X Model
%     
%     % hazard ratio %    
%     for mm = 2:size(hazard_all,2)
%         temp = hazard_all{1};
%         fun_true = temp(:,1)./temp(:,2);
%         if hzd == 1
%             temp = hazard_all{mm};
%             temp_1 = temp{1};
%             temp_2 = temp{2};
%             fun_est = temp_1./temp_2;
%         else
%             temp = PDF_all{mm};
%             temp_1 = temp{1};
%             pp1 = temp_1(:,1);
%             temp_2 = temp{2};
%             pp2 = temp_2(:,1);
%             temp = Survival_all{mm};
%             temp_1 = temp{1};
%             ss1 = temp_1(:,1);
%             temp_2 = temp{2};
%             ss2 = temp_2(:,1);
%             fun_est = (pp1./ss1)./(pp2./ss2);
%         end
%         
%         idx_valid = ~isnan(fun_true);
%         disp(['True-----invalid index: ',num2str(find(~idx_valid))])
%         idx_valid_est = ~isnan(fun_est);
%         disp(['Model ',num2str(mm),'-----invalid index: ',num2str(find(~idx_valid_est))])
%         idd = idx_valid .* idx_valid_est==1;
%         fun_true = fun_true(idd);
%         fun_est = fun_est(idd);
%         percentage(1,mm-1) = 100*mean((fun_true>1).*(fun_est>1));
%         
%     end
%     
%     
%     % survival difference %
%     temp = Survival_all{1};
%     fun_true = temp(:,1)-temp(:,2);
%     for mm = 2:size(Survival_all,2)
%         temp = Survival_all{mm};
%         temp_1 = temp{1};
%         ss1 = temp_1(:,1);
%         temp_2 = temp{2};
%         ss2 = temp_2(:,1);
%         fun_est = ss1 - ss2;
%         percentage(2,mm-1) = 100*mean((fun_true<0).*(fun_est<0));
%     end
%     
%     
%     
%     disp('===================================== Percentage (difference) =====================================')
%     tbl = percentage;
%     tbl = array2table(tbl,'VariableNames',{'HDP-WMM','DP-WMM-All','DP-WMM-Each'});
%     Function = {'hazard ratio';'survival difference'};
%     Function = cell2table(Function);
%     T = [Function,tbl];
%     disp(T)
% end


toc




