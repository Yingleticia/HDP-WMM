%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% --  Hierarchical Dirichlet Process Weibull Mixture Model          -- %%%
%%% ----------    for prognosis analysis (cancer patients)    ---------- %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      IISE - Revision 1     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% What is new:
% 1: Evaluation metrics: modeling (goodness-of-fit) + clustering (ARI)
% 2: Model comparison: DP-WMM-All; DP-WMM-Each
% 3: Censoring cases
% 4: Computation time
% 5: Convergence check
% 6: Sensitivity analysis: J_group, Num_J, r^c (censoring rate)
 
clear;clc;close all;

% dbstop if error
% dbclear if error
% dbstop if infnan
% dbclear if infnan

%=========================== Parameter Setting ===========================%
seed = 2;

label = 'group'; % 'censor' ; 'group';     'both'; 'none'

censor_rate_all = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]; 
J_group_all = 2:5; % 2 - 5
%=========================================================================%

%% Parallel implementation
if strcmp(label,'group') % sensitivity analysis on group number
    
    censor_rate_all = censor_rate_all(1);
    NRep = length(J_group_all);
    for ss = 1:NRep
%     parfor ss = 1:NRep 
        MCMCRep(censor_rate_all,seed,J_group_all(ss),label);
    end
    
elseif strcmp(label,'censor') % sensitivity analysis on censoring rate
    
    J_group_all = J_group_all(1);
    NRep = length(censor_rate_all);
    for ss = 1:NRep
%     parfor ss = 1:NRep 
        MCMCRep(censor_rate_all(ss),seed,J_group_all,label);
    end
    
elseif strcmp(label,'both') % sensitivity analysis on both
    
    [x, y] = meshgrid(censor_rate_all, J_group_all);
    comb = [x(:),y(:)];
    NRep = size(comb, 1);
    for ss = 1:NRep
%     parfor ss = 1:NRep 
        [x, y] = meshgrid(censor_rate_all, J_group_all);
        comb = [x(:),y(:)];
        censor_rate_ss = comb(ss,1);
        J_group_ss = comb(ss,2);
        MCMCRep(censor_rate_ss,seed,J_group_ss,label);
    end
    
    
elseif strcmp(label,'none') % sensitivity analysis on none (one experiment)
    
    censor_rate_all = censor_rate_all(1);
    J_group_all = J_group_all(1);
    MCMCRep(censor_rate_all,seed,J_group_all,label);
    
end


%%
function [] = MCMCRep(censor_rate,seed,J_group,label)
       
rng(seed);
RandStream.getGlobalStream;

%% Parameter setting
iterSize = 10000;
burnin = 5000;
gap = 20;
B = floor((iterSize-burnin)/gap); % posterior sample size

% model parameter
% epsilon = 0.001; % for DP approximation
L = 200; % for grid values in function inference
numDiscrete = 200; % for discretization sampling


%% Data 

%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulation Data %%%
%%%%%%%%%%%%%%%%%%%%%%%

%------- parameter setting for generating data -------%

Num_jj = 200;

mu_all = [0.1, 1.0, 1.4];
sigma2_all = [0.04, 0.02, 0.02];
weight_all = [0.6, 0.4, 0.0;... % group 1
              0.3, 0.0, 0.7;... % group 2
              0.2, 0.8, 0.0;... % group 3
              0.4, 0.2, 0.4;... % group 4
              1.0, 0.0, 0.0];   % group 5

Num_J = Num_jj * ones(1,J_group);
N = sum(Num_J);
K0 = length(mu_all); % number of components

if strcmp(label,'group') % sensitivity analysis on group number
    folder = ['results/SimulationStudy/',label,'/seed_',num2str(seed),'/',...
        label,'_',num2str(J_group),'/'];
elseif strcmp(label,'censor') % sensitivity analysis on censoring rate    
    folder = ['results/SimulationStudy/',label,'/seed_',num2str(seed),'/',...
        label,'_',num2str(censor_rate*10),'/'];
end
%----------------- generating data ------------------%

% 1: survival time (1-2)
% 2: censor indicator (0: not censored; 1: censored)
% 3: group (restaurant)
% 4: original survival time
% 5: dish index
% 6: component parameter alpha
% 7: component parameter beta
% 8: table index

data = nan(N, 4);
data(:,2) = binornd(1,censor_rate,N,1); % censor indicator
mu_i = zeros(N,1);
sigma2_i = zeros(N,1);
idx_true = nan(N, 1);
for j=1:J_group
    data(sum(Num_J(1:j-1))+1:sum(Num_J(1:j)),3) = j; % group
    idx_temp = randsample(K0,Num_J(j),true,weight_all(j,:));
    mu_i(data(:,3)==j) = mu_all(idx_temp);
    sigma2_i(data(:,3)==j) = sigma2_all(idx_temp);
    idx_true(data(:,3)==j) = idx_temp; % true labels
end
sur_orig = lognrnd(mu_i,sqrt(sigma2_i)); % original survival time
censor_idx = data(:,2)==1;
sur_temp = sur_orig(censor_idx);
sur_rand = unifrnd(0,sur_temp);
sur_orig_2 = sur_orig; % with potential right-censored cases
sur_orig_2(censor_idx) = sur_rand; % disp(sum(sur_orig_2~=sur_orig))

data(:,4) = sur_orig_2;
X_max = max(data(:,4));
data(:,1) = 1 + data(:,4)/X_max; % survival time (1-2)


% figure
if censor_rate == 0
    
    X_upper = zeros(1, J_group);
    xx = zeros(L, J_group);
    yy = zeros(L, J_group);

    for j = 1:J_group
        x_max_temp = X_max;
        X_upper(j) = x_max_temp + 0.5;
        xx_temp = linspace(0,X_upper(j),L); % grid values
        xx_temp = xx_temp';
        xx(:,j) = xx_temp;
        yy(:,j) = 1 + xx_temp/X_max; 
    end

    %--------- True ---------%
    PDF_true = zeros(L, J_group);
    CDF_true = zeros(L, J_group);
    Survival_true = zeros(L, J_group);
    hazard_true = zeros(L, J_group);
    Hazard_true = zeros(L, J_group);

    for j = 1:J_group
        weight_j = weight_all(j,:);
        weight_j = weight_j';
        PDF_true(:,j) = lognpdf(xx(:,j),mu_all,sqrt(sigma2_all)) * weight_j; 
        CDF_true(:,j) = logncdf(xx(:,j),mu_all,sqrt(sigma2_all)) * weight_j; 

        Survival_true(:,j) = 1 - CDF_true(:,j);
        hazard_true(:,j) = PDF_true(:,j) ./ (1 - CDF_true(:,j)); % hazard = pdf/(1-cdf)
        Hazard_true(:,j) = -log(1 - CDF_true(:,j)); % Hazard = -log(1-cdf)
    end

    fig_num = 0;
    fig_num = fig_num + 1;
    figure(fig_num);
    tiledlayout(1,4,'TileSpacing','compact');
    color_order = get(gca,'colororder');
    legend_name={'Group 1(empirical)','Group 1(true)','Group 2(empirical)','Group 2(true)'};
    nexttile(1)
    for j=1:J_group
        sur_original = data(data(:,3)==j,4);
        histogram(sur_original,'Normalization','probability','FaceAlpha',0.3, 'FaceColor',color_order(j,:)); hold on
        xlim([0 X_upper(j)])
        title('Density','FontSize',20,'FontName', 'times', 'FontWeight','normal')
        plot(xx(:,j),PDF_true(:,j),'LineWidth', 2,'Color',color_order(j,:))
    end
    set(gca,'FontSize',16,'FontName', 'times')
    title('Density functions','FontSize',20,'FontName', 'times', 'FontWeight','normal')
    legend(legend_name,'FontSize',14,'FontName', 'times') %,'Location','southeast'
    legend('boxoff')

    nexttile(2)
    legend_name={'Group 1','Group 2'};
    for j=1:J_group
        plot(xx(:,j),Survival_true(:,j),'LineWidth', 2,'Color',color_order(j,:)); hold on
        xlim([0 X_upper(j)])
    end
    set(gca,'FontSize',16,'FontName', 'times')
    title('Survival functions','FontSize',20,'FontName', 'times', 'FontWeight','normal')
    legend(legend_name,'FontSize',14,'FontName', 'times') %,'Location','southeast'
    legend('boxoff')

    nexttile(3)
    for j=1:J_group
        plot(xx(:,j), hazard_true(:,j),'LineWidth', 2,'Color',color_order(j,:)); hold on
        xlim([0 X_upper(j)])

    end
    set(gca,'FontSize',16,'FontName', 'times')
    title('Hazard rate functions','FontSize',20,'FontName', 'times', 'FontWeight','normal') 
    legend(legend_name,'FontSize',14,'FontName', 'times','Location','northwest') %,'Location','southeast'
    legend('boxoff') 

    nexttile(4)
    plot(xx(:,1), hazard_true(:,1)./hazard_true(:,2),'LineWidth', 2,'Color',color_order(1,:));
    set(gca,'FontSize',16,'FontName', 'times')
    title('Hazard ratio','FontSize',20,'FontName', 'times', 'FontWeight','normal') 

end
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   Prior inference & Initialization  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[parmHat,~] = wblfit(data(:,1),0.05,data(:,2));
estimated_alpha = parmHat(2); % shape v
estimated_scale = parmHat(1); % scale
estimated_beta = estimated_scale^estimated_alpha; % w

d = 2;
% eta
a_eta = 1;
b_eta = (sqrt(2) - 1) / estimated_beta;
% eta = gamrnd(a_eta, 1/b_eta);

% rho
a_rho = 2;
b_rho = 4 * estimated_alpha / 3;
% rho = gprnd(1/a_rho, b_rho/a_rho, b_rho); % gprnd(k,sigma,theta)  with k=1/a; sigma=b/a, theta=b
% rho = b_rho/((1-unifrnd(0,1))^(1/a_rho));

% gamma
a_gamma = 1;
b_gamma = 0.001;
gamma = gamrnd(a_gamma, 1/b_gamma); % mean = 2/0.1 = 20

% alpha_0
a_alpha_0 = 1;
b_alpha_0 = 0.001;
alpha_0 = gamrnd(a_alpha_0, 1/b_alpha_0); % mean = 3/0.05 = 60

% dish initialization (K-means)
kk = 5; % number of clusters
idx = kmeans(data(:,1),kk); % dish idx
alpha_kk = zeros(size(idx));
beta_kk = zeros(size(idx));
data = [data,idx,alpha_kk,beta_kk];
for i=1:kk
    dat_temp = data(idx==i,1);
    censor_temp = data(idx==i,2);
    [parmHat,~] = wblfit(dat_temp,0.05,censor_temp);
    estimated_alpha = parmHat(2);
    estimated_scale = parmHat(1); % scale
    estimated_beta = estimated_scale^estimated_alpha;
    data(idx==i,6) = estimated_alpha; % shape u
    data(idx==i,7) = estimated_beta; % w
end

% table initialization
TableNum = zeros(1, J_group); % table counts
table = zeros(size(idx));
data = [data,table];
for i=1:J_group
    idx_temp = data(data(:,3)==i,5);
    idx_unique = unique(idx_temp);
    for j=1:length(idx_unique)
        data((data(:,3)==i) .* (data(:,5)==idx_unique(j))==1,8)=j;
    end
    T_j = max(data(data(:,3)==i,8));
    TableNum(i) = T_j;
end

data_original = data; % save the generated data for DP-WMM comparison

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HDP-WMM (proposed) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
%--------------------- iteration start ---------------------%
disp('========================  HDP-WMM ========================')
[data_record_HDP, gamma_record_HDP, alpha_0_record_HDP,...
    rho_record_HDP, eta_record_HDP] = HDP_WMM_MCMC(data,d,a_rho,b_rho,...
    a_eta,b_eta,a_gamma,b_gamma,gamma,a_alpha_0,b_alpha_0,alpha_0,...
    TableNum,J_group,Num_J,N,B,iterSize,burnin,gap,numDiscrete);
toc
computation_time = toc;
computation_time = computation_time/60/60; % hours

%--------------- save data ---------------%
filename = [folder,'HDP_seed_',num2str(seed),'_censor_',num2str(censor_rate*10),'.mat'];
save(filename,'L','X_max','J_group','weight_all','mu_all','sigma2_all','B',...
    'seed','d','censor_rate','data','idx_true','Num_J','N','K0',...
    'data_record_HDP','gamma_record_HDP','alpha_0_record_HDP',...
    'rho_record_HDP','eta_record_HDP','computation_time')


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  DP-WMM (all data)  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1: survival time (1-2)
% 2: censor indicator (0: not censored; 1: censored)
% 3: group (restaurant)
% 4: original survival time
% 5: dish index
% 6: component parameter alpha
% 7: component parameter beta
% 8: table index

data = data_original; % initialize 'data'
tic
%--------------------- iteration start ---------------------%
disp('========================  DP-WMM (All data) ========================')
[data_record_DP_all, gamma_record_DP_all,rho_record_DP_all, ...
    eta_record_DP_all] = DP_WMM_MCMC(data,d,a_rho,b_rho,...
    a_eta,b_eta,a_gamma,b_gamma,gamma,...
    N,B,iterSize,burnin,gap,numDiscrete);
computation_time = toc;
computation_time = computation_time/60/60; % hours

%--------------- save data ---------------%
filename = [folder,'DP_all_seed_',num2str(seed),'_censor_',num2str(censor_rate*10),'.mat'];
save(filename,'L','X_max','J_group','weight_all','mu_all','sigma2_all','B',...
    'seed','d','censor_rate','data','idx_true','Num_J','N','K0',...
    'data_record_DP_all','gamma_record_DP_all','rho_record_DP_all',...
    'eta_record_DP_all','computation_time')



%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DP-WMM (each group) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = data_original; % initialize 'data'
results_all = cell(J_group,4); % group X (data, gamma, rho, eta)
tic
for jj = 1:J_group
    data_jj = data(data(:,3)==jj,:);
    N_jj = Num_J(jj);
    [data_record_DP_jj, gamma_record_DP_jj,rho_record_DP_jj, ...
        eta_record_DP_jj] = DP_WMM_MCMC(data_jj,d,a_rho,b_rho,...
        a_eta,b_eta,a_gamma,b_gamma,gamma,...
        N_jj,B,iterSize,burnin,gap,numDiscrete);
    results_all{jj,1} = data_record_DP_jj;
    results_all{jj,2} = gamma_record_DP_jj;
    results_all{jj,3} = rho_record_DP_jj;
    results_all{jj,4} = eta_record_DP_jj;
end
computation_time = toc;
computation_time = computation_time/60/60; % hours

%--------------- save data ---------------%
filename = [folder,'DP_each_seed_',num2str(seed),'_censor_',num2str(censor_rate*10),'.mat'];
save(filename,'L','X_max','J_group','weight_all','mu_all','sigma2_all','B',...
    'seed','d','censor_rate','data','idx_true','Num_J','N','K0',...
    'results_all','computation_time')


end



