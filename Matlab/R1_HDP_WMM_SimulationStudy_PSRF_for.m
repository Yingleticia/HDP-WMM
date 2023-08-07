%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% --  Hierarchical Dirichlet Process Weibull Mixture Model          -- %%%
%%% ----------    for prognosis analysis (cancer patients)    ---------- %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      IISE - Revision 1     %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% What is new:
% 1: Evaluation metrics: modeling (goodness-of-fit) + clustering (accuracy)
% 2: Model comparison: DP-WMM (using all data to fit)
% 3: Censoring cases
% 4: Computation time
% 5: Convergence check
% 6: Sensitivity analysis: J_group, Num_J, r^c (censoring rate)
 
clear;clc;close all;

dbstop if error
% dbclear if error
% dbstop if infnan
% dbclear if infnan

%=========================== Parameter Setting ===========================%

seed = 3;
label = 'group'; % 'censor' ; 'group';     'both'; 'none'
censor_rate = 0.0; 
J_group_all = 2:5; % 2 - 5

chain_num = 5;
dataGenerate = 1; % 1: generate new MCMC samples; 0: using saved data

ParaInitialFlag = 'Censor'; % Censor; NoCensor
% Censor: including censor indicator for alpha beta initialization 
% NoCensor: not including
%=========================================================================%

%% MCMC implementation
if dataGenerate == 1
    
    NRep = length(J_group_all);
    for ss = 2:NRep
%     parfor ss = 1:NRep 
        MCMCRep(censor_rate,seed,J_group_all(ss),label,chain_num,ParaInitialFlag);
    end
    
end

%% Convergence diagnostics for each J
if dataGenerate == 0
    
    %%% Load MCMC samples %%%
    J_group = 2;
    folder = ['results/SimulationStudy/',label,'/seed_',num2str(seed),'/',...
        label,'_',num2str(J_group),'/'];
    filename = [folder,'PSRF_seed_',num2str(seed),'_censor_',num2str(censor_rate*10),'.mat'];
    load(filename)
    
    %%% PSRF plots %%%
    %------- selected for each J ------%
    chain_select = [2,5]; % 1 - 5; 
    %----------------------------------%
    
    point_start = 5;
    color_order = rand(30,3);
    Ymax = 0; % 0: no setting
    lineSZ = 1;
    fontSZ = 20;
    lgd = 0;% legend or not
    
    %-------  PSRF computation --------%
    % save(filename,'L','X_max','J_group','weight_all','mu_all','sigma2_all','B',...
    %    'seed','d','censor_rate','data','idx_true','Num_J','N','K0',...
    %    'data_storage_all','gamma_storage_all','alpha_0_storage_all',...
    %    'rho_storage_all','eta_storage_all','comp_time','iterSize','burnin','gap')
    
    %--- gamma ---%
    gamma_chain = gamma_storage_all;
    X = zeros(B,1,chain_num);
    for ch = 1:chain_num
        temp = gamma_chain{ch};
        X(:,1,ch) = temp;
    end
    X = X(:,:,chain_select);
    R_gamma = zeros(1,B-point_start);
    for aa = 1:(B-point_start)
        [R,~,~,~,~] = psrf(X(1:(aa+point_start),:,:));
        R_gamma(aa) = R;
    end
    
    %--- alpha_0 ---%
    alpha_0_chain = alpha_0_storage_all;
    X = zeros(B,1,chain_num);
    for ch = 1:chain_num
        temp = alpha_0_chain{ch};
        X(:,1,ch) = temp;
    end
    X = X(:,:,chain_select);
    R_alpha_0 = zeros(1,B-point_start);
    for aa = 1:(B-point_start)
        [R,~,~,~,~] = psrf(X(1:(aa+point_start),:,:));
        R_alpha_0(aa) = R;
    end
    
    %--- rho ---%
    rho_chain = rho_storage_all;
    X = zeros(B,1,chain_num);
    for ch = 1:chain_num
        temp = rho_chain{ch};
        X(:,1,ch) = temp;
    end
    X = X(:,:,chain_select);
    R_rho = zeros(1,B-point_start);
    for aa = 1:(B-point_start)
        [R,~,~,~,~] = psrf(X(1:(aa+point_start),:,:));
        R_rho(aa) = R;
    end
    
    %--- eta ---%
    eta_chain = eta_storage_all;
    X = zeros(B,1,chain_num);
    for ch = 1:chain_num
        temp = eta_chain{ch};
        X(:,1,ch) = temp;
    end
    X = X(:,:,chain_select);
    R_eta = zeros(1,B-point_start);
    for aa = 1:(B-point_start)
        [R,~,~,~,~] = psrf(X(1:(aa+point_start),:,:));
        R_eta(aa) = R;
    end
    
    %-------  PSRF figure --------%
    ll = 1;
    legend_name=cell(1,30);
    figure(1);
    %--- gamma ---%
    plot(burnin+gap*((point_start+1):B),R_gamma,'-.','LineWidth',lineSZ,'Color',color_order(ll,:)');hold on;
    legend_name{ll}='$\gamma$';
    ll = ll + 1;
    %--- alpha_0 ---%
    plot(burnin+gap*((point_start+1):B),R_alpha_0,'-.','LineWidth',lineSZ,'Color',color_order(ll,:)');hold on;
    legend_name{ll}='$\alpha_0$';
    ll = ll + 1;
    %--- rho ---%
    plot(burnin+gap*((point_start+1):B),R_rho,'-.','LineWidth',lineSZ,'Color',color_order(ll,:)');hold on;
    legend_name{ll}='$\rho$';
    ll = ll + 1;
    %--- eta ---%
    plot(burnin+gap*((point_start+1):B),R_eta,'-.','LineWidth',lineSZ,'Color',color_order(ll,:)');hold on;
    legend_name{ll}='$\eta$';
    ll = ll + 1;
    
    line([burnin,iterSize],[1,1],'linestyle','-','LineWidth',lineSZ+2,'Color',[0 0 0]');
    set(gca,'FontSize',fontSZ,'fontname','times')
    ylabel('PSRF','FontSize',fontSZ)
    xlabel('Iterarion','FontSize',fontSZ)
    legend_name{ll} = 'reference';
    if lgd == 1
        legend(legend_name(1:ll),'Interpreter','latex','FontSize',fontSZ,'fontname','times') % ,'Location','northwest'
        legend('boxoff')
    end
    if Ymax~=0
        ylim([0 Ymax])
    end
    
    
    
end
    


%% function
function [] = MCMCRep(censor_rate,seed,J_group,label,chain_num,ParaInitialFlag)
       
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

%----------- save samples for PSRF plots ------------%

folder = ['results/SimulationStudy/',label,'/seed_',num2str(seed),'/PSRF_For/',...
        label,'_',num2str(J_group),'/'];
comp_time = nan(1,chain_num);   
data_storage_all = cell(1,chain_num);
gamma_storage_all = cell(1,chain_num);
alpha_0_storage_all = cell(1,chain_num);
rho_storage_all = cell(1,chain_num);
eta_storage_all = cell(1,chain_num);

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



%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   Prior inference & Initialization  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(ParaInitialFlag,'Censor')
    [parmHat,~] = wblfit(data(:,1),0.05,data(:,2));
else
    [parmHat,~] = wblfit(data(:,1));
end

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
tbl = tabulate(idx);
while sum(tbl(:,2)>1)~=kk
    kk = kk - 1;
    idx = kmeans(data(:,1),kk); % dish idx
    tbl = tabulate(idx);
end
alpha_kk = zeros(size(idx));
beta_kk = zeros(size(idx));
data = [data,idx,alpha_kk,beta_kk];
for i=1:kk
    dat_temp = data(idx==i,1);
    if strcmp(ParaInitialFlag,'Censor')
        censor_temp = data(idx==i,2);
        [parmHat,~] = wblfit(dat_temp,0.05,censor_temp);
    else
        [parmHat,~] = wblfit(dat_temp);
    end
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
for ch = 1:chain_num
    data = data_original;
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
    comp_time(1,ch) = computation_time;
    data_storage_all{1,ch} = data_record_HDP;
    gamma_storage_all{1,ch} = gamma_record_HDP;
    alpha_0_storage_all{1,ch} = alpha_0_record_HDP;
    rho_storage_all{1,ch} = rho_record_HDP;
    eta_storage_all{1,ch} = eta_record_HDP;
end

filename = [folder,'PSRF_seed_',num2str(seed),'_censor_',num2str(censor_rate*10),'.mat'];
save(filename,'L','X_max','J_group','weight_all','mu_all','sigma2_all','B',...
    'seed','d','censor_rate','data','idx_true','Num_J','N','K0',...
    'data_storage_all','gamma_storage_all','alpha_0_storage_all',...
    'rho_storage_all','eta_storage_all','comp_time','iterSize','burnin','gap')

end



