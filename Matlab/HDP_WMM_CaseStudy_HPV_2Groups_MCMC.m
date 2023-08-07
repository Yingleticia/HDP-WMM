%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% --  Hierarchical Dirichlet Process Weibull Mixture Model          -- %%%
%%% ----------    for prognosis analysis (cancer patients)    ---------- %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear;clc;close all;
tic

seed=3;       
rng(seed);
RandStream.getGlobalStream;

%============ Problem Setting ============%
groupVar = 'HPV'; 
SurvivalFlag = 'Overall'; % Specific; Overall
Operation = 'FunEst'; % MCMC; FunEst; All
site = 'Oropharynx'; % All_Site; Oropharynx;
%=========================================%

folder = ['results/CaseStudy/HPV/2Groups/',site,'/seed_',num2str(seed),'/'];
disp(['---> Seed = ',num2str(seed)])
disp(['---> Selected group variable: ',groupVar])
disp(['---> Selected site: ',site])
disp(['---> Selected Survival Flag: ',SurvivalFlag])

%% Parameter setting
iterSize = 10000;
burnin = 5000;
gap = 50;
B = floor((iterSize-burnin)/gap); % posterior sample size

% model parameter
epsilon = 0.001; % for DP approximation
L = 200; % for grid values in function inference
numDiscrete = 200; % for discretization sampling

% function estimation
year_extent = 5; % extented range (years)
medianEst = 0; 
CI_coef = 0.05;
if medianEst == 0 % posterior mean
    CI_method = 1; % 1: t-distribution
elseif medianEst == 1 % posterior median
    CI_method = 2; % 2: quantile
end


%% Data 

%%%%%%%%%%%%%%%%%
%%% Real Data %%%
%%%%%%%%%%%%%%%%%

%------------ Data Extraction ------------%
if strcmp(Operation,'MCMC') || strcmp(Operation,'All')
    
load('CaseStudyData/HNC_HPV_08222022.mat') % 'data_HNC_HPV' 
N_total = size(data_HNC_HPV, 1); % 22888
disp(['The number of all HNC patients with known HPV status: ',num2str(N_total)])
VarNames = data_HNC_HPV.Properties.VariableNames;

%%%%%====== data selection ======%%%%%
if strcmp(site,'Oropharynx')
    ColId = ismember(VarNames,{'Primary Site - labeled'});
    temp = {'C09.9-Tonsil, NOS','C09.0-Tonsillar fossa',...
                'C10.9-Oropharynx, NOS','C09.1-Tonsillar pillar',...
                'C10.8-Overlapping lesion of oropharynx',...
                'C10.0-Vallecula','C10.3-Posterior wall of oropharynx',...
                'C09.8-Overlapping lesion of tonsil',...
                'C10.2-Lateral wall of oropharynx',...
                ''};
    var = data_HNC_HPV(:,ColId);
    var = table2cell(var);
    data_select = data_HNC_HPV(ismember(var,temp),:);
    disp(['The number of Oropharynx patients: ',num2str(size(data_select,1))])
elseif strcmp(site,'All_Site')    
    data_select = data_HNC_HPV;
end


%%%%%====== variables ======%%%%%

% original survival time (months) 
ColId = ismember(VarNames,{'Survival months'});
sur_original = data_select(:,ColId);
sur_original = table2array(sur_original);

% censoring indicator (1) ----> Overall Survival
ColId = ismember(VarNames,{'Vital status recode (study cutoff used)'});
censoring_temp = data_select(:,ColId);
censoring_temp = table2cell(censoring_temp);
censoring_overall = nan(size(censoring_temp));
censoring_overall(ismember(censoring_temp,{'Alive'})) = 1; % 1: censored; 
censoring_overall(ismember(censoring_temp,{'Dead'})) = 0; % 0: uncensored

% censoring indicator (2) ----> Cancer-specific Survival
ColId = ismember(VarNames,{'SEER cause-specific death classification'});
censoring_temp = data_select(:,ColId);
censoring_temp = table2cell(censoring_temp);
censoring_specific = nan(size(censoring_temp));
censoring_specific(ismember(censoring_temp,{'Alive or dead of other cause'})) = 1; % 1: censored; 
censoring_specific(ismember(censoring_temp,{'Dead (attributable to this cancer dx)'})) = 0; % 0: uncensored

% HPV status
ColId = ismember(VarNames,{'HPV recode (2010-2017)'});
HPV_temp = data_select(:,ColId);
HPV_temp = table2cell(HPV_temp);
HPV = nan(size(HPV_temp));
HPV(ismember(HPV_temp,{'HPV Positive'})) = 1; % HPV Positive
HPV(ismember(HPV_temp,{'HPV Negative'})) = 2; % HPV Negative


%%%%%====== dataset ======%%%%%
condition1 = ~isnan(HPV);
condition2 = ~isnan(sur_original);
condition3 = ~isnan(censoring_overall);
condition4 = ~isnan(censoring_specific);
idx_select = (condition1.*condition2.*condition3.*condition4)==1;
censoring1 = censoring_overall(idx_select);
censoring2 = censoring_specific(idx_select);
sur_original = sur_original(idx_select);
HPV = HPV(idx_select);
N_temp = sum(idx_select); %
disp(['The number of included HNC patients (no missing): ',num2str(N_temp)]) % 22795
data_select = data_select(idx_select,:);


% final group
group = HPV;
groupName = {'HPV+','HPV-'};

%--- show counts ---%
tbl1 = tabulate(censoring1);
tbl1 = array2table(tbl1,'VariableNames',{'Overall censor indicator','Count','Frequency'});
disp(tbl1)
tbl2 = tabulate(censoring2);
tbl2 = array2table(tbl2,'VariableNames',{'Disease-specific censor indicator','Count','Frequency'});
disp(tbl2)
tbl3 = tabulate(group);
tbl3 = array2table(tbl3,'VariableNames',{'Group','Count','Frequency'});
disp(tbl3)

%----- Overall survival -----%
data_temp1 = [censoring1, group, sur_original];
%----- Disease-specific survival -----%
data_temp2 = [censoring2, group, sur_original];

J_group = length(unique(group));
pp = [];
for j=1:J_group
    pp = [pp,sum(group==j)];
end
pp = pp/sum(pp);

if 1==0
    N = 50/min(pp); 
else
    N = J_group*200;
end
Num_J = nan(1,J_group);
Num_J(1:(end-1)) = ceil(pp(1:(end-1)) * N);
Num_J(end) = N - sum(Num_J(1:(end-1)));
disp(['The number of analyzed HNC patients (no missing): ',num2str(N)]) % 

idx_all = [];
for j=1:J_group
    temp = find(group==j);
    idx_temp = randsample(temp,Num_J(j));
    idx_all = [idx_all;idx_temp];
end
data_select_final = data_select(idx_all,:);
filename = [folder,'original data selected_seed_',num2str(seed),'.xlsx'];
writetable(data_select_final,filename)

%----- Overall survival -----%
data1 = data_temp1(idx_all,:);
sur = 1 + data1(:,3)/max(data1(:,3));
data1 = [sur, data1];
X_max = max(data1(:,4));
%----- Disease-specific survival -----%
data2 = data_temp2(idx_all,:);
data2 = [sur, data2];

%--- show counts ---%
% censoring
tbl3 = tabulate(data1(:,2));
tbl3 = array2table(tbl3,'VariableNames',{'Overall censor indicator','Count','Frequency'});
disp(tbl3)

tbl3 = tabulate(data2(:,2));
tbl3 = array2table(tbl3,'VariableNames',{'Disease-specific censor indicator','Count','Frequency'});
disp(tbl3)

% group
tbl4 = array2table(Num_J,'VariableNames',groupName);
disp(tbl4)


% 1: survival time (1-2)
% 2: censor indicator (0: not censored; 1: censored)
% 3: group (restaurant)
% 4: original survival time
% 5: dish index
% 6: component parameter alpha
% 7: component parameter beta
% 8: table index

end
%% Function implementation


% MCMC
if strcmp(Operation,'MCMC') || strcmp(Operation,'All')
    
if strcmp(SurvivalFlag,'Overall')
    data = data1;
elseif strcmp(SurvivalFlag,'Specific')
    data = data2;
end      
    
[compute_time,data_record,gamma_record,...
    alpha_0_record,rho_record,eta_record] = HDP_WMM_MCMC(SurvivalFlag,data,...
    iterSize,burnin,gap,numDiscrete,J_group,Num_J,N,seed,folder);
disp(['Sampling computational time (',SurvivalFlag,') is ',num2str(compute_time),' hours'])

end

% Posterior inference
if strcmp(Operation,'FunEst')
    
if strcmp(SurvivalFlag,'Overall')
    filename = [folder,'results_seed_',num2str(seed),'_Overall.mat']; 
elseif strcmp(SurvivalFlag,'Specific')
    filename = [folder,'results_seed_',num2str(seed),'_Specific.mat']; 
end

load(filename,'J_group','B','seed','d','data','Num_J','N',...
    'data_record','gamma_record','alpha_0_record','rho_record','eta_record')
groupName = {'HPV+','HPV-'};

end

if strcmp(Operation,'FunEst') || strcmp(Operation,'All')
    
[FunEstAll,FunEstAll_more,xx,X_upper] = ...
    HDP_WMM_Inference(SurvivalFlag,data_record,gamma_record,alpha_0_record,...
    rho_record,eta_record,epsilon,N,B,Num_J,J_group,groupName,L,folder,seed,...
    year_extent,medianEst,CI_method,CI_coef);

%----------- Function estimates save -----------%
filename = [folder,'FunEstAll_seed_',num2str(seed),'_',SurvivalFlag,'.mat']; 
save(filename,'FunEstAll','FunEstAll_more','xx','X_upper','Num_J','J_group','groupName')

end

%%%%%%%%%%%% Functions %%%%%%%%%%%%
%% HDP_WMM_MCMC
function [computation_time,data_record,gamma_record,...
    alpha_0_record,rho_record,eta_record] = HDP_WMM_MCMC(SurvivalFlag,data,...
    iterSize,burnin,gap,numDiscrete,J_group,Num_J,N,seed,folder)
tic
B = floor((iterSize-burnin)/gap); % posterior sample size

%% Prior inference & Initialization
[parmHat,~] = wblfit(data(:,1),0.05,data(:,2));
estimated_alpha = parmHat(2); % shape v
estimated_scale = parmHat(1); % scale
estimated_beta = estimated_scale^estimated_alpha; % w

d = 2;
% eta
a_eta = 1;
b_eta = (sqrt(2) - 1) / estimated_beta;
eta = gamrnd(a_eta, 1/b_eta);
% rho
a_rho = 2;
b_rho = 4 * estimated_alpha / 3;
% rho = gprnd(1/a_rho, b_rho/a_rho, b_rho); % gprnd(k,sigma,theta)  with k=1/a; sigma=b/a, theta=b
rho = b_rho/((1-unifrnd(0,1))^(1/a_rho));
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
%     censor_temp = data(idx==i,2);
%     [parmHat,~] = wblfit(dat_temp,0.05,censor_temp);
    [parmHat,~] = wblfit(dat_temp);
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


%% MCMC sampling iterations
% storage
data_record = zeros(size(data,1),size(data,2),B);
gamma_record = zeros(1,B);
alpha_0_record = zeros(1,B);
rho_record = zeros(1,B);
eta_record = zeros(1,B);
b = 1;

%--------------------- iteration start ---------------------%
for iter = 1:iterSize %  iterSize
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Step 0: Update hyperparameters (rho & eta & alpha_0 & gamma) %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    K = max(data(:,5)); % number of dishes
    
    %--- update rho ---%
    a_rho_new = a_rho + K;
    b_rho_new = max([b_rho, max(data(:,6))]);
%     rho = gprnd(1/a_rho_new, b_rho_new/a_rho_new, b_rho_new); % gprnd(k,sigma,theta)  with k=1/a; sigma=b/a, theta=b
    rho = b_rho_new/((1-unifrnd(0,1))^(1/a_rho_new));
    
    %--- update eta ---%
    a_eta_new = a_eta + d*K;
    beta_unique = unique(data(:,7));
    b_eta_new = b_eta + sum(1./beta_unique);
    eta = gamrnd(a_eta_new, 1/b_eta_new);
     
    %--- update alpha_0 ---%
    T = sum(TableNum); % total number of tables in all restaurants
    % auxiliary variables
    uu = betarnd((alpha_0+1)*ones(1,J_group),Num_J);
    ss = binornd(1,Num_J./(Num_J+alpha_0));
    a_alpha_0_new = a_alpha_0 + T - sum(ss);
    b_alpha_0_new = b_alpha_0 - sum(log(uu));
    alpha_0 = gamrnd(a_alpha_0_new, 1/b_alpha_0_new);
        
    %--- update gamma ---%
    % auxiliary variables
    uu = betarnd(gamma+1,N);
    pp = (N*(b_gamma-log(uu))) / (N*(b_gamma-log(uu))+a_gamma+K-1);
    ss = binornd(1,pp);
    a_gamma_new = a_gamma + K - ss;
    b_gamma_new = b_gamma - log(uu);
    gamma = gamrnd(a_gamma_new, 1/b_gamma_new);
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Step 1: Update Table Index (t) %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:N
        x_i = data(i, 1); % obs
        delta_i = data(i, 2); % censor indicator
        j = data(i, 3); % restaurant
        data_j = data(data(:,3)==j,:); % all data in restaurant j
        T_j = max(data_j(:,8)); % number of tables in this restaurant
        % sample t_ji
        [t_i_new, prob_new] = tableSampler(data, i, d, eta, rho, alpha_0, J_group, gamma);
        % prob_new: probability vector of choosing a dish k for a new table
        if t_i_new > T_j % new table
            K = max(data(:,5)); % number of dishes
            k_for_tnew = randsample(K+1, 1, true, prob_new);
            if k_for_tnew > K % new dish
                [data(i, 6), data(i, 7)] = NewDishSampler_x_ji(x_i, delta_i, d, eta, rho, numDiscrete);
            else % existing dish
                idx = find(data(:,5)==k_for_tnew,1);
                data(i, 6) = data(idx,6); % alpha
                data(i, 7) = data(idx,7); % beta
            end
            data(i, 5) = k_for_tnew; % dish index 
        else % existing table 
            idx = find(data_j(:,8)==t_i_new,1);
            data(i, 5) = data_j(idx,5); % dish index for table t_i_new
            data(i, 6) = data_j(idx,6); % alpha
            data(i, 7) = data_j(idx,7); % beta
        end
        data(i, 8) = t_i_new; % table index update
    end
    
    %%% re-label the dish index %%%
    K = max(data(:,5)); % number of dishes
    Dish = unique(data(:,5)); % the current set of dish labels
    if length(Dish)<K % some labels are not used -> remove
        for k = 1:length(Dish)
            data(data(:,5)==Dish(k),5)=k;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Step 2: Update Dish Index (k) %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:J_group
        
        %%% re-label the table index %%%
        Table_j = unique(data(data(:,3)==j,8)); % the current set of table labels
        T_j = max(data(data(:,3)==j,8));
        if length(Table_j)<T_j % some labels are not used -> remove
            for t = 1:length(Table_j)
                data((data(:,3)==j) .* (data(:,8)==Table_j(t))==1,8) = t;
            end
            T_j = length(Table_j);
        end
        TableNum(j) = T_j;
        
        %%% k_jt %%%
        for t = 1:T_j
            K = max(data(:,5)); % number of dishes
            dat_jt = data((data(:,3)==j) .* (data(:,8)==t)==1,:);
            % sample k_jt
            k_jt_new = dishSampler(data, j, t, d, eta, rho, J_group, gamma);
            if k_jt_new > K % new dish
                [alpha_new, beta_new] = NewDishSampler_x_jt(dat_jt, d, eta, rho, numDiscrete);
                data((data(:,3)==j) .* (data(:,8)==t)==1,6) = alpha_new;
                data((data(:,3)==j) .* (data(:,8)==t)==1,7) = beta_new;
            else % existing dish
                data((data(:,3)==j) .* (data(:,8)==t)==1, 6) = data(find(data(:,5)==k_jt_new,1),6); % alpha
                data((data(:,3)==j) .* (data(:,8)==t)==1, 7) = data(find(data(:,5)==k_jt_new,1),7); % beta
            end
            data((data(:,3)==j) .* (data(:,8)==t)==1,5) = k_jt_new; % dish index update
        end
        
    end

    %%% re-label the dish index %%%
    K = max(data(:,5)); % number of dishes
    Dish = unique(data(:,5)); % the current set of dish labels
    if length(Dish)<K % some labels are not used -> remove
        for k = 1:length(Dish)
            data(data(:,5)==Dish(k),5)=k;
        end
        K = length(Dish);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Step 3: Update locations (alpha & beta) %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k = 1:K
        alpha_k = data(find(data(:,5)==k,1),6); % current alpha
        dat_k = data(data(:,5)==k,:); % all obs with dish k
        X_k = dat_k(:,1); % observed survival
        n_k = length(X_k); 
        X_k_o = dat_k(dat_k(:,2)==0,1); % not censored
        n_k_o = length(X_k_o); 
                
        %--- update beta_k ---%
        aa = d + n_k_o; % shape
        bb = eta + sum(X_k.^alpha_k); % scale
        beta_k_new = 1/gamrnd(aa, 1/bb);
        data(data(:,5)==k,7) = beta_k_new;
        
        %--- update alpha_k ---%
        % auxiliary variables
        if n_k_o ~= 0
%             g_0 = prod(X_k_o).^alpha_k;
            sum_log = sum(log(X_k_o));
            u_temp1 = unifrnd(0, 1);
            log_u_0 = log(u_temp1) + alpha_k * sum_log;
        else
            log_u_0 = 0;
            sum_log = 1;
        end
        alpha_lower = max([0, log_u_0/sum_log]);
        
%         g_l = exp(-(X_k.^alpha_k)./beta_k_new);
%         u_l = unifrnd(zeros(n_k,1),g_l);
%         min_temp = min(log(-beta_k_new.*log(u_l))./log(X_k));
%         alpha_upper = min([rho, min_temp]);
        
        u_temp = unifrnd(zeros(n_k,1), ones(n_k,1));
        log_u_l_item = log(X_k.^alpha_k - beta_k_new * log(1-u_temp))./log(X_k);
        alpha_upper = min([rho, min(log_u_l_item)]);
        
        % debug 
        if alpha_upper < alpha_lower
            disp('Errors: alpha_upper < alpha_lower')
            disp(['log(u_temp1): ',num2str(log(u_temp1))])
            disp(['alpha_k: ',num2str(alpha_k)])
            disp(['sum_log: ',num2str(sum_log)])
            disp(['alpha_lower: ',num2str(alpha_lower)])
            disp(['rho: ',num2str(rho)])
            disp(['min(log_u_l_item): ',num2str(min(log_u_l_item))])
            
            % return
        end
        
        
        % inverse CDF method
        r = unifrnd(0, 1);
        
        if alpha_lower~=0
            temp1 = r * ((alpha_upper/alpha_lower)^(n_k_o+1));
            temp2 = 1 - r;
            temp = (temp1 + temp2)^(1/(n_k_o+1));
            alpha_k_new = alpha_lower * temp;
        else
            alpha_k_new = (r^(1/(n_k_o+1))) * alpha_upper;
        end
        
        
        % debug 
%         if isnan(alpha_k_new)
%             disp('new alpha is nan')
%             disp(['alpha_lower: ',num2str(alpha_lower)])
%             disp(['alpha_upper: ',num2str(alpha_upper)])
%             disp(['temp1: ',num2str(temp1)])
%             disp(['temp2: ',num2str(temp2)])
%             disp(['temp: ',num2str(temp)])
%         end


%         temp1 = r * (alpha_upper^(n_k_o+1));
%         temp2 = (1 - r) * (alpha_lower^(n_k_o+1));
%         alpha_k_new = (temp1 + temp2)^(1/(n_k_o+1));
        
%         % debug
%         if alpha_k_new == inf
%             disp('new alpha is inf')
%             disp(['alpha_lower: ',num2str(alpha_lower)])
%             disp(['alpha_upper: ',num2str(alpha_upper)])
%             disp(['temp1: ',num2str(temp1)])
%             disp(['temp2: ',num2str(temp2)])
%             disp(['n_k_o: ',num2str(n_k_o)])
%         end
        
        data(data(:,5)==k,6) = alpha_k_new;
    end

 
    %%%%%%%%%%%%%%%%%%% 
    %%%%% Storage %%%%%
    %%%%%%%%%%%%%%%%%%% 
    if iter>burnin && mod(iter,gap)==0
        data_record(:,:,b) = data;
        gamma_record(b) = gamma;
        alpha_0_record(b) = alpha_0;
        rho_record(b) = rho;
        eta_record(b) = eta;
        b = b + 1;
    end
    
    %%%%%%%%%%%%%%%%%% 
    %%% print info %%%
    %%%%%%%%%%%%%%%%%%
    K = max(data(:,5)); % number of dishes
    fprintf('iter = %i, [%i] global dishes, Tables = {',iter,K);
    for j = 1:J_group
        fprintf(' [%i]',TableNum(j));
    end
    fprintf(' }. \n');
    dish1 = unique(data(:,6)); dish1 = dish1';
    dish2 = unique(data(:,7)); dish2 = dish2';
    disp(['alpha','[',num2str(length(dish1)),']: ',num2str(dish1)])
    disp(['beta','[',num2str(length(dish2)),']: ',num2str(dish2)])
end

toc
computation_time = toc;
computation_time = computation_time/60/60; % hours

%--------------- save data ---------------%
if strcmp(SurvivalFlag,'Overall')
    filename = [folder,'results_seed_',num2str(seed),'_Overall.mat']; 
elseif strcmp(SurvivalFlag,'Specific')
    filename = [folder,'results_seed_',num2str(seed),'_Specific.mat']; 
end

save(filename,'J_group','B','seed','d','data','Num_J','N',...
    'data_record','gamma_record','alpha_0_record','rho_record','eta_record')

end


%% HDP_WMM_Inference
function [FunEstAll,FunEstAll_more,xx,X_upper] = ...
    HDP_WMM_Inference(SurvivalFlag,data_record,gamma_record,alpha_0_record,...
    rho_record,eta_record,epsilon,N,B,Num_J,J_group,groupName,L,folder,seed,...
    year_extent,medianEst,CI_method,CI_coef)

% row 1: mean/median; row 2: lower bound; row 3: upper bound; 
% column 1-4: density, suvivial, hazard, cumulative hazard
FunEstAll = cell(3,4); 
% row 1: mean/median; row 2: lower bound; row 3: upper bound; 
% column 1-3: hazard ratio, survival difference; median difference
FunEstAll_more = cell(3,3);

d = 2;
fig_num = 0;
%---------------- Bayesian clustering datasets ----------------%
ClusterDraws = nan(N, B);
PDF_fitted = zeros(N, B); % density
DensityFitted = nan(N, 2); % obs, density estimate
data = data_record(:,:,1);
DensityFitted(:,1) = data(:,4);
DensityPred = nan(J_group*L, 2); % grid values, density estimate

%--------------------------------------------------------------%

X_upper = zeros(1, J_group);
xx = zeros(L, J_group);
yy = zeros(L, J_group);

for j = 1:J_group
%     x_max_temp = max(data(data(:,3)==j,4));
    X_max = max(data(:,4));
    x_max_temp = X_max;
    X_upper(j) = x_max_temp + 12*year_extent; % one more year 
    xx_temp = linspace(0,X_upper(j),L); % grid values
    xx_temp = xx_temp';
    xx(:,j) = xx_temp;
    yy(:,j) = 1 + xx_temp/X_max; 
    DensityPred((j-1)*L+1:j*L,1) = xx_temp;
end


%--------- Estimation ---------%
PDF_est = zeros(L, J_group, B); % density
CDF_est = zeros(L, J_group, B); 
Survival_est = zeros(L, J_group, B);
hazard_est = zeros(L, J_group, B);
if J_group == 2
    hazard_ratio_est = zeros(L, B); % group1/group2
    Survival_diff_est = zeros(L, B); % group1 - group2
    Survival_median_est = nan(1, B); % group1 - group2 > 0
else
    hazard_ratio_est = [];
    Survival_diff_est = [];
    Survival_median_est = [];
end
Hazard_est = zeros(L, J_group, B);  

K_est = zeros(1, B); % number of components
K_common_est = zeros(1, B); % number of components that are shared
T_est = zeros(J_group, B); % number of tables

for b = 1:B
    
    % b-th posterior samples
    data_b = data_record(:,:,b);
%     gamma_b = gamma_record(b);
    gamma_b = max(gamma_record);
    alpha_0_b = alpha_0_record(b);
    rho_b = rho_record(b);
    eta_b = eta_record(b);
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
%         temp1 = wblcdf(1, G_j_beta.^(1./G_j_alpha),G_j_alpha); % Pr(Y<1)
        temp1 = weibull_cdf(1, G_j_alpha, G_j_beta);
%         temp2 = wblcdf(2.05, G_j_beta.^(1./G_j_alpha),G_j_alpha); % Pr(Y<2)
        temp = 1 - temp1; % temp2 - temp1; 1 - temp1
        idd = temp==0;
        if sum(idd)~=0
          disp(['The number of components that have Pr(Y<1)=1: ',num2str(sum(idd))])
          disp(['The sum of the corresponding weights: ',num2str(sum(G_j_omega(idd)))])
        end
        for i = 1:L
            x_i = yy(i,j);
            % pdf
%             pdf_temp = (1/X_max) * wblpdf(x_i, G_j_beta.^(1./G_j_alpha),G_j_alpha)./temp;
            pdf_temp = (1/X_max) * weibull_pdf(x_i, G_j_alpha, G_j_beta)./temp;
            pdf_temp(idd) = 0;
            PDF_est(i,j,b) = sum(G_j_omega .* pdf_temp);
            % cdf
%             cdf_temp = (wblcdf(x_i, G_j_beta.^(1./G_j_alpha),G_j_alpha) - temp1)./temp;
            cdf_temp = (weibull_cdf(x_i, G_j_alpha,G_j_beta) - temp1)./temp;
            cdf_temp(idd) = 0;
            CDF_est(i,j,b) = sum(G_j_omega .* cdf_temp);
            % survival
            Survival_est(i,j,b) = 1 - CDF_est(i,j,b);
            % hazard
%             hazard_est(i,j,b) = sum(G_j_omega .* (pdf_temp./(1-cdf_temp)));
            hazard_est(i,j,b) = PDF_est(i,j,b)/(1-CDF_est(i,j,b));
            % Hazard
            Hazard_est(i,j,b) = -log(1-CDF_est(i,j,b));
        end
        
        for i = 1:Num_J(j)
            x_i = data_b(sum(Num_J(1:(j-1)))+i,1);
            pdf_temp = (1/X_max) * weibull_pdf(x_i, G_j_alpha, G_j_beta)./temp;
            pdf_temp(idd) = 0;
            PDF_fitted(sum(Num_J(1:(j-1)))+i,b) = sum(G_j_omega .* pdf_temp);
        end
        
    end
    
    if J_group == 2
        % hazard ratio
        hazard_ratio_est(:,b) = hazard_est(:,1,b)./hazard_est(:,2,b);
        % survival difference
        Survival_diff_est(:,b) = Survival_est(:,1,b) - Survival_est(:,2,b);
        % median survival difference
        id1 = find(Survival_est(:,1,b)-0.5<0,1,'first');
        id2 = find(Survival_est(:,2,b)-0.5<0,1,'first');
        if  isempty(id1) && ~isempty(id2)
            Survival_median_est(b) = 1;
        elseif ~isempty(id1) && ~isempty(id2)
            Survival_median_est(b) = xx(id1,1) - xx(id2,2)>0;
        end
        
    end
    
    
end

disp(['Posterior mode of the number of all dishes: ',num2str(mode(K_est))])
disp(['Posterior mode of the number of shared dishes: ',num2str(mode(K_common_est))])
disp(['Posterior probability of sharing components: ',num2str(mean(100*K_common_est>0))])


%--------------------- FunEstAll ---------------------%
%%% PDF %%%
PDF_est_mean = mean(PDF_est,3);
PDF_est_median = median(PDF_est,3);
if CI_method==1
    PDF_est_std = std(PDF_est,0,3);
    PDF_est_upper = PDF_est_mean + norminv(1 - CI_coef/2) * PDF_est_std/sqrt(B-1);
    PDF_est_lower = PDF_est_mean - norminv(1 - CI_coef/2) * PDF_est_std/sqrt(B-1);
else
    PDF_est_upper = quantile(PDF_est,CI_coef/2,3);
    PDF_est_lower = quantile(PDF_est,1 - CI_coef/2,3);
end
if medianEst == 1
    FunEstAll{1,1} = PDF_est_median;
else
    FunEstAll{1,1} = PDF_est_mean;
end
FunEstAll{2,1} = PDF_est_lower;
FunEstAll{3,1} = PDF_est_upper;

%%% Survival %%%
Survival_est_mean = mean(Survival_est,3);
Survival_est_median = median(Survival_est,3);
if CI_method==1
    Survival_est_std = std(Survival_est,0,3);
    Survival_est_upper = Survival_est_mean + norminv(1 - CI_coef/2) * Survival_est_std/sqrt(B-1);
    Survival_est_lower = Survival_est_mean - norminv(1 - CI_coef/2) * Survival_est_std/sqrt(B-1);
else
    Survival_est_upper = quantile(Survival_est,CI_coef/2,3);
    Survival_est_lower = quantile(Survival_est,1 - CI_coef/2,3);
end
if medianEst == 1
    FunEstAll{1,2} = Survival_est_median;
else
    FunEstAll{1,2} = Survival_est_mean;
end
FunEstAll{2,2} = Survival_est_lower;
FunEstAll{3,2} = Survival_est_upper;


%%% hazard %%%
hazard_est_mean = mean(hazard_est,3);
hazard_est_median = median(hazard_est,3);
if CI_method==1
    hazard_est_std = std(hazard_est,0,3);
    hazard_est_upper = hazard_est_mean + norminv(1 - CI_coef/2) * hazard_est_std/sqrt(B-1);
    hazard_est_lower = hazard_est_mean - norminv(1 - CI_coef/2) * hazard_est_std/sqrt(B-1);
else
    hazard_est_upper = quantile(hazard_est,CI_coef/2,3);
    hazard_est_lower = quantile(hazard_est,1 - CI_coef/2,3);
end
if medianEst == 1
    FunEstAll{1,3} = hazard_est_median;
else
    FunEstAll{1,3} = hazard_est_mean;
end
FunEstAll{2,3} = hazard_est_lower;
FunEstAll{3,3} = hazard_est_upper;


%%% Hazard %%%
Hazard_est_mean = mean(Hazard_est,3);
Hazard_est_median = median(Hazard_est,3);
if CI_method==1
    Hazard_est_std = std(Hazard_est,0,3);
    Hazard_est_upper = Hazard_est_mean + norminv(1 - CI_coef/2) * Hazard_est_std/sqrt(B-1);
    Hazard_est_lower = Hazard_est_mean - norminv(1 - CI_coef/2) * Hazard_est_std/sqrt(B-1);
else
    Hazard_est_upper = quantile(Hazard_est,CI_coef/2,3);
    Hazard_est_lower = quantile(Hazard_est,1 - CI_coef/2,3);
end
if medianEst == 1
    FunEstAll{1,4} = Hazard_est_median;
else
    FunEstAll{1,4} = Hazard_est_mean;
end
FunEstAll{2,4} = Hazard_est_lower;
FunEstAll{3,4} = Hazard_est_upper;

%--------------------- FunEstAll_more ---------------------%

if J_group == 2
    %%% hazard ratio %%%  
    HR_est_mean = mean(hazard_ratio_est,2);
    HR_est_median = median(hazard_ratio_est,2);
    if CI_method==1
        HR_est_std = std(hazard_ratio_est,0,2);
        HR_est_upper = HR_est_mean + norminv(1 - CI_coef/2) * HR_est_std/sqrt(B-1);
        HR_est_lower = HR_est_mean - norminv(1 - CI_coef/2) * HR_est_std/sqrt(B-1);
    else
        HR_est_upper = quantile(hazard_ratio_est,CI_coef/2,2);
        HR_est_lower = quantile(hazard_ratio_est,1 - CI_coef/2,2);
    end
    if medianEst == 1
        FunEstAll_more{1,1} = HR_est_median;
    else
        FunEstAll_more{1,1} = HR_est_mean;
    end
    FunEstAll_more{2,1} = HR_est_lower;
    FunEstAll_more{3,1} = HR_est_upper;
    
    %%% Survival difference %%%
    Survival_diff_est_mean = mean(Survival_diff_est,2);
    Survival_diff_est_median = median(Survival_diff_est,2);
    if CI_method==1
        Survival_diff_est_std = std(Survival_diff_est,0,2);
        Survival_diff_est_upper = Survival_diff_est_mean + norminv(1 - CI_coef/2) * Survival_diff_est_std/sqrt(B-1);
        Survival_diff_est_lower = Survival_diff_est_mean - norminv(1 - CI_coef/2) * Survival_diff_est_std/sqrt(B-1);
    else
        Survival_diff_est_upper = quantile(Survival_diff_est,CI_coef/2,2);
        Survival_diff_est_lower = quantile(Survival_diff_est,1 - CI_coef/2,2);
    end
    if medianEst == 1
        FunEstAll_more{1,2} = Survival_diff_est_median;
    else
        FunEstAll_more{1,2} = Survival_diff_est_mean;
    end
    FunEstAll_more{2,2} = Survival_diff_est_lower;
    FunEstAll_more{3,2} = Survival_diff_est_upper;
    
    %%% Survival median difference %%%
    FunEstAll_more{1,3} = Survival_median_est;
    
end


%%% Functions %%%
medianEst = 1; % using posterior medians for infernece
PDF_fitted_mean = mean(PDF_fitted,2);
PDF_fitted_median = median(PDF_fitted,2);
if medianEst==1
    
    for j = 1:J_group
        DensityPred(L*(j-1)+1:L*j,2) = PDF_est_median(:,j);
    end
    DensityFitted(:,2) = PDF_fitted_median;
    
else
    
    for j = 1:J_group
        DensityPred(L*(j-1)+1:L*j,2) = PDF_est_mean(:,j);
    end
    DensityFitted(:,2) = PDF_fitted_mean;
end


% check fitted estimates and the predicted estimates %
legend_name = groupName;
for j = 1:J_group
    fig_num = fig_num + 1;
    figure(fig_num);
    color_order = get(gca,'colororder');
    set(figure(fig_num),'visible','off');
    plot(DensityPred(L*(j-1)+1:L*j,1),DensityPred(L*(j-1)+1:L*j,2),'LineWidth', 2,'Color',color_order(1,:))
    xlim([0 X_upper(j)]) % xlim([1 2])
    set(gca,'FontSize',16,'FontName', 'times')
    hold on;
    scatter(DensityFitted(sum(Num_J(1:j-1))+1:sum(Num_J(1:j)),1),...
        DensityFitted(sum(Num_J(1:j-1))+1:sum(Num_J(1:j)),2)); %,'Color',color_order(1,:)
    xlabel('Time (month)')
    title(legend_name{j},'FontSize',20,'FontName', 'times', 'FontWeight','normal')
    if strcmp(SurvivalFlag,'Overall')
        saveas(gcf,[folder,'Overall_fitted_group',int2str(j),'_seed_',num2str(seed),'.png']);
    elseif strcmp(SurvivalFlag,'Specific')
        saveas(gcf,[folder,'Specific_fitted_group',int2str(j),'_seed_',num2str(seed),'.png']);
    end
    
end


ClusterDraws = ClusterDraws';
ClusterDraws = array2table(ClusterDraws);
if strcmp(SurvivalFlag,'Overall')
    filename = [folder,'ClusterDraws_HPV_Overall_seed_',num2str(seed),'.xlsx'];
elseif strcmp(SurvivalFlag,'Specific')
    filename = [folder,'ClusterDraws_HPV_Specific_seed_',num2str(seed),'.xlsx'];
end

writetable(ClusterDraws,filename,'WriteVariableNames',0)

end




