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
% 3: Censoring cases (keep the original one or set to fixed value)?
% 4: Computation time
% 5: only used for MCMC; move FunEst to 'figure' file


clear;clc;close all;

%================== seed =================%
seed = 3;       
rng(seed);
RandStream.getGlobalStream;

%============ Problem Setting ============%
SurvivalFlag = 'Overall'; % Specific; Overall
DataSelectFlag = 'Original'; % Equal (200/group); Original (ratio) among groups
CensorFlag = 'Original'; % Fixed; Original (ratio)

ParaInitialFlag = 'Censor'; % Censor; NoCensor
% Censor: including censor indicator for alpha beta initialization
% NoCensor: not including

if strcmp(CensorFlag,'Original')
    ParaInitialFlag = 'NoCensor'; 
elseif strcmp(CensorFlag,'Fixed')
    censor_rate = 0.3;
end
%=========================================%

folder = ['results/CaseStudy/seed_',num2str(seed),'/'];
disp(['---> Seed = ',num2str(seed)])
disp(['---> Selected Survival Flag: ',SurvivalFlag])

%% Parameter setting

% data
Num_jj = 200;

% MCMC
iterSize = 10000;
burnin = 5000;
gap = 20;
B = floor((iterSize-burnin)/gap); % posterior sample size

% model parameter
epsilon = 0.001; % for DP approximation
L = 200; % for grid values in function inference
numDiscrete = 200; % for discretization sampling


%% Data 

%%%%%%%%%%%%%%%%%
%%% Real Data %%%
%%%%%%%%%%%%%%%%%

%------------ Data Extraction ------------%

load('CaseStudyData/HNC_HPV_08222022.mat') % 'data_HNC_HPV' 
N_total = size(data_HNC_HPV, 1); % 22888
disp(['The number of all HNC patients with known HPV status: ',num2str(N_total)])
VarNames = data_HNC_HPV.Properties.VariableNames;

%%%%%====== data selection ======%%%%%
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
disp(['The number of Oropharynx patients: ',num2str(size(data_select,1))]) % 11336


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
disp(['The number of included HNC patients (no missing): ',num2str(N_temp)]) % 11292
data_select = data_select(idx_select,:);

% final group
group = HPV;
groupName = {'HPV+','HPV-'};

%--- show counts ---%
tbl3 = tabulate(group);
tbl3 = array2table(tbl3(:,2:3),'VariableNames',{'Count','Frequency'});
tbl3 = [cell2table(groupName','VariableNames',{'Group'}),tbl3];
disp(tbl3)
if strcmp(SurvivalFlag,'Overall')
    tbl1 = tabulate(censoring1);
    tbl1 = array2table(tbl1,'VariableNames',{'Overall censor indicator','Count','Frequency'});
    disp(tbl1)
    tbl4 = tabulate(censoring1(group == 1));
    tbl4 = array2table(tbl4,'VariableNames',{'Overall censor indicator (HPV+)','Count','Frequency'});
    disp(tbl4)
    tbl5 = tabulate(censoring1(group == 2));
    tbl5 = array2table(tbl5,'VariableNames',{'Overall censor indicator (HPV-)','Count','Frequency'});
    disp(tbl5)
elseif strcmp(SurvivalFlag,'Specific')
    tbl2 = tabulate(censoring2);
    tbl2 = array2table(tbl2,'VariableNames',{'Disease-specific censor indicator','Count','Frequency'});
    disp(tbl2)
    tbl6 = tabulate(censoring2(group == 1));
    tbl6 = array2table(tbl6,'VariableNames',{'Disease-specific censor indicator (HPV+)','Count','Frequency'});
    disp(tbl6)
    tbl7 = tabulate(censoring2(group == 2));
    tbl7 = array2table(tbl7,'VariableNames',{'Disease-specific censor indicator (HPV-)','Count','Frequency'});
    disp(tbl7)
end 

%----- Overall survival -----%
data_temp1 = [censoring1, group, sur_original];
%----- Disease-specific survival -----%
data_temp2 = [censoring2, group, sur_original];

J_group = length(unique(group));
if strcmp(DataSelectFlag,'Equal') % Equal; Original (ratio)
    Num_J = Num_jj * ones(1,J_group);
    N = sum(Num_J);
else
    N = J_group*Num_jj;
    pp = zeros(1,J_group);
    for j=1:J_group
        pp(j) = sum(group==j);
    end
    pp = pp/sum(pp);
    Num_J = nan(1,J_group);
    Num_J(1:(end-1)) = ceil(pp(1:(end-1)) * N);
    Num_J(end) = N - sum(Num_J(1:(end-1)));
end
disp(['The number of analyzed HNC patients (no missing): ',num2str(N)]) % 

idx_all = [];
if strcmp(CensorFlag,'Original') % Fixed (censor_rate); Original (ratio)
    for j=1:J_group
        temp = find(group==j);
        idx_temp = randsample(temp,Num_J(j));
        idx_all = [idx_all;idx_temp];
    end
else
    if strcmp(SurvivalFlag,'Overall')
        censor = censoring1;
    elseif strcmp(SurvivalFlag,'Specific')
        censor = censoring2;
    end 
    for j=1:J_group
        temp0 = find((group==j).*(censor==0)==1); % 0: uncensored;
        temp1 = find((group==j).*(censor==1)==1); % 1: censored;
        Num1 = floor(Num_J(j)*censor_rate);
        Num0 = Num_J(j) - Num1;
        idx_temp0 = randsample(temp0,Num0);
        idx_temp1 = randsample(temp1,Num1);
        idx_temp = [idx_temp0;idx_temp1];
        idx_all = [idx_all;idx_temp];
    end
end
data_select_final = data_select(idx_all,:);
if strcmp(CensorFlag,'Fixed')
    filename = [folder,'data_censor_',CensorFlag,'_',num2str(censor_rate*10),'_seed_',num2str(seed),'.xlsx'];
else
    filename = [folder,'data_censor_',CensorFlag,'_seed_',num2str(seed),'.xlsx'];
end
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
% group
tbl4 = array2table(Num_J,'VariableNames',groupName);
disp(tbl4)

% censoring
if strcmp(SurvivalFlag,'Overall')
    tbl3 = tabulate(data1(:,2));
    tbl3 = array2table(tbl3,'VariableNames',{'Overall censor indicator','Count','Frequency'});
    disp(tbl3)
elseif strcmp(SurvivalFlag,'Specific')
    tbl3 = tabulate(data2(:,2));
    tbl3 = array2table(tbl3,'VariableNames',{'Disease-specific censor indicator','Count','Frequency'});
    disp(tbl3)
end 

% 1: survival time (1-2)
% 2: censor indicator (0: not censored; 1: censored)
% 3: group (restaurant)
% 4: original survival time
% 5: dish index
% 6: component parameter alpha
% 7: component parameter beta
% 8: table index

%% Function implementation

if strcmp(SurvivalFlag,'Overall')
    data = data1;
elseif strcmp(SurvivalFlag,'Specific')
    data = data2;
end      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   Prior inference & Initialization  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ParaInitialFlag,'Censor')
    [parmHat,~] = wblfit(data(:,1),0.05,data(:,2));
elseif strcmp(ParaInitialFlag,'NoCensor')
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
alpha_kk = zeros(size(idx));
beta_kk = zeros(size(idx));
data = [data,idx,alpha_kk,beta_kk];
for i=1:kk
    dat_temp = data(idx==i,1);
    if strcmp(ParaInitialFlag,'Censor')
        censor_temp = data(idx==i,2);
        [parmHat,~] = wblfit(dat_temp,0.05,censor_temp);
    elseif strcmp(ParaInitialFlag,'NoCensor')
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
if strcmp(CensorFlag,'Fixed')
    filename = [folder,'HDP_seed_',num2str(seed),'_censor_',CensorFlag,'_',num2str(censor_rate*10),'.mat'];
    save(filename,'L','X_max','J_group','B',...
        'seed','d','censor_rate','data','Num_J','N',...
        'data_record_HDP','gamma_record_HDP','alpha_0_record_HDP',...
        'rho_record_HDP','eta_record_HDP','computation_time')
else
    filename = [folder,'HDP_seed_',num2str(seed),'_censor_',CensorFlag,'.mat'];
    save(filename,'L','X_max','J_group','B',...
        'seed','d','data','Num_J','N',...
        'data_record_HDP','gamma_record_HDP','alpha_0_record_HDP',...
        'rho_record_HDP','eta_record_HDP','computation_time')
end



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
if strcmp(CensorFlag,'Fixed')
    filename = [folder,'DP_all_seed_',num2str(seed),'_censor_',CensorFlag,'_',num2str(censor_rate*10),'.mat'];
    save(filename,'L','X_max','J_group','B',...
        'seed','d','censor_rate','data','Num_J','N',...
        'data_record_DP_all','gamma_record_DP_all','rho_record_DP_all',...
        'eta_record_DP_all','computation_time')
else
    filename = [folder,'DP_all_seed_',num2str(seed),'_censor_',CensorFlag,'.mat'];
    save(filename,'L','X_max','J_group','B',...
        'seed','d','data','Num_J','N',...
        'data_record_DP_all','gamma_record_DP_all','rho_record_DP_all',...
        'eta_record_DP_all','computation_time')
end



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
if strcmp(CensorFlag,'Fixed')
    filename = [folder,'DP_each_seed_',num2str(seed),'_censor_',CensorFlag,'_',num2str(censor_rate*10),'.mat'];
    save(filename,'L','X_max','J_group','B',...
        'seed','d','censor_rate','data','Num_J','N',...
        'results_all','computation_time')
else
    filename = [folder,'DP_each_seed_',num2str(seed),'_censor_',CensorFlag,'.mat'];
    save(filename,'L','X_max','J_group','B',...
        'seed','d','data','Num_J','N',...
        'results_all','computation_time')
end






