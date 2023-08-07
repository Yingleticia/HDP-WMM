%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% --  Hierarchical Dirichlet Process Weibull Mixture Model          -- %%%
%%% ----------    for prognosis analysis (cancer patients)    ---------- %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------- Figures for two outcomes (overall & specific) --------%
clear;clc;close all;
tic
%======================= Parameter setting =======================%
seed = 3;

Operation = 'metric'; % figure; metric

SurvivalFlag = 'Overall'; % Specific; Overall; 
CensorFlag = 'Original'; % Fixed; Original (ratio)
if strcmp(CensorFlag,'Original')
    ParaInitialFlag = 'NoCensor'; 
elseif strcmp(CensorFlag,'Fixed')
    censor_rate = 0.3;
end

Fun = {'Survival','Hazard rate'};
% 'Density','Survival','Hazard rate','Cumulative hazard'
FunMore = {'Hazard ratio','Survival difference','Median difference'};
% 'Hazard ratio','Survival difference','Median difference'

EstBound = 1; % confidence interval for HDP-WMM estimation
EstBound_KM = 0; % confidence interval for KM estimation (in Survival)
epsilon = 0.001; % for DP approximation
CI_coef = 0.05;
year_extent = 5; % extented range (years)
medianEst = 0; 
if medianEst == 0 % posterior mean
    CI_method = 1; % 1: t-distribution
elseif medianEst == 1 % posterior median
    CI_method = 2; % 2: quantile
end

%=================================================================%

folder = ['results/CaseStudy/seed_',num2str(seed),'/'];
disp(['---> Seed = ',num2str(seed)])
disp(['---> Selected Survival Flag: ',SurvivalFlag])

%% load data 

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HDP-WMM (proposed) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(CensorFlag,'Fixed')
    filename = [folder,'HDP_seed_',num2str(seed),'_censor_',CensorFlag,'_',num2str(censor_rate*10),'.mat'];
else
    filename = [folder,'HDP_seed_',num2str(seed),'_censor_',CensorFlag,'.mat'];
end
% 'L','X_max','J_group','B',...
% 'seed','d','censor_rate','data','Num_J','N',...
% 'data_record_HDP','gamma_record_HDP','alpha_0_record_HDP',...
% 'rho_record_HDP','eta_record_HDP','computation_time')
load(filename)
disp(['    HDP-WMM (proposed): Sampling computation time is ',num2str(computation_time),' hours'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  DP-WMM (all data)  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(CensorFlag,'Fixed')
    filename = [folder,'DP_all_seed_',num2str(seed),'_censor_',CensorFlag,'_',num2str(censor_rate*10),'.mat'];
else
    filename = [folder,'DP_all_seed_',num2str(seed),'_censor_',CensorFlag,'.mat'];
end
% 'L','X_max','J_group','B',...
% 'seed','d','censor_rate','data','Num_J','N',...
% 'data_record_DP_all','gamma_record_DP_all','rho_record_DP_all',...
% 'eta_record_DP_all','computation_time'
load(filename)
disp(['    DP-WMM (All data): Sampling computation time is ',num2str(computation_time),' hours'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DP-WMM (each group) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(CensorFlag,'Fixed')
    filename = [folder,'DP_each_seed_',num2str(seed),'_censor_',CensorFlag,'_',num2str(censor_rate*10),'.mat'];
else
    filename = [folder,'DP_each_seed_',num2str(seed),'_censor_',CensorFlag,'.mat'];
end
% 'L','X_max','J_group','B',...
% 'seed','d','censor_rate','data','Num_J','N',...
% 'results_all','computation_time'
load(filename)
disp(['    DP-WMM (Each group): Sampling computation time is ',num2str(computation_time),' hours'])

%% function estimate
censor = data(:,2);
sur_original = data(:,4);
group = data(:,3);
disp(['       The maximum survival month in the analyzed data is:',num2str(max(sur_original))])

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HDP-WMM (proposed) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
[HDP_FunEstAll,HDP_FunEstAll_more,HDP_xx,HDP_X_upper] = R1_HDP_WMM_FunEst(SurvivalFlag,...
    data_record_HDP,gamma_record_HDP,alpha_0_record_HDP,...
    rho_record_HDP,eta_record_HDP,epsilon,N,B,J_group,L,folder,...
    year_extent,medianEst,CI_method,CI_coef,Operation);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  DP-WMM (all data)  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
[DP_all_FunEstAll,DP_all_xx,DP_all_X_upper] = ...
    R1_DP_WMM_All_FunEst(SurvivalFlag,data_record_DP_all,gamma_record_DP_all,...
    rho_record_DP_all,eta_record_DP_all,epsilon,N,B,J_group,L,folder,...
    year_extent,medianEst,CI_method,CI_coef,Operation);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DP-WMM (each group) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
[DP_Each_FunEstAll,DP_Each_FunEstAll_more,DP_Each_xx,DP_Each_X_upper] = ...
    R1_DP_WMM_Each_FunEst(SurvivalFlag,data,results_all,...
    epsilon,N,B,J_group,L,folder,Num_J,...
    year_extent,medianEst,CI_method,CI_coef,Operation);



%% figure
if strcmp(Operation,'figure')   

fig_num = 1;
figure(fig_num);
color_order = get(gca,'colororder');
linestyle = {'-','-.',':','--'};
XlabelName = {'(a)','(b)','(c)','(d)','(e)'};


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HDP-WMM (proposed) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

FunNum = length(Fun);
h = cell(1,8);
xx = HDP_xx;
X_upper = HDP_X_upper;
%------------- estimated functions -------------%
for ff = 1:FunNum
    % row 1: mean/median; row 2: lower bound; row 3: upper bound; 
    % column 1-4: density, suvivial, hazard, cumulative hazard
    if strcmp(Fun{ff},'Density') 
        Fun_est = HDP_FunEstAll{1,1};
        Fun_lower = HDP_FunEstAll{2,1};
        Fun_upper = HDP_FunEstAll{3,1}; 
        nexttile(ff)
        for j = 1:J_group
            plot(xx(:,j),Fun_est(:,j),'LineWidth', 2,'Color',...
                color_order(j,:),'linestyle',linestyle{1}); hold on; % 
            if EstBound == 1 % confidence interval
                x = xx(:,j);
                p_low = Fun_lower(:,j);
                p_up = Fun_upper(:,j);
                x = x';
                p_low = p_low';
                p_up = p_up';
                h_temp=fill([x,fliplr(x)],[p_low,fliplr(p_up)],color_order(j,:));
                set(h_temp,'edgealpha',0,'facealpha',0.1)
                legend_name = {'HPV+ (HDP-WMM)','HPV+ (HDP-WMM, 95% CI)',...
                    'HPV- (HDP-WMM)','HPV- (HDP-WMM, 95% CI)'};
            else
                legend_name = {'HPV+ (HDP-WMM)','HPV- (HDP-WMM)'};
            end
            xlim([0 max(X_upper)]); % ylim([0 1]);
        end
        set(gca,'FontSize',16,'FontName', 'times')
        xlabel([XlabelName{ff},' ',Fun{ff},' function'])
        
    elseif strcmp(Fun{ff},'Survival')
        Fun_est = HDP_FunEstAll{1,2};
        Fun_lower = HDP_FunEstAll{2,2};
        Fun_upper = HDP_FunEstAll{3,2};
        for j = 1:J_group
            h{(j-1)*4+1} = plot(xx(:,j),Fun_est(:,j),'LineWidth', 2,'Color',color_order(j,:),...
                'linestyle',linestyle{1}); hold on; % ,'HandleVisibility','off'
            if EstBound == 1 % confidence interval
                x = xx(:,j);
                p_low = Fun_lower(:,j);
                p_up = Fun_upper(:,j);
                x = x';
                p_low = p_low';
                p_up = p_up';
                h{(j-1)*4+2} = fill([x,fliplr(x)],[p_low,fliplr(p_up)],color_order(j,:));
                set(h{(j-1)*4+2},'edgealpha',0,'facealpha',0.1)
                legend_name = {'HPV+ (HDP-WMM)','HPV+ (HDP-WMM, 95% CI)',...
                    'HPV+ (KM)','HPV- (HDP-WMM)',...
                    'HPV- (HDP-WMM, 95% CI)','HPV- (KM)'};
            else
                legend_name = {'HPV+ (HDP-WMM)','HPV+ (KM)','HPV- (HDP-WMM)','HPV- (KM)'};
            end
            xlim([0 max(X_upper)]); ylim([0 1]);
            sur_temp = sur_original(group == j);
            censor_temp = censor(group == j);
            [ss_e,x,flo,fup] = ecdf(sur_temp,'Censoring',censor_temp,'function','survivor');
            h{(j-1)*4+3} = stairs(x, ss_e,'LineWidth', 2,'Color',color_order(j,:),'linestyle',linestyle{2});
            if EstBound_KM == 1
                h{(j-1)*4+4} = stairs(x, flo,'LineWidth', 1.5,'Color',color_order(j,:),'linestyle',linestyle{3});
                stairs(x, fup,'LineWidth', 1.5,'Color',color_order(j,:),'linestyle',linestyle{3},'HandleVisibility','off')
                legend_name = {'HPV+ (HDP-WMM)','HPV+ (HDP-WMM, 95% CI)',...
                    'HPV+ (KM)','HPV+ (KM, 95% CI)',...
                    'HPV- (HDP-WMM)','HPV- (HDP-WMM, 95% CI)',...
                    'HPV- (KM)','HPV- (KM, 95% CI)'};
            end
        end
        set(gca,'FontSize',16,'FontName', 'times')
        xlabel('Time (months)')
        legend(legend_name,'FontSize',14,'FontName', 'times','Location','best')
        legend('boxoff')
        
    elseif strcmp(Fun{ff},'Hazard rate')
        Fun_est = HDP_FunEstAll{1,3};
        Fun_lower = HDP_FunEstAll{2,3};
        Fun_upper = HDP_FunEstAll{3,3};
        fig_num = fig_num + 1;
        figure(fig_num);
        for j = 1:J_group
            plot(xx(:,j),Fun_est(:,j),'LineWidth', 2,'Color',color_order(j,:),...
                'linestyle',linestyle{1}); hold on; % 
            if EstBound == 1 % confidence interval
                x = xx(:,j);
                p_low = Fun_lower(:,j);
                p_up = Fun_upper(:,j);
                x = x';
                p_low = p_low';
                p_up = p_up';
                h_temp=fill([x,fliplr(x)],[p_low,fliplr(p_up)],color_order(j,:));
                set(h_temp,'edgealpha',0,'facealpha',0.1)
                legend_name = {'HPV+ (HDP-WMM)','HPV+ (HDP-WMM, 95% CI)',...
                    'HPV- (HDP-WMM)','HPV- (HDP-WMM, 95% CI)'};
            else
                legend_name = {'HPV+ (HDP-WMM)','HPV- (HDP-WMM)'};
            end
            xlim([0 max(X_upper)]);
        end
        set(gca,'FontSize',16,'FontName', 'times')
        xlabel('Time (months)')
        legend(legend_name,'FontSize',14,'FontName', 'times','Location','best')
        legend('boxoff')
        
    end
end

%------------- group differences -------------%    
FunMoreNum = length(FunMore);
if FunMoreNum>0
    fig_num = fig_num + 1;
    figure(fig_num);
    %------ hazard ratio ------%
    % row 1: mean/median; row 2: lower bound; row 3: upper bound; 
    % column 1-3: hazard ratio, survival difference; median difference (1 row)
    Fun_est = HDP_FunEstAll_more{1,1};
    Fun_lower = HDP_FunEstAll_more{2,1};
    Fun_upper = HDP_FunEstAll_more{3,1};
    plot(xx(:,1),Fun_est,'LineWidth', 2,'Color',color_order(3,:),...
            'linestyle',linestyle{1}); 
    if EstBound == 1 % confidence interval
        hold on;
        x = xx(:,1);
        p_low = Fun_lower;
        p_up = Fun_upper;
        x = x';
        p_low = p_low';
        p_up = p_up';
        h_temp=fill([x,fliplr(x)],[p_low,fliplr(p_up)],color_order(3,:));
        set(h_temp,'edgealpha',0,'facealpha',0.1)
        legend_name = {'Posterior means','95% CI'};
    else
        legend_name = {['(',groupName{1},')','/','(',groupName{2},')']};
    end 
    xlim([0 max(X_upper)]);
    set(gca,'FontSize',16,'FontName', 'times')
    xlabel('Time (months)')
    legend(legend_name,'FontSize',14,'FontName', 'times','Location','southeast')
    legend('boxoff')

    %------ survival probability difference ------%
    Fun_est = HDP_FunEstAll_more{1,2};
    Fun_lower = HDP_FunEstAll_more{2,2};
    Fun_upper = HDP_FunEstAll_more{3,2};
    fig_num = fig_num + 1;
    figure(fig_num);
    plot(xx(:,1),Fun_est,'LineWidth', 2,'Color',color_order(3,:),...
            'linestyle',linestyle{1});
    if EstBound == 1 % confidence interval
        hold on;
        x = xx(:,1);
        p_low = Fun_lower;
        p_up = Fun_upper;
        x = x';
        p_low = p_low';
        p_up = p_up';
        h_temp=fill([x,fliplr(x)],[p_low,fliplr(p_up)],color_order(3,:));
        set(h_temp,'edgealpha',0,'facealpha',0.1)
        legend_name = {'Posterior means','95% CI'};
    else
        legend_name = {['(',groupName{1},')','$-$','(',groupName{2},')']};
    end 
    xlim([0 max(X_upper)]);
    set(gca,'FontSize',16,'FontName', 'times')
    xlabel('Time (months)')
    legend(legend_name,'FontSize',14,'FontName', 'times','Location','southeast','Interpreter','latex')
    legend('boxoff')

    %------ survival median difference ------%
    temp = HDP_FunEstAll_more{1,3};
    temp = temp(~isnan(temp));
    temp = 100 * sum(temp)/length(temp);
    disp(['===> Posterior probability of survival median difference > 0: ',num2str(temp)])
end


%% Metrics: 
elseif strcmp(Operation,'metric')
    
%------------ Model fitting ------------%
% functions: survival
% fitness: KM vs estimation (MSE)
figure;
color_order = get(gca,'colororder');
linestyle = {':','-','-.','--'};
if EstBound_KM ==1
    legend_name = {'HPV+(KM)','HPV+ (KM, 95% CI)',...
        'HPV+(HDP-WMM)','HPV+(DM-WMM-All)','HPV+(DP-WMM-Each)',...
        'HPV-(KM)','HPV- (KM, 95% CI)',...
        'HPV-(HDP-WMM)','HPV-(DM-WMM-All)','HPV-(DP-WMM-Each)'};
else
    legend_name = {'HPV+(KM)','HPV+(HDP-WMM)','HPV+(DM-WMM-All)','HPV+(DP-WMM-Each)',...
        'HPV-(KM)','HPV-(HDP-WMM)','HPV-(DM-WMM-All)','HPV-(DP-WMM-Each)'};
end
fontSZ = 20;

MSE_fit = zeros(J_group, 3); % Survival: Groups X Model 
MAE_fit = zeros(J_group, 3); % Survival: Groups X Model

for j = 1:J_group
    
    %--- KM ---%
    sur_temp = sur_original(group == j);
    censor_temp = censor(group == j);
    [KM_sur,KM_x,KM_lo,KM_up] = ecdf(sur_temp,'Censoring',censor_temp,'function','survivor');
    L_j = length(KM_x);
    
    %--- HDP-WMM ---%
    % row 1: mean/median; row 2: lower bound; row 3: upper bound; 
    % column 1-4: density, suvivial, hazard, cumulative hazard
    Fun_est = HDP_FunEstAll{1,2};
    HDP_sur = Fun_est(:,j);
    HDP_sur = HDP_sur(1:L_j);
    MSE_fit(j,1) = mean((KM_sur-HDP_sur).^2);
    MAE_fit(j,1) = mean(abs(KM_sur-HDP_sur));
    
    %--- DP-WMM-All ---%
    Fun_est = DP_all_FunEstAll{1,2};
    DP_all_sur = Fun_est(:,j);
    DP_all_sur = DP_all_sur(1:L_j);
    MSE_fit(j,2) = mean((KM_sur-DP_all_sur).^2);
    MAE_fit(j,2) = mean(abs(KM_sur-DP_all_sur));
    
    %--- DP-WMM-Each ---%
    Fun_est = DP_Each_FunEstAll{1,2};
    DP_each_sur = Fun_est(:,j);
    DP_each_sur = DP_each_sur(1:L_j);
    MSE_fit(j,3) = mean((KM_sur-DP_each_sur).^2);
    MAE_fit(j,3) = mean(abs(KM_sur-DP_each_sur));
        
    
    if 1 == 1 % figure for survival 
        KM_x = KM_x';
        KM_sur = KM_sur';
        KM_lo = KM_lo';
        KM_up = KM_up';
        HDP_sur = HDP_sur';
        DP_all_sur = DP_all_sur';
        DP_each_sur = DP_each_sur';
        stairs(KM_x, KM_sur,'LineWidth', 2,'Color',color_order(j,:),...
            'linestyle',linestyle{1});  hold on; 
        if EstBound_KM == 1
            h = fill([KM_x(2:end),fliplr(KM_x(2:end))],[KM_lo(2:end),fliplr(KM_up(2:end))],color_order(j,:));
            set(h,'edgealpha',0,'facealpha',0.1)
        end
        plot(KM_x,HDP_sur,'LineWidth', 2,'Color',color_order(j,:),...
                'linestyle',linestyle{2});
        plot(KM_x,DP_all_sur,'LineWidth', 2,'Color',color_order(j,:),...
                'linestyle',linestyle{3});
        plot(KM_x,DP_each_sur,'LineWidth', 2,'Color',color_order(j,:),...
                'linestyle',linestyle{4});
    end
    
end
ylim([0,1])
lg = legend(legend_name,'Orientation','Vertical','FontSize',fontSZ-2,...
    'FontName', 'times','NumColumns',2,'Location','southwest'); % 'southwest';'eastoutside'
legend('boxoff')
set(gca,'FontSize',fontSZ,'FontName', 'times')
xlabel('Time (months)','FontSize',fontSZ,'FontName', 'times')
ylabel('Survival probability','FontSize',fontSZ,'FontName', 'times')



disp('===================================== MSE (fitting) =====================================')
tbl_sum = zeros(1,3);
for j = 1:J_group
    disp(['-------------------------- Group ',num2str(j),' --------------------------'])
    tbl = MSE_fit(j,:);
    tbl_sum = tbl_sum + tbl;
    tbl = array2table(tbl,'VariableNames',{'HDP-WMM','DP-WMM-All','DP-WMM-Each'});
    disp(tbl)
end

disp('-------------------------- All Groups --------------------------')
tbl_sum = array2table(tbl_sum,'VariableNames',{'HDP-WMM','DP-WMM-All','DP-WMM-Each'});
disp(tbl_sum)




disp('===================================== MAE (fitting) =====================================')
tbl_sum = zeros(1,3);
for j = 1:J_group
    disp(['-------------------------- Group ',num2str(j),' --------------------------'])
    tbl = MAE_fit(j,:);
    tbl_sum = tbl_sum + tbl;
    tbl = array2table(tbl,'VariableNames',{'HDP-WMM','DP-WMM-All','DP-WMM-Each'});
    disp(tbl)
end

disp('-------------------------- All Groups --------------------------')
tbl_sum = array2table(tbl_sum,'VariableNames',{'HDP-WMM','DP-WMM-All','DP-WMM-Each'});
disp(tbl_sum)




end


toc




