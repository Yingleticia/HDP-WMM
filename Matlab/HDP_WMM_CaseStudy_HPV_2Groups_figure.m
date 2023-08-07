%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% --  Hierarchical Dirichlet Process Weibull Mixture Model          -- %%%
%%% ----------    for prognosis analysis (cancer patients)    ---------- %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------- Figures for two outcomes (overall & specific) --------%
clear;clc;close all;

seed = 3;
site = 'Oropharynx'; % All_Site; Oropharynx;
operation = 'Overall'; % Specific; Overall; 
Fun = {'Survival','Hazard rate'};
% 'Density','Survival','Hazard rate','Cumulative hazard'
FunMore = {'Hazard ratio','Survival difference','Median difference'};
% 'Hazard ratio','Survival difference','Median difference'

EstBound = 1; % confidence interval for HDP-WMM estimation
EstBound_KM = 1; % confidence interval for KM estimation (in Survival)

folder = ['results/CaseStudy/HPV/2Groups/',site,'/seed_',num2str(seed),'/'];
disp(['---> Seed = ',num2str(seed)])
disp(['---> Selected site: ',site])
disp(['---> Selected Survival Flag: ',operation])
%% load data
xx_est = cell(1,2);
X_upper_est = cell(1,2);

%---- Overall ----%
if strcmp(operation,'Overall') 
    
    SurvivalFlag = 'Overall'; % Specific; Overall
    filename = [folder,'FunEstAll_seed_',num2str(seed),'_',SurvivalFlag,'.mat']; 
    load(filename,'FunEstAll','FunEstAll_more','xx','X_upper','Num_J','J_group','groupName')
    FunEstAll_Overall = FunEstAll;
    FunEstAll_more_Overall = FunEstAll_more;
    xx_est{1} = xx;
    X_upper_est{1} = X_upper;

    filename = [folder,'results_seed_',num2str(seed),'_',SurvivalFlag,'.mat'];
    load(filename,'data')
    censor_overall = data(:,2);
    sur_original = data(:,4);
    group = data(:,3);
    disp(['  The maximum survival month in the analyzed data is:',num2str(max(sur_original))])

end

%---- Specific ----%
if strcmp(operation,'Specific')

    SurvivalFlag = 'Specific'; % Specific; Overall
    filename = [folder,'FunEstAll_seed_',num2str(seed),'_',SurvivalFlag,'.mat']; 
    load(filename,'FunEstAll','FunEstAll_more','xx','X_upper','Num_J','J_group','groupName')
    FunEstAll_Specific = FunEstAll;
    FunEstAll_more_Specific = FunEstAll_more;
    xx_est{2} = xx;
    X_upper_est{2} = X_upper;

    filename = [folder,'results_seed_',num2str(seed),'_',SurvivalFlag,'.mat'];
    load(filename,'data')
    censor_specific = data(:,2);
    sur_original = data(:,4);
    group = data(:,3);
    disp(['  The maximum survival month in the analyzed data is:',num2str(max(sur_original))])

end
%% figure
fig_num = 1;
figure(fig_num);
color_order = get(gca,'colororder');
linestyle = {'-','-.',':','--'};
XlabelName = {'(a)','(b)','(c)','(d)','(e)'};

if strcmp(operation,'Overall') 
    
    FunNum = length(Fun);
%     tiledlayout(1,FunNum,'TileSpacing','loose'); % loose; compact; tight; none
    h = cell(1,8);
%     tiledlayout(2,ceil(FunNum/2),'TileSpacing','compact'); % loose; compact; tight; none
    xx = xx_est{1};
    X_upper = X_upper_est{1};
    
    for ff = 1:FunNum
        % row 1: mean/median; row 2: lower bound; row 3: upper bound; 
        % column 1-4: density, suvivial, hazard, cumulative hazard
        if strcmp(Fun{ff},'Density') 
            
            Fun_est = FunEstAll_Overall{1,1};
            Fun_lower = FunEstAll_Overall{2,1};
            Fun_upper = FunEstAll_Overall{3,1};            
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
%             title([Fun{ff} ' function'],'FontSize',20,'FontName', 'times', 'FontWeight','normal')
%             legend(legend_name,'FontSize',14,'FontName', 'times','Location','best')
%             legend('boxoff')
            
        elseif strcmp(Fun{ff},'Survival')
            
            Fun_est = FunEstAll_Overall{1,2};
            Fun_lower = FunEstAll_Overall{2,2};
            Fun_upper = FunEstAll_Overall{3,2};
%             nexttile(ff)
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
                censor_temp = censor_overall(group == j);
                [ss_e,x,flo,fup] = ecdf(sur_temp,'Censoring',censor_temp,'function','survivor');
                h{(j-1)*4+3} = stairs(x, ss_e,'LineWidth', 2,'Color',color_order(j,:),'linestyle',linestyle{2});
                if EstBound_KM == 1
                    h{(j-1)*4+4} = stairs(x, flo,'LineWidth', 1.5,'Color',color_order(j,:),'linestyle',linestyle{3});
                    stairs(x, fup,'LineWidth', 1.5,'Color',color_order(j,:),'linestyle',linestyle{3},'HandleVisibility','off')
                    legend_name = {'HPV+ (HDP-WMM)','HPV+ (HDP-WMM, 95% CI)',...
                        'HPV+ (KM)','HPV+ (KM, 95% CI)',...
                        'HPV- (HDP-WMM)','HPV- (HDP-WMM, 95% CI)',...
                        'HPV- (KM)','HPV- (KM, 95% CI)'};
%                     legend_name = {'HPV+ (KM)','HPV+ (KM, 95% CI)',...
%                         'HPV- (KM)','HPV- (KM, 95% CI)'};
                end
            end
            set(gca,'FontSize',16,'FontName', 'times')
            xlabel('Time (months)')
%             xlabel([XlabelName{ff},' ',Fun{ff},' function'])
%             title([Fun{ff} ' function'],'FontSize',20,'FontName', 'times', 'FontWeight','normal')
            legend(legend_name,'FontSize',14,'FontName', 'times','Location','best')
            legend('boxoff')
            
        elseif strcmp(Fun{ff},'Hazard rate')  
            
            Fun_est = FunEstAll_Overall{1,3};
            Fun_lower = FunEstAll_Overall{2,3};
            Fun_upper = FunEstAll_Overall{3,3};
            
            fig_num = fig_num + 1;
            figure(fig_num);
%             nexttile(ff)
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
%             xlabel([XlabelName{ff},' ',Fun{ff},' function'])
%             title([Fun{ff} ' function'],'FontSize',20,'FontName', 'times', 'FontWeight','normal')
            legend(legend_name,'FontSize',14,'FontName', 'times','Location','best')
            legend('boxoff')
            
        elseif strcmp(Fun{ff},'Cumulative hazard') 
            
            Fun_est = FunEstAll_Overall{1,4};
            Fun_lower = FunEstAll_Overall{2,4};
            Fun_upper = FunEstAll_Overall{3,4};
            nexttile(ff)
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
            xlabel([XlabelName{ff},' ',Fun{ff},' function'])
%             title([Fun{ff} ' function'],'FontSize',20,'FontName', 'times', 'FontWeight','normal')
%             legend(legend_name,'FontSize',14,'FontName', 'times','Location','best')
%             legend('boxoff')
        end
        
        
    end
%     legend_name = {'HPV+ (HDP-WMM)','HPV+ (HDP-WMM, 95% CI)',...
%         'HPV+ (KM)','HPV+ (KM, 95% CI)',...
%         'HPV- (HDP-WMM)','HPV- (HDP-WMM, 95% CI)',...
%         'HPV- (KM)','HPV- (KM, 95% CI)'};
%     hL = legend([h{1},h{2},h{3},h{4},h{5},h{6},h{7},h{8}]); 
%     hL.Layout.Tile = 3; % 'east'; 1-3
%     hL.String = legend_name;
%     hL.Box = 'off';
%     hL.Location = 'best';
    
    FunMoreNum = length(FunMore);
    if FunMoreNum>0
        fig_num = fig_num + 1;
        figure(fig_num);
%         tiledlayout(1,FunMoreNum-1,'TileSpacing','loose'); % loose; compact; tight; none
        %------ hazard ratio ------%
        % row 1: mean/median; row 2: lower bound; row 3: upper bound; 
        % column 1-3: hazard ratio, survival difference; median difference (1 row)
        Fun_est = FunEstAll_more_Overall{1,1};
        Fun_lower = FunEstAll_more_Overall{2,1};
        Fun_upper = FunEstAll_more_Overall{3,1};
%         nexttile(1)
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
%             legend_name = {['(',groupName{1},')','/','(',groupName{2},')'],'95% CI'};
        else
            legend_name = {['(',groupName{1},')','/','(',groupName{2},')']};
        end 
        xlim([0 max(X_upper)]);
        set(gca,'FontSize',16,'FontName', 'times')
        xlabel('Time (months)')
%         xlabel([XlabelName{1},' Hazard ratio']) % between HPV+ and HPV-
%         title('Hazard ratio','FontSize',20,'FontName', 'times', 'FontWeight','normal')
        legend(legend_name,'FontSize',14,'FontName', 'times','Location','southeast')
        legend('boxoff')
        
        %------ survival probability difference ------%
        Fun_est = FunEstAll_more_Overall{1,2};
        Fun_lower = FunEstAll_more_Overall{2,2};
        Fun_upper = FunEstAll_more_Overall{3,2};
        fig_num = fig_num + 1;
        figure(fig_num);
%         nexttile(2)
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
%             legend_name = {['(',groupName{1},')','$-$','(',groupName{2},')'],'95% CI'};
        else
            legend_name = {['(',groupName{1},')','$-$','(',groupName{2},')']};
        end 
        xlim([0 max(X_upper)]);
        set(gca,'FontSize',16,'FontName', 'times')
        xlabel('Time (months)')
%         xlabel([XlabelName{2},' Difference of survival functions']) % between HPV+ and HPV-
%         title('Difference of survival functions','FontSize',20,'FontName', 'times', 'FontWeight','normal')
        legend(legend_name,'FontSize',14,'FontName', 'times','Location','southeast','Interpreter','latex')
        legend('boxoff')
        
        %------ survival median difference ------%
        temp = FunEstAll_more_Overall{1,3};
        temp = temp(~isnan(temp));
        temp = 100 * sum(temp)/length(temp);
        disp(['===> Posterior probability of survival median difference > 0: ',num2str(temp)])
    end
    
    
    
end

if strcmp(operation,'Specific') 
    
    Survival_est_temp = Survival_est_all{2};
    xx = xx_est{2};
    X_upper = X_upper_est{2};
    for j = 1:J_group
        plot(xx(:,j),Survival_est_temp(:,j),'LineWidth', 2,'Color',color_order(j,:),'linestyle',linestyle{1}); hold on; % 
        xlim([0 max(X_upper)]); ylim([0 1]); 
        sur_temp = sur_original(group == j);
        censor_temp = censor_specific(group == j);
        [ss_e,x,flo,fup] = ecdf(sur_temp,'Censoring',censor_temp,'function','survivor');
        stairs(x, ss_e,'LineWidth', 2,'Color',color_order(j,:),'linestyle',linestyle{2})
%         plot(x, flo,'LineWidth', 1.5,'Color',color_order(j,:),'linestyle',linestyle{3})
%         plot(x, fup,'LineWidth', 1.5,'Color',color_order(j,:),'linestyle',linestyle{3})
    end
    set(gca,'FontSize',16,'FontName', 'times')
    xlabel('Time (month)')
    title('Disease-specific survival function','FontSize',20,'FontName', 'times', 'FontWeight','normal')
    legend(legend_name,'FontSize',14,'FontName', 'times','Location','best')
    legend('boxoff')
    
end



%% figure saving
if 1==0
    
name1 = ['FunEstAll_',operation,'_seed_',num2str(seed),'.png'];
name2 = ['FunEstAll_',operation,'_seed_',num2str(seed),'.fig'];

saveas(gcf,[folder,name1]);
saveas(gcf,[folder,name2]);

end

