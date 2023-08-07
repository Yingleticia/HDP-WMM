%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% --  Hierarchical Dirichlet Process Weibull Mixture Model          -- %%%
%%% ----------    for prognosis analysis (cancer patients)    ---------- %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear;clc;close all;

tic

seed = 3;
operation = 'ClusterSurvivalVariable'; % RawVariable; ClusterSurvivalVariable; All
% RawVariable: create new excel file for Chi-square test
% ClusterSurvivalVariable: cluster-specific data summary & survival figures

CensorFlag = 'Original'; % Fixed; Original (ratio)
if strcmp(CensorFlag,'Original')
    ParaInitialFlag = 'NoCensor'; 
elseif strcmp(CensorFlag,'Fixed')
    censor_rate = 0.3;
end

SurvivalFlag = 'Overall'; % Specific; Overall
AgeFlag = 'Continuous'; % Continuous; Discrete
Est_KM = 0; % show KM or not
EstBound_KM = 0; % confidence interval for KM estimation (in Survival)

folder = ['results/CaseStudy/seed_',num2str(seed),'/'];
L = 200; % for grid values in function estimates
% Unknown = 1; % 1: include Unknown; 0: exclude Unknown
year_extent = 5; % extented range (years)

lgd_cluster = 1;% legend or not for cluster-specific survival
lgd_KM = 1; % legend or not for KM

%% load data
%---- Selected raw data ----%

if strcmp(CensorFlag,'Fixed')
    filename = [folder,'data_censor_',CensorFlag,'_',num2str(censor_rate*10),'_seed_',num2str(seed),'.xlsx'];
    filename2 = [folder,'HDP_seed_',num2str(seed),'_censor_',CensorFlag,'_',num2str(censor_rate*10),'.mat'];

else
    filename = [folder,'data_censor_',CensorFlag,'_seed_',num2str(seed),'.xlsx'];
    filename2 = [folder,'HDP_seed_',num2str(seed),'_censor_',CensorFlag,'.mat'];
end
data_select_final = readtable(filename,'VariableNamingRule','preserve'); 
data_select_final_new = data_select_final; % included variables for Pearson's Chi-squared test
VarNames = data_select_final.Properties.VariableNames;
J_group = 2;
group_name = {'HPV Positive','HPV Negative'}; 

if strcmp(operation,'ClusterSurvivalVariable') || strcmp(operation,'All')
    
%---- Selected analyzed data ----%
load(filename2,'data','J_group','Num_J')
% 1: survival time (1-2)
% 2: censor indicator (0: not censored; 1: censored)
% 3: group (restaurant)
% 4: original survival time
% 5: dish index
% 6: component parameter alpha
% 7: component parameter beta
% 8: table index

%---- Estimated clustering label ----% LabelEst_Overall_HDP
if strcmp(SurvivalFlag,'Overall')
    filename_HDP = [folder,'LabelEst_Overall_HDP.xlsx'];
    filename_DP_All = [folder,'LabelEst_Overall_DP_All.xlsx'];
    filename_DP_Each = [folder,'LabelEst_Overall_DP_Each.xlsx'];
elseif strcmp(SurvivalFlag,'Specific')
    filename_HDP = [folder,'LabelEst_Specific_HDP.xlsx'];
    filename_DP_All = [folder,'LabelEst_Specific_DP_All.xlsx'];
    filename_DP_Each = [folder,'LabelEst_Specific_DP_Each.xlsx'];
end

%%% HDP-WMM (proposed) %%%
ClusterLabelEst = readtable(filename_HDP,'VariableNamingRule','preserve');
LabelEst_HDP = ClusterLabelEst(:,1);
LabelEst_HDP = table2array(LabelEst_HDP);

%%%  DP-WMM (all data)  %%%
ClusterLabelEst = readtable(filename_DP_All,'VariableNamingRule','preserve');
LabelEst_DP_All = ClusterLabelEst(:,1);
LabelEst_DP_All = table2array(LabelEst_DP_All);

%%% DP-WMM (each group) %%%
ClusterLabelEst = readtable(filename_DP_Each,'VariableNamingRule','preserve');
LabelEst_DP_Each = ClusterLabelEst(:,1);
LabelEst_DP_Each = table2array(LabelEst_DP_Each);

end


%% variable summary
if strcmp(operation,'RawVariable') || strcmp(operation,'All')
    
    disp('============================== Variable summary (All) ==============================')
    
    % group (HPV)
    ColId = ismember(VarNames,{'HPV recode (2010-2017)'});
    HPV = data_select_final(:,ColId);
    HPV = table2cell(HPV);
    tbl = tabulate(HPV);
    disp('===> Group summary: ')
    disp(tbl)
    group = HPV;
    
    
    % survival month
    disp('------------------------ survival month ------------------------ ')
    for j=1:J_group
        disp(['---> ',group_name{j},':'])
        idd = ismember(group,group_name{j});
        ColId = ismember(VarNames,{'Survival months'});
        var = data_select_final(idd,ColId);
        var = table2array(var);
        disp(['     mean & std :  ',num2str(mean(var)),', ',num2str(std(var))])
        disp(['     min & max :  ',num2str(min(var)),', ',num2str(max(var))])
    end
    
    % censor indicator
    disp('------------------------ censor indicator ------------------------ ')
    for j=1:J_group
        disp(['---> ',group_name{j},':'])
        idd = ismember(group,group_name{j});
        ColId = ismember(VarNames,{'Vital status recode (study cutoff used)'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        tbl = tabulate(var);
        disp(tbl)
    end
    
    % age
    disp('------------------------ age ------------------------ ')
    for j=1:J_group
        disp(['---> ',group_name{j},':'])
        idd = ismember(group,group_name{j});
        if strcmp(AgeFlag,'Continuous') % (Continuous)
           ColId = ismember(VarNames,{'Age'});
           var = data_select_final(idd,ColId);
           var = table2array(var);
           disp(['     mean & std :  ',num2str(mean(var)),', ',num2str(std(var))])
        elseif strcmp(AgeFlag,'Discrete') % (groups)
           ColId = ismember(VarNames,{'Age Standard for Survival (15-44,45-54,55-64,65-74,75+)'});
           var = data_select_final(idd,ColId);
           var = table2cell(var);
           var(ismember(var,{'00 years','35-39 years','40-44 years',...
               '45-49 years','50-54 years'})) = {'<= 54'};
           var(ismember(var,{'55-59 years','60-64 years'})) = {'55-64'};
           var(ismember(var,{'65-69 years','70-74 years'})) = {'65-74'};
           var(ismember(var,{'75-79 years','80-84 years','85+ years'})) = {'75+'};
           tbl = tabulate(var);
           disp(tbl)
           data_select_final_new(idd,ColId) = var;
        end
        
    end
    
    
    % sex
    disp('------------------------ sex ------------------------ ')
    for j=1:J_group
        disp(['---> ',group_name{j},':'])
        idd = ismember(group,group_name{j});
        ColId = ismember(VarNames,{'Sex'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        tbl = tabulate(var);
        disp(tbl)
    end
    
    
    % marital status
    disp('------------------------ marital status ------------------------ ')
    for j=1:J_group
        disp(['---> ',group_name{j},':'])
        idd = ismember(group,group_name{j});
        ColId = ismember(VarNames,{'Marital status at diagnosis'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        var(ismember(var,{'Married (including common law)'})) = {'Married'};
        var(ismember(var,{'Single (never married)','Divorced','Widowed',...
            'Separated','Unmarried or Domestic Partner'}))... % ,'Unknown'
            = {'Unmarried'}; % /Unknown
        tbl = tabulate(var);
        disp(tbl)
        data_select_final_new(idd,ColId) = var;
    end
    
    
    % race 
    disp('------------------------ race ------------------------')
    for j=1:J_group
        disp(['---> ',group_name{j},':'])
        idd = ismember(group,group_name{j});
        ColId = ismember(VarNames,{'Race recode (White, Black, Other)'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        var(~ismember(var,{'White','Black'})) = {'Other/Unknown'};
        tbl = tabulate(var);
        disp(tbl)
        data_select_final_new(idd,ColId) = var;
    end
 
    
    % cancer site
    disp('------------------------ cancer site ------------------------ ')
    for j=1:J_group
        disp(['---> ',group_name{j},':'])
        idd = ismember(group,group_name{j});
        if 1==1
            ColId = ismember(VarNames,{'Primary Site - labeled'});
            var = data_select_final(idd,ColId);
            var = table2cell(var);
            var(ismember(var,{'C09.9-Tonsil, NOS','C09.0-Tonsillar fossa',...
                'C10.9-Oropharynx, NOS','C09.1-Tonsillar pillar',...
                'C10.8-Overlapping lesion of oropharynx',...
                'C10.0-Vallecula','C10.3-Posterior wall of oropharynx',...
                'C09.8-Overlapping lesion of tonsil',...
                'C10.2-Lateral wall of oropharynx',...
                ''})) = {'Oropharynx'};
            var(~ismember(var,{'Oropharynx'})) = {'Non-Oropharynx'};
            tbl = tabulate(var);
            disp(tbl)
        else
            ColId = ismember(VarNames,{'Site recode ICD-O-3/WHO 2008'});
            var = data_select_final(idd,ColId);
            var = table2cell(var);
            var(~ismember(var,{'Oropharynx'})) = {'Non-Oropharynx'};
            tbl = tabulate(var);
            disp(tbl)
        end
        data_select_final_new(idd,ColId) = var;
    end
    
    % summary stage
    disp('------------------------ summary stage ------------------------ ')
    for j=1:J_group
        disp(['---> ',group_name{j},':'])
        idd = ismember(group,group_name{j});
        ColId = ismember(VarNames,{'Combined Summary Stage (2004+)'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        tbl = tabulate(var);
        disp(tbl)
    end
        
    % Surgery 
    disp('------------------------ Surgery ------------------------ ')
    for j=1:J_group
        disp(['---> ',group_name{j},':'])
        idd = ismember(group,group_name{j});
        ColId = ismember(VarNames,{'Reason no cancer-directed surgery'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        var(~ismember(var,{'Surgery performed'})) = {'No/Unknown'};
        var(ismember(var,{'Surgery performed'})) = {'Yes'};
        tbl = tabulate(var);
        disp(tbl)
        data_select_final_new(idd,ColId) = var;
    end
    
    % Radiation 
    disp('------------------------ Radiation ------------------------ ')
    for j=1:J_group
        disp(['---> ',group_name{j},':'])
        idd = ismember(group,group_name{j});
        ColId = ismember(VarNames,{'Radiation recode'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        var(~ismember(var,{'Beam radiation',...
            'Radiation, NOS  method or source not specified',...
            'Combination of beam with implants or isotopes',...
            'Radioactive implants (includes brachytherapy) (1988+)',...
            'Radioisotopes (1988+)'})) = {'No/Unknown'};
        var(ismember(var,{'Beam radiation',...
            'Radiation, NOS  method or source not specified',...
            'Combination of beam with implants or isotopes',...
            'Radioactive implants (includes brachytherapy) (1988+)',...
            'Radioisotopes (1988+)'})) = {'Yes'};
        tbl = tabulate(var);
        disp(tbl)
        data_select_final_new(idd,ColId) = var;
    end

    % Chemotherapy 
    disp('------------------------ Chemotherapy ------------------------ ')
    for j=1:J_group
        disp(['---> ',group_name{j},':'])
        idd = ismember(group,group_name{j});
        ColId = ismember(VarNames,{'Chemotherapy recode (yes, no/unk)'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        tbl = tabulate(var);
        disp(tbl)
        data_select_final_new(idd,ColId) = var;
    end
    
    
    % Treatment combination
    ColId = ismember(VarNames,{'Reason no cancer-directed surgery'});
    Surgery = data_select_final(:,ColId);
    Surgery = table2cell(Surgery);
    Surgery(ismember(Surgery,{'Surgery performed'})) = {'Yes'};
    ColId = ismember(VarNames,{'Radiation recode'});
    Radiation = data_select_final(:,ColId);
    Radiation = table2cell(Radiation);
    Radiation(ismember(Radiation,{'Beam radiation',...
            'Radiation, NOS  method or source not specified',...
            'Combination of beam with implants or isotopes',...
            'Radioactive implants (includes brachytherapy) (1988+)',...
            'Radioisotopes (1988+)'})) = {'Yes'};
    ColId = ismember(VarNames,{'Chemotherapy recode (yes, no/unk)'});
    Chemotherapy = data_select_final(:,ColId);
    Chemotherapy = table2cell(Chemotherapy);
    Treatment = Chemotherapy;
    idd = (ismember(Surgery,{'Yes'})) .* (~ismember(Radiation,{'Yes'})) .* (~ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Surgery only'};
    idd = (~ismember(Surgery,{'Yes'})) .* (ismember(Radiation,{'Yes'})) .* (~ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Radiation only'};
    idd = (~ismember(Surgery,{'Yes'})) .* (~ismember(Radiation,{'Yes'})) .* (ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Chemotherapy only'};
    idd = (ismember(Surgery,{'Yes'})) .* (ismember(Radiation,{'Yes'})) .* (~ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Surgery + Radiation'};
    idd = (ismember(Surgery,{'Yes'})) .* (~ismember(Radiation,{'Yes'})) .* (ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Surgery + Chemotherapy'};
    idd = (~ismember(Surgery,{'Yes'})) .* (ismember(Radiation,{'Yes'})) .* (ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Radiation + Chemotherapy'};
    idd = (ismember(Surgery,{'Yes'})) .* (ismember(Radiation,{'Yes'})) .* (ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Surgery + Radiation + Chemotherapy'};
    idd = (~ismember(Surgery,{'Yes'})) .* (~ismember(Radiation,{'Yes'})) .* (~ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'No treatment'};
    
    Treatment = cell2table(Treatment);
    data_select_final_new = [data_select_final_new,Treatment];
        
    % Treatment combination 
    disp('------------------------ Treatment combination  ------------------------ ')
    for j=1:J_group
        disp(['---> ',group_name{j},':'])
        idd = ismember(group,group_name{j});
        var = Treatment(idd,1);
        var = table2cell(var);
        tbl = tabulate(var);
        disp(tbl)
    end
    
    %-------------------------- Save new excel file for Chi-square test --------------------------% 
    if strcmp(CensorFlag,'Fixed')
        filename = [folder,'data_censor_',CensorFlag,'_',num2str(censor_rate*10),'_seed_',num2str(seed),'_Chi.xlsx'];

    else
        filename = [folder,'data_censor_',CensorFlag,'_seed_',num2str(seed),'_Chi.xlsx'];
    end
    writetable(data_select_final_new,filename)


       
end


%% Clustering Survival (HDP-WMM)
if strcmp(operation,'ClusterSurvivalVariable') || strcmp(operation,'All')
    % for cluster-specific survival figures
    figure(1);
    tiledlayout(1,3,'TileSpacing','compact')
    
    % survival data
    sur_original = data(:,4);
    X_max = max(sur_original);
    sur = data(:,1);
    group = data(:,3);
    % censoring indicator (1) ----> Overall Survival
    censor = data(:,2);    
 
    disp('==================================== Overall Survival (HDP-WMM) ===================================')
    % Weibull fit for each cluster (Overall Survival)
    ClusterNum = length(unique(LabelEst_HDP));
    ParaEst = nan(ClusterNum,2);
    disp(['--- Total number of estimated clusters: ', num2str(ClusterNum)])
    disp('--------------- Cluster summary (All) --------------------')
    tbl = tabulate(LabelEst_HDP);
    tbl = array2table(tbl,'VariableNames',{'Cluster','Count','Frequency'});
    disp(tbl)
    for j=1:J_group
        disp(['--------------- Cluster summary (Group ',num2str(j),') --------------------'])
        idd = (group==j);
        Label_temp = LabelEst_HDP(idd);
        tbl = tabulate(Label_temp);
        tbl = array2table(tbl,'VariableNames',{'Cluster','Count','Frequency'});
        disp(tbl)
    end
    
    %------------ order clustering labels based survival prob ------------%
    SurProb_1year = nan(1,ClusterNum); % 1 year (12 month); 
    for i=1:ClusterNum
        idd = (LabelEst_HDP==i);
        dat = sur(idd);
        censor_temp = censor(idd);
        [parmHat,~] = wblfit(dat,0.05,censor_temp);
        ParaEst(i,:) = parmHat; % parmHat: a (scale) and b (shape) 
        temp = wblcdf(1, parmHat(1),parmHat(2));
        t = 1 + (1*12)/X_max; % 1-year
        SurProb_1year(i) = 100 * (1 - (wblcdf(t,parmHat(1),parmHat(2))-temp)/(1-temp));
    end
    [~,I] = sort(SurProb_1year,'descend');
    LabelEst_HDP_sort = nan(size(LabelEst_HDP));
    for i=1:ClusterNum
        LabelEst_HDP_sort(LabelEst_HDP==I(i)) = i;
    end
    LabelEst_HDP = LabelEst_HDP_sort;
    
    %------------ figures and data summary based ordered labels ------------%
    SurProb = nan(ClusterNum,3); % 1 year (12 month); 5 year (60 month); 10 year (120 month)
    nexttile
    color_order = get(gca,'colororder');
%     color_order = rand(30,3);
    if Est_KM == 1
        if EstBound_KM == 1
            legend_name = cell(1,3*ClusterNum);
        else
            legend_name = cell(1,2*ClusterNum);
        end  
    else
        legend_name = cell(1,ClusterNum);
    end
    linestyle = {'-','-.',':','--'}; 
    xx = linspace(0,X_max+12*year_extent,L); % grid values
    yy = 1 + xx/X_max; 
    for i=1:ClusterNum
        idd = (LabelEst_HDP==i);
        dat = sur(idd);
        censor_temp = censor(idd);
        [parmHat,~] = wblfit(dat,0.05,censor_temp);
        ParaEst(i,:) = parmHat; % parmHat: a (scale) and b (shape) 
        temp = wblcdf(1, parmHat(1),parmHat(2));
        t = 1 + (1*12)/X_max; % 1-year
        SurProb(i,1) = 100 * (1 - (wblcdf(t,parmHat(1),parmHat(2))-temp)/(1-temp));
        t = 1 + (3*12)/X_max; % 3-year
        SurProb(i,2) = 100 * (1 - (wblcdf(t,parmHat(1),parmHat(2))-temp)/(1-temp));
        t = 1 + (5*12)/X_max; % 5-year
        SurProb(i,3) = 100 * (1 - (wblcdf(t,parmHat(1),parmHat(2))-temp)/(1-temp));
        if Est_KM == 1
            if EstBound_KM == 1
                legend_name{(i-1)*3+1} = ['Cluster ',num2str(i)]; %,'(HDP-WMM)'
                legend_name{(i-1)*3+2} = ['Cluster ',num2str(i),'(KM)'];
                legend_name{(i-1)*3+3} = ['Cluster ',num2str(i),'(KM, 95% CI)'];
            else
                legend_name{(i-1)*2+1} = ['Cluster ',num2str(i)]; %,'(HDP-WMM)'
                legend_name{(i-1)*2+2} = ['Cluster ',num2str(i),'(KM)'];
            end
        else
            legend_name{i} = ['Cluster ',num2str(i)];
        end
        cc = 1 - (wblcdf(yy,parmHat(1),parmHat(2))-temp)/(1-temp);
        plot(xx,cc,'LineWidth', 2,'Color',color_order(i,:),'linestyle',linestyle{1}); hold on % color_order: +J_group+1
        
        if Est_KM == 1
            sur_original_temp = sur_original(idd);
            [ss_e,x,flo,fup] = ecdf(sur_original_temp,'Censoring',censor_temp,'function','survivor');
            stairs(x, ss_e,'LineWidth', 2,'Color',color_order(i,:),'linestyle',linestyle{2});
            if EstBound_KM == 1
                stairs(x, flo,'LineWidth', 1.5,'Color',color_order(i,:),'linestyle',linestyle{3})
                stairs(x, fup,'LineWidth', 1.5,'Color',color_order(i,:),'linestyle',linestyle{3},'HandleVisibility','off')
            end
        end
    end
    
    set(gca,'FontSize',16,'FontName', 'times')
    xlim([0 max(xx)]); ylim([0 1]);
    title('HDP-WMM','FontSize',16,'FontName', 'times','FontWeight','Normal')
    ylabel('Survival probability','FontSize',16,'FontName', 'times')
    xlabel('Time (month)','FontSize',16,'FontName', 'times')
    
%     if lgd_cluster==1
%         legend(legend_name,'FontSize',14,'FontName', 'times','Location','best') % 
%         legend('boxoff')
%     end

    
    disp('--------------- Survival function parameters --------------------')
    ParaEst = array2table(ParaEst,'VariableNames',{'Scale','Shape'});
    disp(ParaEst)
    disp('--------------- Survival Probabilities --------------------')
    SurProb = array2table(SurProb,'VariableNames',{'1-year survival',...
        '3-year survival','5-year survival'});
    disp(SurProb)
    
    
    
    disp('==================================== Overall Survival (DP-WMM-All) ===================================')
    % Weibull fit for each cluster (Overall Survival)
    ClusterNum = length(unique(LabelEst_DP_All));
    ParaEst = nan(ClusterNum,2);
    disp(['--- Total number of estimated clusters: ', num2str(ClusterNum)])
    disp('--------------- Cluster summary (All) --------------------')
    tbl = tabulate(LabelEst_DP_All);
    tbl = array2table(tbl,'VariableNames',{'Cluster','Count','Frequency'});
    disp(tbl)
    for j=1:J_group
        disp(['--------------- Cluster summary (Group ',num2str(j),') --------------------'])
        idd = (group==j);
        Label_temp = LabelEst_DP_All(idd);
        tbl = tabulate(Label_temp);
        tbl = array2table(tbl,'VariableNames',{'Cluster','Count','Frequency'});
        disp(tbl)
    end
    
    %------------ order clustering labels based survival prob ------------%
    SurProb_1year = nan(1,ClusterNum); % 1 year (12 month); 
    for i=1:ClusterNum
        idd = (LabelEst_DP_All==i);
        dat = sur(idd);
        censor_temp = censor(idd);
        [parmHat,~] = wblfit(dat,0.05,censor_temp);
        ParaEst(i,:) = parmHat; % parmHat: a (scale) and b (shape) 
        temp = wblcdf(1, parmHat(1),parmHat(2));
        t = 1 + (1*12)/X_max; % 1-year
        SurProb_1year(i) = 100 * (1 - (wblcdf(t,parmHat(1),parmHat(2))-temp)/(1-temp));
    end
    [~,I] = sort(SurProb_1year,'descend');
    LabelEst_DP_All_sort = nan(size(LabelEst_DP_All));
    for i=1:ClusterNum
        LabelEst_DP_All_sort(LabelEst_DP_All==I(i)) = i;
    end
    LabelEst_DP_All = LabelEst_DP_All_sort;
    
    %------------ figures and data summary based ordered labels ------------%
    SurProb = nan(ClusterNum,3); % 1 year (12 month); 5 year (60 month); 10 year (120 month)
    nexttile
%     color_order = get(gca,'colororder');
    if Est_KM == 1
        if EstBound_KM == 1
            legend_name = cell(1,3*ClusterNum);
        else
            legend_name = cell(1,2*ClusterNum);
        end  
    else
        legend_name = cell(1,ClusterNum);
    end  
    linestyle = {'-','-.',':','--'}; 
    xx = linspace(0,X_max+12*year_extent,L); % grid values
    yy = 1 + xx/X_max; 
    for i=1:ClusterNum
        idd = (LabelEst_DP_All==i);
        dat = sur(idd);
        censor_temp = censor(idd);
        [parmHat,~] = wblfit(dat,0.05,censor_temp);
        ParaEst(i,:) = parmHat; % parmHat: a (scale) and b (shape) 
        temp = wblcdf(1, parmHat(1),parmHat(2));
        t = 1 + (1*12)/X_max; % 1-year
        SurProb(i,1) = 100 * (1 - (wblcdf(t,parmHat(1),parmHat(2))-temp)/(1-temp));
        t = 1 + (3*12)/X_max; % 3-year
        SurProb(i,2) = 100 * (1 - (wblcdf(t,parmHat(1),parmHat(2))-temp)/(1-temp));
        t = 1 + (5*12)/X_max; % 5-year
        SurProb(i,3) = 100 * (1 - (wblcdf(t,parmHat(1),parmHat(2))-temp)/(1-temp));
        if Est_KM == 1
            if EstBound_KM == 1
                legend_name{(i-1)*3+1} = ['Cluster ',num2str(i)];% ,'(DP-WMM-All)'
                legend_name{(i-1)*3+2} = ['Cluster ',num2str(i),'(KM)'];
                legend_name{(i-1)*3+3} = ['Cluster ',num2str(i),'(KM, 95% CI)'];
            else
                legend_name{(i-1)*2+1} = ['Cluster ',num2str(i)]; %,'(DP-WMM-All)'
                legend_name{(i-1)*2+2} = ['Cluster ',num2str(i),'(KM)'];
            end
        else
            legend_name{i} = ['Cluster ',num2str(i)];
        end
        cc = 1 - (wblcdf(yy,parmHat(1),parmHat(2))-temp)/(1-temp);
        plot(xx,cc,'LineWidth', 2,'Color',color_order(i,:),'linestyle',linestyle{1}); hold on
        
        if Est_KM == 1
            sur_original_temp = sur_original(idd);
            [ss_e,x,flo,fup] = ecdf(sur_original_temp,'Censoring',censor_temp,'function','survivor');
            stairs(x, ss_e,'LineWidth', 2,'Color',color_order(i,:),'linestyle',linestyle{2});
            if EstBound_KM == 1
                stairs(x, flo,'LineWidth', 1.5,'Color',color_order(i,:),'linestyle',linestyle{3})
                stairs(x, fup,'LineWidth', 1.5,'Color',color_order(i,:),'linestyle',linestyle{3},'HandleVisibility','off')
            end
        end
    end
    set(gca,'FontSize',16,'FontName', 'times')
    xlim([0 max(xx)]); ylim([0 1]);
    title('DP-WMM-All','FontSize',16,'FontName', 'times','FontWeight','Normal')
%     xlabel('Time (month)')
   
%     if lgd_cluster==1
%         legend(legend_name,'FontSize',14,'FontName', 'times','Location','best') % 
%         legend('boxoff')
%     end

    
    disp('--------------- Survival function parameters --------------------')
    ParaEst = array2table(ParaEst,'VariableNames',{'Scale','Shape'});
    disp(ParaEst)
    disp('--------------- Survival Probabilities --------------------')
    SurProb = array2table(SurProb,'VariableNames',{'1-year survival',...
        '3-year survival','5-year survival'});
    disp(SurProb)
    
    
    
    
    disp('==================================== Overall Survival (DP-WMM-Each) ===================================')
    % Weibull fit for each cluster (Overall Survival)
    ClusterNum = length(unique(LabelEst_DP_Each));
    ParaEst = nan(ClusterNum,2);
    disp(['--- Total number of estimated clusters: ', num2str(ClusterNum)])
    disp('--------------- Cluster summary (All) --------------------')
    tbl = tabulate(LabelEst_DP_Each);
    tbl = array2table(tbl,'VariableNames',{'Cluster','Count','Frequency'});
    disp(tbl)
    for j=1:J_group
        disp(['--------------- Cluster summary (Group ',num2str(j),') --------------------'])
        idd = (group==j);
        Label_temp = LabelEst_DP_Each(idd);
        tbl = tabulate(Label_temp);
        tbl = array2table(tbl,'VariableNames',{'Cluster','Count','Frequency'});
        disp(tbl)
    end
    
    %------------ order clustering labels based survival prob ------------%
    SurProb_1year = nan(1,ClusterNum); % 1 year (12 month); 
    for i=1:ClusterNum
        idd = (LabelEst_DP_Each==i);
        dat = sur(idd);
        censor_temp = censor(idd);
        [parmHat,~] = wblfit(dat,0.05,censor_temp);
        ParaEst(i,:) = parmHat; % parmHat: a (scale) and b (shape) 
        temp = wblcdf(1, parmHat(1),parmHat(2));
        t = 1 + (1*12)/X_max; % 1-year
        SurProb_1year(i) = 100 * (1 - (wblcdf(t,parmHat(1),parmHat(2))-temp)/(1-temp));
    end
    [~,I] = sort(SurProb_1year,'descend');
    LabelEst_DP_Each_sort = nan(size(LabelEst_DP_Each));
    for i=1:ClusterNum
        LabelEst_DP_Each_sort(LabelEst_DP_Each==I(i)) = i;
    end
    LabelEst_DP_Each = LabelEst_DP_Each_sort;
    
    %------------ figures and data summary based ordered labels ------------%
    SurProb = nan(ClusterNum,3); % 1 year (12 month); 5 year (60 month); 10 year (120 month)
    nexttile
%     color_order = get(gca,'colororder');
    if Est_KM == 1
        if EstBound_KM == 1
            legend_name = cell(1,3*ClusterNum);
        else
            legend_name = cell(1,2*ClusterNum);
        end  
    else
        legend_name = cell(1,ClusterNum);
    end
    linestyle = {'-','-.',':','--'}; 
    xx = linspace(0,X_max+12*year_extent,L); % grid values
    yy = 1 + xx/X_max; 
    for i=1:ClusterNum
        idd = (LabelEst_DP_Each==i);
        dat = sur(idd);
        censor_temp = censor(idd);
        [parmHat,~] = wblfit(dat,0.05,censor_temp);
        ParaEst(i,:) = parmHat; % parmHat: a (scale) and b (shape) 
        temp = wblcdf(1, parmHat(1),parmHat(2));
        t = 1 + (1*12)/X_max; % 1-year
        SurProb(i,1) = 100 * (1 - (wblcdf(t,parmHat(1),parmHat(2))-temp)/(1-temp));
        t = 1 + (3*12)/X_max; % 3-year
        SurProb(i,2) = 100 * (1 - (wblcdf(t,parmHat(1),parmHat(2))-temp)/(1-temp));
        t = 1 + (5*12)/X_max; % 5-year
        SurProb(i,3) = 100 * (1 - (wblcdf(t,parmHat(1),parmHat(2))-temp)/(1-temp));
        if Est_KM == 1
            if EstBound_KM == 1
                legend_name{(i-1)*3+1} = ['Cluster ',num2str(i)]; %,'(DP-WMM-Each)'
                legend_name{(i-1)*3+2} = ['Cluster ',num2str(i),'(KM)'];
                legend_name{(i-1)*3+3} = ['Cluster ',num2str(i),'(KM, 95% CI)'];
            else
                legend_name{(i-1)*2+1} = ['Cluster ',num2str(i)]; %,'(DP-WMM-Each)'
                legend_name{(i-1)*2+2} = ['Cluster ',num2str(i),'(KM)'];
            end
        else
            legend_name{i} = ['Cluster ',num2str(i)];
        end
        cc = 1 - (wblcdf(yy,parmHat(1),parmHat(2))-temp)/(1-temp);
        plot(xx,cc,'LineWidth', 2,'Color',color_order(i,:),'linestyle',linestyle{1}); hold on
        
        if Est_KM == 1
            sur_original_temp = sur_original(idd);
            [ss_e,x,flo,fup] = ecdf(sur_original_temp,'Censoring',censor_temp,'function','survivor');
            stairs(x, ss_e,'LineWidth', 2,'Color',color_order(i,:),'linestyle',linestyle{2});
            if EstBound_KM == 1
                stairs(x, flo,'LineWidth', 1.5,'Color',color_order(i,:),'linestyle',linestyle{3})
                stairs(x, fup,'LineWidth', 1.5,'Color',color_order(i,:),'linestyle',linestyle{3},'HandleVisibility','off')
            end
        end
    end
    set(gca,'FontSize',16,'FontName', 'times')
    xlim([0 max(xx)]); ylim([0 1]);
    title('DP-WMM-Each','FontSize',16,'FontName', 'times','FontWeight','Normal')
%     xlabel('Time (month)')
%     
%     if lgd_cluster==1
%         legend(legend_name,'FontSize',14,'FontName', 'times','Location','best') % 
%         legend('boxoff')
%     end
    lg = legend(legend_name,'Orientation','Vertical','FontSize',14,'FontName', 'times'); 
    legend('boxoff')
    lg.Layout.Tile = 'East'; 
    

    disp('--------------- Survival function parameters --------------------')
    ParaEst = array2table(ParaEst,'VariableNames',{'Scale','Shape'});
    disp(ParaEst)
    disp('--------------- Survival Probabilities --------------------')
    SurProb = array2table(SurProb,'VariableNames',{'1-year survival',...
        '3-year survival','5-year survival'});
    disp(SurProb)
    
   
end

%% Variable summary for each cluster
if strcmp(operation,'ClusterSurvivalVariable') || strcmp(operation,'All')

    disp('=========================== Variable summary (each cluster, HDP-WMM) ===========================')
    LabelEst = LabelEst_HDP;
    ClusterNum = length(unique(LabelEst));
   
    
    % survival month
    disp('------------------------ survival month ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Survival months'});
        var = data_select_final(idd,ColId);
        var = table2array(var);
        disp(['     mean & std :  ',num2str(mean(var)),', ',num2str(std(var))])
        disp(['     min & max :  ',num2str(min(var)),', ',num2str(max(var))])
    end
    
    % censor indicator
    disp('------------------------ censor indicator ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Vital status recode (study cutoff used)'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        tbl = tabulate(var);
        disp(tbl)
    end
    
    % HPV status
    disp('------------------------ HPV status ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'HPV recode (2010-2017)'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        tbl = tabulate(var);
        disp(tbl)
    end
    
    
    % age
    disp('------------------------ age ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        if strcmp(AgeFlag,'Continuous') % (Continuous)
           ColId = ismember(VarNames,{'Age'});
           var = data_select_final(idd,ColId);
           var = table2array(var);
           disp(['     mean & std :  ',num2str(mean(var)),', ',num2str(std(var))])
        elseif strcmp(AgeFlag,'Discrete') % (groups)
           ColId = ismember(VarNames,{'Age Standard for Survival (15-44,45-54,55-64,65-74,75+)'});
           var = data_select_final(idd,ColId);
           var = table2cell(var);
           var(ismember(var,{'00 years','35-39 years','40-44 years',...
               '45-49 years','50-54 years'})) = {'<= 54'};
           var(ismember(var,{'55-59 years','60-64 years'})) = {'55-64'};
           var(ismember(var,{'65-69 years','70-74 years'})) = {'65-74'};
           var(ismember(var,{'75-79 years','80-84 years','85+ years'})) = {'75+'};
           tbl = tabulate(var);
           disp(tbl)
        end
    end
    
    % sex
    disp('------------------------ sex ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Sex'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        tbl = tabulate(var);
        disp(tbl)
    end
    
    % marital status
    disp('------------------------ marital status ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Marital status at diagnosis'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        var(ismember(var,{'Married (including common law)'})) = {'Married'};
        var(ismember(var,{'Single (never married)','Divorced','Widowed',...
            'Separated','Unmarried or Domestic Partner'}))... % ,'Unknown'
            = {'Unmarried'}; % /Unknown
        tbl = tabulate(var);
        disp(tbl)
    end
    
    % race 
    disp('------------------------ race ------------------------')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Race recode (White, Black, Other)'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        var(~ismember(var,{'White','Black'})) = {'Other/Unknown'};
        tbl = tabulate(var);
        disp(tbl)
    end
    
    
    % cancer site
    disp('------------------------ cancer site ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Primary Site - labeled'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        var(ismember(var,{'C09.9-Tonsil, NOS','C09.0-Tonsillar fossa',...
                'C10.9-Oropharynx, NOS','C09.1-Tonsillar pillar',...
                'C10.8-Overlapping lesion of oropharynx',...
                'C10.0-Vallecula','C10.3-Posterior wall of oropharynx',...
                'C09.8-Overlapping lesion of tonsil',...
                'C10.2-Lateral wall of oropharynx',...
                ''})) = {'Oropharynx'};
        var(~ismember(var,{'Oropharynx'})) = {'Non-Oropharynx'};
        tbl = tabulate(var);
        disp(tbl)
    end
    
    % summary stage
    disp('------------------------ summary stage ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Combined Summary Stage (2004+)'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        tbl = tabulate(var);
        disp(tbl)
    end
    
%     % AJCC stage, 7th 
%     disp('------------------------ AJCC stage, 7th ------------------------ ')
%     for j=1:ClusterNum
%         disp(['---> Cluster ',num2str(j)])
%         idd = (LabelEst==j);
%         var = data_select_final(idd,15);
%         var = table2cell(var);
%         var(ismember(var,{'IVA','IVB','IVC','IVNOS'})) = {'IV'};
%         var(ismember(var,{'Blank(s)','NA','UNK Stage'})) = {'Unknown'};
%         tbl = tabulate(var);
%         disp(tbl)
%     end
    

    % Surgery 
    disp('------------------------ Surgery ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Reason no cancer-directed surgery'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        var(~ismember(var,{'Surgery performed'})) = {'No/Unknown'};
        var(ismember(var,{'Surgery performed'})) = {'Yes'};
        tbl = tabulate(var);
        disp(tbl)
    end
    
    
    % Radiation 
    disp('------------------------ Radiation ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Radiation recode'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        var(~ismember(var,{'Beam radiation',...
            'Radiation, NOS  method or source not specified',...
            'Combination of beam with implants or isotopes',...
            'Radioactive implants (includes brachytherapy) (1988+)',...
            'Radioisotopes (1988+)'})) = {'No/Unknown'};
        var(ismember(var,{'Beam radiation',...
            'Radiation, NOS  method or source not specified',...
            'Combination of beam with implants or isotopes',...
            'Radioactive implants (includes brachytherapy) (1988+)',...
            'Radioisotopes (1988+)'})) = {'Yes'};
        tbl = tabulate(var);
        disp(tbl)
    end
    
    
    % Chemotherapy 
    disp('------------------------ Chemotherapy ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Chemotherapy recode (yes, no/unk)'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        tbl = tabulate(var);
        disp(tbl)
    end
    
    % Treatment combination
    ColId = ismember(VarNames,{'Reason no cancer-directed surgery'});
    Surgery = data_select_final(:,ColId);
    Surgery = table2cell(Surgery);
    Surgery(ismember(Surgery,{'Surgery performed'})) = {'Yes'};
    ColId = ismember(VarNames,{'Radiation recode'});
    Radiation = data_select_final(:,ColId);
    Radiation = table2cell(Radiation);
    Radiation(ismember(Radiation,{'Beam radiation',...
            'Radiation, NOS  method or source not specified',...
            'Combination of beam with implants or isotopes',...
            'Radioactive implants (includes brachytherapy) (1988+)',...
            'Radioisotopes (1988+)'})) = {'Yes'};
    ColId = ismember(VarNames,{'Chemotherapy recode (yes, no/unk)'});
    Chemotherapy = data_select_final(:,ColId);
    Chemotherapy = table2cell(Chemotherapy);
    Treatment = Chemotherapy;
    idd = (ismember(Surgery,{'Yes'})) .* (~ismember(Radiation,{'Yes'})) .* (~ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Surgery only'};
    idd = (~ismember(Surgery,{'Yes'})) .* (ismember(Radiation,{'Yes'})) .* (~ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Radiation only'};
    idd = (~ismember(Surgery,{'Yes'})) .* (~ismember(Radiation,{'Yes'})) .* (ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Chemotherapy only'};
    idd = (ismember(Surgery,{'Yes'})) .* (ismember(Radiation,{'Yes'})) .* (~ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Surgery + Radiation'};
    idd = (ismember(Surgery,{'Yes'})) .* (~ismember(Radiation,{'Yes'})) .* (ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Surgery + Chemotherapy'};
    idd = (~ismember(Surgery,{'Yes'})) .* (ismember(Radiation,{'Yes'})) .* (ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Radiation + Chemotherapy'};
    idd = (ismember(Surgery,{'Yes'})) .* (ismember(Radiation,{'Yes'})) .* (ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'All'};
    idd = (~ismember(Surgery,{'Yes'})) .* (~ismember(Radiation,{'Yes'})) .* (~ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'No treatment'};
    
    % Treatment combination 
    disp('------------------------ Treatment combination  ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        var = Treatment(idd);
        tbl = tabulate(var);
        disp(tbl)
    end
    
    
    
    disp('=========================== Variable summary (each cluster, DP-All) ===========================')
    LabelEst = LabelEst_DP_All;
    ClusterNum = length(unique(LabelEst));
   
    
    % survival month
    disp('------------------------ survival month ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Survival months'});
        var = data_select_final(idd,ColId);
        var = table2array(var);
        disp(['     mean & std :  ',num2str(mean(var)),', ',num2str(std(var))])
        disp(['     min & max :  ',num2str(min(var)),', ',num2str(max(var))])
    end
    
    % censor indicator
    disp('------------------------ censor indicator ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Vital status recode (study cutoff used)'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        tbl = tabulate(var);
        disp(tbl)
    end
    
    % HPV status
    disp('------------------------ HPV status ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'HPV recode (2010-2017)'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        tbl = tabulate(var);
        disp(tbl)
    end
    
    
    % age
    disp('------------------------ age ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        if strcmp(AgeFlag,'Continuous') % (Continuous)
           ColId = ismember(VarNames,{'Age'});
           var = data_select_final(idd,ColId);
           var = table2array(var);
           disp(['     mean & std :  ',num2str(mean(var)),', ',num2str(std(var))])
        elseif strcmp(AgeFlag,'Discrete') % (groups)
           ColId = ismember(VarNames,{'Age Standard for Survival (15-44,45-54,55-64,65-74,75+)'});
           var = data_select_final(idd,ColId);
           var = table2cell(var);
           var(ismember(var,{'00 years','35-39 years','40-44 years',...
               '45-49 years','50-54 years'})) = {'<= 54'};
           var(ismember(var,{'55-59 years','60-64 years'})) = {'55-64'};
           var(ismember(var,{'65-69 years','70-74 years'})) = {'65-74'};
           var(ismember(var,{'75-79 years','80-84 years','85+ years'})) = {'75+'};
           tbl = tabulate(var);
           disp(tbl)
        end
    end
    
    % sex
    disp('------------------------ sex ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Sex'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        tbl = tabulate(var);
        disp(tbl)
    end
    
    % marital status
    disp('------------------------ marital status ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Marital status at diagnosis'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        var(ismember(var,{'Married (including common law)'})) = {'Married'};
        var(ismember(var,{'Single (never married)','Divorced','Widowed',...
            'Separated','Unmarried or Domestic Partner'}))... % ,'Unknown'
            = {'Unmarried'}; % /Unknown
        tbl = tabulate(var);
        disp(tbl)
    end
    
    % race 
    disp('------------------------ race ------------------------')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Race recode (White, Black, Other)'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        var(~ismember(var,{'White','Black'})) = {'Other/Unknown'};
        tbl = tabulate(var);
        disp(tbl)
    end
    
    
    % cancer site
    disp('------------------------ cancer site ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Primary Site - labeled'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        var(ismember(var,{'C09.9-Tonsil, NOS','C09.0-Tonsillar fossa',...
                'C10.9-Oropharynx, NOS','C09.1-Tonsillar pillar',...
                'C10.8-Overlapping lesion of oropharynx',...
                'C10.0-Vallecula','C10.3-Posterior wall of oropharynx',...
                'C09.8-Overlapping lesion of tonsil',...
                'C10.2-Lateral wall of oropharynx',...
                ''})) = {'Oropharynx'};
        var(~ismember(var,{'Oropharynx'})) = {'Non-Oropharynx'};
        tbl = tabulate(var);
        disp(tbl)
    end
    
    % summary stage
    disp('------------------------ summary stage ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Combined Summary Stage (2004+)'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        tbl = tabulate(var);
        disp(tbl)
    end
    

    % Surgery 
    disp('------------------------ Surgery ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Reason no cancer-directed surgery'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        var(~ismember(var,{'Surgery performed'})) = {'No/Unknown'};
        var(ismember(var,{'Surgery performed'})) = {'Yes'};
        tbl = tabulate(var);
        disp(tbl)
    end
    
    
    % Radiation 
    disp('------------------------ Radiation ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Radiation recode'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        var(~ismember(var,{'Beam radiation',...
            'Radiation, NOS  method or source not specified',...
            'Combination of beam with implants or isotopes',...
            'Radioactive implants (includes brachytherapy) (1988+)',...
            'Radioisotopes (1988+)'})) = {'No/Unknown'};
        var(ismember(var,{'Beam radiation',...
            'Radiation, NOS  method or source not specified',...
            'Combination of beam with implants or isotopes',...
            'Radioactive implants (includes brachytherapy) (1988+)',...
            'Radioisotopes (1988+)'})) = {'Yes'};
        tbl = tabulate(var);
        disp(tbl)
    end
    
    
    % Chemotherapy 
    disp('------------------------ Chemotherapy ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Chemotherapy recode (yes, no/unk)'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        tbl = tabulate(var);
        disp(tbl)
    end
    
    % Treatment combination
    ColId = ismember(VarNames,{'Reason no cancer-directed surgery'});
    Surgery = data_select_final(:,ColId);
    Surgery = table2cell(Surgery);
    Surgery(ismember(Surgery,{'Surgery performed'})) = {'Yes'};
    ColId = ismember(VarNames,{'Radiation recode'});
    Radiation = data_select_final(:,ColId);
    Radiation = table2cell(Radiation);
    Radiation(ismember(Radiation,{'Beam radiation',...
            'Radiation, NOS  method or source not specified',...
            'Combination of beam with implants or isotopes',...
            'Radioactive implants (includes brachytherapy) (1988+)',...
            'Radioisotopes (1988+)'})) = {'Yes'};
    ColId = ismember(VarNames,{'Chemotherapy recode (yes, no/unk)'});
    Chemotherapy = data_select_final(:,ColId);
    Chemotherapy = table2cell(Chemotherapy);
    Treatment = Chemotherapy;
    idd = (ismember(Surgery,{'Yes'})) .* (~ismember(Radiation,{'Yes'})) .* (~ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Surgery only'};
    idd = (~ismember(Surgery,{'Yes'})) .* (ismember(Radiation,{'Yes'})) .* (~ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Radiation only'};
    idd = (~ismember(Surgery,{'Yes'})) .* (~ismember(Radiation,{'Yes'})) .* (ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Chemotherapy only'};
    idd = (ismember(Surgery,{'Yes'})) .* (ismember(Radiation,{'Yes'})) .* (~ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Surgery + Radiation'};
    idd = (ismember(Surgery,{'Yes'})) .* (~ismember(Radiation,{'Yes'})) .* (ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Surgery + Chemotherapy'};
    idd = (~ismember(Surgery,{'Yes'})) .* (ismember(Radiation,{'Yes'})) .* (ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Radiation + Chemotherapy'};
    idd = (ismember(Surgery,{'Yes'})) .* (ismember(Radiation,{'Yes'})) .* (ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'All'};
    idd = (~ismember(Surgery,{'Yes'})) .* (~ismember(Radiation,{'Yes'})) .* (~ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'No treatment'};
    
    % Treatment combination 
    disp('------------------------ Treatment combination  ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        var = Treatment(idd);
        tbl = tabulate(var);
        disp(tbl)
    end
    
    
    
    disp('=========================== Variable summary (each cluster, DP-Each) ===========================')
    LabelEst = LabelEst_DP_Each;
    ClusterNum = length(unique(LabelEst));
   
    
    % survival month
    disp('------------------------ survival month ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Survival months'});
        var = data_select_final(idd,ColId);
        var = table2array(var);
        disp(['     mean & std :  ',num2str(mean(var)),', ',num2str(std(var))])
        disp(['     min & max :  ',num2str(min(var)),', ',num2str(max(var))])
    end
    
    % censor indicator
    disp('------------------------ censor indicator ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Vital status recode (study cutoff used)'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        tbl = tabulate(var);
        disp(tbl)
    end
    
    % HPV status
    disp('------------------------ HPV status ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'HPV recode (2010-2017)'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        tbl = tabulate(var);
        disp(tbl)
    end
    
    
    % age
    disp('------------------------ age ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        if strcmp(AgeFlag,'Continuous') % (Continuous)
           ColId = ismember(VarNames,{'Age'});
           var = data_select_final(idd,ColId);
           var = table2array(var);
           disp(['     mean & std :  ',num2str(mean(var)),', ',num2str(std(var))])
        elseif strcmp(AgeFlag,'Discrete') % (groups)
           ColId = ismember(VarNames,{'Age Standard for Survival (15-44,45-54,55-64,65-74,75+)'});
           var = data_select_final(idd,ColId);
           var = table2cell(var);
           var(ismember(var,{'00 years','35-39 years','40-44 years',...
               '45-49 years','50-54 years'})) = {'<= 54'};
           var(ismember(var,{'55-59 years','60-64 years'})) = {'55-64'};
           var(ismember(var,{'65-69 years','70-74 years'})) = {'65-74'};
           var(ismember(var,{'75-79 years','80-84 years','85+ years'})) = {'75+'};
           tbl = tabulate(var);
           disp(tbl)
        end
    end
    
    % sex
    disp('------------------------ sex ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Sex'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        tbl = tabulate(var);
        disp(tbl)
    end
    
    % marital status
    disp('------------------------ marital status ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Marital status at diagnosis'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        var(ismember(var,{'Married (including common law)'})) = {'Married'};
        var(ismember(var,{'Single (never married)','Divorced','Widowed',...
            'Separated','Unmarried or Domestic Partner'}))... % ,'Unknown'
            = {'Unmarried'}; % /Unknown
        tbl = tabulate(var);
        disp(tbl)
    end
    
    % race 
    disp('------------------------ race ------------------------')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Race recode (White, Black, Other)'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        var(~ismember(var,{'White','Black'})) = {'Other/Unknown'};
        tbl = tabulate(var);
        disp(tbl)
    end
    
    
    % cancer site
    disp('------------------------ cancer site ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Primary Site - labeled'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        var(ismember(var,{'C09.9-Tonsil, NOS','C09.0-Tonsillar fossa',...
                'C10.9-Oropharynx, NOS','C09.1-Tonsillar pillar',...
                'C10.8-Overlapping lesion of oropharynx',...
                'C10.0-Vallecula','C10.3-Posterior wall of oropharynx',...
                'C09.8-Overlapping lesion of tonsil',...
                'C10.2-Lateral wall of oropharynx',...
                ''})) = {'Oropharynx'};
        var(~ismember(var,{'Oropharynx'})) = {'Non-Oropharynx'};
        tbl = tabulate(var);
        disp(tbl)
    end
    
    % summary stage
    disp('------------------------ summary stage ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Combined Summary Stage (2004+)'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        tbl = tabulate(var);
        disp(tbl)
    end
    

    % Surgery 
    disp('------------------------ Surgery ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Reason no cancer-directed surgery'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        var(~ismember(var,{'Surgery performed'})) = {'No/Unknown'};
        var(ismember(var,{'Surgery performed'})) = {'Yes'};
        tbl = tabulate(var);
        disp(tbl)
    end
    
    
    % Radiation 
    disp('------------------------ Radiation ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Radiation recode'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        var(~ismember(var,{'Beam radiation',...
            'Radiation, NOS  method or source not specified',...
            'Combination of beam with implants or isotopes',...
            'Radioactive implants (includes brachytherapy) (1988+)',...
            'Radioisotopes (1988+)'})) = {'No/Unknown'};
        var(ismember(var,{'Beam radiation',...
            'Radiation, NOS  method or source not specified',...
            'Combination of beam with implants or isotopes',...
            'Radioactive implants (includes brachytherapy) (1988+)',...
            'Radioisotopes (1988+)'})) = {'Yes'};
        tbl = tabulate(var);
        disp(tbl)
    end
    
    
    % Chemotherapy 
    disp('------------------------ Chemotherapy ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        ColId = ismember(VarNames,{'Chemotherapy recode (yes, no/unk)'});
        var = data_select_final(idd,ColId);
        var = table2cell(var);
        tbl = tabulate(var);
        disp(tbl)
    end
    
    % Treatment combination
    ColId = ismember(VarNames,{'Reason no cancer-directed surgery'});
    Surgery = data_select_final(:,ColId);
    Surgery = table2cell(Surgery);
    Surgery(ismember(Surgery,{'Surgery performed'})) = {'Yes'};
    ColId = ismember(VarNames,{'Radiation recode'});
    Radiation = data_select_final(:,ColId);
    Radiation = table2cell(Radiation);
    Radiation(ismember(Radiation,{'Beam radiation',...
            'Radiation, NOS  method or source not specified',...
            'Combination of beam with implants or isotopes',...
            'Radioactive implants (includes brachytherapy) (1988+)',...
            'Radioisotopes (1988+)'})) = {'Yes'};
    ColId = ismember(VarNames,{'Chemotherapy recode (yes, no/unk)'});
    Chemotherapy = data_select_final(:,ColId);
    Chemotherapy = table2cell(Chemotherapy);
    Treatment = Chemotherapy;
    idd = (ismember(Surgery,{'Yes'})) .* (~ismember(Radiation,{'Yes'})) .* (~ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Surgery only'};
    idd = (~ismember(Surgery,{'Yes'})) .* (ismember(Radiation,{'Yes'})) .* (~ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Radiation only'};
    idd = (~ismember(Surgery,{'Yes'})) .* (~ismember(Radiation,{'Yes'})) .* (ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Chemotherapy only'};
    idd = (ismember(Surgery,{'Yes'})) .* (ismember(Radiation,{'Yes'})) .* (~ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Surgery + Radiation'};
    idd = (ismember(Surgery,{'Yes'})) .* (~ismember(Radiation,{'Yes'})) .* (ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Surgery + Chemotherapy'};
    idd = (~ismember(Surgery,{'Yes'})) .* (ismember(Radiation,{'Yes'})) .* (ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'Radiation + Chemotherapy'};
    idd = (ismember(Surgery,{'Yes'})) .* (ismember(Radiation,{'Yes'})) .* (ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'All'};
    idd = (~ismember(Surgery,{'Yes'})) .* (~ismember(Radiation,{'Yes'})) .* (~ismember(Chemotherapy,{'Yes'})) ==1;
    Treatment(idd) = {'No treatment'};
    
    % Treatment combination 
    disp('------------------------ Treatment combination  ------------------------ ')
    for j=1:ClusterNum
        disp(['---> Cluster ',num2str(j)])
        idd = (LabelEst==j);
        var = Treatment(idd);
        tbl = tabulate(var);
        disp(tbl)
    end
    
    
    
end

toc


