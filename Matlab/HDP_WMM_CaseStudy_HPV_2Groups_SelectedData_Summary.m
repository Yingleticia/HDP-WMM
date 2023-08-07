%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% --  Hierarchical Dirichlet Process Weibull Mixture Model          -- %%%
%%% ----------    for prognosis analysis (cancer patients)    ---------- %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear;clc;close all;

seed = 3;
site = 'Oropharynx'; % All_Site; Oropharynx;
operation = 'RawVariable'; % RawVariable; ClusterSurvivalVariable; All
SurvivalFlag = 'Overall'; % Specific; Overall
AgeFlag = 'Discrete'; % Continuous; Discrete
EstBound = 1; % confidence interval for HDP-WMM estimation
EstBound_KM = 1; % confidence interval for KM estimation (in Survival)

folder = ['results/CaseStudy/HPV/2Groups/',site,'/seed_',num2str(seed),'/'];
L = 200; % for grid values in function estimates
% Unknown = 1; % 1: include Unknown; 0: exclude Unknown
year_extent = 5; % extented range (years)

%% load data
%---- Selected raw data ----%

filename = [folder,'original data selected_seed_',num2str(seed),'.xlsx'];
data_select_final = readtable(filename,'VariableNamingRule','preserve'); 
data_select_final_new = data_select_final; % included variables for Pearson's Chi-squared test
VarNames = data_select_final.Properties.VariableNames;
J_group = 2;
group_name = {'HPV Positive','HPV Negative'}; 

if strcmp(operation,'ClusterSurvivalVariable') || strcmp(operation,'All')
    
%---- Selected analyzed data ----%
filename = [folder,'results_seed_',num2str(seed),'_',SurvivalFlag,'.mat'];
load(filename,'data','J_group','Num_J')
% 1: survival time (1-2)
% 2: censor indicator (0: not censored; 1: censored)
% 3: group (restaurant)
% 4: original survival time
% 5: dish index
% 6: component parameter alpha
% 7: component parameter beta
% 8: table index

%---- Estimated clustering label ----%
filename = [folder,'ClusterLabelEst_HPV_',SurvivalFlag,'_seed_',num2str(seed),'.xlsx'];
ClusterLabelEst = readtable(filename,'VariableNamingRule','preserve');
LabelEst = ClusterLabelEst(:,1);
LabelEst = table2array(LabelEst);

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
            'Separated','Unmarried or Domestic Partner','Unknown'}))...
            = {'Unmarried/Unknown'};
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
    
    % AJCC stage, 7th 
%     disp('------------------------ AJCC stage, 7th ------------------------ ')
%     for j=1:J_group
%         disp(['---> ',group_name{j},':'])
%         idd = ismember(HPV,group_name{j});
%         ColId = ismember(VarNames,{'Derived AJCC Stage Group, 7th ed (2010-2015)'});
%         var = data_select_final(idd,ColId);
%         var = table2cell(var);
%         var(ismember(var,{'IVA','IVB','IVC','IVNOS'})) = {'IV'};
%         var(ismember(var,{'Blank(s)','NA','UNK Stage'})) = {'Unknown'};
%         tbl = tabulate(var);
%         disp(tbl)
%     end
    
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
    
    filename = [folder,'original data selected_seed_',num2str(seed),'_new.xlsx'];
    writetable(data_select_final_new,filename)
    
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

       
end


%% Clustering Survival
if strcmp(operation,'ClusterSurvivalVariable') || strcmp(operation,'All')
    
    % survival data
    sur_original = data(:,4);
    X_max = max(sur_original);
    sur = data(:,1);
    group = data(:,3);
    
    % censoring indicator (1) ----> Overall Survival
    censor = data(:,2);    
   
    % Weibull fit for each cluster (Overall Survival)
    disp('============================== Overall Survival ==============================')
    ClusterNum = length(unique(LabelEst));
    ParaEst = nan(ClusterNum,2);
    disp(['--- Total number of estimated clusters: ', num2str(ClusterNum)])
    disp('--------------- Cluster summary (All) --------------------')
    tbl = tabulate(LabelEst);
    tbl = array2table(tbl,'VariableNames',{'Cluster','Count','Frequency'});
    disp(tbl)
    for j=1:J_group
        disp(['--------------- Cluster summary (Group ',num2str(j),') --------------------'])
        idd = (group==j);
        Label_temp = LabelEst(idd);
        tbl = tabulate(Label_temp);
        tbl = array2table(tbl,'VariableNames',{'Cluster','Count','Frequency'});
        disp(tbl)
    end
    SurProb = nan(ClusterNum,3); % 1 year (12 month); 5 year (60 month); 10 year (120 month)
    fig_num = 1;
    figure(fig_num);
    color_order = get(gca,'colororder');
    if EstBound_KM == 1
        legend_name = cell(1,3*ClusterNum);
    else
        legend_name = cell(1,2*ClusterNum);
    end   
    linestyle = {'-','-.',':','--'}; 
    xx = linspace(0,X_max+12*year_extent,L); % grid values
    yy = 1 + xx/X_max; 
    for i=1:ClusterNum
        idd = (LabelEst==i);
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
        if EstBound_KM == 1
            legend_name{(i-1)*3+1} = ['Cluster ',num2str(i),'(HDP-WMM)'];
            legend_name{(i-1)*3+2} = ['Cluster ',num2str(i),'(KM)'];
            legend_name{(i-1)*3+3} = ['Cluster ',num2str(i),'(KM, 95% CI)'];
        else
            legend_name{(i-1)*2+1} = ['Cluster ',num2str(i),'(HDP-WMM)'];
            legend_name{(i-1)*2+2} = ['Cluster ',num2str(i),'(KM)'];
        end
        cc = 1 - (wblcdf(yy,parmHat(1),parmHat(2))-temp)/(1-temp);
        plot(xx,cc,'LineWidth', 2,'Color',color_order(i+J_group+1,:),'linestyle',linestyle{1}); hold on
        
        sur_original_temp = sur_original(idd);
        [ss_e,x,flo,fup] = ecdf(sur_original_temp,'Censoring',censor_temp,'function','survivor');
        stairs(x, ss_e,'LineWidth', 2,'Color',color_order(i+J_group+1,:),'linestyle',linestyle{2});
        if EstBound_KM == 1
            stairs(x, flo,'LineWidth', 1.5,'Color',color_order(i+J_group+1,:),'linestyle',linestyle{3})
            stairs(x, fup,'LineWidth', 1.5,'Color',color_order(i+J_group+1,:),'linestyle',linestyle{3},'HandleVisibility','off')
        end
    end
    set(gca,'FontSize',16,'FontName', 'times')
    xlim([0 max(xx)]); ylim([0 1]);
%     xlabel('Time (month)')
    if strcmp(SurvivalFlag,'Overall')
        name = 'Overall survival function';
    elseif strcmp(SurvivalFlag,'Specific')
        name = 'Disease-specific survival function';
    end
%     title(name,'FontSize',20,'FontName', 'times', 'FontWeight','normal')
    legend(legend_name,'FontSize',14,'FontName', 'times','Location','best') % 
    legend('boxoff')
    
    if strcmp(SurvivalFlag,'Overall') 
        name1 = 'Survival_Functions_Overall_Cluster.png';
        name2 = 'Survival_Functions_Overall_Cluster.fig';
    elseif strcmp(SurvivalFlag,'Specific') 
        name1 = 'Survival_Functions_Specific_Cluster.png';
        name2 = 'Survival_Functions_Specific_Cluster.fig';
    end
    saveas(gcf,[folder,name1]);
    saveas(gcf,[folder,name2]);
    
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
    
    disp('=========================== Variable summary (each cluster) ===========================')
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
            'Separated','Unmarried or Domestic Partner','Unknown'}))...
            = {'Unmarried/Unknown'};
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
    Treatment(idd) = {'Surgery + Radiation + Chemotherapy'};
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




