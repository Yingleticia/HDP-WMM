function [FunEstAll,xx,X_upper,KM_sur] = ...
    R1_DP_WMM_All_FunEst(SurvivalFlag,data_record,gamma_record,...
    rho_record,eta_record,epsilon,N,B,J_group,L,folder,...
    year_extent,medianEst,CI_method,CI_coef,Operation)

% row 1: mean/median; row 2: lower bound; row 3: upper bound; 
% column 1-4: density, suvivial, hazard, cumulative hazard
FunEstAll = cell(3,4); 

d = 2;
data = data_record(:,:,1);

%---------------- Bayesian clustering datasets ----------------%
ClusterDraws = nan(N, B);

%------------------------ Grid values -------------------------%
X_max = max(data(:,4));
if strcmp(Operation,'figure')
    
    X_upper = zeros(1, J_group);
    xx = zeros(L, J_group);
    yy = zeros(L, J_group);
    for j = 1:J_group
        X_upper(j) = X_max + 12*year_extent; % one more year 
        xx_temp = linspace(0,X_upper(j),L); % grid values
        xx_temp = xx_temp';
        xx(:,j) = xx_temp;
        yy(:,j) = 1 + xx_temp/X_max; 
    end 
    KM_sur = [];

elseif strcmp(Operation,'metric')
    
    X_upper = [];
    xx = cell(1, J_group);
    yy = cell(1, J_group);
    KM_sur = cell(1, J_group);
    sur_original = data(:,4);
    group = data(:,3);
    censor = data(:,2);
    for j = 1:J_group
        %--- KM ---%
        sur_temp = sur_original(group == j);
        censor_temp = censor(group == j);
        [KM_sur_j,KM_x_j,~,~] = ecdf(sur_temp,'Censoring',censor_temp,'function','survivor');
        KM_sur{j} = KM_sur_j;
        xx{j} = KM_x_j;
        yy{j} = 1 + KM_x_j/X_max; 
    end
    
end


%--------- Estimation ---------%
PDF_est = nan(L, J_group, B); % density
CDF_est = nan(L, J_group, B); 
Survival_est = nan(L, J_group, B);
hazard_est = nan(L, J_group, B);
Hazard_est = nan(L, J_group, B);  


for b = 1:B
    
    % b-th posterior samples
    data_b = data_record(:,:,b);
    gamma_b = max(gamma_record);
    rho_b = rho_record(b);
    eta_b = eta_record(b);
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
        
        if strcmp(Operation,'figure') 
            yy_j = yy(:,j);
            L_j = L;
        elseif strcmp(Operation,'metric')
            yy_j = yy{j};
            L_j = length(yy_j);
        end
        
        for i = 1:L_j
            x_i = yy_j(i);
            % pdf
            pdf_temp = (1/X_max) * weibull_pdf(x_i, G_alpha, G_beta)./temp;
            pdf_temp(idd) = 0;
            PDF_est(i,j,b) = sum(G_omega .* pdf_temp);
            % cdf
            cdf_temp = (weibull_cdf(x_i, G_alpha,G_beta) - temp1)./temp;
            cdf_temp(idd) = 0;
            CDF_est(i,j,b) = sum(G_omega .* cdf_temp);
            % survival
            Survival_est(i,j,b) = 1 - CDF_est(i,j,b);
            % hazard
            hazard_est(i,j,b) = PDF_est(i,j,b)/Survival_est(i,j,b);
            % Hazard
            Hazard_est(i,j,b) = -log(Survival_est(i,j,b));
        end
    end
    
       
end


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




ClusterDraws = ClusterDraws';
ClusterDraws = array2table(ClusterDraws);
if strcmp(SurvivalFlag,'Overall')
    filename = [folder,'Overall_ClusterDraws_DP_WMM_All.xlsx'];
elseif strcmp(SurvivalFlag,'Specific')
    filename = [folder,'Specific_ClusterDraws_DP_WMM_All.xlsx'];
end
writetable(ClusterDraws,filename,'WriteVariableNames',0)

end