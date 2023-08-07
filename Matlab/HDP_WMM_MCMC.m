
%%% performing MCMC sampling for HDP-WMM

function [data_record_HDP, gamma_record_HDP, alpha_0_record_HDP,...
    rho_record_HDP, eta_record_HDP] = HDP_WMM_MCMC(data,d,a_rho,b_rho,...
    a_eta,b_eta,a_gamma,b_gamma,gamma,a_alpha_0,b_alpha_0,alpha_0,...
    TableNum,J_group,Num_J,N,B,iterSize,burnin,gap,numDiscrete)

% storage
data_record_HDP = zeros(size(data,1),size(data,2),B);
gamma_record_HDP = zeros(1,B);
alpha_0_record_HDP = zeros(1,B);
rho_record_HDP = zeros(1,B);
eta_record_HDP = zeros(1,B);
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
        data_record_HDP(:,:,b) = data;
        gamma_record_HDP(b) = gamma;
        alpha_0_record_HDP(b) = alpha_0;
        rho_record_HDP(b) = rho;
        eta_record_HDP(b) = eta;
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


end