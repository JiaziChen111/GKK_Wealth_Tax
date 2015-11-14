% Graphs GKK Wealth Tax
% Sergio Ocampo Diaz


 javaaddpath('jxl.jar');
 javaaddpath('MXL.jar');

 import mymxl.*;
 import jxl.*;   

%%  Parameters

% Grids
    n_a = 201 ; 
    n_z = 7   ;
    n_l = 5   ;
    n_e = 5   ;
    
    Max_Age = 81 ;
    Ret_Age = 45 ;

% Utility and technology
    sigma = 4.00  ;
    gamma = 0.4494; 
    mu    = 0.9   ;
    delta = 0     ;
    
% Taxes
    Threshold_Factor = 0.00 ;
	tauC    = 0.075         ;
    tauK    = 0.25          ;
	
% Set Tax Regime
    % tauPL = 0.185     ;
    % psi   = 0.77      ;
    tauPL = 0.0         ;
    psi   = 0.776       ;
    tauW  = 0.017072675596579098 ;
        
% Age brackets 
    age = [5 , 15, 25, 35, 45, 55, Max_Age ];
    %%%    24, 34, 44, 54, 64, 74, 80
    n_age = numel(age) ;
    
% Percentiles
prctl = [10, 25, 50, 75, 90, 99, 99.9];
    
%% Result_Folders

    if ((tauPL==0.0)&&(sigma~=1.0))  
        Result_Folder = strcat('../NSU_LT_Results/Factor_',num2str(Threshold_Factor,'%.2f'),'/') ;
        Bench_Folder  = '../NSU_LT_Results/Bench_Files/' ;
    elseif ((tauPL~=0.0)&&(sigma~=1.0))  
        Result_Folder = strcat('../NSU_PT_Results/Factor_',num2str(Threshold_Factor,'%.2f'),'/') ;
        Bench_Folder  = '../NSU_PT_Results/Bench_Files/' ;
    elseif ((tauPL==0.0)&&(sigma==1.0))  
        Result_Folder = strcat('../SU_LT_Results/Factor_',num2str(Threshold_Factor,'%.2f'),'/') ;
        Bench_Folder  = '../SU_LT_Results/Bench_Files/' ;
    elseif ((tauPL~=0.0)&&(sigma==1.0))  
        Result_Folder = strcat('../SU_PT_Results/Factor_',num2str(Threshold_Factor,'%.2f'),'/') ;
        Bench_Folder  = '../SU_PT_Results/Bench_Files/' ;
    end
    
    xls_file = 'AZ_Tables.xls' ;

%% Grids 
   
% A grid
    eval(['load ',Result_Folder,'agrid']);
    
    A_mat = repmat(agrid,[Max_Age,1,n_z,n_l,n_e]);
    
% Z grid
    eval(['load ',Result_Folder,'zgrid']);
    
    Z_mat = repmat(reshape(zgrid,[1,1,n_z,1,1]),[Max_Age,n_a,1,n_l,n_e]);
    
%% Read files - Benchmark Economy
% Distribution
    eval(['load ',Bench_Folder,'DBN'])
    DBN_bench = reshape(DBN,[Max_Age,n_a,n_z,n_l,n_e]) ;
    clear DBN
    
% Wage
    eval(['load ',Bench_Folder,'wage'])
    W_bench = wage ;
    clear wage
    
% Interst Rate
    eval(['load ',Bench_Folder,'rr'])
    R_bench = rr ;
    clear rr
    
% Interst Rate
    eval(['load ',Bench_Folder,'EBAR'])
    E_bench = EBAR ;
    clear EBAR
    
    
%% Read files - Experimental Economy

% Distribution
    eval(['load ',Result_Folder,'Exp_results_DBN']) ;
    DBN_exp = reshape(Exp_results_DBN,[Max_Age,n_a,n_z,n_l,n_e]) ;
    clear Exp_results_DBN

% Wage
    eval(['load ',Result_Folder,'Exp_results_wage'])
    W_exp = Exp_results_wage ;
    clear Exp_results_wage
    
% Interst Rate
    eval(['load ',Result_Folder,'Exp_results_rr'])
    R_exp = Exp_results_rr ;
    clear Exp_results_rr
    
% Wealth Taxes
    Threshold = Threshold_Factor*E_bench ;
    
    
%% Simulate the economy
% Load simulation resutls 
    eval(['load ',Result_Folder,'Sim_age_ben'])
    eval(['load ',Result_Folder,'Sim_A_ben'])
    eval(['load ',Result_Folder,'Sim_C_ben'])
    eval(['load ',Result_Folder,'Sim_H_ben'])
    eval(['load ',Result_Folder,'Sim_R_ben'])
    eval(['load ',Result_Folder,'Sim_R_at_ben'])
    eval(['load ',Result_Folder,'Sim_Z_ben'])
    eval(['load ',Result_Folder,'Sim_Ap_ben'])
    eval(['load ',Result_Folder,'Sim_Yh_ben'])
    Sim_S_ben = 100*(Sim_Ap_ben./Sim_A_ben-1) ;

    eval(['load ',Result_Folder,'Sim_age_exp'])
    eval(['load ',Result_Folder,'Sim_A_exp'])
    eval(['load ',Result_Folder,'Sim_C_exp'])
    eval(['load ',Result_Folder,'Sim_H_exp'])
    eval(['load ',Result_Folder,'Sim_R_exp'])
    eval(['load ',Result_Folder,'Sim_R_at_exp'])
    eval(['load ',Result_Folder,'Sim_Z_exp'])
    eval(['load ',Result_Folder,'Sim_Ap_exp'])
    eval(['load ',Result_Folder,'Sim_Yh_exp'])
    Sim_S_exp = 100*(Sim_Ap_exp./Sim_A_exp-1) ;
    
%%  Mean and std by age and age-z
    for i=1:Max_Age
        ind = (Sim_age_ben==i) ;
        n   = sum(ind)         ;
        sim_mean_A_ben(i)  = sum(Sim_A_ben(ind))/n  ;
        sim_mean_C_ben(i)  = sum(Sim_C_ben(ind))/n  ;
        sim_mean_H_ben(i)  = sum(Sim_H_ben(ind))/n  ;
        sim_mean_Ap_ben(i) = sum(Sim_Ap_ben(ind))/n ;
        sim_mean_S_ben(i)  = 100*(sim_mean_Ap_ben(i)/sim_mean_A_ben(i)-1) ;
        sim_mean_R_ben(i)     = sum(Sim_R_ben(ind))/n     ;
        sim_mean_R_at_ben(i)  = sum(Sim_R_at_ben(ind))/n  ;
        sim_mean_wR_ben(i)    = sum(Sim_R_ben(ind).*Sim_A_ben(ind))/sum(Sim_A_ben(ind))    ;
        sim_mean_wR_at_ben(i) = sum(Sim_R_at_ben(ind).*Sim_A_ben(ind))/sum(Sim_A_ben(ind)) ;
        sim_mean_S_ben_aux(i) = sum(Sim_S_ben(ind))/n ;
        
        sim_std_A_ben(i)  = ( sum( (Sim_A_ben(ind) - sim_mean_A_ben(i)).^2 )/(n-1) )^0.5 ;
        sim_std_C_ben(i)  = ( sum( (Sim_C_ben(ind) - sim_mean_C_ben(i)).^2 )/(n-1) )^0.5 ;
        sim_std_H_ben(i)  = ( sum( (Sim_H_ben(ind) - sim_mean_H_ben(i)).^2 )/(n-1) )^0.5 ;
        sim_std_Ap_ben(i) = ( sum( (Sim_Ap_ben(ind) - sim_mean_Ap_ben(i)).^2 )/(n-1) )^0.5 ;
        sim_std_R_ben(i)     = ( sum( (Sim_R_ben(ind) - sim_mean_R_ben(i)).^2 )/(n-1) )^0.5 ;
        sim_std_R_at_ben(i)  = ( sum( (Sim_R_at_ben(ind) - sim_mean_R_at_ben(i)).^2 )/(n-1) )^0.5 ;
        sim_std_wR_ben(i)    = ( sum( (Sim_R_ben(ind) - sim_mean_wR_ben(i)).^2 .*Sim_A_ben(ind) )/sum(Sim_A_ben(ind)) )^0.5 ;
        sim_std_wR_at_ben(i) = ( sum( (Sim_R_at_ben(ind) - sim_mean_wR_at_ben(i)).^2 .*Sim_A_ben(ind) )/sum(Sim_A_ben(ind)) )^0.5 ;
        
        
        ind = (Sim_age_exp==i) ;
        n   = sum(ind)         ;
        sim_mean_A_exp(i)  = sum(Sim_A_exp(ind))/n  ;
        sim_mean_C_exp(i)  = sum(Sim_C_exp(ind))/n  ;
        sim_mean_H_exp(i)  = sum(Sim_H_exp(ind))/n  ;
        sim_mean_Ap_exp(i) = sum(Sim_Ap_exp(ind))/n ;
        sim_mean_S_exp(i)  = 100*(sim_mean_Ap_exp(i)/sim_mean_A_exp(i)-1) ;
        sim_mean_R_exp(i)     = sum(Sim_R_exp(ind))/n     ;
        sim_mean_R_at_exp(i)  = sum(Sim_R_at_exp(ind))/n  ;
        sim_mean_wR_exp(i)    = sum(Sim_R_exp(ind).*Sim_A_exp(ind))/sum(Sim_A_exp(ind))    ;
        sim_mean_wR_at_exp(i) = sum(Sim_R_at_exp(ind).*Sim_A_exp(ind))/sum(Sim_A_exp(ind)) ;
        sim_mean_S_exp_aux(i) = sum(Sim_S_exp(ind))/n ;
        
        sim_std_A_exp(i)  = ( sum( (Sim_A_exp(ind) - sim_mean_A_exp(i)).^2 )/(n-1) )^0.5 ;
        sim_std_C_exp(i)  = ( sum( (Sim_C_exp(ind) - sim_mean_C_exp(i)).^2 )/(n-1) )^0.5 ;
        sim_std_H_exp(i)  = ( sum( (Sim_H_exp(ind) - sim_mean_H_exp(i)).^2 )/(n-1) )^0.5 ;
        sim_std_Ap_exp(i) = ( sum( (Sim_Ap_exp(ind) - sim_mean_Ap_exp(i)).^2 )/(n-1) )^0.5 ;
        sim_std_R_exp(i)     = ( sum( (Sim_R_exp(ind) - sim_mean_R_exp(i)).^2 )/(n-1) )^0.5 ;
        sim_std_R_at_exp(i)  = ( sum( (Sim_R_at_exp(ind) - sim_mean_R_at_exp(i)).^2 )/(n-1) )^0.5 ;
        sim_std_wR_exp(i)    = ( sum( (Sim_R_exp(ind) - sim_mean_wR_exp(i)).^2 .*Sim_A_exp(ind) )/sum(Sim_A_exp(ind)) )^0.5 ;
        sim_std_wR_at_exp(i) = ( sum( (Sim_R_at_exp(ind) - sim_mean_wR_at_exp(i)).^2 .*Sim_A_exp(ind) )/sum(Sim_A_exp(ind)) )^0.5 ;
        
    
        for j=1:n_z
            ind = ((Sim_age_ben==i).*(Sim_Z_ben==j))==1 ;
            n   = sum(ind)                         ;
            sim_mean_A_ben_AZ(i,j)  = sum(Sim_A_ben(ind))/n  ;
            sim_mean_C_ben_AZ(i,j)  = sum(Sim_C_ben(ind))/n  ;
            sim_mean_H_ben_AZ(i,j)  = sum(Sim_H_ben(ind))/n  ;
            sim_mean_Ap_ben_AZ(i,j) = sum(Sim_Ap_ben(ind))/n ;
            sim_mean_S_ben_AZ(i,j)  = 100*(sim_mean_Ap_ben_AZ(i,j)/sim_mean_A_ben_AZ(i,j)-1) ;
            sim_mean_R_ben_AZ(i,j)     = sum(Sim_R_ben(ind))/n     ;
            sim_mean_R_at_ben_AZ(i,j)  = sum(Sim_R_at_ben(ind))/n  ;
            sim_mean_wR_ben_AZ(i,j)    = sum(Sim_R_ben(ind).*Sim_A_ben(ind))/sum(Sim_A_ben(ind))    ;
            sim_mean_wR_at_ben_AZ(i,j) = sum(Sim_R_at_ben(ind).*Sim_A_ben(ind))/sum(Sim_A_ben(ind)) ;

            sim_std_A_ben_AZ(i,j)  = ( sum( (Sim_A_ben(ind) - sim_mean_A_ben_AZ(i,j)).^2 )/(n-1) )^0.5 ;
            sim_std_C_ben_AZ(i,j)  = ( sum( (Sim_C_ben(ind) - sim_mean_C_ben_AZ(i,j)).^2 )/(n-1) )^0.5 ;
            sim_std_H_ben_AZ(i,j)  = ( sum( (Sim_H_ben(ind) - sim_mean_H_ben_AZ(i,j)).^2 )/(n-1) )^0.5 ;
            sim_std_Ap_ben_AZ(i,j) = ( sum( (Sim_Ap_ben(ind) - sim_mean_Ap_ben_AZ(i,j)).^2 )/(n-1) )^0.5 ;
            sim_std_R_ben_AZ(i,j)     = ( sum( (Sim_R_ben(ind) - sim_mean_R_ben_AZ(i,j)).^2 )/(n-1) )^0.5 ;
            sim_std_R_at_ben_AZ(i,j)  = ( sum( (Sim_R_at_ben(ind) - sim_mean_R_at_ben_AZ(i,j)).^2 )/(n-1) )^0.5 ;
            sim_std_wR_ben_AZ(i,j)    = ( sum( (Sim_R_ben(ind) - sim_mean_wR_ben_AZ(i,j)).^2 .*Sim_A_ben(ind) )/sum(Sim_A_ben(ind)) )^0.5 ;
            sim_std_wR_at_ben_AZ(i,j) = ( sum( (Sim_R_at_ben(ind) - sim_mean_wR_at_ben_AZ(i,j)).^2 .*Sim_A_ben(ind) )/sum(Sim_A_ben(ind)) )^0.5 ;

            ind = ((Sim_age_exp==i).*(Sim_Z_exp==j))==1 ;
            n   = sum(ind)                         ;
            sim_mean_A_exp_AZ(i,j)  = sum(Sim_A_exp(ind))/n  ;
            sim_mean_C_exp_AZ(i,j)  = sum(Sim_C_exp(ind))/n  ;
            sim_mean_H_exp_AZ(i,j)  = sum(Sim_H_exp(ind))/n  ;
            sim_mean_Ap_exp_AZ(i,j) = sum(Sim_Ap_exp(ind))/n ;
            sim_mean_S_exp_AZ(i,j)  = 100*(sim_mean_Ap_exp_AZ(i,j)/sim_mean_A_exp_AZ(i,j)-1) ;
            sim_mean_R_exp_AZ(i,j)     = sum(Sim_R_exp(ind))/n     ;
            sim_mean_R_at_exp_AZ(i,j)  = sum(Sim_R_at_exp(ind))/n  ;
            sim_mean_wR_exp_AZ(i,j)    = sum(Sim_R_exp(ind).*Sim_A_exp(ind))/sum(Sim_A_exp(ind))    ;
            sim_mean_wR_at_exp_AZ(i,j) = sum(Sim_R_at_exp(ind).*Sim_A_exp(ind))/sum(Sim_A_exp(ind)) ;

            sim_std_A_exp_AZ(i,j)  = ( sum( (Sim_A_exp(ind) - sim_mean_A_exp_AZ(i,j)).^2 )/(n-1) )^0.5 ;
            sim_std_C_exp_AZ(i,j)  = ( sum( (Sim_C_exp(ind) - sim_mean_C_exp_AZ(i,j)).^2 )/(n-1) )^0.5 ;
            sim_std_H_exp_AZ(i,j)  = ( sum( (Sim_H_exp(ind) - sim_mean_H_exp_AZ(i,j)).^2 )/(n-1) )^0.5 ;
            sim_std_Ap_exp_AZ(i,j) = ( sum( (Sim_Ap_exp(ind) - sim_mean_Ap_exp_AZ(i,j)).^2 )/(n-1) )^0.5 ;
            sim_std_R_exp_AZ(i,j)     = ( sum( (Sim_R_exp(ind) - sim_mean_R_exp_AZ(i,j)).^2 )/(n-1) )^0.5 ;
            sim_std_R_at_exp_AZ(i,j)  = ( sum( (Sim_R_at_exp(ind) - sim_mean_R_at_exp_AZ(i,j)).^2 )/(n-1) )^0.5 ;
            sim_std_wR_exp_AZ(i,j)    = ( sum( (Sim_R_exp(ind) - sim_mean_wR_exp_AZ(i,j)).^2 .*Sim_A_exp(ind) )/sum(Sim_A_exp(ind)) )^0.5 ;
            sim_std_wR_at_exp_AZ(i,j) = ( sum( (Sim_R_at_exp(ind) - sim_mean_wR_at_exp_AZ(i,j)).^2 .*Sim_A_exp(ind) )/sum(Sim_A_exp(ind)) )^0.5 ; 
        end
    end

    
%% 
    for j=1:n_age
        
        
        for i=1:numel(prctl)
            if j==1
                ind = (Sim_age_ben<age(j))  ;
            elseif j<n_age && j>1
                ind = ((Sim_age_ben>=age(j-1)).*(Sim_age_ben<age(j)))==1  ;
            else
                ind = (Sim_age_ben>=age(j-1))  ;
            end
            prc_A_ben(j,i)    = prctile(Sim_A_ben(ind),prctl(i))  ;
            prc_C_ben(j,i)    = prctile(Sim_C_ben(ind),prctl(i))  ;
            prc_H_ben(j,i)    = prctile(Sim_H_ben(ind),prctl(i))  ;
            prc_Ap_ben(j,i)   = prctile(Sim_Ap_ben(ind),prctl(i)) ;
            prc_R_ben(j,i)    = prctile(Sim_R_ben(ind),prctl(i))    ;
            prc_R_at_ben(j,i) = prctile(Sim_R_at_ben(ind),prctl(i)) ;
            
            if j==1
                ind = (Sim_age_exp<age(j))  ;
            elseif j<n_age && j>1
                ind = ((Sim_age_exp>=age(j-1)).*(Sim_age_exp<age(j)))==1  ;
            else
                ind = (Sim_age_exp>=age(j-1))  ;
            end
            prc_A_exp(j,i)    = prctile(Sim_A_exp(ind),prctl(i))  ;
            prc_C_exp(j,i)    = prctile(Sim_C_exp(ind),prctl(i))  ;
            prc_H_exp(j,i)    = prctile(Sim_H_exp(ind),prctl(i))  ;
            prc_Ap_exp(j,i)   = prctile(Sim_Ap_exp(ind),prctl(i)) ;
            prc_R_exp(j,i)    = prctile(Sim_R_exp(ind),prctl(i))    ;
            prc_R_at_exp(j,i) = prctile(Sim_R_at_exp(ind),prctl(i)) ;
        end
        
        for i=1:n_z
            if j==1
                ind = ((Sim_age_ben<age(j)).*(Sim_Z_ben==i))==1  ;
            elseif j<n_age && j>1
                ind = ((Sim_age_ben>=age(j-1)).*(Sim_age_ben<age(j)).*(Sim_Z_ben==i))==1  ;
            else
                ind = ((Sim_age_ben>=age(j-1)).*(Sim_Z_ben==i))==1  ;
            end
            n   = sum(ind)                              ;
            A_AZ_ben(j,i)   = sum(Sim_A_ben(ind))/n  ;
            C_AZ_ben(j,i)   = sum(Sim_C_ben(ind))/n  ;
            H_AZ_ben(j,i)   = sum(Sim_H_ben(ind))/n  ;
            Ap_AZ_ben(j,i)  = sum(Sim_Ap_ben(ind))/n ;
            S_AZ_ben(j,i)   = 100*(Ap_AZ_ben(j,i)/A_AZ_ben(j,i)-1) ;
            R_AZ_ben(j,i)     = sum(Sim_R_ben(ind))/n     ;
            R_at_AZ_ben(j,i)  = sum(Sim_R_at_ben(ind))/n  ;
            wR_AZ_ben(j,i)    = sum(Sim_R_ben(ind).*Sim_A_ben(ind))/sum(Sim_A_ben(ind))    ;
            wR_at_AZ_ben(j,i) = sum(Sim_R_at_ben(ind).*Sim_A_ben(ind))/sum(Sim_A_ben(ind)) ;
        
            if j==1
                ind = ((Sim_age_exp<age(j)).*(Sim_Z_exp==i))==1  ;
            elseif j<n_age && j>1
                ind = ((Sim_age_exp>=age(j-1)).*(Sim_age_exp<age(j)).*(Sim_Z_exp==i))==1  ;
            else
                ind = ((Sim_age_exp>=age(j-1)).*(Sim_Z_exp==i))==1  ;
            end
            n   = sum(ind)                              ;
            A_AZ_exp(j,i)   = sum(Sim_A_exp(ind))/n  ;
            C_AZ_exp(j,i)   = sum(Sim_C_exp(ind))/n  ;
            H_AZ_exp(j,i)   = sum(Sim_H_exp(ind))/n  ;
            Ap_AZ_exp(j,i)  = sum(Sim_Ap_exp(ind))/n ;
            S_AZ_exp(j,i)   = 100*(Ap_AZ_exp(j,i)/A_AZ_exp(j,i)-1) ;
            R_AZ_exp(j,i)     = sum(Sim_R_exp(ind))/n     ;
            R_at_AZ_exp(j,i)  = sum(Sim_R_at_exp(ind))/n  ;
            wR_AZ_exp(j,i)    = sum(Sim_R_exp(ind).*Sim_A_exp(ind))/sum(Sim_A_exp(ind))    ;
            wR_at_AZ_exp(j,i) = sum(Sim_R_at_exp(ind).*Sim_A_exp(ind))/sum(Sim_A_exp(ind)) ;
        end
    end 

%% Wealth Distribution
    n_w = 1000 ;
    wealth_grid = linspace(log(min(Sim_A_ben)),log(max(Sim_A_ben)),n_w)' ;
    Pareto_ben = NaN(n_w,1) ;
    Pareto_exp = NaN(n_w,1) ;
    N = numel(Sim_A_ben) ;
    log_A_ben = log(Sim_A_ben) ;
    log_A_exp = log(Sim_A_exp) ;
    for i=1:n_w
       Pareto_ben(i) = log(sum(log_A_ben>=wealth_grid(i))/N) ;
       Pareto_exp(i) = log(sum(log_A_exp>=wealth_grid(i))/N) ;
    end
    
    figure;
    plot(wealth_grid,[Pareto_ben Pareto_exp]); xlim([wealth_grid(1),wealth_grid(end)])
    xlabel('Log(Wealth)'); legend('Bench','Exp','location','northeast')
    print('-dpdf','./Simulation/Wealth_Pareto.pdf') ;
    
    prctl_vec = linspace(1,99,99) ;
    prctl_W_ben = NaN(99,1) ;
    prctl_W_exp = NaN(99,1) ;
    for i=1:99
        prctl_W_ben(i) = prctile(log_A_ben,prctl_vec(i)) ;
        prctl_W_exp(i) = prctile(log_A_exp,prctl_vec(i)) ;
    end
    
    figure;
    subplot(2,1,1); plot(prctl_vec,[prctl_W_ben prctl_W_exp]);
    legend('Bench','Exp','location','northeast'); title('Percentiles of Wealth Dist.')
    subplot(2,1,2); plot(prctl_vec,(prctl_W_exp-prctl_W_ben));
    xlabel('percentiles'); title('Change in Percentiles of Wealth Dist.')
    print('-dpdf','./Simulation/Wealth_Prctl.pdf') ;

%% Simulation Graphs  
    z_vec = [2 4 6] ;

    figure; 
        subplot(1,4,1); plot(20:100,[sim_mean_A_ben' sim_mean_A_exp']); xlim([19+1,19+Max_Age]); title('Mean Assets')
        subplot(1,4,2); plot(20:80,[sim_mean_S_ben(1:end-20)' sim_mean_S_exp(1:end-20)']); xlim([19+1,19+Max_Age-20]); title('Mean Saving Rate')
        subplot(1,4,3); plot(20:100,[sim_mean_C_ben' sim_mean_C_exp']); xlim([19+1,19+Max_Age]); title('Mean Consumption')
        subplot(1,4,4); plot(20:100,[sim_mean_H_ben' sim_mean_H_exp']); xlim([19+1,19+Max_Age]); title('Mean Hours')
        legend('Bench','Exp','location','southeast')
    print('-dpdf','./Simulation/Mean_Variables_Life_Cycle.pdf') ;
    
    figure; 
        subplot(2,4,1); plot(20:100,sim_mean_A_ben_AZ(:,z_vec)'); xlim([19+1,19+Max_Age]); title('Mean Assets - bench')
        subplot(2,4,2); plot(20:80,sim_mean_S_ben_AZ(1:end-20,z_vec)'); xlim([19+1,19+Max_Age-20]); title('Mean Saving Rate - bench')
        subplot(2,4,3); plot(20:100,sim_mean_C_ben_AZ(:,z_vec)'); xlim([19+1,19+Max_Age]); title('Mean Consumption - bench')
        subplot(2,4,4); plot(20:100,sim_mean_H_ben_AZ(:,z_vec)'); xlim([19+1,19+Max_Age]); title('Mean Hours - bench')
        legend('z2','z4','z6','location','southeast')
        subplot(2,4,5); plot(20:100,sim_mean_A_exp_AZ(:,z_vec)'); xlim([19+1,19+Max_Age]); title('Mean Assets - exp')
        subplot(2,4,6); plot(20:80,sim_mean_S_exp_AZ(1:end-20,z_vec)'); xlim([19+1,19+Max_Age-20]); title('Mean Saving Rate - exp')
        subplot(2,4,7); plot(20:100,sim_mean_C_exp_AZ(:,z_vec)'); xlim([19+1,19+Max_Age]); title('Mean Consumption - exp')
        subplot(2,4,8); plot(20:100,sim_mean_H_exp_AZ(:,z_vec)'); xlim([19+1,19+Max_Age]); title('Mean Hours - exp')
        legend('z2','z4','z6','location','southeast')
    print('-dpdf','./Simulation/Mean_Variables_Life_Cycle_by_Z.pdf') ;
    
    figure; 
        subplot(1,3,1); plot(20:100,[sim_std_A_ben' sim_std_A_exp']); xlim([19+1,19+Max_Age]); title('Std Assets')
        subplot(1,3,2); plot(20:100,[sim_std_C_ben' sim_std_C_exp']); xlim([19+1,19+Max_Age]); title('Std Consumption')
        subplot(1,3,3); plot(20:100,[sim_std_H_ben' sim_std_H_exp']); xlim([19+1,19+Max_Age]); title('Std Hours')
        legend('Bench','Exp','location','southeast')
    print('-dpdf','./Simulation/Std_Variables_Life_Cycle.pdf') ;
    
    figure; 
        subplot(2,3,1); plot(20:100,sim_std_A_ben_AZ(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Assets - bench')
        subplot(2,3,2); plot(20:100,sim_std_C_ben_AZ(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Consumption - bench')
        subplot(2,3,3); plot(20:100,sim_std_H_ben_AZ(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Hours - bench')
        legend('z2','z4','z6','location','southeast')
        subplot(2,3,4); plot(20:100,sim_std_A_exp_AZ(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Assets - exp')
        subplot(2,3,5); plot(20:100,sim_std_C_exp_AZ(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Consumption - exp')
        subplot(2,3,6); plot(20:100,sim_std_H_exp_AZ(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Hours - exp')
        legend('z2','z4','z6','location','southeast')
    print('-dpdf','./Simulation/Std_Variables_Life_Cycle_by_Z.pdf') ;
    
    figure; 
        subplot(2,2,1); plot(20:100,[100*(sim_mean_wR_ben-1)' 100*(sim_mean_wR_exp-1)']); xlim([19+1,19+Max_Age]); title('Mean Return')
        subplot(2,2,2); plot(20:100,[sim_std_wR_ben' sim_std_wR_exp']); xlim([19+1,19+Max_Age]); title('Std Return')
        subplot(2,2,3); plot(20:100,[100*(sim_mean_wR_at_ben-1)' 100*(sim_mean_wR_at_exp-1)']); xlim([19+1,19+Max_Age]); title('Mean Return - After tax')
        subplot(2,2,4); plot(20:100,[sim_std_wR_at_ben' sim_std_wR_at_exp']); xlim([19+1,19+Max_Age]); title('Std Return - After tax')
        legend('Bench','Exp','location','southeast')
    print('-dpdf','./Simulation/Return_Life_Cycle.pdf') ;
    
    figure; 
        subplot(2,2,1); plot(20:100,100*(sim_mean_wR_ben_AZ(:,z_vec)-1)); xlim([19+1,19+Max_Age]); title('Mean Return - bench')
        subplot(2,2,2); plot(20:100,sim_std_wR_ben_AZ(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Return - bench')
        legend('z2','z4','z6','location','southeast')
        subplot(2,2,3); plot(20:100,100*(sim_mean_wR_exp_AZ(:,z_vec)-1)); xlim([19+1,19+Max_Age]); title('Mean Return - exp')
        subplot(2,2,4); plot(20:100,sim_std_wR_exp_AZ(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Return - exp')
        legend('z2','z4','z6','location','southeast')
    print('-dpdf','./Simulation/Return_Life_Cycle_by_Z.pdf') ;

    
%% Simulation Tables
    prc_file = './Simulation/prc_Tables.xls' ;
    AZ_file  = './Simulation/AZ_Tables.xls'  ;

% Titles 
    % Weights by z
        Weights_Z_bench = NaN(1,n_z) ;
        Weights_Z_exp   = NaN(1,n_z) ;
        for i=1:n_z
            Weights_Z_bench(i) = sum(sum(sum(sum(DBN_bench(:,:,i,:,:))))) ;
            Weights_Z_exp(i)   = sum(sum(sum(sum(DBN_exp(:,:,i,:,:)))))   ;
        end 

    % Weights by age
        Weights_Age_bench = zeros(n_age,1) ;
        Weights_Age_exp   = zeros(n_age,1) ;
        j = 1;
        for i=1:Max_Age
            if i<=age(j)
            Weights_Age_bench(j) = Weights_Age_bench(j) + sum(sum(sum(sum(DBN_bench(i,:,:,:,:))))) ;
            Weights_Age_exp(j)   = Weights_Age_exp(j)   + sum(sum(sum(sum(DBN_exp(i,:,:,:,:))))) ;
            end 
            if i==age(j); j=j+1; end 
        end 

    % Weights by age-Z
        Weights_AZ_bench = zeros(n_age,n_z) ;
        Weights_AZ_exp   = zeros(n_age,n_z) ;
        for i=1:n_z
            j=1;
            for l=1:Max_Age
                if l<=age(j)
                Weights_AZ_bench(j,i) = Weights_AZ_bench(j,i) + sum(sum(sum(sum(DBN_bench(l,:,i,:,:))))) ;
                Weights_AZ_exp(j,i)   = Weights_AZ_exp(j,i)   + sum(sum(sum(sum(DBN_exp(l,:,i,:,:))))) ;
                end 
                if l==age(j); j=j+1; end 
            end
        end 

    Weights_Z   = 100*Weights_Z_bench   ; 
    Weights_Age = 100*Weights_Age_bench ; 
    Weights_AZ  = 100*Weights_AZ_bench  ;
    clear Weights_Z_bench Weights_Z_exp Weights_Age_bench Weights_Age_exp Weights_AZ_bench Weights_AZ_exp
    
    z_title{1} = ' ' ;
    for i=1:n_z
        z_title{i+1} = strcat('z',int2str(i),'(',num2str(Weights_Z(i),'%2.1f'),')') ;
    end 
    
    
    age_title{1,1} = strcat('<25 (',num2str(Weights_Age(1),'%2.1f'),')') ;
    age_title{2,1} = strcat('25-34 (',num2str(Weights_Age(2),'%2.1f'),')') ;
    age_title{3,1} = strcat('35-44 (',num2str(Weights_Age(3),'%2.1f'),')') ;
    age_title{4,1} = strcat('45-54 (',num2str(Weights_Age(4),'%2.1f'),')') ;
    age_title{5,1} = strcat('55-64 (',num2str(Weights_Age(5),'%2.1f'),')') ;
    age_title{6,1} = strcat('65-74 (',num2str(Weights_Age(6),'%2.1f'),')') ;
    age_title{7,1} = strcat('>75 (',num2str(Weights_Age(7),'%2.1f'),')') ;
    
% Percentiles 
    p_title{1} = ' ' ;
    for i=1:n_z
        p_title{i+1} = strcat('p',num2str(prctl(i))) ;
    end 
    prc_R_ben = 100*(prc_R_ben-1) ;
    prc_R_exp = 100*(prc_R_exp-1) ;
    prc_R_at_ben = 100*(prc_R_at_ben-1) ;
    prc_R_at_exp = 100*(prc_R_at_exp-1) ;
    
    Mat_A_ben = [{'prc_A_ben',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(prc_A_ben)] ;
    Mat_C_ben = [{'prc_C_ben',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(prc_C_ben)] ;
    Mat_H_ben = [{'prc_H_ben',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(prc_H_ben)] ;
    Mat_R_ben = [{'prc_R_ben',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(prc_R_ben)] ;
    Mat_R_at_ben = [{'prc_R_at_ben',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(prc_R_at_ben)] ;
    
    Mat_A_exp = [{'prc_A_exp',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(prc_A_exp)] ;
    Mat_C_exp = [{'prc_C_exp',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(prc_C_exp)] ;
    Mat_H_exp = [{'prc_H_exp',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(prc_H_exp)] ;
    Mat_R_exp = [{'prc_R_exp',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(prc_R_exp)] ;
    Mat_R_at_exp = [{'prc_R_at_exp',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(prc_R_at_exp)] ;
    
    Mat_A_diff = [{'prc_A_diff',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(100*(prc_A_exp./prc_A_ben-1))] ;
    Mat_C_diff = [{'prc_C_diff',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(100*(prc_C_exp./prc_C_ben-1))] ;
    Mat_H_diff = [{'prc_H_diff',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(100*(prc_H_exp./prc_H_ben-1))] ;
    Mat_R_diff = [{'prc_R_diff',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(prc_R_exp-prc_R_ben)] ;
    Mat_R_at_diff = [{'prc_R_at_diff',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(prc_R_at_exp-prc_R_at_ben)] ;
    
   
    Mat_7 = cell(1,8*3+3) ;
    Mat_8 = cell(9,1) ;
    
    Mat = [Mat_7;Mat_8 Mat_A_ben Mat_8 Mat_A_exp Mat_8 Mat_A_diff] ;
    status = xlwrite(prc_file,Mat,'Assets') ;
    
    Mat = [Mat_7;Mat_8 Mat_C_ben Mat_8 Mat_C_exp Mat_8 Mat_C_diff] ;
    status = xlwrite(prc_file,Mat,'Cons') ;
    
    Mat = [Mat_7;Mat_8 Mat_H_ben Mat_8 Mat_H_exp Mat_8 Mat_H_diff] ;
    status = xlwrite(prc_file,Mat,'Hours') ;
    
    Mat = [Mat_7;Mat_8 Mat_R_ben Mat_8 Mat_R_exp Mat_8 Mat_R_diff] ;
    status = xlwrite(prc_file,Mat,'Return') ;
    
    Mat = [Mat_7;Mat_8 Mat_R_at_ben Mat_8 Mat_R_at_exp Mat_8 Mat_R_at_diff] ;
    status = xlwrite(prc_file,Mat,'Return - AT') ;
    
% AZ Tables
    Mat_A_ben = [{'AZ_A_ben',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(A_AZ_ben)] ;
    Mat_C_ben = [{'AZ_C_ben',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(C_AZ_ben)] ;
    Mat_H_ben = [{'AZ_H_ben',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(H_AZ_ben)] ;
    Mat_S_ben = [{'AZ_S_ben',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(S_AZ_ben)] ;
    Mat_wR_ben = [{'AZ_wR_ben',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(100*(wR_AZ_ben-1))] ;
    Mat_wR_at_ben = [{'AZ_wR_at_ben',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(100*(wR_at_AZ_ben-1))] ;
    
    Mat_A_exp = [{'AZ_A_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(A_AZ_exp)] ;
    Mat_C_exp = [{'AZ_C_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(C_AZ_exp)] ;
    Mat_H_exp = [{'AZ_H_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(H_AZ_exp)] ;
    Mat_S_exp = [{'AZ_S_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(S_AZ_exp)] ;
    Mat_wR_exp = [{'AZ_wR_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(100*(wR_AZ_exp-1))] ;
    Mat_wR_at_exp = [{'AZ_wR_at_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(100*(wR_at_AZ_exp-1))] ;
    
    Mat_A_diff = [{'AZ_A_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(100*(A_AZ_exp./A_AZ_ben-1))] ;
    Mat_C_diff = [{'AZ_C_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(100*(C_AZ_exp./C_AZ_ben-1))] ;
    Mat_H_diff = [{'AZ_H_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(100*(H_AZ_exp./H_AZ_ben-1))] ;
    Mat_S_diff = [{'AZ_S_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(S_AZ_exp-S_AZ_ben)] ;
    Mat_wR_diff = [{'AZ_wR_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(wR_AZ_exp-wR_AZ_ben)] ;
    Mat_wR_at_diff = [{'AZ_wR_at_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(wR_at_AZ_exp-wR_at_AZ_ben)] ;
    
    Mat = [Mat_7;Mat_8 Mat_A_ben Mat_8 Mat_A_exp Mat_8 Mat_A_diff] ;
    status = xlwrite(AZ_file,Mat,'Assets') ;
    
    Mat = [Mat_7;Mat_8 Mat_S_ben Mat_8 Mat_S_exp Mat_8 Mat_S_diff] ;
    status = xlwrite(AZ_file,Mat,'Saving Rate') ;
    
    Mat = [Mat_7;Mat_8 Mat_C_ben Mat_8 Mat_C_exp Mat_8 Mat_C_diff] ;
    status = xlwrite(AZ_file,Mat,'Cons') ;
    
    Mat = [Mat_7;Mat_8 Mat_H_ben Mat_8 Mat_H_exp Mat_8 Mat_H_diff] ;
    status = xlwrite(AZ_file,Mat,'Hours') ;
    
    Mat = [Mat_7;Mat_8 Mat_wR_ben Mat_8 Mat_wR_exp Mat_8 Mat_wR_diff] ;
    status = xlwrite(AZ_file,Mat,'Return') ;
    
    Mat = [Mat_7;Mat_8 Mat_wR_at_ben Mat_8 Mat_wR_at_exp Mat_8 Mat_wR_at_diff] ;
    status = xlwrite(AZ_file,Mat,'Return - AT') ;

    
%% Wealth stats 

N = numel(Sim_A_ben) ;
% Wealth of people with more tna 10.000 units of assets 
    n = sum(Sim_A_ben>=5000) ;
    Av_Wealth = sum(Sim_A_ben.*(Sim_A_ben>=5000))/n                 ;
    Wealth_Share = sum(Sim_A_ben.*(Sim_A_ben>=5000))/sum(Sim_A_ben) ;
    col={'Wealth','Share of Pop','Av. Wealth','Wealth Share'} ;
    Mat=[col;num2cell([5000,100*n/N,Av_Wealth,100*Wealth_Share])] ;
    disp(Mat)

% Wealth of people with more tna 10.000 units of assets 
    n = sum(Sim_A_ben>=10000) ;
    Av_Wealth = sum(Sim_A_ben.*(Sim_A_ben>=10000))/n                 ;
    Wealth_Share = sum(Sim_A_ben.*(Sim_A_ben>=10000))/sum(Sim_A_ben) ;
    col={'Wealth','Share of Pop','Av. Wealth','Wealth Share'} ;
    Mat=[col;num2cell([10000,100*n/N,Av_Wealth,100*Wealth_Share])] ;
    disp(Mat)
    
% Wealth of people with more tna 50.000 units of assets 
    n = sum(Sim_A_ben>=50000) ;
    Av_Wealth = sum(Sim_A_ben.*(Sim_A_ben>=50000))/n                 ;
    Wealth_Share = sum(Sim_A_ben.*(Sim_A_ben>=50000))/sum(Sim_A_ben) ;
    col={'Wealth','Share of Pop','Av. Wealth','Wealth Share'} ;
    Mat=[col;num2cell([50000,100*n/N,Av_Wealth,100*Wealth_Share])] ;
    disp(Mat)
    
% Wealth of top ten wealth holders
    Wealth_top_10 = sort(Sim_A_ben);
    Wealth_top_10 = sum(Wealth_top_10(end-10:end)) ;
    Wealth_top_10_Share = 100* Wealth_top_10 / sum(Sim_A_ben) ; 
    
% Average Labor income
    Av_labor_income_no_ret = sum(Sim_Yh_ben.*Sim_age_ben<Ret_Age.*Sim_H_ben>0)/sum(Sim_age_ben<Ret_Age.*Sim_H_ben>0) ;
    
% Ratio of wealth of top ten earners to labor income:
    Ratio_top10Wealth_Labor_income = (Wealth_top_10/10)/Av_labor_income_no_ret ;
    

    
%% Return by asset distribution 

    for j=1:n_age
            
        for i=1:numel(prctl)
            if j==1
                ind = (Sim_age_ben<age(j))  ;
            elseif j<n_age && j>1
                ind = ((Sim_age_ben>=age(j-1)).*(Sim_age_ben<age(j)))==1  ;
            else
                ind = (Sim_age_ben>=age(j-1))  ;
            end
            Z_aux = Sim_Z_ben(ind) ;
            prc_A_ben(j,i)    = prctile(Sim_A_ben(ind),prctl(i))  ;
            a_ind             = find(Sim_A_ben(ind)==prc_A_ben(j,i),1) ;
%             prc_Z_ben(j,i)    = Z_aux(a_ind) ;
            %prc_R_ben(j,i)    = 1 + ( R_bench * (Z_aux(a_ind))^mu * prc_A_ben(j,i)^(mu-1)  - delta) ;
            %prc_R_ben(j,i)    = 1 + ( R_bench * (Z_aux(a_ind))^mu * prc_A_ben(j,i)^(mu-1)  - delta)*(1-tauK) ;
            
            if j==1
                ind = (Sim_age_exp<age(j))  ;
            elseif j<n_age && j>1
                ind = ((Sim_age_exp>=age(j-1)).*(Sim_age_exp<age(j)))==1  ;
            else
                ind = (Sim_age_exp>=age(j-1))  ;
            end
            Z_aux = Sim_Z_exp(ind) ;
            prc_A_exp(j,i)    = prctile(Sim_A_ben(ind),prctl(i))  ;
            a_ind             = find(Sim_A_ben(ind)==prc_A_ben(j,i),1) ;
%             prc_Z_exp(j,i)    = Z_aux(a_ind) ;
%             prc_R_exp(j,i)    = 1 + ( R_exp * (Z_aux(a_ind))^mu * prc_A_exp(j,i)^(mu-1)  - delta) ;
%             prc_R_exp(j,i)    = (1 + ( R_exp * (Z_aux(a_ind))^mu * prc_A_exp(j,i)^(mu-1)  - delta))*(1-tauW) ;
        end
        
    end 


    