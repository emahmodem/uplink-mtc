function [ output_args ] = gen_coverage_results_mMTC_old( input_args )
%GEN_COVERAGE_RESULTS_MMTC-D  Summary of this function goes here
%   Author : Mahmoud Kamel
%   Date :  12-Dec-2017 12:01:39
%   This function generates results of simulating the effect of clustering
%   on the load of cells and coverage
%/////////////////////////////////////////////////////////////////////////
clc;close all;clear all;
%hold on;
% Simulation Parameters
params.Pm = 100e-3;                          % Uplink Transmit power of M2M nodes in watts 20 dBm
params.Pmin = 1e-6;                          % Min. received power at the base station in watts > BS sensitivity
params.simulation_area_side = [-500 500];    % simulation area side to side
params.space_realizations = 100;
params.time_slots = 10;
params.BW = 20e6;
params.N_RB = 100;
%params.N = 10;    % Orthogonal carriers allocated to the RACH transmission of M2M traffic
params.RB_BW = 180e3; %Hz
params.No = 10^(-17.4) * params.RB_BW;
params.SEPL.alpha = 0.94;
params.SEPL.beta = 0.5;

%/////////////////////////////////////////////////////////////////////////
% Simulation Switches
%sim = 'Debug';
%sim = 'CoverageThreshold_Pmin';
%sim = 'CoverageThreshold_SEPL';
%sim = 'CoverageThreshold_Pm';
%sim = 'CellDensity_Pmin';

%sim = 'MaxTransmitPower';
%sim = 'PowerTruncationThreshold';
%sim = 'ResourceBlocks'
%sim = 'Validate_PDF_OF_P';

%sim = 'Expectation_Log_P';
%sim = 'Expectation_Log__K_P';
%sim = 'CellDensity';
%sim = 'M2MDensity';
%sim = 'ClusterHeadDensity';

switch(sim)
    case 'Debug'
        params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 10;
        params.time_slots = 10;
        Pm_dB = 20;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        params.rho_m = 0.1;                                    % of M2M nodes traffic activity
        params.aggregation_mode = 'DIRECT';
        Threshold_dB = 0 ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        params.rho_m = 0.1;
        params.LA_B = [100 300 500 800 1000 3000 5000 8000 10000]*1e-6   ;       % Small Cell Denisty (cells/m^2)
        
        
        Pmin_dBm = -90 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_smallcells_density(params);
        end

        subplot(1,2,2);
        f1 = plot(params.LA_B*1e6  ,Pcov_analy , 'b-');
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        xlabel('Small cell density (cells/$km^2$)  ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$P_{min} = -90$ dBm  ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
    case 'CoverageThreshold_Pmin'
        params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 10;
        params.time_slots = 10;
        params.LA_B = 1000e-6   ;       % Small Cell Denisty (cells/m^2)
        %params.Beta_Distribution = [3 4];   % Beta Distribution parameters                                   % of M2M nodes traffic activity
        Threshold_dB = (-20:1:20) ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        LA_M = [1e-2 5e-2 1e-1];
        params.rho_m = 0.1;
        Pmin_dBm = -70 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            [Pcov_approx_analy(i,:) , Pcov_exact_analy(i,:)] = compute_uplink_coverage_with_coverage_threshold(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_coverage_threshold(params);
        end
        figure(1)
        subplot(1,2,1);
        f1 = plot(Threshold_dB  ,Pcov_exact_analy , 'k-',Threshold_dB,Pcov_simul,'ro' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Coverage threshold ($\tau$)  ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$P_{min} = -70$ dBm  ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
        Pmin_dBm = -80 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            [Pcov_approx_analy(i,:) , Pcov_exact_analy(i,:)] = compute_uplink_coverage_with_coverage_threshold(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_coverage_threshold(params);
        end
        subplot(1,2,2);
        f1 = plot(Threshold_dB  ,Pcov_exact_analy , 'k-',Threshold_dB,Pcov_simul,'ro' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Coverage threshold ($\tau$)  ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$P_{min} = -80$ dBm  ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'CoverageThreshold_SEPL'
        params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 10;
        params.time_slots = 10;
        params.LA_B = 1000e-6   ;       % Small Cell Denisty (cells/m^2)
        %params.Beta_Distribution = [3 4];   % Beta Distribution parameters                                   % of M2M nodes traffic activity
        Threshold_dB = (-20:1:20) ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        LA_M = [1e-2 5e-2 1e-1];
        params.rho_m = 0.1;
        Pmin_dBm = -70 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            [Pcov_approx_analy(i,:) , Pcov_exact_analy(i,:)] = compute_uplink_coverage_with_coverage_threshold(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_coverage_threshold(params);
        end
        figure(1)
        subplot(1,2,1);
        f1 = plot(Threshold_dB  ,Pcov_exact_analy , 'k-',Threshold_dB,Pcov_simul,'ro' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Coverage threshold ($\tau$)  ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$\alpha = 0.94, \beta = \frac{1}{2}$','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
        params.SEPL.alpha = 0.3;
        params.SEPL.beta = 2/3;
        
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            [Pcov_approx_analy(i,:) , Pcov_exact_analy(i,:)] = compute_uplink_coverage_with_coverage_threshold(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_coverage_threshold(params);
        end
        subplot(1,2,2);
        f1 = plot(Threshold_dB  ,Pcov_exact_analy , 'k-',Threshold_dB,Pcov_simul,'ro' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Coverage threshold ($\tau$)  ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$\alpha = 0.3, \beta = \frac{2}{3}$ ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'CoverageThreshold_Pm'
        params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 10;
        params.time_slots = 10;
        params.LA_B = 1000e-6   ;       % Small Cell Denisty (cells/m^2)
        %params.Beta_Distribution = [3 4];   % Beta Distribution parameters                                   % of M2M nodes traffic activity
        Threshold_dB = (-20:1:20) ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        
        params.rho_m = 0.1;
        Pmin_dBm = -80 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        
        Pm_dB = 20;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            [Pcov_approx_analy(i,:) , Pcov_exact_analy(i,:)] = compute_uplink_coverage_with_coverage_threshold(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_coverage_threshold(params);
        end
        figure(1)
        subplot(1,2,1);
        f1 = plot(Threshold_dB  ,Pcov_exact_analy , 'k-',Threshold_dB,Pcov_simul,'ro' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Coverage threshold ($\tau$)  ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$P_m = 20$ dBm','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
        Pm_dB = 33;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            [Pcov_approx_analy(i,:) , Pcov_exact_analy(i,:)] = compute_uplink_coverage_with_coverage_threshold(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_coverage_threshold(params);
        end
        subplot(1,2,2);
        f1 = plot(Threshold_dB  ,Pcov_exact_analy , 'k-',Threshold_dB,Pcov_simul,'ro' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Coverage threshold ($\tau$)  ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$P_m = 33$ dBm ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'CellDensity_Pmin'
        params.simulation_area_side = [-50 50];    % simulation area side to side
        params.space_realizations = 10;
        params.time_slots = 10;
        Pm_dB = 20;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        params.rho_m = 0.1;                                    % of M2M nodes traffic activity
        params.aggregation_mode = 'DIRECT';
        Threshold_dB = 0 ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        
        %params.LA_B = [100 300 500 800 1000 3000 5000 8000 10000]*1e-6   ;       % Small Cell Denisty (cells/m^2)
        params.LA_B = [100 300 500]*1e-6;
        %
        Pmin_dBm = -70 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_smallcells_density(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_small_cell_density(params);
        end
        figure(1)
        subplot(1,2,1);
        f1 = plot(params.LA_B*1e6  ,Pcov_analy , 'k-',params.LA_B*1e6,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Small cell density (cells/$km^2$)  ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$P_{min} = -70$ dBm  ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
        Pmin_dBm = -90 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_smallcells_density(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_small_cell_density(params);
        end
        figure(1)
        subplot(1,2,2);
        f1 = plot(params.LA_B*1e6  ,Pcov_analy , 'k-',params.LA_B*1e6,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Small cell density (cells/$km^2$)  ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$P_{min} = -90$ dBm  ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'ResourceBlocks'
        params.N_RB = 10:10:100;
        params.LA_B = 1000e-6   ;       % Small Cell Denisty (cells/m^2)
        params.Beta_Distribution = [3 4];   % Beta Distribution parameters                                   % of M2M nodes traffic activity
        Threshold_dB = 0 ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        LA_M = [1e-2 5e-2 1e-1];
        params.rho_m = 0.1;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_resource_blocks(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_resource_blocks(params);
        end
        figure(1)
        %subplot(1,2,1);
        f1 = plot(params.N_RB  ,Pcov_analy , 'k-',params.N_RB,Pcov_simul,'ro' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Number of orthogonal resource blocks($N_{RB}$)  ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
    case 'PowerTruncationThreshold'
        Pmin = -100:10:0; % in dBm
        params.Pmin = 10.^(Pmin/10) * 1e-3;
        params.LA_B = 1000e-6   ;       % Small Cell Denisty (cells/m^2
        Threshold_dB = 0 ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        LA_M = [1e-2 5e-2 1e-1];
        params.rho_m = 0.1;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_power_truncation_threshold(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_power_truncation_threshold(params);
        end
        figure(1)
        %subplot(1,2,1);
        f1 = plot(Pmin  ,Pcov_analy , 'k-',Pmin,Pcov_simul,'ro' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel(' Power truncation threshold $P_{min}$(dBm) ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
    case 'MaxTransmitPower'
        Pm = 0:5:33; % in dBm
        params.Pm = 10.^(Pm/10) * 1e-3;
        params.LA_B = 1000e-6   ;       % Small Cell Denisty (cells/m^2
        Threshold_dB = 0 ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        LA_M = [1e-2 5e-2 1e-1];
        params.rho_m = 0.1;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_max_transmit_power(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_max_transmit_power(params);
        end
        figure(1)
        %subplot(1,2,1);
        f1 = plot(Pm  ,Pcov_analy , 'k-',Pm,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Maximum transmit power $P_m$(dBm)','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
    case 'Expectation_Log__K_P'
        params.LA_B = [10:10:100 200:100:1000 2000:1000:10000] .* 1e-6    ;
        K = 1;
        [Ep] = compute_Exep_of_log_K_P(params,K);
        Ep_ = compute_Exep_of_log_P(params);
        f1 = semilogx((params.LA_B * 1e6) , Ep_ , 'b^--' , (params.LA_B * 1e6) , Ep , 'k-' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend({'$E_p(\ln(P)) $' '$E_p(\ln^k(P) , k = 1)$' , '$\ln(P_o)$'},'FontSize',25,'FontWeight','bold','Location','northeast','Interpreter','LaTex');
        
        xlabel('Small cells density ($\lambda_s)~ cells/km^2$   ','Interpreter','LaTex');
        ylabel('$E_P(\ln^k(P))$ ','Interpreter','LaTex');
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
    case 'Expectation_Log_P'
        params.LA_B = [10:10:100 200:100:1000 2000:1000:10000] .* 1e-6    ;
        Ep = compute_Exep_of_log_P(params);
        f1 = semilogx((params.LA_B * 1e6) , Ep , 'k-' , (params.LA_B * 1e6) , log(params.Pmin)*ones(1,numel(params.LA_B)), 'r--');
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend({'$E_p(\ln(P))$' , '$\ln(P_o)$'},'FontSize',25,'FontWeight','bold','Location','northeast','Interpreter','LaTex');
        
        xlabel('Small cells density ($\lambda_s)~ cells/km^2$   ','Interpreter','LaTex');
        ylabel('$E_P(\ln(P))$ ','Interpreter','LaTex');
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
    case 'Validate_PDF_OF_P'
        params.LA_B = 1000e-6   ;
        params.LA_M = 1;
        params.Beta_Distribution = [3 4];   % Beta Distribution parameters
        [f_analy ,  x] = compute_pdf_of_P(params);
        [f_simul , y , P ]= simulate_pdf_of_P(params);
        f1 = plot(x,f_analy, 'k-' ,y,f_simul , 'rx' );
        hold on;
        
        xlim([0 0.0005])
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1(1:2),{'$Analysis$' , '$Simulation$'},'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('$x$  ','Interpreter','LaTex');
        ylabel('$f_P(x)$ ','Interpreter','LaTex');
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
        %edges = params.Pmin:params.Pmin/10:params.Pm;
        %          [n , x] = hist(P,1000);
        %          n_ = n./numel(P) * 100;
        %           bar(x,n_);
        
        
        
    case 'Loading'
        params.LA_B = [100e-6 1000e-6 5000e-6 10000e-6]   ;       % Small Cell Denisty (cells/m^2)
        params.LA_H = 500e-6;     % H2H Users Density(users/m^2)
        %params.LA_M = 5e-1;       % M2M Nodes Density(nodes/m^2)
        params.Beta_Distribution = [3 4];   % Beta Distribution parameters
        % of M2M nodes traffic activity
        params.aggregation_mode = 'DIRECT';
        LA_M = [1e-1 5e-1 1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            [loads_dir(i,:),~] = compute_average_load_with_small_cell_density(params);
        end
        
        params.aggregation_mode = 'AGGREGATION';
        params.fraction_of_cluster_heads = 0.05;
        params.LA_C = params.fraction_of_cluster_heads * params.LA_M;       % M2M Nodes Density(nodes/m^2)
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            params.fraction_of_cluster_heads = 0.01;
            params.LA_C = params.fraction_of_cluster_heads * params.LA_M;       % M2M Nodes Density(nodes/m^2)
            [loads_agg(i,:),loads_h2h(i,:),loads_ch(i,:)] = compute_average_load_with_small_cell_density(params);
        end
        
        figure(1)
        subplot(1,2,1);
        f1 = semilogx(params.LA_B ,loads_dir);
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1(1:3),{'$\lambda_m  = 10^5~nodes/km^2$' , '$\lambda_m  = 5 \times 10^5~nodes/km^2$' , '$\lambda_m  = 10^6~nodes/km^2$'},'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
        ylabel('Cell load (M2M nodes / cell)','Interpreter','LaTex');
        title('No aggregation')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
        subplot(1,2,2);
        f2 = semilogx(params.LA_B ,loads_agg);
        set(f2,'MarkerSize',15);
        set(f2,'LineWidth',4);
        legend(f2(1:3),{'$\lambda_m  = 10^5~nodes/km^2$' , '$\lambda_m  = 5 \times 10^5~nodes/km^2$' , '$\lambda_m  = 10^6~nodes/km^2$'},'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
        ylabel('Cell load (M2M nodes / cell)','Interpreter','LaTex');
        title('Aggregation')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
        
        figure(2)
        subplot(1,2,1);
        f3 = semilogx(params.LA_B ,loads_h2h);
        set(f3,'MarkerSize',15);
        set(f3,'LineWidth',4);
        legend(f3(1:3),{'$\lambda_m  = 10^5~nodes/km^2$' , '$\lambda_m  = 5 \times 10^5~nodes/km^2$' , '$\lambda_m  = 10^6~nodes/km^2$'},'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
        ylabel('H2H aggregation load (M2M nodes / H2H user)','Interpreter','LaTex');
        title('H2H aggregator')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
        subplot(1,2,2);
        f4 = semilogx(params.LA_B ,loads_ch);
        set(f4,'MarkerSize',15);
        set(f4,'LineWidth',4);
        legend(f4(1:3),{'$\lambda_m  = 10^5~nodes/km^2$' , '$\lambda_m  = 5 \times 10^5~nodes/km^2$' , '$\lambda_m  = 10^6~nodes/km^2$'},'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
        ylabel('CH aggregation load (M2M nodes / CH user)','Interpreter','LaTex');
        title('CH aggregator')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
    case 'NetworkVisualize'
        params.LA_B = 1000e-6;       % Small Cell Denisty (cells/m^2)
        params.LA_H = 500e-6;     % H2H Users Density(users/m^2)
        params.LA_M = 5e-1;       % M2M Nodes Density(nodes/m^2)
        params.Beta_Distribution = [3 4];   % Beta Distribution parameters
        % of M2M nodes traffic activity
        params.fraction_of_cluster_heads = 0.01;
        params.LA_C = params.fraction_of_cluster_heads * params.LA_M;       % M2M Nodes Density(nodes/m^2)
        params.aggregation_mode = 'AGGREGATION';
        visualize_the_network(params)
        
end

end

%/////////////////////////////////////////////////////////////////////////
function [Ep] = compute_Exep_of_log_K_P(params,K)
Ep = zeros(K,numel(params.LA_B));
a = params.SEPL.alpha;
b = params.SEPL.beta;
Pu= params.Pm;
Po = params.Pmin;
theta = log((Pu/Po).^(1/a)).^(1/b);
n = (2/b) - 1;
A = (1 -  exp(- pi .* params.LA_B .* theta^2) ).^-1;
for k = 0:K
    Ep(k+1,:) =  a^k .* gamma (k/(n+1) + 1) .* gammainc( pi .* params.LA_B .* theta^2 , k/(n+1) + 1 ,'lower') ./ ((pi .* params.LA_B) .^(k/(n+1)) .* (1 - exp(- pi .* params.LA_B .* theta^2)));
end
end
%/////////////////////////////////////////////////////////////////////////
function [Ep ] = compute_Exep_of_log_P(params)
a = params.SEPL.alpha;
b = params.SEPL.beta;
Pu= params.Pm;
Po = params.Pmin;
theta = log((Pu/Po).^(1/a)).^(1/b);
% Ep = log(params.Pmin) +  (a .* gamma (b/2 + 1) .* gammainc( pi .* params.LA_B .* theta^2 , b/2 + 1 ,'lower')) ./ ((pi .* params.LA_B) .^(b/2) .*(1 - exp(- pi .* params.LA_B .* theta^2)));
Ep = (a .* gamma (b/2 + 1) .* gammainc( pi .* params.LA_B .* theta^2 , b/2 + 1 ,'lower')) ./ (((pi .* params.LA_B) .^(b/2)) .* (1 - exp(- pi .* params.LA_B .* theta^2)));
end
%/////////////////////////////////////////////////////////////////////////
function [f_analy ,  x]= compute_pdf_of_P(params)
a = params.SEPL.alpha;
b = params.SEPL.beta;
Pu= params.Pm;
Po = params.Pmin;
LA_B = params.LA_B;
theta = log((Pu/Po).^(1/a)).^(1/b);
x = params.Pmin:params.Pmin:params.Pm;
f_analy = (2 * pi * LA_B) / (a * b * (1 - exp(- pi * LA_B .* theta.^2))) * (1 ./ x) .* (log((x/Po).^(1/a)).^(2/b - 1)) .* exp(- pi * LA_B .* log((x/Po).^(1/a)).^(2/b));
end
%/////////////////////////////////////////////////////////////////////////
function [f_simul , x , P] = simulate_pdf_of_P(params)
simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.LA_B);
P = [];
for p = 1:points
    
    fprintf('\n')
    disp(['Smallcell Density: ' , num2str(params.LA_B(p))]);
    disp(['M2M Density: ' , num2str(params.LA_M)]);
    disp(['Path Loss Parameters: ' ,'alpha=', num2str(params.SEPL.alpha) ,' beta=', num2str(params.SEPL.beta)]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_b = params.LA_B(p) * simulation_area;
        rho_m = betarnd(params.Beta_Distribution(1),params.Beta_Distribution(2));
        mu_m = params.LA_M * simulation_area;
        
        N_cells = poissrnd(mu_b);
        while(N_cells < 1)
            N_cells = poissrnd(mu_b);
        end
        N_users_M2M = round(rho_m * poissrnd(mu_m));
        
        locations.BS = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        locations.M2M = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users_M2M, 2);
        
        % Association
        distances_M2M_to_BS = pdist2(locations.BS,locations.M2M,'euclidean') ;
        [distance_to_server , server_m2m] = min(distances_M2M_to_BS);
        
        
        server_m2m_P = params.Pmin * exp(params.SEPL.alpha .* (distance_to_server .^ params.SEPL.beta));
        P = [P server_m2m_P];
        
        
        
        
    end
    
end
pts = params.Pmin:params.Pmin:params.Pm;
%pts = params.Pmin:params.Pmin/20:params.Pm;
factor = params.space_realizations * params.time_slots;
[f_simul ,x] = ksdensity(P,pts);
%f_simul = f_simul./factor;
end
%/////////////////////////////////////////////////////////////////////////
function P_cov = compute_uplink_coverage_with_coverage_threshold(params)
SNR = params.Pmin  / params.No;
noise_term = exp(- params.Threshold / SNR);
a = params.SEPL.alpha;
b = params.SEPL.beta;
t = params.Threshold;
n = (2/b) - 1;
F = @(y,v,tau) log(tau ./ y).^(v)  ./ ( y + 1);
Pu= params.Pm;
Po = params.Pmin;
theta = log((Pu/Po).^(1/a)).^(1/b);
A = (n+1) * (pi / a^(n+1));
B = params.rho_m * (params.LA_M / params.N_RB);
for k = 0:n
    Ep(k+1) =  a^k .* gamma (k/(n+1) + 1) .* gammainc( pi .* params.LA_B .* theta^2 , k/(n+1) + 1 , 'lower') ./ ((pi .* params.LA_B) .^(k/(n+1)) .* (1 - exp(- pi .* params.LA_B .* theta^2)));
    NK(k+1) = nchoosek(n,k);
    for p = 1:numel(t)
        J_v(k+1,p) = integral(@(y)F(y,n - k,t(p)),0,t(p));
    end
    T(k+1,:) = NK(k+1) * Ep(k+1) * J_v(k+1,:);
end

S = sum(T,1);
interference_term = exp(- A * B * S);
P_cov = noise_term .* interference_term;

end
%/////////////////////////////////////////////////////////////////////////
function Pcov = simulate_uplink_coverage_with_coverage_threshold(params)
simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.Threshold);
PoCoverage = zeros(points,  params.space_realizations , params.time_slots );

fprintf('\n')
%disp(['coverage Threshold: ' , num2str(params.Threshold(p))]);
disp(['Smallcell Density: ' , num2str(params.LA_B)]);
disp(['M2M Density: ' , num2str(params.LA_M)]);
disp(['Path Loss Parameters: ' ,'alpha=', num2str(params.SEPL.alpha) ,' beta=', num2str(params.SEPL.beta)]);
disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
for m = 1:params.space_realizations;
    if(mod(m,params.space_realizations/100) == 0)
        fprintf('|');
    end
    mu_b = params.LA_B * simulation_area;
    %rho_m = betarnd(params.Beta_Distribution(1),params.Beta_Distribution(2));
    mu_m = params.LA_M * simulation_area;
    
    N_cells = poissrnd(mu_b);
    while(N_cells < 1)
        N_cells = poissrnd(mu_b);
    end
    N_users_M2M = round(params.rho_m * poissrnd(mu_m));
    
    locations.BS = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
    locations.M2M = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users_M2M, 2);
    
    % Association
    distances_M2M_to_BS = pdist2(locations.BS,locations.M2M,'euclidean') ;
    [distance_to_server , server_m2m] = min(distances_M2M_to_BS);
    
    server_m2m(2,:)= randi(params.N_RB,1,N_users_M2M);   % each M2M node selects a carrier from N carriers randomly (uniform distribution)
    server_m2m_P = params.Pmin * exp(params.SEPL.alpha .* distance_to_server .^ params.SEPL.beta);
    server_m2m_P(server_m2m_P > params.Pm) = 0;
    %server_m2m(3,:) = server_m2m_P;
    
    S_m = zeros(N_users_M2M,1);
    I_m = zeros(N_users_M2M,1);
    SINR_m = zeros(N_users_M2M,1);
    
    for t = 1:params.time_slots
        Hm = exprnd(1,N_cells,N_users_M2M);
        for i = 1:N_users_M2M
            ser = server_m2m(1,i);
            rb = server_m2m(2,i);
            same_ser =  (server_m2m(1,:) == ser);
            another_ser =  (server_m2m(1,:) ~= ser);
            same_rb = (server_m2m(2,:) == rb);
            intra_interferers = find(same_ser + same_rb == 2);
            inter_interferers = find(another_ser + same_rb == 2);
            S_m(i) = params.Pmin * Hm(ser,i);
            I_intra =  sum(params.Pmin * Hm(ser,intra_interferers));
            distance_to_interferers = pdist2(locations.BS(ser,:),locations.M2M(inter_interferers,:),'euclidean') ;
            I_inter = sum( server_m2m_P(inter_interferers) .* Hm(ser,inter_interferers).* exp(-params.SEPL.alpha .* distance_to_interferers .^ params.SEPL.beta) );
            %I_m(i) = I_intra + I_inter;
            I_m(i) =  I_inter;  % special case #1 no intra interference
        end
        SINR_m = S_m./ (I_m + params.No);
        SINR_m_dB = 10*log10(SINR_m);
        for p = 1:points
            PoCoverage(p,m,t) = sum(SINR_m > params.Threshold(p)) / N_users_M2M;
        end
    end
    
end

normfact = params.space_realizations * params.time_slots ;% * simulation_area;
Pcov = sum(sum(PoCoverage,3),2) / normfact;


end
%/////////////////////////////////////////////////////////////////////////
function P_exact = compute_uplink_coverage_with_smallcells_density(params)
SNR = params.Pmin  / params.No;
noise_term = exp(- params.Threshold / SNR);
a = params.SEPL.alpha;
b = params.SEPL.beta;
Pu= params.Pm;
Po = params.Pmin;
n = (2/b) - 1;
t = params.Threshold;

theta = log((Pu./Po).^(1/a)).^(1/b);
F = @(y,v,tau) log(tau ./ y).^(v)  ./ ( y + 1);
% compute E_P(ln^k(P))
%A = (1 -  exp(- pi .* params.LA_B .* theta^2) ).^-1;
A = (n+1) * (pi / a^(n+1));
B = params.rho_m * (params.LA_M / params.N_RB);
for k = 0:n
    Ep(k+1,:) =  a^k .* gamma (k/(n+1) + 1) .* gammainc( pi .* params.LA_B .* theta.^2 , k/(n+1) + 1 , 'lower') ./ ((pi .* params.LA_B) .^(k/(n+1)) .* (1 - exp(- pi .* params.LA_B .* theta.^2)));
    NK(k+1) = nchoosek(n,k);
    J_v(k+1) = integral(@(y)F(y,n - k,t),0,t);
    T(k+1,:) = NK(k+1) *  J_v(k+1) * Ep(k+1,:) ;
end

S = sum(T,1);
interference_term = exp(- A * B * S);
P_exact = noise_term .* interference_term;
end

function Pcov = simulate_uplink_coverage_with_small_cell_density(params)
simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.LA_B);
%PoCoverage = zeros(points,  params.space_realizations , params.time_slots );
PoCoverage = zeros(params.space_realizations , params.time_slots );


for p = 1:points
    if(params.LA_B(p) > 5000)
        params.simulation_area_side = [-350 350];
        simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
    end
    fprintf('\n')
    %disp(['Done: ' , num2str(params.LA_B(p))]);
    disp(['Smallcell Density: ' , num2str(params.LA_B(p))]);
    disp(['M2M Density: ' , num2str(params.LA_M)]);
    disp(['Path Loss Parameters: ' ,'alpha=', num2str(params.SEPL.alpha) ,' beta=', num2str(params.SEPL.beta)]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_b = params.LA_B(p) * simulation_area;
        %rho_m = betarnd(params.Beta_Distribution(1),params.Beta_Distribution(2));
        mu_m = params.LA_M * simulation_area;
        
        N_cells = poissrnd(mu_b);
        while(N_cells < 1)
            N_cells = poissrnd(mu_b);
        end
        N_users_M2M = round(params.rho_m * poissrnd(mu_m));
        
        locations.BS = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        locations.M2M = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users_M2M, 2);
        
        % Association
        distances_M2M_to_BS = pdist2(locations.BS,locations.M2M,'euclidean') ;
        [distance_to_server , server_m2m] = min(distances_M2M_to_BS);
        
        server_m2m(2,:)= randi(params.N_RB,1,size(server_m2m,2));  % each M2M node selects a carrier from N carriers randomly (uniform distribution)
        server_m2m_P = params.Pmin * exp(params.SEPL.alpha .* distance_to_server .^ params.SEPL.beta);
        server_m2m_P(server_m2m_P > params.Pm) = 0;
        %server_m2m(3,:) = server_m2m_P;
        
        S_m = zeros(N_users_M2M,1);
        I_m = zeros(N_users_M2M,1);
        SIR_m = zeros(N_users_M2M,1);
        
        for t = 1:params.time_slots
            Hm = exprnd(1,N_cells,N_users_M2M);
            for i = 1:N_users_M2M
                ser = server_m2m(1,i);
                rb = server_m2m(2,i);
                same_ser =  (server_m2m(1,:) == ser);
                another_ser =  (server_m2m(1,:) ~= ser);
                same_rb = (server_m2m(2,:) == rb);
                %intra_interferers = find(same_ser + same_rb == 2);
                inter_interferers = find(another_ser + same_rb == 2);
                S_m(i) = params.Pmin * Hm(ser,i);
                %I_intra =  sum(params.Pmin * Hm(ser,intra_interferers));
                distance_to_interferers = pdist2(locations.BS(ser,:),locations.M2M(inter_interferers,:),'euclidean') ;
                I_inter = sum( server_m2m_P(inter_interferers) .* Hm(ser,inter_interferers).* exp(-params.SEPL.alpha .* distance_to_interferers .^ params.SEPL.beta) );
                %I_m(i) = I_intra + I_inter;
                I_m(i) =  I_inter;  % special case #1 on intra interference
            end
            SINR_m = S_m./ (I_m + params.No);
            SIR_m_dB = 10*log10(SIR_m);
            PoCoverage(m,t) = sum(SIR_m_dB > params.Threshold) / N_users_M2M;
        end
        
    end
    
    normfact = params.space_realizations * params.time_slots ;% * simulation_area;
    Pcov(p) = sum(sum(PoCoverage)) / normfact
    %     fprintf(resultfile,'%s\n',char(datetime));
    %     fprintf(resultfile,'%s\n',num2str(params.LA_M));
    %     fprintf(resultfile,'%s\n',num2str(params.LA_B));
    %     fprintf(resultfile,'%s\n',num2str(Pcov));
end

end
%/////////////////////////////////////////////////////////////////////////
function P_exact = compute_uplink_coverage_with_max_transmit_power(params)
SNR = params.Pmin  / params.No;
noise_term = exp(- params.Threshold / SNR);
a = params.SEPL.alpha;
b = params.SEPL.beta;
Pu= params.Pm;
Po = params.Pmin;
n = (2/b) - 1;
t = params.Threshold;

theta = log((Pu./Po).^(1/a)).^(1/b);
F = @(y,v,tau) log(tau ./ y).^(v)  ./ ( y + 1);
% compute E_P(ln^k(P))
%A = (1 -  exp(- pi .* params.LA_B .* theta^2) ).^-1;
A = (n+1) * (pi / a^(n+1));
B = params.rho_m * (params.LA_M / params.N_RB);
for k = 0:n
    Ep(k+1,:) =  a^k .* gamma (k/(n+1) + 1) .* gammainc( pi .* params.LA_B .* theta.^2 , k/(n+1) + 1 , 'lower') ./ ((pi .* params.LA_B) .^(k/(n+1)) .* (1 - exp(- pi .* params.LA_B .* theta.^2)));
    NK(k+1) = nchoosek(n,k);
    J_v(k+1) = integral(@(y)F(y,n - k,t),0,t);
    T(k+1,:) = NK(k+1) *  J_v(k+1) * Ep(k+1,:) ;
end

S = sum(T,1);
interference_term = exp(- A * B * S);
P_exact = noise_term .* interference_term;

end

function Pcov = simulate_uplink_coverage_with_max_transmit_power(params)
simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.Pm);
PoCoverage = zeros(points,  params.space_realizations , params.time_slots );


for p = 1:points
    
    fprintf('\n')
    disp(['Max. Tranmit Power: ' , num2str(params.Pm(p))]);
    disp(['Coverage Threshold: ' , num2str(params.Threshold)]);
    disp(['Smallcell Density: ' , num2str(params.LA_B)]);
    disp(['M2M Density: ' , num2str(params.LA_M)]);
    disp(['Path Loss Parameters: ' ,'alpha=', num2str(params.SEPL.alpha) ,' beta=', num2str(params.SEPL.beta)]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_b = params.LA_B * simulation_area;
        %rho_m = betarnd(params.Beta_Distribution(1),params.Beta_Distribution(2));
        mu_m = params.LA_M * simulation_area;
        
        N_cells = poissrnd(mu_b);
        while(N_cells < 1)
            N_cells = poissrnd(mu_b);
        end
        N_users_M2M = round(params.rho_m * poissrnd(mu_m));
        
        locations.BS = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        locations.M2M = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users_M2M, 2);
        
        % Association
        distances_M2M_to_BS = pdist2(locations.BS,locations.M2M,'euclidean') ;
        [distance_to_server , server_m2m] = min(distances_M2M_to_BS);
        
        server_m2m(2,:)= randi(params.N_RB,1,N_users_M2M);   % each M2M node selects a carrier from N carriers randomly (uniform distribution)
        server_m2m_P = params.Pmin * exp(params.SEPL.alpha .* distance_to_server .^ params.SEPL.beta);
        server_m2m_P(server_m2m_P > params.Pm(p)) = 0;
        %server_m2m(3,:) = server_m2m_P;
        
        S_m = zeros(N_users_M2M,1);
        I_m = zeros(N_users_M2M,1);
        SINR_m = zeros(N_users_M2M,1);
        
        for t = 1:params.time_slots
            Hm = exprnd(1,N_cells,N_users_M2M);
            for i = 1:N_users_M2M
                ser = server_m2m(1,i);
                rb = server_m2m(2,i);
                %same_ser =  (server_m2m(1,:) == ser);
                another_ser =  (server_m2m(1,:) ~= ser);
                same_rb = (server_m2m(2,:) == rb);
                %intra_interferers = find(same_ser + same_rb == 2);
                inter_interferers = find(another_ser + same_rb == 2);
                S_m(i) = params.Pmin * Hm(ser,i);
                %I_intra =  sum(params.Pmin * Hm(ser,intra_interferers));
                distance_to_interferers = pdist2(locations.BS(ser,:),locations.M2M(inter_interferers,:),'euclidean') ;
                I_inter = sum( server_m2m_P(inter_interferers) .* Hm(ser,inter_interferers).* exp(-params.SEPL.alpha .* distance_to_interferers .^ params.SEPL.beta) );
                %I_m(i) = I_intra + I_inter;
                I_m(i) =  I_inter;  % special case #1 no intra interference
            end
            SINR_m = S_m./ (I_m + params.No);
            %SINR_m_dB = 10*log10(SINR_m);
            
            PoCoverage(p,m,t) = sum(SINR_m > params.Threshold) / N_users_M2M;
            
        end
        
    end
end
normfact = params.space_realizations * params.time_slots ;% * simulation_area;
Pcov = sum(sum(PoCoverage,3),2) / normfact;
end
%/////////////////////////////////////////////////////////////////////////
function P_exact = compute_uplink_coverage_with_resource_blocks(params)
SNR = params.Pmin  / params.No;
noise_term = exp(- params.Threshold ./ SNR);
a = params.SEPL.alpha;
b = params.SEPL.beta;
Pu= params.Pm;
Po = params.Pmin;
n = (2/b) - 1;
t = params.Threshold;

theta = log((Pu./Po).^(1/a)).^(1/b);
F = @(y,v,tau) log(tau ./ y).^(v)  ./ ( y + 1);
% compute E_P(ln^k(P))
%A = (1 -  exp(- pi .* params.LA_B .* theta^2) ).^-1;
A = (n+1) * (pi / a^(n+1));
B = params.rho_m * (params.LA_M ./ params.N_RB);
for k = 0:n
    Ep(k+1,:) =  a^k .* gamma (k/(n+1) + 1) .* gammainc( pi .* params.LA_B .* theta.^2 , k/(n+1) + 1 , 'lower') ./ ((pi .* params.LA_B) .^(k/(n+1)) .* (1 - exp(- pi .* params.LA_B .* theta.^2)));
    NK(k+1) = nchoosek(n,k);
    J_v(k+1) = integral(@(y)F(y,n - k,t),0,t);
    T(k+1,:) = NK(k+1) *  J_v(k+1) * Ep(k+1,:) ;
end

S = sum(T,1);
interference_term = exp(- A * B * S);
P_exact = noise_term .* interference_term;

end

function Pcov = simulate_uplink_coverage_with_resource_blocks(params)
simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.N_RB);
PoCoverage = zeros(points,  params.space_realizations , params.time_slots );


for p = 1:points
    
    fprintf('\n')
    disp(['Number of Resource Blocks: ' , num2str(params.N_RB(p))]);
    disp(['Coverage Threshold: ' , num2str(params.Threshold)]);
    disp(['Smallcell Density: ' , num2str(params.LA_B)]);
    disp(['M2M Density: ' , num2str(params.LA_M)]);
    disp(['Path Loss Parameters: ' ,'alpha=', num2str(params.SEPL.alpha) ,' beta=', num2str(params.SEPL.beta)]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_b = params.LA_B * simulation_area;
        %rho_m = betarnd(params.Beta_Distribution(1),params.Beta_Distribution(2));
        mu_m = params.LA_M * simulation_area;
        
        N_cells = poissrnd(mu_b);
        while(N_cells < 1)
            N_cells = poissrnd(mu_b);
        end
        N_users_M2M = round(params.rho_m * poissrnd(mu_m));
        
        locations.BS = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        locations.M2M = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users_M2M, 2);
        
        % Association
        distances_M2M_to_BS = pdist2(locations.BS,locations.M2M,'euclidean') ;
        [distance_to_server , server_m2m] = min(distances_M2M_to_BS);
        
        server_m2m(2,:)= randi(params.N_RB(p),1,N_users_M2M);   % each M2M node selects a carrier from N carriers randomly (uniform distribution)
        server_m2m_P = params.Pmin * exp(params.SEPL.alpha .* distance_to_server .^ params.SEPL.beta);
        server_m2m_P(server_m2m_P > params.Pm) = 0;
        %server_m2m(3,:) = server_m2m_P;
        
        S_m = zeros(N_users_M2M,1);
        I_m = zeros(N_users_M2M,1);
        SINR_m = zeros(N_users_M2M,1);
        
        for t = 1:params.time_slots
            Hm = exprnd(1,N_cells,N_users_M2M);
            for i = 1:N_users_M2M
                ser = server_m2m(1,i);
                rb = server_m2m(2,i);
                %same_ser =  (server_m2m(1,:) == ser);
                another_ser =  (server_m2m(1,:) ~= ser);
                same_rb = (server_m2m(2,:) == rb);
                %intra_interferers = find(same_ser + same_rb == 2);
                inter_interferers = find(another_ser + same_rb == 2);
                S_m(i) = params.Pmin * Hm(ser,i);
                %I_intra =  sum(params.Pmin * Hm(ser,intra_interferers));
                distance_to_interferers = pdist2(locations.BS(ser,:),locations.M2M(inter_interferers,:),'euclidean') ;
                I_inter = sum( server_m2m_P(inter_interferers) .* Hm(ser,inter_interferers).* exp(-params.SEPL.alpha .* distance_to_interferers .^ params.SEPL.beta) );
                %I_m(i) = I_intra + I_inter;
                I_m(i) =  I_inter;  % special case #1 no intra interference
            end
            SINR_m = S_m./ (I_m + params.No);
            %SINR_m_dB = 10*log10(SINR_m);
            
            PoCoverage(p,m,t) = sum(SINR_m > params.Threshold) / N_users_M2M;
            
        end
        
    end
end
normfact = params.space_realizations * params.time_slots ;% * simulation_area;
Pcov = sum(sum(PoCoverage,3),2) / normfact;
end
%/////////////////////////////////////////////////////////////////////////
function P_exact = compute_uplink_coverage_with_power_truncation_threshold(params)
SNR = params.Pmin  / params.No;
noise_term = exp(- params.Threshold ./ SNR);
a = params.SEPL.alpha;
b = params.SEPL.beta;
Pu= params.Pm;
Po = params.Pmin;
n = (2/b) - 1;
t = params.Threshold;

theta = log((Pu./Po).^(1/a)).^(1/b);
F = @(y,v,tau) log(tau ./ y).^(v)  ./ ( y + 1);
% compute E_P(ln^k(P))
%A = (1 -  exp(- pi .* params.LA_B .* theta^2) ).^-1;
A = (n+1) * (pi / a^(n+1));
B = params.rho_m * (params.LA_M / params.N_RB);
for k = 0:n
    Ep(k+1,:) =  a^k .* gamma (k/(n+1) + 1) .* gammainc( pi .* params.LA_B .* theta.^2 , k/(n+1) + 1 , 'lower') ./ ((pi .* params.LA_B) .^(k/(n+1)) .* (1 - exp(- pi .* params.LA_B .* theta.^2)));
    NK(k+1) = nchoosek(n,k);
    J_v(k+1) = integral(@(y)F(y,n - k,t),0,t);
    T(k+1,:) = NK(k+1) *  J_v(k+1) * Ep(k+1,:) ;
end

S = sum(T,1);
interference_term = exp(- A * B * S);
P_exact = noise_term .* interference_term;

end

function Pcov = simulate_uplink_coverage_with_power_truncation_threshold(params)
simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.Pmin);
PoCoverage = zeros(points,  params.space_realizations , params.time_slots );


for p = 1:points
    
    fprintf('\n')
    disp(['Power Truncation Threshold: ' , num2str(params.Pmin(p))]);
    disp(['Coverage Threshold: ' , num2str(params.Threshold)]);
    disp(['Smallcell Density: ' , num2str(params.LA_B)]);
    disp(['M2M Density: ' , num2str(params.LA_M)]);
    disp(['Path Loss Parameters: ' ,'alpha=', num2str(params.SEPL.alpha) ,' beta=', num2str(params.SEPL.beta)]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_b = params.LA_B * simulation_area;
        %rho_m = betarnd(params.Beta_Distribution(1),params.Beta_Distribution(2));
        mu_m = params.LA_M * simulation_area;
        
        N_cells = poissrnd(mu_b);
        while(N_cells < 1)
            N_cells = poissrnd(mu_b);
        end
        N_users_M2M = round(params.rho_m * poissrnd(mu_m));
        
        locations.BS = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        locations.M2M = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users_M2M, 2);
        
        % Association
        distances_M2M_to_BS = pdist2(locations.BS,locations.M2M,'euclidean') ;
        [distance_to_server , server_m2m] = min(distances_M2M_to_BS);
        
        server_m2m(2,:)= randi(params.N_RB,1,N_users_M2M);   % each M2M node selects a carrier from N carriers randomly (uniform distribution)
        server_m2m_P = params.Pmin(p) * exp(params.SEPL.alpha .* distance_to_server .^ params.SEPL.beta);
        server_m2m_P(server_m2m_P > params.Pm) = 0;
        %server_m2m(3,:) = server_m2m_P;
        
        S_m = zeros(N_users_M2M,1);
        I_m = zeros(N_users_M2M,1);
        SINR_m = zeros(N_users_M2M,1);
        
        for t = 1:params.time_slots
            Hm = exprnd(1,N_cells,N_users_M2M);
            for i = 1:N_users_M2M
                ser = server_m2m(1,i);
                rb = server_m2m(2,i);
                %same_ser =  (server_m2m(1,:) == ser);
                another_ser =  (server_m2m(1,:) ~= ser);
                same_rb = (server_m2m(2,:) == rb);
                %intra_interferers = find(same_ser + same_rb == 2);
                inter_interferers = find(another_ser + same_rb == 2);
                S_m(i) = params.Pmin(p) * Hm(ser,i);
                %I_intra =  sum(params.Pmin * Hm(ser,intra_interferers));
                distance_to_interferers = pdist2(locations.BS(ser,:),locations.M2M(inter_interferers,:),'euclidean') ;
                I_inter = sum( server_m2m_P(inter_interferers) .* Hm(ser,inter_interferers).* exp(-params.SEPL.alpha .* distance_to_interferers .^ params.SEPL.beta) );
                %I_m(i) = I_intra + I_inter;
                I_m(i) =  I_inter;  % special case #1 no intra interference
            end
            SINR_m = S_m./ (I_m + params.No);
            %SINR_m_dB = 10*log10(SINR_m);
            
            PoCoverage(p,m,t) = sum(SINR_m > params.Threshold) / N_users_M2M;
            
        end
        
    end
end
normfact = params.space_realizations * params.time_slots ;% * simulation_area;
Pcov = sum(sum(PoCoverage,3),2) / normfact;
end
%/////////////////////////////////////////////////////////////////////////
function [cell_loads h2h_loads ch_loads] = compute_average_load_with_small_cell_density(params)
simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.LA_B);
AVERAGE_LOAD= zeros(points,  params.space_realizations);
AVERAGE_LOAD_AGG_H2H= zeros(points,  params.space_realizations); % load of H2H and CH aggregators
AVERAGE_LOAD_AGG_CH= zeros(points,  params.space_realizations); % load of H2H and CH aggregators
for p = 1:points
    fprintf('\n')
    disp(['Smallcell Density: ' , num2str(params.LA_B(p))]);
    disp(['M2M Density: ' , num2str(params.LA_M)]);
    %disp(['Activiation Probability rho_m: ' , num2str(rho_m)]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_b = params.LA_B(p) * simulation_area;
        rho_m = betarnd(params.Beta_Distribution(1),params.Beta_Distribution(2));
        mu_h = params.LA_H * simulation_area;
        mu_m = params.LA_M * simulation_area;
        N_cells = poissrnd(mu_b);
        while(N_cells < 1)
            N_cells = poissrnd(mu_b);
        end
        N_users_H2H = poissrnd(mu_h);
        N_users_M2M = round(rho_m * poissrnd(mu_m));
        
        
        locations.BS = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        locations.H2H = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users_H2H, 2);
        locations.M2M = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users_M2M, 2);
        if(strcmp(params.aggregation_mode , 'AGGREGATION'))
            mu_c = params.LA_C * simulation_area;
            N_nodes_CH = round(mu_c);    % Fixed Privilged Cluster Heads
            locations.CH= params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_nodes_CH, 2);
        end
        [~ , LOADS] = associate_nodes(locations,params.aggregation_mode);
        
        if(strcmp(params.aggregation_mode , 'DIRECT'))
            AVERAGE_LOAD(p,m) = sum(LOADS.BS) / sum(LOADS.BS ~= 0);
        else
            AVERAGE_LOAD(p,m) = mean(LOADS.UE_H2H) + mean(LOADS.UE_CH);
            AVERAGE_LOAD_AGG_H2H(p,m) =  mean(LOADS.H2H);
            AVERAGE_LOAD_AGG_CH(p,m) =  mean(LOADS.CH);
            save('results','LOADS')
        end
    end
end
cell_loads  = ceil(mean(AVERAGE_LOAD,2));
h2h_loads   = ceil(mean(AVERAGE_LOAD_AGG_H2H,2));
ch_loads   = ceil(mean(AVERAGE_LOAD_AGG_CH,2));
end
%/////////////////////////////////////////////////////////////////////////
function [association loads] = associate_nodes(locations,mode)
switch(mode)
    case 'AGGREGATION'
        distances_H2H_to_BS = pdist2(locations.BS,locations.H2H,'euclidean') ;
        [~ , server_h2h] = min(distances_H2H_to_BS);
        association.H2H_to_BS = server_h2h;
        
        distances_CH_to_BS = pdist2(locations.BS,locations.CH,'euclidean') ;
        [~ , server_ch] = min(distances_CH_to_BS);
        association.CH_to_BS = server_ch;
        
        aggregators = [locations.H2H ; locations.CH];
        N_H2H = size(locations.H2H ,1);
        N_CH = size(locations.CH ,1);
        distances_M2M_to_H2H = pdist2(aggregators,locations.M2M,'euclidean') ;
        [~ , server_m2m] = min(distances_M2M_to_H2H);
        
        Im2m_h2h = find(server_m2m <= N_H2H);
        Im2m_ch = find(server_m2m > N_H2H);
        Ih2h = server_m2m(Im2m_h2h);
        Ich = server_m2m(Im2m_ch) - N_H2H  ;
        
        association.M2M_to_H2H = [Im2m_h2h' Ih2h'];
        association.M2M_to_CH = [Im2m_ch' Ich'];
        
        for i = 1:N_H2H
            loads.H2H(i) = sum(association.M2M_to_H2H(:,2) == i);
        end
        
        for i = 1:N_CH
            loads.CH(i) = sum(association.M2M_to_CH(:,2) == i);
        end
        
        for i = 1:size(locations.BS ,1)
            loads.UE_H2H(i) = sum(association.H2H_to_BS == i);
            loads.UE_CH(i) = sum(association.CH_to_BS == i);
        end
        
    case 'DIRECT'
        
        distances_H2H_to_BS = pdist2(locations.BS,locations.H2H,'euclidean') ;
        [~ , server_h2h] = min(distances_H2H_to_BS);
        association.H2H_to_BS = server_h2h;
        
        distances_M2M_to_BS = pdist2(locations.BS,locations.M2M,'euclidean') ;
        [~ , server_m2m] = min(distances_M2M_to_BS);
        association.M2M_to_BS = server_m2m;
        
        N_BS = size(locations.BS ,1);
        for i = 1:N_BS
            loads.BS(i) = sum(association.M2M_to_BS == i);
        end
end
end
%/////////////////////////////////////////////////////////////////////////
function visualize_the_network(params)
simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;

mu_b = params.LA_B * simulation_area;
rho_m = betarnd(params.Beta_Distribution(1),params.Beta_Distribution(2));
mu_h = params.LA_H * simulation_area;
mu_m = params.LA_M * simulation_area;
mu_c = params.LA_C * simulation_area;
N_cells = poissrnd(mu_b);
while(N_cells < 1)
    N_cells = poissrnd(mu_b);
end
N_users_H2H = poissrnd(mu_h);
N_users_M2M = round(rho_m * poissrnd(mu_m));
N_nodes_CH = round(mu_c);    % Fixed Privilged Cluster Heads

locations.BS = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
locations.H2H = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users_H2H, 2);
locations.M2M = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users_M2M, 2);
locations.CH= params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_nodes_CH, 2);

[association,loads] = associate_nodes(locations,params.aggregation_mode);

%plot (0,0,'X','MarkerSize',20,'LineWidth',4,'Color','black')
hold on;

M2M = plot_nodes(locations.M2M(:,1),locations.M2M(:,2),'.r',30);
H2H = plot_nodes(locations.H2H(:,1),locations.H2H(:,2),'sg',15);
cells = plot_nodes(locations.BS(:,1),locations.BS(:,2),'^b',15);
annotate_node_id(locations.BS,15)
annotate_node_id(locations.H2H,15)
annotate_node_id(locations.M2M,10)
if(strcmp(params.aggregation_mode,'AGGREGATION') )
    CH = plot_nodes(locations.CH(:,1),locations.CH(:,2),'*m',15);
    annotate_node_id(locations.CH,15)
    legend ({'M2M'  , 'H2H' ,'BS','CH' })
else
    legend ({'M2M'  , 'H2H' ,'BS' })
end
connect_nodes_to_servers(locations,association,params.aggregation_mode)

axis equal;
X_lim = params.simulation_area_side;
Y_lim = params.simulation_area_side;
axis([X_lim Y_lim]);
xlabel('x (meters)');
ylabel('y (meters)');
set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 20);
set(gca, 'FontWeight', 'Bold');
grid on;

end
%/////////////////////////////////////////////////////////////////////////
function y = plot_nodes(x,y,marker,size)
v = plot(x,y,marker);
set(v,'MarkerSize',size);
set(v,'LineWidth',4);
end
%/////////////////////////////////////////////////////////////////////////
function annotate_node_id(nodes,fontsize)
for i = 1:size(nodes,1)
    text(nodes(i,1) + 1,nodes(i,2) + 1,num2str(i),'Color','black','FontSize',fontsize);
end
end
%/////////////////////////////////////////////////////////////////////////
function connect_nodes_to_servers(locations,association,mode)
switch(mode)
    case 'AGGREGATION'
        for i = 1:size(association.M2M_to_H2H,1);
            x = [locations.M2M(association.M2M_to_H2H(i,1),1) , locations.H2H(association.M2M_to_H2H(i,2),1) ];
            y = [locations.M2M(association.M2M_to_H2H(i,1),2) , locations.H2H(association.M2M_to_H2H(i,2),2) ];
            line(x ,y,'Color','green','LineStyle','--','LineWidth',1)
        end
        
        for i = 1:size(association.M2M_to_CH,1);
            x = [locations.M2M(association.M2M_to_CH(i,1),1) , locations.CH(association.M2M_to_CH(i,2),1) ];
            y = [locations.M2M(association.M2M_to_CH(i,1),2) , locations.CH(association.M2M_to_CH(i,2),2) ];
            line(x ,y,'Color','magenta','LineStyle','--','LineWidth',1)
        end
        
        for i = 1:size(association.H2H_to_BS,2);
            x = [locations.H2H(i,1) , locations.BS(association.H2H_to_BS(i),1) ];
            y = [locations.H2H(i,2) , locations.BS(association.H2H_to_BS(i),2) ];
            line(x ,y,'Color','green','LineStyle','-','LineWidth',3)
        end
        
        for i = 1:size(association.CH_to_BS,2);
            x = [locations.CH(i,1) , locations.BS(association.CH_to_BS(i),1) ];
            y = [locations.CH(i,2) , locations.BS(association.CH_to_BS(i),2) ];
            line(x ,y,'Color','magenta','LineStyle','-','LineWidth',3)
        end
    case 'DIRECT'
        for i = 1:size(association.M2M_to_BS,2);
            x = [locations.M2M(i,1) , locations.BS(association.M2M_to_BS(i),1) ];
            y = [locations.M2M(i,2) , locations.BS(association.M2M_to_BS(i),2) ];
            line(x ,y,'Color','red','LineStyle','--','LineWidth',1)
        end
        
        for i = 1:size(association.H2H_to_BS,2);
            x = [locations.H2H(i,1) , locations.BS(association.H2H_to_BS(i),1) ];
            y = [locations.H2H(i,2) , locations.BS(association.H2H_to_BS(i),2) ];
            line(x ,y,'Color','green','LineStyle','-','LineWidth',3)
        end
end
end
