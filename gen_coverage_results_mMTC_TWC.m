function [ output_args ] = gen_coverage_results_mMTC(sim)
%GEN_COVERAGE_RESULTS_MMTC-D  Summary of this function goes here
%   Author : Mahmoud Kamel
%   Date :  12-Dec-2017 12:01:39
%   This function generates results of simulating the effect of clustering
%   on the load of cells and coverage
%/////////////////////////////////////////////////////////////////////////
%close all;
hold on;
% Simulation Parameters
Pm_dB = 20;
params.Pm = 10.^(Pm_dB/10) * 1e-3;                         % Uplink Transmit power of M2M nodes in watts 20 dBm
Pmin_dBm = -100 ;
params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;                         % Min. received power at the base station in watts > BS sensitivity
params.simulation_area_side = [-250 250];    % simulation area side to side
params.space_realizations = 100;
params.time_slots = 1;
params.BW = 20e6;
params.N_RB = 100;
%params.N = 10;    % Orthogonal carriers allocated to the RACH transmission of M2M traffic
params.RB_BW = 180e3; %Hz
params.No = 10^(-17.4) * params.RB_BW * 1e-3;
params.SEPL.alpha = 0.94;
params.SEPL.beta = 1/2;
params.LA_M = 0.1;
params.rho_m = 0.1;
params.sim_style = {'b^','g*' ,'rx'};
params.ana_style = {'b:','g--','r-'};
params.ana_only_style = {'b:^','g*--','rx-'};
%/////////////////////////////////////////////////////////////////////////
% Simulation Switches


%sim = 'PC_with_SMALL_SCALE_FADING';

%sim = 'General_Beta';

%sim = 'AveragedCoverage_tau_Pmin';
%sim = 'AveragedCoverage_tau_Pm';
%sim = 'AveragedCoverage_tau_N_RB';
%sim = 'AveragedCoverage_tau_SEPL';


%sim = 'MaxTransmitPower_SEPL';
%sim = 'MaxTransmitPower_Pmin';
%sim = 'MaxTransmitPower_N_RB';


%sim = 'CoverageThreshold_N_RB';
%sim = 'CoverageThreshold_Pmin';
%sim = 'CoverageThreshold_SEPL';
%sim = 'CoverageThreshold_Pm';

%sim = 'CellDensity_Pmin';
%sim = 'CellDensity_SEPL';
%sim = 'CellDensity_Pm';
%sim = 'CellDensity_N_RB';
%sim = 'CellDensity_Scaling_N_RB';


%sim = 'ResourceBlocks_Pmin';
%sim = 'ResourceBlocks_SEPL';
%sim = 'ResourceBlocks_Pm';

%sim = 'PowerTruncationThreshold_SEPL';
%sim = 'PowerTruncationThreshold_Pm';
%sim = 'PowerTruncationThreshold_N_RB';
%sim = 'PowerTruncationThreshold_Cell_Density';

%sim = 'M2MDensity_Pmin';
%sim = 'M2MDensity_Pm';
%sim = 'M2MDensity_SEPL';
%sim = 'M2MDensity_N_RB';



%sim = 'Log_P_Po_Pm';
%sim = 'Log_P_Po_SEPL';

%sim = 'Truncation_Outage_Pmin';
%sim = 'Truncation_Outage_SEPL';
%sim = 'Truncation_Outage_CellDensity_SEPL';
%sim = 'Truncation_Outage_Pmax_SEPL';



%sim = 'Validate_PDF_OF_P';

%sim = 'Expectation_Log_P';
%sim = 'Expectation_Log__K_P';
%sim = 'CellDensity';

%sim = 'ClusterHeadDensity';

%sim = 'Closed_Form';


switch(sim)
    
    case 'Closed_Form'
        params.simulation_area_side = [-100 100];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 10;
        params.LA_B = 1000e-6   ;       % Small Cell Denisty (cells/m^2)
        Threshold_dB = (-20:1:20) ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        LA_M = [1e-2 5e-2 1e-1];
        params.rho_m = 0.1;
        Pmin_dBm = -100 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        params.SEPL.alpha = 5e-2;
        params.SEPL.beta = 2;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            [Pcov_analy(i,:) Pcov_analy_cf(i,:)] = compute_uplink_coverage_with_coverage_threshold_closed_form(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_coverage_threshold(params);
        end
        figure;
        subplot(1,2,1);
        hold on;
        
        for plt = 1:numel(LA_M)
            plot(Threshold_dB  ,Pcov_analy_cf(plt,:),'r.',Threshold_dB  ,Pcov_analy(plt,:), params.ana_style{plt},Threshold_dB,Pcov_simul(plt,:), params.sim_style{plt});
        end
        
        h = get(gca,'Children') ;
        deco.xlabel = 'Coverage threshold ($\tau$ dB)  ';
        deco.ylabel = 'Uplink coverage probability ';
        deco.title = '$P_o = -100$ dBm';
        deco.legend.items = 1:2*numel(LA_M) ;
        ind = 0;
        for i = 1:2:2*numel(LA_M)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);
        %####################################################
        Pmin_dBm = -110 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            [Pcov_analy(i,:) Pcov_analy_cf(i,:)] = compute_uplink_coverage_with_coverage_threshold_closed_form(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_coverage_threshold(params);
        end
        subplot(1,2,2);
        hold on;
        
        for plt = 1:numel(LA_M)
            plot(Threshold_dB  ,Pcov_analy_cf(plt,:),'r.',Threshold_dB  ,Pcov_analy(plt,:), params.ana_style{plt},Threshold_dB,Pcov_simul(plt,:), params.sim_style{plt});
        end
        
        h = get(gca,'Children') ;
        deco.title = '$P_o = -110$ dBm';
        deco.legend.items = 1:2*numel(LA_M) ;
        ind = 0;
        for i = 1:2:2*numel(LA_M)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);
        
        savefig('TWCResults/coverage_tau_Po_cf');
        
        
    case 'PC_with_SMALL_SCALE_FADING';
        params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 10;
        params.time_slots = 10;
        params.LA_B = 1000e-6   ;       % Small Cell Denisty (cells/m^2)
        Threshold_dB = (-20:1:20) ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        LA_M = [1e-2 5e-2 1e-1];
        params.rho_m = 0.1;
        Pmin_dBm = -100 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_coverage_threshold_PC_SmalelScale(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_coverage_threshold_PC_SmalelScale(params);
        end
        figure;
        subplot(1,2,1);
        hold on;
        for plt = 1:numel(LA_M)
            plot(Threshold_dB  ,Pcov_analy_cf(plt,:),'r.',Threshold_dB  ,Pcov_analy(plt,:), params.ana_style{plt},Threshold_dB,Pcov_simul(plt,:), params.sim_style{plt});
        end
        
        h = get(gca,'Children') ;
        deco.xlabel = 'Coverage threshold ($\tau$ dB)  ';
        deco.ylabel = 'Uplink coverage probability';
        deco.title = '$\alpha = 0.94, \beta = \frac{1}{2}$';
        deco.legend.items = 1:2*numel(LA_M) ;
        ind = 0;
        for i = 1:2:2*numel(LA_M)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);
        
        params.SEPL.alpha = 0.3;
        params.SEPL.beta = 2/3;
        
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:)= compute_uplink_coverage_with_coverage_threshold(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_coverage_threshold(params);
        end
        subplot(1,2,2);
        hold on;
        
        for plt = 1:numel(LA_M)
            plot(Threshold_dB  ,Pcov_analy(plt,:), params.ana_style{plt},Threshold_dB,Pcov_simul(plt,:), params.sim_style{plt});
        end
        
        h = get(gca,'Children') ;
        deco.title = '$\alpha = 0.3, \beta = \frac{2}{3}$';
        deco.legend.items = 1:2*numel(LA_M) ;
        ind = 0;
        for i = 1:2:2*numel(LA_M)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);
        
        savefig('TWCResults/coverage_tau_SEPL_PC');
    case 'General_Beta'
        params.Nterms = 7;
        params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 10;
        params.LA_B = 1000e-6   ;       % Small Cell Denisty (cells/m^2)
        params.Beta_Distribution = [3 4];   % Beta Distribution parameters                                   % of M2M nodes traffic activity
        Threshold_dB = (-20:1:20) ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        LA_M = [1e-2 5e-2 1e-1];
        params.rho_m = 0.1;
        Pmin_dBm = -70 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,: )= compute_uplink_coverage_with_coverage_threshold_general(params);
            %Pcov_simul(i,:) = simulate_uplink_coverage_with_coverage_threshold_general(params);
        end
        hold on;
        %figure;
        %subplot(1,2,1);
        %f1 = plot(Threshold_dB  ,Pcov_analy , 'k-',Threshold_dB,Pcov_simul,'ro' );
        f1 = plot(Threshold_dB  ,Pcov_analy , 'b-' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        %legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Coverage threshold ($\tau$)  ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$P_{min} = -70$ dBm  ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
    case 'MaxTransmitPower_SEPL'
        disp('MaxTransmitPower_SEPL');
        params.simulation_area_side = [-100 100];    % simulation area side to side
        params.space_realizations = 10;
        params.time_slots = 10;
        params.rho_m = 0.1;                                    % of M2M nodes traffic activity
        params.aggregation_mode = 'DIRECT';
        Threshold_dB = 0 ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        params.rho_m = 0.1;
        params.LA_B = 1000*1e-6 ;       % Small Cell Denisty (cells/m^2)
        Pmin_dBm = -100 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        params.N_RB = 100;
        
        Pm_dBm = 0:5:20; % in dBm
        params.Pm = 10.^(Pm_dBm/10) * 1e-3;
        
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_max_transmit_power(params);
        end
        figure
        subplot(1,2,1);
        
        hold on;
        for plt = 1:numel(LA_M)
            plot(Pm_dBm  ,Pcov_analy(plt,:), params.ana_only_style{plt});
        end
        
        h = get(gca,'Children') ;
        deco.xlabel = 'Maximum transmit power ($P_{m}$ dBm)';
        deco.ylabel = 'Uplink coverage probability ';
        deco.title = '$\alpha = 0.94, \beta = \frac{1}{2}$';
        deco.legend.items = 1:2*numel(LA_M) ;
        for i = 1:numel(LA_M)
            deco.legend.labels{i} = strcat('$(\lambda_m =', num2str(LA_M(i)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);
        %#################################################################
        params.SEPL.alpha = 0.3;
        params.SEPL.beta = 2/3;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_max_transmit_power(params);
        end
        subplot(1,2,2);
        hold on;
        
        for plt = 1:numel(LA_M)
            plot(Pm_dBm  ,Pcov_analy(plt,:), params.ana_only_style{plt});
        end
        
        h = get(gca,'Children') ;
        deco.title = '$\alpha = 0.3, \beta = \frac{2}{3}$';
        
        decorate_plot(h,deco);
        
        savefig('TWCResults/coverage_Pm_SEPL');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
    case 'MaxTransmitPower_Pmin'
        disp('MaxTransmitPower_Pmin');
        params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 10;
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        params.rho_m = 0.1;                                    % of M2M nodes traffic activity
        params.aggregation_mode = 'DIRECT';
        Threshold_dB = 0 ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        params.rho_m = 0.1;
        params.LA_B = 1000*1e-6 ;       % Small Cell Denisty (cells/m^2)
        
        params.N_RB = 100;
        
        Pm_dBm = 0:2:33; % in dBm
        params.Pm = 10.^(Pm_dBm/10) * 1e-3;
        
        Pmin_dBm = -70 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_max_transmit_power(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_max_transmit_power(params);
        end
        figure
        subplot(1,2,1);
        f1 = plot(Pm_dBm ,Pcov_analy , 'k-',Pm_dBm,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel(' Maximum transmit power $P_{m}$(dBm) ' ,'Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$P_{min} = -70$ dBm  ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        set(gca, 'LineWidth', 2);
        set(gca, 'GridAlpha', 0.5);
        set(gca, 'MinorGridAlpha', 0.5);
        
        Pmin_dBm = -80 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_max_transmit_power(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_max_transmit_power(params);
        end
        subplot(1,2,2);
        f1 = plot(Pm_dBm,Pcov_analy , 'k-',Pm_dBm,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel(' Maximum transmit power $P_{m}$(dBm) ' ,'Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$P_{min} = -80$ dBm  ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        set(gca, 'LineWidth', 2);
        set(gca, 'GridAlpha', 0.5);
        set(gca, 'MinorGridAlpha', 0.5);
        
        savefig('FinalResults/coverage_Pm_Pmin');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
    case 'MaxTransmitPower_N_RB'
        disp('MaxTransmitPower_N_RB');
        params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 10;
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        params.rho_m = 0.1;                                    % of M2M nodes traffic activity
        params.aggregation_mode = 'DIRECT';
        Threshold_dB = 0 ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        params.rho_m = 0.1;
        params.LA_B = 1000*1e-6 ;       % Small Cell Denisty (cells/m^2)
        
        params.N_RB = 100;
        
        Pm_dBm = 0:2:33; % in dBm
        params.Pm = 10.^(Pm_dBm/10) * 1e-3;
        
        Pmin_dBm = -70 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        
        
        params.N_RB = 50;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_max_transmit_power(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_max_transmit_power(params);
        end
        figure
        subplot(1,2,1);
        f1 = plot(Pm_dBm ,Pcov_analy , 'k-',Pm_dBm,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel(' Maximum transmit power $P_{m}$(dBm) ' ,'Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$N_{RB} = 50$ ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        set(gca, 'LineWidth', 2);
        set(gca, 'GridAlpha', 0.5);
        set(gca, 'MinorGridAlpha', 0.5);
        
        params.N_RB = 100;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_max_transmit_power(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_max_transmit_power(params);
        end
        subplot(1,2,2);
        f1 = plot(Pm_dBm,Pcov_analy , 'k-',Pm_dBm,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel(' Maximum transmit power $P_{m}$(dBm) ' ,'Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$N_{RB} = 100$ ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        set(gca, 'LineWidth', 2);
        set(gca, 'GridAlpha', 0.5);
        set(gca, 'MinorGridAlpha', 0.5);
        
        savefig('FinalResults/coverage_Pm_N_RB');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
        
        
    case 'AveragedCoverage_tau_Pmin'
        params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 10;
        params.LA_B = 1000e-6   ;       % Small Cell Denisty (cells/m^2)
        params.Beta_Distribution = [3 4];   % Beta Distribution parameters
        % of M2M nodes traffic activity
        
        
        Threshold_dB = (-10:1:20) ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        
        LA_M = [1e-2 5e-2 1e-1];
        
        Pmin_dBm = -100 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov__analy(i,:) = compute_average_coverage_with_coverage_threshold(params);
            Pcov_simul(i,:) = simulate_average_coverage_with_coverage_threshold(params);
        end
        figure;
        subplot(1,2,1);
        params.sim_style = {'rx','g*' ,'b^'};
        params.ana_style = {'r-','g--','b:'};
        hold on;
        
        for plt = 1:numel(LA_B)
            plot(Threshold_dB  ,Pcov__analy(plt,:) , params.ana_style{plt},Threshold_dB,Pcov_simul(plt,:),params.sim_style{plt} );
        end
        
        h = get(gca,'Children') ;
        deco.xlabel = 'Coverage threshold ($\tau$)  ';
        deco.ylabel = 'Uplink coverage probability ';
        deco.title = 'HTC';
        deco.legend.items = 1:2*numel(LA_B) ;
        ind = 0;
        for i = 1:2:2*numel(LA_B)
            ind = ind + 1;
            deco.legend.labels{i+1} = strcat('Analysis $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
            deco.legend.labels{i} = strcat('Simulation $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
        end
        deco.legend.location = 'northeast';
        decorate_plot(h,deco);
        
        
        
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Coverage threshold ($\tau$)  ','Interpreter','LaTex');
        ylabel('Averaged uplink coverage probability ','Interpreter','LaTex');
        title('$P_{min} = -70$ dBm  ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        set(gca, 'LineWidth', 2);
        set(gca, 'GridAlpha', 0.5);
        set(gca, 'MinorGridAlpha', 0.5);
        
        Pmin_dBm = -110 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov__analy(i,:) = compute_average_coverage_with_coverage_threshold(params);
            Pcov_simul(i,:) = simulate_average_coverage_with_coverage_threshold(params);
        end
        subplot(1,2,2);
        f1 = plot(Threshold_dB  ,Pcov__analy , 'k-',Threshold_dB,Pcov_simul,'ro' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Coverage threshold ($\tau$)  ','Interpreter','LaTex');
        ylabel('Averaged uplink coverage probability ','Interpreter','LaTex');
        title('$P_{min} = -80$ dBm  ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        set(gca, 'LineWidth', 2);
        set(gca, 'GridAlpha', 0.5);
        set(gca, 'MinorGridAlpha', 0.5);
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'AveragedCoverage_tau_Pm'
        disp('AveragedCoverage_tau_Pm');
        params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 10;
        params.LA_B = 1000e-6   ;       % Small Cell Denisty (cells/m^2)
        params.Beta_Distribution = [3 4];   % Beta Distribution parameters
        % of M2M nodes traffic activity
        Pmin_dBm = -70 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        
        Threshold_dB = (-10:1:20) ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        
        LA_M = [1e-2 5e-2 1e-1];
        
        Pm_dB = 20;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov__analy(i,:) = compute_average_coverage_with_coverage_threshold(params);
            Pcov_simul(i,:) = simulate_average_coverage_with_coverage_threshold(params);
        end
        figure;
        subplot(1,2,1);
        f1 = plot(Threshold_dB  ,Pcov__analy , 'k-',Threshold_dB,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Coverage threshold ($\tau$)  ','Interpreter','LaTex');
        ylabel('Averaged uplink coverage probability ','Interpreter','LaTex');
        title('$P_{m} = 20$ dBm  ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        set(gca, 'LineWidth', 2);
        set(gca, 'GridAlpha', 0.5);
        set(gca, 'MinorGridAlpha', 0.5);
        
        Pm_dB = 33;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov__analy(i,:) = compute_average_coverage_with_coverage_threshold(params);
            Pcov_simul(i,:) = simulate_average_coverage_with_coverage_threshold(params);
        end
        subplot(1,2,2);
        f1 = plot(Threshold_dB  ,Pcov__analy , 'k-',Threshold_dB,Pcov_simul,'ro' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Coverage threshold ($\tau$)  ','Interpreter','LaTex');
        ylabel('Averaged uplink coverage probability ','Interpreter','LaTex');
        title('$P_{m} = 33$ dBm  ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        set(gca, 'LineWidth', 2);
        set(gca, 'GridAlpha', 0.5);
        set(gca, 'MinorGridAlpha', 0.5);
        savefig('FinalResults/average_coverage_tau_Pm');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'AveragedCoverage_tau_SEPL'
        disp('AveragedCoverage_tau_SEPL');
        params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 10;
        params.LA_B = 1000e-6   ;       % Small Cell Denisty (cells/m^2)
        params.Beta_Distribution = [3 4];   % Beta Distribution parameters
        % of M2M nodes traffic activity
        Pmin_dBm = -70 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        Pm_dB = 20;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        
        Threshold_dB = (-10:1:20) ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        
        LA_M = [1e-2 5e-2 1e-1];
        
        
        
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov__analy(i,:) = compute_average_coverage_with_coverage_threshold(params);
            Pcov_simul(i,:) = simulate_average_coverage_with_coverage_threshold(params);
        end
        figure;
        subplot(1,2,1);
        f1 = plot(Threshold_dB  ,Pcov__analy , 'k-',Threshold_dB,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Coverage threshold ($\tau$)  ','Interpreter','LaTex');
        ylabel('Averaged uplink coverage probability ','Interpreter','LaTex');
        title('$\alpha = 0.94, \beta = \frac{1}{2}$','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        set(gca, 'LineWidth', 2);
        set(gca, 'GridAlpha', 0.5);
        set(gca, 'MinorGridAlpha', 0.5);
        
        params.SEPL.alpha = 0.3;
        params.SEPL.beta = 2/3;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov__analy(i,:) = compute_average_coverage_with_coverage_threshold(params);
            Pcov_simul(i,:) = simulate_average_coverage_with_coverage_threshold(params);
        end
        subplot(1,2,2);
        f1 = plot(Threshold_dB  ,Pcov__analy , 'k-',Threshold_dB,Pcov_simul,'ro' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Coverage threshold ($\tau$)  ','Interpreter','LaTex');
        ylabel('Averaged uplink coverage probability ','Interpreter','LaTex');
        title('$\alpha = 0.94, \beta = \frac{1}{2}$','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        set(gca, 'LineWidth', 2);
        set(gca, 'GridAlpha', 0.5);
        set(gca, 'MinorGridAlpha', 0.5);
        savefig('FinalResults/average_coverage_tau_SEPL');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'AveragedCoverage_tau_N_RB'
        disp('AveragedCoverage_tau_N_RB');
        params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 10;
        params.LA_B = 1000e-6   ;       % Small Cell Denisty (cells/m^2)
        params.Beta_Distribution = [3 4];   % Beta Distribution parameters
        % of M2M nodes traffic activity
        Pmin_dBm = -70 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        Pm_dB = 20;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        
        Threshold_dB = (-10:1:20) ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        
        LA_M = [1e-2 5e-2 1e-1];
        
        
        
        params.N_RB = 50;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov__analy(i,:) = compute_average_coverage_with_coverage_threshold(params);
            Pcov_simul(i,:) = simulate_average_coverage_with_coverage_threshold(params);
        end
        figure;
        subplot(1,2,1);
        f1 = plot(Threshold_dB  ,Pcov__analy , 'k-',Threshold_dB,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Coverage threshold ($\tau$)  ','Interpreter','LaTex');
        ylabel('Averaged uplink coverage probability ','Interpreter','LaTex');
        title('$N_{RB} = 50$ ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        set(gca, 'LineWidth', 2);
        set(gca, 'GridAlpha', 0.5);
        set(gca, 'MinorGridAlpha', 0.5);
        
        params.N_RB = 100;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov__analy(i,:) = compute_average_coverage_with_coverage_threshold(params);
            Pcov_simul(i,:) = simulate_average_coverage_with_coverage_threshold(params);
        end
        subplot(1,2,2);
        f1 = plot(Threshold_dB  ,Pcov__analy , 'k-',Threshold_dB,Pcov_simul,'ro' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Coverage threshold ($\tau$)  ','Interpreter','LaTex');
        ylabel('Averaged uplink coverage probability ','Interpreter','LaTex');
        title('$N_{RB} = 100$ ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        set(gca, 'LineWidth', 2);
        set(gca, 'GridAlpha', 0.5);
        set(gca, 'MinorGridAlpha', 0.5);
        savefig('FinalResults/average_coverage_tau_N_RB');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
    case 'Truncation_Outage_Pmin'
        disp('Truncation_Outage_Pmin');
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        % of M2M nodes traffic activity
        params.aggregation_mode = 'DIRECT';
        
        Pmin_dBm = -100:2:0; % in dBm
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        
        Pm_dB = 20;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        LA_B = [10 100 1000]*1e-6 ;
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            outage_analy(i,:) = compute_truncation_outage_with_truncation_power(params);
            %outage_simul(i,:) = simulate_truncation_outage_with_truncation_power(params);
        end
        figure;
        subplot(1,2,1);
        hold on;
        for plt = 1:numel(LA_B)
            semilogx(Pmin_dBm  ,outage_analy(plt,:), params.ana_style{plt});
        end
        
        h = get(gca,'Children') ;
        deco.xlabel = 'Power truncation threshold ($P_o$~dBm)';
        deco.ylabel = 'Power truncation outage';
        deco.title = '$\alpha = 0.94, \beta = \frac{1}{2}$';
        deco.legend.items = 1:numel(LA_B) ;
        for i = 1:numel(LA_B)
            deco.legend.labels{i} = strcat(' $(\lambda_s =', num2str(LA_B(i)*1e6),'$)');
        end
        deco.legend.location = 'northwest';
        decorate_plot(h,deco);
        
        
        params.SEPL.alpha = 0.3;
        params.SEPL.beta = 2/3;
        
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            outage_analy(i,:) = compute_truncation_outage_with_Pmax(params);
            %outage_simul(i,:) = simulate_truncation_outage_with_cell_density(params);
        end
        
        subplot(1,2,2);
        hold on;
        for plt = 1:numel(LA_B)
            semilogx(Pmin_dBm ,outage_analy(plt,:), params.ana_style{plt});
        end
        h = get(gca,'Children') ;
        deco.title = '$\alpha = 0.3, \beta = \frac{2}{3}$';
        decorate_plot(h,deco);
        
        savefig('TWCResults/outage_Po_SEPL');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'Truncation_Outage_CellDensity_SEPL'
        params.simulation_area_side = [-5000 5000];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 1;
        Pm_dB = 20;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        %params.LA_B = [0.1:0.1:1 2:10]*1e-6 ;
        params.LA_B = [1:10]*1e-6 ;
        % of M2M nodes traffic activity
        params.aggregation_mode = 'DIRECT';
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        Pmin_dBm = [-120 -100 -80]; % in dBm
        Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        for i = 1:numel(Pmin)
            params.Pmin = Pmin(i);
            outage_analy(i,:) = compute_truncation_outage_with_cell_density(params);
            outage_simul(i,:) = simulate_truncation_outage_with_cell_density(params);
        end
        figure;
        subplot(1,2,1);
        hold on;
        for plt = 1:numel(Pmin)
            plot(params.LA_B*1e6  ,outage_analy(plt,:), params.ana_style{plt},params.LA_B*1e6,outage_simul(plt,:), params.sim_style{plt});
        end
        
        h = get(gca,'Children') ;
        deco.xlabel = 'Density of small cells( $\lambda_s$ cells/km$^2$)';
        deco.ylabel = 'Power truncation outage';
        deco.title = '$\alpha = 0.94, \beta = \frac{1}{2}$';
        deco.legend.items = 1:2*numel(Pmin) ;
        ind = 0;
        for i = 1:2:2*numel(Pmin)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(P_o =', num2str(Pmin_dBm(ind)),'$ dBm)');
            deco.legend.labels{i+1} = strcat('Simulation $(P_o =', num2str(Pmin_dBm(ind)),'$ dBm)');
        end
        deco.legend.location = 'northeast';
        decorate_plot(h,deco);
        
        
        params.SEPL.alpha = 0.3;
        params.SEPL.beta = 2/3;
        
        for i = 1:numel(Pmin)
            params.Pmin = Pmin(i);
            outage_analy(i,:) = compute_truncation_outage_with_cell_density(params);
            outage_simul(i,:) = simulate_truncation_outage_with_cell_density(params);
        end
        
        subplot(1,2,2);
        hold on;
        for plt = 1:numel(Pmin)
            plot(params.LA_B*1e6  ,outage_analy(plt,:), params.ana_style{plt},params.LA_B*1e6,outage_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.title = '$\alpha = 0.3, \beta = \frac{2}{3}$';
        decorate_plot(h,deco);
        
        savefig('TWCResults/outage_LA_B_SEPL');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'Truncation_Outage_Pmax_SEPL'
        params.simulation_area_side = [-5000 5000];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 1;

        
        % of M2M nodes traffic activity
        params.aggregation_mode = 'DIRECT';

        Pmin_dBm = -100; % in dBm
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        
        params.Pm = [10:10:90 100:100:900 1e3:1e3:1e4 1e4:1e4:1e5]* params.Pmin;
        
        LA_B = [10 100 1000]*1e-6 ;
        
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            outage_analy(i,:) = compute_truncation_outage_with_Pmax(params);
            %outage_simul(i,:) = simulate_truncation_outage_with_cell_density(params);
        end
        figure;
        subplot(1,2,1);
        hold on;
        for plt = 1:numel(LA_B)
            semilogx(params.Pm/ params.Pmin  ,outage_analy(plt,:), params.ana_style{plt});
        end
        
        h = get(gca,'Children') ;
        deco.xlabel = '$P_m / P_o$';
        deco.ylabel = 'Power truncation outage';
        deco.title = '$\alpha = 0.94, \beta = \frac{1}{2}$';
        deco.legend.items = 1:numel(LA_B) ;
        for i = 1:numel(LA_B)
            deco.legend.labels{i} = strcat(' $(\lambda_s =', num2str(LA_B(i)*1e6),'$)');
        end
        deco.legend.location = 'northeast';
        decorate_plot(h,deco);
        
        
        params.SEPL.alpha = 0.3;
        params.SEPL.beta = 2/3;
        
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            outage_analy(i,:) = compute_truncation_outage_with_Pmax(params);
            %outage_simul(i,:) = simulate_truncation_outage_with_cell_density(params);
        end
        
        subplot(1,2,2);
        hold on;
        for plt = 1:numel(LA_B)
            semilogx(params.Pm/ params.Pmin  ,outage_analy(plt,:), params.ana_style{plt});
        end
        h = get(gca,'Children') ;
        deco.title = '$\alpha = 0.3, \beta = \frac{2}{3}$';
        decorate_plot(h,deco);
        
        savefig('TWCResults/outage_Pmax_SEPL');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
    case 'Log_P_Po_Pm'
        disp('Log_P_Po_Pm');
        
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        % of M2M nodes traffic activity
        params.aggregation_mode = 'DIRECT';
        
        Pmin_dBm = -100:2:0; % in dBm
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        
        Pm_dB = 20;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        LA_B = [100 1000 10000]*1e-6 ;
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            ElogPPo_analy(i,:) = compute_expected_log_relative_power_with_truncation_power(params);
        end
        figure;
        subplot(1,2,1);
        f1 = plot(Pmin_dBm ,ElogPPo_analy , 'k-' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        %legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','southeast','Interpreter','LaTex');
        xlabel(' Power truncation threshold ($P_{min}$~dBm) ' ,'Interpreter','LaTex');
        ylabel('Expected Log Relative Transmit Power','Interpreter','LaTex');
        title('$P_{m} = 20$ dBm  ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        set(gca, 'LineWidth', 2);
        set(gca, 'GridAlpha', 0.5);
        set(gca, 'MinorGridAlpha', 0.5);
        
        Pm_dB = 33;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        LA_B = [100 1000 10000]*1e-6 ;
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            ElogPPo_analy(i,:) = compute_expected_log_relative_power_with_truncation_power(params);
        end
        subplot(1,2,2);
        f1 = plot(Pmin_dBm  ,ElogPPo_analy , 'k-' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        %legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','southeast','Interpreter','LaTex');
        xlabel(' Power truncation threshold ($P_{min}$~dBm) ' ,'Interpreter','LaTex');
        ylabel('Expected Log Relative Transmit Power','Interpreter','LaTex');
        title('$P_{m} = 33$ dBm  ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        set(gca, 'LineWidth', 2);
        set(gca, 'LineWidth', 2);
        set(gca, 'GridAlpha', 0.5);
        set(gca, 'MinorGridAlpha', 0.5);
        % savefig('FinalResults/coverage_LA_M_Pmin');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
    case 'Log_P_Po_SEPL'
        disp('Log_P_Po_SEPL');
        Pm_dB = 20;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        
        % of M2M nodes traffic activity
        params.aggregation_mode = 'DIRECT';
        
        Pmin_dBm = -100:2:0; % in dBm
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        LA_B = [100 1000 10000]*1e-6 ;
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            ElogPPo_analy(i,:) = compute_expected_log_relative_power_with_truncation_power(params);
        end
        figure;
        subplot(1,2,1);
        f1 = plot(Pmin_dBm ,ElogPPo_analy , 'k-' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        %legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','southeast','Interpreter','LaTex');
        xlabel(' Power truncation threshold ($P_{min}$~dBm) ' ,'Interpreter','LaTex');
        ylabel('Expected Log Relative Transmit Power','Interpreter','LaTex');
        title('$\alpha = 0.94, \beta = \frac{1}{2}$','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        set(gca, 'LineWidth', 2);
        set(gca, 'GridAlpha', 0.5);
        set(gca, 'MinorGridAlpha', 0.5);
        
        params.SEPL.alpha = 0.3;
        params.SEPL.beta = 2/3;
        LA_B = [100 1000 10000]*1e-6 ;
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            ElogPPo_analy(i,:) = compute_expected_log_relative_power_with_truncation_power(params);
        end
        subplot(1,2,2);
        f1 = plot(Pmin_dBm  ,ElogPPo_analy , 'k-' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        %legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','southeast','Interpreter','LaTex');
        xlabel(' Power truncation threshold ($P_{min}$~dBm) ' ,'Interpreter','LaTex');
        ylabel('Expected Log Relative Transmit Power','Interpreter','LaTex');
        title('$\alpha = 0.3, \beta = \frac{2}{3}$','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        set(gca, 'LineWidth', 2);
        set(gca, 'LineWidth', 2);
        set(gca, 'GridAlpha', 0.5);
        set(gca, 'MinorGridAlpha', 0.5);
        %savefig('FinalResults/coverage_LA_M_Pmin');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
    case 'CoverageThreshold_Pmin'
        params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 10;
        params.LA_B = 1000e-6   ;       % Small Cell Denisty (cells/m^2)
        Threshold_dB = (-20:1:20) ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        LA_M = [1e-2 5e-2 1e-1];
        params.rho_m = 0.1;
        Pmin_dBm = -100 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:)= compute_uplink_coverage_with_coverage_threshold(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_coverage_threshold(params);
        end
        figure;
        subplot(1,2,1);
        hold on;
        
        for plt = 1:numel(LA_M)
            plot(Threshold_dB  ,Pcov_analy(plt,:), params.ana_style{plt},Threshold_dB,Pcov_simul(plt,:), params.sim_style{plt});
        end
        
        h = get(gca,'Children') ;
        deco.xlabel = 'Coverage threshold ($\tau$ dB)  ';
        deco.ylabel = 'Uplink coverage probability ';
        deco.title = '$P_o = -100$ dBm';
        deco.legend.items = 1:2*numel(LA_M) ;
        ind = 0;
        for i = 1:2:2*numel(LA_M)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);
        %####################################################
        Pmin_dBm = -110 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:)= compute_uplink_coverage_with_coverage_threshold(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_coverage_threshold(params);
        end
        subplot(1,2,2);
        hold on;
        
        for plt = 1:numel(LA_M)
            plot(Threshold_dB  ,Pcov_analy(plt,:), params.ana_style{plt},Threshold_dB,Pcov_simul(plt,:), params.sim_style{plt});
        end
        
        h = get(gca,'Children') ;
        deco.title = '$P_o = -110$ dBm';
        deco.legend.items = 1:2*numel(LA_M) ;
        ind = 0;
        for i = 1:2:2*numel(LA_M)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);
        
        savefig('TWCResults/coverage_tau_Po');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'CoverageThreshold_SEPL'
        params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 10;
        params.LA_B = 1000e-6   ;       % Small Cell Denisty (cells/m^2)
        Threshold_dB = (-20:1:20) ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        LA_M = [1e-2 5e-2 1e-1];
        params.rho_m = 0.1;
        Pmin_dBm = -100 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_coverage_threshold(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_coverage_threshold(params);
        end
        figure;
        subplot(1,2,1);
        hold on;
        for plt = 1:numel(LA_M)
            plot(Threshold_dB  ,Pcov_analy(plt,:), params.ana_style{plt},Threshold_dB,Pcov_simul(plt,:), params.sim_style{plt});
        end
        
        h = get(gca,'Children') ;
        deco.xlabel = 'Coverage threshold ($\tau$ dB)  ';
        deco.ylabel = 'Uplink coverage probability';
        deco.title = '$\alpha = 0.94, \beta = \frac{1}{2}$';
        deco.legend.items = 1:2*numel(LA_M) ;
        ind = 0;
        for i = 1:2:2*numel(LA_M)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);
        
        params.SEPL.alpha = 0.3;
        params.SEPL.beta = 2/3;
        
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:)= compute_uplink_coverage_with_coverage_threshold(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_coverage_threshold(params);
        end
        subplot(1,2,2);
        hold on;
        
        for plt = 1:numel(LA_M)
            plot(Threshold_dB  ,Pcov_analy(plt,:), params.ana_style{plt},Threshold_dB,Pcov_simul(plt,:), params.sim_style{plt});
        end
        
        h = get(gca,'Children') ;
        deco.title = '$\alpha = 0.3, \beta = \frac{2}{3}$';
        deco.legend.items = 1:2*numel(LA_M) ;
        ind = 0;
        for i = 1:2:2*numel(LA_M)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);
        
        savefig('TWCResults/coverage_tau_SEPL');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'CoverageThreshold_N_RB'
        params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 100;
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
        
        
        
        params.N_RB = 50;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_coverage_threshold(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_coverage_threshold(params);
        end
        figure;
        subplot(1,2,1);
        f1 = plot(Threshold_dB  ,Pcov_analy , 'k-',Threshold_dB,Pcov_simul,'ro' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Coverage threshold ($\tau$)  ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$N_{RB} = 50$ ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
        params.N_RB = 100;
        
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_coverage_threshold(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_coverage_threshold(params);
        end
        subplot(1,2,2);
        f1 = plot(Threshold_dB  ,Pcov_analy , 'k-',Threshold_dB,Pcov_simul,'ro' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Coverage threshold ($\tau$)  ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$N_{RB} = 100$ ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
        
    case 'CoverageThreshold_Pm'
        params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 100;
        params.time_slots = 10;
        params.LA_B = 1000e-6   ;       % Small Cell Denisty (cells/m^2)
        params.Beta_Distribution = [3 4];   % Beta Distribution parameters                                   % of M2M nodes traffic activity
        Threshold_dB = (-20:1:20) ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        LA_M = [1e-2 5e-2 1e-1];
        params.rho_m = 0.1;
        Pmin_dBm = -80 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        
        Pm_dB = 20;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            [Pcov_approx_analy(i,:) , Pcov_exact_analy(i,:)] = compute_uplink_coverage_with_coverage_threshold(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_coverage_threshold(params);
        end
        figure;
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
    case 'CellDensity_Scaling_N_RB'
        disp('CellDensity_Scaling_N_RB');
        params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 10;
        params.time_slots = 10;
        Pm_dB = 20;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
                                  % of M2M nodes traffic activity
        params.aggregation_mode = 'DIRECT';
        Threshold_dB = 0 ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        params.LA_B = [100 300 500 800 1000 1500 2000 2500 3000 4000 5000 6000 7000 8000 9000 10000]*1e-6 ;       % Small Cell Denisty (cells/m^2)
        %params.LA_B = [100 300 500]*1e-6 ;
        Pmin_dBm = -100 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        
        params.rho_m = 0.1;  
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Nmtc_active(i,:) = compute_Nmtc_active_with_smallcells_density(params);
            Nmtc_active(i,:) = simulate_Nmtc_active_with_small_cell_density(params);
        end
        figure;
        subplot(1,2,1);
        hold on;
        for plt = 1:numel(LA_M)
            plot( params.LA_B*1e6  ,Pcov_analy(plt,:), params.ana_style{plt},params.LA_B*1e6,Pcov_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.xlabel = 'Density of small cells( $\lambda_s$ cells/km$^2$)';
        deco.ylabel = 'Number of active MTC nodes';
        deco.title = '$\rho_m = 0.1$ ';
        deco.legend.items = 1:2*numel(LA_M) ;
        ind = 0;
        for i = 1:2:2*numel(LA_M)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);
        
        Pmin_dBm = -110 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Nmtc_active(i,:) = compute_Nmtc_active_with_smallcells_density(params);
            Nmtc_active(i,:) = simulate_Nmtc_active_with_small_cell_density(params);
        end
        subplot(1,2,2);
        hold on;
        for plt = 1:numel(LA_M)
            plot(params.LA_B*1e6  ,Pcov_analy(plt,:), params.ana_style{plt},params.LA_B*1e6,Pcov_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.title = '$\rho_m = 0.3$ ';
        decorate_plot(h,deco);
        
        savefig('TWCResults/coverage_LA_B_Po');     
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'CellDensity_Pmin'
        disp('CellDensity_Pmin');
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
        params.LA_B = [100 300 500 800 1000 1500 2000 2500 3000 4000 5000 6000 7000 8000 9000 10000]*1e-6 ;       % Small Cell Denisty (cells/m^2)
        %params.LA_B = [100 300 500]*1e-6 ;
        Pmin_dBm = -100 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_smallcells_density(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_small_cell_density(params);
        end
        figure;
        subplot(1,2,1);
        hold on;
        for plt = 1:numel(LA_M)
            plot( params.LA_B*1e6  ,Pcov_analy(plt,:), params.ana_style{plt},params.LA_B*1e6,Pcov_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.xlabel = 'Density of small cells( $\lambda_s$ cells/km$^2$)';
        deco.ylabel = 'Uplink coverage probability';
        deco.title = '$P_o = -100$ dBm';
        deco.legend.items = 1:2*numel(LA_M) ;
        ind = 0;
        for i = 1:2:2*numel(LA_M)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);
        
        Pmin_dBm = -110 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_smallcells_density(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_small_cell_density(params);
        end
        subplot(1,2,2);
        hold on;
        for plt = 1:numel(LA_M)
            plot(params.LA_B*1e6  ,Pcov_analy(plt,:), params.ana_style{plt},params.LA_B*1e6,Pcov_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.title = '$P_o = -110$ dBm';
        decorate_plot(h,deco);
        
        savefig('TWCResults/coverage_LA_B_Po');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'CellDensity_N_RB'
        disp('CellDensity_N_RB');
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
        params.LA_B = [100 300 500 800 1000 1500 2000 2500 3000 4000 5000 6000 7000 8000 9000 10000]*1e-6 ;       % Small Cell Denisty (cells/m^2)
        %params.LA_B = [1000 3000 5000]*1e-6 ;
        Pmin_dBm = -70 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        
        params.N_RB = 50;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_smallcells_density(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_small_cell_density(params);
        end
        figure;
        subplot(1,2,1);
        f1 = plot(params.LA_B*1e6  ,Pcov_analy , 'k-',params.LA_B*1e6,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Small cell density (cells/$km^2$)  ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$N_{RB} = 50$ ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        set(gca, 'LineWidth', 2);
        set(gca, 'GridAlpha', 0.5);
        set(gca, 'MinorGridAlpha', 0.5);
        
        params.N_RB = 100;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_smallcells_density(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_small_cell_density(params);
        end
        subplot(1,2,2);
        f1 = plot(params.LA_B*1e6  ,Pcov_analy , 'k-',params.LA_B*1e6,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Small cell density (cells/$km^2$)  ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$N_{RB} = 100$ ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        set(gca, 'LineWidth', 2);
        set(gca, 'GridAlpha', 0.5);
        set(gca, 'MinorGridAlpha', 0.5);
        print('FinalResults/coverage_LA_B_N_RB', '-depsc');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'CellDensity_SEPL'
        disp('CellDensity_SEPL');
        params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 10;
        params.time_slots = 10;
        Pm_dB = 20;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        params.rho_m = 0.1;                                    % of M2M nodes traffic activity
        params.aggregation_mode = 'DIRECT';
        Threshold_dB = 0 ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        
        
        params.LA_B = [100 300 500 800 1000 1500 2000 2500 3000 4000 5000 6000 7000 8000 9000 10000]*1e-6   ;       % Small Cell Denisty (cells/m^2)
        %params.LA_B = [1000 3000 5000]*1e-6   ;
        Pmin_dBm = -70 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_smallcells_density(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_small_cell_density(params);
        end
        figure;
        subplot(1,2,1);
        hold on;
        for plt = 1:numel(LA_M)
            plot( params.LA_B*1e6  ,Pcov_analy(plt,:), params.ana_style{plt},params.LA_B*1e6,Pcov_simul(plt,:), params.sim_style{plt});
        end
        
        h = get(gca,'Children') ;
        deco.xlabel = 'Density of small cells( $\lambda_s$ cells/km$^2$)';
        deco.ylabel = 'Uplink coverage probability';
        deco.title = '$\alpha = 0.94, \beta = \frac{1}{2}$';
        deco.legend.items = 1:2*numel(LA_M) ;
        ind = 0;
        for i = 1:2:2*numel(LA_M)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);        subplot(1,2,1);
        
        params.SEPL.alpha = 0.3;
        params.SEPL.beta = 2/3;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_smallcells_density(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_small_cell_density(params);
        end
        subplot(1,2,2);
        hold on;
        for plt = 1:numel(LA_M)
            plot(params.LA_B*1e6  ,Pcov_analy(plt,:), params.ana_style{plt},params.LA_B*1e6,Pcov_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.title = '$\alpha = 0.3, \beta = \frac{2}{3}$';
        decorate_plot(h,deco);
        
        savefig('TWCResults/coverage_LA_B_SEPL');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'CellDensity_Pm'
        disp('CellDensity_Pm');
        params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 10;
        params.time_slots = 10;
        params.rho_m = 0.1;                                    % of M2M nodes traffic activity
        params.aggregation_mode = 'DIRECT';
        Threshold_dB = 0 ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        
        params.LA_B = [100 300 500 800 1000 1500 2000 2500 3000 4000 5000 6000 7000 8000 9000 10000]*1e-6   ;       % Small Cell Denisty (cells/m^2)
        %params.LA_B = [1000 3000 5000]*1e-6   ;
        Pmin_dBm = -80 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        
        Pm_dB = 20;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_smallcells_density(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_small_cell_density(params);
        end
        figure
        subplot(1,2,1);
        f1 = plot(params.LA_B*1e6  ,Pcov_analy , 'k-',params.LA_B*1e6,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Small cell density (cells/$km^2$)  ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$P_m = 20$ dBm ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
        Pm_dB = 33;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_smallcells_density(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_small_cell_density(params);
        end
        subplot(1,2,2);
        f1 = plot(params.LA_B*1e6  ,Pcov_analy , 'k-',params.LA_B*1e6,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Density of small cells(cells/$km^2$)  ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$P_m = 33$ dBm ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        savefig('FinalResults/coverage_LA_B_Pm');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'M2MDensity_Pmin'
        disp('M2MDensity_Pmin');
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
        params.LA_M = [10000:5000:50000 60000:10000:100000]*1e-6 ;       % Small Cell Denisty (cells/m^2)
        
        Pmin_dBm = -100 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        LA_B = [1000 5000 10000]*1e-6 ;
        %LA_B = 100*1e-6 ;
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_m2m_density(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_m2m_density(params);
        end
        figure;
        subplot(1,2,1);
        hold on;
        for plt = 1:numel(LA_B)
            semilogx( params.LA_M*1e6  ,Pcov_analy(plt,:), params.ana_style{plt},params.LA_M*1e6,Pcov_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.xlabel = 'Density of mMTC nodes($\lambda_m$ nodes/km$^2$)';
        deco.ylabel = 'Uplink coverage probability';
        deco.title = '$P_o = -100$ dBm';
        deco.legend.items = 1:2*numel(LA_B) ;
        ind = 0;
        for i = 1:2:2*numel(LA_B)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);
        
        Pmin_dBm = -110 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_m2m_density(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_m2m_density(params);
        end
        subplot(1,2,2);
        hold on;
        for plt = 1:numel(LA_B)
            semilogx(params.LA_M*1e6  ,Pcov_analy(plt,:), params.ana_style{plt},params.LA_M*1e6,Pcov_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.title = '$P_o = -110$ dBm';
        decorate_plot(h,deco);
        
        savefig('TWCResults/coverage_LA_M_Po');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'M2MDensity_SEPL'
        disp('M2MDensity_SEPL');
        params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 10;
        params.time_slots = 10;
        Pm_dB = 20;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        params.rho_m = 0.1;                                    % of M2M nodes traffic activity
        params.aggregation_mode = 'DIRECT';
        Threshold_dB = 0 ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        
        
        params.LA_M = [10000:5000:50000 60000:10000:100000]*1e-6 ;       % Small Cell Denisty (cells/m^2)
        %params.LA_B = [1000 3000 5000]*1e-6   ;
        Pmin_dBm = -100 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        
        LA_B = [1000 5000 10000]*1e-6 ;
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_m2m_density(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_m2m_density(params);
        end
        figure;
        subplot(1,2,1);
        hold on;
        for plt = 1:numel(LA_B)
            semilogx( params.LA_M*1e6  ,Pcov_analy(plt,:), params.ana_style{plt},params.LA_M*1e6,Pcov_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.xlabel = 'Density of mMTC nodes($\lambda_m$ nodes/km$^2$)';
        deco.ylabel = 'Uplink coverage probability';
        deco.title = '$\alpha = 0.94, \beta = \frac{1}{2}$';
        deco.legend.items = 1:2*numel(LA_B) ;
        ind = 0;
        for i = 1:2:2*numel(LA_B)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_s =', num2str(LA_B(ind)*1e6),'$)');
        end
        deco.legend.location = 'southwest';
        decorate_plot(h,deco);
        
        
        params.SEPL.alpha = 0.3;
        params.SEPL.beta = 2/3;
        
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_m2m_density(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_m2m_density(params);
        end
        subplot(1,2,2);
        hold on;
        for plt = 1:numel(LA_B)
            semilogx(params.LA_M*1e6  ,Pcov_analy(plt,:), params.ana_style{plt},params.LA_M*1e6,Pcov_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.title = '$\alpha = 0.3, \beta = \frac{2}{3}$';
        decorate_plot(h,deco);
        
        savefig('TWCResults/coverage_LA_M_SEPL');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'M2MDensity_Pm'
        disp('M2MDensity_Pm');
        params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 10;
        params.time_slots = 10;
        params.rho_m = 0.1;                                    % of M2M nodes traffic activity
        params.aggregation_mode = 'DIRECT';
        Threshold_dB = 0 ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        
        params.LA_M = [10000:5000:50000 60000:10000:100000]*1e-6 ;
        %params.LA_B = [1000 3000 5000]*1e-6   ;
        Pmin_dBm = -80 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        
        Pm_dB = 20;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        
        LA_B = [100 1000 10000]*1e-6 ;
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_m2m_density(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_m2m_density(params);
        end
        figure
        subplot(1,2,1);
        f1 = plot(params.LA_M*1e6  ,Pcov_analy , 'k-',params.LA_M*1e6,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','southeast','Interpreter','LaTex');
        xlabel('mMTC nodes density (nodes/$km^2$)  ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$P_m = 20$ dBm ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        set(gca, 'LineWidth', 2);
        set(gca, 'GridAlpha', 0.5);
        set(gca, 'MinorGridAlpha', 0.5);
        
        Pm_dB = 33;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        
        
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_m2m_density(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_m2m_density(params);
        end
        subplot(1,2,2);
        f1 = plot(params.LA_M*1e6  ,Pcov_analy , 'k-',params.LA_M*1e6,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','southeast','Interpreter','LaTex');
        xlabel('mMTC nodes density (nodes/$km^2$)  ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$P_m = 33$ dBm ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        set(gca, 'LineWidth', 2);
        set(gca, 'GridAlpha', 0.5);
        set(gca, 'MinorGridAlpha', 0.5);
        savefig('FinalResults/coverage_LA_M_Pm');
        
        
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'M2MDensity_N_RB'
        disp('M2MDensity_N_RB');
        params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 10;
        params.time_slots = 10;
        params.rho_m = 0.1;                                    % of M2M nodes traffic activity
        params.aggregation_mode = 'DIRECT';
        Threshold_dB = 0 ; % dB
        params.Threshold = 10.^(Threshold_dB/10);
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        
        params.LA_M = [10000:5000:50000 60000:10000:100000]*1e-6 ;
        %params.LA_B = [1000 3000 5000]*1e-6   ;
        Pmin_dBm = -80 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        
        Pm_dB = 20;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        
        params.N_RB = 50;
        
        LA_B = [100 1000 10000]*1e-6 ;
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_m2m_density(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_m2m_density(params);
        end
        figure
        subplot(1,2,1);
        f1 = plot(params.LA_M*1e6  ,Pcov_analy , 'k-',params.LA_M*1e6,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','southeast','Interpreter','LaTex');
        xlabel('mMTC nodes density (nodes/$km^2$)  ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$N_{RB} = 50$  ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        set(gca, 'LineWidth', 2);
        set(gca, 'GridAlpha', 0.5);
        set(gca, 'MinorGridAlpha', 0.5);
        
        params.N_RB = 100;
        
        
        for i = 1:numel(LA_B)
            params.LA_B = LA_B(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_m2m_density(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_m2m_density(params);
        end
        subplot(1,2,2);
        f1 = plot(params.LA_M*1e6  ,Pcov_analy , 'k-',params.LA_M*1e6,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','southeast','Interpreter','LaTex');
        xlabel('mMTC nodes density (nodes/$km^2$)  ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$N_{RB} = 100$','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        set(gca, 'LineWidth', 2);
        set(gca, 'GridAlpha', 0.5);
        set(gca, 'MinorGridAlpha', 0.5);
        savefig('FinalResults/coverage_LA_M_N_RB');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
        
    case 'ResourceBlocks_Pmin'
        disp('ResourceBlocks_Pmin');
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
        params.LA_B = 1000*1e-6 ;       % Small Cell Denisty (cells/m^2)
        
        params.N_RB = 10:10:200;
        Pmin_dBm = -100 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_resource_blocks(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_resource_blocks(params);
        end
        figure
        subplot(1,2,1);
        hold on;
        for plt = 1:numel(LA_M)
            plot(params.N_RB  ,Pcov_analy(plt,:), params.ana_style{plt},params.N_RB,Pcov_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.xlabel = 'Number of orthogonal channels($N_{RB}$)';
        deco.ylabel = 'Uplink coverage probability';
        deco.title = '$P_o = -100$ dBm';
        deco.legend.items = 1:2*numel(LA_M) ;
        ind = 0;
        for i = 1:2:2*numel(LA_M)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
        end
        deco.legend.location = 'southeast';
        decorate_plot(h,deco);
        
        Pmin_dBm = -110 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_resource_blocks(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_resource_blocks(params);
        end
        subplot(1,2,2);
        hold on;
        for plt = 1:numel(LA_M)
            plot(params.N_RB  ,Pcov_analy(plt,:), params.ana_style{plt},params.N_RB,Pcov_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.title = '$P_o = -110$ dBm';
        decorate_plot(h,deco);
        
        savefig('TWCResults/coverage_N_RB_Po');
        
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'ResourceBlocks_SEPL'
        disp('ResourceBlocks_SEPL');
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
        params.LA_B = 1000*1e-6 ;       % Small Cell Denisty (cells/m^2)
        
        params.N_RB = 10:10:200;
        Pmin_dBm = -70 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_resource_blocks(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_resource_blocks(params);
        end
        figure
        subplot(1,2,1);
        f1 = plot(params.N_RB ,Pcov_analy , 'k-',params.N_RB,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Number of orthogonal resource blocks($N_{RB}$)' ,'Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$\alpha = 0.94, \beta = \frac{1}{2}$','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
        params.SEPL.alpha = 0.3;
        params.SEPL.beta = 2/3;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_resource_blocks(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_resource_blocks(params);
        end
        subplot(1,2,2);
        f1 = plot(params.N_RB,Pcov_analy , 'k-',params.N_RB,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Number of orthogonal resource blocks($N_{RB}$) ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$\alpha = 0.3, \beta = \frac{2}{3}$','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        savefig('FinalResults/coverage_N_RB_SEPL');
        
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'ResourceBlocks_Pm'
        disp('ResourceBlocks_Pm');
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
        params.LA_B = 1000*1e-6 ;       % Small Cell Denisty (cells/m^2)
        
        params.N_RB = 10:10:200;
        Pmin_dBm = -70 ;
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        
        Pm_dB = 20;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_resource_blocks(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_resource_blocks(params);
        end
        figure
        subplot(1,2,1);
        f1 = plot(params.N_RB ,Pcov_analy , 'k-',params.N_RB,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Number of orthogonal resource blocks($N_{RB}$)' ,'Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$P_m = 20$ dBm ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
        Pm_dB = 33;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_resource_blocks(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_resource_blocks(params);
        end
        subplot(1,2,2);
        f1 = plot(params.N_RB,Pcov_analy , 'k-',params.N_RB,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Number of orthogonal resource blocks($N_{RB}$) ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$P_m = 33$ dBm ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        savefig('FinalResults/coverage_N_RB_Pm');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'PowerTruncationThreshold_SEPL'
        disp('PowerTruncationThreshold_SEPL');
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
        params.LA_B = 1000*1e-6 ;       % Small Cell Denisty (cells/m^2)
        
        params.N_RB = 100;
        Pmin_dBm = -100:2:-40; % in dBm
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        
        params.SEPL.alpha = 0.94;
        params.SEPL.beta = 0.5;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_power_truncation_threshold(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_power_truncation_threshold(params);
        end
        figure
        subplot(1,2,1);
        f1 = plot(Pmin_dBm ,Pcov_analy , 'k-',Pmin_dBm,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel(' Power truncation threshold $P_{min}$(dBm) ' ,'Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$\alpha = 0.94, \beta = \frac{1}{2}$','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
        params.SEPL.alpha = 0.3;
        params.SEPL.beta = 2/3;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_power_truncation_threshold(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_power_truncation_threshold(params);
        end
        subplot(1,2,2);
        f1 = plot(Pmin_dBm,Pcov_analy , 'k-',Pmin_dBm,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel(' Power truncation threshold $P_{min}$(dBm) ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$\alpha = 0.3, \beta = \frac{2}{3}$','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        savefig('FinalResults/coverage_Pmin_SEPL');
        
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'PowerTruncationThreshold_Pm'
        disp('PowerTruncationThreshold_Pm');
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
        params.LA_B = 1000*1e-6 ;       % Small Cell Denisty (cells/m^2)
        
        params.N_RB = 100;
        
        Pmin_dBm = -100:2:-40; % in dBm
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        
        Pm_dB = 20;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_power_truncation_threshold(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_power_truncation_threshold(params);
        end
        figure
        subplot(1,2,1);
        f1 = plot(Pmin_dBm ,Pcov_analy , 'k-',Pmin_dBm,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel(' Power truncation threshold $P_{min}$(dBm) ' ,'Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$P_m = 20$ dBm','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
        Pm_dB = 33;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_power_truncation_threshold(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_power_truncation_threshold(params);
        end
        subplot(1,2,2);
        f1 = plot(Pmin_dBm,Pcov_analy , 'k-',Pmin_dBm,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel(' Power truncation threshold $P_{min}$(dBm) ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$P_m = 33$ dBm','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        savefig('FinalResults/coverage_Pmin_Pm');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
    case 'PowerTruncationThreshold_N_RB'
        disp('PowerTruncationThreshold_N_RB');
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
        params.LA_B = 1000*1e-6 ;       % Small Cell Denisty (cells/m^2)
        
        
        
        Pmin_dBm = -100:2:-40; % in dBm
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        
        params.N_RB = 50;
        
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_power_truncation_threshold(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_power_truncation_threshold(params);
        end
        figure
        subplot(1,2,1);
        f1 = plot(Pmin_dBm ,Pcov_analy , 'k-',Pmin_dBm,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel(' Power truncation threshold $P_{min}$(dBm) ' ,'Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$N_{RB} = 50$ ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
        params.N_RB = 100;
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_power_truncation_threshold(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_power_truncation_threshold(params);
        end
        subplot(1,2,2);
        f1 = plot(Pmin_dBm,Pcov_analy , 'k-',Pmin_dBm,Pcov_simul,'r.' );
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel(' Power truncation threshold $P_{min}$(dBm) ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('$N_{RB} = 100$ ','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        savefig('FinalResults/coverage_Pmin_N_RB');
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
    case 'PowerTruncationThreshold_Cell_Density'
        disp('PowerTruncationThreshold_Cell_Density');
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
        
        
        
        
        Pmin_dBm = [-130:2:-100 -95:5:-70] ; % in dBm
        params.Pmin = 10.^(Pmin_dBm/10) * 1e-3;
        
        params.LA_B = 500*1e-6 ;       % Small Cell Denisty (cells/m^2)
        
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_power_truncation_threshold(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_power_truncation_threshold(params);
        end
        figure
        subplot(1,2,1);
        hold on;
        for plt = 1:numel(LA_M)
            plot( Pmin_dBm  ,Pcov_analy(plt,:), params.ana_style{plt},Pmin_dBm,Pcov_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.xlabel = 'Power truncation threshold ($P_o$ dBm) ';
        deco.ylabel = 'Uplink coverage probability';
        deco.title = '$\lambda_s = 500$ cells/km$^2$';
        deco.legend.items = 1:2*numel(LA_M) ;
        ind = 0;
        for i = 1:2:2*numel(LA_M)
            ind = ind + 1;
            deco.legend.labels{i} = strcat('Analysis $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
            deco.legend.labels{i+1} = strcat('Simulation $(\lambda_m =', num2str(LA_M(ind)*1e6),'$)');
        end
        deco.legend.location = 'southeast';
        decorate_plot(h,deco);
        
        
        params.LA_B = 5000*1e-6 ;       % Small Cell Denisty (cells/m^2)
        LA_M = [1e-2 5e-2 1e-1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            Pcov_analy(i,:) = compute_uplink_coverage_with_power_truncation_threshold(params);
            Pcov_simul(i,:) = simulate_uplink_coverage_with_power_truncation_threshold(params);
        end
        subplot(1,2,2);
        hold on;
        for plt = 1:numel(LA_M)
            plot( Pmin_dBm  ,Pcov_analy(plt,:), params.ana_style{plt},Pmin_dBm,Pcov_simul(plt,:), params.sim_style{plt});
        end
        h = get(gca,'Children') ;
        deco.title = '$\lambda_s = 5000$ cells/km$^2$';
        decorate_plot(h,deco);
        
        savefig('TWCResults/coverage_P_o_LA_B');
        
        
        
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
        
        
        
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
        params.simulation_area_side = [-500 500];    % simulation area side to side
        params.space_realizations = 10;
        params.time_slots = 10;
        Pm_dB = -20;
        params.Pm = 10.^(Pm_dB/10) * 1e-3;
        
        Pmin_dB = -110;
        params.Pmin = 10.^(Pmin_dB/10) * 1e-3;
        
        
        params.LA_B = 1000e-6   ;
        params.LA_M = 1;
        params.Beta_Distribution = [3 4];   % Beta Distribution parameters
        [f_analy ,  x] = compute_pdf_of_P(params);
        [f_simul , y , P ]= simulate_pdf_of_P(params);
        f1 = plot(x,f_analy , 'k-' ,y,f_simul , 'rx' );
        hold on;
        
        xlim([0 1e-11])
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