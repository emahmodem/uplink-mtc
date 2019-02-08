function [ output_args ] = gen_clustering_coverage_results_mMTC( input_args )
%GEN_CLUSTERING_COVERAGE_RESULTS_MMTC Summary of this function goes here
%   Author : Mahmoud Kamel
%   Date :  21-Nov-2017 10:02:13
%   This function generates results of simulating the effect of clustering
%   on the load of cells and coverage
%/////////////////////////////////////////////////////////////////////////
clc;close all;clear all;
% Simulation Parameters
params.Pm = 100e-3;                          % Uplink Transmit power of M2M nodes in watts 20 dBm
params.Pmin = 1e-6;                          % Min. received power at the base station in watts > BS sensitivity 
params.simulation_area_side = [-50 50];    % simulation area side to side
params.space_realizations = 10;
params.time_slots = 10;
params.BW = 20e6;
params.N_RB = 100;
params.N = 10;    % Orthogonal carriers allocated to the RACH transmission of M2M traffic 
params.RB_BW = 180e3; %Hz
params.No = 10^(-17.4) * params.RB_BW;
params.SEPL.alpha = 0.94;
params.SEPL.beta = 0.5;

%/////////////////////////////////////////////////////////////////////////
% Simulation Switches


sim = 'UplinkCoverage_SmallCellsDensity';
%sim = 'NetworkVisualize';
%sim = 'Loading';
%sim = 'ProbOfH2H_ClusterHeads_H2HDensity';
%sim = 'CellDensity';
%sim = 'M2MDensity';
%sim = 'ClusterHeadDensity';
%sim = 'CoverageThreshold';
switch(sim)
    case 'UplinkCoverage_SmallCellsDensity'
        params.LA_B = [1000e-6 5000e-6 10000e-6]   ;       % Small Cell Denisty (cells/m^2)
        params.Beta_Distribution = [3 4];   % Beta Distribution parameters
                                            % of M2M nodes traffic activity
        params.aggregation_mode = 'DIRECT';
         params.Threshold = -5 ; % dB
        LA_M = [1e-1 5e-1 1];
        for i = 1:numel(LA_M)
            params.LA_M = LA_M(i);
            [Pcov(i,:)] =simulate_uplink_coverage_with_small_cell_density(params);
        end
        figure(1)
        subplot(1,2,1);
        f1 = semilogx(params.LA_B ,Pcov);
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1(1:3),{'$\lambda_m  = 10^5~nodes/km^2$' , '$\lambda_m  = 5 \times 10^5~nodes/km^2$' , '$\lambda_m  = 10^6~nodes/km^2$'},'FontSize',25,'FontWeight','bold','Location','northwest','Interpreter','LaTex');
        xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
        ylabel('Uplink coverage probability ','Interpreter','LaTex');
        title('No aggregation')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
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
function Pcov = simulate_uplink_coverage_with_small_cell_density(params)
simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.LA_B);
PoCoverage = zeros(points,  params.space_realizations , params.time_slots );

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
   
    server_m2m(2,:)= randi(params.N,1,N_users_M2M);   % each M2M node selects a carrier from N carriers randomly (uniform distribution)
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
            intra_interferers = find(same_ser + same_rb == 2);
            inter_interferers = find(another_ser + same_rb == 2);
            S_m(i) = params.Pmin * Hm(ser,i);
            I_intra =  sum(params.Pmin * Hm(ser,intra_interferers));
            distance_to_interferers = pdist2(locations.BS(ser,:),locations.M2M(inter_interferers,:),'euclidean') ;
            I_inter = sum( server_m2m_P(inter_interferers) .* Hm(ser,inter_interferers).* exp(-params.SEPL.alpha .* distance_to_interferers .^ params.SEPL.beta) );
            %I_m(i) = I_intra + I_inter; 
            I_m(i) =  I_inter;  % special case #1 on intra interference 
        end   
        SIR_m = S_m./ I_m;
        SIR_m_dB = 10*log10(SIR_m);
        PoCoverage(p,m,t) = sum(SIR_m_dB > params.Threshold) / N_users_M2M;
        end
        
    end
    
    normfact = params.space_realizations * params.time_slots ;% * simulation_area;
    Pcov = sum(sum(PoCoverage,3),2) / normfact;
end

end
%/////////////////////////////////////////////////////////////////////////
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
