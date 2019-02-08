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