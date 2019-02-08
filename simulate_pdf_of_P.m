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
pts = [linspace(params.Pmin,9*params.Pmin,10) linspace(10*params.Pmin,90*params.Pmin,10) linspace(100*params.Pmin,900*params.Pmin,10) linspace(1000*params.Pmin,params.Pm,10)  ];
%pts = params.Pmin:params.Pmin/20:params.Pm;
%factor = params.space_realizations * params.time_slots;
[f_simul ,x] = ksdensity(P,pts);
%f_simul = f_simul./factor;
end
