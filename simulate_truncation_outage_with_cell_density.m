function [Op_Simul ]  = simulate_truncation_outage_with_cell_density(params)
simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.LA_B);
%PoCoverage = zeros(points,  params.space_realizations , params.time_slots );
Outage = zeros(1,params.space_realizations);
for p = 1:points
    simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
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
        Outage(m) = sum(server_m2m_P == 0) / N_users_M2M;
        %server_m2m(3,:) = server_m2m_P; 
    end
    
    normfact = params.space_realizations ;% * simulation_area;
    Op_Simul(p) = sum(Outage) / normfact
    %     fprintf(resultfile,'%s\n',char(datetime));
    %     fprintf(resultfile,'%s\n',num2str(params.LA_M));
    %     fprintf(resultfile,'%s\n',num2str(params.LA_B));
    %     fprintf(resultfile,'%s\n',num2str(Pcov));
end

end