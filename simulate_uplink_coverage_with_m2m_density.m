function Pcov = simulate_uplink_coverage_with_m2m_density(params)
simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.LA_M);
%PoCoverage = zeros(points,  params.space_realizations , params.time_slots );
PoCoverage = zeros(params.space_realizations , params.time_slots );


for p = 1:points
    %     if(params.LA_B(p) > 5000)
    %         params.simulation_area_side = [-350 350];
    %         simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
    %     end
    fprintf('\n')
    %disp(['Done: ' , num2str(params.LA_B(p))]);
    disp(['Smallcell Density: ' , num2str(params.LA_B)]);
    disp(['M2M Density: ' , num2str(params.LA_M(p))]);
    disp(['Path Loss Parameters: ' ,'alpha=', num2str(params.SEPL.alpha) ,' beta=', num2str(params.SEPL.beta)]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_b = params.LA_B * simulation_area;
        %rho_m = betarnd(params.Beta_Distribution(1),params.Beta_Distribution(2));
        mu_m = params.LA_M(p) * simulation_area;
        
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
                I_m(i) =  I_inter;  % special case #1 on intra interference
            end
            SINR_m = S_m./ (I_m + params.No);
            %SIR_m_dB = 10*log10(SIR_m);
            PoCoverage(m,t) = sum(SINR_m > params.Threshold) / N_users_M2M;
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