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