function [Op_Analy ]  = compute_truncation_outage_with_cell_density(params)
a = params.SEPL.alpha;
b = params.SEPL.beta;
Pu= params.Pm;
Po = params.Pmin;
theta = log((Pu./Po).^(1/a)).^(1/b);

Op_Analy =  exp(- pi .* params.LA_B .* theta.^2);
end