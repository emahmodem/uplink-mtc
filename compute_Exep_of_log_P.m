function [Ep ] = compute_Exep_of_log_P(params)
a = params.SEPL.alpha;
b = params.SEPL.beta;
Pu= params.Pm;
Po = params.Pmin;
theta = log((Pu/Po).^(1/a)).^(1/b);
% Ep = log(params.Pmin) +  (a .* gamma (b/2 + 1) .* gammainc( pi .* params.LA_B .* theta^2 , b/2 + 1 ,'lower')) ./ ((pi .* params.LA_B) .^(b/2) .*(1 - exp(- pi .* params.LA_B .* theta^2)));
Ep = (a .* gamma (b/2 + 1) .* gammainc( pi .* params.LA_B .* theta^2 , b/2 + 1 ,'lower')) ./ (((pi .* params.LA_B) .^(b/2)) .* (1 - exp(- pi .* params.LA_B .* theta^2)));
end