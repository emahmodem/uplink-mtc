function [Ep ]  = compute_expected_log_relative_power_with_truncation_power(params)
a = params.SEPL.alpha;
b = params.SEPL.beta;
Pu= params.Pm;
Po = params.Pmin;
theta = log((Pu./Po).^(1/a)).^(1/b);
k = 1;
Ep =  a^k .* gamma (b*k/2 + 1) .* gammainc( pi .* params.LA_B .* theta.^2 , b*k/2 + 1 , 'lower') ./ ((pi .* params.LA_B) .^(b*k/2) .* (1 - exp(- pi .* params.LA_B .* theta.^2)));
end