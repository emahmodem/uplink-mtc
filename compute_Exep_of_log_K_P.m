function [Ep] = compute_Exep_of_log_K_P(params,K)
Ep = zeros(K,numel(params.LA_B));
a = params.SEPL.alpha;
b = params.SEPL.beta;
Pu= params.Pm;
Po = params.Pmin;
theta = log((Pu/Po).^(1/a)).^(1/b);
n = (2/b) - 1;
A = (1 -  exp(- pi .* params.LA_B .* theta^2) ).^-1;
for k = 0:K
    Ep(k+1,:) =  a^k .* gamma (k/(n+1) + 1) .* gammainc( pi .* params.LA_B .* theta^2 , k/(n+1) + 1 ,'lower') ./ ((pi .* params.LA_B) .^(k/(n+1)) .* (1 - exp(- pi .* params.LA_B .* theta^2)));
end
end