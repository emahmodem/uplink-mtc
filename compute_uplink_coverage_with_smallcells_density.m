function P_exact = compute_uplink_coverage_with_smallcells_density(params)
SNR = params.Pmin  / params.No;
noise_term = exp(- params.Threshold / SNR);
a = params.SEPL.alpha;
b = params.SEPL.beta;
Pu= params.Pm;
Po = params.Pmin;
n = (2/b) - 1;
t = params.Threshold;

theta = log((Pu./Po).^(1/a)).^(1/b);
F = @(y,v,tau) log(tau ./ y).^(v)  ./ ( y + 1);
% compute E_P(ln^k(P))
%A = (1 -  exp(- pi .* params.LA_B .* theta^2) ).^-1;
A = (n+1) * (pi / a^(n+1));
B = params.rho_m * (params.LA_M / params.N_RB);
for k = 0:n
    Ep(k+1,:) =  a^k .* gamma (k/(n+1) + 1) .* gammainc( pi .* params.LA_B .* theta.^2 , k/(n+1) + 1 , 'lower') ./ ((pi .* params.LA_B) .^(k/(n+1)) .* (1 - exp(- pi .* params.LA_B .* theta.^2)));
    NK(k+1) = nchoosek(n,k);
    J_v(k+1) = integral(@(y)F(y,n - k,t),0,t);
    T(k+1,:) = NK(k+1) *  J_v(k+1) * Ep(k+1,:) ;
end

S = sum(T,1);
interference_term = exp(- A * B * S);
P_exact = noise_term .* interference_term;
end