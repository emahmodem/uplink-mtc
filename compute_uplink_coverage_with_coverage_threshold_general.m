function P_cov = compute_uplink_coverage_with_coverage_threshold_general(params)
SNR = params.Pmin  / params.No;
noise_term = exp(- params.Threshold / SNR);
a = params.SEPL.alpha;
b = params.SEPL.beta;
t = params.Threshold;
r = (2/b) - 1;
F = @(y,v,tau) log(tau ./ y).^(v)  ./ ( y + 1);
Pu= params.Pm;
Po = params.Pmin;
theta = log((Pu/Po).^(1/a)).^(1/b);
A = (2 * pi) / (b * a^(2/b));
B = params.rho_m * (params.LA_M / params.N_RB);
for k = 0:params.Nterms 
    Ep(k+1) =  a^k .* gamma (b*k/2 + 1) .* gammainc( pi .* params.LA_B .* theta^2 , b*k/2 + 1 , 'lower') ./ ((pi .* params.LA_B) .^(b*k/2) .* (1 - exp(- pi .* params.LA_B .* theta^2)));  
    RK(k+1) = pochhammer(r,k) / factorial(k);
    for p = 1:numel(t)
        J_v(k+1,p) = integral(@(y)F(y,r-k,t(p)),0,t(p));
    end
    T(k+1,:) = RK(k+1) * Ep(k+1) * J_v(k+1,:);
end

S = sum(T,1);
interference_term = exp(- A * B * S);
P_cov = noise_term .* interference_term;

end