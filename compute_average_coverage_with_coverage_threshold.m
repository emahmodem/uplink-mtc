function P_average = compute_average_coverage_with_coverage_threshold(params)
SNR = params.Pmin  / params.No;
noise_term = exp(- params.Threshold / SNR);

a = params.SEPL.alpha;
b = params.SEPL.beta;
n = (2/b) - 1;
t = params.Threshold;
Pu= params.Pm;
Po = params.Pmin;
theta = log((Pu/Po).^(1/a)).^(1/b);
F = @(y,v,tau) log(tau ./ y).^(v)  ./ ( y + 1);
A = (n+1) * (pi / a^(n+1)) * (params.LA_M / params.N_RB);
for k = 0:n
    Ep(k+1,:) =  a^k .* gamma (k/(n+1) + 1) .* gammainc( pi .* params.LA_B .* theta^2 , k/(n+1) + 1 , 'lower') ./ ((pi .* params.LA_B) .^(k/(n+1)) .* (1 - exp(- pi .* params.LA_B .* theta^2)));
    NK(k+1) = nchoosek(n,k);
    for p = 1:numel(t)
        J_v(k+1,p) = integral(@(y)F(y,n - k,t(p)),0,t(p));
    end
    T(k+1,:) = NK(k+1) * Ep(k+1) * J_v(k+1,:);
end

S = sum(T,1);
G = A * S;
interference_term = 1./(beta(3,4) .* G.^6) .* ( (6 * G.^2 + 48 * G + 120 ).* exp(-G) + (2 .* G.^3 - 18 .* G.^2 + 72 * G -120)) ;

P_average = noise_term .* interference_term;

end