function [P_cov  P_cov_cf]= compute_uplink_coverage_with_coverage_threshold_closed_form(params)
SNR = params.Pmin  / params.No;
noise_term = exp(- params.Threshold / SNR);
a = params.SEPL.alpha;
b = params.SEPL.beta;
t = params.Threshold;
n = (2/b) - 1;
F = @(y,v,tau) log(tau ./ y).^(v)  ./ ( y + 1);
Pu= params.Pm;
Po = params.Pmin;
theta = log((Pu/Po).^(1/a)).^(1/b);
A = (n+1) * (pi / a^(n+1));
B = params.rho_m * (params.LA_M / params.N_RB);
for k = 0:n
    Ep(k+1) =  a^k .* gamma (k/(n+1) + 1) .* gammainc( pi .* params.LA_B .* theta^2 , k/(n+1) + 1 , 'lower') ./ ((pi .* params.LA_B) .^(k/(n+1)) .* (1 - exp(- pi .* params.LA_B .* theta^2)));
    NK(k+1) = nchoosek(n,k);
    for p = 1:numel(t)
        J_v(k+1,p) = integral(@(y)F(y,n - k,t(p)),0,t(p));
    end
    T(k+1,:) = NK(k+1) * Ep(k+1) * J_v(k+1,:);
end

S = sum(T,1);
interference_term = exp(- A * B * S);
P_cov = noise_term .* interference_term;


P_cov_cf = exp(- t/SNR + A*B*log(1+t)) ;

semilogx(t,P_cov,'r-',t,P_cov_cf,'r^');
end