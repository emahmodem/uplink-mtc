Pu = 500e-3;
Po = 1e-3;
x = 0:0.0001:Pu;
LA_B =1e-3;
a = 4;
d = 2/a;
fx =  2 * pi * LA_B * x.^(d - 1).* exp(- pi * LA_B .* (x/Po).^d) / (a * Po^d .* (1 - exp(- pi * LA_B .* (Pu/Po).^d)));

x = Po:0.0001:Pu;
a = 0.94;
b = 0.5;
theta = log((Pu/Po).^(1/a)).^(1/b);
gx = (2 * pi * LA_B) / (a * b * (1 - exp(- pi * LA_B .* theta.^2))) * (1 ./ x) .* (log((x/Po).^(1/a)).^(2/b - 1)) .* exp(- pi * LA_B .* log((x/Po).^(1/a)).^(2/b));
F = @(x) (2 * pi * LA_B) / (a * b * (1 - exp(- pi * LA_B .* theta.^2))) * (1 ./ x) .* (log((x/Po).^(1/a)).^(2/b - 1)) .* exp(- pi * LA_B .* log((x/Po).^(1/a)).^(2/b));
I = integral(@(x)F(x),Po,Pu)
plot(x,gx)