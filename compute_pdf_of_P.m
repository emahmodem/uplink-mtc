function [f_analy ,  x]= compute_pdf_of_P(params)
a = params.SEPL.alpha;
b = params.SEPL.beta;
Pu= params.Pm;
Po = params.Pmin;
LA_B = params.LA_B;
theta = log((Pu/Po).^(1/a)).^(1/b);
x = [linspace(params.Pmin,9*params.Pmin,10) linspace(10*params.Pmin,90*params.Pmin,10) linspace(100*params.Pmin,900*params.Pmin,10) linspace(1000*params.Pmin,params.Pm,10)  ];
%x = params.Pmin:params.Pmin:params.Pm;
f_analy = (2 * pi * LA_B) / (a * b * (1 - exp(- pi * LA_B .* theta.^2))) * (1 ./ x) .* (log((x/Po).^(1/a)).^(2/b - 1)) .* exp(- pi * LA_B .* log((x/Po).^(1/a)).^(2/b));
end