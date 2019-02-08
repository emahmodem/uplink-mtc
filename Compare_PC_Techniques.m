clear all; close all;
F = @(r,la,Po,Pm,a,b) 2 * pi * la .* r .* exp(- pi .* la .* r.^2) .* exp(- Po/Pm * exp(a * r .^ b));
%la = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1];
la = 1e-7:1e-7:1e-5;
la_ = [0.1 1 10 100 1000]*1e-6;
Po = 1e-13;   % -100 dBm
Pm = 0.1;     % 20 dbm
a = 0.3;
b = 2/3;

rm = log((Pm/Po).^(1/a)).^(1/b);
O0 = exp(- pi * rm^2 * la)
O0_ =exp(- pi * rm^2 * la_)
for p = 1:numel(la)
    Omin(p) = 1 - integral(@(r)F(r,la(p),Po,Pm,a,b),0,inf);
end


for p = 1:numel(la_)
    Omin_(p) = 1 - integral(@(r)F(r,la_(p),Po,Pm,a,b),0,inf) 
end
Omin_

f1 = semilogx(la*1e6,Omin,'-',la*1e6,O0,'--');
grid on;
set(f1,'MarkerSize',15);
set(f1,'LineWidth',4);
%legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','southeast','Interpreter','LaTex');
xlabel(' Small cell''s density  ($\lambda_s$ cells/km$^2$) ' ,'Interpreter','LaTex');
ylabel('Power truncation outage','Interpreter','LaTex');
title('$\alpha = 0.3, \beta = \frac{2}{3}$','Interpreter','LaTex')
set(gca, 'FontSize', 30);
set(gca, 'FontWeight', 'Bold');
set(gca, 'LineWidth', 2);
set(gca, 'GridAlpha', 0.5);
set(gca, 'MinorGridAlpha', 0.5);
legend('Power control technique 1' , 'Power control technique 2')