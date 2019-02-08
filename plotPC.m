Po = 10e-13;
Pm = 0.1;
a = 0.3;
b = 2/3;

r = 0:0.1:1e3;


P1 = Po*exp(a *r .^b);
P1(P1 > Pm) = 0;
figure;
f1 = plot(r,P1);

grid on;
set(f1,'MarkerSize',15);
set(f1,'LineWidth',4);
%legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','southeast','Interpreter','LaTex');
xlabel(' Distance to serving cell  ($m$) ' ,'Interpreter','LaTex');
ylabel('Transmit power','Interpreter','LaTex');
title('$\alpha = 0.3, \beta = \frac{2}{3}$','Interpreter','LaTex')
set(gca, 'FontSize', 30);
set(gca, 'FontWeight', 'Bold');
set(gca, 'LineWidth', 2);
set(gca, 'GridAlpha', 0.5);
set(gca, 'MinorGridAlpha', 0.5);


P2 = min(Po*exp(a *r .^b),Pm);
figure;
f2 = plot(r,P2);

grid on;
set(f2,'MarkerSize',15);
set(f2,'LineWidth',4);
%legend(f1([1 4]),{'Analysis' , 'Simulation' },'FontSize',25,'FontWeight','bold','Location','southeast','Interpreter','LaTex');
xlabel(' Distance to serving cell  ($m$) ' ,'Interpreter','LaTex');
ylabel('Transmit power','Interpreter','LaTex');
title('$\alpha = 0.3, \beta = \frac{2}{3}$','Interpreter','LaTex')
set(gca, 'FontSize', 30);
set(gca, 'FontWeight', 'Bold');
set(gca, 'LineWidth', 2);
set(gca, 'GridAlpha', 0.5);
set(gca, 'MinorGridAlpha', 0.5);