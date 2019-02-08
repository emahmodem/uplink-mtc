hold on;
F = @(r,la,Po,Pm,a,b) 2 * pi * la .* r .* exp(- pi .* la .* r.^2) .* exp(- Po/Pm * exp(a * r .^ b));
%la = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1];
la = (1:10).*1e-6;
Po = 1e-11;   % -100 dBm
Pm = 0.1;     % 20 dbm
a = 0.94;
b = 1/2;


for p = 1:numel(la)
    Omin(p) = 1 - integral(@(r)F(r,la(p),Po,Pm,a,b),0,inf);
end



subplot(1,2,1);
hold on;
f1 = plot(la*1e6,Omin,'rs-');
grid on;
set(f1,'MarkerSize',10);
set(f1,'LineWidth',4);
