function  test_dilogarithm()
clc;close all;clear all;

x = -10:0.1:10;
F_acc =  dilog(1 + x) + dilog(1 + x.^-1);
F_appr = - ((pi^2)/6 + 1/2 *( log(x)).^2);

plot(x,F_acc,'-r',x,F_appr,'.k');

legend({'Accurate','Approximate'})
end

