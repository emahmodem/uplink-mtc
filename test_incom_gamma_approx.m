function  test_incom_gamma_approx()
clc;close all;clear all;
a = 1;
x = 0:0.01:10;
F_acc =  gammainc(x,a);
F_appr = 1 - exp(-x).* expo_partial(a,x);

plot(x,F_acc,'-r',x,F_appr,'.k');

legend({'Accurate','Approximate'})
end

function y = expo_partial(n,x)
y= zeros(1,numel(x));
for i = 0:n-1
    term = x.^i/factorial(i);
    y = y + term;
end
end