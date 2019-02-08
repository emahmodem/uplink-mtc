t = 1000;
y = 0:0.001:t;
Po = 1e-7;
P = 100*Po;

f = log(t./y) ./ log(P/Po)

figure,
plot(y,f)