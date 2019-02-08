data = randn(1000000,1);

figure(1);
[n , x] = hist(data,1000);
n_ = n./numel(data) * 100;
bar(x,n_);
hold on;
ksdensity(data);

