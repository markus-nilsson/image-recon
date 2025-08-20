function m = imop_mask(d)
% assume we're getting dwi data

tmp = mean(d.w,2);

th = linspace(quantile(tmp(:), 0.01), quantile(tmp(:), 0.99), 100);

% Compute the fraction of voxels left at various thresholds
[n,x] = hist(tmp(:), th);
n = cumsum(n);
n = n / max(n);
n = 1 - n;


ind = x > th(round(0.5 * end));
p = polyfit(x(ind), n(ind), 1);

th_ind = find(n > (polyval(p, x) + 0.05), 1, 'last');

m = d.new(tmp > th(th_ind));
