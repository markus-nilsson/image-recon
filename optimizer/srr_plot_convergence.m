function srr_plot_convergence(cgr)

y = [];
for c = 1:numel(cgr)
    y = cat(2, y, abs(cgr{c}));
end

semilogy(y, 'k.-');